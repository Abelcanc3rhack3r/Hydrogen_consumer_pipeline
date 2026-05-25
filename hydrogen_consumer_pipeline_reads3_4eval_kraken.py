#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Hydrogen consumer pipeline:
  • Hydrogenases (HydDB): use external mapper to get post-filter read->subject, then join with Kraken2 taxonomy,
    aggregate by species/family/phylum × hydrogenase group, compute counts & RPKM.
  • Marker genes (user-specified FASTAs): do the same (DIAMOND + filters) per marker, then join & aggregate by
    species/family/phylum × marker, compute counts & RPKM.

Hyd-group extraction uses subject IDs containing "_-_[GROUP]" (e.g. FeFe_Group_A3):contentReference[oaicite:2]{index=2}.
Hydrogenase read→subject mapping CSV is produced by Diamond_blast_hyDB_reads_mapping.py:contentReference[oaicite:3]{index=3}.

Outputs:
  <out_dir>/taxon3_hydgroup_counts.csv
  <out_dir>/taxon3_hydgroup_rpkm.csv
  <out_dir>/taxon3_marker_counts.csv
  <out_dir>/taxon3_marker_rpkm.csv
"""

import os
import re
import csv
import sys
import argparse
import subprocess
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

# ---------- Defaults ----------
DEFAULT_MAPPING_SCRIPT = "/oceanstor/scratch/tllseedorf/e1103389/hyDB/Diamond_blast_hyDB_reads_mapping.py"
DEFAULT_DIAMOND_EXE    = "/oceanstor/home/e1103389/DIAMOND/diamond"
DEFAULT_HYDB_FASTA_CANDIDATES = [
    "/tllhome/abel/hydrogen_consumer_pipeline/hyDB/HydDB_reformated.fasta",
    "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/hyDB/HydDB_reformated.fasta",
    r"C:\Users\abel\Documents\hydrogenases\hyDB\HydDB_reformated.fasta",
]
DEFAULT_KRAKEN_EXE = "kraken2"
DEFAULT_THREADS    = 32

# Regex to pull hyd group like FeFe_Group_A3 from subject IDs formatted with "_-_[FeFe_Group_A3]"
GROUP_PAT = re.compile(r'_-_\[(.+?)\]')

# ---------- NCBI taxonomy (from Kraken DB: taxonomy/nodes.dmp, names.dmp) ----------
def load_ncbi_taxonomy(kraken_db: str):
    tax_dir = os.path.join(kraken_db, "taxonomy")
    nodes_path = os.path.join(tax_dir, "nodes.dmp")
    names_path = os.path.join(tax_dir, "names.dmp")
    if not (os.path.exists(nodes_path) and os.path.exists(names_path)):
        raise FileNotFoundError(f"Missing nodes.dmp/names.dmp under {tax_dir}")

    parent = {}
    rank = {}
    with open(nodes_path, "r") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            tax_id, parent_tax, r = parts[0], parts[1], parts[2]
            parent[tax_id] = parent_tax
            rank[tax_id]   = r

    sci_name = {}
    with open(names_path, "r") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            tax_id, name_txt, _, name_class = parts[0], parts[1], parts[2], parts[3]
            if name_class == "scientific name":
                sci_name[tax_id] = name_txt
    return parent, rank, sci_name

def lineage_levels(taxid: str, parent: dict, rank: dict, sci_name: dict):
    target = {"species": None, "family": None, "phylum": None}
    cur = taxid
    seen = set()
    while cur and cur not in seen:
        seen.add(cur)
        r = rank.get(cur, "")
        if r in target and target[r] is None:
            target[r] = (cur, sci_name.get(cur, ""))
            if all(target[k] is not None for k in target):
                break
        if cur == parent.get(cur):
            break
        cur = parent.get(cur)
    for k in target:
        if target[k] is None:
            target[k] = (None, "")
    return target

# ---------- Kraken2 ----------
def run_kraken2(reads1: str, kraken_db: str, output_path: str,
                kraken_exe: str, threads: int, paired: bool, reads2: str = None):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    report_path = os.path.splitext(output_path)[0] + ".report"
    cmd = [
        kraken_exe, "--db", kraken_db, "--threads", str(threads),
        "--report", report_path, "--output", output_path
    ]
    if paired:
        if not reads2:
            raise ValueError("paired=True but --reads2 is missing")
        cmd += ["--paired", reads1, reads2]
    else:
        cmd += ["--fasta-input", reads1]
    print("RUN KRAKEN2:", " ".join(cmd))
    subprocess.run(cmd, check=False)
    if not os.path.exists(output_path):
        raise FileNotFoundError(f"Kraken output not found: {output_path}")
    return output_path, report_path

def parse_kraken_output(kraken_output_path: str):
    read_to_tax = {}
    with open(kraken_output_path, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:  # status, read_id, taxid, ...
                continue
            if parts[0] == "C":
                read_to_tax[parts[1]] = parts[2]
    return read_to_tax

def fastq_count_reads(*fastq_paths):
    total = 0
    for fq in fastq_paths:
        if not fq:
            continue
        with open(fq, "r") as f:
            total += sum(1 for _ in f) // 4
    return total

# ---------- Hydrogenases (use external mapping script) ----------
def run_hydrogenase_mapper(reads_path: str, blast_dir: str, diamond_exe: str, mapping_script: str,
                           min_qcov: float, min_pid: float, min_scov: float, train_percent: float):
    os.makedirs(blast_dir, exist_ok=True)
    cmd = [
        sys.executable, mapping_script,
        reads_path, blast_dir, diamond_exe,
        "--min_qcov", str(min_qcov),
        "--min_pid",  str(min_pid),
        "--min_scov", str(min_scov),
        "--train",    str(train_percent),
    ]
    print("RUN HYD-MAP:", " ".join(cmd))
    subprocess.run(cmd, check=False)
    mapping_csv = os.path.join(blast_dir, "Filtered", "read_to_subject.csv")
    if not os.path.exists(mapping_csv):
        raise FileNotFoundError(f"Hydrogenase mapping CSV not found: {mapping_csv}")
    return mapping_csv

def extract_hyd_group(subject_id: str):
    m = GROUP_PAT.search(subject_id)
    return m.group(1) if m else None

def load_subject_lengths_nt_from_fasta(fasta_path: str):
    subj_len_nt = {}
    with open(fasta_path, "r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            subj_len_nt[rec.id] = len(rec.seq) * 3  # AA -> nt
    return subj_len_nt

def aggregate_hyd(mapping_csv: str, read_to_tax: dict, subject_len_nt: dict,
                  total_reads: int, parent: dict, rank: dict, sci_name: dict):
    total_reads = max(1, int(total_reads))
    counts_tax_subject = defaultdict(int)
    counts_lineage_group = defaultdict(int)

    with open(mapping_csv, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2:
                continue
            read_id, subject_id = row[0], row[1]
            taxid = read_to_tax.get(read_id)
            if not taxid:
                continue
            hyd = extract_hyd_group(subject_id)
            if not hyd:
                continue
            lv = lineage_levels(taxid, parent, rank, sci_name)
            sp_id, sp_nm = lv["species"]
            fa_id, fa_nm = lv["family"]
            ph_id, ph_nm = lv["phylum"]
            key = (sp_id, sp_nm, fa_id, fa_nm, ph_id, ph_nm, hyd)
            counts_lineage_group[key] += 1
            counts_tax_subject[(taxid, subject_id)] += 1

    # subject-level RPKM then roll up
    rpkm_tax_subject = {}
    for (taxid, subject_id), n in counts_tax_subject.items():
        glen = subject_len_nt.get(subject_id)
        if not glen or glen <= 0:
            continue
        rpkm_tax_subject[(taxid, subject_id)] = 1e9 * float(n) / (float(glen) * float(total_reads))

    rpkm_lineage_group = defaultdict(float)
    for (taxid, subject_id), val in rpkm_tax_subject.items():
        hyd = extract_hyd_group(subject_id)
        if not hyd:
            continue
        lv = lineage_levels(taxid, parent, rank, sci_name)
        sp_id, sp_nm = lv["species"]
        fa_id, fa_nm = lv["family"]
        ph_id, ph_nm = lv["phylum"]
        key = (sp_id, sp_nm, fa_id, fa_nm, ph_id, ph_nm, hyd)
        rpkm_lineage_group[key] += val

    return counts_lineage_group, rpkm_lineage_group

# ---------- Markers (run DIAMOND directly for each marker FASTA) ----------
def build_dmnd_if_needed(db_fasta: str, out_dmnd: str, diamond_exe: str):
    if os.path.exists(out_dmnd):
        return out_dmnd
    cmd = [diamond_exe, "makedb", "--in", db_fasta, "--db", out_dmnd]
    print("MAKE DMND:", " ".join(cmd))
    subprocess.run(cmd, check=False)
    if not os.path.exists(out_dmnd + ".dmnd"):
        # diamond writes "<out>.dmnd"
        raise FileNotFoundError(f"DIAMOND DB not created: {out_dmnd}.dmnd")
    return out_dmnd + ".dmnd"

def run_marker_diamond(reads_for_diamond: str, marker_name: str, marker_fasta: str,
                       out_dir: str, diamond_exe: str,
                       min_qcov: float, min_scov: float, min_pid: float):
    """
    Run DIAMOND blastx against marker FASTA, filter by coverage & identity,
    produce mapping CSV (read,subject) for that marker, and return:
      (mapping_csv_path, subject_len_nt_dict)
    """
    os.makedirs(out_dir, exist_ok=True)
    # Query lengths
    qext = "fasta" if reads_for_diamond.endswith((".fa", ".fna", ".fasta")) else "fastq"
    query_lengths = {}
    with open(reads_for_diamond, "r") as f:
        for rec in SeqIO.parse(f, qext):
            query_lengths[rec.id] = len(rec.seq)

    # Subject lengths (AA->nt)
    subj_len_nt = {}
    with open(marker_fasta, "r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            subj_len_nt[rec.id] = len(rec.seq) * 3

    # Ensure DMND
    dmnd_base = os.path.join(out_dir, f"{marker_name}")
    dmnd_db = build_dmnd_if_needed(marker_fasta, dmnd_base, diamond_exe)

    # Raw DIAMOND
    raw_tsv = os.path.join(out_dir, f"{marker_name}_raw.tsv")
    cmd = [
        diamond_exe, "blastx",
        "--db", dmnd_db,
        "--query", reads_for_diamond,
        "--out", raw_tsv,
        "--max-hsps", "1",
        "--max-target-seqs", "1",
        "--threads", "16",
        "--outfmt", "6"  # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    ]
    print("RUN MARKER:", " ".join(cmd))
    subprocess.run(cmd, check=False)

    # Prefilter by coverages
    pref_tsv = os.path.join(out_dir, f"{marker_name}_prefiltered.tsv")
    with open(raw_tsv, "r") as infile, open(pref_tsv, "w") as outfile:
        for line in infile:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            qid, sid = cols[0], cols[1]
            qstart, qend = int(cols[6]), int(cols[7])
            sstart, send = int(cols[8]), int(cols[9])
            qlen = query_lengths.get(qid)
            slen_nt = subj_len_nt.get(sid)  # nt
            if not qlen or not slen_nt:
                continue
            # Subject length is protein*3; DIAMOND coords are AA → compare in AA
            slen_aa = slen_nt // 3
            qcov = (abs(qend - qstart) + 1) / max(qlen, 1)
            scov = (abs(send - sstart) + 1) / max(slen_aa, 1)
            if scov >= float(min_scov) and qcov >= float(min_qcov):
                outfile.write(line)

    # Identity filter
    filt_tsv = os.path.join(out_dir, f"{marker_name}_filtered.tsv")
    with open(pref_tsv, "r") as infile, open(filt_tsv, "w") as outfile:
        for line in infile:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            try:
                pid = float(cols[2])
            except ValueError:
                continue
            if pid >= float(min_pid):
                outfile.write(line)

    # Mapping CSV (read,subject) — post-filter only
    mapping_csv = os.path.join(out_dir, f"{marker_name}_read_to_subject.csv")
    with open(filt_tsv, "r") as infile, open(mapping_csv, "w", newline="") as outcsv:
        w = csv.writer(outcsv)
        for line in infile:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                w.writerow([cols[0], cols[1]])

    return mapping_csv, subj_len_nt

def aggregate_markers(marker_map_list, read_to_tax: dict,
                      total_reads: int, parent: dict, rank: dict, sci_name: dict):
    """
    marker_map_list: list of tuples (marker_name, mapping_csv, subject_len_nt_dict)
    Returns two dicts keyed by (species_taxid,species,family_taxid,family,phylum_taxid,phylum,marker):
      counts_lineage_marker, rpkm_lineage_marker
    """
    total_reads = max(1, int(total_reads))
    counts_tax_subject = defaultdict(int)  # (taxid, subject_id) -> n
    subj_len_nt_all = {}
    # we will also keep marker label per subject_id
    subject_to_marker = {}

    # Read all marker mappings
    for marker_name, mapping_csv, subj_len_nt in marker_map_list:
        subj_len_nt_all.update(subj_len_nt)
        with open(mapping_csv, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) < 2:
                    continue
                read_id, subject_id = row[0], row[1]
                taxid = read_to_tax.get(read_id)
                if not taxid:
                    continue
                subject_to_marker[subject_id] = marker_name
                counts_tax_subject[(taxid, subject_id)] += 1

    # subject-level RPKM
    rpkm_tax_subject = {}
    for (taxid, subject_id), n in counts_tax_subject.items():
        glen = subj_len_nt_all.get(subject_id)
        if not glen or glen <= 0:
            continue
        rpkm_tax_subject[(taxid, subject_id)] = 1e9 * float(n) / (float(glen) * float(total_reads))

    # roll up to lineage × marker
    counts_lineage_marker = defaultdict(int)
    rpkm_lineage_marker   = defaultdict(float)

    for (taxid, subject_id), n in counts_tax_subject.items():
        marker_name = subject_to_marker.get(subject_id)
        if not marker_name:
            continue
        lv = lineage_levels(taxid, parent, rank, sci_name)
        sp_id, sp_nm = lv["species"]
        fa_id, fa_nm = lv["family"]
        ph_id, ph_nm = lv["phylum"]
        key = (sp_id, sp_nm, fa_id, fa_nm, ph_id, ph_nm, marker_name)
        counts_lineage_marker[key] += n

    for (taxid, subject_id), r in rpkm_tax_subject.items():
        marker_name = subject_to_marker.get(subject_id)
        if not marker_name:
            continue
        lv = lineage_levels(taxid, parent, rank, sci_name)
        sp_id, sp_nm = lv["species"]
        fa_id, fa_nm = lv["family"]
        ph_id, ph_nm = lv["phylum"]
        key = (sp_id, sp_nm, fa_id, fa_nm, ph_id, ph_nm, marker_name)
        rpkm_lineage_marker[key] += r

    return counts_lineage_marker, rpkm_lineage_marker

# ---------- Writers ----------
def write_counts_csv(path: str, rows: dict, last_col_name: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["species_taxid","species","family_taxid","family","phylum_taxid","phylum", last_col_name])
        for key, val in sorted(rows.items(),
                               key=lambda kv: (kv[0][2] or '', kv[0][4] or '', kv[0][0] or '', kv[0][-1] or '')):
            sp_id, sp_nm, fa_id, fa_nm, ph_id, ph_nm, label = key
            w.writerow([sp_id or "", sp_nm, fa_id or "", fa_nm, ph_id or "", ph_nm, f"{label}:{val}" if last_col_name=="hyd_group" else val])

def write_counts_csv_splitcol(path: str, rows: dict, label_col: str, value_col: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["species_taxid","species","family_taxid","family","phylum_taxid","phylum", label_col, value_col])
        for key, val in sorted(rows.items(),
                               key=lambda kv: (kv[0][2] or '', kv[0][4] or '', kv[0][0] or '', kv[0][-1] or '')):
            sp_id, sp_nm, fa_id, fa_nm, ph_id, ph_nm, label = key
            w.writerow([sp_id or "", sp_nm, fa_id or "", fa_nm, ph_id or "", ph_nm, label, val])

def write_rpkm_csv_splitcol(path: str, rows: dict, label_col: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["species_taxid","species","family_taxid","family","phylum_taxid","phylum", label_col, "rpkm"])
        for key, val in sorted(rows.items(),
                               key=lambda kv: (kv[0][2] or '', kv[0][4] or '', kv[0][0] or '', kv[0][-1] or '')):
            sp_id, sp_nm, fa_id, fa_nm, ph_id, ph_nm, label = key
            w.writerow([sp_id or "", sp_nm, fa_id or "", fa_nm, ph_id or "", ph_nm, label, f"{val:.6f}"])

# ---------- Main pipeline ----------
def pipeline(args):
    # Resolve HydDB FASTA (for hydrogenase subject lengths if needed later)
    hydb_fasta = args.hydb_fasta
    if not hydb_fasta:
        for cand in DEFAULT_HYDB_FASTA_CANDIDATES:
            if os.path.exists(cand):
                hydb_fasta = cand
                break
    if not hydb_fasta or not os.path.exists(hydb_fasta):
        raise FileNotFoundError("HydDB FASTA not found; provide --hydb_fasta")

    os.makedirs(args.tmp_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    # 0) Taxonomy
    parent, rank, sci_name = load_ncbi_taxonomy(args.kraken_db)

    # 1) Kraken2 (per-read taxonomy)
    kraken_out = os.path.join(args.tmp_dir, "kraken", "classified.tsv")
    run_kraken2(
        reads1=args.reads1, reads2=args.reads2,
        paired=args.paired, kraken_db=args.kraken_db,
        output_path=kraken_out, kraken_exe=args.kraken_exe,
        threads=args.kraken_threads
    )
    read_to_tax = parse_kraken_output(kraken_out)
    total_reads = fastq_count_reads(args.reads1, args.reads2 if args.paired else None)
    print(f"TOTAL READS (RPKM denominator): {total_reads}")

    # 2) HYDROGENASES: call external mapping script (produces read_to_subject.csv)
    hyd_blast_dir = os.path.join(args.tmp_dir, "hydblast")
    hyd_mapping_csv = run_hydrogenase_mapper(
        reads_path=args.reads_for_diamond, blast_dir=hyd_blast_dir,
        diamond_exe=args.diamond_exe, mapping_script=args.mapping_script,
        min_qcov=args.min_qcov, min_pid=args.min_pid, min_scov=args.min_scov,
        train_percent=args.train_percent
    )
    # Subject lengths (from HydDB FASTA)
    hyd_subject_len_nt = load_subject_lengths_nt_from_fasta(hydb_fasta)
    hyd_counts, hyd_rpkm = aggregate_hyd(
        mapping_csv=hyd_mapping_csv,
        read_to_tax=read_to_tax,
        subject_len_nt=hyd_subject_len_nt,
        total_reads=total_reads,
        parent=parent, rank=rank, sci_name=sci_name
    )
    # Write hydrogenase outputs
    hyd_counts_csv = os.path.join(args.out_dir, "taxon3_hydgroup_counts.csv")
    hyd_rpkm_csv   = os.path.join(args.out_dir, "taxon3_hydgroup_rpkm.csv")
    write_counts_csv_splitcol(hyd_counts_csv, hyd_counts, label_col="hyd_group", value_col="read_count")
    write_rpkm_csv_splitcol(hyd_rpkm_csv, hyd_rpkm, label_col="hyd_group")
    print("WROTE:", hyd_counts_csv)
    print("WROTE:", hyd_rpkm_csv)

    # 3) MARKERS: for each marker provided
    marker_map_list = []  # list of (marker_name, mapping_csv, subject_len_nt_dict)
    markers_dir = os.path.join(args.tmp_dir, "markers")
    if args.marker:
        for item in args.marker:
            if "=" not in item:
                raise ValueError(f"--marker expects name=path, got: {item}")
            name, fasta = item.split("=", 1)
            name = name.strip()
            fasta = fasta.strip()
            if not (name and os.path.exists(fasta)):
                raise FileNotFoundError(f"Marker missing or FASTA not found: {item}")
            out_dir = os.path.join(markers_dir, name)
            mapping_csv, subj_len_nt = run_marker_diamond(
                reads_for_diamond=args.reads_for_diamond,
                marker_name=name,
                marker_fasta=fasta,
                out_dir=out_dir,
                diamond_exe=args.diamond_exe,
                min_qcov=args.marker_min_qcov,
                min_scov=args.marker_min_scov,
                min_pid=args.marker_min_pid,
            )
            marker_map_list.append((name, mapping_csv, subj_len_nt))

    if marker_map_list:
        m_counts, m_rpkm = aggregate_markers(
            marker_map_list=marker_map_list,
            read_to_tax=read_to_tax,
            total_reads=total_reads,
            parent=parent, rank=rank, sci_name=sci_name
        )
        m_counts_csv = os.path.join(args.out_dir, "taxon3_marker_counts.csv")
        m_rpkm_csv   = os.path.join(args.out_dir, "taxon3_marker_rpkm.csv")
        write_counts_csv_splitcol(m_counts_csv, m_counts, label_col="marker", value_col="read_count")
        write_rpkm_csv_splitcol(m_rpkm_csv, m_rpkm, label_col="marker")
        print("WROTE:", m_counts_csv)
        print("WROTE:", m_rpkm_csv)
    else:
        print("No markers provided; skipping marker outputs.")

    return 0

# ---------- CLI ----------
def build_arg_parser():
    p = argparse.ArgumentParser(
        description="Hydrogenases (HydDB) + Kraken2 taxonomy + species/family/phylum hyd-group RPKM; and the same for user-specified marker genes."
    )
    # Kraken inputs
    p.add_argument("--reads1", required=True, help="FASTQ R1 for Kraken2.")
    p.add_argument("--reads2", default=None, help="FASTQ R2 for Kraken2 (if paired).")
    p.add_argument("--paired", action="store_true", help="Set if reads1/reads2 are paired.")
    # DIAMOND read input (could reuse reads1)
    p.add_argument("--reads_for_diamond", required=True, help="Reads file (FASTQ or FASTA) for DIAMOND.")
    # DBs & tools
    p.add_argument("--kraken_db", required=True, help="Kraken2 DB (must contain taxonomy/nodes.dmp & names.dmp).")
    p.add_argument("--mapping_script", default=DEFAULT_MAPPING_SCRIPT, help="Path to Diamond_blast_hyDB_reads_mapping.py")
    p.add_argument("--diamond_exe", default=DEFAULT_DIAMOND_EXE, help="Path to DIAMOND executable.")
    p.add_argument("--kraken_exe", default=DEFAULT_KRAKEN_EXE, help="Kraken2 executable.")
    p.add_argument("--kraken_threads", type=int, default=DEFAULT_THREADS, help="Threads for Kraken2.")
    p.add_argument("--hydb_fasta", default=None, help="HydDB reformatted FASTA (protein). If omitted, tries common locations.")
    # Output dirs
    p.add_argument("--out_dir", required=True, help="Final output directory.")
    p.add_argument("--tmp_dir", required=True, help="Working directory.")
    # HYD thresholds
    p.add_argument("--min_qcov", type=float, default=0.8, help="Hydrogenase DIAMOND min query coverage.")
    p.add_argument("--min_scov", type=float, default=0.0, help="Hydrogenase DIAMOND min subject coverage.")
    p.add_argument("--min_pid",  type=float, default=80.0, help="Hydrogenase DIAMOND min percent identity.")
    p.add_argument("--train_percent", type=float, default=100.0, help="Hydrogenase DB subsample percent.")
    # MARKER inputs: repeatable --marker name=/path/to/marker.fasta
    p.add_argument("--marker", action="append", help="Marker spec as name=/path/to/marker.fasta (repeatable).")
    # MARKER thresholds (separate knobs if you want them different from hyd)
    p.add_argument("--marker_min_qcov", type=float, default=0.8, help="Marker DIAMOND min query coverage.")
    p.add_argument("--marker_min_scov", type=float, default=0.0, help="Marker DIAMOND min subject coverage.")
    p.add_argument("--marker_min_pid",  type=float, default=80.0, help="Marker DIAMOND min percent identity.")
    return p

if __name__ == "__main__":
    parser = build_arg_parser()
    args = parser.parse_args()

    # Echo key args
    print(
        "ARGS:",
        f"reads1={args.reads1}",
        f"reads2={args.reads2}",
        f"paired={args.paired}",
        f"reads_for_diamond={args.reads_for_diamond}",
        f"kraken_db={args.kraken_db}",
        f"out_dir={args.out_dir}",
        f"tmp_dir={args.tmp_dir}",
        f"mapping_script={args.mapping_script}",
        f"diamond_exe={args.diamond_exe}",
        f"kraken_exe={args.kraken_exe}",
        f"hydb_fasta={args.hydb_fasta}",
        f"min_qcov={args.min_qcov}",
        f"min_scov={args.min_scov}",
        f"min_pid={args.min_pid}",
        f"train_percent={args.train_percent}",
        f"markers={args.marker}",
        f"marker_min_qcov={args.marker_min_qcov}",
        f"marker_min_scov={args.marker_min_scov}",
        f"marker_min_pid={args.marker_min_pid}",
        sep="\n  "
    )

    sys.exit(pipeline(args))
