#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Hydrogen consumer pipeline:
  • Hydrogenases (HydDB): use external mapper to get post-filter read->subject, then join with Kraken2 taxonomy,
    aggregate by species/family/phylum × hydrogenase group, compute counts & RPKM.
  • Marker genes (user-specified FASTAs): do the same (DIAMOND + filters) per marker, then join & aggregate by
    species/family/phylum × marker, compute counts & RPKM.

Hyd-group extraction uses subject IDs containing "_-_[GROUP]" (e.g. FeFe_Group_A3).
Hydrogenase read→subject mapping CSV is produced by Diamond_blast_hyDB_reads_mapping.py.

Outputs:
  <out_dir>/taxon3_hydgroup_counts.csv
  <out_dir>/taxon3_hydgroup_rpkm.csv
  <out_dir>/taxon3_marker_counts.csv
  <out_dir>/taxon3_marker_rpkm.csv
"""

from datetime import datetime
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

# ---------- NCBI taxonomy ----------
def _find_taxonomy_dir(kraken_db: str, explicit_tax_dir: str = None):
    """
    Resolve a directory that contains nodes.dmp and names.dmp.
    Precedence:
      1) explicit_tax_dir (if provided, and contains both files)
      2) <kraken_db>/taxonomy
      3) <kraken_db> itself (in case user points directly to a taxonomy dump dir)
    """
    candidates = []
    if explicit_tax_dir:
        candidates.append(explicit_tax_dir)
    if kraken_db:
        candidates.append(os.path.join(kraken_db, "taxonomy"))
        candidates.append(kraken_db)

    for d in candidates:
        nodes = os.path.join(d, "nodes.dmp")
        names = os.path.join(d, "names.dmp")
        if os.path.exists(nodes) and os.path.exists(names):
            return d
    # If nothing found, return first checked for error messaging
    return candidates[0] if candidates else None
def lineage_levels(taxid: str, parent: dict, rank: dict, sci_name: dict):
    """
    Walk up the taxonomy from a taxid to collect the nearest labels for
    species, family, and phylum. Returns a dict:
      {
        "species": (taxid or None, scientific_name or ""),
        "family" : (taxid or None, scientific_name or ""),
        "phylum" : (taxid or None, scientific_name or "")
      }
    Missing levels are returned as (None, "").
    """
    target = {"species": None, "family": None, "phylum": None}
    cur = taxid
    seen = set()
    while cur and cur not in seen:
        seen.add(cur)
        r = rank.get(cur, "")
        if r in target and target[r] is None:
            target[r] = (cur, sci_name.get(cur, ""))
            # stop early if we've filled all targets
            if all(target[k] is not None for k in target):
                break
        # reached root or self-loop
        if cur == parent.get(cur):
            break
        # step upward
        cur = parent.get(cur)

    # normalize any missing levels
    for k in target:
        if target[k] is None:
            target[k] = (None, "")
    return target

def load_ncbi_taxonomy_from_dir(tax_dir: str):
    """
    Load minimal taxonomy from a directory containing nodes.dmp and names.dmp.
    Returns:
      parent[taxid], rank[taxid], sci_name[taxid]
    """
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

# ---------- Kraken2 ----------
def run_kraken2(reads1: str, kraken_db: str, output_path: str,
                kraken_exe: str, threads: int):
    '''
    kraken2 \
  --db /oceanstor/scratch/tllseedorf/e1103389/krakendbb/db2/ \
  --output ERR2855939_kraken.names.out \
  --report ERR2855939_kraken.report \
  --use-names \
  ERR2855939_2.subset.fastq
    '''
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    report_path = os.path.splitext(output_path)[0] + ".report"
    cmd = [
        kraken_exe, "--db", "/oceanstor/scratch/tllseedorf/e1103389/krakendbb/db2/", "--threads", str(threads),
        "--report", report_path, "--output", output_path, reads1
    ]
    
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

def _detect_format(path: str):
    lower = path.lower()
    if lower.endswith((".fa", ".fna", ".fasta")):
        return "fasta"
    elif lower.endswith((".fq", ".fastq")):
        return "fastq"
    # Fallback: attempt FASTA first
    return "fasta"

def count_reads_generic(*paths):
    """
    Count sequences in FASTA or FASTQ using Biopython SeqIO.
    If both R1/R2 provided and paired, counts are summed (mates counted separately).
    """
    total = 0
    for p in paths:
        if not p:
            continue
        fmt = _detect_format(p)
        with open(p, "r") as handle:
            total += sum(1 for _ in SeqIO.parse(handle, fmt))
    return total

# ---------- Hydrogenases (use external mapping script) ----------
def run_hydrogenase_mapper(reads_path: str, blast_dir: str, diamond_exe: str, mapping_script: str,
                           min_qcov: float, min_pid: float, min_scov: float, train_percent: float, num_threads=16):
    os.makedirs(blast_dir, exist_ok=True)
    cmd = [
        sys.executable, mapping_script,
        reads_path, blast_dir, diamond_exe,
        "--min_qcov", str(min_qcov),
        "--min_pid",  str(min_pid),
        "--min_scov", str(min_scov),
        "--train",    str(train_percent),
        "--threads",  str(num_threads),
    ]
    print("RUN HYD-MAP:", " ".join(cmd))
    subprocess.run(cmd, check=False)
    mapping_csv = os.path.join(blast_dir, "Filtered", "read_to_subject.csv")
    if not os.path.exists(mapping_csv):
        raise FileNotFoundError(f"Hydrogenase mapping CSV not found: {mapping_csv}")
    return mapping_csv

# Replace the old GROUP_PAT/extract_hyd_group with these:
GROUP_BRACKET = re.compile(r'_-_\[(.*?)\]')
FAMILY_GROUP  = re.compile(r'(FeFe|NiFe)[-_ ]*(?:Group[_ -]*)?([A-Za-z0-9]+)')

def extract_hyd_group(subject_id: str):
    """
    Return labels like 'FeFe A3' / 'NiFe 1h' from HydDB subject IDs.
    Tries the bracket payload first (if present), then searches the whole ID.
    Falls back to plain 'FeFe'/'NiFe' only if no subgroup is detectable.
    """
    regex = r'\[.*'
    m = re.search(regex, subject_id)
    if m:
        #remove the _Group_ and []
        p=m.group(0).replace("_-_[", "").replace("]", "").replace("_Group_", " ").replace("-", " ").replace("_", " ")
        p=p.replace("[", "")
        return p
    return None
#gen_consumer_pipeline_reads3_4eval_kraken_fixed.py", line 253, in aggregate_hyd
#    hyd = extract_hyd_group(subject_id)
#  File "/oceanstor/scratch/tllseedorf/e1103389/hydrogen_consumer_pipeline/hydrogen_consumer_pipeline_reads3_4eval_kraken_fixed.py", line 226, in extract_hyd_group
##           ~~~~~~~^^^
#IndexError: no such group


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
            print("subject_id:", subject_id)
            taxid = read_to_tax.get(read_id)
            if not taxid:
                continue
            hyd = extract_hyd_group(subject_id)
            print("hyd:", hyd)
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
                       min_qcov: float, min_scov: float, min_pid: float, num_threads=16):
    """
    Run DIAMOND blastx against marker FASTA, filter by coverage & identity,
    produce mapping CSV (read,subject) for that marker, and return:
      (mapping_csv_path, subject_len_nt_dict)

    If thresholds_file is provided (TSV), it overrides min_pid/min_qcov/min_scov
    for this marker_name. Expected columns per row:
        name<TAB>pid<TAB>qcov[<TAB>scov]
    """
    import csv

    os.makedirs(out_dir, exist_ok=True)

    # (A) Query lengths
    qext = "fasta" if reads_for_diamond.endswith((".fa", ".fna", ".fasta")) else "fastq"
    query_lengths = {}
    with open(reads_for_diamond, "r") as f:
        for rec in SeqIO.parse(f, qext):
            query_lengths[rec.id] = len(rec.seq)

    # (B) Subject lengths (AA->nt)
    subj_len_nt = {}
    with open(marker_fasta, "r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            subj_len_nt[rec.id] = len(rec.seq) * 3
    thresholds_file = '/oceanstor/scratch/tllseedorf/e1103389/hydrogen_consumer_pipeline/marker_thresholds.csv'
    # (C) Optional per-marker thresholds override
    if thresholds_file and os.path.exists(thresholds_file):
        try:
            with open(thresholds_file, "r") as tf:
                reader = csv.reader(tf, delimiter="\t")
                for row in reader:
                    if not row:
                        continue
                    # accept name match (case-sensitive by default; tweak if you prefer .lower())
                    if row[0] == marker_name:
                        # row: name  pid  qcov  [scov]
                        if len(row) >= 2 and row[1]:
                            min_pid = float(row[1])
                        if len(row) >= 3 and row[2]:
                            min_qcov = float(row[2])
                        if len(row) >= 4 and row[3]:
                            min_scov = float(row[3])
                        print(f"[marker thresholds] {marker_name}: pid={min_pid} qcov={min_qcov} scov={min_scov}")
                        break
        except Exception as e:
            print(f"WARNING: could not read thresholds file '{thresholds_file}': {e}")

    # (D) Ensure DIAMOND DB
    dmnd_base = os.path.join(out_dir, f"{marker_name}")
    dmnd_db = build_dmnd_if_needed(marker_fasta, dmnd_base, diamond_exe)

    # (E) Raw DIAMOND
    raw_tsv = os.path.join(out_dir, f"{marker_name}_raw.tsv")
    cmd = [
        diamond_exe, "blastx",
        "--db", dmnd_db,
        "--query", reads_for_diamond,
        "--out", raw_tsv,
        "--max-hsps", "1",
        "--max-target-seqs", "1",
        "--threads", str(num_threads),
        "--outfmt", "6"  # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    ]
    print("RUN MARKER:", " ".join(cmd))
    subprocess.run(cmd, check=False)

    # (F) Prefilter by coverages (qcov & scov)
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
            # DIAMOND coords are AA for subject → compare in AA
            slen_aa = slen_nt // 3
            qcov = (abs(qend - qstart) + 1) / max(qlen, 1)
            scov = (abs(send - sstart) + 1) / max(slen_aa, 1)
            if scov >= float(min_scov) and qcov >= float(min_qcov):
                outfile.write(line)

    # (G) Identity filter
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

    # (H) Mapping CSV (read,subject) — post-filter only
    mapping_csv = os.path.join(out_dir, f"{marker_name}_read_to_subject.csv")
    with open(filt_tsv, "r") as infile, open(mapping_csv, "w", newline="") as outcsv:
        w = csv.writer(outcsv)
        for line in infile:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                w.writerow([cols[0], cols[1]])

    return mapping_csv, subj_len_nt


# ---------- Combined Markers (single DB + single blast, then split per marker) ----------
def read_fasta_subject_lengths_and_map(marker_name: str, marker_fasta: str):
    """
    Returns (subj_len_nt, subject_to_marker) for this single marker FASTA.
    subject_to_marker maps original subject ids -> marker_name (no renaming).
    """
    subj_len_nt = {}
    subj_to_marker = {}
    with open(marker_fasta, "r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            subj_len_nt[rec.id] = len(rec.seq) * 3
            subj_to_marker[rec.id] = marker_name
    return subj_len_nt, subj_to_marker

def build_combined_marker_db(marker_specs: list, combined_dir: str, diamond_exe: str):
    """
    marker_specs: list of (name, fasta_path)
    Builds one combined FASTA and a single DIAMOND DB.
    Returns:
      dmnd_db_path, subj_len_nt_all, subject_to_marker
    """
    os.makedirs(combined_dir, exist_ok=True)
    combined_fa = os.path.join(combined_dir, "all_markers.faa")

    subj_len_nt_all = {}
    subject_to_marker = {}

    # Write combined FASTA without changing IDs; keep an index of subject->marker.
    with open(combined_fa, "w") as outfa:
        for name, fasta in marker_specs:
            if not os.path.exists(fasta):
                raise FileNotFoundError(f"Marker FASTA not found: {fasta}")
            slen_nt, s2m = read_fasta_subject_lengths_and_map(name, fasta)
            subj_len_nt_all.update(slen_nt)
            subject_to_marker.update(s2m)
            for rec in SeqIO.parse(fasta, "fasta"):
                SeqIO.write(rec, outfa, "fasta")

    dmnd_base = os.path.join(combined_dir, "all_markers")
    dmnd_db = build_dmnd_if_needed(combined_fa, dmnd_base, diamond_exe)
    return dmnd_db, subj_len_nt_all, subject_to_marker

def run_single_markers_blast_and_split(reads_for_diamond: str,
                                       dmnd_db: str,
                                       out_dir: str,
                                       diamond_exe: str,
                                       subject_to_marker: dict,
                                       subj_len_nt_all: dict,
                                       default_min_qcov: float,
                                       default_min_scov: float,
                                       default_min_pid: float,
                                       threads: int):
    """
    One blastx against the combined DB, then apply per-marker thresholds and
    split to per-marker mapping CSVs. Returns marker_map_list compatible with
    aggregate_markers(): [(marker_name, mapping_csv, subj_len_nt_dict), ...]
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1) Query lengths (FASTA/FASTQ)
    qext = "fasta" if reads_for_diamond.endswith((".fa",".fna",".fasta")) else "fastq"
    query_lengths = {}
    with open(reads_for_diamond, "r") as f:
        for rec in SeqIO.parse(f, qext):
            query_lengths[rec.id] = len(rec.seq)

    # 2) Single DIAMOND run
    raw_tsv = os.path.join(out_dir, "all_markers_raw.tsv")
    cmd = [
        diamond_exe, "blastx",
        "--db", dmnd_db,
        "--query", reads_for_diamond,
        "--out", raw_tsv,
        "--max-hsps", "1",
        "--max-target-seqs", "1",
        "--threads", str(threads),
        "--outfmt", "6"
    ]
    #add the timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"RUN MARKERS (single blast) at {timestamp}:", " ".join(cmd))
    subprocess.run(cmd, check=False)

    # 3) Optional thresholds TSV (same location as original code)
    thresholds_file = '/oceanstor/scratch/tllseedorf/e1103389/hydrogen_consumer_pipeline/marker_thresholds.csv'
    per_marker_thr = {}  # marker -> (pid, qcov, scov)
    if os.path.exists(thresholds_file):
        try:
            with open(thresholds_file, "r") as tf:
                r = csv.reader(tf, delimiter="\t")
                #skip the header
                next(r, None)
                for row in r:
                    if not row:
                        continue
                    m = row[0]
                    pid = float(row[1]) if len(row)>=2 and row[1] else default_min_pid
                    qcv = float(row[2]) if len(row)>=3 and row[2] else default_min_qcov
                    scv = float(row[3]) if len(row)>=4 and row[3] else default_min_scov
                    per_marker_thr[m] = (pid, qcv, scv)
        except Exception as e:
            print(f"WARNING: could not read thresholds file '{thresholds_file}': {e}")

    # 4) Prepare writers per marker
    per_marker_files = {}  # marker -> mapping_csv_path
    per_marker_writers = {}
    per_marker_tsv = {}    # marker -> filtered TSV handle

    def _ensure_marker_io(marker):
        if marker in per_marker_writers:
            return
        mdir = os.path.join(out_dir, marker)
        os.makedirs(mdir, exist_ok=True)
        filt_tsv = os.path.join(mdir, f"{marker}_filtered.tsv")
        map_csv  = os.path.join(mdir, f"{marker}_read_to_subject.csv")
        per_marker_tsv[marker] = open(filt_tsv, "w")
        per_marker_files[marker] = map_csv
        per_marker_writers[marker] = csv.writer(open(map_csv, "w", newline=""))

    # 5) Stream raw hits, compute coverage/identity filters per-marker, and write
    with open(raw_tsv, "r") as infile:
        for line in infile:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            qid, sid = cols[0], cols[1]
            marker = subject_to_marker.get(sid)
            if not marker:
                continue
            pid, qcv, scv = per_marker_thr.get(marker, (default_min_pid, default_min_qcov, default_min_scov))

            # coverage calcs
            try:
                qstart, qend = int(cols[6]), int(cols[7])
                sstart, send = int(cols[8]), int(cols[9])
            except ValueError:
                continue
            qlen = query_lengths.get(qid)
            slen_nt = subj_len_nt_all.get(sid)
            if not qlen or not slen_nt:
                continue
            slen_aa = slen_nt // 3
            qcov = (abs(qend - qstart) + 1) / max(qlen, 1)
            scov = (abs(send - sstart) + 1) / max(slen_aa, 1)

            # identity
            try:
                pident = float(cols[2])
            except ValueError:
                continue

            if (scov >= scv) and (qcov >= qcv) and (pident >= pid):
                _ensure_marker_io(marker)
                per_marker_tsv[marker].write(line)
                per_marker_writers[marker].writerow([qid, sid])

    # close TSV handles
    for fh in per_marker_tsv.values():
        try:
            fh.close()
        except Exception:
            pass

    # 6) Return list expected by aggregate_markers()
    out = []
    for marker, map_csv in per_marker_files.items():
        slen_nt_subset = {sid: ln for sid, ln in subj_len_nt_all.items()
                          if subject_to_marker.get(sid) == marker}
        out.append((marker, map_csv, slen_nt_subset))
    return out


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

# ---------- Writers with taxon-name mode ----------
def _format_lineage_columns(key, mode: str):
    """
    key = (species_taxid, species_name, family_taxid, family_name, phylum_taxid, phylum_name, label)
    mode: 'both' | 'names' | 'ids'
    Returns: (cols, label) list for the lineage columns and the label (hyd_group/marker)
    """
    sp_id, sp_nm, fa_id, fa_nm, ph_id, ph_nm, label = key
    if mode == "names":
        return [sp_nm, fa_nm, ph_nm], label
    elif mode == "ids":
        return [sp_id or "", fa_id or "", ph_id or ""], label
    else:
        return [sp_id or "", sp_nm, fa_id or "", fa_nm, ph_id or "", ph_nm], label

def write_counts_csv_splitcol(path: str, rows: dict, label_col: str, value_col: str, mode: str = "both"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        if mode == "names":
            header = ["species","family","phylum", label_col, value_col]
        elif mode == "ids":
            header = ["species_taxid","family_taxid","phylum_taxid", label_col, value_col]
        else:
            header = ["species_taxid","species","family_taxid","family","phylum_taxid","phylum", label_col, value_col]
        w.writerow(header)

        for key, val in sorted(rows.items(),
                               key=lambda kv: (kv[0][2] or '', kv[0][4] or '', kv[0][0] or '', kv[0][-1] or '')):
            cols, label = _format_lineage_columns(key, mode)
            w.writerow(cols + [label, val])

def write_rpkm_csv_splitcol(path: str, rows: dict, label_col: str, mode: str = "both"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        if mode == "names":
            header = ["species","family","phylum", label_col, "rpkm"]
        elif mode == "ids":
            header = ["species_taxid","family_taxid","phylum_taxid", label_col, "rpkm"]
        else:
            header = ["species_taxid","species","family_taxid","family","phylum_taxid","phylum", label_col, "rpkm"]
        w.writerow(header)

        for key, val in sorted(rows.items(),
                               key=lambda kv: (kv[0][2] or '', kv[0][4] or '', kv[0][0] or '', kv[0][-1] or '')):
            cols, label = _format_lineage_columns(key, mode)
            w.writerow(cols + [label, f"{val:.6f}"])

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

    # 0) Taxonomy (allow explicit taxonomy dir override)
    tax_dir = _find_taxonomy_dir(args.kraken_db, args.kraken_taxonomy_dir)
    if not tax_dir:
        raise FileNotFoundError("Unable to resolve a taxonomy directory (nodes.dmp/names.dmp). "
                                "Provide --kraken_taxonomy_dir or ensure <kraken_db>/taxonomy exists.")
    print(f"Using taxonomy dir: {tax_dir}")
    parent, rank, sci_name = load_ncbi_taxonomy_from_dir(tax_dir)

    # 1) Kraken2 (per-read taxonomy) - using FASTA input explicitly
    kraken_out = os.path.join(args.tmp_dir, "kraken", "classified.tsv")
    run_kraken2(
        reads1=args.reads1,  kraken_db=args.kraken_db,
        output_path=kraken_out, kraken_exe=args.kraken_exe,
        threads=args.kraken_threads
    )
    read_to_tax = parse_kraken_output(kraken_out)

    # 2) Total reads for RPKM denominator (support FASTA or FASTQ)
    total_reads = count_reads_generic(args.reads1)
    print(f"TOTAL READS (RPKM denominator): {total_reads}")

    # 3) HYDROGENASES: call external mapping script (produces read_to_subject.csv)
    hyd_blast_dir = os.path.join(args.tmp_dir, "hydblast")
    hyd_mapping_csv = run_hydrogenase_mapper(
        reads_path=args.reads_for_diamond, blast_dir=hyd_blast_dir,
        diamond_exe=args.diamond_exe, mapping_script=args.mapping_script,
        min_qcov=args.min_qcov, min_pid=args.min_pid, min_scov=args.min_scov,
        train_percent=args.train_percent, num_threads=args.kraken_threads
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
    write_counts_csv_splitcol(hyd_counts_csv, hyd_counts, label_col="hyd_group", value_col="read_count", mode=args.taxon_name_mode)
    write_rpkm_csv_splitcol(hyd_rpkm_csv, hyd_rpkm, label_col="hyd_group", mode=args.taxon_name_mode)
    print("WROTE:", hyd_counts_csv)
    print("WROTE:", hyd_rpkm_csv)
    #add the timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"Completed hydrogenase aggregation at {timestamp}")

    
    # 4) MARKERS: single combined DB & one DIAMOND run, then split/aggregate
    marker_specs = [
        ("mcrA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/methanogen/mcrA.fasta"),
        ("acsB",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/acetogen/acsB.fasta"),
        ("dsrA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/sulphate_reducer/dsrA.fasta"),
        ("aprA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/sulphate_reducer/aprA.fasta"),
        ("asrA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/sulphate_reducer/asrA.fasta"),
        ("dmsA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/trimethylamine_N_oxide_reduction/DmsA.fasta"),
        ("napA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/nitrate_reducer/napA.fasta"),
        ("narG",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/nitrate_reducer/narG.fasta"),
        ("nrfA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/nitrate_reducer/nrfA.fasta"),
        ("frdA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/fumarate_reducer/frdA.fasta"),
        ("cydA",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/cydA.fasta"),
        ("fhl",        "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/hyd_production/FHL.fasta"),
        ("pfor",       "/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/hyd_production/PFOR.fasta"),
        ("nitrogenase","/oceanstor/scratch/tllseedorf/e1103389/hydrogenases/hyd_production/nifH.fasta"),
    ]
    markers_dir = os.path.join(args.tmp_dir, "markers_singleblast")
    os.makedirs(markers_dir, exist_ok=True)
    dmnd_db, subj_len_nt_all, subject_to_marker = build_combined_marker_db(
        marker_specs, markers_dir, args.diamond_exe
    )
    marker_map_list = run_single_markers_blast_and_split(
        reads_for_diamond=args.reads_for_diamond,
        dmnd_db=dmnd_db,
        out_dir=markers_dir,
        diamond_exe=args.diamond_exe,
        subject_to_marker=subject_to_marker,
        subj_len_nt_all=subj_len_nt_all,
        default_min_qcov=args.marker_min_qcov,
        default_min_scov=args.marker_min_scov,
        default_min_pid=args.marker_min_pid,
        threads=args.kraken_threads
    )
    if marker_map_list:
        m_counts, m_rpkm = aggregate_markers(
            marker_map_list=marker_map_list,
            read_to_tax=read_to_tax,
            total_reads=total_reads,
            parent=parent, rank=rank, sci_name=sci_name
        )
        m_counts_csv = os.path.join(args.out_dir, "taxon3_marker_counts.csv")
        m_rpkm_csv   = os.path.join(args.out_dir, "taxon3_marker_rpkm.csv")
        write_counts_csv_splitcol(m_counts_csv, m_counts, label_col="marker", value_col="read_count", mode=args.taxon_name_mode)
        write_rpkm_csv_splitcol(m_rpkm_csv, m_rpkm, label_col="marker", mode=args.taxon_name_mode)
        print("WROTE:", m_counts_csv)
        print("WROTE:", m_rpkm_csv)
        #add the timestamp
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"Completed marker aggregation at {timestamp}")
    else:
        print("No markers provided; skipping marker outputs.")

    return 0

# ---------- CLI ----------
def build_arg_parser():
    p = argparse.ArgumentParser(
        description="Hydrogenases (HydDB) + Kraken2 taxonomy + species/family/phylum hyd-group RPKM; and the same for user-specified marker genes."
    )
    # Kraken inputs
    p.add_argument("--reads1", required=True, help="FASTA/FASTQ R1 for Kraken2.")
    # DIAMOND read input (could reuse reads1)
    p.add_argument("--reads_for_diamond", required=True, help="Reads file (FASTQ or FASTA) for DIAMOND.")
    # DBs & tools
    p.add_argument("--kraken_db", required=True, help="Kraken2 DB (point this at the DB root used by kraken2 --db).",default='/oceanstor/scratch/tllseedorf/e1103389/krakendbb/')
    p.add_argument("--kraken_taxonomy_dir", default='/oceanstor/scratch/tllseedorf/e1103389/krakendbb/taxonomy',
                   help="(Optional) Directory that contains nodes.dmp and names.dmp. "
                        "Use if different from <kraken_db>/taxonomy or if taxonomy lives elsewhere.")
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
    #p.add_argument("--marker", action="append", help="Marker spec as name=/path/to/marker.fasta (repeatable).")
    # MARKER thresholds (separate knobs if you want them different from hyd)
    p.add_argument("--marker_min_qcov", type=float, default=0.8, help="Marker DIAMOND min query coverage.")
    p.add_argument("--marker_min_scov", type=float, default=0.0, help="Marker DIAMOND min subject coverage.")
    p.add_argument("--marker_min_pid",  type=float, default=80.0, help="Marker DIAMOND min percent identity.")
    # Taxon-name expansion mode
    p.add_argument("--taxon_name_mode", choices=["both","names","ids"], default="both",
                   help="How to print taxon columns in outputs: both (default), names, or ids.")
    return p

if __name__ == "__main__":
    parser = build_arg_parser()
    args = parser.parse_args()

    # Echo key args
    print(
        "ARGS:",
        f"reads1={args.reads1}",
        f"reads_for_diamond={args.reads_for_diamond}",
        f"kraken_db={args.kraken_db}",
        f"kraken_taxonomy_dir={args.kraken_taxonomy_dir}",
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
        f"marker_min_qcov={args.marker_min_qcov}",
        f"marker_min_scov={args.marker_min_scov}",
        f"marker_min_pid={args.marker_min_pid}",
        f"taxon_name_mode={args.taxon_name_mode}",
        sep="\n  "
    )

    sys.exit(pipeline(args))
