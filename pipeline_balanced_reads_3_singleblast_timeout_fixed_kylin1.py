#!/usr/bin/env python3
import os
import subprocess
import shutil
from pathlib import Path
import zipfile
import re
import sys
import random
import pickle
from collections import Counter
from datetime import datetime
#conda activate shotgun first!
# -------- NEW: central config for hydrogen pipeline --------
HYDROGEN_SCRIPT = "/oceanstor/scratch/tllseedorf/e1103389/hydrogen_consumer_pipeline/hydrogen_consumer_pipeline_reads_3_4_eval_kraken_singleblast_unclassified1.py"
#hydrogen_consumer_pipeline_reads3_4eval_kraken_singleblast.py"
# Use the UHGG v2.0.2 Kraken2 DB
KRAKEN_DB = "/oceanstor/scratch/tllseedorf/e1103389/human_gut/kraken_db"

HYDB_FASTA      = "/oceanstor/scratch/tllseedorf/e1103389/hyDB/HydDB_reformated.fasta"
DIAMOND_EXE     = "diamond"  # assume in PATH
MAPPING_SCRIPT  = "/oceanstor/scratch/tllseedorf/e1103389/hyDB/Diamond_blast_hyDB_reads_mapping.py"
KRAKEN_THREADS  = 25
MIN_QCOV        = 0.8
MIN_SCOV        = 0.0
MIN_PID         = 80
TRAIN_PERCENT   = 100
MARKERS                 = []
MARKER_THRESHOLDS_FILE  = None
TIMEOUT_SEC = 1800  # example: 30 minutes
# Male and female directories
MALE_DIR = Path("/oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/reads3/human_male_healthy/successful/fasta")
FEMALE_DIR = Path("/oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/reads3/human_female_healthy")

# Quotas
TARGET_MALE   = 50
TARGET_FEMALE = 50

from datetime import datetime  # (already imported below; safe to keep once)

RUN_TAG = "run2"  # e.g., uhgg_20251016_1234

BASE_WORKING_DIR = Path("/oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/reads3/qc_and_trim")
BASE_RESULT_DIR  = Path(f"/oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/reads3/hydrogen_results/{RUN_TAG}")
BASE_TMP_DIR     = Path(f"/oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/reads3/tmp/hydrogen_pipe/{RUN_TAG}")

BASE_WORKING_DIR.mkdir(exist_ok=True)
BASE_RESULT_DIR.mkdir(exist_ok=True)
BASE_TMP_DIR.mkdir(exist_ok=True)


#read the PROBLEM_LIST from a file
PROBLEM_LIST_FILE = "/oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/reads3/problem_list.txt"
if os.path.exists(PROBLEM_LIST_FILE):
    with open(PROBLEM_LIST_FILE, "r") as f:
        PROBLEM_LIST = [line.strip() for line in f if line.strip()]
else:
    PROBLEM_LIST = []
print("Loaded PROBLEM_LIST with", len(PROBLEM_LIST), "entries.")
print(PROBLEM_LIST)
def log_progress(logfile, message):
    from datetime import datetime
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(logfile, "a") as f:
        f.write(f"[{timestamp}] {message}\n")



def run_command(cmd, dry_run=False, capture_output=False, timeout_sec=None, output_file=None, offending_file=None):
    """
    Run a shell command with an overall timeout.

    Args:
        cmd (str): Command to execute.
        dry_run (bool): If True, only print the command.
        capture_output (bool): If True, return stdout as text.
        timeout_sec (int|float|None): Timeout in seconds. If None, uses TIMEOUT_SEC.
        output_file (str|Path|None): If provided and capture_output=False, stdout is written here.

    Returns:
        str | int: stdout (str) if capture_output=True, else the return code (int).

    Raises:
        subprocess.CalledProcessError: Non-zero exit code.
        subprocess.TimeoutExpired: Command timed out.
    """
    print(f"[CMD] {cmd}")
    if dry_run:
        return "" if capture_output else 0

    if timeout_sec is None:
        timeout_sec = TIMEOUT_SEC

    try:
        if capture_output:
            # Return stdout as text
            return subprocess.check_output(
                cmd, shell=True, text=True, timeout=timeout_sec
            )
        else:
            # Optionally capture to file
            if output_file is not None:
                proc = subprocess.run(
                    cmd, shell=True, text=True, capture_output=True, timeout=timeout_sec
                )
                # Write stdout to file
                Path(output_file).write_text(proc.stdout or "", encoding="utf-8")
                if proc.returncode != 0:
                    raise subprocess.CalledProcessError(
                        proc.returncode, cmd, output=proc.stdout, stderr=proc.stderr
                    )
                return proc.returncode
            else:
                proc = subprocess.run(cmd, shell=True, timeout=timeout_sec)
                if proc.returncode != 0:
                    raise subprocess.CalledProcessError(proc.returncode, cmd)
                return proc.returncode

    except subprocess.TimeoutExpired as e:
        # Surface a clear timeout error with partial outputs if any
        msg = f"Command timed out after {timeout_sec}s: {cmd}"
        #print
        #add the timestamp
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {msg}")
        #add the offending file to the problem list if provided
        if offending_file is not None:
            PROBLEM_LIST.append(offending_file)
            #log_progress("problem_list.txt", f"Timeout on {offending_file}")
            #save the problem list to a file
            with open(PROBLEM_LIST_FILE, "w") as f:
                for item in PROBLEM_LIST:
                    f.write(f"{item}\n")
        return None
    
# --------------------------------------------------------------------------------------
# Keep all the existing helper functions from your previous script (not reprinted here):
# - log_progress, custom_print, parse_fastqc_and_adapters, run_command, get_fastqc_info,
# - parse_wget, download_fastq, output_file_exists, process_one_fastq, etc.
# --------------------------------------------------------------------------------------
def parse_fastqc_and_adapters(filename_prefix, cwd=".", phred_thresh=30, adapter_threshold=0.01):
    zip_path = Path(cwd) / f"{filename_prefix}_fastqc.zip"
    fastqc_data_path = f"{filename_prefix}_fastqc/fastqc_data.txt"

    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(fastqc_data_path) as fastqc_file:
            content = fastqc_file.read().decode(errors="ignore")

    # Per-base Quality
    mean_quals = []
    parsing = False
    for line in content.splitlines():
        if line.startswith(">>Per base sequence quality"):
            parsing = True
            continue
        if parsing and line.startswith(">>END_MODULE"):
            break
        if parsing and not (line.startswith("#") or not line.strip()):
            fields = re.split(r"\s+", line.strip())
            # Support for base ranges
            if "-" in fields[0]:
                start, end = map(int, fields[0].split("-"))
            else:
                start = end = int(fields[0])
            mean = float(fields[1])
            for _ in range(start, end+1):
                mean_quals.append(mean)

    # Trim length = where phred drops below thresh, and most after are < thresh
    trim_len = None
    for i, q in enumerate(mean_quals):
        rest = mean_quals[i:]
        if q < phred_thresh and sum(rq < phred_thresh for rq in rest) >= len(rest) / 2:
            trim_len = i
            break
    if trim_len is None:
        trim_len = len(mean_quals)

    # Trunc length = position where all previous are < thresh
    trunc_len = 0
    for i in range(len(mean_quals)):
        if all(q < phred_thresh for q in mean_quals[:i+1]):
            trunc_len = i+1

    # Adapter content parsing
    adapter_content = {}  # {adapter_name: [fraction,...]}
    adapters = []
    parsing = False
    for line in content.splitlines():
        if line.startswith(">>Adapter Content"):
            parsing = True
            continue
        if parsing and line.startswith(">>END_MODULE"):
            break
        if parsing:
            if line.startswith("#Position"):
                fields = re.split(r"\s+",line.strip())[1:]
                adapters = ["Illumina Universal Adapter", "Illumina Small RNA 3' Adapter",
                            "Illumina Small RNA 5' Adapter", "Nextera Transposase Sequence", "PolyA", "PolyG"]
                for ad in adapters:
                    adapter_content[ad] = []
                continue
            if re.match(r"^\d", line):
                vals = re.split(r"\s+", line.strip())
                for i, ad in enumerate(adapters):
                    adapter_content[ad].append(float(vals[i+1]))

    # Adapters with >1% at any pos
    bad_adapters = []
    for adapter, abunds in adapter_content.items():
        if any(a > adapter_threshold for a in abunds):
            bad_adapters.append(adapter)

    # Output as tsv
    print(f"filename\ttrim_len\ttrunc_len\tbad_adapters")
    print(f"{filename_prefix}\t{trim_len}\t{trunc_len}\t{','.join(bad_adapters) if bad_adapters else 'None'}")

    # Also return for use in pipeline
    return {
        "trim_len": trim_len,
        "trunc_len": trunc_len,
        "adapter_content": adapter_content,
        "bad_adapters": bad_adapters
    }
def process_one_fastq(
    fq_path, working_dir, result_dir, tmp_dir, do_cleanup=True, cache=None,):

    fq = Path(fq_path)
    #if fq is a fasta file, skip to the pipeline
    if not fq.suffix in ['.fasta', '.fa', '.fna', '.fasta.gz']:
        
        fq_basename = fq.stem.replace('.fastq', '').replace('.fq', '').replace('.gz','')
        #if the abort counter is > 5, skip the file
        if cache is not None and cache["abort_counter"][fq_basename] > 5:
            print(f"[SKIP] {fq_basename} skipped due to too many previous errors.")
            log_progress("fastqc_error.log", f"{fq_basename} skipped due to too many previous errors.")
            #add to the errors
            cache["errors"].append(fq_basename)
            return

        fastqc_zip = Path(working_dir) / f"{fq_basename}_fastqc.zip"
        trimmed_fq = Path(working_dir) / f"{fq_basename}.trimmed.fastq.gz"
        fasta_file = Path(working_dir) / f"{fq_basename}.fasta"
        tmp_fasta = Path(tmp_dir) / fasta_file.name
        bad_adapters= None
        # 1. FASTQC
        for tries in range(3):
            print("tries:", tries)
            if not fastqc_zip.exists():
                #try multiple times to run fastqc
                print(f"[INFO] FastQC not produced for {fq_basename}, running now...")
                try:
                    outp=run_command(f"fastqc -o {working_dir} {fq}",capture_output=True)
                    #save output to a file
                    with open(f"{working_dir}/fastqc_output.log", "w") as f:
                        f.write(outp)
                except subprocess.CalledProcessError as e:
                    print(f"[ERROR] FastQC failed: {e}")
                    #add 1 to the cache abort counter
                    if cache is not None:
                        cache["abort_counter"][fq_basename] += 1
                        print(f"[ERROR] Cache abort counter for {fq_basename} incremented. It is now {cache['abort_counter'][fq_basename]}")
                    return
            else:
                print("[SKIP] FastQC already produced.")
                log_progress("fastqc.log", f"FastQC already produced for {fq_basename}, skipping.")
                #import pdb; pdb.set_trace()
                print("FASTQC output:", fastqc_zip)
                # 2. Parse FastQC
            try:
                qc_res = parse_fastqc_and_adapters(fq_basename, cwd=working_dir)
                trim_len, trunc_len, bad_adapters = qc_res['trim_len'], qc_res['trunc_len'], qc_res['bad_adapters']
                #add the timestamp
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[INFO] FastQC results for {fq_basename} at {timestamp}:")
                print(f"[INFO] trim_len={trim_len}, trunc_len={trunc_len}, bad_adapters={bad_adapters}")
            except Exception as e:
                print(f"[ERROR] Failed to parse FastQC: {e}, deleting fastqc zip file.")
                log_progress("fastqc_error.log", f"Failed to parse FastQC for {fq_basename}: {e}")
                #delete the fastqc zip file
                #add 1 to the cache abort counter
                #if cache is not None:
                    #cache["abort_counter"][fq_basename] += 1
                os.remove(fastqc_zip)
            break
        else:
            print(f"[ERROR] FastQC failed after 3 attempts for {fq_basename}.")
            log_progress("fastqc_error.log", f"FastQC failed after 3 attempts for {fq_basename}.")
            #add 1 to the cache abort counter
            if cache is not None:
                cache["abort_counter"][fq_basename] += 1
            return
        # 3. Cutadapt
        if not trimmed_fq.exists():
            #print info, cutadapt not produced
            print(f"[INFO] fastp not produced for {fq_basename}, running now...")
            fastp_cmd = (
                f"fastp --in1 {fq} --out1 {trimmed_fq} "
                f"-l 40 "  # Minimum length of read to keep after trimming
                f"--cut_front --cut_front_mean_quality 30 "
                f"--cut_tail --cut_tail_mean_quality 30 "
                f"-w 2 --dont_overwrite "
            )
            #if bad adapters is not defined bad_adapters
            # Add adapter sequences if any bad adapters were found
            if not bad_adapters:
                #add the timestamp
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[INFO] No bad adapters found, skipping adapter trimming at {timestamp}.")
            elif bad_adapters is None:
                #use auto adapter detection
                #add the timestamp
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[INFO] No bad adapters defined, using auto adapter detection at {timestamp}.")
                fastp_cmd += "--detect_adapter_for_pe "
            for adapter in bad_adapters:
                if adapter == "Illumina Universal Adapter":
                    fastp_cmd += '--adapter_sequence AGATCGGAAGAGC '
                elif adapter == "Illumina Small RNA 3\' Adapter":
                    fastp_cmd += "--adapter_sequence TGGAATTCTCGGGTGCCAAGG "
                elif adapter == "Illumina Small RNA 5\' Adapter":
                    fastp_cmd += "--adapter_sequence_r2 GTTCAGAGTTCTACAGTCCGACGATC "
                elif adapter == "Nextera Transposase Sequence":
                    fastp_cmd += '--adapter_sequence CTGTCTCTTATACACATCT '
                # etc.

            fastp_cmd += f"--json {trimmed_fq}.fastp.json --html {trimmed_fq}.fastp.html"
            try:
                run_command(fastp_cmd)
            except Exception as e:
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[ERROR] FastP failed at {timestamp}: {e}")
                #add 1 to the cache abort counter
                if cache is not None:
                    cache["abort_counter"][fq_basename] += 1
                return
        else:
            print("[SKIP] fastp already produced.")

        # 4. seqtk FASTA conversion
        if not fasta_file.exists():
            print(f"[INFO] seqtk FASTA not produced for {fq_basename}, running now...")
            run_command(f"seqtk seq -a {trimmed_fq} > {fasta_file}")
        else:
            print("[SKIP] seqtk FASTA already produced.")
    fasta_file=fq
    fq_basename = fasta_file.stem.replace('.fasta', '').replace('.fa', '').replace('.fna','').replace('.fasta.gz','')
    tmp_fasta = Path(tmp_dir) / fasta_file.name
    # 5. Copy to tmp dir
    if not tmp_fasta.exists():
        print(f"[INFO] Copying FASTA {fasta_file} to tmp dir: {tmp_fasta}")
        shutil.copy(fasta_file, tmp_fasta)
    else:
        print("[SKIP] FASTA already in tmpdir.")

    # 6. Run hydrogen pipeline (check with new output search)
    print(result_dir, fq_basename,"is exists?", Path(result_dir) / fq_basename, "?", output_file_exists(result_dir, fq_basename))
    if not output_file_exists(result_dir, fq_basename):
        # Make all paths absolute
        tmp_dir_abs     = str(Path(tmp_dir).resolve())
        #add a folder with the same basename as the sample, resolve it
        tmp_dir_abs   = str((Path(tmp_dir_abs) / fq_basename).resolve())
        result_dir_abs  = str(Path(result_dir).resolve())
        #add a folder with the same basename as the sample, resolve it
        out_dir_abs   = str((Path(result_dir_abs) / fq_basename).resolve())

        working_dir_abs = str(Path(working_dir).resolve())
        fasta_guess1    = Path(tmp_dir_abs) / f"{fq_basename}.fasta"
        fasta_guess2    = Path(working_dir_abs) / f"{fq_basename}.fasta"
        fasta_path_abs  = str(fasta_guess1 if fasta_guess1.exists() else fasta_guess2)
        '''
        python /oceanstor/scratch/tllseedorf/e1103389/hydrogen_consumer_pipeline/hydrogen_consumer_pipeline_reads3_4eval_kraken_fixed.py \
  --reads1 /oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/human_diabetes/reads/ERR2855939/ERR2855939_2.subset.fasta \
  --reads_for_diamond /oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/human_diabetes/reads/ERR2855939/ERR2855939_2.subset.fasta \
  --kraken_db /oceanstor/scratch/tllseedorf/e1103389/krakendbb \
  --kraken_taxonomy_dir /oceanstor/scratch/tllseedorf/e1103389/krakendbb/taxonomy \
  --hydb_fasta /oceanstor/scratch/tllseedorf/e1103389/hydrogenases/hyDB/HydDB_reformated.fasta \
  --diamond_exe /oceanstor/home/e1103389/DIAMOND/diamond \
  --mapping_script /oceanstor/scratch/tllseedorf/e1103389/hyDB/Diamond_blast_hyDB_reads_mapping.py \
  --out_dir /oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/human_diabetes/reads/ERR2855939/ERR2855939_2.subset_out \
  --tmp_dir /oceanstor/scratch/tllseedorf/e1103389/metagenomics_disease/human_diabetes/reads/ERR2855939/ERR2855939_2.subset_tmp \
  --min_qcov 0.8 --min_scov 0.0 --min_pid 80 --train_percent 100 \
  --kraken_threads 16 \
  --taxon_name_mode both

        
        '''
        #make the tmp and output dirs
        Path(tmp_dir_abs).mkdir(parents=True, exist_ok=True)
        Path(out_dir_abs).mkdir(parents=True, exist_ok=True)
        #if the input filename is in PROBLEM_LIST, then skip
        print("fasta_path_abs:", fasta_path_abs, "PROBLEM_LIST:", PROBLEM_LIST)
        if fasta_path_abs in PROBLEM_LIST:
            print(f"[SKIP] {fasta_path_abs} is in the problem list.")
            return
        # Build the new hydrogen pipeline command
        hcmd = (
            f"python {HYDROGEN_SCRIPT} "
            f"--reads1 {fasta_path_abs} "
            f"--reads_for_diamond {fasta_path_abs} "
            f"--kraken_db {KRAKEN_DB} "
            f"--hydb_fasta {HYDB_FASTA} "
            f"--diamond_exe {DIAMOND_EXE} "
            f"--mapping_script {MAPPING_SCRIPT} "
            f"--out_dir {out_dir_abs} "
            f"--tmp_dir {tmp_dir_abs} "
            f"--min_qcov {MIN_QCOV} --min_scov {MIN_SCOV} --min_pid {MIN_PID} "
            f"--train_percent {TRAIN_PERCENT} "
            f"--kraken_threads {KRAKEN_THREADS}"
        )
       
        #add the timestamp
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[INFO] Running hydrogen pipeline for {fq_basename} at {timestamp}: {hcmd}")
        #
       # print("[INFO] Running hydrogen pipeline:", hcmd)
        run_command(hcmd,offending_file=fasta_path_abs)
        if output_file_exists(out_dir_abs, fq_basename):
            #add the timestamp
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[DONE] Hydrogen output CSV for {fq_basename} created at {timestamp}.")
        else:
            print(f"[WARNING] result for {fq_basename} not found; check hydrogen pipeline logs.")
    else:
        print("[SKIP] hydrogen pipeline already produced output.")

    if do_cleanup:
        #delete the \tmp folder and /working_dir
        #rmdir the tmp dir
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        if os.path.exists(working_dir):
            shutil.rmtree(working_dir)
        print(f"[CLEANUP] Removed temporary directories: {tmp_dir} and {working_dir}")
        #make the empty tmp_dir again
        Path(tmp_dir).mkdir(exist_ok=True)
        # Optionally, remove the trimmed FASTQ
        #make the empty working_dir again
        Path(working_dir).mkdir(exist_ok=True)
def output_file_exists(result_dir, fq_basename):
    """Return True if the hydrogen pipeline outputs exist in result_dir."""
    rd = Path(result_dir) / Path(fq_basename)
    core = ["taxon3_hydgroup_counts.csv", "taxon3_hydgroup_rpkm.csv","taxon3_marker_counts.csv", "taxon3_marker_rpkm.csv"]
    if all((rd / x).exists() for x in core):
        return True
    return False

# --- Main balanced driver ---
def balanced_main():
    male_files   = list(MALE_DIR.glob("*.fasta"))
    female_files = list(FEMALE_DIR.glob("*.fastq.gz"))

    random.shuffle(male_files)
    random.shuffle(female_files)

    n_male_done = 0
    n_female_done = 0

    while n_male_done < TARGET_MALE or n_female_done < TARGET_FEMALE:
        if n_male_done < TARGET_MALE and male_files:
            fq = male_files.pop()
            sample_id = fq.stem
            result_dir = BASE_RESULT_DIR / sample_id
            tmp_dir    = BASE_TMP_DIR / sample_id
            work_dir   = BASE_WORKING_DIR / sample_id
            result_dir.mkdir(exist_ok=True, parents=True)
            tmp_dir.mkdir(exist_ok=True, parents=True)
            work_dir.mkdir(exist_ok=True, parents=True)
            try:
                if not output_file_exists(BASE_RESULT_DIR, sample_id):
                    shutil.copy(fq, work_dir / fq.name)
                    #add the timestamp
                    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    print(f"[INFO] Copied {fq} to working dir: {work_dir / fq.name} at {timestamp}")
                    process_one_fastq(str(work_dir / fq.name), str(work_dir), str(result_dir), str(tmp_dir), do_cleanup=True, cache=None)
                n_male_done += 1
                #add the timestamp
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[MALE DONE] {n_male_done}/{TARGET_MALE} at {timestamp}")
            except Exception as e:
                #add the timestamp
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[ERROR] Male {fq}: {e} at {timestamp}")

        if n_female_done < TARGET_FEMALE and female_files:
            fq = female_files.pop()
            sample_id = fq.stem.replace(".fastq", "").replace(".fq","").replace(".gz","")
            result_dir = BASE_RESULT_DIR / sample_id
            tmp_dir    = BASE_TMP_DIR / sample_id
            work_dir   = BASE_WORKING_DIR / sample_id
            result_dir.mkdir(exist_ok=True, parents=True)
            tmp_dir.mkdir(exist_ok=True, parents=True)
            work_dir.mkdir(exist_ok=True, parents=True)
            try:
                if not output_file_exists(BASE_RESULT_DIR, sample_id):
                    process_one_fastq(str(fq), str(work_dir), str(result_dir), str(tmp_dir), do_cleanup=True, cache=None)
                n_female_done += 1
                #add the timestamp
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[FEMALE DONE] {n_female_done}/{TARGET_FEMALE} at {timestamp}")
            except Exception as e:
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[ERROR] Female {fq}: {e} at {timestamp}")

        if n_male_done >= TARGET_MALE and n_female_done >= TARGET_FEMALE:
            break

    print(f"✅ Finished: {n_male_done} males, {n_female_done} females")

if __name__ == "__main__":
    balanced_main()
