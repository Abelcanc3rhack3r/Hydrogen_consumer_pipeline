import argparse
import subprocess
import os
from pathlib import Path

# Argument parsing setup
parser = argparse.ArgumentParser(description='Use Diamond to search against hydrogenase database.')
parser.add_argument('Bin_DIR', help='Enter the genome file pathways')
parser.add_argument('BLAST_DIR', help='Enter the output directory')

#add  the path to the diamond executable
parser.add_argument('Diamond_EXE', help='Enter the path to the diamond executable')
#add the path to the diamond db dir
parser.add_argument('DiaDB_DIR', help='Enter the path to the diamond database')
args = parser.parse_args()

# Set up your directories
DiaDB_DIR = args.DiaDB_DIR
#DiaDB_DIR = 'C:\\Users\\abel\\Documents\\hydrogenases\\hyDB'

# No need to change into directories, specify paths in commands instead
db_path = os.path.join(DiaDB_DIR, 'HydDB.dmnd')

# Create necessary directories
for sub_dir in ['Raw', 'Prefiltered', 'Filtered', 'Summary']:
    os.makedirs(os.path.join(args.BLAST_DIR, sub_dir), exist_ok=True)
diamond_exe = args.Diamond_EXE
# Perform diamond blast
for bin_file in Path(args.Bin_DIR).glob('*.fasta'):
    base = bin_file.stem
    output_path = os.path.join(args.BLAST_DIR, 'Raw', f"{base}_hydrogenase.txt")
    subprocess.run([
        diamond_exe, 'blastp',
        '--db', db_path,
        '--query', str(bin_file),
        '--out', output_path,
        '--max-hsps', '1',
        '--max-target-seqs', '1',
        '--threads', '10',
        '--outfmt','tab'
        #'--outfmt', '6 qtitle stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp full_qseq'
    ])

# Prefilter based on query and subject cover
for hit_path in Path(args.BLAST_DIR, 'Raw').glob("*_hydrogenase.txt"):
    prefiltered_path = os.path.join(args.BLAST_DIR, 'Prefiltered', f"{hit_path.stem}_prefiltered.txt")
    with open(hit_path) as infile, open(prefiltered_path, 'w') as outfile:
        for line in infile:
            fields = line.split('\t')
            if float(fields[10]) >= 80 or float(fields[11]) >= 80:
                outfile.write(line)

# Filter raw hits
for prefiltered_path in Path(args.BLAST_DIR, 'Prefiltered').glob("*_prefiltered.txt"):
    filtered_path = os.path.join(args.BLAST_DIR, 'Filtered', f"{prefiltered_path.stem.replace('_prefiltered', '')}_filtered.txt")
    with open(prefiltered_path) as infile, open(filtered_path, 'w') as outfile:
        for line in infile:
            fields = line.split('\t')
            percentage_identity = float(fields[2])
            title = fields[1]
            if "FeFe-" in title and percentage_identity >= 60:
                outfile.write(line)
            elif "[Fe]" in title and percentage_identity >= 50:
                outfile.write(line)
            elif "NiFe-" in title:
                if "Group 4" not in title or percentage_identity >= 60:
                    outfile.write(line)

# Concatenate all summary files for a dataset
with open(os.path.join(args.BLAST_DIR, 'Filtered', 'all_hydrogenase.txt'), 'w') as outfile:
    for filtered_file in Path(args.BLAST_DIR, 'Filtered').glob("*_filtered.txt"):
        with open(filtered_file) as infile:
            for line in infile:
                outfile.write(f"{line.strip()}\t{filtered_file.name}\n")