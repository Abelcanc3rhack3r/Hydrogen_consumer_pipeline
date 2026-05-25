import argparse
import random
import re
import subprocess
import os
from pathlib import Path
from Bio import SeqIO
# Argument parsing setup
def main(args):
    
    print("running")

    # Set up your directories
    DiaDB_DIR = "/tllhome/abel/hydrogen_consumer_pipeline/hyDB"
    #if DiaDB_DIR does not exist, use args.DIA_DB_DIR
    if not os.path.exists(DiaDB_DIR):
        DiaDB_DIR = '/oceanstor/home/e1103389/hyDB'#HydDB_reformated.fasta
    #if no exist, use C:\Users\abel\Documents\hydrogenases\hyDB
    if not os.path.exists(DiaDB_DIR):
        DiaDB_DIR = r'C:\Users\abel\Documents\hydrogenases\hyDB'
    #r"C:\Users\abel\Documents\hydrogenases\hyDB"
    hydrogenase_file=r"C:\Users\abel\Documents\hydrogenases\hyDB\HydDB_reformated.fasta"
    #if hydrogenase file not exist, use the alternate : sftp://abel@hpc.tll.org.sg/tllhome/abel/hydrogen_consumer_pipeline/hyDB/HydDB_reformated.fasta
    if not os.path.exists(hydrogenase_file):
        hydrogenase_file = r'/tllhome/abel/hydrogen_consumer_pipeline/hyDB/HydDB_reformated.fasta'
    #if still not exist, use /oceanstor/home/e1103389/hyDB/HydDB_reformated.fasta
    if not os.path.exists(hydrogenase_file):
        hydrogenase_file = r'/oceanstor/home/e1103389/hyDB/HydDB_reformated.fasta'
    #get the read lengths
    subj_lens={}
    for record in SeqIO.parse(hydrogenase_file, "fasta"):
        subj_lens[record.id]=len(record.seq)
    
    # No need to change into directories, specify paths in commands instead
    db_path = os.path.join(DiaDB_DIR, 'HydDB.dmnd')
    #if train is not 100, then make a new diamond db
    if args.train != 100:
            #open the hydrogenase file and randomly select x% of the sequences
            #write the sequences    
            with open(hydrogenase_file) as file:
                records = list(SeqIO.parse(file, "fasta"))
                #calculate the number of sequences to select
                num_records = int(len(records) * (float(args.train)/100))
                #randomly select the sequences
                selected_records = random.sample(records, num_records)
                #write the sequences to a new file
                with open(os.path.join(DiaDB_DIR, 'HydDB_train.fasta'), 'w') as outfile:
                    SeqIO.write(selected_records, outfile, "fasta")
            #create a new diamond database
            subprocess.run([
                args.Diamond_EXE, 'makedb',
                '--in', os.path.join(DiaDB_DIR, 'HydDB_train.fasta'),
                '--db', os.path.join(DiaDB_DIR, 'HydDB_train.dmnd')
            ])
            print("made a new db with", num_records, "sequences out of", len(records))
            db_path = os.path.join(DiaDB_DIR, 'HydDB_train.dmnd')
    # Create necessary directories
    for sub_dir in ['Raw', 'Prefiltered', 'Filtered', 'Summary']:
        os.makedirs(os.path.join(args.BLAST_DIR, sub_dir), exist_ok=True)
    diamond_exe = r"D:\diamond-windows\diamond.exe"
    #if diamond no exist, use sftp://e1103389@kylin.cbis.nus.edu.sg/oceanstor/home/e1103389/DIAMOND/diamond
    if not os.path.exists(diamond_exe):
        diamond_exe = '/oceanstor/home/e1103389/DIAMOND/diamond'
    # Perform diamond blast

    def detect_sequences(input_folder):
        # Define regular expressions for DNA and amino acid sequences
        dna_regex = re.compile(r'^[ACGTacgt]+$')
        amino_acid_regex = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+$')
        
        # List all files in the input folder
        files = os.listdir(input_folder)
        
        # Ensure there is at least one file
        if not files:
            return "No files in the input folder"
        
        # Choose the first file
        first_file = files[0]
        filepath = os.path.join(input_folder, first_file)
        
        # Ensure the path is a file
        if os.path.isfile(filepath):
            sequence_data = []
            sequence_started = False
            
            with open(filepath, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        if sequence_started:
                            break
                        sequence_started = True
                    elif sequence_started:
                        sequence_data.append(line)
                        if sum(len(seq) for seq in sequence_data) >= 100:
                            break
            
            # Concatenate and truncate the sequence data to the first 100 characters
            content = ''.join(sequence_data)[:100]
            
            # Check if the content matches DNA or amino acid sequences
            if amino_acid_regex.fullmatch(content):
                return 'aa'
            elif dna_regex.fullmatch(content):
                return 'dna'
            else:
                return 'Unknown'
        else:
            return "First item is not a file"


    blasttype='blastx'

    for bin_file in  os.listdir(args.Bin_DIR):
        #if file ends with fa, fasta or fna
        if not bin_file.endswith(('.fa', '.fasta', '.fna')):
            continue
        query_lengths={}
        #open the file
        with open(os.path.join(args.Bin_DIR, bin_file)) as file:
            #parse the file
            for record in SeqIO.parse(file, "fasta"):
                query_lengths[record.id]=len(record.seq)

        filepath = os.path.join(args.Bin_DIR, bin_file)
        base = os.path.basename(bin_file).split('.')[0]
        output_path = os.path.join(args.BLAST_DIR, 'Raw', f"{base}_hydrogenase.txt")
        subprocess.run([
            diamond_exe, blasttype,
            '--db', db_path,
            '--query', str(filepath),
            '--out', output_path,
            '--max-hsps', '1',
            '--max-target-seqs', '1',
            '--threads', '10',
            #'--outfmt','tab'
            '--outfmt', '6'
        ])
    #'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    min_scov=args.min_scov
    min_qcov=args.min_qcov
    # Prefilter based on query and subject cover
    for hit_path in Path(args.BLAST_DIR, 'Raw').glob("*_hydrogenase.txt"):
        prefiltered_path = os.path.join(args.BLAST_DIR, 'Prefiltered', f"{hit_path.stem}_prefiltered.txt")
        with open(hit_path) as infile, open(prefiltered_path, 'w') as outfile:
            for line in infile:
                fields = line.split('\t')
                query_len=query_lengths[fields[0]]
                subj_len=subj_lens[fields[1]]
                qextent= abs(int(fields[6])-int(fields[7]))
                sextent= abs(int(fields[8])-int(fields[9]))
                qcov= qextent/query_len
                scov= sextent/subj_len
                #if qcov >= 0.8 or scov >= 0.8:
                if scov >= float(min_scov) and qcov >= float(min_qcov):
                    outfile.write(line)
                #if float(fields[10]) >= 80 or float(fields[11]) >= 80:
                    #outfile.write(line)

    # Filter raw hits
    min_pid=args.min_pid
    
    for prefiltered_path in Path(args.BLAST_DIR, 'Prefiltered').glob("*_prefiltered.txt"):
        filtered_path = os.path.join(args.BLAST_DIR, 'Filtered', f"{prefiltered_path.stem.replace('_prefiltered', '')}_filtered.txt")
        with open(prefiltered_path) as infile, open(filtered_path, 'w') as outfile:
            for line in infile:
                fields = line.split('\t')
                percentage_identity = float(fields[2])
                title = fields[1]
                '''if "FeFe-" in title and percentage_identity >= 60:
                    outfile.write(line)
                elif "[Fe]" in title and percentage_identity >= 50:
                    outfile.write(line)
                elif "NiFe-" in title:
                    if "Group 4" not in title or percentage_identity >= 60:
                        outfile.write(line)'''
                if percentage_identity >= float(min_pid):
                    outfile.write(line)

    # Concatenate all summary files for a dataset
    with open(os.path.join(args.BLAST_DIR, 'Filtered', 'all_hydrogenase.txt'), 'w') as outfile:
        for filtered_file in Path(args.BLAST_DIR, 'Filtered').glob("*_filtered.txt"):
            with open(filtered_file) as infile:
                for line in infile:
                    outfile.write(f"{line.strip()}\t{filtered_file.name}\n")


def tes():
    input_file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\hyinput"
    blast_dir=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\blast"
    diamond_exe=r"C:\Users\abel\Documents\hydrogenases\hyDB"
    args=argparse.Namespace(Bin_DIR=input_file, BLAST_DIR=blast_dir, Diamond_EXE=diamond_exe)
    main(args)
rc_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
def rc(string):
    return ''.join([rc_dict.get(base, "N") for base in string[::-1]])

def get_files(root_dir):
    return {"hydrogen_producer":os.path.join(root_dir,"hydrogen_producer","fasta"),
            "hydrogen_consumer":os.path.join(root_dir,"hydrogen_consumer","fasta"),
            "non_hydrogen_producer":os.path.join(root_dir,"non_hydrogen_producer","fasta"),
            "non_hydrogen_consumer":os.path.join(root_dir,"non_hydrogen_consumer","fasta")}
from Bio import SeqIO
def make_reads(files_list, read_len=100, frag_length=(1500,3000),  num_reads=1000):
    reads_list={}
    for typed, path in files_list.items():
        
        files_dict=[]
        with open(os.path.join(path, f"{typed}.fasta"), 'w') as outfile:
            for fasta_file in os.listdir(path):
                #skip if not fasta
                if not fasta_file.endswith(('.fa', '.fasta', '.fna')):
                    continue
                files_dict.append(fasta_file)
            #allocate the number of reads to each file, genrate a multinomial distribution
            probabilities=[1]*len(files_dict)
            multi=random.multinomial(num_reads, probabilities)

            x=0
            for fasta_file in files_dict:
                file_id=fasta_file.split('.')[0]
                print("writing reads for", fasta_file)
                with open(os.path.join(path, fasta_file)) as file:
                    records = list(SeqIO.parse(file, "fasta"))
                    
                    
                    for i in range(multi[files_dict.index(fasta_file)]):

                        #choose a random record
                        rec = random.choice(records)
                        start = random.randint(0, len(rec.seq)-frag_length[1])
                        end = start + random.randint(frag_length[0], frag_length[1])
                        frag = rec.seq[start:end]
                        #get the first read_len characters
                        read = frag[:read_len]
                        #get the last read_len characters
                        read2 = frag[-read_len:]
                        #reverse complement the second read
                        read2 = rc(read2)
                        #write the read
                        id_fwd=f"{file_id}_{rec.id}_fwd_{start}_{end}"
                        id_rev=f"{file_id}_{rec.id}_rev_{start}_{end}"
                        outfile.write(f">{id_fwd}\n{read}\n")
                        outfile.write(f">{id_rev}\n{read2}\n")
                        x+=1
                        if x%10==0:
                            print("writing read", id_fwd)
                            print("writing read", id_rev)
        reads_list[typed]=path
    return reads_list

                        
def run_diamond(hyd_groups_folder,reads_list,root_dir,diamond_exe,min_pids=[80],min_qcovs=[0.8],train=100,paired=True):
    hyd_groups=[f for f in os.listdir(hyd_groups_folder) if f.endswith('.fasta')]
    output_files={}
    for grp in hyd_groups:
        for read_type, path in reads_list.items():
            print("running diamond for", grp, read_type)
            grp_name=grp.split('.')[0]
            output_dir=os.path.join(root_dir, "blast", grp_name, read_type)
            os.makedirs(output_dir, exist_ok=True)

            db_path = os.path.join(hyd_groups_folder, grp)
            #make a diamond db with a random sample of the hydrogenase group
            if train != 100:
                with open(db_path) as file:
                    records = list(SeqIO.parse(file, "fasta"))
                    num_records = int(len(records) * (float(train)/100))
                    selected_records = random.sample(records, num_records)
                    with open(os.path.join(hyd_groups_folder, f"{grp_name}_train.fasta"), 'w') as outfile:
                        SeqIO.write(selected_records, outfile, "fasta")
                subprocess.run([
                    diamond_exe, 'makedb',
                    '--in', os.path.join(hyd_groups_folder, f"{grp_name}_train.fasta"),
                    '--db', os.path.join(hyd_groups_folder, f"{grp_name}_train.dmnd")
                ])
                print("made a new db with", num_records, "sequences out of", len(records))
                db_path = os.path.join(hyd_groups_folder, f"{grp_name}_train.dmnd")

            output_file=os.path.join(output_dir, f"{grp_name}_{read_type}_{train}_hydrogenase.txt")
            print("running diamond for", grp, read_type)
            #  os.system(f'{diamond_exe} {blasttype} --db {db_path} --query "{str(filepath)}" --out "{output_path}" --max-hsps 1 --max-target-seqs 1 --threads 10 --outfmt 6 qtitle stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp  --quiet')
            command=f'{diamond_exe} blastx --db {db_path} --query "{path}" --out "{output_file}" --max-hsps 1 --max-target-seqs 1 --threads 10 --outfmt 6 qtitle stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp  --quiet'
            print(command)
            os.system(command)
            read_pairs=defaultdict(lambda : [False, False])
            hits={}
            with open (output_file , 'r') as infile:
                for min_pid in min_pids:
                    for min_qcov in min_qcovs:
                        for line in infile:
                            fields = line.split('\t')
                            
                            qcov= fields[-2]
                            #scov= fields[-1]
                            pident=fields[2]
                            if qcov >= min_qcov and pident >= min_pid:
                                read_pairs_id=fields[0].replace("_fwd", "").replace("_rev", "")
                                #if fwd in the id, set the first element to True
                                if "_fwd" in fields[0]:
                                    read_pairs[read_pairs_id][0]=True
                                    print("read_pairs_id hit", read_pairs_id)
                                #if rev in the id, set the second element to True
                                if "_rev" in fields[0]:
                                    read_pairs[read_pairs_id][1]=True
                                    print("read_pairs_id hit", read_pairs_id)
                                    hits[read_pairs_id]=line
                    #if paired is True, write the read only if both the forward and reverse reads have hits
                    for read_id, hit in hits.items():
                        if read_pairs[read_id][0] and read_pairs[read_id][1] or not paired:
                            with open(os.path.join(output_dir, f"{train}_{min_pid}_{min_qcov}_{read_type}_hits.txt"), 'w') as outfile:
                                outfile.write(f"{hit}\n")
                    output_files[(grp_name, read_type, train, min_pid, min_qcov)]=output_file
    return output_files
import csv
def evaluate(output_files, num_reads, hyd_direction_file, min_pids=[80], min_qcovs=[0.8], trains=[100]):
    hyd_direction={"consuming":[], "producing":[]}
    with open(hyd_direction_file) as infile:
        reader=csv.reader(infile, delimiter='\t')
        for line in reader:
            cls=line[0]
            grp=line[1]
            #if grp is 3c, add to consuming and continue
            if grp=="3c":
                hyd_direction["consuming"].append(grp)
                continue
            #replace Hmd with ""
            grp=grp.replace("Hmd", "")
            #replace NiFeSe class with NiFe
            grp=grp.replace("NiFeSe", "NiFe")
            #if the 3rd col contains Bidirectional, Bifurcation or uptake, set the direction to consuming
            dire=line[2]
            if "bidirectional" in dire.lower() or "bifurcation" in dire.lower() or "uptake" in dire.lower():
                hyd_direction["consuming"].append(grp)
            #if the 3rd col contains evolution or bidirectional, or bifurcation, set the direction to producing
            if "evolution" in dire.lower() or "bidirectional" in dire.lower() or "bifurcation" in dire.lower():
                hyd_direction["producing"].append(grp)

        for grp in hyd_direction["consuming"]:
            print("evaluating", grp, read_type)
            with open(output_files[grp]) as infile:
                for line in infile:
                    fields = line.split('\t')
                    read_id=fields[0]
                    if read_id not in num_reads[read_type]:
                        num_reads[read_type][read_id]=0
                    num_reads[read_type][read_id]+=1

                    
        

            
            





def evaluate():
    pass






if __name__ == '__main__':
    #tes()
    parser = argparse.ArgumentParser(description='Use Diamond to search against hydrogenase database.')
    parser.add_argument('Bin_DIR', help='Enter the genome file pathways')
    parser.add_argument('BLAST_DIR', help='Enter the output directory')
    #add  the path to the diamond executable
    parser.add_argument('Diamond_EXE', help='Enter the path to the diamond executable', default='/oceanstor/home/e1103389/DIAMOND/diamond')
    #add an argument for minimum subject coverage
    parser.add_argument('--min_scov', help='Enter the minimum subject coverage', default=0)
    #add an argument for minimum query coverage
    parser.add_argument('--min_qcov', help='Enter the minimum query coverage', default=0.8)
    #add an argument for minimum percent identity
    parser.add_argument('--min_pid', help='Enter the minimum percent identity', default=80)
    #add an argument for train test, reduces the database to x% size
    parser.add_argument('--train', help='Enter the percentage of the database to use', default=100)
    #print the args
    
    args = parser.parse_args()
    print(f"diamond_blast hydb reads.py", f"--min_scov {args.min_scov}", f"--min_qcov {args.min_qcov}", f"--min_pid {args.min_pid}", f"--train {args.train}")
    main(args)
