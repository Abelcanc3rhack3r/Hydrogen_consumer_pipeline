import os
import csv
import re
import numpy as np
from Bio import SeqIO
import pandas as pd
debug=True
parent_dir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
default_diamond_exe_path='/tllhome/abel/DIAMOND/diamond'
#refactor to use parent dir as tllhome/abel
default_diamond_exe_path=os.path.join(parent_dir,'DIAMOND','diamond')
default_DIADB_dir=os.path.join(parent_dir,'hyDB')
#default_prodigal_path=r'/tllhome/abel/hydrogen_cons
#r"D:\diamond-windows\diamond.exe"
#default_hydrogenase_script_path=r"C:\Users\abel\Documents\hydrogenases\hyDB\Diamond_blastp_hyDB.sh"
#now the hydrogenase script is a python script
default_hydrogenase_script_path=r'/tllhome/abel/hydrogen_consumer_pipeline/hyDB/Diamond_blast_hyDB.py'
#refactor to use parent dir as /tllhome/abel/hydrogen_consumer_pipeline
def create_diamond_database(fasta_file, tmp_folder,diamond_path,train_percent=100):
    #format the fasta file so that each record has a len
    #open the fasta file
    with open(fasta_file, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        #create a new file with the formatted records
        formatted_records = []
        for record in records:
            formatted_records.append(record)
    #write the formatted records to a new file
    formatted_fasta_file = os.path.join(tmp_folder, os.path.basename(fasta_file))
    #if train is not 100, sample the records
    if train_percent!=100:
        formatted_records=random.sample(formatted_records,int(train_percent*len(formatted_records)/100))
    SeqIO.write(formatted_records, formatted_fasta_file, 'fasta')

    #create a diamond database from the fasta file
    #get the name of the fasta file
    name=os.path.basename(formatted_fasta_file)
    #get the name without the extension
    name=os.path.splitext(name)[0]+".dmnd"
    #add the name to the tmp folder
    name=os.path.join(tmp_folder,name)
    #create the database
    command=diamond_path+' makedb --in '+formatted_fasta_file+' -d '+name
    os.system(command)
    #return the path to the database
    return name
def batch_diamond(joined_records, diamond_exe_path, database_fasta, type, tmp_dir, min_percent_identity=90, min_qcov=0.9, resume=False,paired=False,train_percent=100):
    database_path=database_fasta.replace('.fasta','.dmnd')
    #if database path not exists, create a diamond database from the fasta file
    if not os.path.exists(database_path) or train_percent!=100:
        database_path=create_diamond_database(database_fasta,tmp_dir,diamond_exe_path,train_percent)
    #if database path is a fasta file, create a diamond database from it
    #if database_path.endswith('.fasta'):
        #create a diamond database from the fasta file
    #    database_path=create_diamond_database(database_path,tmp_dir,diamond_exe_path)
    #join the records in the input dir into one fasta file, add an identifier to the start of each record representing the file it came from
    #joined_records=join_and_add_labels_to_records(input_dir,os.path.join(tmp_dir,type+"_joined_records.fasta"))
    #run diamond on the joined records
    output_file=os.path.join(tmp_dir,type+"_diamond_output_filtered.tsv")
    
    #skip the diamond if output file exists and resume is true
    if not (resume and os.path.exists(output_file)):
        command=diamond_exe_path+' blastx -d '+database_path+' -q ' +joined_records+' -o '+os.path.join(tmp_dir,type+"_diamond_output.tsv")+' --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen --fast --quiet'
        print("COMMAND:",command)
        os.system(command)
    # open the database fasta file
    if database_fasta is None:
        raise Exception("Database fasta file is None")
    # get the lengths of the sequences in the database
    lengths = {}
    with open(database_fasta, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            lengths[record.id] = len(record.seq)*3
    # filter the rows by min percent identity and min scov
    with open(os.path.join(tmp_dir,type+"_diamond_output.tsv"), 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = [row for row in reader]
    # filter the rows by min percent identity and min scov
    filtered_rows = []
    for row in rows:
        pid = False
        qcov = False
        if float(row[2]) >= min_percent_identity:
            pid = True
        qstart=int(row[6])
        qend=int(row[7])
        #sstart = int(row[8])
        #send = int(row[9])
        #align = (send - sstart + 1)
        #scov1 = align / lengths[row[1]]
        qlen=int(row[12])
        qcov1=(qend-qstart+1)/qlen

        if qcov1 >= min_qcov:
            qcov = True
        if pid and qcov:
            filtered_rows.append(row)
   #write the filtered rows to a new file
    with open(os.path.join(tmp_dir,type+"_diamond_output_filtered.tsv"), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(filtered_rows)
    hist,total=parse_reads(os.path.join(tmp_dir,type+"_diamond_output_filtered.tsv"),joined_records, lengths,paired=paired)
    return hist, total
from collections import defaultdict
def extract_hyd_group(subj):
    grp=subj.split("_-_[")[1]
    grp=grp.replace('[','').replace(']','')
    return grp
def group_paired_reads(diamond_output_file,reads_file, lengths,hyddb=False):
    #read the diamond output file
    hits=defaultdict(lambda : defaultdict(int))
    with open(diamond_output_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            #if len of query is 0, continue
            if len(row)==0:
                continue
            sbj=row[1]
            if hyddb:
                subject=extract_hyd_group(row[1])
                
            else:
                subject="match"

            query=row[0]
            evalue=row[10]
            hits[query][subject]=(evalue,sbj)
    #find pairs of reads in the hits
    new_rows=[]
    #rids=set()
    #the subject of a read pair is for now, the hit of the FWD read
    for query in hits:
            #other query , find 1:N or 1:Y, and replace it with  2:N or 2:Y
            regex=re.compile(r'1:[NY]')
            #check if the query contains this regex
            if not regex.search(query):
                continue
            oquery=regex.sub('2:N',query)
            #find the other query in the hits
            if oquery in hits:
                #find the common subjects
                common_subjects=set(hits[query].keys()).intersection(set(hits[oquery].keys()))
                for subject in common_subjects:
                    #add the forward and reverse read id to the rids
                    #rids.add((query, subject))
                    #rids.add((oquery, subject))
                    #combine the evalues
                    #E_combined = 1 - (1 - E1) * (1 - E2)
                    evalue1=hits[query][subject][0]
                    evalue2=hits[oquery][subject][0]
                    log_evalue1=np.log10(float(evalue1))
                    log_evalue2=np.log10(float(evalue2))
                    evalue_combined=(log_evalue1+log_evalue2)
                    sbj=hits[query][subject][1]
                    #evalue_combined=1-(1-evalue1)*(1-evalue2)
                    #make a row with 12 columns, 10 is evalue, 11 is bitscore
                    new_row=[query,sbj]+[""]*8+[evalue_combined]+[""]
                    new_rows.append(new_row)
    #return the new rows
    return new_rows


consuming_hydrogenases= [
    "FeFe_group_A2",
    "FeFe_group_A3",
    "FeFe_group_A4",
    "FeFe",
    "NiFe_group_1a",
    "NiFeSe_group_1a",
    "NiFe_group_1b",
    "NiFe_group_1c",
    "NiFe_group_1d",
    "NiFe_group_1e",
    "NiFe_group_1f",
    "NiFe_group_1g",
    "NiFe_group_1h",
    "NiFe_group_1i",
    "NiFe_group_1j",
    "NiFe_group_1k",
    "NiFe_group_2a",
    "NiFe_group_2d",
    "NiFe_group_2e",
    "NiFe_group_3a",
    "NiFe_group_3b",
    "NiFe_group_3c",
    "NiFe_group_3d",
    "NiFe_group_4h",
    "NiFe_group_4i"
]
producing_hydrogenases = [
    "FeFe_group_A1",
    "FeFe_group_A3",
    "FeFe_group_A4",
    "FeFe_group_B",
    "FeFe",
    "NiFe_group_3a",
    "NiFe_group_3b",
    "NiFe_group_3c",
    "NiFe_group_3d",
    "NiFe_group_4a",
    "NiFe_group_4b",
    "NiFe_group_4c",
    "NiFe_group_4d",
    "NiFe_group_4e",
    "NiFe_group_4f",
    "NiFe_group_4g"
]

def parse_reads(diamond_output_file,reads_file, lengths,hyddb=False,paired=False):
    #open the reads file and count the number of reads
    #if read file is fastq
    ext=os.path.splitext(reads_file)[1]
    print("EXTENSION:",ext)
    if ext in [".fq",".fastq"]:
        num_lines=0
        with open(reads_file, 'r') as f:
            for line in f:
                num_lines+=1
        num_reads= num_lines/4
    else:
        #count the number of >
        num_reads=0
        with open(reads_file, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    num_reads+=1
    print("There are",num_reads,"reads in the reads file",reads_file)
    #open the diamond output file
    with open(diamond_output_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = [row for row in reader if len(row)>0]
    if paired:
        new_rows=group_paired_reads(diamond_output_file,reads_file, lengths,hyddb)
        nr_0=[new_row[0] for new_row in new_rows]
        oquery_0=[new_row[0].replace('1:N','2:N') for new_row in new_rows]
        #remove the rows with the same query as new_rows
        for row in rows:

            if row[0] in nr_0 or row[0] in oquery_0:
                #print("get rid of",row  )
                continue
            new_rows.append(row)

        rows=new_rows
    #for each read, find the rows
    read_hits=defaultdict(list)
    for row in rows:
            if len(row)==0:
                continue
            #if row[0] contains "prod" and row[1] contains one of the producing hydrogenases, then skip
            if "prod" in row[0] and any([prod in row[1] for prod in producing_hydrogenases]):
                print("skipping",row,"coz it is a producing hydrogenase")
                continue
            #if row[0] contains "cons" and row[1] contains one of the consuming hydrogenases, then skip
            if "cons" in row[0] and any([cons in row[1] for cons in consuming_hydrogenases]):
                print("skipping",row,"coz it is a consuming hydrogenase")
                continue
            read_hits[row[0]].append(row)
    #the e value is row K, get the lowest e value for each read
    best_hits={}
    for read in read_hits:
        best_hits[read]=read_hits[read][0]
        log_evalue_best=np.log10(float(best_hits[read][10]))
        for hit in read_hits[read]:
            log_evalue=np.log10(float(hit[10]))
            #if float(hit[10])>float(best_hits[read][10]):
            if log_evalue<log_evalue_best:
                best_hits[read]=hit
    #make a mapping file of read to subject
    if debug:
        with open(r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\read_to_subject_mapping.txt", 'w') as f:
            f.write("read,subject\n")
            for read in best_hits:
                f.write(read+","+best_hits[read][1]+"\n")
    #make a histogram of subjects and their reads
    subject_reads={}
    for read in best_hits:
        subject=best_hits[read][1]
        if subject not in subject_reads:
            subject_reads[subject]=[]
        subject_reads[subject].append(read)
    #extract the lengths of the subjects
    subject_lengths=lengths
    #print the lengths to a file
    if debug:
        with open(R"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\subject_lengths.txt", 'w') as f:
            #print the total number of reads
            f.write("Total number of reads,"+str(num_reads)+"\n")
            for subj in subject_lengths:
                f.write(subj+","+str(subject_lengths[subj])+"\n")
            
    #calculate the rkpm
    #numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
    total_reads=num_reads
    rkpm={}
    num_reads1={}
    if debug:
            with open (r'C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\rkpm_formula.txt','w') as f:
                f.write(f"rkpm formula\n")
    #relative_abundances={}
    for subj in subject_reads:
        num_reads=len(subject_reads[subj])
        num_reads1[subj]=num_reads
        gene_length=subject_lengths[subj]
        #relative_abundances[subj]=num_reads/total_reads
        rkpm[subj]=1e9*num_reads/(gene_length*total_reads)
        #write the formula to a file
        if debug:
            with open (r'C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\rkpm_formula.txt','a') as f:
                f.write(f"rkpm{subj}=1e9*{num_reads}/{gene_length}*{total_reads} = {rkpm[subj]}\n")
    #print the rkpm
    for subj in rkpm:
        print(subj,"has an rkpm of",rkpm[subj])
    #sum the rkpm
    total_rkpm=0
    num_reads=0
    for subj in rkpm:
        total_rkpm+=rkpm[subj]
        num_reads+=num_reads1[subj]
    if hyddb:
        grp_sum=sum_hyd_groups(rkpm)
        num_reads_sum=sum_hyd_groups(num_reads1)
        return num_reads_sum,grp_sum
    return num_reads,total_rkpm
def sum_hyd_groups(rkpm):
    grp_sum=defaultdict(float)
    for key in rkpm:
        group=key.split("_-_[")[1]
        #remove the square brackets
        group=group.replace('[','').replace(']','')
        grp_sum[group]+=rkpm[key]
    return grp_sum



import shutil
def batch_run_hydrogenase_script(ijoined_records,hydrogenase_script_path,tmp_dir,resume=False,
                                 diamond_exe_path=default_diamond_exe_path,hydb_fasta_file=None,paired=False,train_percent=100, min_pident=80, min_qcov=0.8):
    print("Running hydrogenase script")
    #joined_records = join_and_add_labels_to_records(input_dir, os.path.join(tmp_dir, "hyd_joined_records.fasta"))
    '''command:
    bash {hydrogenase_script_path} {hyinput} {hyoutput}
    '''
    #copy the joined_records to tmp_dir/hyDinput
    #make the hyinput dir, if hyinput dir exists, and resume is false, then remove the hyinput dir

    output_filepath = os.path.join(tmp_dir, "hyDB_output")
    lengths = {}
    #open the database fasta file aan get the lengths of the sequences
    print("HYDB FILE:",hydb_fasta_file)
    with open(hydb_fasta_file, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            lengths[record.id] = len(record.seq)*3

    # if in resume mode, and the output filepath exists, then return the output file
    
    #print("Output file exists, skipping hydrogenase script")
        

    '''if os.path.exists(os.path.join(tmp_dir,'hyinput')) and not resume:
        shutil.rmtree(os.path.join(tmp_dir,'hyinput'))
    os.makedirs(os.path.join(tmp_dir,'hyinput'),exist_ok=True)
    hyinput_dir=os.path.join(tmp_dir,'hyinput')
    shutil.copy(ijoined_records,os.path.join(tmp_dir,'hyinput'))
    #rename the file to .fasta
    #remove the extension
    new_name=os.path.splitext(os.path.basename(ijoined_records))[0]
    new_name= new_name+'.fasta'
    #if new name file already exists, remove it
    if os.path.exists(os.path.join(hyinput_dir,new_name)):
        os.remove(os.path.join(hyinput_dir,new_name))
    os.rename(os.path.join(hyinput_dir,os.path.basename(ijoined_records)),os.path.join(hyinput_dir,new_name))'''
    #cp the file to a tmp dir
    #hyinput_dir=os.path.join(tmp_dir,"hyinput")
    #shutil.copy(ijoined_records,os.path.join(hyinput_dir,os.path.basename(ijoined_records)))
    hyinput_dir=os.path.dirname(ijoined_records)
    #hyinput_dir=os.path.join(tmp_dir,"hyinput")
    #mkdir the output dir if it does not exist, else if not resume, remove the dir and make it
    if os.path.exists(output_filepath) and not resume:
        shutil.rmtree(output_filepath)
    os.makedirs(output_filepath,exist_ok=True)
    print("MAKE DIR:",output_filepath)
    if not (os.path.exists(os.path.join(output_filepath, "Filtered", "all_hydrogenase.txt")) and resume):
    #the last argument is the diamond exe path
        command='python '+hydrogenase_script_path+' '+ijoined_records+' '+output_filepath+' '+diamond_exe_path+" --min_pid "+str(min_pident)+" --min_qcov "+str(min_qcov)
        print("COMMAND:",command)
        os.system(command)
    #delete all the files in the hyinput dir
    #shutil.rmtree(hyinput_dir)
    #get the output file
    ijoined_records_bn=os.path.basename(ijoined_records)
    #remove the extension
    ijoined_records_bn=os.path.splitext(ijoined_records_bn)[0]
    #add the min_scov and min_pident to the output file
    #multiply the min_qcov by 100 and round
    min_qcov_name=int(min_qcov*100)
    ijoined_records_bn1=ijoined_records_bn+""+str(min_qcov_name)+""+str(min_pident)

    output_file=os.path.join(output_filepath,"Filtered",f"{ijoined_records_bn}_hydrogenase_filtered.txt")
    output_file1=os.path.join(output_filepath,"Filtered",f"{ijoined_records_bn1}_hydrogenase_filtered.txt")
    #move the output file to the tmp dir
    shutil.copy(output_file,os.path.join(tmp_dir,os.path.basename(output_file1)))
    #return the output file
    #output_file=os.path.join(tmp_dir,"all_hydrogenase.txt")
    hist,total=parse_reads(output_file,ijoined_records,lengths,hyddb=True,paired=paired)
    return hist,total









def make_reads_file(fasta_file, nreads, read_len,output_file=None,prefix="",paired=False):
    #sample reads of length read_len from the fasta file
    #open the fasta file
    with open(fasta_file, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        #create a new file with the formatted records
        formatted_records = []
        for record in records:
            formatted_records.append(record)
    reads=[]
    for record in formatted_records:
        for i in range(0,len(record.seq)-read_len,100):
            recordid=record.id+"_1:N"+str(i)+prefix
            read=record.seq[i:i+read_len]
            #reverse translate the read
            read=reverse_translate(str(read))
            reads.append(">"+recordid+"\n"+str(read))
            #if paired, add a second read
            if paired:
                recordid=record.id+"_2:N"+str(i)+prefix
                read=record.seq[i:i+read_len]
                #reverse translate the read
                read=reverse_translate(str(read))
                reads.append(">"+recordid+"\n"+str(read))
            if len(reads)==nreads:
                break
    #write the reads to a new file
    if output_file is  None:
        reads_file=fasta_file.replace('.fasta','_reads.fasta')
    else:
        reads_file=output_file
    with open(reads_file, 'a') as f:
        #write a new line
        #if first line is not empty, write a new line
        f.write("\n".join(reads))
        f.write("\n")
    return reads_file

    # Codon table dictionary
codon_table = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
    'Q': 'CAA', 'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTT', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
    'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
    'O': 'TAA'  ,'X':'TGC'# Stop codon
}

def reverse_translate(protein_sequence):
    dna_sequence = ''
    for amino_acid in protein_sequence:
        if amino_acid in codon_table:
            dna_sequence += codon_table[amino_acid]
        else:
            #add a gap ---
            dna_sequence += 'NNN'
    return dna_sequence








#batch_diamond(joined_records, diamond_exe_path, database_fasta, type, tmp_dir, min_percent_identity=90, min_scov=0.9, resume=True):

def pipeline(reads_file,args):
    tmp_folder=args.tmp_folder
    #make a folder for the reads file in the tmp folder
    reads_folder=os.path.join(tmp_folder,os.path.basename(reads_file)+"_reads")
    os.makedirs(reads_folder,exist_ok=True)
    #run the mcra diamond script
    output_filepath=os.path.join(args.output_folder,os.path.basename(reads_file)+"_output.csv")
    #make an empty with the header, gene, total_rkpm
    #if output filepath no exists, create it
    print("OUTPUT FILEPATH:",output_filepath)
    if not os.path.exists(output_filepath):
        with open(output_filepath, 'w') as f:
            f.write("gene,total_counts,total_rkpm\n")
    #else , if output filepath contains FINISHED, then return
    if "nitrogenase" in open(output_filepath).read():
        return output_filepath
    if args.hydb_fasta_file is not None:
        hist,total=batch_run_hydrogenase_script(reads_file,args.hydrogenase_script_path,args.tmp_folder,False,
                                                args.diamond_exe_path,args.hydb_fasta_file,args.resume)
        print("hyd hist",hist)
        print("hyd total",total)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            #for each key in the hist, write the key and the value
            print("write",hist)
            for key in hist:
                pass
                f.write(key+","+str(hist[key])+","+str(total[key])+"\n")
    if args.mcrA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.mcrA_database_fasta_file,"mcrA",reads_folder,
                                       args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("mcrA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the acsB diamond script
    if args.acsb_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.acsb_database_fasta_file,"acsB",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("acsB,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the dsrA diamond script
    if args.dsrA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.dsrA_database_fasta_file,"dsrA",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("dsrA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the aprA diamond script
    if args.aprA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.aprA_database_fasta_file,"aprA",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("aprA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the asrA diamond script
    if args.asrA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.asrA_database_fasta_file,"asrA",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("asrA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the dmsA diamond script
    if args.dmsA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.dmsA_database_fasta_file,"dmsA",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("dmsA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the napA diamond script
    if args.napA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.napA_database_fasta_file,"napA",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("napA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the narG diamond script
    if args.narG_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.narG_database_fasta_file,"narG",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("narG,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the nrfA diamond script
    if args.nrfA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.nrfA_database_fasta_file,"nrfA",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("nrfA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the frdA diamond script
    if args.frdA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.frdA_database_fasta_file,"frdA",reads_folder,
                                        args.resume)
        #write the output to the output file
        print("frdA rkpm",rkpm)
        print("frdA total rkpm",total_rkpm)
        with open(output_filepath, 'a') as f:
            f.write("frdA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the cydA diamond script
    if args.cydA_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.cydA_database_fasta_file,"cydA",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("cydA,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the fhl diamond script
    if args.fhl_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.fhl_database_fasta_file,"fhl",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("fhl,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the pfor diamond script
    if args.pfor_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,args.pfor_database_fasta_file,"pfor",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("pfor,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the nitrogenase diamond script
    if args.nitrogenase_database_fasta_file is not None:
        rkpm, total_rkpm=batch_diamond(reads_file,args.diamond_exe_path,
                                       args.nitrogenase_database_fasta_file,"nitrogenase",reads_folder,
                                        args.resume)
        #write the output to the output file
        with open(output_filepath, 'a') as f:
            f.write("nitrogenase,"+str(rkpm)+","+str(total_rkpm)+"\n")
    #run the hydrogenase script
    #def batch_run_hydrogenase_script(ijoined_records,hydrogenase_script_path,tmp_dir,resume=False,
    #                             diamond_exe_path=default_diamond_exe_path,hydb_fasta_file=None,paired=False,train_percent=100, min_pident=80, min_qcov=0.8):
    
    return output_filepath
        
import random
def reverse_complement(seq):
    new_seq=""
    for i in range(len(seq)-1,-1,-1):
        if seq[i]=="A":
            new_seq+="T"
        elif seq[i]=="T":
            new_seq+="A"
        elif seq[i]=="G":
            new_seq+="C"
        elif seq[i]=="C":
            new_seq+="G"
    #reverse the sequence
    return new_seq

def evaluation_script(sim=True,repeat=1, max_reads=100,percent_junk=0.9,train_test_split=0.8,junk=False,marker=True,hyd=True,min_hyd_pidents=[80],min_hyd_qcovs=[0.8]):
    root_folder=r'C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads'
    
    default_root_folder=r'/tllhome/abel/hydrogen_consumer_pipeline/test_reads'
    #if root folder not exists, use the default root folder
    if not os.path.exists(root_folder):
        root_folder=default_root_folder
    #train test split records will be 
    #max_reads=100
    #percent_junk=0.9
    reads_file=os.path.join(root_folder,"reads.fasta")
    if sim:
        default_hyddb_fasta_file=r'/tllhome/abel/hydrogen_consumer_pipeline/hyDB/HydDB_reformated.fasta'
        hyddb_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyDB\HydDB_reformated.fasta"
        #if hyddb_fasta_file does not exist, use the default hyddb_fasta_file
        if not os.path.exists(hyddb_fasta_file):
            hyddb_fasta_file=default_hyddb_fasta_file
        hyd_groups=defaultdict(list)
        rec_dict={}
        #set the random_seed
        random.seed(0)
        #hyd=False
        all_reads=[]
        gene_stats=defaultdict(lambda : [0,0,''])  
        hist=defaultdict(int)
        if hyd:
            #open the hyddb fasta file and extract the hyd groups
            with open(hyddb_fasta_file, 'r') as f:
                records = list(SeqIO.parse(f, 'fasta'))
                for record in records:
                    #get the hyd group
                    hyd_group=record.id.split("_-_[")[1]
                    hyd_group=hyd_group.replace('[','').replace(']','')
                    hyd_groups[hyd_group].append(record.id)
                    rec_dict[record.id]=record
            #make a random multinomial distribution of reads
            #for each hyd group, make a random number of reads
            
            #rkpm[subj]=num_reads/(gene_length/1000*total_reads/1000000)
            
            for hyd_grp in hyd_groups:
                #generate a random number between 0 nad max number of reads
                num_reads=random.randint(0,max_reads)
                hist[hyd_grp]=num_reads
            
            for hyd_grp in hist:
                print(hyd_grp,"has",hist[hyd_grp],"reads")
                #sample reads
                reads=[]
                for i in range(hist[hyd_grp]):
                    #sample a random read
                    random_rec=random.choice(hyd_groups[hyd_grp])
                    #sample a random read from the record
                    record=rec_dict[random_rec]
                    gene_stats[record.id][0]=len(record.seq)*3
                    gene_stats[record.id][1]+=1
                    gene_stats[record.id][2]=hyd_grp
                    recseq=reverse_translate(str(record.seq))
                    rand_start=random.randint(0,len(recseq)-100)
                    read=recseq[rand_start:rand_start+100]
                    reads.append(">"+random_rec+"1:N"+str(rand_start)+"\n"+str(read))
                    
                    #make a second read
                    rand_start2=random.randint(rand_start,len(recseq)-100)
                    read= recseq[rand_start2:rand_start2+100]
                    #reverse complement the read
                    read=reverse_complement(read)

                    reads.append(">"+random_rec+"2:N"+str(rand_start)+"\n"+str(read))
                    
                all_reads+=reads
        #add reads from the marker database, they are in the /tllhome/abel/hydrogen_consumer_pipeline/DIAMOND dir
        
                              
        if marker:
            marker_files={"mcrA":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\mcrA.fasta",
                        "acsB":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\acsB.fasta",
                            "dsrA":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\dsrA.fasta",
                            "aprA":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\aprA.fasta",
                            "asrA":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\asrA.fasta",
                            "dmsA":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\dmsA.fasta",
                            "napA":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\napA.fasta",
                            "narG":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\narG.fasta",
                            "nrfA":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\nrfA.fasta",
                            "frdA":r"c:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\frdA.fasta",
                            "cydA":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\cydA.fasta",
                            "fhl":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\fhl.fasta",
                            "pfor":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\pfor.fasta",
                            "nitrogenase":r"C:\Users\abel\Documents\hydrogenases\curated_marker\final_curated_marker\nifH.fasta"}
            #if the marker files do not exist, use the default marker files
            #if any ([not os.path.exists(marker_files[marker]) for marker in marker_files]):
            default_marker_files={"mcrA":r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/mcrA.fasta",
                                  "acsB": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/acsB.fasta",
                                    "dsrA": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/dsrA.fasta",
                                    "aprA": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/aprA.fasta",
                                    "asrA": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/asrA.fasta",
                                    "dmsA": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/dmsA.fasta",
                                    "napA": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/napA.fasta",
                                    "narG": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/narG.fasta",
                                    "nrfA": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/nrfA.fasta",
                                    "frdA": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/frdA.fasta",
                                    "cydA": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/cydA.fasta",
                                    "fhl": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/fhl.fasta",
                                    "pfor": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/pfor.fasta",
                                    "nitrogenase": r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/nitrogenase.fasta"}
            
            for marker in marker_files:
                #if no exists, use the default marker files
                if not os.path.exists(marker_files[marker]):
                    print("Using default marker files", marker_files[marker])
                    marker_files[marker]=default_marker_files[marker]

            marker_gene_stats=defaultdict(lambda : [0,0,''])
            for marker in marker_files:
                    #make a random number of reads
                    hist[marker]=random.randint(1,max_reads)
                    print(marker,"has",hist[marker],"reads")
                    #open the fasta file
                    with open(marker_files[marker], 'r') as f:
                        records = SeqIO.parse(f, 'fasta')
                        #create a new file with the formatted records
                        formatted_records = []
                        for record in records:
                            formatted_records.append(record)
                    reads=[]
                    #randomly choose a record from the formatted records
                    while(len(reads)<hist[marker]*2):
                        choice=random.choice(formatted_records)
                        record=choice
                        recseq=reverse_translate(str(record.seq))

                        rand_start=random.randint(0,len(recseq)-100)
                        marker_gene_stats[record.id][0]=len(recseq)*3
                        read=recseq[rand_start:rand_start+100]
                        reads.append(">"+record.id+"1:N"+str(rand_start)+"_"+marker+"\n"+str(read))
                        #make a second read
                        rand_start2=random.randint(rand_start,len(recseq)-100)
                        read= recseq[rand_start2:rand_start2+100]
                        #reverse complement the read
                        read=reverse_complement(read)
                        reads.append(">"+record.id+"2:N"+str(rand_start)+"_"+marker+"\n"+str(read))
                        #update the gene stats
                        marker_gene_stats[record.id][1]+=1
                        marker_gene_stats[record.id][2]=marker
                    all_reads+=reads



            

        if junk:
            #junk_genomes_dir=r'C:\Users\abel\Documents\mouse_gut\mouse_gut_catalogue\mgnify_,mouse_catalogue\accepted_genomes'
            #default_junk_genomes_dir=r'/tllhome/abel/mgnify_mouse_catalogue/hq_genomes'
            junk_folders=[("non_consumer",r"C:\Users\abel\Documents\literature_search\literature_review\genomes\non_hydrogen_consumer"),
                          ("non_producer",r"C:\Users\abel\Documents\literature_search\literature_review\genomes\non_hydrogen_producer")]
            #if junk genome drive not exist, use  /tllhome/abel/mgnify_mouse_catalogue/hq_genomes
            #if not os.path.exists(junk_genomes_dir):
            #    junk_genomes_dir=r'/tllhome/abel/mgnify_mouse_catalogue/hq_genomes'
            junk_genome_files=[]
            for folder in junk_folders:
                #if not os.path.exists(folder[1]):
                #    folder[1]=os.path.join(default_junk_genomes_dir,folder[0])
                junk_genome_files+=[os.path.join(folder[1],f) for f in os.listdir(folder[1]) if f.endswith('.fna')]
            #choose a random fasta file from the junk genomes dir and sample reads from it
            nu_junk=(percent_junk*(1/(1-percent_junk)))*len(all_reads)
            for junk in range(0,int(nu_junk),100):
                #choose a random fasta file from the junk genomes dir
                junk_genome=random.choice(junk_genome_files)
                #open the fasta file
                with open(junk_genome, 'r') as f:
                    #if not fasta or fa, skip
                    if not any([junk_genome.endswith(ext) for ext in ['fasta','fa',"fna"]]):
                        continue
                    records = list(SeqIO.parse(f, 'fasta'))
                    #if file name contains producer, tag is Prod, else tag is Cons
                    if "producer" in junk_genome:
                        tag="Prod"
                    else:
                        tag="Cons"
                    for batch in range(100):
                        record=random.choice(records)
                        rand_start=random.randint(0,len(record.seq)-100)
                        read=record.seq[rand_start:rand_start+100]
                        recordid=record.id.split("_-_")[0]
                        all_reads.append(">"+recordid+"1_N"+str(rand_start)+"_"+tag+"\n"+str(read))
                        #make a second read
                        rand_start2=random.randint(rand_start,len(record.seq)-100)
                        read= record.seq[rand_start2:rand_start2+100]
                        #reverse complement the read
                        read=reverse_complement(read)
                        all_reads.append(">"+recordid+"2_N"+str(rand_start2)+"_"+tag+"\n"+str(read))
                        if batch+junk>=nu_junk:
                            break
        #with open(r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\rkpm_formula_gold.txt", 'w') as f:
        #    f.write("rkpm formula\n")
        hyd_grp_rkpms=defaultdict(float)
        for gene in gene_stats:
            #rkpm = 10^9 * read_mapped/(gene_length*total_reads)
            hyd_grp_rkpms[gene_stats[gene][2]]+=10**9*(gene_stats[gene][1]/(gene_stats[gene][0]*len(all_reads)))
            #write the rkpm formula to a file
            #with open(r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\rkpm_formula_gold.txt", 'a') as f:
            #    f.write(f"rkpm{gene}=10^9*{gene_stats[gene][1]}/({gene_stats[gene][0]}*{len(all_reads)} = {10**9*(gene_stats[gene][1]/(gene_stats[gene][0]*len(all_reads)))})\n")
        if marker:
            marker_rkpms=defaultdict(float)
            for gene in marker_gene_stats:
                #rkpm = 10^9 * read_mapped/(gene_length*total_reads)
                marker_rkpms[marker_gene_stats[gene][2]]+=10**9*(marker_gene_stats[gene][1]/(marker_gene_stats[gene][0]*len(all_reads)))
                #write the rkpm formula to a file
                continue
            #with open(r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\rkpm_formula_gold.txt", 'a') as f:
            #    f.write(f"rkpm{gene}=10^9*{marker_gene_stats[gene][1]}/({marker_gene_stats[gene][0]}*{len(all_reads)} = {10**9*(marker_gene_stats[gene][1]/(marker_gene_stats[gene][0]*len(all_reads)))})\n")
        #write the reads to a file
        #reads_file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\hyinput\pipeline_hyd_reads.fasta"
        #if hyinput does not exist, make it
        if not os.path.exists(os.path.join(root_folder,"hyinput")):
            os.makedirs(os.path.join(root_folder,"hyinput"),exist_ok=True)
        reads_file=os.path.join(root_folder,"hyinput","pipeline_hyd_reads.fasta")
        #copy it to the hyinput
        
        #run the diamon dblast script
        with open(reads_file, 'w') as f:
            f.write("\n".join(all_reads))
            f.write("\n")
        #shutil.copy(reads_file,r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\hyinput")
    #print the gene stats to a file
    #with open(r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\gene_stats.txt", 'w') as f:
        #write the header
        #PRINT THE TOTAL NUMBER OF READS
        #f.write("Total number of reads,"+str(len(all_reads))+"\n")
        #f.write("gene,length,num_reads,group\n")
        #for gene in gene_stats:
        #    f.write(gene+","+",".join([str(x) for x in gene_stats[gene]])+"\n")
    
    #run the diamond blast
    hyd_script_path=r"C:\Users\abel\Documents\hydrogenases\hyDB\Diamond_blast_hyDB_reads.py"
    #if the hyd_script_path does not exist, use the default hyd_script_path
    if not os.path.exists(hyd_script_path):
        hyd_script_path=r"/tllhome/abel/hydrogen_consumer_pipeline/hyDB/Diamond_blast_hyDB_reads.py"
    tmp_dir=os.path.join(root_folder,"tmp")
    #if tmp dir exist, remove it
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir,exist_ok=True)
    diamond_exe_path=r"D:\diamond-windows\diamond.exe"
    #if the diamond_exe_path does not exist, use the default diamond_exe_path
    #sftp://abel@hpc.tll.org.sg/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/diamond
    if not os.path.exists(diamond_exe_path):
        diamond_exe_path=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/diamond"
    hyddb_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyDB\HydDB_reformated.fasta"
    #if the hyddb_fasta_file does not exist, use the default hyddb_fasta_file
    if not os.path.exists(hyddb_fasta_file):
        hyddb_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/hyDB/HydDB_reformated.fasta"
    if hyd:
        for min_pident_hyd in min_hyd_pidents:
            for min_qcov_hyd in min_hyd_qcovs:
                

                hist1,total=batch_run_hydrogenase_script(reads_file, hyd_script_path, tmp_dir, resume=False, diamond_exe_path=diamond_exe_path, hydb_fasta_file=hyddb_fasta_file,
                                                         paired=True,train_percent=train_test_split*100,min_pident=min_pident_hyd,min_qcov=min_qcov_hyd)
                
                #get the rmsd of total and hyd_group_rkpms
                group_counts={}
                for gs in gene_stats:
                    group=gene_stats[gs][2]
                    if group not in group_counts:
                        group_counts[group]=0
                    group_counts[group]+=gene_stats[gs][1]
                rmsd=0
                devs=defaultdict(list)
                for key in hyd_grp_rkpms:
                    gold_count=group_counts[key]
                    gold_rkpm=hyd_grp_rkpms[key]
                    observed_count=hist1[key]
                    observed_rkpm=total[key]
                    dev=gold_count-observed_count
                    dev_rkpm=gold_rkpm-observed_rkpm
                    devs[key]=[gold_count,observed_count,gold_rkpm,observed_rkpm,dev,dev_rkpm]
                rmsd=(rmsd/len(devs))**0.5
                #save the devs to a file
                #with open(fr"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\_{percent_junk}_{len(all_reads)}_{repeat}_devs.txt", 'w') as f:
                with open (os.path.join(root_folder,f"_{percent_junk}_{len(all_reads)}_{repeat}_{min_pident_hyd}_{min_qcov_hyd}_devs.csv"), 'w') as f:
                    #write the header
                    f.write("hyd_group,gold_count,observed_count,gold_rkpm,observed_rkpm,dev,dev_rkpm\n")
                    for key in devs:
                        f.write(key+","+",".join([str(x) for x in devs[key]])+"\n")
    #for each marker, run the diamond script
    if (marker):
        devs=defaultdict(list)
        group_counts={}
        for gs in marker_gene_stats:
                group=marker_gene_stats[gs][2]
                if group not in group_counts:
                    group_counts[group]=0
                group_counts[group]+=marker_gene_stats[gs][1]
        for marker in marker_files:
            count1,rkpm1=batch_diamond(reads_file,diamond_exe_path,marker_files[marker],marker,tmp_dir,paired=True,train_percent=train_test_split*100)
            #get the rmsd of total and hyd_group_rkpms
            
            
            rmsd=0
            
            key=marker
            gold_count=group_counts[key]
            gold_rkpm=marker_rkpms[key]
            observed_count=count1
            observed_rkpm= rkpm1
            dev=gold_count-observed_count
            dev_rkpm=gold_rkpm-observed_rkpm
            devs[key]=[gold_count,observed_count,gold_rkpm,observed_rkpm,dev,dev_rkpm]
            rmsd=(rmsd/len(devs))**0.5
            #save the devs to a file
        #with open(fr"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\marker_devs_{percent_junk}_{len(all_reads)}_{repeat}.txt", 'w') as f:
        with open (os.path.join(root_folder,f"marker_devs_{percent_junk}_{len(all_reads)}_{repeat}.csv"), 'w') as f:
            #write the header
            f.write("marker,gold_count,observed_count,gold_rkpm,observed_rkpm,dev,dev_rkpm\n")
            for key in devs:
                f.write(key+","+",".join([str(x) for x in devs[key]])+"\n")
    return 
    #open the row rkpm formula file as df
    #rkpm_observed=pd.read_csv(r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\rkpm_formula.txt",sep=",")
    #rkpm_gold=pd.read_csv(r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\rkpm_formula_gold.txt",sep=",")
    #for each row in rkpm observed, find the row that sarts with the same header
    headers_observed={}
    new_rows=[]
    for index,row in rkpm_observed.iterrows():
        header=row[0].split("_-_")[1]
        headers_observed[header]=row[0]
    for index ,row in rkpm_gold.iterrows():
        header=row[0].split("_-_")[1]
        if header not in headers_observed:
            continue
        observed_row=headers_observed[header]
        #get the rkpm values
        rkpm_gold= row[0].split("=")[-1]
        rkpm_observed=observed_row.split("=")[-1]
        re_gold=re.compile(r"10\^9\*\d+")
        count_gold=re_gold.findall(row[0])[0]
        #remove the 10^9*
        count_gold=count_gold.replace("10^9*","")
        re_observed=re.compile(r"1e9\*\d+")
        count_observed=re_observed.findall(observed_row)[0]
        count_observed=count_observed.replace("1e9*","")

        #if both rkpm are have the same integer value then they are equal
        #remove brackets and strip
        rkpm_gold=rkpm_gold.replace("(","").replace(")","").strip()
        rkpm_observed=rkpm_observed.replace("(","").replace(")","").strip()
        equals= int(float(rkpm_gold))==int(float(rkpm_observed))
        new_row=[header,count_gold,count_observed,equals]
        new_rows.append(new_row)
    #sort by equals, not equals first
    new_rows=sorted(new_rows,key=lambda x: x[3])
    #write new rows to a file
    with open(r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\rkpm_comparison.txt", 'w') as f:
        f.write("header,count_gold,count_observed,equals\n")
        for row in new_rows:
            f.write(",".join([str(x) for x in row])+"\n")


    return rmsd








def tes():
    diamond_exe_path=r"D:\diamond-windows\diamond.exe"
    database_fasta=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\dsrA.fasta"
    tmp_dir=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads"
    reads_file=make_reads_file(database_fasta, 100, 100)
    joined_records=reads_file
    batch_diamond(joined_records, diamond_exe_path, database_fasta, "test", tmp_dir, min_percent_identity=90, min_scov=0.9, resume=False)
    parse_reads(os.path.join(tmp_dir,"test_diamond_output_filtered.tsv"),reads_file)
def tes2():
    hydb_fasta_file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\HydDB_reformated.fasta"
    diamond_exe_path=r"D:\diamond-windows\diamond.exe"
    hyd_script_path=r"C:\Users\abel\Documents\hydrogenases\hyDB\Diamond_blast_hyDB_reads.py"
    tmp_dir=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads"
    reads_file=make_reads_file(hydb_fasta_file, 100, 100)
    batch_run_hydrogenase_script(reads_file, hyd_script_path, tmp_dir, resume=False, diamond_exe_path=diamond_exe_path, hydb_fasta_file=hydb_fasta_file)

import argparse
def tes4():
    args=argparse.Namespace()
    args.tmp_folder=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\tmp"
    args.output_folder=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads"
    args.diamond_exe_path=r"D:\diamond-windows\diamond.exe"
    args.mcrA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\methanogen\mcrA.fasta"
    args.acsb_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\acetogen\acsb.fasta"
    args.dsrA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\dsrA.fasta"
    args.aprA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\aprA.fasta"
    args.asrA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\asrA.fasta"
    args.dmsA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\trimethylamine_N_oxide_reduction\DmsA.fasta"
    args.napA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\napA.fasta"
    args.narG_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\narG.fasta"
    args.nrfA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\nrfA.fasta"
    args.frdA_database_fasta_file=r"c:\Users\abel\Documents\hydrogenases\fumarate_reducer\frdA.fasta"
    args.cydA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\cydA.fasta"
    args.fhl_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyd_production\FHL.fasta"
    args.pfor_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyd_production\PFOR.fasta"
    args.nitrogenase_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyd_production\nifH.fasta"
    args.hydb_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyDB\HydDB_reformated.fasta"
    #make a reads file
    reads_file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\pipeline_reads.fasta"
    reads_file=make_reads_file(args.hydb_fasta_file, 100, 100,reads_file)
    reads_file=make_reads_file(args.mcrA_database_fasta_file, 100, 100,reads_file,prefix="_mcrA")
    reads_file=make_reads_file(args.acsb_database_fasta_file, 100, 100,reads_file,prefix="_acsb")
    reads_file=make_reads_file(args.dsrA_database_fasta_file, 100, 100,reads_file,prefix="_dsrA")
    reads_file=make_reads_file(args.aprA_database_fasta_file, 100, 100,reads_file,prefix="_aprA")
    reads_file=make_reads_file(args.asrA_database_fasta_file, 100, 100,reads_file,prefix="_asrA")
    reads_file=make_reads_file(args.dmsA_database_fasta_file, 100, 100,reads_file,prefix="_dmsA")
    reads_file=make_reads_file(args.napA_database_fasta_file, 100, 100,reads_file,prefix="_napA")
    reads_file=make_reads_file(args.narG_database_fasta_file, 100, 100,reads_file,prefix="_narG")
    reads_file=make_reads_file(args.nrfA_database_fasta_file, 100, 100,reads_file,prefix="_nrfA")
    reads_file=make_reads_file(args.frdA_database_fasta_file, 100, 100,reads_file,prefix="_frdA")
    reads_file=make_reads_file(args.cydA_database_fasta_file, 100, 100,reads_file,prefix="_cydA")
    reads_file=make_reads_file(args.fhl_database_fasta_file, 100, 100,reads_file,prefix="_fhl")
    reads_file=make_reads_file(args.pfor_database_fasta_file, 100, 100,reads_file,prefix="_pfor")
    reads_file=make_reads_file(args.nitrogenase_database_fasta_file, 100, 100,reads_file,prefix="_nitrogenase")
    pipeline(reads_file,args)

def tes_paired():
    args=argparse.Namespace()
    args.tmp_folder=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\tmp"
    args.output_folder=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads"
    args.diamond_exe_path=r"D:\diamond-windows\diamond.exe"
    args.mcrA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\methanogen\mcrA.fasta"
    args.acsb_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\acetogen\acsb.fasta"
    args.dsrA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\dsrA.fasta"
    args.aprA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\aprA.fasta"
    args.asrA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\asrA.fasta"
    args.dmsA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\trimethylamine_N_oxide_reduction\DmsA.fasta"
    args.napA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\napA.fasta"
    args.narG_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\narG.fasta"
    args.nrfA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\nrfA.fasta"
    args.frdA_database_fasta_file=r"c:\Users\abel\Documents\hydrogenases\fumarate_reducer\frdA.fasta"
    args.cydA_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\cydA.fasta"
    args.fhl_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyd_production\FHL.fasta"
    args.pfor_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyd_production\PFOR.fasta"
    args.nitrogenase_database_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyd_production\nifH.fasta"
    args.hydb_fasta_file=r"C:\Users\abel\Documents\hydrogenases\hyDB\HydDB_reformated.fasta"
    args.hydrogenase_script_path=r"C:\Users\abel\Documents\hydrogenases\hyDB\Diamond_blast_hyDB_reads.py"
    #make a reads file
    reads_file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test_reads\pipeline_reads.fasta"
    reads_file=make_reads_file(args.hydb_fasta_file, 100, 100,reads_file,paired=True)
    #reads_file=make_reads_file(args.mcrA_database_fasta_file, 100, 100,reads_file,prefix="_mcrA",paired=True)
    #batch diamond the reads
    #batch_diamond(reads_file, args.diamond_exe_path, args.mcrA_database_fasta_file, "mcrA",
    #               args.tmp_folder, min_percent_identity=90, min_scov=0.9, resume=True,paired=True)
    #run the hydrogenase script
    diamond_path=r"D:\diamond-windows\diamond.exe"
    batch_run_hydrogenase_script(reads_file, args.hydrogenase_script_path,args.tmp_folder, resume=True, diamond_exe_path=diamond_path, hydb_fasta_file=args.hydb_fasta_file,paired=True)

    


import csv
def main():
    #this is meant to be run on kylin
    #the database dir is at /oceanstor/home/e1103389/hydrogenases
    args=argparse.ArgumentParser()
    args.add_argument('tmp_folder', help='Enter the tmp folder')
    args.add_argument('output_folder', help='Enter the output folder')
    #add argument for the input folder
    args.add_argument('input_folder', help='Enter the input folder')
    #add a resume argument, it is store as a boolean
    args.add_argument('--resume', help='Enter True if you want to resume the pipeline', default=False,required=False,action='store_true')
    args.add_argument('--diamond_exe_path', help='Enter the path to the diamond executable', default='/oceanstor/home/e1103389/DIAMOND/diamond',required=False)
    args.add_argument('--mcrA_database_fasta_file', help='Enter the path to the mcrA database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/methanogen/mcrA.fasta')
    args.add_argument('--acsb_database_fasta_file', help='Enter the path to the acsb database fasta file',required=False,
                       default='/oceanstor/home/e1103389/hydrogenases/acetogen/acsB.fasta')
    args.add_argument('--dsrA_database_fasta_file', help='Enter the path to the dsrA database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/sulphate_reducer/dsrA.fasta')
    args.add_argument('--aprA_database_fasta_file', help='Enter the path to the aprA database fasta file',required=False,
                       default='/oceanstor/home/e1103389/hydrogenases/sulphate_reducer/aprA.fasta')
    args.add_argument('--asrA_database_fasta_file', help='Enter the path to the asrA database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/sulphate_reducer/asrA.fasta')
    args.add_argument('--dmsA_database_fasta_file', help='Enter the path to the dmsA database fasta file',required=False,
                       default='/oceanstor/home/e1103389/hydrogenases/trimethylamine_N_oxide_reduction/DmsA.fasta')
    args.add_argument('--napA_database_fasta_file', help='Enter the path to the napA database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/nitrate_reducer/napA.fasta')
    args.add_argument('--narG_database_fasta_file', help='Enter the path to the narG database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/nitrate_reducer/narG.fasta')
    args.add_argument('--nrfA_database_fasta_file', help='Enter the path to the nrfA database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/nitrate_reducer/nrfA.fasta')
    args.add_argument('--frdA_database_fasta_file', help='Enter the path to the frdA database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/fumarate_reducer/frdA.fasta')
    args.add_argument('--cydA_database_fasta_file', help='Enter the path to the cydA database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/cydA.fasta')
    args.add_argument('--fhl_database_fasta_file', help='Enter the path to the fhl database fasta file', required=False,
                      default='/oceanstor/home/e1103389/hydrogenases/hyd_production/FHL.fasta')
    args.add_argument('--pfor_database_fasta_file', help='Enter the path to the pfor database fasta file', default='/oceanstor/home/e1103389/hydrogenases/hyd_production/PFOR.fasta')
    args.add_argument('--nitrogenase_database_fasta_file', help='Enter the path to the nitrogenase database fasta file', default='/oceanstor/home/e1103389/hydrogenases/hyd_production/nifH.fasta')
    args.add_argument('--hydb_fasta_file', help='Enter the path to the hydb database fasta file', default='/oceanstor/home/e1103389/hyDB/HydDB_reformated.fasta')
    #add the hydrogenase script path
    args.add_argument('--hydrogenase_script_path', help='Enter the path to the hydrogenase script', 
                      default='/oceanstor/home/e1103389/hydrogenases/hyDB/Diamond_blast_hyDB_reads.py')
    #add an argument for log file
    args.add_argument('--log_file', help='Enter the path to the log file', default=None)
    args=args.parse_args()
    #if log file is none, set it to input folder
    if args.log_file is None:
        args.log_file=os.path.join(args.input_folder,"log.txt")
    #make the output folder if it does not exist
    os.makedirs(args.output_folder,exist_ok=True)
    #make the tmp folder if it does not exist
    os.makedirs(args.tmp_folder,exist_ok=True)
    #open the hydrogenase fasta file and load all the hyd groups
    with open(args.hydb_fasta_file, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        hyd_groups=[]
        for record in records:
            ide=record.description
            print(ide)
            grp=ide.split("_-_[")[1]
            grp=grp.replace('[','').replace(']','')
            hyd_groups.append(grp)
    hyd_groups=set(hyd_groups)
    all_keys=hyd_groups
    #log_file=os.path.join(args.output_folder,"log.txt")
    for file in os.listdir(args.input_folder):
        ext=file.split(".")[-1]
        if ext in["fasta","fa","fna","fastq","fq"]:
            reads_file=os.path.join(args.input_folder,file)
            #remove internal dots
            name=os.path.basename(reads_file)
            ext=name.split(".")[-1]
            name=name.replace(".","_")
            new_name=name+"."+ext
            #copy the reads file to the tmp folder
            shutil.copy(reads_file,os.path.join(args.tmp_folder,new_name))
            print("copied",reads_file)
    for file in os.listdir(args.tmp_folder):
        ext=file.split(".")[-1]
        if ext in["fasta","fa","fna","fastq","fq"]:
            reads_file=os.path.join(args.tmp_folder,file)
            #if resume and the file is in the output folder, skip

            output_filepath=pipeline(reads_file,args)
            #append output filepath to log file
            with open(args.log_file, 'a') as f:
                f.write(output_filepath+"\n")
            if len(all_keys)==0:
                with open(output_filepath, 'r') as f:
                    for line in f:
                        if line.startswith("gene"):
                            continue
                        all_keys.add(line.split(",")[0])
                        print("added",line.split(",")[0])
    #concatenate all the output files in output folder
    output_file=os.path.join(args.output_folder,"all_output.csv")
    all_samples_dict=[]
    all_samples_counts=[]
    all_keys=[]
    for output_filen in os.listdir(args.output_folder):
        sample_name=output_filen.split(".")[0]
        row={"sample":sample_name}
        row2={"sample":sample_name}
        if output_filen.endswith('.csv'):
            with open(os.path.join(args.output_folder,output_filen), 'r') as f:
                print("reading",output_filen)
                next(f)
                for line in f:
                    gene=line.split(",")[0]
                    count=line.split(",")[1]
                    rkpm=line.split(",")[2]
                    row[gene]=rkpm
                    row2[gene]=count
                    all_keys.append(gene)
            all_samples_dict.append(row)
            all_samples_counts.append(row2)
    all_keys=["sample"]+list(set(all_keys))
    with open(output_file, 'w') as f:
        #write the header
        writer=csv.DictWriter(f,all_keys)
        writer.writeheader()
        pad_dict={}
        for sample_row in all_samples_counts:
            for key in all_keys:
                if key not in sample_row:
                    pad_dict[key]=0
            sample_row.update(pad_dict)
            writer.writerow(sample_row)
    #make the rkpm out
    rkpm_filename=output_file.replace(".csv","_rkpm.csv")
    with open(rkpm_filename, 'w') as f:
        #write the header
        writer=csv.DictWriter(f,all_keys)
        writer.writeheader()
        pad_dict={}
        for sample_row in all_samples_dict:
            for key in all_keys:
                if key not in sample_row:
                    pad_dict[key]=0
            sample_row.update(pad_dict)
            writer.writerow(sample_row)

    #pipeline(args)



if __name__=="__main__":
    #evaluation_script(sim=True,train_test_split=1,repeat=1,percent_junk=0.9,hyd=True,marker=False,junk=True)
    main()