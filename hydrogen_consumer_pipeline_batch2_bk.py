


import argparse
import concurrent
import pandas as pd
import csv
import os
import random
import shutil
from collections import defaultdict
import create_feature_table as CFT
import fix_feature_table as FFT
import numpy as np
from Bio import Entrez,SeqIO
import csv
import re
import os

from hydrogen_consumer_pipeline_deuplicate_names import deduplicate_names

acetobase_dir=r'/home/ec2-user/acetoscan/AcetoBase.fasta'
acetoscan_input_dir=r'/home/ec2-user/acetoscan/genome'
#hyDb_input_dir=r'/home/ec2-user/hyDB/genome'
#hyDb_output_dir=r'/home/ec2-user/hydrogenae_Scan/genome_output'
# diamond_exe_path,hydrogenase_script_path,
#             dsr_database_path, mcr_database_path
default_diamond_exe_path=r'/home/ec2-user/hyDB/diamond'
default_hydrogenase_script_path=r'/home/ec2-user/hyDB/Diamond_blastp_hyDB.sh'
default_dsr_database_path=r'/home/ec2-user/hyDB/dsrA.dmnd'
default_mcr_database_path=r'/home/ec2-user/hyDB/mcr1.dmnd'
default_prodigal_path=r'/home/ec2-user/prodigal/Prodigal/prodigal'
dsr_database_fasta_file=r'/home/ec2-user/dsr/dsrA_curated.fasta'
mcr_database_fasta_file=r'/home/ec2-user/mcr/mcr1_curated2.fasta'
MCR_MIN_PERCENT_IDENTITY=90
MCR_MIN_SUBJECT_COVERAGE=0.9
DSR_MIN_PERCENT_IDENTITY=50
DSR_MIN_SUBJECT_COVERAGE=0.9


CAT_taxonomy_path="/home/ec2-user/nr_dmnd/CAT_data/CAT_prepare_20210107/2021-01-07_taxonomy"
CAT_database_path="/home/ec2-user/nr_dmnd/CAT_data/CAT_prepare_20210107/2021-01-07_CAT_database"
#the create feature script is in the same directory as this script

DEFAULT_HYDROGENASE_ACTIVITY_FILE='/home/ec2-user/hyDB/hydrogenase_classes_activity2t.csv'
class BatchFetcher:
    #fetches the batch of accessions from NCBI
    def __init__(self,database):
        self.accessions=set()
        self.database=database
        self.descriptions={}
    def add_accession(self,accession):
        self.accessions.add(accession)
        uid="{{"+accession+"}}"
        return uid
    def fetch(self):
        #fetch the accessions
        Entrez.email = "abel@tll.org.sg"
        descriptions={}
        #fetch all the accessions in batches of 100
        for i in range(0,len(self.accessions),100):
            batch_accessions=list(self.accessions)[i:i+100]
            print("Fetching from NCBI for batch", i)
            print(f"Entrez.efetch(db={self.database}, id={batch_accessions}, rettype='fasta', retmode='text')")
            #try 5 times, if uncessful, raise the error
            tries=0
            while tries<5:
                try:
                    print("FETCH try:",tries)
                    handle = Entrez.efetch(db=self.database, id=batch_accessions, rettype="fasta", retmode="text")
                    records = SeqIO.parse(handle, "fasta")

                    for record in records:
                        descriptions[record.id] = record.description
                    handle.close()
                    break
                except Exception as e:

                    tries+=1
            if tries==5:
                print("Error occurred for batch", i)
                print(e)
                continue

        #look for the accessions that were not found
        for accession in self.accessions:
            if accession not in descriptions:
                descriptions[accession]=accession

        self.descriptions=descriptions
        return descriptions
    def replace(self, string):
        #replaces all the uids with the descriptions
        #print("Parsing",string)
        for accession in self.accessions:
            uid="{{"+accession+"}}"
            repl=self.descriptions.get(accession,accession)
            string=string.replace(uid,repl)
        return string

Batch_Fetcher=BatchFetcher("protein")
def fetch_from_ncbi(accession,database):
    Entrez.email = "abel@tll.org.sg"
    tries=0
    uid=Batch_Fetcher.add_accession(accession)
    return uid
    #return "Unknown"
    while tries<3:
        try:
            # get the accession record from the database
            print("Fetching from NCBI for accession", accession)
            handle = Entrez.efetch(db=database, id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            #get the description
            desc=record.description
            handle.close()
            return desc
        except Exception as e:
            print("Error occurred for record", accession)
            tries+=1
            if tries==3:
                raise e
            print(e)
            return None

def is_sulphate_reducer(dsr_diamond_output):
    with open(dsr_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(dsr_diamond_output)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}


def is_methanogen(mcr_diamond_output):
    with open(mcr_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(mcr_diamond_output)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}



def process_diamond_output(diamond_output_file,bit_score_index=-1,inject_desc=True):
    #open the file,for each query, find the subject with the highest bit score.
    #the bit score is the last column
    unique_queries=defaultdict(list)
    with open(diamond_output_file) as f:
        reader=csv.reader(f,delimiter='\t')
        for row in reader:
            #if len of row is less than 3, skip
            if len(row)<3:
                continue
            query=row[0]
            subject=row[1]
            bit_score=float(row[bit_score_index])
            unique_queries[query].append((subject,bit_score))
    #for each query, sort the subjects by bit score
    for query,subject_list in unique_queries.items():
        subject_list.sort(key=lambda x:x[1],reverse=True)
        unique_queries[query]=subject_list[0][0]
    #fetch the description from ncbi
    for query,subject in unique_queries.items():
        #if subject contains " - ", split it and take the first part
        print("subject:",subject)
        subject1=subject
        if " - " in subject:
            subject1=subject.split(" - ")[0]
            #split the subject by - and take the 2nd part
            if "-" in subject1:
                subject1=subject.split("-")[1]

        desc=fetch_from_ncbi(subject1,'protein')
        if inject_desc==False:
            unique_queries[query]=subject
        else:
            unique_queries[query]=desc
    return unique_queries


#use the acetoscan out.csv file
def process_acetoscan_output(acetoscan_output):
    '''the output looks like this, retrieve the species name
    Subject_Accession	Kingdom	Phylum	Class	Order	Family	Genus	Species	Percentage_identity	Evalue	Query_length	Alignment_length[blastx]	Bitscore	Query_seq
NP_0000006450	Bacteria	Bacteroidetes	Bacteroidia	Bacteroidales	Prevotellaceae	Prevotella	Prevotella_sp	45.833	0.00000461	226	48	40.8	KEAKMSKAAKITAICNQKGGVGKTVTTVNLGIGLAREGKKVLLVDVDP
NP_0000016492	Bacteria	Actinobacteria	Actinobacteria	Propionibacteriales	Propionibacteriaceae	Acidipropionibacterium	Acidipropionibacterium_thoenii	33.824	0.0000124	400	68	42	HCPITVASFEAGKHVMCEKPMAHNTEDAQKMLDAWKKSGKKFTIGYQNRLRDDTQTLHASCEAGELGD
NP_0000006450	Bacteria	Bacteroidetes	Bacteroidia	Bacteroidales	Prevotellaceae	Prevotella	Prevotella_sp	47.368	0.000727	436	38	37	VIALANQKGGTGKTTTAVNLGVGLANEGKKVLLVDADP
'''
    unique_queries= defaultdict(list)
    print("processing acetoscan output",acetoscan_output)
    with open(acetoscan_output) as f:
        #ITS A COMMA SEPARATED FILE
        reader=list(csv.reader(f,delimiter=','))
        #if length of reader is 0, return empty dict
        print("reader",reader)
        if len(reader)==0:
            print("LENGTH OF READER IS 0 FOR ACETOSCAN OUTPUT",acetoscan_output)
            return {}
        #skip the first line
        #next(reader)
        #for each line, get the species name and the bit score
        for row in reader:
            #if len of row is less than 3, skip
            print("row",row,len(row))
            if len(row)<3:
                continue
            query=row[0]
            subject=row[7]
            bit_score=float(row[-2])
            print("Found acetoscan match:",subject,bit_score)
            unique_queries[query].append((subject,bit_score))
    #for each query, sort the subjects by bit score
    for query,subject_list in unique_queries.items():
        subject_list.sort(key=lambda x:x[1],reverse=True)
        unique_queries[query]=subject_list[0][0]
    #the description is simply species_name_FTFHS
    for query,subject in unique_queries.items():
        unique_queries[query]=subject+"_FTFHS"
        print("Best acetoscan match:",subject)
    print("unique queries",str(unique_queries))
    return unique_queries

    
    pass
def is_acetogen(acetoscan_output):
    with open(acetoscan_output) as f:
        reader=csv.reader(f,delimiter=',')
        readerlist=process_acetoscan_output(acetoscan_output)
        if len(readerlist)>0:
            print("READER LIST IS",readerlist)
            return True,readerlist
        else:
            return False,{}
hyd_data = {
    "[FeFe] Group A1": "H2-evolution",
    "[FeFe] Group A2": "H2-uptake ?",
    "[FeFe] Group A3": "Electron-bifurcation",
    "[FeFe] Group A4": "H2-evolution or electron-bifurcation",
    "[FeFe] Group B" : "H2-evolution?",
    "[FeFe] Group C1": "H2-sensing?",
    "[FeFe] Group C2": "H2-sensing?",
    "[FeFe] Group C3": "H2-sensing?",
    "[Fe] Group Hmd": "Bidirectional",
    "[NiFe] Group 1a": "H2-uptake",
    "[NiFeSe] Group 1a": "H2-uptake",
    "[NiFe] Group 1b": "H2-uptake",
    "[NiFe] Group 1c": "H2-uptake",
    "[NiFe] Group 1d": "H2-uptake",
    "[NiFe] Group 1e": "Bidirectional",
    "[NiFe] Group 1f": "H2-uptake ?",
    "[NiFe] Group 1g": "H2-uptake",
    "[NiFe] Group 1h": "H2-uptake",
    "[NiFe] Group 1i": "H2-uptake",
    "[NiFe] Group 1j": "H2-uptake",
    "[NiFe] Group 1k": "H2-uptake",
    "[NiFe] Group 2a": "H2-uptake",
    "[NiFe] Group 2b": "H2-sensing",
    "[NiFe] Group 2c": "H2-sensing?",
    "[NiFe] Group 2d": "H2-uptake?",
    "[NiFe] Group 2e": "H2-uptake?",
    "[NiFe] Group 3a": "Bidirectional",
    "[NiFe] Group 3b": "Bidirectional",
    "[NiFe] Group 3c": "Electron-bifurcation",
    "[NiFe] Group 3d": "Bidirectional",
    "[NiFe] Group 4a": "H2-evolution",
    "[NiFe] Group 4b": "H2-evolution",
    "[NiFe] Group 4c": "H2-evolution",
    "[NiFe] Group 4d": "H2-evolution",
    "[NiFe] Group 4e": "Bidirectional",
    "[NiFe] Group 4f": "H2-evolution?",
    "[NiFe] Group 4g": "H2-evolution?",
    "[NiFe] Group 4h": "H2-uptake",
    "[NiFe] Group 4i": "H2-uptake",
}
def is_hydrogenase2(hyDB_output):
    processed_output=process_diamond_output(hyDB_output,bit_score_index=-4,inject_desc=False)
    unique_queries={}
    at_least_one_uptake = False
    at_least_one_bifurcating_or_bidirectional = False
    for query in processed_output:
        hyd_type=processed_output[query]
        hyd_types = hyd_type.split(' - ')
        # get the last one part and remove the fina
        hyd_type = hyd_types[-1]
        # find the hyd type in the dict
        if hyd_type in hyd_data:
            value = hyd_data[hyd_type]
            if "uptake" in value:
                at_least_one_uptake = True
                unique_queries[query] = processed_output[query]
            if "bifurcation" in value or "bidirectional" in value:
                at_least_one_bifurcating_or_bidirectional = True
                unique_queries[query] = processed_output[query]
            else:
                unique_queries[query] = processed_output[query]
    return at_least_one_uptake,at_least_one_bifurcating_or_bidirectional,unique_queries


def is_hydrogenase(hyDB_output):
    processed_output=process_diamond_output(hyDB_output,bit_score_index=-4)
    with open(hyDB_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=list(reader)
        at_least_one_uptake=False
        at_least_one_bifurcating_or_bidirectional=False
        for row in readerlist:
            #example of a row: reverse translation of NiFe-WP_011973026.1 - Methanococcus aeolicus	NiFe-WP_011973026.1 - Methanococcus aeolicus - [NiFe] Group 4i	100	374	1	1122	1	374
            #get the 2nd col
            query=row[0]
            if len(row)>1:
                hyd_type=row[1]
                #split by -
                hyd_type=hyd_type.split('-')
                #get the 2nd part
                hyd_type=hyd_type[1]
                #find the hyd type in the dict
                if hyd_type in hyd_data:
                    value=hyd_data[hyd_type]
                    #if value contains uptake, then return true
                    if 'uptake' in value:
                        at_least_one_uptake=True
                    #if value contains bifurcation or bidirectional, then return true
                    if 'bifurcation' in value or 'Bidirectional' in value:
                        at_least_one_bifurcating_or_bidirectional=True
        if at_least_one_uptake:
            return True
        elif at_least_one_bifurcating_or_bidirectional:
            return "Maybe"
        else:
            return False

def is_hydrogen_consumer(query,dsr_diamond_output,mcr_diamond_output,acetoscan_output,hyDB_output):
    #check if it is a sulphate reducer
    sulphate_reducer=False
    methanogen=False
    acetogen=False
    hydrogenase=False
    is_sulphate, unique_dictd = is_sulphate_reducer(dsr_diamond_output)
    is_methano, unique_dictm = is_methanogen(mcr_diamond_output)
    is_aceto, unique_dicta = is_acetogen(acetoscan_output)
    is_consume, is_bifurcating, unique_dicth = is_hydrogenase2(hyDB_output)
    is_hydrogen_consumer=False
    if is_sulphate and (is_consume or is_bifurcating):
        print("is sulphate reducer")
        is_hydrogen_consumer=True
        #
    #check if it is a methanogen
    elif is_methano and (is_consume or is_bifurcating):
        print("is methanogen")
        is_hydrogen_consumer=True
    #check if it is an acetogen
    elif is_aceto and (is_consume or is_bifurcating):
        print("is acetogen")
        is_hydrogen_consumer=True
    #check if it is a hydrogenase
    elif is_consume:
        print("contains a consumption hydrogenase")
        is_hydrogen_consumer=True
    front_row=[query,is_hydrogen_consumer,is_sulphate,is_methano,is_aceto,is_consume,is_bifurcating]
    #make the evidence rows
    evidence_rows=make_evidence_rows(unique_dictd,unique_dictm,unique_dicta,unique_dicth)
    #pad the evidence rows with empty strings from the front
    evidence_rows=[[""]*len(front_row)+row for row in evidence_rows]
    #pad the front rows with empty strings from the back
    #if no evidence rows, then add a empty row
    if len(evidence_rows)==0:
        evidence_rows=[[""]*len(front_row)]
    front_row=front_row+[""]*(len(evidence_rows[0]))
    #add the front row to the evidence rows
    all_rows=[front_row]+evidence_rows
    return all_rows


def make_evidence_rows(unique_dictd,unique_dictm,unique_dicta,unique_dicth):
    #get the dict with the max rows
    lengths=[len(unique_dictd),len(unique_dictm),len(unique_dicta),len(unique_dicth)]
    lisd=list(unique_dictd.values())
    lism=list(unique_dictm.values())
    lisa=list(unique_dicta.values())
    print("LISA:",lisa)
    lish=list(unique_dicth.values())
    max_length=max(lengths)
    rows=[]
    for r in range(max_length):
        row=[]
        if r<len(lisd):
            row.append(lisd[r])
        else:
            row.append('')
        if r<len(lism):
            row.append(lism[r])
        else:
            row.append('')
        if r<len(lisa):
            row.append(lisa[r])
        else:
            row.append('')
        if r<len(lish):
            row.append(lish[r])
        else:
            row.append('')
        rows.append(row)
    return rows

def run_diamond(input_file,diamond_exe_path,database_path,type,tmp_dir,min_percent_identity=90,min_scov=0.9,database_fasta=None):
    diamond_output=os.path.join(os.path.dirname(input_file),type+"_diamond_output.tsv")
    command=diamond_exe_path+' blastp -d '+database_path+' -q '+input_file+' -o '+diamond_output+' --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen'
    os.system(command)
    #move the output file to the tmp dir
    #shutil.copy(diamond_output,os.path.join(tmp_dir,os.path.basename(diamond_output)))
    #new_damond_output=os.path.join(tmp_dir,os.path.basename(diamond_output))
    #open the diamond output and sort the rows by bitscore
    #open the database fasta file
    if database_fasta is None:
        raise Exception("Database fasta file is None")
    #get the lengths of the sequences in the database
    lengths={}
    with open(database_fasta,'r') as f:
        records=SeqIO.parse(f,'fasta')
        for record in records:
            lengths[record.id]=len(record.seq)
    #filter the rows by min percent identity and min scov
    with open(diamond_output,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        rows=[row for row in reader]
    #filter the rows by min percent identity and min scov
    filtered_rows=[]
    for row in rows:
        pid=False
        scov=False
        if float(row[2])>=min_percent_identity:
            pid=True
        sstart=int(row[8])
        send=int(row[9])
        align=(send-sstart+1)
        scov=align/lengths[row[1]]
        if scov>=min_scov:
            scov=True
        if pid and scov:
            filtered_rows.append(row)
    #sort the rows by bitscore
    sorted_rows=sorted(filtered_rows,key=lambda x: float(x[11]),reverse=True)
    #write the sorted rows to a new file
    diamond_output_dir=os.path.dirname(diamond_output)
    new_diamond_output=os.path.join(diamond_output_dir,type+"_diamond_output_filtered.tsv")
    with open(new_diamond_output,'w') as f:
        writer=csv.writer(f,delimiter='\t')
        writer.writerows(sorted_rows)


    return new_diamond_output
def join_and_add_labels_to_records(input_dir,output_fasta):
    #join the records in the input dir into one fasta file, add an identifier to the start of each record representing the file it came from
    #get the files in the input dir
    sep="_file_"
    files=os.listdir(input_dir)
    #make a dict to store the records
    records={}
    #iterate over the files
    for file in files:
        #if file not ends with .faa, or .fasta or .fna, skip
        if not file.endswith(".faa") and not file.endswith(".fasta") and not file.endswith(".fna"):
            continue
        #open the file
        with open(os.path.join(input_dir,file),'r') as f:
            #get the records
            records[file]=list(SeqIO.parse(f,'fasta'))
    #make a new dict to store the joined records
    joined_records={}
    #iterate over the records
    for file in records:
        #iterate over the records in the file
        for record in records[file]:
            #add the file identifier to the start of the record id
            record.id=file+sep+record.id
            #add the record to the joined records dict
            joined_records[record.id]=record
    #write the joined records to the output file
    with open(output_fasta,'w') as f:
        SeqIO.write(joined_records.values(),f,'fasta')
    return output_fasta
def batch_diamond(joined_records, diamond_exe_path, database_path, type, tmp_dir, min_percent_identity=90, min_scov=0.9, database_fasta=None,resume=True):
    #join the records in the input dir into one fasta file, add an identifier to the start of each record representing the file it came from
    #joined_records=join_and_add_labels_to_records(input_dir,os.path.join(tmp_dir,type+"_joined_records.fasta"))
    #run diamond on the joined records
    output_file=os.path.join(tmp_dir,type+"_diamond_output_filtered.tsv")
    #if resume and output file exists, skip
    if os.path.exists(output_file) and resume:
        print("Diamond output file exists",type," skipping")
        return output_file
    command=diamond_exe_path+' blastp -d '+database_path+' -q ' +joined_records+' -o '+os.path.join(tmp_dir,type+"_diamond_output.tsv")+' --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen --fast'
    os.system(command)
    # open the database fasta file
    if database_fasta is None:
        raise Exception("Database fasta file is None")
    # get the lengths of the sequences in the database
    lengths = {}
    with open(database_fasta, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            lengths[record.id] = len(record.seq)
    # filter the rows by min percent identity and min scov
    with open(os.path.join(tmp_dir,type+"_diamond_output.tsv"), 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = [row for row in reader]
    # filter the rows by min percent identity and min scov
    filtered_rows = []
    for row in rows:
        pid = False
        scov = False
        if float(row[2]) >= min_percent_identity:
            pid = True
        sstart = int(row[8])
        send = int(row[9])
        align = (send - sstart + 1)
        scov = align / lengths[row[1]]
        if scov >= min_scov:
            scov = True
        if pid and scov:
            filtered_rows.append(row)
   #write the filtered rows to a new file
    with open(os.path.join(tmp_dir,type+"_diamond_output_filtered.tsv"), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(filtered_rows)
    return os.path.join(tmp_dir,type+"_diamond_output_filtered.tsv")




def batch_acetoscan(joined_records,tmp_dir,resume=True,run_acetoscan=True):
    #make the acetoscan inout dir in the tmp dir
    acetoscan_input_dir=os.path.join(tmp_dir,"acetoscan_input")
    if not os.path.exists(acetoscan_input_dir):
        os.mkdir(acetoscan_input_dir)
    #run acetoscan on the joined records
    output_file = os.path.join(acetoscan_input_dir, 'aceto_out.csv')
    new_output_file = os.path.join(tmp_dir, os.path.basename(output_file))
    #if in resume mode and output file exists, return
    if resume and os.path.exists(new_output_file):
        print("Output file exists, skipping acetoscan")
        return new_output_file
    '''command:
    sudo docker run --rm -v {input_dir}:/acetoscan/input_dir --entrypoint acetotax -it abhijeetsingh1704/acetoscan -o /acetoscan/input_dir/aceto_out  -i /acetoscan/input_dir/{genome_file}
    '''
    #move the input file to the acetoscan input dir
    shutil.copy(joined_records,os.path.join(acetoscan_input_dir,os.path.basename(joined_records)))
    #run acetoscan
    command='sudo docker run --rm -v '+acetoscan_input_dir+':/acetoscan/input_dir --entrypoint acetotax -it abhijeetsingh1704/acetoscan -o /acetoscan/input_dir/aceto_out  -i /acetoscan/input_dir/'+os.path.basename(joined_records)
    #if run the acetoscan command
    if run_acetoscan:
        os.system(command)
    #get the output file

    #if output file not exists, then make a blank file
    if not os.path.exists(output_file):
        open(output_file,'w').close()
    #remove the input file from the acetoscan input dir if it exists
    if os.path.exists(os.path.join(acetoscan_input_dir,os.path.basename(joined_records))):
        os.remove(os.path.join(acetoscan_input_dir,os.path.basename(joined_records)))
    #move the output file to the tmp dir
    shutil.copy(output_file,os.path.join(tmp_dir,os.path.basename(output_file)))
    new_output_file=os.path.join(tmp_dir,os.path.basename(output_file))
    return new_output_file


def batch_run_hydrogenase_script(ijoined_records,hydrogenase_script_path,tmp_dir,resume=True):
    #joined_records = join_and_add_labels_to_records(input_dir, os.path.join(tmp_dir, "hyd_joined_records.fasta"))
    '''command:
    bash {hydrogenase_script_path} {hyinput} {hyoutput}
    '''
    #copy the joined_records to tmp_dir/hyDinput
    #make the hyinput dir, if hyinput dir exists, and resume is false, then remove the hyinput dir

    output_filepath = os.path.join(tmp_dir, "hyDB_output")

    # if in resume mode, and the output filepath exists, then return the output file
    if os.path.exists(os.path.join(output_filepath, "Filtered", "all_hydrogenase.txt")) and resume:
        print("Output file exists, skipping hydrogenase script")
        return os.path.join(output_filepath, "Filtered", "all_hydrogenase.txt")

    if os.path.exists(os.path.join(tmp_dir,'hyinput')) and not resume:
        shutil.rmtree(os.path.join(tmp_dir,'hyinput'))
    os.makedirs(os.path.join(tmp_dir,'hyinput'),exist_ok=True)
    hyinput_dir=os.path.join(tmp_dir,'hyinput')
    shutil.copy(ijoined_records,os.path.join(tmp_dir,'hyinput'))
    #rename the file to .fasta
    #remove the extension
    new_name=os.path.splitext(os.path.basename(ijoined_records))[0]
    new_name= new_name+'.fasta'
    os.rename(os.path.join(hyinput_dir,os.path.basename(ijoined_records)),os.path.join(hyinput_dir,new_name))

    #mkdir the output dir if it does not exist, else if not resume, remove the dir and make it
    if os.path.exists(output_filepath) and not resume:
        shutil.rmtree(output_filepath)
    os.makedirs(output_filepath,exist_ok=True)
    print("MAKE DIR:",output_filepath)


    command='bash '+hydrogenase_script_path+' '+hyinput_dir+' '+output_filepath
    print("COMMAND:",command)
    os.system(command)
    #get the output file
    output_file=os.path.join(output_filepath,"Filtered","all_hydrogenase.txt")
    #move the output file to the tmp dir
    shutil.copy(output_file,os.path.join(tmp_dir,os.path.basename(output_file)))
    #return the output file
    output_file=os.path.join(tmp_dir,"all_hydrogenase.txt")
    return output_file
    #remove the input file from the hyDB input dir
    #os.remove(os.path.join(hyDb_input_dir,os.path.basename(input_file))+'.fasta')
    #os.remove(os.path.join(hyDb_input_dir,os.path.basename(input_file)))
    #return output_file
def run_prodigal(input_file,prodigal_exe_path,out_dir=None,resume=True):
    '''command:
    {prodigal_exe_path} -i {input_file} -a {output_file} -o {output_file2}
    '''
    basename=os.path.basename(input_file)
    #remove the extension
    basename=os.path.splitext(basename)[0]
    if out_dir is None:
        output_file=os.path.join(os.path.dirname(input_file),basename+'_prodigal.faa')
    else:
        output_file=os.path.join(out_dir,basename+'_prodigal.faa')
    #if prodigal output file exists, then return the prodigal output file
    if os.path.exists(output_file) and resume:
        return output_file
    command=prodigal_exe_path+' -i '+input_file+' -a '+output_file +' -o '+output_file+'.txt'
    os.system(command)
    return output_file



def group_by_genome(input_dir, dsr_output, mcr_output, acetoscan_output, hyDB_output,tmp_dir,extensions):
    #if tmp dir not exists, make it, else, remove all files in it
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir,exist_ok=True)
    #else:
        #shutil.rmtree(tmp_dir)
        #os.makedirs(tmp_dir,exist_ok=True)
    #for each genome, get the dsr, mcr, acetoscan and hyDB results
    groups={}
    #open the dsr output file
    dsr_rows=[]
    with open(dsr_output) as f:
        reader=csv.reader(f,delimiter='\t')
        #skip the header
        #next(reader)
        for row in reader:
            dsr_rows.append(row)
    #open the mcr output file
    mcr_rows=[]
    with open(mcr_output) as f:
        reader=csv.reader(f,delimiter='\t')
        #if
        #skip the header
        #next(reader)
        for row in reader:
            mcr_rows.append(row)
    #open the acetoscan output file
    acetoscan_rows=[]
    print("acetoscan output:",acetoscan_output)
    header=None
    with open(acetoscan_output) as f:
        reader=csv.reader(f,delimiter=',')
        #skip the header, if header is not None

        #next(reader)
        for row in reader:
            if len(row)==0:
                continue
            if header is None:
                header=row
                continue
            acetoscan_rows.append(row)
    #open the hyDB output file
    hyDB_rows=[]
    with open(hyDB_output) as f:
        reader=csv.reader(f,delimiter='\t')
        #skip the header
        #next(reader)
        for row in reader:
            hyDB_rows.append(row)
    #get all the files in the input dir
    files=os.listdir(input_dir)
    #for each file, get the results
    for file in files:
        #if file is not fna, faa fasta, skip
        #if not file.endswith(".fna") and not file.endswith(".faa") and not file.endswith(".fasta"):
        if not any([file.endswith(ext) for ext in extensions]):
            continue
        file_no_ext=os.path.splitext(file)[0]
        #for all dsr rows, find the row[0] that starts with the file name
        dsr_res=[]
        for row in dsr_rows:
            if row[0].startswith(file_no_ext):
                dsr_res.append(row)
        #for all mcr rows, find the row[0] that starts with the file name
        mcr_res=[]
        for row in mcr_rows:
            if row[0].startswith(file_no_ext):
                mcr_res.append(row)
        #for all acetoscan rows, find the row[0] that starts with the file name
        acetoscan_res=[]
        for row in acetoscan_rows:
            if row[0].startswith(file_no_ext):
                acetoscan_res.append(row)
        #for all hyDB rows, find the row[0] that starts with the file name
        hyDB_res=[]
        for row in hyDB_rows:
            if row[0].startswith(file_no_ext):
                hyDB_res.append(row)
        #add the results to the groups
        groups[file]=(dsr_res,mcr_res,acetoscan_res,hyDB_res)
    #for each group, write each result to a separate file
    for i,group in enumerate(groups):
        print("GENOME:",i,"/",len(groups)," ",group)
        dsr_rows=groups[group][0]
        #write the dsr results to a file
        with open(os.path.join(tmp_dir,group+"_dsr_matches.tsv"),'w') as f:
            writer=csv.writer(f,delimiter='\t')
            writer.writerows(dsr_rows)
        mcr_rows=groups[group][1]
        #write the mcr results to a file
        with open(os.path.join(tmp_dir,group+"_mcr_matches.tsv"),'w') as f:
            writer=csv.writer(f,delimiter='\t')
            writer.writerows(mcr_rows)
        acetoscan_rows=groups[group][2]
        #write the acetoscan results to a file
        with open(os.path.join(tmp_dir,group+"_acetoscan_matches.tsv"),'w') as f:
            writer=csv.writer(f,delimiter=',')
            writer.writerows(acetoscan_rows)
        hyDB_rows=groups[group][3]
        #write the hyDB results to a file
        with open(os.path.join(tmp_dir,group+"_hyDB_matches.tsv"),'w') as f:
            writer=csv.writer(f,delimiter='\t')
            writer.writerows(hyDB_rows)
        #run the is_hydrogen_consumer function on the group
        res=is_hydrogen_consumer(group,os.path.join(tmp_dir,group+"_dsr_matches.tsv"),os.path.join(tmp_dir,group+"_mcr_matches.tsv"),
                                 os.path.join(tmp_dir,group+"_acetoscan_matches.tsv"),os.path.join(tmp_dir,group+"_hyDB_matches.tsv"))
        yield res

    #return groups
def batch_run2(input_dir,output_file,tmp_dir,args):
    print(f"batch_run2 (input_dir={input_dir},output_file={output_file},tmp_dir={tmp_dir},args={args})")
    #if tmp dir is None, then tmp dir is input_dir/tmp
    if tmp_dir is None:
        tmp_dir=os.path.join(input_dir,"tmp")
    #if tmp dir not exists, make it, else, if not in resume mode, delete it and make it again
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    elif not args.resume:
        shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir)
    files=os.listdir(input_dir)
    #create a list to store the results
    results= []
    #for each file, run the pipeline
    prodigal_tmp_dir=os.path.join(tmp_dir,"prodigal")
    #if the tmp dir does not exist, create it. if it exists, delete it and create it again, if not in resume mode
    if os.path.exists(prodigal_tmp_dir):
        if not args.resume:
            shutil.rmtree(prodigal_tmp_dir)
    os.makedirs(prodigal_tmp_dir,exist_ok=True)
    #run the prodigal
    for file in files:
        #get only files with .faa, fasta or .fna
        #if file.endswith(".faa") or file.endswith(".fasta") or file.endswith(".fna") or file.endswith(".fa"):
        if any([file.endswith(ext) for ext in args.file_extensions]):
            #get the full path
            file=os.path.join(input_dir,file)
            #run the prodigal
            #if args.prodigal_path is None, set it to the default
            if args.prodigal_exe_path is None:
                prodigal_path=default_prodigal_path
            else:
                prodigal_path=args.prodigal_exe_path


            prodigal_output=run_prodigal(file,prodigal_path,prodigal_tmp_dir,resume=args.resume)
    #join the fasta files into one
    dna_joined=join_and_add_labels_to_records(input_dir,os.path.join(tmp_dir,"dna_joined.fna"))
    #join the records in the prodigal output
    prodigal_joined=join_and_add_labels_to_records(prodigal_tmp_dir,os.path.join(tmp_dir,"prodigal_joined.faa"))
    #run the batch dsr
    #def batch_diamond(joined_records, diamond_exe_path, database_path, type, tmp_dir, min_percent_identity=90, min_scov=0.9, database_fasta=None):
    #if args.diamond_exe_path is None, set it to the default
    if args.diamond_exe_path is None:
        diamond_exe_path=default_diamond_exe_path
    else:
        diamond_exe_path=args.diamond_exe_path
    #if args.dsr_database_path is None, set it to the default
    if args.dsr_database_path is None:
        dsr_database_path=default_dsr_database_path
    else:
        dsr_database_path=args.dsr_database_path
    dsr_output=batch_diamond(prodigal_joined,diamond_exe_path,dsr_database_path,"dsr",tmp_dir,min_percent_identity=90,min_scov=0.9,
                             database_fasta=dsr_database_fasta_file,resume=args.resume)
    #run the batch mcr
    #if args.mcr_database_path is None, set it to the default
    if args.mcr_database_path is None:
        mcr_database_path=default_mcr_database_path
    else:
        mcr_database_path=args.mcr_database_path
    mcr_output=batch_diamond(prodigal_joined,default_diamond_exe_path,mcr_database_path,"mcr",tmp_dir,min_percent_identity=90,min_scov=0.9,
                             database_fasta=mcr_database_fasta_file,resume=args.resume)
    #run the batch hyDB
    #batch_run_hydrogenase_script(ijoined_records,hydrogenase_script_path,tmp_dir)
    #if args.hydrogenase_script_path is None, set it to the default
    if args.hydrogenase_script_path is None:
        hydrogenase_script_path=default_hydrogenase_script_path
    else:
        hydrogenase_script_path=args.hydrogenase_script_path
    hyDB_output=batch_run_hydrogenase_script(prodigal_joined,hydrogenase_script_path,tmp_dir,resume=args.resume)
    #run the batch acetoscan
    #batch_acetoscan(joined_records,tmp_dir)
    #
    acetoscan_output=batch_acetoscan(dna_joined,tmp_dir, run_acetoscan=args.acetoscan,resume=args.resume)
    results=[]
    #group_by_genome(input_dir, dsr_output, mcr_output, acetoscan_output, hyDB_output,tmp_dir)
    for res in group_by_genome(input_dir, dsr_output, mcr_output, acetoscan_output, hyDB_output,tmp_dir,args.file_extensions):
        results.extend(res)
    # check if the output file already exists, if exists, append to it
    header = False
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            str1 = f.read()
            # if the first row is a header, then set header=True
            if str1.startswith("filename"):
                header = True
    # write the results to the output file
    with open(output_file, 'a') as f:
        writer = csv.writer(f, delimiter='\t')
        # if header is not there, write the header
        if not header:
            # the headers are the filename, hydrogen_consumer, sulphate_reducer, methanogen, acetogen, hydrogenase
            writer.writerow(["filename", "hydrogen_consumer", "sulphate_reducer", "methanogen", "acetogen",
                             "consuming hydrogenase", "bifurcating_or_bidirectional_hydrogenase",
                             "dsr_matches", "mcr_matches", "acetoscan_matches", "hyDB_matches"])
        for row in results:
            writer.writerow(row)
    # open the file again and replace the uids with descriptions
    with open(output_file, 'r') as f:
        str1 = f.read()
    # replace the uids with descriptions
    Batch_Fetcher.fetch()
    repl = Batch_Fetcher.replace(str1)
    # write the replaced string to the file
    with open(output_file, 'w') as f:
        f.write(repl)
    return results


def generate_dummy_CAT_data(filenames, out_file):

        scores = [0.99, 0.85, 0.72, 0.65]


        taxonomy_levels = ['root', 'cellular organisms', 'Bacteria', 'clade', 'phylum',
                           'class', 'order', 'family', 'genus', 'species']
        names = {
            'root': 'root',
            'cellular organisms': 'cellular organisms',
            'Bacteria': 'Bacteria',
            'clade': ['Terrabacteria group', 'FCB group'],
            'phylum': ['Firmicutes', 'Bacteroidetes/Chlorobi group'],
            'class': ['Clostridia', 'Bacteroidia'],
            'order': ['Clostridiales', 'Bacteroidales'],
            'family': ['Lachnospiraceae', 'Bacteroidaceae'],
            'genus': ['Lachnospiraceae bacterium', 'Bacteroides'],
            'species': ['unclassified Lachnospiraceae', 'Bacteroides genus']
        }

        # Generate the random data
        cat_data = []
        for filename in filenames:
            #skip files not containing .fna, .fa, .fas, .fasta or .faa
            if not filename.endswith(".fna") and not filename.endswith(".fa") and not filename.endswith(".fas") and not filename.endswith(".fasta") and not filename.endswith(".faa"):
                continue
            taxid = 'taxid assigned'
            orfs = f"based on {random.randint(1000, 6000)}/{random.randint(6000, 10000)} ORFs"
            lineage = ['1', '131567', '2'] + [str(random.randint(100000, 999999)) for _ in range(7)]
            lineage_scores = [f"{random.uniform(0.3, 1):.2f}" for _ in taxonomy_levels]
            full_lineage_names = [
                f"{random.choice(names[taxonomy_levels[i]])} ({taxonomy_levels[i]}): {lineage_scores[i]}"
                #f"{taxonomy_levels[i]} ({taxonomy_levels[i]}): {lineage_scores[i]}"
                for i in range(len(taxonomy_levels))
            ]

            cat_line = [filename, taxid, orfs, ';'.join(lineage), ';'.join(lineage_scores)] + full_lineage_names
            cat_data.append('\t'.join(cat_line))


        #return data
        # save the data to the out file
        with open(out_file, 'w') as f:
            #writer=csv.writer(f,delimiter='\t')
            for line in cat_data:
                f.write(line + '\n')
        return cat_data













def run_CAT(input_dir, tmp_folder, resume_mode=True, dummy_data=False):
    print(f"run_CAT (input_dir={input_dir},tmp_folder={tmp_folder},resume_mode={resume_mode},dummy_data={dummy_data})")
    # create a CAT folder in the tmp folder, if it does not exist
    print("CAT folder:",os.path.join(tmp_folder, "CAT"), os.path.exists(os.path.join(tmp_folder, "CAT")))
    if not os.path.exists(os.path.join(tmp_folder, "CAT")):
        print("CAT file not exists, creating")
        os.makedirs(os.path.join(tmp_folder, "CAT"))
    else:
        # if it exists, and resume mode if false, delete it and create it again
        if not resume_mode:

            shutil.rmtree(os.path.join(tmp_folder, "CAT"))
            os.makedirs(os.path.join(tmp_folder, "CAT"))

    cat_folder = os.path.join(tmp_folder, "CAT")
    current_wd = os.getcwd()
    # make a input folder in the CAT folder
    if not os.path.exists(os.path.join(cat_folder, "input")):
        os.makedirs(os.path.join(cat_folder, "input"))

    # copy the input files to the input folder
    for file in os.listdir(input_dir):
        #only copy the fna, fa, fas, fasta and faa files
        ext=os.path.splitext(file)[1]
        if ext in [".fna",".fa",".fas",".fasta",".faa"]:
            shutil.copy(os.path.join(input_dir, file), os.path.join(cat_folder, "input"))
    #deduplicate the headers
    deduplicate_names(os.path.join(cat_folder, "input"))
    # if dummy data is False
    if not dummy_data:
        # change dir to the CAT folder

        os.chdir(cat_folder)
        # run the CAT command
        # CAT bins -b {bin folder} -d {database folder} -t {taxonomy folder}
        classification_file = os.path.join(cat_folder, "out.BAT.bin2classification.txt")
        #if classification file exists, then dont run CAT
        if not os.path.exists(classification_file):
            command = f"CAT bins -b {os.path.join(cat_folder, 'input')} -d {CAT_database_path} -t {CAT_taxonomy_path} --force"
            print("COMMAND:",command)
            os.system(command)
        else:
            print("CAT output file exists, skipping CAT bins")
        # add the names: CAT add_names -i {ORF2LCA / classification file} -o {output file} -t {taxonomy folder}
        #if output file exists, then dont run CAT
        output_file = os.path.join(cat_folder, "out.BAT.bin2classification.named.txt")
        if not os.path.exists(output_file):

            command = f"CAT add_names -i {classification_file} -o {output_file} -t {CAT_taxonomy_path}"
            print("COMMAND:",command)
            os.system(command)
        else:
            print("Named Output file exists, skipping CAT add names")
        # cd back to the original wd
        os.chdir(current_wd)
    # if dummy data is True
    else:
        # create a dummy classification file
        dummy_classification_file = os.path.join(cat_folder, "out.CAT.contig2classification.txt")
        dummy_classification_data = generate_dummy_CAT_data(os.listdir(input_dir), dummy_classification_file)
        output_file= os.path.join(cat_folder, "out.CAT.contig2classification.named.txt")
    # return the classification file
    return output_file


def create_feature_table (tmp_filename,hydrogen_classes_activity_file,hydrogen_consumer_pipeline_output_folder,out_file):
    hydrogen_consumer_pipeline_output_folder=os.path.join(hydrogen_consumer_pipeline_output_folder,"batches")
    #extract the basename of tmp_filename
    tmp_filename=os.path.basename(tmp_filename)
    #COMBINE THE CAT OUTPUT
    # combine the out.CAT.contig2classification.named.txt from the /tmp/CAT from each folder in the hydrogen_consumer_pipeline_output_folder
    # get the folders in the hydrogen_consumer_pipeline_output_folder
    folders = os.listdir(hydrogen_consumer_pipeline_output_folder)
    # get the folders that are dirs
    folders = [folder for folder in folders if
               os.path.isdir(os.path.join(hydrogen_consumer_pipeline_output_folder, folder))]
    combined_classification_file = os.path.join(hydrogen_consumer_pipeline_output_folder,
                                                'combined_classification_file.txt')
    # combine the classification files
    with open(combined_classification_file, 'w') as f:

        for folder in folders:
            header=None
            # get the classification file
            classification_file = os.path.join(hydrogen_consumer_pipeline_output_folder, folder, tmp_filename,"CAT",
            'out.BAT.bin2classification.named.txt')
            # open the file
            with open(classification_file, 'r') as f1:
                #skip the header, if header  not  None
                if header is not None:
                    next(f1)
                else:
                    header=next(f1)
                    #write the header to the combined file
                    f.write(header)
                #read the rest of the file
                str1 = f1.read()
                #write the file to the combined file
                f.write(str1)
    taxon_output=combined_classification_file




    #get the arguments
    #def main(hydrogen_consumer_pipeline_output_folder, hydrogen_classes_activity_file, taxon_output, out_file):
    out_file=CFT.main(hydrogen_consumer_pipeline_output_folder,hydrogen_classes_activity_file,taxon_output,out_file)
    #create_feature_table_script_path=CREATE_FEATURE_TABLE_SCRIPT_PATH
    #command=f'python {create_feature_table_script_path} --hydrogen_consumer_pipeline_output_folder {hydrogen_consumer_pipeline_output_folder} --hydrogen_classes_activity_file {hydrogen_classes_activity_file} --taxon_output {taxon_output} --out_file {out_file}'
    #run the command
    #print("COMMAND:",command)
    #os.system(command)
    #return the output file
    #fix the feature table
    #main(input_file,out_file)
    out_file=FFT.main(out_file,out_file)
    #fix_feature_table_script_path=FIX_FEATURE_TABLE_SCRIPT_PATH
    #command=f'python {fix_feature_table_script_path} --input_file {out_file} --out_file {out_file}'
    #print("COMMAND:",command)
    #os.system(command)
    print("DONE! The feature table is in:",out_file)
    return out_file

def batch_the_files(input_dir, batch_root_dir, batch_size=100,resume=True):
    print(f"batch_the_files(input_dir={input_dir}, batch_root_dir={batch_root_dir}, batch_size={batch_size},resume={resume})")
    #get all the files that end in .fna, .fasta , .fa or .fas
    files=os.listdir(input_dir)
    #get only the files that end in .fna, .fasta , .fa or .fas
    files=[file for file in files if file.endswith(".fna") or file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fas")]
    #sort the files
    files.sort()
    #make the root batch dir if it does not exist
    if not os.path.exists(batch_root_dir):
        os.makedirs(batch_root_dir)
    #if resume is true and the batch root dir is not empty, then return the batch root dir
    if resume and len(os.listdir(batch_root_dir))>0:
        return batch_root_dir
    for i in range(0,len(files),batch_size):
        #get the batch
        batch=files[i:i+batch_size]
        #make a folder for each batch
        batch_dir=os.path.join(batch_root_dir,"batch_"+str(i))
        if not os.path.exists(batch_dir):
            os.makedirs(batch_dir)
        #copy the files to the batch dir
        for file in batch:
            shutil.copy(os.path.join(input_dir,file),os.path.join(batch_dir,file))
    return batch_root_dir
def create_dummy_results(input_dir,output_file):
    #create a dummy df with the index as the filenames
    files=os.listdir(input_dir)
    #filter for fa fas fasta fna
    files=[file for file in files if file.endswith(".fa") or file.endswith(".fas") or file.endswith(".fasta") or file.endswith(".fna")]
    files.sort()
    #create a dummy df
    df=pd.DataFrame(index=files)
    #add the columns, the columns are filename	hydrogen_consumer	sulphate_reducer	methanogen	acetogen	consuming hydrogenase	bifurcating_or_bidirectional_hydrogenase	dsr_matches	mcr_matches	acetoscan_matches	hyDB_matches
    #make random, True/Falses for hydrogen consumer, sulphate reducer, methanogen, acetogen, consuming hydrogenase, bifurcating_or_bidirectional_hydrogenase
    df["hydrogen_consumer"]=np.random.choice([True,False],size=len(files))
    df["sulphate_reducer"]=np.random.choice([True,False],size=len(files))
    df["methanogen"]=np.random.choice([True,False],size=len(files))
    df["acetogen"]=np.random.choice([True,False],size=len(files))
    df["consuming hydrogenase"]=np.random.choice([True,False],size=len(files))
    df["bifurcating_or_bidirectional_hydrogenase"]=np.random.choice([True,False],size=len(files))
    #make random names for dsr_matches, mcr_matches, acetoscan_matches, hyDB_matches
    dsr_choices='''1. WP_745832199.2 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
2. WP_485720431.3 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
3. WP_994885126.1 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
4. WP_853217504.4 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
5. WP_326501879.5 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
6. WP_213776598.7 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
7. WP_140988632.8 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
8. WP_667090285.6 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
9. WP_571930482.9 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]
10. WP_422377150.0 dissimilatory-type sulfite reductase subunit alpha [uncultured Desulfovibrio sp.]'''
    dsr_choices=dsr_choices.split("\n")
    #remove the numbers
    #dsr_choices=[choice.split(".")[1:] for choice in dsr_choices]
    df["dsr_matches"]=np.random.choice(dsr_choices,size=len(files))
    mcr_choices='''1. WP_201994832.1 methyl-coenzyme M reductase subunit alpha [uncultured Methanogen sp.]
2. WP_252073189.2 methyl-coenzyme M reductase subunit beta [uncultured Methanosarcina sp.]
3. WP_325077650.3 methyl-coenzyme M reductase subunit gamma [uncultured Methanobacterium sp.]
4. WP_489002871.4 methyl-coenzyme M reductase subunit alpha [uncultured Methanococcus sp.]
5. WP_568731645.5 methyl-coenzyme M reductase subunit beta [uncultured Methanopyrus sp.]
6. WP_611372090.6 methyl-coenzyme M reductase subunit gamma [uncultured Methanoculleus sp.]
7. WP_734005812.7 methyl-coenzyme M reductase subunit alpha [uncultured Methanospirillum sp.]
8. WP_850192476.8 methyl-coenzyme M reductase subunit beta [uncultured Methanomicrobium sp.]
9. WP_967345028.9 methyl-coenzyme M reductase subunit gamma [uncultured Methanocorpusculum sp.]
10. WP_1033827590.0 methyl-coenzyme M reductase subunit alpha [uncultured Methanothrix sp.]'''
    mcr_choices=mcr_choices.split("\n")
    #remove the numbers
    #mcr_choices=[choice.split(".")[1] for choice in mcr_choices]
    df["mcr_matches"]=np.random.choice(mcr_choices,size=len(files))
    acetoscan_choices='''1. PV_100392871.1 fibronectin-binding protein FTFHS [uncultured Prevotella sp.]
2. PV_210583210.2 hemagglutinin FTFHS subunit [uncultured Prevotella ruminicola]
3. PV_320664532.3 adhesin complex protein FTFHS [uncultured Prevotella bryantii]
4. PV_430755493.4 collagen-binding surface protein FTFHS [uncultured Prevotella melaninogenica]
5. PV_540846774.5 iron acquisition FTFHS receptor [uncultured Prevotella copri]
6. PV_650937865.6 carbohydrate metabolism enzyme FTFHS [uncultured Prevotella dentalis]
7. PV_761028946.7 protease inhibitor FTFHS [uncultured Prevotella intermedia]
8. PV_872119027.8 lipase FTFHS component [uncultured Prevotella loescheii]
9. PV_983210108.9 sialidase FTFHS protein [uncultured Prevotella paludivivens]
10. PV_1094391189.0 FTFHS structural protein [uncultured Prevotella oryzae]'''
    acetoscan_choices=acetoscan_choices.split("\n")
    #remove the numbers
    #acetoscan_choices=[choice.split(".")[1] for choice in acetoscan_choices]
    df["acetoscan_matches"]=np.random.choice(acetoscan_choices,size=len(files))
    hyDB_choices='''1. FeFe-WP_038211902.2 - Caldicellulosiruptor kronotskyensis - [FeFe] Group A2
2. FeFe-WP_049532285.1 - Thermotoga lettingae - [FeFe] Group B1
3. FeFe-WP_060755466.3 - Clostridium clariflavum - [FeFe] Group A4
4. FeFe-WP_071998577.4 - Herbinix luporum - [FeFe] Group B2
5. NiFe-WP_082219688.1 - Bacteroides fragilis - [NiFe] Group 3c
6. NiFe-WP_093440799.2 - Methanosarcina barkeri - [NiFe] Group 1d
7. NiFe-WP_104661910.3 - Desulfovibrio vulgaris - [NiFe] Group 2a
8. FeFe-WP_115882021.5 - Paenibacillus polymyxa - [FeFe] Group A5
9. NiFe-WP_126003132.1 - Syntrophobacter fumaroxidans - [NiFe] Group 1b
10. NiFe-WP_137124243.4 - Geobacter metallireducens - [NiFe] Group 4a'''
    hyDB_choices=hyDB_choices.split("\n")
    #remove the numbers
    #hyDB_choices=[choice.split(".")[1] for choice in hyDB_choices]
    df["hyDB_matches"]=np.random.choice(hyDB_choices,size=len(files))
    #write the df to the output file
    df.to_csv(output_file,sep="\t")
    return output_file
def run_batches(batch_root_dir,output_file,tmp_dir,args):
    print(f"run_batches(batch_root_dir={batch_root_dir},output_file={output_file},tmp_dir={tmp_dir},args={args})")
    #get all the batch dirs
    batch_dirs=[os.path.join(batch_root_dir,batch) for batch in os.listdir(batch_root_dir)]
    #only get the folders
    batch_dirs=[batch for batch in batch_dirs if os.path.isdir(batch)]
    #sort the batch dirs
    batch_dirs.sort()
    for batch_dir in batch_dirs:
        #for each batch dir, run the pipeline
        #extract the filename of the output file
        output_file=os.path.basename(output_file)
        #add the batch folder name to the output file
        #get the parent dir of the batch dir
        parent_dir=os.path.dirname(batch_dir)

        output_file=os.path.join(parent_dir,batch_dir+"results.tsv")
        #get the foldername of the tmp dir
        tmp_dir=os.path.basename(tmp_dir)
        #add the batch folder name to the tmp dir
        tmp_dir1=os.path.join(batch_dir,tmp_dir)
        batch_run2(batch_dirs[0],output_file,tmp_dir1,args)
        #generate the dummy data
        #create_dummy_results(batch_dir,output_file)
        #run the CAT
        run_CAT(batch_dir, tmp_dir1, resume_mode=True, dummy_data=False)



import concurrent.futures
import threading


def run_batches_multithreaded(batch_root_dir, output_file, tmp_dir, args, max_batches=None):
    print(f"run_batches_multithreaded(batch_root_dir={batch_root_dir},output_file={output_file},tmp_dir={tmp_dir},args={args},max_batches={max_batches})")
    # Get all the batch dirs
    batch_dirs = [os.path.join(batch_root_dir, batch) for batch in os.listdir(batch_root_dir)]
    # Only get the folders
    batch_dirs = [batch for batch in batch_dirs if os.path.isdir(batch)]
    # Sort the batch dirs
    batch_dirs.sort()

    if max_batches is None or max_batches > len(batch_dirs):
        max_batches = len(batch_dirs)

    semaphore = threading.Semaphore(max_batches)

    def process_batch(batch_dir):
        try:
            # For each batch dir, run the pipeline
            # Extract the filename of the output file
            #output_filename = os.path.basename(output_file)
            parent_dir = os.path.dirname(batch_dir)

            output_file = os.path.join(parent_dir, batch_dir + "results.tsv")
            # Add the batch folder name to the output file
            #output_filename = os.path.join(batch_dir, output_filename)
            # Get the foldername of the tmp dir
            tmp_dirname = os.path.basename(tmp_dir)
            # Add the batch folder name to the tmp dir
            tmp_dir1 = os.path.join(batch_dir, tmp_dirname)

            # Run the pipeline
            batch_run2(batch_dir, output_file, tmp_dir1, args)
            #generate the dummy data
            #create_dummy_results(batch_dir,output_file)
            #run the CAT
            run_CAT(batch_dir, tmp_dir1, resume_mode=True, dummy_data=False)
            #
        finally:
            semaphore.release()  # Release the semaphore permit even if an exception occurs

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_batches) as executor:
        futures = []
        for batch_dir in batch_dirs:
            # Acquire a semaphore permit
            semaphore.acquire()

            # Submit the task to the executor
            future = executor.submit(process_batch, batch_dir)
            futures.append(future)

        # Wait for all the tasks to complete
        concurrent.futures.wait(futures)




def main():
    #get the arguments
    print("RUNNING")
    args=argparse.ArgumentParser()
    args.add_argument("-i","--input_dir",help="the input directory containing the genomes to be analyzed",required=True)
    args.add_argument("-o","--output_file",help="the output file to save the results",required=True)
    args.add_argument("-p","--prodigal_exe_path",help="the path to the prodigal executable",default=default_prodigal_path)
    args.add_argument("-d","--diamond_exe_path",help="the path to the diamond executable",default=default_diamond_exe_path)
    args.add_argument("-hy","--hydrogenase_script_path",help="the path to the hydrogenase script",default=default_hydrogenase_script_path)
    args.add_argument("-dsr","--dsr_database_path",help="the path to the dsr database",default=default_dsr_database_path)
    args.add_argument("-mcr","--mcr_database_path",help="the path to the mcr database",default=default_mcr_database_path)
    #add a flag for max threads
    args.add_argument("-th","--max_threads",help="the maximum number of threads to use",default=1,type=int)
    #add a flag for batch size
    args.add_argument("-b","--batch_size",help="the batch size",default=100,type=int)
    args.add_argument('-t',"--tmp_dir",help="the path to the temporary directory",default="./tmp")
    #add a flag to run the acetoscan
    args.add_argument("-a","--acetoscan",help="run acetoscan",action="store_true",default=False)
    #add an argument for file extensions to be considered
    args.add_argument("-x","--file_extensions",help="the file extensions to be considered",default=[".fa",".fna",".fasta"],nargs="+")
    #add the path to the hydrogenase activity file
    args.add_argument("-hyd","--hydrogen_classes_activity_file",help="the path to the hydrogen classes activity file",default=DEFAULT_HYDROGENASE_ACTIVITY_FILE)
    #add a flag to enable resume mode, to continue in case of failure,its default is true
    args.add_argument("-r","--resume",help="resume mode, in case the previous run crash, use cached output from previous run",action="store_true",default=True)
    args=args.parse_args()
    #batch the files
    batch_root_dir=batch_the_files(args.input_dir,os.path.join(args.input_dir,"batches"),batch_size=args.batch_size,resume=args.resume)
    #run the batches
    #run the pipeline
    #batch_run2(args.input_dir,args.output_file,tmp_dir=args.tmp_dir, args=args)
    #run the batches
    #if max threads is 1, use the normal batch run
    if args.max_threads==1:
        run_batches(batch_root_dir,args.output_file,args.tmp_dir,args)
    else:
        run_batches_multithreaded(batch_root_dir,args.output_file,args.tmp_dir,args,max_batches=args.max_threads)
    #run the CAT command
    #tax_output=run_CAT(args.input_dir,args.tmp_dir,resume_mode=args.resume)
    # create the feature table
    create_feature_table(args.tmp_dir, args.hydrogen_classes_activity_file, args.input_dir
                         , args.output_file)

def tes():
    # fetch from ncbi
    acc_list = ['TCO66124.1', 'SJZ83014.1', 'EMS81424.1', 'SHJ62183.1', 'SHO53957.1', 'SFB35260.1', 'SET75422.1',
                'CRL34893.1', 'SEU02316.1', 'SHM11210.1', 'SEW18904.1', 'SHN81160.1', 'SHH66267.1', 'SEP28343.1',
                'SHJ21341.1', 'SDN39900.1', 'SJZ87696.1', 'SFM16154.1',  'SKB00378.1', 'SEH22604.1',
                'SEE18335.1', 'SKI21864.1', 'SEE18034.1', 'SEE17036.1', 'SEE17393.1', 'SES16927.1', 'SEE18845.1',
                'SEE19522.1', 'SEE13657.1', 'SKA07415.1', 'SEE10790.1', 'SEH30585.1', 'SEE06718.1', 'SEE08601.1',
                'SEQ39612.1']
    text=""" The records found were as follows:
    {{TCO66124.1}}: 1-2 of 2
    {{SJZ83014.1}}: 1-2 of 2
    {{EMS81424.1}}: 1-2 of 2
    {{SHJ62183.1}}: 1-2 of 2
    {{SHO53957.1}}: 1-2 of 2
    {{SFB35260.1}}: 1-2 of 2
    {{SET75422.1}}: 1-2 of 2
    """
    batch_fetcher=BatchFetcher("protein")
    for acc in acc_list:
        batch_fetcher.add_accession(acc)
    batch_fetcher.fetch()
    new_text=batch_fetcher.replace(text)
    print(new_text)
def tes2():
    acetoscan_out_file=r"C:\Users\abel\Desktop\scratch\blautia\aceto_out.csv"
    dsr_out_file=r"C:\Users\abel\PycharmProjects\hydrogenases\test_pipeline\blautia\dsr_diamond_output.tsv"
    mcr_out_file=r"C:\Users\abel\PycharmProjects\hydrogenases\test_pipeline\blautia\mcr_diamond_output.tsv"
    hydrogenase_out_file=r"C:\Users\abel\PycharmProjects\hydrogenases\test_pipeline\blautia\all_hydrogenase.txt"
    rows=is_hydrogen_consumer("Blautia",dsr_out_file,mcr_out_file,acetoscan_out_file,hydrogenase_out_file)
    #save the rows to an output file
    dirname=os.path.dirname(dsr_out_file)
    output_file=os.path.join(dirname,"blautia_output.csv")
    with open(output_file,'w') as f:
        writer=csv.writer(f,delimiter='\t')
        #write the headers
        writer.writerow(["filename","hydrogen_consumer","sulphate_reducer","methanogen","acetogen","consuming hydrogenase","bifurcating_or_bidirectional_hydrogenase",
                         "dsr_matches","mcr_matches","acetoscan_matches","hyDB_matches"])

        writer.writerows(rows)
    #reopen the file and replace the uids with descriptions
    with open(output_file,'r') as f:
        str1=f.read()
    #replace the uids with descriptions
    Batch_Fetcher.fetch()
    repl=Batch_Fetcher.replace(str1)
    #write the replaced string to the file
    with open(output_file,'w') as f:
        f.write(repl)
    print(rows)
def tes_diamond_batch():
    in_dir=r'/home/ec2-user/mouse_catalogues/test_batch'
    out_dsr=r'/home/ec2-user/mouse_catalogues/test_batch_out/dsr_diamond_output.tsv'
    out_mcr=r'/home/ec2-user/mouse_catalogues/test_batch_out/mcr_diamond_output.tsv'
    tmp_dir=r'/home/ec2-user/mouse_catalogues/test_tmp'
    #batch_diamond(input_dir, diamond_exe_path, database_path, type, tmp_dir, min_percent_identity=90, min_scov=0.9,
    #              database_fasta=None):
    '''
    dsr_diamond_output=run_diamond(input_file,diamond_exe_path,dsr_database_path,"dsr",tmp_dir,min_percent_identity=DSR_MIN_PERCENT_IDENTITY
                                   ,min_scov=DSR_MIN_SUBJECT_COVERAGE,database_fasta=dsr_database_fasta_file)

    #run diamond on the input file using the mcr database
    mcr_diamond_output=run_diamond(input_file,diamond_exe_path,mcr_database_path,"mcr",tmp_dir,min_percent_identity=MCR_MIN_PERCENT_IDENTITY
                                      ,min_scov=MCR_MIN_SUBJECT_COVERAGE,database_fasta=mcr_database_fasta_file)'''
    joined_records = join_and_add_labels_to_records(in_dir, os.path.join(tmp_dir, "hyd_joined_records.fasta"))
    batch_diamond(joined_records, default_diamond_exe_path, default_dsr_database_path, "dsr",tmp_dir ,
                  min_percent_identity=DSR_MIN_PERCENT_IDENTITY,
                  min_scov=DSR_MIN_SUBJECT_COVERAGE, database_fasta=dsr_database_fasta_file)
    batch_diamond(joined_records, default_diamond_exe_path, default_mcr_database_path, "mcr", r'/home/ec2-user/mouse_catalogues/test_tmp',
                    min_percent_identity=MCR_MIN_PERCENT_IDENTITY,
                    min_scov=MCR_MIN_SUBJECT_COVERAGE, database_fasta=mcr_database_fasta_file)
def test_acetoscan_batch():
    in_dir=r'/home/ec2-user/mouse_catalogues/test_batch/test_acetoscan'
    tmp_dir = r'/home/ec2-user/mouse_catalogues/test_tmp'
    joined_records = join_and_add_labels_to_records(in_dir, os.path.join(tmp_dir, "ace_joined_records.fasta"))
    out_ace=batch_acetoscan(joined_records,r'/home/ec2-user/mouse_catalogues/test_tmp')
def test_hydrogenase_batch():
    in_dir=r'/home/ec2-user/mouse_catalogues/test_batch'
    tmp_dir = r'/home/ec2-user/mouse_catalogues/test_tmp'
    joined_records = join_and_add_labels_to_records(in_dir, os.path.join(tmp_dir, "hyd_joined_records.fasta"))
    out_hyd=batch_run_hydrogenase_script(joined_records,default_hydrogenase_script_path,r'/home/ec2-user/mouse_catalogues/test_tmp')

def test_results_combine():
    results=[]
    input_dir=r'C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test'
    dsr_output=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\dsr_diamond_output_filtered.tsv"
    mcr_output=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\mcr_diamond_output_filtered.tsv"
    acetoscan_output=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\aceto_out.csv"
    hyDB_output=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\all_hydrogenase.txt"
    tmp_dir=r'C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\tmp1'
    try:
        for res in group_by_genome(input_dir, dsr_output, mcr_output, acetoscan_output, hyDB_output, tmp_dir):
            results.extend(res)
    except StopIteration:
        pass
    #write the results to a file
    with open(r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\results.csv",'w') as f:
        writer=csv.writer(f,delimiter='\t')
        writer.writerow(["filename","hydrogen_consumer","sulphate_reducer","methanogen","acetogen","consuming hydrogenase","bifurcating_or_bidirectional_hydrogenase",
                         "dsr_matches","mcr_matches","acetoscan_matches","hyDB_matches"])
        writer.writerows(results)
    output_file=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\results.csv"
    # open the file again and replace the uids with descriptions
    with open(output_file, 'r') as f:
        str1 = f.read()
    # replace the uids with descriptions
    Batch_Fetcher.fetch()
    repl = Batch_Fetcher.replace(str1)
    # write the replaced string to the file
    with open(output_file, 'w') as f:
        f.write(repl)
    return results


def test_batch_run2():
    test_data=r'C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\test_data'
    output_file=r'C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\test_data\results.csv'
    tmp_dir=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\test_data\tmp"
    args = argparse.ArgumentParser()
    args.add_argument("-i", "--input_dir", help="the input directory containing the genomes to be analyzed",default=test_data,
                      required=False)
    args.add_argument("-o", "--output_file", help="the output file to save the results", required=False,default=output_file)
    args.add_argument("-p", "--prodigal_exe_path", help="the path to the prodigal executable",
                      default=default_prodigal_path)
    args.add_argument("-d", "--diamond_exe_path", help="the path to the diamond executable",
                      default=default_diamond_exe_path)
    args.add_argument("-hy", "--hydrogenase_script_path", help="the path to the hydrogenase script",
                      default=default_hydrogenase_script_path)
    args.add_argument("-dsr", "--dsr_database_path", help="the path to the dsr database",
                      default=default_dsr_database_path)
    args.add_argument("-mcr", "--mcr_database_path", help="the path to the mcr database",
                      default=default_mcr_database_path)
    # add a flag for max threads
    args.add_argument("-th", "--max_threads", help="the maximum number of threads to use", default=10, type=int)
    # add a flag for batch size
    args.add_argument("-b", "--batch_size", help="the batch size", default=100, type=int)
    args.add_argument('-t', "--tmp_dir", help="the path to the temporary directory", default="./tmp")
    # add a flag to run the acetoscan
    args.add_argument("-a", "--acetoscan", help="run acetoscan", action="store_true", default=False)
    # add an argument for file extensions to be considered
    args.add_argument("-x", "--file_extensions", help="the file extensions to be considered",
                      default=[".fa", ".fna", ".fasta"], nargs="+")
    # add the path to the hydrogenase activity file
    args.add_argument("-hyd", "--hydrogen_classes_activity_file", help="the path to the hydrogen classes activity file",
                      default=DEFAULT_HYDROGENASE_ACTIVITY_FILE)
    # add a flag to enable resume mode, to continue in case of failure,its default is true
    args.add_argument("-r", "--resume",
                      help="resume mode, in case the previous run crash, use cached output from previous run",
                      action="store_true", default=True)
    args = args.parse_args()
    #parse the args
    #args = args.parse_args()
    batch_root_dir = batch_the_files(test_data, os.path.join(test_data, "batches"), batch_size=3)
    # run the batches
    # run the pipeline
    # batch_run2(args.input_dir,args.output_file,tmp_dir=args.tmp_dir, args=args)
    # run the batches
    # if max threads is 1, use the normal batch run
    #if args.max_threads == 1:

    #run_batches(batch_root_dir, output_file, tmp_dir,args)
    #else:
    run_batches_multithreaded(batch_root_dir, args.output_file, args.tmp_dir, args, max_batches=args.max_threads)
    # create the feature table
    #create_feature_table (tmp_filename,hydrogen_classes_activity_file,hydrogen_consumer_pipeline_output_folder,out_file)
    create_feature_table(tmp_dir, r"C:\Users\abel\Documents\hydrogen_consumers\hydrogenase\hydrogenase_classes_activity2t.csv",
                         test_data, output_file)



if __name__=="__main__":
    #test_batch_run2()
    #exit()
    main()
    #tes()
    #tes2()
    #tes_diamond_batch()
    #test_acetoscan_batch()
    #test_hydrogenase_batch()
    #test_results_combine()
    #test_batch_run2()
