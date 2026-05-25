


import argparse
import concurrent
import sys
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




default_cdd_search_path=r'tllhome/abel/hydrogen_consumer_pipeline/hydrogen_consumer_pipeline/cdd_search.py'
#if default cdd search path not exist, use the alternative: /oceanstor/home/e1103389/hydrogen_consumer_pipeline
if not os.path.exists(default_cdd_search_path):
    default_cdd_search_path=r'/oceanstor/home/e1103389/hydrogen_consumer_pipeline/cdd_search.py'


#get the parent dir of the script
parent_dir= os.path.dirname(os.path.realpath(__file__))
#get the parent of the parent dir
best_results_file=os.path.join(parent_dir,'all_best_df.csv')
parent_dir=os.path.dirname(parent_dir)


def parse_best_results(best_results_file):
    df=pd.read_csv(best_results_file,sep=',')
    #filter for train_percent =0.2
    df=df[df['train_percent']==0.2]
    threshold_dict={}
    for index,row in df.iterrows():
        threshold_dict[row['gene']]=(row['pid'],row['scov'])
    return threshold_dict


thresholds=parse_best_results(best_results_file)




#r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\cdd_search.py"
#Refactor to use parent dir as /home/ec2-user/
#acetobase_dir=r'/home/ec2-user/acetoscan/AcetoBase.fasta'
acetobase_dir=os.path.join(parent_dir,'acetoscan','AcetoBase.fasta')
#acetoscan_input_dir=r'/home/ec2-user/acetoscan/genome'
acetoscan_input_dir=os.path.join(parent_dir,'acetoscan','genome')
#hyDb_input_dir=r'/home/ec2-user/hyDB/genome'
#hyDb_output_dir=r'/home/ec2-user/hydrogenae_Scan/genome_output'
# diamond_exe_path,hydrogenase_script_path,
#             dsr_database_path, mcr_database_path

default_diamond_exe_path='/tllhome/abel/DIAMOND/diamond'
#refactor to use parent dir as tllhome/abel
default_diamond_exe_path=os.path.join(parent_dir,'DIAMOND','diamond')

#r"D:\diamond-windows\diamond.exe"
#default_hydrogenase_script_path=r"C:\Users\abel\Documents\hydrogenases\hyDB\Diamond_blastp_hyDB.sh"
#now the hydrogenase script is a python script
default_hydrogenase_script_path=r'/tllhome/abel/hydrogen_consumer_pipeline/hyDB/Diamond_blast_hyDB2.py'

#refactor to use parent dir as /tllhome/abel/hydrogen_consumer_pipeline
default_hydrogenase_script_path=os.path.join(parent_dir,'hyDB','Diamond_blast_hyDB2.py')
#r"C:\Users\abel\Documents\hydrogenases\hyDB\Diamond_blast_hyDB.py"
#default_DIADB_dir=r'/tllhome/abel/hydrogen_consumer_pipeline/hyDB'
#refactor to use parent dir as /tllhome/abel/hydrogen_consumer_pipeline
default_DIADB_dir=os.path.join(parent_dir,'hyDB')
#default_prodigal_path=r'/tllhome/abel/hydrogen_consumer_pipeline/prodigal.linux'
#refactor to use parent dir as /tllhome/abel/hydrogen_consumer_pipeline
default_prodigal_path=os.path.join(parent_dir,'prodigal.linux')
#if the prodigal path does not exist,just use prodigal
if not os.path.exists(default_prodigal_path):
    default_prodigal_path='prodigal'

#r"D:\prodigal\prodigal.windows.exe"
default_dsr_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/dsrA_curated.fasta"

default_mcr_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/mcr1_curated2.fasta"
default_dsr_database_path=r'/home/ec2-user/hyDB/dsrA.dmnd'
default_mcr_database_path=r'/home/ec2-user/hyDB/mcr1.dmnd'
#refactor to use parent dir as /tllhome/abel/hydrogen_consumer_pipeline
default_dsr_database_fasta_file=os.path.join(parent_dir,'dsrA_curated.fasta')
default_mcr_database_fasta_file=os.path.join(parent_dir,'mcr1_curated2.fasta')
default_dsr_database_path=os.path.join(parent_dir,'hyDB','dsrA.dmnd')
default_mcr_database_path=os.path.join(parent_dir,'hyDB','mcr1.dmnd')
default_acsb_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/acsB.fasta"
#if the acsB database fasta file does not exist, use the default
if not os.path.exists(default_acsb_database_fasta_file):
    default_acsb_database_fasta_file=os.path.join(parent_dir,'DIAMOND','acsB.fasta')
#aprA
default_aprA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/aprA.fasta"
#if the aprA database fasta file does not exist, use the default
if not os.path.exists(default_aprA_database_fasta_file):
    default_aprA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','aprA.fasta')
#asrA
default_asrA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/asrA.fasta"
#if the asrA database fasta file does not exist, use the default
if not os.path.exists(default_asrA_database_fasta_file):
    default_asrA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','asrA.fasta')
#cydA
default_cydA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/cydA.fasta"
#if the cydA database fasta file does not exist, use the default
if not os.path.exists(default_cydA_database_fasta_file):
    default_cydA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','cydA.fasta')

#DmsA
default_dmsA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/DmsA.fasta"
#if the dmsA database fasta file does not exist, use the default
if not os.path.exists(default_dmsA_database_fasta_file):
    default_dmsA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','DmsA.fasta')

#dsrA
default_dsrA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/dsrA.fasta"
#if the dsrA database fasta file does not exist, use the default
if not os.path.exists(default_dsrA_database_fasta_file):
    default_dsrA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','dsrA.fasta')

#frdA
default_frdA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/frdA.fasta"
#if the frdA database fasta file does not exist, use the default
if not os.path.exists(default_frdA_database_fasta_file):
    default_frdA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','frdA.fasta')
#mcrA
default_mcrA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/mcrA.fasta"
#if the mcrA database fasta file does not exist, use the default
if not os.path.exists(default_mcrA_database_fasta_file):
    default_mcrA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','mcrA.fasta')
#napA
default_napA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/napA.fasta"
#if the napA database fasta file does not exist, use the default
if not os.path.exists(default_napA_database_fasta_file):
    default_napA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','napA.fasta')
#narG
default_narG_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/narG.fasta"
#if the narG database fasta file does not exist, use the default
if not os.path.exists(default_narG_database_fasta_file):
    default_narG_database_fasta_file=os.path.join(parent_dir,'DIAMOND','narG.fasta')
#nrfA
default_nrfA_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/nrfA.fasta"
#if the nrfA database fasta file does not exist, use the default
if not os.path.exists(default_nrfA_database_fasta_file):
    default_nrfA_database_fasta_file=os.path.join(parent_dir,'DIAMOND','nrfA.fasta')
#fhl
default_fhl_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/FHL.fasta"
#if the fhl database fasta file does not exist, use the default
if not os.path.exists(default_fhl_database_fasta_file):
    default_fhl_database_fasta_file=os.path.join(parent_dir,'DIAMOND','FHL.fasta')
#pfor
default_pfor_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/pfor.fasta"
#if the pfor database fasta file does not exist, use the default
if not os.path.exists(default_pfor_database_fasta_file):
    default_pfor_database_fasta_file=os.path.join(parent_dir,'DIAMOND','pfor.fasta')
#nitrogenase
default_nitrogenase_database_fasta_file=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/nifH.fasta"
#if the nitrogenase database fasta file does not exist, use the default
if not os.path.exists(default_nitrogenase_database_fasta_file):
    default_nitrogenase_database_fasta_file=os.path.join(parent_dir,'DIAMOND','nifH.fasta')
#default_dsr_database_path=r'/tllhome/abel/hydrogen_consumer_pipeline/hyDB/dsrA.dmnd'
MCR_MIN_PERCENT_IDENTITY=90
MCR_MIN_SUBJECT_COVERAGE=0.9
DSR_MIN_PERCENT_IDENTITY=50
DSR_MIN_SUBJECT_COVERAGE=0.9

default_cat_path="/tllhome/abel/CAT/CAT/CAT_pack/CAT"
#refactor to use parent dir as  /tllhome/abel/hydrogen_consumer_pipeline/CAT/CAT_pack/CAT_pack/CAT_pack
default_cat_path=os.path.join(parent_dir,'CAT','CAT_pack','CAT_pack','CAT_pack')
default_CAT_taxonomy_path="/tlldata/HenningLab/abel/CAT/CAT_prepare_20210107/2021-01-07_taxonomy"
default_CAT_database_path="/tlldata/HenningLab/abel/CAT/CAT_prepare_20210107/2021-01-07_CAT_database"
#refactor to use parent dir as /tlldata/HenningLab/abel/CAT/CAT_prepare_20210107/2021-01-07_taxonomy
CAT_taxonomy_path=os.path.join(parent_dir,'CAT','CAT_prepare_20210107','2021-01-07_taxonomy')
#refactor to use parent dir as /tlldata/HenningLab/abel/CAT/CAT_prepare_20210107/2021-01-07_CAT_database
CAT_database_path=os.path.join(parent_dir,'CAT','CAT_prepare_20210107','2021-01-07_CAT_database')
#the create feature script is in the same directory as this script

DEFAULT_HYDROGENASE_ACTIVITY_FILE=r"/tllhome/abel/hydrogenases/hydrogenase_classes_activity.csv"

#refactor to use parent dir as /tllhome/abel/hydrogenases/hydrogenase_classes_activity.csv
DEFAULT_HYDROGENASE_ACTIVITY_FILE=os.path.join(parent_dir,'hydrogenase_classes_activity.csv')
#if DEFAULT_HYDROGENASE_ACTIVITY_FILE does not exist, use the alternative
if not os.path.exists(DEFAULT_HYDROGENASE_ACTIVITY_FILE):
    DEFAULT_HYDROGENASE_ACTIVITY_FILE=r'/oceanstor/home/e1103389/hyDB/hydrogenase_classes_activity.csv'

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
                    print("Error occurred for batch", i)
                    print(e)
                    tries+=1
            if tries==5:
                print("Error occurred for batch", i)
                #print(e)
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

'''def is_sulphate_reducer(dsr_diamond_output):
    with open(dsr_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(dsr_diamond_output)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}'''

def is_sulphate_reducer(non_hyd_outputs,args):
    dsrA_diamond_output=non_hyd_outputs["dsrA"]
    aprA_diamond_output=non_hyd_outputs["aprA"]
    asrA_diamond_output=non_hyd_outputs["asrA"]
    #for output in [dsrA_diamond_output,aprA_diamond_output,asrA_diamond_output]:
    readerlist={}
    with open(dsrA_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist.update(process_diamond_output(dsrA_diamond_output,database_file=args.dsrA_database_fasta_file))
    with open(aprA_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist.update(process_diamond_output(aprA_diamond_output,database_file=args.aprA_database_fasta_file))
    with open(asrA_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist.update(process_diamond_output(asrA_diamond_output,database_file=args.asrA_database_fasta_file))
    if len(readerlist)>0:
        return True,readerlist
    return False,{}

def is_cydA(non_hyd_outputs,args):
    cydA_diamond_output=non_hyd_outputs["cydA"]
    with open(cydA_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(cydA_diamond_output,database_file=args.cydA_database_fasta_file)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}


def is_methanogen(non_hyd_outputs,args):
    mcrA_diamond_output=non_hyd_outputs["mcrA"]
    with open(mcrA_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(mcrA_diamond_output,database_file=args.mcrA_database_fasta_file)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}
def is_acetogen(non_hyd_outputs,args):
    acsB_diamond_output=non_hyd_outputs["acsB"]

    with open(acsB_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(acsB_diamond_output,database_file=args.acsB_database_fasta_file)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}
def is_nitrate_reducer(non_hyd_outputs,args):
    narg_diamond_output=non_hyd_outputs["narG"]
    napA_diamond_output=non_hyd_outputs["napA"]
    nrfA_diamond_output=non_hyd_outputs["nrfA"]
    readerlist={}
    #for output,db_file in [narg_diamond_output,napA_diamond_output,nrfA_diamond_output]:
    for output,db_file in [(narg_diamond_output,args.narG_database_fasta_file),(napA_diamond_output,args.napA_database_fasta_file),(nrfA_diamond_output,args.nrfA_database_fasta_file)]:
        with open(output) as f:
            reader=csv.reader(f,delimiter='\t')
            readerlist.update(process_diamond_output(output,database_file=db_file))
    if len(readerlist)>0:
        return True,readerlist
    return False,{}
def is_fumarate_reducer(non_hyd_outputs,args):
    frdA_diamond_output=non_hyd_outputs["frdA"]
    with open(frdA_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(frdA_diamond_output,database_file=args.frdA_database_fasta_file)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}
def is_dms_reducer(non_hyd_outputs,args):
    dmsA_diamond_output=non_hyd_outputs["dmsA"]
    with open(dmsA_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(dmsA_diamond_output,database_file=args.dmsA_database_fasta_file)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}
def is_pforl(non_hyd_outputs,args):
    pfor_diamond_output=non_hyd_outputs["pfor"]
    with open(pfor_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(pfor_diamond_output,database_file=args.pfor_database_fasta_file)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}
def is_fhll(non_hyd_outputs,args):
    fhl_diamond_output=non_hyd_outputs["fhl"]
    with open(fhl_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(fhl_diamond_output,database_file=args.fhl_database_fasta_file)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}
def is_nitrogenase(non_hyd_outputs,args):
    nifH_diamond_output=non_hyd_outputs["nifH"]
    with open(nifH_diamond_output) as f:
        reader=csv.reader(f,delimiter='\t')
        readerlist=process_diamond_output(nifH_diamond_output,database_file=args.nitrogenase_database_fasta_file)
        if len(readerlist)>0:
            return True,readerlist
        else:
            return False,{}




def process_diamond_output(diamond_output_file,bit_score_index=-1,inject_desc=True,database_file=None):
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
        #if the database file is not none, then fetch the description from the database file
        if database_file is not None:
            with open(database_file) as f:
                recs=SeqIO.parse(f,'fasta')
                for rec in recs:
                    if rec.id==subject:
                        desc=rec.description
                        break
        #if subject contains " - ", split it and take the first part
        else:
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

    


'''hyd_data = {
    "[FeFe] Group A1": "H2-evolution",
    "[FeFe] Group A2": "H2-uptake ?",
    "[FeFe] Group A3": "Electron-bifurcation",
    "[FeFe] Group A4": "H2-evolution or electron-bifurcation",
    "[FeFe] Group B" : "H2-evolution?",
    "[FeFe] Group C1": "H2-sensing?",
    "[FeFe] Group C2": "H2-sensing?",
    "[FeFe] Group C3": "H2-sensing?",
    '[FeFe] NOT FOUND':"???",
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
}'''
def parse_hydrogenase_data(hyd_data_file):
    with open(hyd_data_file) as f:
        reader=csv.reader(f,delimiter=',')
        next(reader)
        hyd_data={}
        for row in reader:
            if row[0]=="Fe":
                index=f"[{row[0]}]"
            else:
                index=f"[{row[0]}] Group {row[1]}"
            hyd_data[index]=row[2]
    return hyd_data
def is_hydrogenase2(hyDB_output,hyd_data_file=DEFAULT_HYDROGENASE_ACTIVITY_FILE):
    hyd_data=parse_hydrogenase_data(hyd_data_file)
    processed_output=process_diamond_output(hyDB_output,bit_score_index=-4,inject_desc=False)
    unique_queries={}
    at_least_one_uptake = False
    at_least_one_bifurcating_or_bidirectional = False
    at_least_one_producing=False
    for query in processed_output:
        hyd_type=processed_output[query]
        hyd_types = hyd_type.split('_-_')
        # get the last one part and remove the fina
        hyd_type = hyd_types[-1]
        # find the hyd type in the dict
        #replace underscores with spaces
        hyd_type=hyd_type.replace("_"," ")
        if hyd_type in hyd_data:
            value = hyd_data[hyd_type]
            if "uptake" in value:
                at_least_one_uptake = True
                unique_queries[query] = processed_output[query]
            if "bifurcation" in value or "bidirectional" in value:
                at_least_one_bifurcating_or_bidirectional = True
                unique_queries[query] = processed_output[query]
            else:
                at_least_one_producing=True
                unique_queries[query] = processed_output[query]
    return at_least_one_uptake,at_least_one_bifurcating_or_bidirectional,at_least_one_producing,unique_queries


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

def is_hydrogen_consumer(query,non_hyd_output, hyDB_output,args):
    #check if it is a sulphate reducer
    sulphate_reducer=False
    methanogen=False
    acetogen=False
    cyda=False
    nitrate_reducer=False
    fumarate_reducer=False
    dms_reducer=False
    pfor=False
    fhl=False

    hydrogenase=False
    is_sulphate, unique_dictd = is_sulphate_reducer(non_hyd_output,args)
    is_methano, unique_dictm = is_methanogen(non_hyd_output,args)
    is_aceto, unique_dicta = is_acetogen(non_hyd_output,args)
    is_cyd, unique_dictcyd = is_cydA(non_hyd_output,args)
    is_nitrate, unique_dictn = is_nitrate_reducer(non_hyd_output,args)
    is_fumarate, unique_dictf = is_fumarate_reducer(non_hyd_output,args)
    is_dms, unique_dictdm = is_dms_reducer(non_hyd_output,args)
    is_pfor, unique_dictp = is_pforl(non_hyd_output,args)
    is_fhl, unique_dictfhl = is_fhll(non_hyd_output,args)
    is_nifh, unique_dictnif = is_nitrogenase(non_hyd_output,args)
    #TODO: uncomment this after debugging
    is_consume, is_bifurcating,is_producing, unique_dicth = is_hydrogenase2(hyDB_output)
    is_hydrogen_consumer=False
    if is_consume:
        print("contains a consumption hydrogenase")
        is_hydrogen_consumer=True
    front_row={"query":query,"sulphate_reducer":is_sulphate, "methanogen":is_methano, "acetogen":is_aceto,
                                           "cydA":is_cyd, "nitrate_reducer":is_nitrate, "fumarate_reducer":is_fumarate, "dms_reducer":is_dms,
                                           "pfor":is_pfor, "fhl":is_fhl, "hydrogen_consumer":is_hydrogen_consumer,
                                           "hydrogen_producer":is_producing,"nitrogenase":is_nifh}
    return front_row,{"sulphate_reducer_matches":unique_dictd, "methanogen_matches":unique_dictm,
                       "acetogen_matches":unique_dicta,
                         "cydA_matches":unique_dictcyd, "nitrate_reducer_matches":unique_dictn,
                         "fumarate_reducer_matches":unique_dictf, "dms_reducer_matches":unique_dictdm,
                         "pfor_matches":unique_dictp, "fhl_matches":unique_dictfhl, "hydrogenase_matches":unique_dicth,"nitrogenase_matches":unique_dictnif}
    #make the evidence rows
    evidence_rows=make_evidence_rows(unique_dictd,unique_dictm,unique_dicta,unique_dicth,
                                     unique_dictcyd,unique_dictn,unique_dictf,unique_dictdm,unique_dictp,unique_dictfhl)
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


def make_evidence_rows(unique_dictd,unique_dictm,unique_dicta,unique_dicth,
                       unique_dictcyd,unique_dictn,unique_dictf,unique_dictdm,unique_dictp,unique_dictfhl):
    #get the dict with the max rows
    #lengths=[len(unique_dictd),len(unique_dictm),len(unique_dicta),len(unique_dicth),
    #            len(unique_dictcyd),len(unique_dictn),len(unique_dictf),len(unique_dictdm),len(unique_dictp),len(unique_dictfhl)]
    dicts=[unique_dictd,unique_dictm,unique_dicta,unique_dicth,
           unique_dictcyd,unique_dictn,unique_dictf,unique_dictdm,unique_dictp,unique_dictfhl]
    lengths=[len(d) for d in dicts]
    lisd=list(unique_dictd.values())
    lism=list(unique_dictm.values())
    lisa=list(unique_dicta.values())
    lish=list(unique_dicth.values())

    max_length=max(lengths)

    print("LISA:",lisa)
    lish=list(unique_dicth.values())
    max_length=max(lengths)
    rows=[]
    for r in range(max_length):
        row=[]
        for d in dicts:
            if r<len(d):
                row.append(d[r])
            else:
                row.append('')
        '''if r<len(lisd):
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
            row.append('')'''
        rows.append(row)
    return rows

def run_diamond(input_file,diamond_exe_path,database_path,type,tmp_dir,min_percent_identity=70,min_scov=0.5,database_fasta=None,tag=None):

    #if tag is not none, use the threshold file
    if tag is not None:
        thresh=thresholds[tag]
        min_percent_identity=thresh[0]
        min_scov=thresh[1]

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
        print("FILEF:",os.path.join(input_dir,file))
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
def create_diamond_database(fasta_file, tmp_folder,diamond_path):
    #create a diamond database from the fasta file
    #get the name of the fasta file
    name=os.path.basename(fasta_file)
    #get the name without the extension
    name=os.path.splitext(name)[0]+".dmnd"
    #add the name to the tmp folder
    name=os.path.join(tmp_folder,name)
    #create the database
    command=diamond_path+' makedb --in '+fasta_file+' -d '+name
    os.system(command)
    #return the path to the database
    return name
def batch_diamond(joined_records, diamond_exe_path, database_fasta, type, tmp_dir, min_percent_identity=70, min_scov=0.5, resume=True,TAG=None):
    #if tag is not none, use the threshold file
    if TAG is not None:
        thresh=thresholds[TAG]
        min_percent_identity=thresh[0]
        min_scov=thresh[1]
    database_path=database_fasta.replace('.fasta','.dmnd')
    #if database path not exists, create a diamond database from the fasta file
    if not os.path.exists(database_path):
        database_path=create_diamond_database(database_fasta,tmp_dir,diamond_exe_path)
    #if database path is a fasta file, create a diamond database from it
    #if database_path.endswith('.fasta'):
        #create a diamond database from the fasta file
    #    database_path=create_diamond_database(database_path,tmp_dir,diamond_exe_path)
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

def batch_acetoscan(input_dir,tmp_dir,resume=True,run_acetoscan=True):
    #if run acetoscan is false, then return None
    if not run_acetoscan:
        return "None"
    #if theres a file called _all_acetoscan.csv in the tmp dir, use it
    if resume:
        #find the file in the tmp dir
        for file in os.listdir(tmp_dir):
            if file.endswith("_all_acetoscan.csv"):
                print("Using cached in tmp dir:",file)
                return os.path.join(tmp_dir,file)
    #run the script

    parent_dir= os.path.dirname(__file__)
    #the parent dir is the parent of the parent dir
    parent_dir=os.path.dirname(parent_dir)
    script_path=os.path.join(parent_dir,'Acetoscan','run_acetoscan_HPC.py')
    #add the script path to the import path
    sys.path.append(os.path.dirname(script_path))
    print("sys path:",sys.path)
    #import the script
    import run_acetoscan_HPC
    output_file=run_acetoscan_HPC.run_acetoscan(input_dir,tmp_dir,None)
    return output_file
    

def batch_run_hydrogenase_script(ijoined_records,hydrogenase_script_path,tmp_dir,resume=True,diamond_exe_path=default_diamond_exe_path,diadb_dir=default_DIADB_dir):
    print("Running hydrogenase script")
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
        print("Return ",os.path.join(output_filepath, "Filtered", "all_hydrogenase.txt"))
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
    #if new name file already exists, remove it
    if os.path.exists(os.path.join(hyinput_dir,new_name)):
        os.remove(os.path.join(hyinput_dir,new_name))
    os.rename(os.path.join(hyinput_dir,os.path.basename(ijoined_records)),os.path.join(hyinput_dir,new_name))

    #mkdir the output dir if it does not exist, else if not resume, remove the dir and make it
    if os.path.exists(output_filepath) and not resume:
        shutil.rmtree(output_filepath)
    os.makedirs(output_filepath,exist_ok=True)
    print("MAKE DIR:",output_filepath)

    #the last argument is the diamond exe path
    command='python '+hydrogenase_script_path+' '+hyinput_dir+' '+output_filepath+' '+diamond_exe_path
    print("COMMAND:",command)
    os.system(command)
    #get the output file
    output_file=os.path.join(output_filepath,"Filtered","all_hydrogenase.txt")
    #move the output file to the tmp dir
    shutil.copy(output_file,os.path.join(tmp_dir,os.path.basename(output_file)))
    #return the output file
    output_file=os.path.join(tmp_dir,"all_hydrogenase.txt")
    print("Return ",output_file)
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



def group_by_genome(input_dir, non_hyd_outputs, hyDB_output,tmp_dir,extensions,args):
    print(f"group by genome ({input_dir}, {non_hyd_outputs}, {hyDB_output}, {tmp_dir})")
    #if tmp dir not exists, make it, else, remove all files in it
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir,exist_ok=True)
    #else:
        #shutil.rmtree(tmp_dir)
        #os.makedirs(tmp_dir,exist_ok=True)
    #for each genome, get the dsr, mcr, acetoscan and hyDB results
    groups={}
    non_hyds={}
    for key, non_hyd_output in non_hyd_outputs.items():
        non_hyds[key]=[]
        with open(non_hyd_output) as f:
            reader=csv.reader(f,delimiter='\t')
            for row in reader:
                non_hyds[key].append(row)

    #open the hyDB output file
    hyDB_rows=[]
    print("hyDB output:",hyDB_output)
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
        non_hyd_res={}
        for key in non_hyds:
            non_hyd_res[key]=[]
            for row in non_hyds[key]:
                #if len of row is 0, skip
                if len(row)==0:
                    continue
                if row[0].startswith(file_no_ext):
                    non_hyd_res[key].append(row)
        #for all mcr rows, find the row[0] that starts with the file name
        '''mcr_res=[]
        for row in mcr_rows:
            #skip if row is empty
            if len(row)==0:
                continue
            if row[0].startswith(file_no_ext):
                mcr_res.append(row)'''
        #for all hyDB rows, find the row[0] that starts with the file name
        hyDB_res=[]
        for row in hyDB_rows:
            #skip if row is empty
            if len(row)==0:
                continue
            if row[0].startswith(file_no_ext):
                hyDB_res.append(row)
        #add the results to the groups
        groups[file]=(non_hyd_res,hyDB_res)
        #groups[file]=(dsr_res,mcr_res,acetoscan_res,hyDB_res)
    #for each group, write each result to a separate file
    for i,group in enumerate(groups):
        file_dict={}
        print("GENOME:",i,"/",len(groups)," ",group)
        #dsr_rows=groups[group][0]
        '''#write the dsr results to a file
        with open(os.path.join(tmp_dir,group+"_dsr_matches.tsv"),'w') as f:
            writer=csv.writer(f,delimiter='\t')
            writer.writerows(dsr_rows)'''
        for key in groups[group][0]:
            with open(os.path.join (tmp_dir,group+"_"+key+"_matches.tsv"),'w') as f:
                writer=csv.writer(f,delimiter='\t')
                writer.writerows(groups[group][0][key])
            file_dict[key]=os.path.join(tmp_dir,group+"_"+key+"_matches.tsv")
        hyDB_rows=groups[group][1]
        #write the hyDB results to a file
        with open(os.path.join(tmp_dir,group+"_hyDB_matches.tsv"),'w') as f:
            writer=csv.writer(f,delimiter='\t')
            writer.writerows(hyDB_rows)
        #run the is_hydrogen_consumer function on the group
        res,dicts=is_hydrogen_consumer(group,file_dict,os.path.join(tmp_dir,group+"_hyDB_matches.tsv"),args)
        print("DICTS:",dicts)
        #header=list(res.keys())+["protein_match"]+ list(dicts.keys())
        
        res_rows=[]
        header= ['query', 'sulphate_reducer', 'methanogen', 'acetogen', 'cydA', 'nitrate_reducer', 'fumarate_reducer', 'dms_reducer', 'pfor', 'fhl', 
                 'hydrogen_consumer', 'hydrogen_producer', 'nitrogenase', 'protein_match',
                   'sulphate_reducer_matches', 'methanogen_matches', 'acetogen_matches', 'cydA_matches', 'nitrate_reducer_matches', 
                   'fumarate_reducer_matches', 'dms_reducer_matches', 'pfor_matches', 'fhl_matches', 'hydrogenase_matches', 'nitrogenase_matches']
        headerp1=['query', 'sulphate_reducer', 'methanogen', 'acetogen', 'cydA', 'nitrate_reducer', 'fumarate_reducer', 'dms_reducer', 'pfor', 'fhl', 
                 'hydrogen_consumer', 'hydrogen_producer', 'nitrogenase']
        headerp2=['sulphate_reducer_matches', 'methanogen_matches', 'acetogen_matches', 'cydA_matches', 'nitrate_reducer_matches',
                    'fumarate_reducer_matches', 'dms_reducer_matches', 'pfor_matches', 'fhl_matches', 'hydrogenase_matches', 'nitrogenase_matches']
        res_row1=[]
        for key in headerp1:
            res_row1.append(res[key])
        #add the protein match
        #res_row1.append(res['protein_match'])
        row=res_row1
        for ni,key in enumerate(headerp2):
            for ele in dicts[key].values():
                row=[ele]+[""]*(ni)+[ele]+[""]*(len(dicts)-ni-1)
                rr=res_row1+row
                for i in range(len(rr)):
                    pass
                    #print('RESS:',i,header[i],rr[i])
                res_rows.append(res_row1+row)
        #write the results to a file
        with open(os.path.join(tmp_dir,group+"_results.tsv"),'w') as f:
            print(os.path.join(tmp_dir,group+"_results.tsv"))
            writer=csv.writer(f,delimiter='\t')
            #writer.writerow(header)
            writer.writerows(res_rows)
        #print("RES ROWS",res_rows)
        yield res_rows,header




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
                print("setting prodigal path to default:",prodigal_path)
            else:
                prodigal_path=args.prodigal_exe_path
                print("setting prodigal path to args:",prodigal_path)

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
    #mcr_output=batch_diamond(prodigal_joined,default_diamond_exe_path,mcr_database_fasta_file,"mcr",tmp_dir,min_percent_identity=90,min_scov=0.9,
    #                         resume=args.resume)
    #run the batch hyDB
    #batch_run_hydrogenase_script(ijoined_records,hydrogenase_script_path,tmp_dir)
    #if args.hydrogenase_script_path is None, set it to the default
    if args.hydrogenase_script_path is None:
        hydrogenase_script_path=default_hydrogenase_script_path
    else:
        hydrogenase_script_path=args.hydrogenase_script_path
    hyDB_output=batch_run_hydrogenase_script(prodigal_joined,hydrogenase_script_path,tmp_dir,resume=args.resume,diamond_exe_path=args.diamond_exe_path,diadb_dir=args.diadb_path)
    #run the batch acetoscan
    
    
    default_diamond_exe_path=diamond_exe_path
    #acsB
    if args.acsB_database_fasta_file is None:
        acsB_database_fasta_file=default_acsb_database_fasta_file
    else:
        acsB_database_fasta_file=args.acsB_database_fasta_file
    acsB_output= batch_diamond(prodigal_joined,default_diamond_exe_path,acsB_database_fasta_file,"acsB",tmp_dir,min_percent_identity=90,min_scov=0.9,
                             resume=args.resume,TAG="AcsB")
    
    #aprA
    if args.aprA_database_fasta_file is None:
        aprA_database_fasta_file=default_aprA_database_fasta_file
    else:
        aprA_database_fasta_file=args.aprA_database_fasta_file
    aprA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,aprA_database_fasta_file,"aprA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                             resume=args.resume,TAG="AprA")
    #asrA
    if args.asrA_database_fasta_file is None:
        asrA_database_fasta_file=default_asrA_database_fasta_file
    else:
        asrA_database_fasta_file=args.asrA_database_fasta_file
    asrA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,asrA_database_fasta_file,"asrA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                             resume=args.resume,TAG=None)
    #cydA
    if args.cydA_database_fasta_file is None:
        cydA_database_fasta_file=default_cydA_database_fasta_file
    else:
        cydA_database_fasta_file=args.cydA_database_fasta_file
    cydA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,cydA_database_fasta_file,"cydA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="CydA")
    #DmsA
    if args.dmsA_database_fasta_file is None:
        dmsA_database_fasta_file=default_dmsA_database_fasta_file
    else:
        dmsA_database_fasta_file=args.dmsA_database_fasta_file
    dmsA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,dmsA_database_fasta_file,"dmsA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="DmsA")
    #dsrA
    if args.dsrA_database_fasta_file is None:
        dsrA_database_fasta_file=default_dsrA_database_fasta_file
    else:
        dsrA_database_fasta_file=args.dsrA_database_fasta_file
    dsrA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,dsrA_database_fasta_file,"dsrA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="DsrA")
    #frdA
    if args.frdA_database_fasta_file is None:
        frdA_database_fasta_file=default_frdA_database_fasta_file
    else:
        frdA_database_fasta_file=args.frdA_database_fasta_file
    frdA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,frdA_database_fasta_file,"frdA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="FrdA")
    
    #mcrA
    if args.mcrA_database_fasta_file is None:
        mcrA_database_fasta_file=default_mcrA_database_fasta_file
    else:
        mcrA_database_fasta_file=args.mcrA_database_fasta_file
    mcrA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,mcrA_database_fasta_file,"mcrA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="McrA")
    #napA
    if args.napA_database_fasta_file is None:
        napA_database_fasta_file=default_napA_database_fasta_file
    else:
        napA_database_fasta_file=args.napA_database_fasta_file
    napA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,napA_database_fasta_file,"napA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="NapA")
    #narG
    if args.narG_database_fasta_file is None:
        narG_database_fasta_file=default_narG_database_fasta_file
    else:
        narG_database_fasta_file=args.narG_database_fasta_file
    narG_output= batch_diamond(prodigal_joined,default_diamond_exe_path,narG_database_fasta_file,"narG",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="NarG")
    #nrfA
    if args.nrfA_database_fasta_file is None:
        nrfA_database_fasta_file=default_nrfA_database_fasta_file
    else:
        nrfA_database_fasta_file=args.nrfA_database_fasta_file
    nrfA_output= batch_diamond(prodigal_joined,default_diamond_exe_path,nrfA_database_fasta_file,"nrfA",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="NrfA")
    #fhl
    if args.fhl_database_fasta_file is None:
        fhl_database_fasta_file=default_fhl_database_fasta_file
    else:
        fhl_database_fasta_file=args.fhl_database_fasta_file
    fhl_output= batch_diamond(prodigal_joined,default_diamond_exe_path,fhl_database_fasta_file,"fhl",tmp_dir,min_percent_identity=90,min_scov=0.9,
                             resume=args.resume,TAG="FHL")
    #pfor
    if args.pfor_database_fasta_file is None:
        pfor_database_fasta_file=default_pfor_database_fasta_file
    else:
        pfor_database_fasta_file=args.pfor_database_fasta_file
    pfor_output= batch_diamond(prodigal_joined,default_diamond_exe_path,pfor_database_fasta_file,"pfor",tmp_dir,min_percent_identity=90,min_scov=0.9,
                             resume=args.resume,TAG="PFOR")
    #nitrogenase
    if args.nitrogenase_database_fasta_file is None:
        nitrogenase_database_fasta_file=default_nitrogenase_database_fasta_file
    else:
        nitrogenase_database_fasta_file=args.nitrogenase_database_fasta_file
    nitrogenase_output= batch_diamond(prodigal_joined,default_diamond_exe_path,nitrogenase_database_fasta_file,"nitrogenase",tmp_dir,min_percent_identity=90,min_scov=0.9,
                                resume=args.resume,TAG="NifH")
    
    #acetoscan is now deprecated in favor of acsb    




    non_hyd_outputs=[acsB_output,aprA_output,asrA_output,cydA_output,dmsA_output,
                     dsrA_output,frdA_output,mcrA_output,napA_output,narG_output,nrfA_output,fhl_output,pfor_output,nitrogenase_output]
    #make it into a dict
    non_hyd_outputs={"acsB":acsB_output,"aprA":aprA_output,"asrA":asrA_output,"cydA":cydA_output,"dmsA":dmsA_output,
                     "dsrA":dsrA_output,"frdA":frdA_output,"mcrA":mcrA_output,"napA":napA_output,"narG":narG_output,"nrfA":nrfA_output,"fhl":fhl_output,"pfor":pfor_output,"nifH":nitrogenase_output}
    #batch_acetoscan(dna_joined,tmp_dir, run_acetoscan=args.acetoscan,resume=args.resume, results=args.acetoscan_results_dir)
    results=[]
    #group_by_genome(input_dir, dsr_output, mcr_output, acetoscan_output, hyDB_output,tmp_dir)

    #if args run the neighbour hood search,
    if args.neighbourhood_search:
        print("neighbourhood search")
        res=run_neighbourhood_search(input_dir,args)
        #set the hydDB output to the neighbourhood search output
        hyDB_output=res
    header1= ['query', 'sulphate_reducer', 'methanogen', 'acetogen', 'cydA', 'nitrate_reducer', 'fumarate_reducer', 'dms_reducer', 'pfor', 'fhl', 
                 'hydrogen_consumer', 'hydrogen_producer', 'nitrogenase', 'protein_match',
                   'sulphate_reducer_matches', 'methanogen_matches', 'acetogen_matches', 'cydA_matches', 'nitrate_reducer_matches', 
                   'fumarate_reducer_matches', 'dms_reducer_matches', 'pfor_matches', 'fhl_matches', 'hydrogenase_matches', 'nitrogenase_matches']

        
    for res, header1 in group_by_genome(input_dir, non_hyd_outputs,  hyDB_output,tmp_dir,args.file_extensions,args):
        results.extend(res)
    
    print("HEADER1:",header1)
    # check if the output file already exists, if exists, append to it
    header = False
    
    #print("RESULTS[0]",results[0])
    '''is_sulphate, unique_dictd = is_sulphate_reducer(non_hyd_output)
    is_methano, unique_dictm = is_methanogen(non_hyd_output)
    is_aceto, unique_dicta = is_acetogen(non_hyd_output)
    is_cyd, unique_dictcyd = is_cydA(non_hyd_output)
    is_nitrate, unique_dictn = is_nitrate_reducer(non_hyd_output)
    is_fumarate, unique_dictf = is_fumarate_reducer(non_hyd_output)
    is_dms, unique_dictdm = is_dms_reducer(non_hyd_output)
    is_pfor, unique_dictp = is_pfor(non_hyd_output)
    is_fhl, unique_dictfhl = is_fhl(non_hyd_output)'''
    #check if any file in input dir is not in results, if not, add it to results, add an empty row
    for file in os.listdir(input_dir):
        #skip if file is not fna, faa, fasta
        if not any([file.endswith(ext) for ext in args.file_extensions]):
            continue
        #res_row1 	hydrogen_consumer	sulphate_reducer	methanogen	acetogen	cydA	nitrate_reducer	fumarate_reducer	DMS_reducer	Pyruvate oxidoreductase	Formate hydrogen lyase	consuming hydrogenase	bifurcating_or_bidirectional_hydrogenase	dsr_matches	mcrA_matches	acsB_matches	cydA_matches	nitrate reductase matches	fumarate reductase matches	DMS_reductase matches	pfor matches	fhl matches	nitrogenase matches
        RESULTS_0= ["filename", "hydrogen_consumer", "sulphate_reducer", "methanogen", 
                    "acetogen","cydA","nitrate_reducer","fumarate_reducer","DMS_reducer","Pyruvate oxidoreductase","Formate hydrogen lyase",
                    "consuming hydrogenase","bifurcating_or_bidirectional_hydrogenase",
                    "protein_match","mcrA_matches","acsB_matches","cydA_matches","nitrate reductase matches","fumarate reductase matches",
                    "DMS reductase matches","pfor matches","fhl matches","nitrogenase matches"]
        
        
        RESULTS_0=["string"]+[True]*12+["string"]*12
        #['Bacteroides_caccae.fna', False, False, False, True, True, True, False, False, True, False, True, False, '{{tr|A0A174III8|A0A174III8_9BACE}}']
        if not any([file.startswith(row[0]) for row in results]):
            lenro=len(RESULTS_0)
            empty_row=[file]+[""]*(lenro-1)
            for i in range(1,lenro):
                #if its a boolean, then set it to false
                if RESULTS_0[i] in [True,False]:
                    empty_row[i]=False
                else:
                    empty_row[i]=""
            print("Adding empty row for ",file,empty_row)
            results.append(empty_row)
                
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            str1 = f.read()
            # if the first row is a header, then set header=True
            if str1.startswith("filename"):
                header = True
            
    # write the results to the output file
    #if output file exists, read the rows
    #current_rows=[]
    #if os.path.exists(output_file):
    #    with open(output_file,'r') as f:
    #        reader=csv.reader(f,delimiter='\t')
    #        current_rows=[row for row in reader]
    with open(output_file, 'a') as f:
        writer = csv.writer(f, delimiter='\t')
        # if header is not there, write the header
        if not header:
            #filename,is_hydrogen_consumer,is_sulphate,is_methano,is_aceto,is_cydA,is_nitrate,is_fumarate,is_dms,is_pfor,is_fhl,
               #is_consume,is_bifurcating
            #
            writer.writerow(header1)

        for row in results:
            # if the row is not in the current rows, then write it
            #if row not in current_rows:

            writer.writerow(row)
    '''# open the file again and replace the uids with descriptions
    with open(output_file, 'r') as f:
        str1 = f.read()
    # replace the uids with descriptions
    Batch_Fetcher.fetch()
    repl = Batch_Fetcher.replace(str1)
    # write the replaced string to the file
    with open(output_file, 'w') as f:
        f.write(repl)'''
    return results
import extract_gene_neighbourhood
def run_neighbourhood_search(input_dir,args):
    #print(f"run_neighbourhood_search (input_dir={input_dir},args={args})")
    res=extract_gene_neighbourhood.main(input_dir,args)
    return res
    



def generate_dummy_CAT_data(filenames, out_file):

        scores = [0.99, 0.85, 0.72, 0.65]


        taxonomy_levels = ['root', 'cellular organisms', 'Bacteria', 'clade', 'phylum',
                           'class', 'order', 'family', 'genus', 'species']
        names = {
            'root': 'root',
            'cellular organisms': 'cellular organisms',
            'Bacteria': 'Bacteria',
            'clade': ['Terrabacteria group', 'FCB group'],
            'phylum': ['DUM_Firmicutes', 'DUM_Bacteroidetes/Chlorobi group'],
            'class': ['DUM_Clostridia', 'DUM_Bacteroidia'],
            'order': ['DUM_Clostridiales', 'DUM_Bacteroidales'],
            'family': ['DUM_Lachnospiraceae', 'DUM_Bacteroidaceae'],
            'genus': ['Lachnospiraceae bacterium', 'Bacteroides'],
            'species': ['DUM_unclassified Lachnospiraceae', 'DUM_Bacteroides genus']
        }
        #make a new dict with all unknown
        names={
            'root':'root',
            'cellular organisms':'cellular organisms',
            'Bacteria':'Bacteria',
            'clade':['NOT AVAILABLE','NOT AVAILABLE'],
            'phylum':['NOT AVAILABLE','NOT AVAILABLE'],
            'class':['NOT AVAILABLE','NOT AVAILABLE'],
            'order':['NOT AVAILABLE','NOT AVAILABLE'],
            'family':['NOT AVAILABLE','NOT AVAILABLE'],
            'genus':['NOT AVAILABLE','NOT AVAILABLE'],
            'species':['NOT AVAILABLE','NOT AVAILABLE']

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
            #write the header
            f.write("filename\ttaxid\tORFs\tlineage\tlineage_scores\tlineage_names\n")
            #writer=csv.writer(f,delimiter='\t')
            for line in cat_data:
                f.write(line + '\n')
        return cat_data






def make_CAT_table(genomes_dir,taxonomy_table,CAT_table_out_file):

    #map the genomes to the ids in the taxonomy table
    filenames={}
    for file in os.listdir(genomes_dir):
        #the id is just the filename minus the extension
        id=os.path.splitext(file)[0]
        filenames[id]=file

    #open the taxonomy table, and find the column containing d__ 
    with open(taxonomy_table,'r') as f:
        reader=csv.reader(f,delimiter=',')
        header=next(reader)
        #get the index of the column containing d__ and p__
        d_index=None
        first_row=next(reader)
        for i,col in enumerate(first_row):
            if col.startswith("d__"):
                d_index=i
                break
    taxons={}
    with open(taxonomy_table,'r') as f:
        #skip the header
        next(f)
        #the first col is the MAG id
        #the d_index is the taxon id
        for row in csv.reader(f,delimiter=','):
            taxon_row=row[d_index]
            #d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__MGBC100174;s__MGBC100174 sp910576895
            #split into the taxons
            taxon_row=taxon_row.split(";")
            domain= taxon_row[0].split("__")[1]
            phylum= taxon_row[1].split("__")[1]
            class_= taxon_row[2].split("__")[1]
            order= taxon_row[3].split("__")[1]
            family= taxon_row[4].split("__")[1]
            genus= taxon_row[5].split("__")[1]
            species= taxon_row[6].split("__")[1]
            taxons[row[0]]=(domain,phylum,class_,order,family,genus,species)
    #make the fake CAT table, it looks like this
    '''# bin	classification	reason	lineage	lineage scores (f: 0.30)	full lineage names										
MGG09353.fna	taxid assigned	based on 1448/1591 ORFs	1;131567;2;1783272;1239;186801;186802;39779;1898207	1.00;1.00;0.99;0.94;0.93;0.53;0.50;0.46;0.46	root (no rank): 1.00	cellular organisms (no rank): 1.00	Bacteria (superkingdom): 0.99	Terrabacteria group (clade): 0.94	Firmicutes (phylum): 0.93	Clostridia (class): 0.53	Clostridiales (order): 0.50	unclassified Clostridiales (no rank): 0.46	Clostridiales bacterium (species): 0.46		
MGG09621.fna	taxid assigned	based on 1118/1255 ORFs	1;131567;2;1783272;1239;186801;186802;31979;1485;59619;1262819	1.00;0.99;0.99;0.86;0.85;0.77;0.64;0.32;0.32;0.31;0.30	root (no rank): 1.00	cellular organisms (no rank): 0.99	Bacteria (superkingdom): 0.99	Terrabacteria group (clade): 0.86	Firmicutes (phylum): 0.85	Clostridia (class): 0.77	Clostridiales (order): 0.64	Clostridiaceae (family): 0.32	Clostridium (genus): 0.32	environmental samples (no rank): 0.31	Clostridium sp. CAG:557 (species): 0.30
    '''
    header="bin\tclassification\treason\tlineage\tlineage scores\tfull lineage names"
    for bin in taxons:
        #if bin is not in the filenames, then skip
        if bin not in filenames:
            continue
        taxon_str=taxons[bin][0]+"(superkingdom)\t" \
                   +taxons[bin][1]+"(phylum)\t" \
                   +taxons[bin][2]+"(class)\t" \
                   +taxons[bin][3]+"(order)\t" \
                   +taxons[bin][4]+"(family)\t" \
                   +taxons[bin][5]+"(genus)\t" \
                   +taxons[bin][6]+"(species)"
        taxon_row={"bin": filenames[bin], "classification": "taxid assigned", "reason": f"based on {random.randint(1000, 6000)}/{random.randint(6000, 10000)} ORFs",
                   "lineage": "1;131567;2;"+";".join([str(random.randint(100000, 999999)) for _ in range(7)]),
                   "lineage scores": ";".join([f"{random.uniform(0.3, 1):.2f}" for _ in range(9)]),
                   "full lineage names": taxon_str}
        taxon_row_str="\t".join([taxon_row[col] for col in header.split("\t")])
        with open(CAT_table_out_file,'a') as f:
            f.write(taxon_row_str+'\n')
    return CAT_table_out_file, taxons


        







def run_CAT1(input_dir, tmp_folder, resume_mode=True, dummy_data=False, CAT_executable_path=None,diamond_exe_path=None,
            CAT_database_path=None,CAT_taxonomy_path=None,taxonomy_table=None):
    '''fills in the taxon information with the taxonomy table if the MAG is in the table, else runs CAT on the MAGs'''
    #if the taxonomy_table is not None, then make the CAT table
    CAT_table_out_file=None
    taxons={}
    
    #if CAT executable path is None, set it to the default
    if CAT_executable_path is None:
        CAT_executable_path=default_cat_path
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
    if taxonomy_table is not None:
        CAT_table_out_file, taxons=make_CAT_table(input_dir,taxonomy_table,os.path.join(tmp_folder,"CAT","out.BAT.bin2classification.named2.txt"))
        for filet in taxons:
            print("taxon alrady in table:",filet)
    # copy the input files to the input folder
    for file in os.listdir(input_dir):
        #only copy the fna, fa, fas, fasta and faa files
        ext=os.path.splitext(file)[1]
        filename=os.path.splitext(file)[0]
        #if filename exists in the taxons, then skip
        
        if ext in [".fna",".fa",".fas",".fasta",".faa"]:
            if filename in taxons:
                print("Filename exists in taxons, skipping:",filename)
                continue
            shutil.copy(os.path.join(input_dir, file), os.path.join(cat_folder, "input"))
    #deduplicate the headers
    deduplicate_names(os.path.join(cat_folder, "input"))
    first=False
    #if there are no files in the input folder, then skip
    output_file=os.path.join(cat_folder, "out.CAT.contig2classification.named.txt")
    if len(os.listdir(os.path.join(cat_folder, "input"))) >0:
        # if dummy data is False
        if not dummy_data:
            # change dir to the CAT folder

            os.chdir(cat_folder)
            # run the CAT command
            # CAT bins -b {bin folder} -d {database folder} -t {taxonomy folder}
            classification_file = os.path.join(cat_folder, "out.BAT.bin2classification.txt")
            #if classification file exists, then dont run CAT
            if not os.path.exists(classification_file):
                command = f"{CAT_executable_path} bins -b {os.path.join(cat_folder, 'input')} -d {CAT_database_path} -t {CAT_taxonomy_path} --force --path_to_diamond {diamond_exe_path}"
                print("COMMAND:",command)
                os.system(command)
            else:
                print("CAT output file exists, skipping CAT bins")
            # add the names: CAT add_names -i {ORF2LCA / classification file} -o {output file} -t {taxonomy folder}
            #if output file exists, then dont run CAT
            output_file = os.path.join(cat_folder, "out.BAT.bin2classification.named.txt")
            if not os.path.exists(output_file):

                command = f"{CAT_executable_path} add_names -i {classification_file} -o {output_file} -t {CAT_taxonomy_path}"
                print("COMMAND:",command)
                os.system(command)
            else:
                print("Named Output file exists, skipping CAT add names")
            # cd back to the original wd
            os.chdir(current_wd)
        # if dummy data is True
        else:
            # create a dummy classification file
            dummy_classification_file = os.path.join(cat_folder, "out.CAT.contig2classification.named.txt")
            dummy_classification_data = generate_dummy_CAT_data(os.listdir(os.path.join(cat_folder, 'input')), dummy_classification_file)
            output_file= os.path.join(cat_folder, "out.CAT.contig2classification.named.txt")
            if CAT_table_out_file is None:
                CAT_table_out_file=output_file
                first=True
        # return the classification file
    #if both the CAT output and the table output exist, concatenate them by appending the table output to the CAT output
    if os.path.exists(output_file) and CAT_table_out_file!=None and os.path.exists(CAT_table_out_file) and not first:
        with open(output_file,'a') as f:
            with open(CAT_table_out_file,'r') as f1:
                f.write(f1.read())
        output_file=output_file
    #else if only the table output exists, then set the output file to the table output
    elif CAT_table_out_file!=None and os.path.exists(CAT_table_out_file):
        #rename the table output to the output file
        os.rename(CAT_table_out_file,output_file)
        output_file=output_file
    #if none exist,raise an error
    else:
        raise FileNotFoundError("No CAT output file found")
    return output_file

def run_CAT(input_dir, tmp_folder, resume_mode=True, dummy_data=False, CAT_executable_path=None,diamond_exe_path=None,
            CAT_database_path=None,CAT_taxonomy_path=None,taxonomy_table=None):
    print(f"run_CAT (input_dir={input_dir},tmp_folder={tmp_folder},resume_mode={resume_mode},dummy_data={dummy_data})")
    '''fills in the taxon information with the taxonomy table if the MAG is in the table, else runs CAT on the MAGs'''
    if taxonomy_table is None:
        taxonomy_table_rows={}
    else:
        taxonomy_table_rows=parse_classification(taxonomy_table)
        for key in taxonomy_table_rows:
            print("MAG in taxon table:",key)
    to_run=[]
    for file in os.listdir(input_dir):
        #print("CHECK FILE:",file)
        #check if file is a fasta file
        if not file.endswith(".fna") and not file.endswith(".fa") and not file.endswith(".fas") and not file.endswith(".fasta") and not file.endswith(".faa"):
            continue
        #get the MAG id
        mag_id=file#os.path.splitext(file)[0]
        #if inperfect_d\+ is in the mag id, then replace it with _
        regex=re.compile(r"_inperfect_\d+")
        mag_id=regex.sub("",mag_id)
        #if the mag id is in the taxonomy table, then skip
        print("MAG ID:",mag_id)
        if mag_id in taxonomy_table_rows:
            print("MAG in taxon table, skipping:",mag_id)
            continue
        to_run.append(file)
    print("to run:",to_run)
    #if there are no files to run, then return the taxonomy table
    if len(to_run)==0:
        return taxonomy_table_rows
    #run CAT on the files
    #copy the files to the CAT folder
    cat_folder = os.path.join(tmp_folder, "CAT")
    #if CAT folder does not exist, create it
    if not os.path.exists(cat_folder):
        os.makedirs(cat_folder)
    #if input folder not exists, create it
    if not os.path.exists(os.path.join(cat_folder, "input")):
        os.makedirs(os.path.join(cat_folder, "input"))
    for file in to_run:
        shutil.copy(os.path.join(input_dir, file), os.path.join(cat_folder, "input"))
    current_wd = os.getcwd()
    deduplicate_names(os.path.join(cat_folder, "input"))
    #if there are no files in the input folder, then skip
    output_file=os.path.join(cat_folder, "out.CAT.contig2classification.named.txt")
    if len(os.listdir(os.path.join(cat_folder, "input"))) >0:
        # if dummy data is False
        #if the cat path is not exist, set dummy data to True
        if CAT_executable_path is None:
            dummy_data=True
        if not dummy_data:
        
            # change dir to the CAT folder

            os.chdir(cat_folder)
            # run the CAT command
            # CAT bins -b {bin folder} -d {database folder} -t {taxonomy folder}
            classification_file = os.path.join(cat_folder, "out.BAT.bin2classification.txt")
            #if classification file exists, then dont run CAT
            if not os.path.exists(classification_file):
                command = f"{CAT_executable_path} bins -b {os.path.join(cat_folder, 'input')} -d {CAT_database_path} -t {CAT_taxonomy_path} --force --path_to_diamond {diamond_exe_path}"
                print("COMMAND:",command)
                os.system(command)
            else:
                print("CAT output file exists, skipping CAT bins")
            # add the names: CAT add_names -i {ORF2LCA / classification file} -o {output file} -t {taxonomy folder}
            #if output file exists, then dont run CAT
            output_file = os.path.join(cat_folder, "out.BAT.bin2classification.named.txt")
            if not os.path.exists(output_file):

                command = f"{CAT_executable_path} add_names -i {classification_file} -o {output_file} -t {CAT_taxonomy_path}"
                print("COMMAND:",command)
                os.system(command)
            else:
                print("Named Output file exists, skipping CAT add names")
            # cd back to the original wd
            os.chdir(current_wd)
        else:
            # create a dummy classification file
            dummy_classification_file = os.path.join(cat_folder, "out.CAT.contig2classification.named.txt")
            dummy_classification_data = generate_dummy_CAT_data(os.listdir(os.path.join(cat_folder, 'input')), dummy_classification_file)
            output_file= os.path.join(cat_folder, "out.CAT.contig2classification.named.txt")
    #parse the output file
    if output_file is not None and os.path.exists(output_file):
        taxons=parse_classification_file(output_file)
    else:
        taxons={}
    #update the taxonomy table with the new taxons
    taxonomy_table_rows.update(taxons)
    return taxonomy_table_rows





def parse_classification(taxon_output):
    print("TAXON OUTPUT:",taxon_output)
    tax_dict={}
    with open(taxon_output,'r') as f:
         if "\t" in f.read():
            delimiter='\t'
         else:
            delimiter=','
    with open(taxon_output,'r') as f:
        reader=csv.reader(f,delimiter=delimiter)
        header=next(reader)
        #get the index of the column containing d__ and p__
        d_index=None
        first_row=next(reader)
        print("FIRST ROW:",first_row)
        for i,col in enumerate(first_row):
            if col.startswith("d__"):
                d_index=i
                break
    
    with open(taxon_output,'r') as f:
        #if there are tabs, then use the tab delimiter
       
        reader=csv.reader(f,delimiter=delimiter)
        header=next(reader)
        #get the taxon dict
        for row in reader:
            #print("tax row:",row)
            tax_data=row[d_index].split(";")
            domain=None
            #find the ele in tax_data that starts with d__
            for ele in tax_data:
                if ele.startswith("d__"):
                    domain=ele.split("__")[1]
                    break
            #find the ele in tax_data that starts with p__
            phylum=None
            for ele in tax_data:
                if ele.startswith("p__"):
                    phylum=ele.split("__")[1]
                    break
            #find the ele in tax_data that starts with c__
            class_=None
            for ele in tax_data:
                if ele.startswith("c__"):
                    class_=ele.split("__")[1]
                    break
            #find the ele in tax_data that starts with o__
            order=None
            for ele in tax_data:
                if ele.startswith("o__"):
                    order=ele.split("__")[1]
                    break
            #find the ele in tax_data that starts with f__
            family=None
            for ele in tax_data:
                if ele.startswith("f__"):
                    family=ele.split("__")[1]
                    break
            #find the ele in tax_data that starts with g__
            genus=None
            for ele in tax_data:
                if ele.startswith("g__"):
                    genus=ele.split("__")[1]
                    break
            #find the ele in tax_data that starts with s__
            species=None
            for ele in tax_data:
                if ele.startswith("s__"):
                    species=ele.split("__")[1]
                    break
            #print("taxons:",domain,phylum,class_,order,family,genus,species)
            taxons={"domain":domain,"phylum":phylum,"class":class_,"order":order,"family":family,"genus":genus,"species":species}
            tax_dict[row[0]]=taxons
    return tax_dict
def parse_classification_file(taxon_output):
    taxon_dict={  }
    with open(taxon_output,'r') as f:
            reader=csv.reader(f,delimiter='\t')
            #skip the header
            next(reader)

            for row in reader:
                #if len of row is 0, skip
                if len(row)==0:
                    continue

                #the taxon data is in the 6th column onwards
                taxon_dict1={}
                for col in range(0,len(row)):
                    #get the taxon level, its in the brackets
                    #use re to find the words in brackets
                    taxon_level=re.findall(r'\((.*?)\)',row[col])
                    #if there is a taxon level
                    if taxon_level:
                        #get the taxon name
                        taxon_name=row[col].split(':')[0]
                        #remove the words in the brackets using re
                        taxon_name=re.sub(r'\(.*?\)','',taxon_name)
                        #strip
                        taxon_name=taxon_name.strip()
                        #strip the taxon level
                        taxon_level=taxon_level[0].strip()
                        #if taxon level is superkingdom, taxon dict1 kingdom is superkingdom
                        if taxon_level=='superkingdom':
                            taxon_dict1['kingdom']=taxon_name
                        #if taxon level is phylum, taxon dict1 phylum is phylum
                        elif taxon_level=='phylum':
                            taxon_dict1['phylum']=taxon_name
                        #if taxon level is class, taxon dict1 class is class
                        elif taxon_level=='class':
                            taxon_dict1['class']=taxon_name
                        #if taxon level is order, taxon dict1 order is order
                        elif taxon_level=='order':
                            taxon_dict1['order']=taxon_name
                        #if taxon level is family, taxon dict1 family is family
                        elif taxon_level=='family':
                            taxon_dict1['family']=taxon_name
                        #if taxon level is genus, taxon dict1 genus is genus
                        elif taxon_level=='genus':
                            taxon_dict1['genus']=taxon_name
                        #if taxon level is species, taxon dict1 species is species
                        elif taxon_level=='species':
                            taxon_dict1['species']=taxon_name
                #add the taxon dict1 to the taxon dict
                taxon_dict[row[0]]=taxon_dict1
    return taxon_dict
import get_protein_names
def create_feature_table1 (out_file,hydrogen_classes_activity_file):
    #read the out_file using pandas
    df=pd.read_csv(out_file,delimiter='\t')
    #read the hydrogen_classes_activity_file using pandas
    df1=pd.read_csv(hydrogen_classes_activity_file,delimiter=',')
    #FeFe-WP_005679640.1_-_Bacteroides_caccae_-_[FeFe]_Group_B
    group_regex=r'(\[.+\])_Group_(.{2})'
    #get the groups
    #compile the regex
    group_regex=re.compile(group_regex)
    def apply_row(row):
        #if the row["protein_match"] is not nan
        if not pd.isna(row["protein_match"]):
            #get the group
            group=group_regex.search(row["protein_match"])
            #if group is not None
            if group is not None:
                #get the group
                class1=group.group(1).strip()
                #remove the []
                class1=class1.replace("[","").replace("]","")
                group1=group.group(2).strip()
                #find the rows in the df1 where Class is equal to class1 and Group is equal to group1
                rows=df1[(df1["Class"]==class1) & (df1["Group"]==group1)]
                #if there are rows
                if len(rows)>0:
                    #get the first row
                    hyd_direction=rows.iloc[0]["direction"].strip()
                    hyd_localization=rows.iloc[0]["localization"].strip()
                    oyxgen_tolerance=rows.iloc[0]["oxygen tolerance"].strip()
                    return pd.Series([hyd_direction,hyd_localization,oyxgen_tolerance])
        return pd.Series([None,None,None])
    if 'oxygen_tolerance' not in df.columns:
        df['oxygen_tolerance'] = None
    if 'hyd_direction' not in df.columns:
        df['hyd_direction'] = None
    if 'hyd_localization' not in df.columns:
        df['hyd_localization'] = None
    df[["hyd_direction","hyd_localization","oxygen_tolerance"]] = df.apply(apply_row, axis=1)
    #group by query and count the number of rows with hyd_direction =H2-uptake 
    queries=df[df["hyd_direction"]=="H2-uptake"]["query"].unique()
    #set all rows where query is in queries to H2-uptake
    df["contains_uptake_hydrogenase"]=False
    df.loc[df["query"].isin(queries),"contains_uptake_hydrogenase"]=True
    #group by query and count the number of rows with hyd_direction =Electron-bifurcation, H2-evolution or electron-bifurcation,  Bidirectional
    queries=df[df["hyd_direction"].isin(["Electron-bifurcation","Bidirectional","H2-evolution or electron-bifurcation"]) ]["query"].unique()
    #set all rows where query is in queries to H2-evolution
    df["contains bidirectional or bifurcating hydrogenase"]=False
    df.loc[df["query"].isin(queries),"contains bidirectional or bifurcating hydrogenase"]=True
    #group  by query and count the number of rows with hyd_direction =H2-evolution or H2-evolution or electron-bifurcation
    queries=df[df["hyd_direction"] .isin (["H2-evolution","H2-evolution or electron-bifurcation"])]["query"].unique()
    #set all rows where query is in queries to H2-evolution
    df["contains_evolution_hydrogenase"]=False
    df.loc[df["query"].isin(queries),"contains_evolution_hydrogenase"]=True
    mapping = {
    "FeFe A1": "H2-evolution",
    "FeFe A2": "H2-uptake",
    "FeFe A3": "Electron-bifurcation",
    "FeFe A4": "H2-evolution or electron-bifurcation",
    "FeFe B": "H2-evolution",
    "FeFe C1": "H2-sensing",
    "FeFe C2": "H2-sensing",
    "FeFe C3": "H2-sensing",
    "Fe Hmd": "Bidirectional",
    "NiFe 1a": "H2-uptake",
    "NiFeSe 1a": "H2-uptake",
    "NiFe 1b": "H2-uptake",
    "NiFe 1c": "H2-uptake",
    "NiFe 1d": "H2-uptake",
    "NiFe 1e": "Bidirectional",
    "NiFe 1f": "H2-uptake",
    "NiFe 1g": "H2-uptake",
    "NiFe 1h": "H2-uptake",
    "NiFe 1i": "H2-uptake",
    "NiFe 1j": "H2-uptake",
    "NiFe 1k": "H2-uptake",
    "NiFe 2a": "H2-uptake",
    "NiFe 2b": "H2-sensing",
    "NiFe 2c": "H2-sensing",
    "NiFe 2d": "H2-uptake",
    "NiFe 2e": "H2-uptake",
    "NiFe 3a": "Bidirectional",
    "NiFe 3b": "Bidirectional",
    "NiFe 3c": "Electron-bifurcation",
    "NiFe 3d": "Bidirectional",
    "NiFe 4a": "H2-evolution",
    "NiFe 4b": "H2-evolution",
    "NiFe 4c": "H2-evolution",
    "NiFe 4d": "H2-evolution",
    "NiFe 4e": "Bidirectional",
    "NiFe 4f": "H2-evolution",
    "NiFe 4g": "H2-evolution",
    "NiFe 4h": "H2-uptake",
    "NiFe 4i": "H2-uptake"
}
    def extract_hydrogenase_group(row):
        #the group is in the protein match, after the word "Group"
        try:
            spl=row["protein_match"].split("Group")
            #get the group number
            group_number=spl[1].strip()
            #get the class
            if "FeFe" in row["protein_match"]:
                hyd_class="FeFe"
            elif "NiFe" in row["protein_match"]:
                hyd_class="NiFe"
            elif "Fe" in row["protein_match"]:
                hyd_class="Fe"
            #strip and rmeove underscores
            hyd_class=hyd_class.strip().replace("_","")
            group_number=group_number.strip().replace("_","")
            return hyd_class+" "+group_number
        except Exception as e:
            #print("ERROR:",row["protein_match"],str(e))
            return "NA"
    #apply the function
    #df_hydd=df.copy()
    #print tyhe columns of df_hydd
    #print("COLUMNS:",df_hydd.columns )
    #df_hydd["hyd_group"]=df_hydd.apply(extract_hydrogenase_group, axis=1)
    #df_hydd["hyd_direction"]=df_hydd["hyd_group"].map(mapping)
    #drop all columns except the genome name , protein match, hyd group and hyd direction
    #df_hydd=df_hydd[["genome_name","protein_match","hyd_group","hyd_direction"]]

    #group by genome name 
    #df_hydd_grouped=df_hydd.groupby('genome_name').agg(lambda x: x.tolist())
    #df_hydd_grouped.reset_index(inplace=True)
    #check if the hyd_direction contains H2-evolution, Bidirectional , Electron-bifurcation,or H2-evolution or electron-bifurcation, set the is_hydrogen_producer to True
    #renmame the query column to genome_name
    df.rename(columns={"query":"genome_name"},inplace=True)
    for i,row in df.iterrows():
        genome_name=row["genome_name"]
        df.loc[df["genome_name"]==genome_name,"is_hydrogen_producer"]=False
        df.loc[df["genome_name"]==genome_name,"is_hydrogen_consumer"]=False
    #check if the hyd_direction contains H2-evolution, Bidirectional , Electron-bifurcation,or H2-evolution or electron-bifurcation, set the is_hydrogen_producer to True
    for i,row in df.iterrows():
        #skip NA 
        if pd.isna(row["hyd_direction"]):
            continue
        if "H2-evolution" in row["hyd_direction"] or "Bidirectional" in row["hyd_direction"] or "Electron-bifurcation" in row["hyd_direction"] or "H2-evolution or electron-bifurcation" in row["hyd_direction"]:
            genome_name=row["genome_name"]
            df.loc[df["genome_name"]==genome_name,"is_hydrogen_producer"]=True
        
    #check if the hyd_direction contains H2-uptake, Bidirectional , Electron-bifurcation,or H2-evolution or electron-bifurcation, set the is_hydrogen_consumer to True
    for i,row in df.iterrows():
        if pd.isna(row["hyd_direction"]):
            continue
        if "H2-uptake" in row["hyd_direction"] or "Bidirectional" in row["hyd_direction"] or "Electron-bifurcation" in row["hyd_direction"] or "H2-evolution or electron-bifurcation" in row["hyd_direction"]:
            genome_name=row["genome_name"]
            df.loc[df["genome_name"]==genome_name,"is_hydrogen_consumer"]=True

    df["hyd_group"]=df.apply(extract_hydrogenase_group, axis=1)

    #rename the query to genome_name, methanogen to is_methanogen, acetogen to is_acetogen, sulfate reducer to is_sulfate_reducer, cydA to aerobic_respirer
    df.rename(columns={"methanogen":"is_methanogen","acetogen":"is_acetogen","sulfate reducer":"is_sulfate_reducer","cydA":"aerobic_respirer"},inplace=True)
    #save it to a csv file
    df.to_csv(out_file,index=False)
    get_protein_names.process(out_file)




def create_feature_table (tmp_filename,hydrogen_classes_activity_file,hydrogen_consumer_pipeline_output_folder,out_file, taxon_table=None):
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
    if taxon_table is None:
        combined_classification_file = os.path.join(hydrogen_consumer_pipeline_output_folder,
                                                    'combined_classification_file.txt')
        #remove folders that contain "tmp"
        folders=[folder for folder in folders if "tmp" not in folder]
        # combine the classification files
    
        with open(combined_classification_file, 'w') as f:

            for folder in folders:
                header=None
                # get the classification file
                #TODO: Reset the file back to out.BAT.bin2classification.named.txt after debugging
                classification_file = os.path.join(hydrogen_consumer_pipeline_output_folder, folder, tmp_filename,"CAT",
                                                'out.CAT.contig2classification.txt')
                #check if the classification file exists, if not, then use out.BAT.bin2classification.named.txt
                if not os.path.exists(classification_file):
                    classification_file=os.path.join(hydrogen_consumer_pipeline_output_folder, folder, tmp_filename,"CAT",
                                                'out.BAT.bin2classification.named.txt')
                #'out.BAT.bin2classification.named.txt')
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
    else:
        taxon_output=taxon_table

    #IF TAXON TABLE IS NOT none
    if taxon_table is not None:
        tax_dict=parse_classification(taxon_table)
    else:
        tax_dict=parse_classification_file(taxon_output)
    
    

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
        print("Batch root dir is not empty, resuming")
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
            #rename the file to .fna (THIS IS UNTESTED)
            #remove the extension, 
            new_name=os.path.splitext(os.path.basename(file))[0]
            new_name= new_name+'.fna'
            os.rename(os.path.join(batch_dir,os.path.basename(file)),os.path.join(batch_dir,new_name))
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
def run_batches(batch_root_dir,f_output_file,tmp_dir,args):
    print(f"run_batches(batch_root_dir={batch_root_dir},output_file={f_output_file},tmp_dir={tmp_dir},args={args})")
    #get all the batch dirs
    batch_dirs=[os.path.join(batch_root_dir,batch) for batch in os.listdir(batch_root_dir)]
    #only get the folders
    batch_dirs=[batch for batch in batch_dirs if os.path.isdir(batch)]
    #remove the dirs that dont start woth batch
    batch_dirs=[batch for batch in batch_dirs if os.path.basename(batch).startswith("batch")]
    #sort the batch dirs
    batch_dirs.sort()
    #if not in resume mode, remove the files in root dir matching batch_xresults.tsv, combined_classification_file.txt, and all_output.csv
    if not args.resume:
        #remove the files
        for file in os.listdir(batch_root_dir):
            if file.endswith("results.tsv") or file=="combined_classification_file.txt" or file=="all_output.csv":
                os.remove(os.path.join(batch_root_dir,file))
    print("batch dir run",batch_root_dir,batch_dirs)
    for batch_dir in batch_dirs:

        print("BATCH_DIR",batch_dir)
        #for each batch dir, run the pipeline
        #extract the filename of the output file
        output_file=os.path.basename(f_output_file)
        #add the batch folder name to the output file
        #get the parent dir of the batch dir
        parent_dir=os.path.dirname(batch_dir)

        output_file=os.path.join(parent_dir,batch_dir+"results.tsv")
        #get the foldername of the tmp dir
        tmp_dir=os.path.basename(tmp_dir)
        #add the batch folder name to the tmp dir
        tmp_dir1=os.path.join(batch_dir,tmp_dir)
        results=batch_run2(batch_dir,output_file,tmp_dir1,args)
        #make the results into a dataframe
        #"query":query,"sulphate_reducer":is_sulphate, "methanogen":is_methano, "acetogen":is_aceto,
        #                                   "cydA":is_cyd, "nitrate_reducer":is_nitrate, "fumarate_reducer":is_fumarate, "dms_reducer":is_dms,
        #                                  "pfor":is_pfor, "fhl":is_fhl, "hydrogen_consumer":is_hydrogen_consumer,"hydrogen_producer":is_producing}
        #header=['query', 'sulphate_reducer', 'methanogen', 'acetogen', 'cydA', 'nitrate_reducer', 'fumarate_reducer', 'dms_reducer', 'pfor', 'fhl', 'hydrogen_consumer', 'hydrogen_producer', 'protein_match']
        #header=['query', 'sulphate_reducer', 'methanogen', 'acetogen', 'cydA', 'nitrate_reducer', 'fumarate_reducer', 'dms_reducer', 'pfor', 'fhl', 'hydrogen_consumer', 'hydrogen_producer', 'nitrogenase', 'protein_match']
        header= ['query', 'sulphate_reducer', 'methanogen', 'acetogen', 'cydA', 'nitrate_reducer', 'fumarate_reducer', 'dms_reducer', 'pfor', 'fhl', 
                 'hydrogen_consumer', 'hydrogen_producer', 'nitrogenase', 'protein_match',
                   'sulphate_reducer_matches', 'methanogen_matches', 'acetogen_matches', 'cydA_matches', 'nitrate_reducer_matches', 
                   'fumarate_reducer_matches', 'dms_reducer_matches', 'pfor_matches', 'fhl_matches', 'hydrogenase_matches', 'nitrogenase_matches']
        #results_transposed = list(map(list, zip(*results)))
        results=pd.DataFrame(results,columns=header)
        #set the index to the query
        results.set_index("query",inplace=True)
        #generate the dummy data
        #create_dummy_results(batch_dir,output_file)
        #run the CAT
        #TODO: Set the dummy data to False in Production
        print("SINGLE BATCH CAT:")
        #if args contains taxonomy_table_path, then dont do the CAT table
        #if args.taxonomy_table_path is  None:
        taxon_dict=run_CAT(batch_dir, tmp_dir1, resume_mode=args.resume, dummy_data=True, CAT_executable_path=args.cat_path,diamond_exe_path=args.diamond_exe_path,
                    CAT_database_path=args.cat_database_path,CAT_taxonomy_path=args.cat_taxonomy_path,taxonomy_table=args.taxonomy_table_path)
        print("END SINGLE BATCH CAT")
        for index, row in results.iterrows():
            #get the taxon data
            if index in taxon_dict:
                taxon_data=taxon_dict[index]
                #update the row with the taxon data
                for key in taxon_data:
                    results.loc[index,key]=taxon_data[key]
            else:
                #put unknown for the taxon data
                for key in taxon_data:
                    results.loc[index,key]="Unknown"
        #write the results to the output file
        print("final results written to:",output_file)
        results.to_csv(output_file,sep="\t")
    #combine all the batch files into the final output file
    #get all the batch files
    all_rows=[]
    for i,batch_dir in enumerate(batch_dirs):
        batch_results_file=os.path.join(batch_root_dir,f"batch_{i}results.tsv")
        print("reading file:",batch_results_file)
        #read the file
        with open(batch_results_file,'r') as f:
            reader=csv.reader(f,delimiter='\t')
            header=next(reader)
            for row in reader:
                all_rows.append(row)
    #write the rows to the output file
    with open(f_output_file,'w') as f:
        writer=csv.writer(f,delimiter='\t')
        writer.writerow(header)
        for row in all_rows:
            writer.writerow(row)
    
    return f_output_file



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
    #if not in resume mode, remove the files in root dir matching batch_xresults.tsv, combined_classification_file.txt, and all_output.csv
    if not args.resume:
        #remove the files
        for file in os.listdir(batch_root_dir):
            if file.endswith("results.tsv") or file=="combined_classification_file.txt" or file=="all_output.csv":
                os.remove(os.path.join(batch_root_dir,file))
    if max_batches is None or max_batches > len(batch_dirs):
        max_batches = len(batch_dirs)

    semaphore = threading.Semaphore(max_batches)

    def process_batch(batch_dir):
        try:
            # For each batch dir, run the pipeline
            # Extract the filename of the output file
            #output_filename = os.path.basename(output_file)
            parent_dir = os.path.dirname(batch_dir)

            output_file1 = os.path.join(parent_dir, batch_dir + "results.tsv")
            # Add the batch folder name to the output file
            #output_filename = os.path.join(batch_dir, output_filename)
            # Get the foldername of the tmp dir
            tmp_dirname = os.path.basename(tmp_dir)
            # Add the batch folder name to the tmp dir
            tmp_dir1 = os.path.join(batch_dir, tmp_dirname)
            #taxon_dict=run_CAT(batch_dir, tmp_dir1, resume_mode=args.resume, dummy_data=False, CAT_executable_path=args.cat_path,diamond_exe_path=args.diamond_exe_path,
            #            CAT_database_path=args.cat_database_path,CAT_taxonomy_path=args.cat_taxonomy_path,taxonomy_table=args.taxonomy_table_path)
            
            # Run the pipeline
            print("BATCH RUN:")
            results=batch_run2(batch_dir, output_file1, tmp_dir1, args)
            #header=['query', 'sulphate_reducer', 'methanogen', 'acetogen', 'cydA', 'nitrate_reducer', 'fumarate_reducer', 'dms_reducer', 'pfor', 'fhl', 'hydrogen_consumer', 'hydrogen_producer', 'nitrogenase', 'protein_match']
            header= ['query', 'sulphate_reducer', 'methanogen', 'acetogen', 'cydA', 'nitrate_reducer', 'fumarate_reducer', 'dms_reducer', 'pfor', 'fhl', 
                 'hydrogen_consumer', 'hydrogen_producer', 'nitrogenase', 'protein_match',
                   'sulphate_reducer_matches', 'methanogen_matches', 'acetogen_matches', 'cydA_matches', 'nitrate_reducer_matches', 
                   'fumarate_reducer_matches', 'dms_reducer_matches', 'pfor_matches', 'fhl_matches', 'hydrogenase_matches', 'nitrogenase_matches']
            #results_transposed = list(map(list, zip(*results)))
            #print all the results
            for row in results:
                print(len(row),row)
            results=pd.DataFrame(results,columns=header)
            #set the index to the query
            results.set_index("query",inplace=True)
            #generate the dummy data
            #create_dummy_results(batch_dir,output_file)
            #run the CAT
            print("BATCH CAT:")
            #if args contains taxonomy_table_path, then dont do the CAT table
            #if args.taxonomy_table_path is None:
            taxon_dict=run_CAT(batch_dir, tmp_dir1, resume_mode=args.resume, dummy_data=False, CAT_executable_path=args.cat_path,diamond_exe_path=args.diamond_exe_path,
                        CAT_database_path=args.cat_database_path,CAT_taxonomy_path=args.cat_taxonomy_path,taxonomy_table=args.taxonomy_table_path)
            
            for index, row in results.iterrows():
                #get the taxon data
                if index in taxon_dict:
                    taxon_data=taxon_dict[index]
                    #update the row with the taxon data
                    for key in taxon_data:
                        results.loc[index,key]=taxon_data[key]
                else:
                    #put unknown for the taxon data
                    for key in taxon_data:
                        results.loc[index,key]="Unknown"
            #write the results to the output file
            print("results.tocsv:",output_file1)
            results.to_csv(output_file1,sep="\t")
            #
        except Exception as e:
            print("THREAD EXCEPTION")
            import traceback
            #print(e)
            print(traceback.format_exc())
            raise e
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
    #get all the batch files
    print("writing output file:",output_file)
    all_rows=[]
    for i,batch_dir in enumerate(batch_dirs):
        batch_results_file=os.path.join(batch_root_dir,f"{batch_dir}results.tsv")
        #read the file
        with open(batch_results_file,'r') as f:
            reader=csv.reader(f,delimiter='\t')
            header=next(reader)
            for row in reader:
                all_rows.append(row)
    #write the rows to the output file
    with open(output_file,'w') as f:
        writer=csv.writer(f,delimiter='\t')
        writer.writerow(header)
        for row in all_rows:
            writer.writerow(row)
    return output_file




def main():
    #get the arguments
    print("RUNNING")
    args=argparse.ArgumentParser()
    args.add_argument("-i","--input_dir",help="the input directory containing the genomes to be analyzed",required=True)
    args.add_argument("-o","--output_file",help="the output file to save the results",required=True)
    args.add_argument("-p","--prodigal_exe_path",help="the path to the prodigal executable",default=default_prodigal_path)
    args.add_argument("-d","--diamond_exe_path",help="the path to the diamond executable",default=default_diamond_exe_path)
    args.add_argument("-hy","--hydrogenase_script_path",help="the path to the hydrogenase script",default=default_hydrogenase_script_path)
    args.add_argument("-dsrA","--dsrA_database_fasta_file",help="the path to the dsr database",default=None)
    #apra
    args.add_argument("-apra","--aprA_database_fasta_file",help="the path to the apra database",default=None)
    #asrA
    args.add_argument("-asrA","--asrA_database_fasta_file",help="the path to the asrA database",default=None)
    #"C:\Users\abel\Documents\hydrogenases\trimethylamine_N_oxide_reduction\DmsA.fasta"
    args.add_argument("-dmsA","--dmsA_database_fasta_file",help="the path to the dmsA database",default=None)
    #"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\napA.fasta"
    args.add_argument("-napA","--napA_database_fasta_file",help="the path to the napA database",default=None)
    #"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\narG.fasta"
    args.add_argument("-narG","--narG_database_fasta_file",help="the path to the narG database",default=None)
    #"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\nrfA.fasta"
    args.add_argument("-nrfA","--nrfA_database_fasta_file",help="the path to the nrfA database",default=None)
    #c:\Users\abel\Documents\hydrogenases\fumarate_reducer\frdA.fasta
    args.add_argument("-frdA","--frdA_database_fasta_file",help="the path to the frdA database",default=None)
    #"C:\Users\abel\Documents\hydrogenases\cydA.fasta"
    args.add_argument("-cydA","--cydA_database_fasta_file",help="the path to the cydA database",default=None)
    #"C:\Users\abel\Documents\hydrogenases\hyd_production\fhl.fasta"
    args.add_argument("-fhl","--fhl_database_fasta_file",help="the path to the fhl database",default=None)
    #"C:\Users\abel\Documents\hydrogenases\hyd_production\PFOR.fasta"
    args.add_argument("-pfor","--pfor_database_fasta_file",help="the path to the pfor database",default=None)
    #"C:\Users\abel\Documents\hydrogenases\nitrogenase\nitrogenase.fasta"
    args.add_argument("-nifH","--nitrogenase_database_fasta_file",help="the path to the nitrogenase database",default=None)

    args.add_argument("-acsb","--acsB_database_fasta_file",help="the path to the acetate kinase database",default=None)
    args
    args.add_argument("-mcrA","--mcrA_database_fasta_file",help="the path to the mcr database",default=None)
    #add the argument to run the neighbourhood search
    args.add_argument("-n","--neighbourhood_search",help="search the neighbourhood of FeFe A hydrogenases to determine the subgroup",action="store_false",default=True)
    #add the argument for the diadb path
    args.add_argument("-diadb","--diadb_path",help="the path to the diadb",default=default_DIADB_dir)
    #add a flag for max threads
    #add an argument for the CAT path
    args.add_argument("-c","--cat_path",help="the path to the CAT executable",default=default_cat_path)
    #add an argument for the max threads
    args.add_argument("-th","--max_threads",help="the maximum number of threads to use",default=1,type=int)
    #add a flag for batch size
    args.add_argument("-b","--batch_size",help="the batch size",default=100,type=int)
    args.add_argument('-t',"--tmp_dir",help="the name of the temporary directory (in the batch files)",default="./tmp")
    #add a flag to run the acetoscan
    args.add_argument("-a","--acetoscan",help="run acetoscan",action="store_true",default=False)
    #add an argument for file extensions to be considered
    args.add_argument("-x","--file_extensions",help="the file extensions to be considered",default=[".fa",".fna",".fasta"],nargs="+")
    #add the path to the hydrogenase activity file
    args.add_argument("-hyd","--hydrogen_classes_activity_file",help="the path to the hydrogen classes activity file",default=DEFAULT_HYDROGENASE_ACTIVITY_FILE)
    #add a flag to enable resume mode, to continue in case of failure,its default is true
    args.add_argument("-r","--resume",help="resume mode, in case the previous run crash, use cached output from previous run",action="store_true",default=False)
    #add the flag for cat database path
    args.add_argument("-catdb","--cat_database_path",help="the path to the CAT database",default=default_CAT_database_path)
    #add the flag for the taxonomy path
    args.add_argument("-tax","--cat_taxonomy_path",help="the path to the CAT taxonomy",default=default_CAT_taxonomy_path)
    #add a flag of the taxonomy table path
    args.add_argument("-taxt","--taxonomy_table_path",help="the path to the taxonomy table",default=None)
    #add a flag for the acetoscan results directory
    args.add_argument("-ar","--acetoscan_results_dir",help="the path to the acetoscan results directory",default=None)
    args=args.parse_args()
    #batch the files
    batch_root_dir=batch_the_files(args.input_dir,os.path.join(args.input_dir,"batches"),batch_size=args.batch_size,resume=args.resume)
    #run the batches
    #run the pipeline
    #batch_run2(args.input_dir,args.output_file,tmp_dir=args.tmp_dir, args=args)
    #run the batches
    #if max threads is 1, use the normal batch run
    if args.max_threads==1:
        print("RUN BATCHES")
        output_file=run_batches(batch_root_dir,args.output_file,args.tmp_dir,args)
    else:
        print("RUN BATCHES MULTITHREADED")
        output_file=run_batches_multithreaded(batch_root_dir,args.output_file,args.tmp_dir,args,max_batches=args.max_threads)
    #run the CAT command
    #tax_output=run_CAT(args.input_dir,args.tmp_dir,resume_mode=args.resume,dummy_data=False,CAT_executable_path=args.cat_path,diamond_exe_path=args.diamond_exe_path,CAT_database_path=args.cat_database_path,CAT_taxonomy_path=args.cat_taxonomy_path)
    # create the feature table
    print("CREATING FEATURE TABLE",output_file)
    create_feature_table1(output_file,args.hydrogen_classes_activity_file)

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


def tese1():
    input_folder=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\test_gene_organization"
    tmp_folder=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\test_gene_organization\tmp"
    output_file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\test_gene_organization\results.csv"
    print("RUNNING")
    args=argparse.ArgumentParser()
    args.add_argument("-i","--input_dir",help="the input directory containing the genomes to be analyzed",default=input_folder,
                      required=False)
    args.add_argument("-o","--output_file",help="the output file to save the results",required=False,default=output_file)
    args.add_argument("-p","--prodigal_exe_path",help="the path to the prodigal executable",default=default_prodigal_path)
    args.add_argument("-d","--diamond_exe_path",help="the path to the diamond executable",default=default_diamond_exe_path)
    args.add_argument("-hy","--hydrogenase_script_path",help="the path to the hydrogenase script",default=default_hydrogenase_script_path)
    args.add_argument("-dsr","--dsr_database_path",help="the path to the dsr database",default=dsr_database_fasta_file)
    args.add_argument("-mcr","--mcr_database_path",help="the path to the mcr database",default=mcr_database_fasta_file)
    #add a flag for max threads
    args.add_argument("-th","--max_threads",help="the maximum number of threads to use",default=1,type=int)
    #add a flag for batch size
    args.add_argument("-b","--batch_size",help="the batch size",default=100,type=int)
    args.add_argument('-t',"--tmp_dir",help="the path to the temporary directory",default=tmp_folder)
    #add a flag to run the acetoscan
    args.add_argument("-a","--acetoscan",help="run acetoscan",action="store_true",default=False)
    #add a flag for the acetoscan results directory
    args.add_argument("-ar","--acetoscan_results_dir",help="the path to the acetoscan results directory",default=None)
    #add an argument for file extensions to be considered
    args.add_argument("-x","--file_extensions",help="the file extensions to be considered",default=[".fa",".fna",".fasta"],nargs="+")
    #add the path to the hydrogenase activity file
    args.add_argument("-hyd","--hydrogen_classes_activity_file",help="the path to the hydrogen classes activity file",default=DEFAULT_HYDROGENASE_ACTIVITY_FILE)
    #add a flag to enable resume mode, to continue in case of failure,its default is true
    args.add_argument("-r","--resume",help="resume mode, in case the previous run crash, use cached output from previous run",action="store_true",default=True)
    args=args.parse_args()
    run_batches(input_folder,output_file,tmp_folder,args)

def tese3():
    input_folder='/tllhome/abel/hydrogen_consumer_pipeline/test/test_genomes'
    tmp_folder='/tllhome/abel/hydrogen_consumer_pipeline/test/test_genomes/tmp'
    output_file='/tllhome/abel/hydrogen_consumer_pipeline/test/test_genomes/results.csv'
    #r"C:\Users\abel\Documents\hydrogenases\hydrogenase_subunits\test_genomes\genomes"
    #tmp_folder=r"C:\Users\abel\Documents\hydrogenases\hydrogenase_subunits\test_genomes\genomes\tmp"
    #output_file=r"C:\Users\abel\Documents\hydrogenases\hydrogenase_subunits\test_genomes\genomes\results.csv"
    print("RUNNING")
    args=argparse.ArgumentParser()
    args.add_argument("-i","--input_dir",help="the input directory containing the genomes to be analyzed",default=input_folder,
                      required=False)
    args.add_argument("-o","--output_file",help="the output file to save the results",required=False,default=output_file)
    args.add_argument("-p","--prodigal_exe_path",help="the path to the prodigal executable",default=default_prodigal_path)
    args.add_argument("-d","--diamond_exe_path",help="the path to the diamond executable",default=default_diamond_exe_path)
    args.add_argument("-hy","--hydrogenase_script_path",help="the path to the hydrogenase script",default=default_hydrogenase_script_path)
    args.add_argument("-dsr","--dsr_database_fasta_file",help="the path to the dsr database",default=default_dsr_database_fasta_file)
    args.add_argument("-mcr","--mcr_database_fasta_file",help="the path to the mcr database",default=default_mcr_database_fasta_file)
    #add an argument for the DIADB path
    args.add_argument("-di","--diadb_path",help="the path to the DIADB database",default=default_DIADB_dir)
    #add the argument to run the neighbourhood search
    args.add_argument("-n","--neighbourhood_search",help="search the neighbourhood of FeFe A hydrogenases to determine the subgroup",action="store_true",default=True)
    #add a flag for max threads
    #add an argument for the CAT path
    args.add_argument("-c","--cat_path",help="the path to the CAT executable",default=default_cat_path)
    args.add_argument("-th","--max_threads",help="the maximum number of threads to use",default=1,type=int)
    #add a flag for batch size
    args.add_argument("-b","--batch_size",help="the batch size",default=100,type=int)
    args.add_argument('-t',"--tmp_dir",help="the path to the temporary directory",default=tmp_folder)
    #add a flag to run the acetoscan
    args.add_argument("-a","--acetoscan",help="run acetoscan",action="store_true",default=False)
    #add an argument for file extensions to be considered
    args.add_argument("-x","--file_extensions",help="the file extensions to be considered",default=[".fa",".fna",".fasta"],nargs="+")
    #add the path to the hydrogenase activity file
    args.add_argument("-hyd","--hydrogen_classes_activity_file",help="the path to the hydrogen classes activity file",default=DEFAULT_HYDROGENASE_ACTIVITY_FILE)
    #add a flag to enable resume mode, to continue in case of failure,its default is true
    args.add_argument("-r","--resume",help="resume mode, in case the previous run crash, use cached output from previous run",action="store_true",default=False)
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
    #create the feature table
    create_feature_table(args.tmp_dir, args.hydrogen_classes_activity_file, args.input_dir
                         , args.output_file)
    
def tese4():
    input_folder=r'C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test\test_genomes'
    tmp_folder='/tllhome/abel/hydrogen_consumer_pipeline/test/test_genomes/tmp'
    output_file=r'C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test\test_results.csv'
    #r"C:\Users\abel\Documents\hydrogenases\hydrogenase_subunits\test_genomes\genomes"
    #tmp_folder=r"C:\Users\abel\Documents\hydrogenases\hydrogenase_subunits\test_genomes\genomes\tmp"
    #output_file=r"C:\Users\abel\Documents\hydrogenases\hydrogenase_subunits\test_genomes\genomes\results.csv"
    print("RUNNING")
    args=argparse.ArgumentParser()
    args.add_argument("-i","--input_dir",help="the input directory containing the genomes to be analyzed",default=input_folder,
                      required=False)
    args.add_argument("-b","--batch_size",help="the batch size",default=1,type=int)
    args.add_argument("-o","--output_file",help="the output file to save the results",required=False,default=output_file)
    args.add_argument("-p","--prodigal_exe_path",help="the path to the prodigal executable",default=r"D:\prodigal\prodigal.windows.exe" )
    args.add_argument("-d","--diamond_exe_path",help="the path to the diamond executable",default=r"D:\diamond-windows\diamond.exe")
    args.add_argument("-hy","--hydrogenase_script_path",help="the path to the hydrogenase script",default=r"C:\Users\abel\Documents\hydrogenases\hyDB\Diamond_blast_hyDB.py")
    args.add_argument("-dsrA","--dsrA_database_fasta_file",help="the path to the dsr database",default=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\dsrA.fasta")
    #apra
    args.add_argument("-apra","--aprA_database_fasta_file",help="the path to the apra database",default=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\aprA.fasta")
    #asrA
    args.add_argument("-asrA","--asrA_database_fasta_file",help="the path to the asrA database",default=r"C:\Users\abel\Documents\hydrogenases\sulphate_reducer\asrA.fasta")
    #"C:\Users\abel\Documents\hydrogenases\trimethylamine_N_oxide_reduction\DmsA.fasta"
    args.add_argument("-dmsA","--dmsA_database_fasta_file",help="the path to the dmsA database",default=r"C:\Users\abel\Documents\hydrogenases\trimethylamine_N_oxide_reduction\DmsA.fasta")
    #"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\napA.fasta"
    args.add_argument("-napA","--napA_database_fasta_file",help="the path to the napA database",default=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\napA.fasta")
    #"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\narG.fasta"
    args.add_argument("-narG","--narG_database_fasta_file",help="the path to the narG database",default=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\narG.fasta")
    #"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\nrfA.fasta"
    args.add_argument("-nrfA","--nrfA_database_fasta_file",help="the path to the nrfA database",default=r"C:\Users\abel\Documents\hydrogenases\nitrate_reducer\nrfA.fasta")
    #c:\Users\abel\Documents\hydrogenases\fumarate_reducer\frdA.fasta
    args.add_argument("-frdA","--frdA_database_fasta_file",help="the path to the frdA database",default=r"c:\Users\abel\Documents\hydrogenases\fumarate_reducer\frdA.fasta")
    #"C:\Users\abel\Documents\hydrogenases\cydA.fasta"
    args.add_argument("-cydA","--cydA_database_fasta_file",help="the path to the cydA database",default=r"C:\Users\abel\Documents\hydrogenases\cydA.fasta")
    #"C:\Users\abel\Documents\hydrogenases\hyd_production\fhl.fasta"
    args.add_argument("-fhl","--fhl_database_fasta_file",help="the path to the fhl database",default=r"C:\Users\abel\Documents\hydrogenases\hyd_production\FHL.fasta")
    #"C:\Users\abel\Documents\hydrogenases\hyd_production\PFOR.fasta"
    args.add_argument("-pfor","--pfor_database_fasta_file",help="the path to the pfor database",default=r"C:\Users\abel\Documents\hydrogenases\hyd_production\PFOR.fasta")
    args.add_argument("-nifH","--nitrogenase_database_fasta_file",help="the path to the nifh database",default=r"C:\Users\abel\Documents\hydrogenases\hyd_production\nifH.fasta")
    args.add_argument("-acsb","--acsB_database_fasta_file",help="the path to the acetate kinase database",default=r"C:\Users\abel\Documents\hydrogenases\acetogen\acsb.fasta")
    args
    args.add_argument("-mcrA","--mcrA_database_fasta_file",help="the path to the mcr database",default=r"C:\Users\abel\Documents\hydrogenases\methanogen\mcrA.fasta")
    #add an argument for the DIADB path
    args.add_argument("-di","--diadb_path",help="the path to the DIADB database",default=r'C:\Users\abel\Documents\hydrogenases\hyDB')
    #add the argument to run the neighbourhood search
    args.add_argument("-x","--file_extensions",help="the file extensions to be considered",default=[".fa",".fna",".fasta"],nargs="+")
    #add the path to the hydrogenase activity file
    args.add_argument("-hyd","--hydrogen_classes_activity_file",help="the path to the hydrogen classes activity file",default=r"C:\Users\abel\Documents\hydrogenases\hydrogenase_classes_activity.csv")
    #add a flag to enable resume mode, to continue in case of failure,its default is true
    args.add_argument("-r","--resume",help="resume mode, in case the previous run crash, use cached output from previous run",action="store_true",default=True)
    args.add_argument('-t',"--tmp_dir",help="the path to the temporary directory",default="tmp")
    #add the argument to run the neighbourhood search
    args.add_argument("-n","--neighbourhood_search",help="search the neighbourhood of FeFe A hydrogenases to determine the subgroup",action="store_true",default=True)
    #add the taxonomic table path
    #args.add_argument("-tax","--taxonomy_table_path",help="the path to the taxonomy table",default=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline_HPC\hydrogen_consumer_pipeline\test\test_metadata2.csv")
    args.add_argument("-taxt","--taxonomy_table_path",help="the path to the taxonomy table",default=None)
    args.add_argument("-c","--cat_path",help="the path to the CAT executable",default=default_cat_path)
    #add the flag for cat database path
    args.add_argument("-catdb","--cat_database_path",help="the path to the CAT database",default=default_CAT_database_path)
    #add the flag for the taxonomy path
    args.add_argument("-tax","--cat_taxonomy_path",help="the path to the CAT taxonomy",default=default_CAT_taxonomy_path)
    args=args.parse_args()
    #batch the files
    batch_root_dir=batch_the_files(args.input_dir,os.path.join(args.input_dir,"batches"),batch_size=args.batch_size,resume=args.resume)

    #run the batches
    run_batches(batch_root_dir,args.output_file,args.tmp_dir,args)
    #create_feature_table(args.tmp_dir, args.hydrogen_classes_activity_file, args.input_dir
    #                     , args.output_file,taxon_table=args.taxonomy_table_path)
def test_is_hydrogen():
    hyd_output=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\test_run\batches\batch_0\tmp\all_hydrogenase_fixed.txt"
    is_hydrogenase2(hyd_output)




def test_CAT_table_input():
    table_file=r"C:\Users\abel\Documents\mouse_gut\mouse_gut_catalogue\mgnify_,mouse_catalogue\genomes-all_metadata.csv"
    input_folder=r"C:\Users\abel\Documents\mouse_gut\mouse_gut_catalogue\mgnify_,mouse_catalogue\test_table"
    run_CAT(input_folder,r"C:\Users\abel\Documents\mouse_gut\mouse_gut_catalogue\mgnify_,mouse_catalogue\tmp",
            resume_mode=False,dummy_data=True,CAT_executable_path=default_cat_path,diamond_exe_path=default_diamond_exe_path,
                CAT_database_path=default_CAT_database_path,CAT_taxonomy_path=default_CAT_taxonomy_path,taxonomy_table=table_file)




if __name__=="__main__":
    #tese4()
    main()
    #test_CAT_table_input()
    #test_is_hydrogen()
    #tese3()
    #test_batch_run2()
    #exit()
    #main()
    #tes()
    #tes2()
    #tes_diamond_batch()
    #test_acetoscan_batch()
    #test_hydrogenase_batch()
    #test_results_combine()
    #test_batch_run2()
