import pandas as pd
from Bio import SeqIO
import random
import subprocess
import os
mcra_file=r"C:\Users\abel\Documents\marker_sens_spec\mcra_table.tsv"
all_file=r"C:\Users\abel\Documents\marker_sens_spec\all_reviewed.tsv"

mcra_fasta=r"C:\Users\abel\Documents\marker_sens_spec\mcra.fasta"
all_fasta=r"C:\Users\abel\Documents\marker_sens_spec\all_reviewed.fasta"
diamond_exe=r"/tllhome/abel/hydrogen_consumer_pipeline/DIAMOND/diamond
hyddb_file=r"C:\Users\abel\Documents\hydrogenases\hyDB\HydDB_reformated.fasta"
hyddb_table=r"C:\Users\abel\Documents\hydrogenases\hyDB\hydrogenase_classes_activity2t.csv"
def make_negative_fasta(all_fasta,  all_file, mcra_file,cwd_folder):
    #read the all reviewed file
    all_df=pd.read_csv(all_file, sep='\t')
    mcr_df=pd.read_csv(mcra_file, sep='\t')

    mcra_entries=mcr_df['Entry'].tolist()
    all_entries=all_df['Entry'].tolist()
    #get the all entries that are not in the mcra entries
    not_in_mcra=[entry for entry in all_entries if entry not in mcra_entries]
    in_mcra=[entry for entry in all_entries if entry in mcra_entries]
    import time
    negative_records=[]
    num_recs_parsed=0
    time1=time.time()
    #open the all reviewed fasta file
    with open(all_fasta, "r") as handle:
        #count number of >
        if False:
            for record in SeqIO.parse(handle, "fasta"):
                #print("id:", record.id)
                num_recs_parsed+=1
                if num_recs_parsed%1000==0:
                    print("Records parsed:", num_recs_parsed)
                    time2=time.time()
                    print("Time taken:", time2-time1)
                    #estimate time remaining
                    time_per_rec=(time2-time1)/num_recs_parsed
                    remaining_recs=len(all_entries)-num_recs_parsed
                    remaining_time=remaining_recs*time_per_rec
                    print("Remaining time:", remaining_time)
                #split by | and get the 2nd element
                ide=record.id.split("|")[1]
                if ide in not_in_mcra:
                    negative_records.append(record)
        #split by >
        recs=handle.read().split("\n>")
        for rec in recs:
            if True:
                rec=rec.strip()
                if len(rec)==0:
                    continue
                #split by \n
                lines=rec.split("\n")
                #get the first line
                first_line=lines[0]
                #split by | and get the 2nd element
                ide=first_line.split("|")[1]
                if ide not in in_mcra:
                    negative_records.append(rec)
                #    pass
                    #negative_records.append(rec)
            num_recs_parsed+=1
            if num_recs_parsed%1000==0:
                print("Records parsed:", num_recs_parsed)
                time2=time.time()
                print("Time taken:", time2-time1)
                #estimate time remaining
                time_per_rec=(time2-time1)/num_recs_parsed
                remaining_recs=len(all_entries)-num_recs_parsed
                remaining_time=remaining_recs*time_per_rec
                print("Remaining time:", remaining_time)
    #write the negative records to a fasta file
    num_recs_written=0
    with open(os.path.join(cwd_folder, "negative.fasta"), "w") as handle:
        for record in negative_records:
            num_recs_written+=1
            if num_recs_written%1000==0:
                print("Records written:", num_recs_written)

            handle.write(">"+record+"\n")



    #writing the negative records to a fasta file
    #SeqIO.write(negative_records, os.path.join(cwd_folder, "negative.fasta"), "fasta")
    return os.path.join(cwd_folder, "negative.fasta")


def parse_results(train_file,test_file, test_results, negative_file, negative_results, min_pid, min_scov):
    #parse the train file and get the length of the sequences
    train_lengths={}
    with open(train_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            train_lengths[record.id]=len(record.seq)
    #count the number of test and negative records
    test_count=0
    with open(test_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            test_count+=1
    negative_count=0
    with open(negative_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            negative_count+=1
    #read the test results
    def get_scov(row):
        return row['aln_len']/train_lengths[row['subject']]
    #if test df is empty, return None
    if os.path.getsize(test_results)==0:
        true_positives=0
    else:
        test_df=pd.read_csv(test_results, sep='\t', header=None)
        test_df.columns=["query", "subject", "pid", "aln_len", "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
        test_df['scov']=test_df.apply(get_scov, axis=1)
        test_df=test_df[(test_df['pid']>=min_pid) & (test_df['scov']>=min_scov)]
        true_positives=test_df.shape[0]
    #read the negative results
    #if negative df is empty, return None
    
    if os.path.getsize(negative_results)==0:
        false_positives=0
    else:
        negative_df=pd.read_csv(negative_results, sep='\t', header=None)
        negative_df.columns=["query", "subject", "pid", "aln_len", "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
        #filter the test results
        negative_df['scov']=negative_df.apply(get_scov, axis=1)
        #filter the negative results
        negative_df=negative_df[(negative_df['pid']>=min_pid) & (negative_df['scov']>=min_scov)]
        false_positives=negative_df.shape[0]
    #make a olumn called scov
    
    #if both true positives and false positives are 0, return None
    if true_positives==0 and false_positives==0:
        row={"gene":"mcrA", "pid":min_pid, "scov":min_scov, "sensitivity":None, "specificity":None, "precision":None, "accuracy":None}
        return row
    
    
    #get the true positives
    
    #get the false positives
    
    #get the true negatives
    true_negatives=negative_count-false_positives
    #get the false negatives
    false_negatives=test_count-true_positives
    #get the sensitivity
    sensitivity=true_positives/(true_positives+false_negatives)
    #get the specificity
    specificity=true_negatives/(true_negatives+false_positives)
    #get the accuracy
    accuracy=(true_positives+true_negatives)/(true_positives+true_negatives+false_positives+false_negatives)
    #get the precision
    precision=true_positives/(true_positives+false_positives)
    #make the row, gene name (mcrA), pid, scov, sensitivity, specificity, precision, accuracy
    row={"gene":"mcrA", "pid":min_pid, "scov":min_scov, "sensitivity":sensitivity, "specificity":specificity, "precision":precision, "accuracy":accuracy}
    return row
#a

def parse_results_batched(train_file,test_file, test_results, negative_file, negative_results, min_pids, min_scovs,gene_name,filtered_folder,train_percent,trial):
    #if filtered folder does not exist, make it
    if not os.path.exists(filtered_folder):
        os.mkdir(filtered_folder)
    #parse the train file and get the length of the sequences
    train_lengths={}
    with open(train_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            train_lengths[record.id]=len(record.seq)
    #count the number of test and negative records
    test_count=0
    with open(test_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            test_count+=1
    negative_count=0
    with open(negative_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            negative_count+=1
    #read the test results
    def get_scov(row):
        return row['aln_len']/train_lengths[row['subject']]
    #if test df is empty, return None
    rows=[]
    for min_pid in min_pids:
        for min_scov in min_scovs:
            if os.path.getsize(test_results)==0:
                true_positives=0
            else:
                test_df=pd.read_csv(test_results, sep='\t', header=None)
                test_df.columns=["query", "subject", "pid", "aln_len", "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
                test_df['scov']=test_df.apply(get_scov, axis=1)
                test_df=test_df[(test_df['pid']>=min_pid) & (test_df['scov']>=min_scov)]
                tag=str(train_percent)+"_"+str(min_pid)+"_"+str(min_scov)+"_"+str(trial)
                #save the filtered results to a file
                test_df.to_csv(os.path.join(filtered_folder, gene_name+"_"+tag+"_test_results.tsv"), sep='\t', index=False)
                true_positives=test_df.shape[0]
            #read the negative results
            #if negative df is empty, return None
            
            if os.path.getsize(negative_results)==0:
                false_positives=0
            else:
                negative_df=pd.read_csv(negative_results, sep='\t', header=None)
                negative_df.columns=["query", "subject", "pid", "aln_len", "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
                #filter the test results
                negative_df['scov']=negative_df.apply(get_scov, axis=1)
                tag=str(train_percent)+"_"+str(min_pid)+"_"+str(min_scov)+"_"+str(trial)
                #filter the negative results
                negative_df=negative_df[(negative_df['pid']>=min_pid) & (negative_df['scov']>=min_scov)]
                false_positives=negative_df.shape[0]
                #save the filtered results to a file
                negative_df.to_csv(os.path.join(filtered_folder, gene_name+"_"+tag+"_negative_results.tsv"), sep='\t', index=False)
            #make a olumn called scov
            
            #if both true positives and false positives are 0, return None
            if true_positives==0 and false_positives==0:
                row={"gene":gene_name, "pid":min_pid, "scov":min_scov, "trial":trial, "train_percent":train_percent, 
                     "sensitivity":None, "specificity":None, "precision":None, "accuracy":None}
                rows.append(row)
                continue
        
    
            #get the true positives
            
            #get the false positives
            
            #get the true negatives
            true_negatives=negative_count-false_positives
            #get the false negatives
            false_negatives=test_count-true_positives
            #get the sensitivity
            sensitivity=true_positives/(true_positives+false_negatives)
            #get the specificity
            specificity=true_negatives/(true_negatives+false_positives)
            #get the accuracy
            accuracy=(true_positives+true_negatives)/(true_positives+true_negatives+false_positives+false_negatives)
            #get the precision
            precision=true_positives/(true_positives+false_positives)
            #make the row, gene name (mcrA), pid, scov, sensitivity, specificity, precision, accuracy
            row={"gene":gene_name, "pid":min_pid, "scov":min_scov, "trial":trial, "train_percent":train_percent, 
                 "sensitivity":sensitivity, "specificity":specificity, "precision":precision, "accuracy":accuracy}
            rows.append(row)
    return rows
#a
def diamond_blast(train_file, test_file, negative_recs_file,diamond_exe,cwd_folder,tag=""):
    #make the diamond database
    command=[diamond_exe, "makedb", "--in", train_file, "-d", os.path.join(cwd_folder, "train_db")]
    subprocess.run(command, check=True)
    #run the diamond blast
    command=[diamond_exe, "blastp", "-d", os.path.join(cwd_folder, "train_db"), "-q", test_file, "-o", os.path.join(cwd_folder, tag+"_test_results.tsv")
                ,'--max-target-seqs', '1']
    subprocess.run(command, check=True)
    #run the diamond blast for the negative records
    command=[diamond_exe, "blastp", "-d", os.path.join(cwd_folder, "train_db"), "-q", negative_recs_file, 
             "-o", os.path.join(cwd_folder, tag+"_negative_results.tsv"),'--max-target-seqs', '1']
    subprocess.run(command, check=True)
    return os.path.join(cwd_folder, tag+"_test_results.tsv"), os.path.join(cwd_folder,tag+ "_negative_results.tsv")
def training_split(positive_file, train_percent,cwd_folder,tag=""):
    positive_dir=os.path.dirname(positive_file)
    #read the positive file
    positive_records=[]
    with open(positive_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            positive_records.append(record)
    #shuffle the positive records
    random.shuffle(positive_records)
    #split the positive records
    train_size=int(len(positive_records)*train_percent)
    train_records=positive_records[:train_size]
    test_records=positive_records[train_size:]
    #write the train records to a fasta file
    SeqIO.write(train_records, os.path.join(cwd_folder, tag+"_train.fasta"), "fasta")
    SeqIO.write(test_records, os.path.join(cwd_folder, tag+"_test.fasta"), "fasta")
    return os.path.join(cwd_folder, tag+"_train.fasta"), os.path.join(cwd_folder, tag+"_test.fasta")

def run (positive_file, negative_file, diamond_exe, cwd_folder, train_percent, min_pid, min_scov,trial):
    #split the training file
    tag=str(train_percent)+"_"+str(min_pid)+"_"+str(min_scov)+"_"+str(trial)
    train_file, test_file=training_split(positive_file, train_percent,cwd_folder,tag)
    #run the diamond blast
    #if tagged file exists
    if os.path.exists(os.path.join(cwd_folder, tag+"_test_results.tsv")):
        test_results=os.path.join(cwd_folder, tag+"_test_results.tsv")
        negative_results=os.path.join(cwd_folder, tag+"_negative_results.tsv")
        print("tagged file exists", test_results)
        row=parse_results(train_file, test_file, test_results, negative_file, negative_results, min_pid, min_scov)
        return row
    else:
        test_results, negative_results=diamond_blast(train_file, test_file, negative_file, diamond_exe, cwd_folder, tag)
    #parse the results
    row=parse_results(train_file, test_file, test_results, negative_file, negative_results, min_pid, min_scov)
    return row



def run_multiple(positive_file, negative_file, diamond_exe, cwd_folder,train_percents, min_pids, min_scovs,trial_count,marker_name):
    rows=[]
    for train_percent in train_percents:
        for min_pid in min_pids[0:1]:
            for min_scov in min_scovs[0:1]:
                for trial in range(trial_count):
                    row=run(positive_file, negative_file, diamond_exe, cwd_folder, train_percent, min_pid, min_scov,trial)
                    row["trial"]=trial
                rows.append(row)
    #make a dataframe
    df=pd.DataFrame(rows)
    #save the dataframe
    df.to_csv(os.path.join(cwd_folder, marker_name+"_results.tsv"), sep='\t', index=False)
    return df



cwd_folder=r"C:\Users\abel\Documents\marker_sens_spec"


def run_marker(marker_file,table_file,cwd_folder):
    positive_file=marker_file
    #make the cwd folder
    if not os.path.exists(cwd_folder):
        os.mkdir(cwd_folder)
    #check if positive file exists
    if not os.path.exists(positive_file):
        raise Exception("File does not exist")
    #check if table file exists
    if not os.path.exists(table_file):
        raise Exception("Table file does not exist")
    #make the negative file
    negative_file=os.path.join(cwd_folder, "negative.fasta")
    if not os.path.exists(negative_file):
        negative_file=make_negative_fasta(all_fasta,  all_file, marker_file,cwd_folder)
    train_percent=0.8
    min_pids=[40,60,80,90]
    min_scovs=[0.5,0.7,0.9]
    train_percents=[0.8,0.2]
    trial_count=3
    marker_name=os.path.basename(table_file).split(".")[0]
    marker_name=marker_name.replace("_table", "")
    #Run the multiple
    df=run_multiple(positive_file, negative_file, diamond_exe, cwd_folder,train_percents, min_pids, min_scovs,trial_count,marker_name)
    #Save the results
    df.to_csv(os.path.join(cwd_folder, marker_name+"_results.tsv"), sep='\t', index=False)


from Bio import SeqIO

negative_file=r"/tllhome/abel\marker_sens_spec\negative.fasta"
def make_hydrogenase_groups(hydb_file, hyddb_table,out_folder):
    hydb_df=pd.read_csv(hyddb_table)
    hyd_groups={}
    for index, row in hydb_df.iterrows():
        class1= row["Class"]
        group=row["Group"]
        #if group is A11, ignore
        if group=="A11":
            continue
        #if class is NiFeSe,ignore
        if class1=="NiFeSe":
            continue
        #if class == Fe, hyd group is just [Fe]
        if class1=="Fe":
            hyd_groups["["+class1+"]"]=[]
        else:
            hyd_groups["["+class1+"]_Group_"+group]=[]
    recs=SeqIO.parse(hydb_file, "fasta")
    for rec in recs:
        #get the description
        desc=rec.description
        #get the class
        for key in hyd_groups:
            if key in desc:
                hyd_groups[key].append((rec.id, rec))
                print("appended", rec.id,"in group", key)
                break
    all_recs=[(rec.id, rec) for rec in SeqIO.parse(hydb_file, "fasta")]
    #parse the negative records
    negative_recs=[]
    with open(negative_file, "r") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            negative_recs.append(rec)
    #write the hydrogenase groups to a folder
    for key in hyd_groups:
        recss=[rec[1] for rec in hyd_groups[key]]
        SeqIO.write(recss, os.path.join(out_folder, key+".fasta"), "fasta")
        #get the recs which are not in the hyd groups
        not_in_group=[rec[1] for rec in all_recs if rec[0] not in [rec[0] for rec in hyd_groups[key]]]
        #add the negative recs to the not in group
        not_in_group+=negative_recs
        SeqIO.write(not_in_group, os.path.join(out_folder, key+"_negative.fasta"), "fasta")

    return hyd_groups
    
def run_hydrogenase(pos_file, neg_file, diamond_exe, cwd_folder,train_percents, min_pids, min_scovs,trial_count,marker_name,hyd_group):
    rows=[]
    for train_percent in train_percents:
        for min_pid in min_pids[0:1]:
            for min_scov in min_scovs[0:1]:
                for trial in range(trial_count):
                    
                        positive_file=os.path.join(cwd_folder, hyd_group+".fasta")
                        negative_file=os.path.join(cwd_folder, hyd_group+"_negative.fasta")
                        #make a folder for the hydrogenase group
                        if not os.path.exists(os.path.join(cwd_folder, hyd_group)):
                            os.mkdir(os.path.join(cwd_folder, hyd_group))
                        hyd_folder=os.path.join(cwd_folder, hyd_group)
                        row=run(positive_file, negative_file, diamond_exe, hyd_folder, train_percent, min_pid, min_scov,trial)
                        row["trial"]=trial
                        row["group"]=hyd_group
                        rows.append(row)
    #make a dataframe
    df=pd.DataFrame(rows)
    #save the dataframe
    df.to_csv(os.path.join(cwd_folder, marker_name+"_results.tsv"), sep='\t', index=False)
    return df
def run_all_hydrogenases( cwd_folder):
    train_percents=[0.8,0.2]
    min_pids=[40,60,80,90]
    min_scovs=[0.5,0.7,0.9]
    trial_count=3
    for file in os.listdir(cwd_folder):
        if not file.endswith(".fasta"):
            continue
        #skip if negative file
        if "negative" in file:
            continue
        marker_name=file.split(".")[0]
        run_hydrogenase(os.path.join(cwd_folder, file), os.path.join(cwd_folder, file.replace(".fasta", "_negative.fasta")), 
                        diamond_exe, cwd_folder,train_percents, min_pids, min_scovs,trial_count,marker_name,hyd_group=file.split(".")[0])





def aggregate_hyd_files(cwd_folder):
    #get the hyd groups form the filenames
    hyd_groups={}
    for file in os.listdir(cwd_folder):
        if not file.endswith(".fasta"):
            continue
        if "negative" in file:
            continue
        key=file.split(".")[0]
        hyd_groups[key]=[]
    #read the files
    for key in hyd_groups:
        df=pd.read_csv(os.path.join(cwd_folder, key+"_results.tsv"), sep='\t')
        df["group"]=key
        hyd_groups[key]=df
    #concatenate the dataframes
    df=pd.concat(hyd_groups.values())
    #save the dataframe
    df.to_csv(os.path.join(cwd_folder, "all_hydrogenase_results.tsv"), sep='\t', index=False)
    return df





def main():
    marker_dir=r"C:\Users\abel\Documents\marker_sens_spec\final_curated_marker"
    table_folder=r"C:\Users\abel\Documents\marker_sens_spec\tables"
    all_table=r"C:\Users\abel\Documents\marker_sens_spec\all_reviewed.tsv"
    root_dir=r"C:\Users\abel\Documents\marker_sens_spec"
    for file in os.listdir(marker_dir):
        #skip if not fasta file
        if not file.endswith(".fasta"):
            continue
        #find the table file, replace fasta with table
        table_file=file.replace(".fasta", "_table.tsv")
        table_file=os.path.join(table_folder, table_file)
        #if table file no exist, raise exception
        if not os.path.exists(table_file):
            print(f"Table file does not exist: {table_file}")   
            #raise Exception(f"Table file does not exist: {table_file}")
            continue
        #make the negative fasta
        #make a cwd folder called marker_name   
        cwd_folder=os.path.join(root_dir, os.path.basename(file).split(".")[0])
        if not os.path.exists(cwd_folder):
            os.mkdir(cwd_folder)
        #check if negative file exists, if not make it
        negative_file=os.path.join(cwd_folder, "negative.fasta")
        if not os.path.exists(negative_file):
            negative_file=make_negative_fasta(all_fasta,  all_file, os.path.join(marker_dir, file),cwd_folder)
        #negative_file=make_negative_fasta(all_fasta,  all_file, table_file,cwd_folder)
            print("made negative file", negative_file)
        #run the marker
        run_marker(os.path.join(marker_dir, file), table_file,cwd_folder)
        print("ran marker", file)


def aggregate_files(root_dir):
    min_pids=[40,60,80,90]
    min_scovs=[0.5,0.7,0.9]
    train_percents=[0.8,0.2]
    trial_count=3
    markers=["acsB","aprA","cydA","DmsA","dsrA","FHL","frdA","mcrA","narG","nifH","nrfA","PFOR","hydrogenase"]
    for fold in markers:
        rows=[]
        for min_pid in min_pids:
            for min_scov in min_scovs:
                for train_percent in train_percents:
                    
                    # parse_results(train_file,test_file, test_results, negative_file, negative_results, min_pid, min_scov)
                    for trial in range(trial_count):
                        tag=str(train_percent)+"_"+str(min_pid)+"_"+str(min_scov)+"_"+str(trial)
                        print(fold,tag)
                        train_file=os.path.join(root_dir, fold, tag+"_train.fasta")
                        test_file=os.path.join(root_dir, fold, tag+"_test.fasta")
                        test_results=os.path.join(root_dir, fold, tag+"_test_results.tsv")
                        negative_file=os.path.join(root_dir, fold, "negative.fasta")
                        negative_results=os.path.join(root_dir, fold, tag+"_negative_results.tsv")
                        row=parse_results(train_file, test_file, test_results, negative_file, negative_results, min_pid, min_scov)
                        row["trial"]=trial
                        row["gene"]=fold
                        rows.append(row)
        #make a dataframe
        df=pd.DataFrame(rows)
        #save the dataframe in root dir
        df.to_csv(os.path.join(root_dir, fold+"_results.tsv"), sep='\t', index=False)
        print("saved", fold)


def aggregate_files_batched(root_dir):
    min_pids=[40,60,80,90]
    min_scovs=[0.5,0.7,0.9]
    train_percents=[0.8,0.2]
    trial_count=3
    markers=["acsB","aprA","cydA","DmsA","dsrA","FHL","frdA","mcrA","narG","nifH","nrfA","PFOR","hydrogenase"]
    for fold in markers[-2:-1]:
        rows=[]
        
        for train_percent in train_percents:
                    
                    # parse_results(train_file,test_file, test_results, negative_file, negative_results, min_pid, min_scov)
                    for trial in range(trial_count):
                        tag=str(train_percent)+"_"+str(min_pids[0])+"_"+str(min_scovs[0])+"_"+str(trial)
                        print(fold,tag)
                        train_file=os.path.join(root_dir, fold, tag+"_train.fasta")
                        test_file=os.path.join(root_dir, fold, tag+"_test.fasta")
                        test_results=os.path.join(root_dir, fold, tag+"_test_results.tsv")
                        negative_file=os.path.join(root_dir, fold, "negative.fasta")
                        negative_results=os.path.join(root_dir, fold, tag+"_negative_results.tsv")
                        rows1=parse_results_batched(train_file, test_file, test_results, negative_file, negative_results,
                                          min_pids=min_pids, min_scovs=min_scovs,gene_name=fold,
                                          filtered_folder=os.path.join(root_dir, fold,"filtered"),train_percent=train_percent,trial=trial)
                        rows+=rows1
        #make a dataframe
        df=pd.DataFrame(rows)
        #save the dataframe in root dir
        df.to_csv(os.path.join(root_dir, fold+"_results.tsv"), sep='\t', index=False)
        print("saved", fold)

#run the PFOR
marker_dir=r"C:\Users\abel\Documents\marker_sens_spec\final_curated_marker"
table_folder=r"C:\Users\abel\Documents\marker_sens_spec\tables"
table_file=r"C:\Users\abel\Documents\marker_sens_spec\tables\PFOR_table.tsv"
cwd_folder=r"C:\Users\abel\Documents\marker_sens_spec\PFOR"
file="PFOR.fasta"

#run_marker(os.path.join(marker_dir, file), table_file,cwd_folder)


out_folder=r'/tllhome/abel/marker_sens_spec'
hydb_file=r"C:/tllhome/abel/hydrogen_consumer_pipeline/hyDB/HydDB_reformated.fasta"
hyddb_table=r"/tllhome/abel/hydrogen_consumer_pipeline/hyDB/hydrogenase_classes_activity2t.csv"
#make_hydrogenase_groups(hydb_file, hyddb_table,out_folder)

#hyd_groups=make_hydrogenase_groups(hydb_file, hyddb_table,out_folder)
cwd_folder=r"/tllhome/abel/marker_sens_spec/hyd"
run_all_hydrogenases( cwd_folder)
aggregate_hyd_files(cwd_folder)
#aggregate_files_batched(r"C:\Users\abel\Documents\marker_sens_spec")






