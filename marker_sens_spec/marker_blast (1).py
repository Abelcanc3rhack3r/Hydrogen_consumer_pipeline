import pandas as pd
from Bio import SeqIO
import random
import subprocess
import os
mcra_file=r"C:\Users\abel\Documents\marker_sens_spec\mcra_table.tsv"
all_file=r"C:\Users\abel\Documents\marker_sens_spec\all_reviewed.tsv"

mcra_fasta=r"C:\Users\abel\Documents\marker_sens_spec\mcra.fasta"
all_fasta=r"C:\Users\abel\Documents\marker_sens_spec\all_reviewed.fasta"
diamond_exe=r"C:\Users\abel\Documents\marker_sens_spec\diamond.exe"
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


def parse_results_batched_hydrogenases(train_file,test_file, test_results, min_pids, min_scovs,gene_name,filtered_folder,train_percent,trial):
    #the test file also contains the negative records
    #if filtered folder does not exist, make it

    def extract_group(row,col):
        query=row[col]
        #split by " - " and get the last element
        groups=query.split("_-_")
        if len(groups)==1:
            return None
        return groups[-1]
    def extract_group1(strr):
        query=strr
        #split by " - " and get the last element
        groups=query.split("_-_")
        if len(groups)==1:
            return None
        return groups[-1]
    if not os.path.exists(filtered_folder):
        os.mkdir(filtered_folder)
    #parse the train file and get the length of the sequences
    train_lengths={}
    with open(train_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            train_lengths[record.id]=len(record.seq)
    #count the number of test and negative records
    test_count={}
    with open(test_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            #extract the group
            group=extract_group1(record.id+"_"+record.description)
            if group not in test_count:
                test_count[group]=0
            test_count[group]+=1

    negative_count=0
    
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
            
            #read the filtered results
            if os.path.getsize(os.path.join(filtered_folder, gene_name+"_"+tag+"_test_results.tsv"))==0:
                false_positives=0
            else:
                hyd_dict={}
                for index, row in test_df.iterrows():
                    query_group=extract_group(row,"query")
                    subject_group=extract_group(row,"subject")
                    if subject_group is not None:
                        if subject_group not in hyd_dict:
                            hyd_dict[subject_group]={"TP":0,"Non hyd FP":0}

                    if query_group is not None and subject_group is not None:
                        if query_group==subject_group:
                            hyd_dict[subject_group]["TP"]+=1
                        else:
                            q_group=extract_group(row,"query")
                            if q_group is not None:
                                if q_group not in hyd_dict[subject_group]:
                                    hyd_dict[subject_group][q_group]=0
                                hyd_dict[subject_group][q_group]+=1
                    if query_group is None and subject_group is not None:
                        hyd_dict[subject_group]["Non hyd FP"]+=1
                
                non_hyd_false_positives=0
                for group in hyd_dict:
                    false_positives=0
                    false_positives+=hyd_dict[group]["Non hyd FP"]
                    true_positives=hyd_dict[group]["TP"]
                    for key in hyd_dict[group]:
                        if key!="TP" and key!="Non hyd FP":
                            false_positives+=hyd_dict[group][key]
                    non_hyd_false_positives=hyd_dict[group]["Non hyd FP"]
                    non_hyd_specificity=1-(non_hyd_false_positives/test_count[None])
                    #get the true negatives
                    non_group_counts=0
                    for group1 in test_count:
                        if  group1!=group:
                            non_group_counts+=test_count[group1]
                    true_negatives=non_group_counts-false_positives
                    #get the false negatives
                    false_negatives=test_count[group]-true_positives
                    #get the sensitivity
                    #if true positives is 0, set sensitivity to 0
                    if true_positives==0:
                        sensitivity=0
                    else:
                        sensitivity=true_positives/(true_positives+false_negatives)
                    #sensitivity=true_positives/(true_positives+false_negatives)
                    #get the specificity
                    #if true negatives is 0, set specificity to 0
                    if true_negatives==0:
                        specificity=0
                    else:
                        specificity=true_negatives/(true_negatives+false_positives)
                    #specificity=true_negatives/(true_negatives+false_positives)
                    #if true positives and false positives are 0, set accuracy to 0
                    if true_positives==0 and false_positives==0:
                        accuracy=0
                    else:
                        accuracy=(true_positives+true_negatives)/(true_positives+true_negatives+false_positives+false_negatives)
                    #accuracy=(true_positives+true_negatives)/(true_positives+true_negatives+false_positives+false_negatives)
                    #get the precision
                    #if true positives and false positives are 0, set precision to 0
                    if true_positives==0 and false_positives==0:
                        precision=0
                    else:
                        precision=true_positives/(true_positives+false_positives)
                    #make the row, gene name (mcrA), pid, scov, sensitivity, specificity, precision, accuracy
                    row={"gene":group, "pid":min_pid, "scov":min_scov, "trial":trial, "train_percent":train_percent,
                        "sensitivity":sensitivity, "specificity":specificity, "precision":precision, "accuracy":accuracy,
                        "non_hyd_false_positives":non_hyd_false_positives, "non_hyd_specificity":non_hyd_specificity}
                    #add the hyd false positives
                    for key in hyd_dict[group]:
                        if key!="TP" and key!="Non hyd FP":
                            row[key+"_false positives"]=hyd_dict[group][key]
                            row[key+"_specificity"]=1-(hyd_dict[group][key]/test_count[key])
                    rows.append(row)
                    
    #save the filtered results to a file
    df=pd.DataFrame(rows)
    tag1=gene_name+"_"+str(train_percent)+"_"+str(trial)
    df.to_csv(os.path.join(filtered_folder, tag1+"test_results.tsv"), sep='\t', index=False)
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
    if negative_recs_file is not None:
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




def training_split_hyd(positive_file, negative_file,train_percent,cwd_folder,tag=""):
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
    #add the negative records to the test records
    with open(negative_file, "r") as handle:
        with open(os.path.join(cwd_folder, tag+"_test.fasta"), "a") as handle2:
            for record in SeqIO.parse(handle, "fasta"):
                SeqIO.write(record, handle2, "fasta")
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

#def train_test

hyddb=r"C:\Users\abel\Documents\hydrogenases\hyDB\HydDB_reformated.fasta"
negative_file=r"C:\Users\abel\Documents\hydrogenases\marker_sens_spec\negative.fasta"
def run_hydrogenases(hyd_file, negative_file, diamond_exe, cwd_folder,train_percents, min_pids, min_scovs,trial_count):
    #make a cwd folder
    #add the hyd file to the negative file
    test_file=hyd_file.replace(".fasta", "_test.fasta")
    #if test file no exists
    '''if not os.path.exists(test_file):
        with open(hyd_file, "r") as handle:
            with open(test_file, "w") as handle2:
                for record in SeqIO.parse(handle, "fasta"):
                    SeqIO.write(record, handle2, "fasta")
            with open(negative_file, "a") as handle3:
                for record in SeqIO.parse(handle, "fasta"):
                    SeqIO.write(record, handle3, "fasta")'''
    #cwd_folder=os.path.join(os.path.dirname(hyd_file), "results")
    if not os.path.exists(cwd_folder):
        os.mkdir(cwd_folder)
    #run the hydrogenases
    all_rows=[]
    for train_percent in train_percents:
        for trial in range(trial_count):
            tag="hydb_"+str(train_percent)+"_"+str(trial)
            #if test file and train file no exist
            if not os.path.exists(os.path.join(cwd_folder, tag+"_train.fasta")):
                train_file, test_file=training_split_hyd(hyd_file, negative_file, train_percent,cwd_folder,tag)
            #run the diamond blast
            #if tagged file exists
            if os.path.exists(os.path.join(cwd_folder, tag+"_test_results.tsv")):
                train_file=os.path.join(cwd_folder, tag+"_train.fasta")
                test_file=os.path.join(cwd_folder, tag+"_test.fasta")
                test_results=os.path.join(cwd_folder, tag+"_test_results.tsv")
                negative_results=os.path.join(cwd_folder, tag+"_negative_results.tsv")
                print("tagged file exists", test_results)
                #parse_results_batched_hydrogenases(train_file,test_file, test_results, min_pids, min_scovs,gene_name,filtered_folder,train_percent,trial):
                rows=parse_results_batched_hydrogenases(train_file, test_file,
                                                         test_results, min_pids, min_scovs,"hyd", os.path.join(cwd_folder, "filtered"),train_percent,trial)
                #return rows
                all_rows+=rows
            else:
                train_file=os.path.join(cwd_folder, tag+"_train.fasta")
                test_file=os.path.join(cwd_folder, tag+"_test.fasta")
                test_results, negative_results=diamond_blast(train_file, test_file, None, diamond_exe, cwd_folder, tag)
            #parse the results
            rows=parse_results_batched_hydrogenases(train_file, test_file, 
                                                   test_results, min_pids, min_scovs,"hyd", os.path.join(cwd_folder, "filtered"),train_percent,trial)
            all_rows+=rows
    #make a dataframe
    df=pd.DataFrame(all_rows)
    #save the dataframe
    df.to_csv(os.path.join(cwd_folder, "hyd_results.tsv"), sep='\t', index=False)
    #aggregate the files





cwd_folder=r"C:\Users\abel\Documents\marker_sens_spec"
run_hydrogenases(hyddb, negative_file, diamond_exe, cwd_folder,[0.8,0.2], [10,20,40,60,80,90], [0.1, 0.3,0.5,0.7,0.9],3)

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
from Bio import SeqIO
def count_groups(root_dir):
    group_counts={}
    for file in os.listdir(root_dir):
        #if file has [ and not negative .fasta
        if "[" in file and "negative" not in file:
            if file.endswith(".fasta"):
                with open(os.path.join(root_dir, file), "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        group=file.replace(".fasta", "")
                        if group not in group_counts:
                            group_counts[group]=0
                        group_counts[group]+=1
    return group_counts



def extract_group(row):
    query=row['query']
    #split by " - " and get the last element
    groups=query.split("_-_")
    if len(groups)==1:
        return None
    return groups[-1]
def parse_average_file(average_file):
    df=pd.read_csv(average_file, sep=',')
    group_probs={}
    av_copy={}
    for index, row in df.iterrows():
        group=row['group']
        splitw=group.split(" ")
        #if len of slpitw is less than 2, continue
        if len(splitw)<2:
            continue
        if group=="Fe":
            new_name="[Fe]"
        else:
            new_name="["+splitw[0]+"]_Group_"+splitw[1]
        prob=row['average']
        group_probs[new_name]=prob
        mouse_copy=row['Mouse_copy']
        human_copy=row['Human_copy']
        av_copy[new_name]=(mouse_copy+human_copy)/2

    return group_probs, av_copy
average_file=r"C:\Users\abel\Documents\marker_sens_spec\average.csv"
import pandas as pd
av_mouse_genes=2684
av_human_genes=2419
neg_file=r"C:\Users\abel\Documents\marker_sens_spec\hydrogenase\negative.fasta"
def count_neg_genes(neg_file):
    count=0
    with open(neg_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            count+=1
    return count
#neg_genes= count_neg_genes(neg_file)
#av_genes=(av_mouse_genes+av_human_genes)/2
import re


#for each false positive hyd group, calculate the probability of at least 1 hit in the average genome
# the probability is the probability that at least 1 copy of hyd group is present in the average genome * the probability of false detection
#the probability of false positive detection = num FP hits/num genes in group
#for non hyd genes, the probability of false positive detection = num FP hits/num genes in negative file(non hydrogenases)
#the probability of at least 1 non hydrogenase FP hit in the average genome = 1- (1-probability of false positive detection)^average genome size
#the probability of at least 1 FP hit in the average genome is weighted by the av copy number of each hyd group in genome , or
#for non hyd genes, the average number of genes in the genome
#the average number of genes is calculated by the unweighted average of the number of genes in mouse catalogue and human catalogue
#the number of hyd genes is taken from the hydDB.
#the probability of true positive detection = num TP hits/num genes in test set (1- split ratio * num genes in group)
#the probability of at least 1 TP hit in the average genome = TP detection prob * prob of at least 1 hyd gene in average genome

def parse_filtered(filtered_dir, group_countsd, average_dict,copy_dict,neg_genes):
    res={}
    for fold in os.listdir(filtered_dir):
        #if fold is not a directory, continue
        if not os.path.isdir(os.path.join(filtered_dir, fold)):
            continue
        group=fold
        #if fold is not in group countsd, set it to 0
        if group not in group_countsd:
            group_countsd[group]=0
        #if fold is not in average_dict, set it to 0
        if group not in average_dict:
            average_dict[group]=0
        for file in os.listdir(os.path.join(filtered_dir, fold)):
            print("file", file)
            tag=file.replace("_test_results.tsv", "")
            tag=tag.replace("_negative_results.tsv", "")
            if tag not in res:
                res[tag]={}
                res[tag]["group"]=fold
            #if file contains negative
            if "negative" in file:
                group_counts={}
                df=pd.read_csv(os.path.join(filtered_dir, fold,file), sep='\t')
                for index, row in df.iterrows():
                    group=extract_group(row)

                    if group not in group_counts:
                        group_counts[group]=0
                    group_counts[group]+=1
                    #get the average
                probs={}
                for group in group_counts:
                    if group is not None:
                        prop=group_counts.get(group,0)/group_countsd.get(group,1)
                        prob=prop*average_dict.get(group,0)
                        probs[group]=prob
                        print("Probability of at least 1 FP",group, prob)
                        res[tag][group+" prob of FP"]=prob
                    else:
                        prop=group_counts[group]/neg_genes
                        prob=1-pow(1-prop,av_genes) 
                        probs[group]=prob
                        print("Probability of at least 1 FP non hydrogenase", prob)
                        res[tag]["non hydrogenase prob of FP"]=prob
                #get the weighted average
                total=0
                total_weight=0
                if len(probs)==0:
                    FP_prob=0
                    print("Probability of at least 1 FP", FP_prob)
                    res[tag]["prob of FP"]=FP_prob
                    continue
                for group in probs:
                    
                    if group is not None:
                        weight= copy_dict.get(group, 0)
                        total+=probs.get(group,0)*weight
                        total_weight+=weight
                    else:
                        weight=av_genes
                        total+=probs[group]*weight
                        total_weight+=weight
                #if total weight is 0, set FP prob to 0
                if total_weight==0:
                    FP_prob=0
                    print("Probability of at least 1 FP", FP_prob)
                    res[tag]["prob of FP"]=FP_prob
                else:
                    FP_prob=total/total_weight
                    print("Probability of at least 1 FP", FP_prob)
                    res[tag]["prob of FP"]=FP_prob
            else:
                #"C:\Users\abel\Documents\marker_sens_spec\filtered\[FeFe]_Group_B\[FeFe]_Group_B_0.2_40_0.9_0_negative_results.tsv"
                #extract the database split, it is 0.2
                pattern = r'\d+\.\d+'
                matches = re.findall(pattern, file)

                # Extract the first match
                if matches:
                    first_float = float(matches[0])
                else:
                    first_float = 0
                test_percent=1-first_float
                

                pos_counts=0
                df=pd.read_csv(os.path.join(filtered_dir, fold,file), sep='\t')
                pos_counts=df.shape[0]
                prop=pos_counts/(group_countsd[fold]*test_percent)
                prob=prop*average_dict[fold]
                print("Probability of at least 1 TP", prob)
                res[tag]["prob of TP"]=prob
        #make the df
    all_keys=set()
    for key in res:
        all_keys.update(res[key].keys())
    rows=[]
    for key in res:
        row1=res[key]
        group=row1['group']
        tage=key.replace(group, "")
        spl=tage.split("_")[1:]
        pid=spl[1]
        database_split=spl[0]
        scov=spl[2]
        trial=spl[3]
        row={}
        row['pid']=pid
        row['database_split']=database_split
        row['scov']=scov
        row['trial']=trial
        row['group']=group
        row["FP"]=row1.get("prob of FP", 0)
        row["TP"]=row1.get("prob of TP", 0)
        row["non hydrogenase FP"]=row1.get("non hydrogenase prob of FP", 0)
        for key in all_keys:
            if key in row1:
                if key=="group":
                    continue
                row[key]=row1[key]
            else:
                row[key]=0
        rows.append(row)
    df=pd.DataFrame(rows)
    df["TP/FP"]=df["prob of TP"]/df["prob of FP"]
    #save the df
    df.to_csv(os.path.join(filtered_dir, "filtered_results.tsv"), sep='\t', index=False)


            
    return probs

root_dir_marker=r"C:\Users\abel\Documents\hydrogenases\marker_sens_spec\marker_results"
def get_marker_group_counts(root_dir):
    res={}
    for file in os.listdir(root_dir):   
        #if file contains results
        if "results" in file and file.endswith(".tsv"):
            df=pd.read_csv(os.path.join(root_dir, file), sep='\t')
            for index, row in df.iterrows():
                pid=row['pid']
                scov=row['scov']
                trial=row['trial']
                train_percent=row['train_percent']
                gene=row['gene']
                res[(gene,train_percent,pid,scov,trial)]=row["sensitivity"]
    return res
mouse_genomes=1676
human_genomes=3338
def make_marker_copies_genomes(root_dir):
    human_marker_genomes=os.path.join(root_dir, "marker_copies_genomes_human.tsv")
    mouse_marker_genomes=os.path.join(root_dir, "marker_copies_genomes_mouse.tsv")
    df_hum=pd.read_csv(human_marker_genomes, sep='\t')
    #divide each column by the number of human genomes
    
    df_hum["genomes"]=df_hum["genomes"]/human_genomes
    df_mouse=pd.read_csv(mouse_marker_genomes, sep='\t')
    #divide each column by the number of mouse genomes
    df_mouse["genomes"]=df_mouse["genomes"]/mouse_genomes
    #get the average of the two
    av = (df_hum["genomes"] + df_mouse["genomes"]) / 2
    # make the dict
    res = {}
    for index, value in av.iteritems():
        res[df_hum.at[index, 'gene']] = value
    return res
    #make the dict
    res={}
    for index, row in av.iterrows():
        res[row['gene']]=row['genomes']
    return res

def make_marker_copies(root_dir):
    human_marker_genomes=os.path.join(root_dir, "marker_copies_human.tsv")
    mouse_marker_genomes=os.path.join(root_dir, "marker_copies_mouse.tsv")
    df_hum=pd.read_csv(human_marker_genomes, sep='\t')
    df_mouse=pd.read_csv(mouse_marker_genomes, sep='\t')
    #divide by the number of human genomes
    df_hum["copies"]=df_hum["copies"]/human_genomes
    df_mouse["copies"]=df_mouse["copies"]/mouse_genomes
    #get the average of the two
    av=(df_hum["copies"]+df_mouse["copies"])/2
    #make the dict
    res={}
    #res = {}
    for index, value in av.iteritems():
        res[df_hum.at[index, 'gene']] = value
    return res




def parse_filtered2(filtered_dir, sens_dict,genomes_dict, copy_dict,neg_genes):
    res={}
    for fold in os.listdir(filtered_dir):
        #if fold is not a directory, continue   
        if not os.path.isdir(os.path.join(filtered_dir, fold)):
            continue
        cont_folder= os.path.join(filtered_dir, fold,"filtered")
        for file in os.listdir(cont_folder):
            #if file is not a tsv file, continue    
            if not file.endswith(".tsv"):
                continue
            print("file", file)
            tag=file.replace("_test_results.tsv", "")
            tag=tag.replace("_negative_results.tsv", "")
            if tag not in res:
                res[tag]={}
                res[tag]["group"]=fold.replace("filtered", "")
            #if file contains negative
            if "negative" in file:
                df=pd.read_csv(os.path.join(cont_folder,file), sep='\t')
                fp=len(df)
                
                prop=fp/neg_genes
                prob=1-pow(1-prop,copy_dict[None])
                
                print("Probability of at least 1 FP non marker", prob)
                
                res[tag]["prob of FP"]=prob
            else:
                #"C:\Users\abel\Documents\marker_sens_spec\filtered\[FeFe]_Group_B\[FeFe]_Group_B_0.2_40_0.9_0_negative_results.tsv"
                #extract the database split, it is 0.2
                pattern = r'\d+\.\d+'
                matches = re.findall(pattern, file)

                # Extract the first match
                if matches:
                    first_float = float(matches[0])
                else:
                    first_float = 0
                test_percent=1-first_float
                
                gene,train_percent,pid,scov,trial=tag.split("_")
                entry=(gene,float(train_percent),int(pid),float(scov),int(trial))
                if entry not in sens_dict:
                    raise Exception("Entry not in sensitivity dict {}".format(entry))
                pos_counts=0
                df=pd.read_csv(os.path.join(cont_folder,file), sep='\t')
                pos_counts=df.shape[0]
                prop=sens_dict[entry]*test_percent
                prob=prop*genomes_dict[gene]
                print("Probability of at least 1 TP", prob)
                res[tag]["prob of TP"]=prob
        #make the df
    all_keys=set()
    for key in res:
        all_keys.update(res[key].keys())
    rows=[]
    for key in res:
        row1=res[key]
        group=row1['group']
        if group=="mcra_":
            group="mcrA_"
        if group=="nfrA_":
            group="nrfA_"
        tage=key.replace(group, "")
        spl=tage.split("_")[:]
        #remove the first item if its empty
        if spl[0]=="":
            spl=spl[1:]
        pid=spl[1]
        database_split=spl[0]
        scov=spl[2]
        trial=spl[3]
        row={}
        row['pid']=pid
        row['database_split']=database_split
        row['scov']=scov
        row['trial']=trial
        row['group']=group
        row["prob of FP"]=row1.get("prob of FP", 0)
        row["prob of TP"]=row1.get("prob of TP", 0)
        for key in all_keys:
            if key in row1:
                if key=="group":
                    continue
                row[key]=row1[key]
            else:
                row[key]=0
        rows.append(row)
    df=pd.DataFrame(rows)
    df["TP/FP"]=df["prob of TP"]/df["prob of FP"]
    #save the df
    df.to_csv(os.path.join(filtered_dir, "filtered_results.tsv"), sep='\t', index=False)
            
import math
import numpy as np
def calculate_misclassification_probability(misclassification_probability, average_copies, odds_of_detection):
    probability_of_hydrogenase = odds_of_detection / (1 + odds_of_detection)
    probability_of_misclassification_given_hydrogenase = 1 - math.exp(-average_copies * misclassification_probability)
    probability_of_misclassification = (probability_of_hydrogenase *
                                        probability_of_misclassification_given_hydrogenase)
    return probability_of_misclassification
def calc_probs():
    human_table=r"C:\Users\abel\Downloads\human_1\human_\compare\human_results1_with_protein_names1.csv"
    mouse_table=r"C:\Users\abel\Downloads\human_1\human_\compare\mouse_results2_with_protein_names1.csv"
    human=pd.read_csv(human_table)
    mouse=pd.read_csv(mouse_table)
    h_copy=human.groupby("genome_name")["hyd_group"].value_counts().unstack()
    #fill nan with 0
    h_copy=h_copy.fillna(0)
    #get the proportion of non zero rows
    h_copy.loc["prop"]=h_copy.apply(lambda x: x[x>0].count()/len(x),axis=0)
    #get the average
    h_copy.loc["average"]=h_copy[:-1].mean()
    h_copy
    m_copy=mouse.groupby("genome_name")["hyd_group"].value_counts().unstack()
    m_copy=m_copy.fillna(0)
    m_copy.loc["prop"]=m_copy.apply(lambda x: x[x>0].count()/len(x),axis=0)
    m_copy.loc["average"]=m_copy[:-1].mean()
    gp=(h_copy.loc["prop"]+m_copy.loc["prop"])/2
    group_probs={}
    for index,row in gp.iteritems():
        sp=index.split(" ")
        if len(sp)==1:
            group_probs["["+sp[0]+"]"]=row
        group_probs["["+sp[0]+"]_Group_"+sp[1]]=row

    group_probs
    group_copies={}
    for index,row in h_copy.loc["average"].iteritems():
        sp=index.split(" ")
        if len(sp)==1:
            group_copies["["+sp[0]+"]"]=row
        group_copies["["+sp[0]+"]_Group_"+sp[1]]=row

    #parse the average file
    #group_probs, av_copy=parse_average_file(average_file)
    #group_counts=count_groups(r"C:\Users\abel\Documents\marker_sens_spec")
    #parse the negative file
    neg_genes= 570525#count_neg_genes(neg_file)
    print("Negative genes", neg_genes)
    #parse the filtered files
    #parse_filtered(r"C:\Users\abel\Documents\marker_sens_spec\filtered", group_counts, av_copy,group_probs,neg_genes)
    filt_file=r"C:\Users\abel\Documents\marker_sens_spec\hyd_results.tsv"
    def get_hyd_FP_prob(row):
        group=row['gene']
        sume=0
        for col in row.index:
            if "specificity" in col and "[" in col:
                spe=col.replace("specificity", "")
                #if the last char is a _, remove it
                if spe[-1]=="_":
                    spe=spe[:-1]
                group1=spe
                spec=row[col]
                prob=calculate_misclassification_probability(1-spec, group_copies.get(group1, 0), group_probs.get(group1, 0))
                #if nan, then zero
                if np.isnan(prob):
                    prob=0
                print(group1, prob)
                sume+=prob
        return sume
    av_genes=(av_mouse_genes+av_human_genes)/2
    def get_FP_prob(row):
        non_hyd=row['non_hyd_specificity']
        prob=calculate_misclassification_probability(1-non_hyd, av_genes, 1)
        print("Non hydrogenase", prob)
        return prob
    def get_TP_prob(row):
        group=row['gene']
        res=row['sensitivity']* group_probs.get(group, 0)
        #if is nan, set to 0
        if np.isnan(res):
            res=0
        return res
    df=pd.read_csv(filt_file, sep='\t')
    df['FP_hydrogenase']=df.apply(get_hyd_FP_prob, axis=1)
    df['FP_non_hydrogenase']=df.apply(get_FP_prob, axis=1)
    df["FP"]=df['FP_hydrogenase']+df['FP_non_hydrogenase']
    df['TP']=df.apply(get_TP_prob, axis=1)
    #save the df
    df.to_csv(filt_file.replace(".tsv", "_probs.tsv"), sep='\t', index=False)


    

    
def calc_probs_marker():
    
    root_dir=r'C:\Users\abel\Documents\hydrogenases\marker_sens_spec'
    marker_genomes=make_marker_copies_genomes(root_dir_marker)
    marker_copies=make_marker_copies(root_dir_marker)
    #8075099/3338 +4499713/1676
    marker_copies[None]=(8075099/3338 +4499713/1676)/2
    neg_genes=570525
    sens_dict=get_marker_group_counts(root_dir_marker)
    #parse_filtered2(filtered_dir, sens_dict,genocopy_dict,neg_genes)
    filtered_dir=r"C:\Users\abel\Documents\hydrogenases\marker_sens_spec\marker_filtered"
    parse_filtered2(filtered_dir, sens_dict,marker_genomes,marker_copies,neg_genes)


calc_probs()



#run the PFOR
marker_dir=r"C:\Users\abel\Documents\marker_sens_spec\final_curated_marker"
table_folder=r"C:\Users\abel\Documents\marker_sens_spec\tables"
table_file=r"C:\Users\abel\Documents\marker_sens_spec\tables\PFOR_table.tsv"
cwd_folder=r"C:\Users\abel\Documents\marker_sens_spec\PFOR"
file="PFOR.fasta"

#run_marker(os.path.join(marker_dir, file), table_file,cwd_folder)




aggregate_files_batched(r"C:\Users\abel\Documents\marker_sens_spec")






