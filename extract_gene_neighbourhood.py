from Bio import SeqIO
import re
#import the SeqRecord class from the Bio.SeqIO module
from Bio.SeqRecord import SeqRecord
import os
import re
default_diamond_exe_path=r"D:\diamond-windows\diamond.exe"
default_cdd_search_script_path=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\cdd_search.py"
manifest_file=r"C:\Users\abel\Documents\hydrogenases\hydrogenase_subunits\final_manifest.csv"
homolog_database_file=r'/tllhome/abel/hydrogen_consumer_pipeline/subunits/combined.fasta'
import csv
homolog_folder=r'C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\gene_organisation_dataset\subunits'
def build_homolog_database(manifest_file, out_dir):
    #if combined file exists and is not empty, return the combined file
    if os.path.exists(os.path.join(out_dir, "combined.fasta")) and os.path.getsize(os.path.join(out_dir, "combined.fasta"))>0:
        return os.path.join(out_dir, "combined.fasta")
    #build the homolog of FeFe subunit of hydrogenase in the manifest file
    with  open(manifest_file,"r",errors="ignore") as f:
        reader=csv.DictReader(f, delimiter="\t")
        for row in reader:
            file=row["file"]
            #remove quotes and whitespace
            file=file.replace('"', "").strip()
            gene_symbol=row["relevant"]
            group=row["group"]
            #if group contains FeFe_A, then copy the file to the out_dir, and rename it to the gene symbol + the group
            if "FeFe A" in group:
                new_name=os.path.join(out_dir, gene_symbol+"_"+group+".fasta")
                with open(file, "r") as f1:
                    content=f1.read()
                    with open(new_name, "w") as f2:
                        f2.write(content)
    #make a combined file of all the files in the out_dir
    combined_file=os.path.join(out_dir, "combined.fasta")
    with open(combined_file, "w") as f:
        for file in os.listdir(out_dir):
            for record in SeqIO.parse(os.path.join(out_dir, file), "fasta"):
                #change the name to the file name
                fn=file.split(".")[0]
                record.id=fn
                record.description=fn
                #write the record to the combined file
                SeqIO.write(record, f, "fasta")
    return combined_file

#if the homolog database dont exist, or is empty, build it
#if not os.path.exists(homolog_folder) or len(os.listdir(homolog_folder))==0:
#    os.makedirs(homolog_folder,exist_ok=True)
#    build_homolog_database(manifest_file, homolog_folder)
def sort_fn(x):
    #split the x by "_"
    x=x.split("_")
    #multiply the first element by 1 billion and add the second element
    return int(x[0])*1000000000+int(x[1])
def parse_genbank_file(file:str):
    with open(file, "r") as f:
        content=f.read()
        contigs={}
        contigs_str=content.split("//")
        for contig_str in contigs_str:
            #extract the seqnum=1
            seqnum=re.findall(r'seqnum=(\d+)', contig_str)
            #if len(seqnum)==0, continue
            if len(seqnum)==0:
                continue
            #add the seqnum to the contigs dictionary
            seqhdr=seqnum[0]
            #remove quotes
            seqhdr=seqhdr.replace('"', "")
            contigs[seqhdr]=[]
            #split by CDS
            cdss=contig_str.split("CDS")[1:]
            for cds in cdss:
                #extract the ID=2_1
                ids=re.findall(r'ID=(\d+_\d+)', cds)
                #split the id by _ and take the last element

                #if len(id)==0, continue
                if len(ids)==0:
                    continue
                idd=ids[0]
                idd1=idd.split("_")[-1]
                #add the id to the contigs dictionary
                contigs[seqhdr].append(idd1)
            #sort the contigs by the id
            contigs[seqhdr].sort(key=lambda x:int(x))
    return contigs

def get_neighbourhood(contigs, contig,id1, window_size=5):
    id1=id1.split("_")[-1]
    #get the index of the id1 in the contigs
    index=int(id1)-1
    #get the start and end of the window
    start=index-window_size
    end=index+window_size
    #if start<0, set start to 0
    if start<0:
        start=0
    #if end>len(contigs[contig]), set end to len(contigs[contig])
    if end>len(contigs[contig]):
        end=len(contigs[contig])
    #return the contigs from start to end
    return contigs[contig][start:end]

def parse_cds_file(cds_file):
    #open the file
    with open(cds_file, "r") as f:
        #parse the file using the SeqIO.parse function
        records = {}
        for record in SeqIO.parse(f, "fasta"):
            desc=record.description
            #extract the ID_x_x from the description
            regex=re.compile(r'ID=(\d+_\d+)')
            matches=regex.findall(desc)
            if len(matches)==0:
                raise ValueError(f"ID not found in {desc}")
            id12=matches[0]
            contig=id12.split("_")[0]
            id1=id12.split("_")[1]

            if contig not in records:
                records[contig]={}
            records[contig][id1]=record
    return records
def extract_contig_id(idd):
    #split idd by _
    idd=idd.split("_")
    #return the first element of the split
    contig="_".join(idd[:-1])
    id1=idd[-1]
    return contig,id1
def extract_gene_neighbourhood(row,tmp_folder="./genomic_neighbourhood"):
    file=row[0]
    contig_num=row[1]
    id1=row[2]
    #get the prodigal dir, it is the parent of the file
    prodigal_folder=os.path.dirname(file)
    
    cds_file=os.path.join(prodigal_folder, file)
    #parse the cds file
    records=parse_cds_file(cds_file)
    #open the file,parse it using gb (genbank) format
    gb_file=file+".txt"
    contigs={}
    contigs=parse_genbank_file(gb_file)
    #get the neighbourhood of the gene
    neighbourhood=get_neighbourhood(contigs, contig_num,id1)
    #get the records of the neighbourhood
    neighbourhood_records={}
    for ide in neighbourhood:

        neighbourhood_records[ide]=records[contig_num][ide]
    #write to the tmp file
    #if temp folder does not exist, create it
    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)
    #write the neighbourhood to the tmp folder, write one fasta file for all the records
    #get the name of the file from the filepath
    fn=os.path.basename(file)
    hydrogenase_gene_id=fn+"_"+contig_num+"_"+id1
    with open(os.path.join(tmp_folder, hydrogenase_gene_id+".fasta"), "w") as f:
        for record in neighbourhood_records:
            SeqIO.write(neighbourhood_records[record], f, "fasta")
            
    return os.path.join(tmp_folder, hydrogenase_gene_id+".fasta")


def run_diamond(diamond_exe_path, neighbourhood_file, homolog_folder, out_dir):
    #get the combined file
    combined_file=homolog_database_file
    #build_homolog_database(manifest_file, homolog_folder)
    #build the diamond database
    diamond_db=os.path.join(out_dir, "combined.dmnd")
    #if the diamond database does not exist, build it
    if not os.path.exists(diamond_db):
        os.system(f"{diamond_exe_path} makedb --in {combined_file} -d {diamond_db}")
    #run the diamond blast
    diamond_output=os.path.join(out_dir, "diamond_output.txt")
    os.system(f"{diamond_exe_path} blastp -d {diamond_db} -q {neighbourhood_file} -o {diamond_output} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore")
    return diamond_output
def get_prodigal_mappings(prodigal_folder):
    print("get_prodigal_mappings( ",prodigal_folder,")")
    id_mappings=defaultdict(dict)
    #get the faa file
    for file in os.listdir(prodigal_folder):
        #if file does not end with .faa, continue
        if not file.endswith(".faa"):
            continue
        #open the file
        with open(os.path.join(prodigal_folder, file), "r") as f:
            org_name=file
            recs=list(SeqIO.parse(f, "fasta"))
            for rec in recs:
                id1=rec.id
                #extratc the prodigal id , it is ID=X_X
                regex=re.compile(r'ID=(\d+_\d+)')
                matches=regex.findall(rec.description)
                if len(matches)==0:
                    print("ERROR", rec.description, "does not contain ID")
                    raise ValueError(f"ID not found in {rec.description}")
                id2=matches[0]
                id_mappings[org_name][id1]=id2
    return id_mappings
def compile_batch_hydrogenase_file(batch_folder2,tmpe_folder="./tmp",tmp_prefix="tmp"):
    #NOW it runs only on one batch folder
    compiled_rows=[]
    print("batch_folder2", batch_folder2)
    for batch_folder in [batch_folder2]:
        #batch root folder is the parent of the batch folder
        batch_root_dir=os.path.dirname(batch_folder)
        batch_folder=os.path.basename(batch_folder)
        print("batch_folder", batch_folder)
        print("batch_root_dir", batch_root_dir)
    #for batch_folder in os.listdir(batch_root_dir):
        #if folder does not start with batch, skip, and is not a directory, skip
        #if not batch_folder.startswith("batch") or not os.path.isdir(os.path.join(batch_root_dir, batch_folder)):
        #    print("skipping", batch_folder,"coz it does not start with batch or is not a directory")
        #    continue
        #get the tmp folder
        tmp_folder=os.path.join(batch_root_dir, batch_folder, tmp_prefix)
        print("tmp_folder", tmp_folder)
        prodigal_file=os.path.join(tmp_folder, "prodigal")
        print("prodigal_file", prodigal_file)
        #if prodigal folder does not exist, continue
        if not os.path.exists(prodigal_file):
            print("ERROR", prodigal_file, "does not exist")
            raise ValueError(f"Prodigal folder {prodigal_file} does not exist")
        id_mappings=get_prodigal_mappings(prodigal_file)
        #print("id_mappings", id_mappings)
        #get the all_hydrogenase.txt in the tmp folder
        if os.path.exists(os.path.join(tmp_folder, "all_hydrogenase.txt")):
            with open(os.path.join(tmp_folder, "all_hydrogenase.txt"), "r") as f:
                reader=csv.reader(f, delimiter="\t")
                for row in reader:
                    #if row[1] contains [FeFe] Group A, then add it to the compiled_rows
                    if "[FeFe]_Group_A" in row[1]:
                        id2=row[0].split("_file_")[1]
                        id2=id2.split(" ")[0]
                        prodigal_id=id_mappings[row[0].split("_file_")[0]][id2]
                        contig=prodigal_id.split("_")[0]
                        cds_no=prodigal_id.split("_")[1]
                        #add the path to the prodigal folder to the row
                        prodigal_folder=os.path.join(tmp_folder, "prodigal")
                        filepath= row[0].split("_file_")[0]
                        full_file=os.path.join(prodigal_folder, filepath)
                        row=[full_file, contig, cds_no,id2]
                        compiled_rows.append(row)
    '''all_recs={}
    new_compiled_rows=[]
    for row in compiled_rows:
        #extract the filename from the row[1]
        filename=row[1].split("_file_")[0]
        #add the prodigal folder to the filename
        filename=os.path.join(row[0], filename)
        sequence_id=row[1].split("_file_")[1].split(" ")[0]
        #open the prodigal file
        #if file is in all_recs, use it
        if filename in all_recs:
            records=all_recs[filename]
        else:
            with open(filename, "r") as f:
                records=list(SeqIO.parse(f, "fasta"))
                all_recs[filename]=records
        sequence_id2=None
        for record in records:
            id_desc=record.id+" "+record.description
            if sequence_id in id_desc:
                #extract the id from id_desc, it looks like ID=289_3
                matches=re.findall(r'ID=(\d+_\d+)', id_desc)
                if len(matches)>0:
                    sequence_id2=matches[0]
        #the new row is the filename, contig_num and sequence_id
        contig_num=sequence_id2.split("_")[0]
        sequence_id1=sequence_id2.split("_")[1]
        new_row=[filename, contig_num, sequence_id1]
        new_compiled_rows.append(new_row)
    compiled_rows=new_compiled_rows'''

            

    #write the compiled rows to a file
    #make a tmp folder in the batch_root_dir
    if not os.path.exists(tmpe_folder):
        os. makedirs(tmpe_folder)
    #write the compiled rows to a file
    with open(os.path.join(tmpe_folder, "compiled_hydrogenase.txt"), "w") as f:
        writer=csv.writer(f, delimiter="\t")
        writer.writerows(compiled_rows)
    return os.path.join(tmpe_folder, "compiled_hydrogenase.txt") 
from Bio import SeqIO
def batch_fileslist(fileslist,tmp_folder, batch_size=200):
    #batch the records in the fileslist such that each batch contains batch_size records
    batches=[]
    current_batch=[]
    for file in fileslist:
        #parse the file
        with open(file, "r") as f:
            records=list(SeqIO.parse(f, "fasta"))
        for record in records:
            current_batch.append(record)
            if len(current_batch)>=batch_size:
                batches.append(current_batch)
                current_batch=[]
    if len(current_batch)>0:
        batches.append(current_batch)
    #write the batches to the tmp_folder
    for i, batch in enumerate(batches):
        with open(os.path.join(tmp_folder, f"batch{i}.fasta"), "w") as f:
            SeqIO.write(batch, f, "fasta")
    return [os.path.join(tmp_folder, f"batch{i}.fasta") for i in range(len(batches))]
from collections import defaultdict
def batch_cdd_search(compiled_hydrogenase_file, tmp_folder, cdd_search_script_path,output_file=None):
    fileslist=[]
    #if tmp folder does not exist, create it
    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)
    #open the compiled_hydrogenase_file
    #make the manifests folder in the tmp_folder
    if not os.path.exists(os.path.join(tmp_folder, "manifests")):
        os.makedirs(os.path.join(tmp_folder, "manifests"))
    with open(compiled_hydrogenase_file, "r") as f:
        
        reader=csv.reader(f, delimiter="\t")
        for row in reader:

            #if len(row)<2, continue
            if len(row)<2:
                continue
            #extract the gene neighbourhood
            #filename=row[0]
            filename=os.path.basename(row[0])
            full_name=filename+"_"+row[1]+"_"+row[2]
            manifest_file=os.path.join(tmp_folder, "manifests",full_name+".txt")
            gene_neighbourhood_file=extract_gene_neighbourhood(row, tmp_folder)
            #add the filename to every record of the gene_neighbourhood_file
            #write the records_id to the manifest_file
            with open(manifest_file, "w") as f:
                for record in SeqIO.parse(gene_neighbourhood_file, "fasta"):
                    f.write(record.id+"\n")
            with open(gene_neighbourhood_file, "r") as f:
                records=list(SeqIO.parse(f, "fasta"))
            for record in records:
                filename=os.path.basename(row[0])
                record.id=filename+"@"+record.id
                #record.description=filename+"_"+record.description
            #write the records to the gene_neighbourhood_file
            with open(gene_neighbourhood_file, "w") as f:
                SeqIO.write(records, f, "fasta")
            #add the gene_neighbourhood_file to the fileslist
            fileslist.append(gene_neighbourhood_file)
    #batch the fileslist
    batches=batch_fileslist(fileslist,tmp_folder)
    rows=[]
    for file in batches:
        basename=os.path.basename(file)
        #remove the extension
        basename=basename.split(".")[0]
        #run the cdd, here is an example command
        #python batch_cd_search.py input.fasta temp_results.txt output_results.txt --dummy
        tmp_file1=os.path.join(tmp_folder, f"tmp_results_{basename}.csv")
        out_file=os.path.join(tmp_folder, f"results_{basename}.csv")
        os.system(f"python {cdd_search_script_path} {file} {tmp_file1} {out_file} ")
        #parse the out file to determine the FeFe type
        '''RULES:
         A sequence is classified as [FeFe] Group A2 if one of the domains “GltA”, “GltD”,
           “glutamate synthase small subunit” or “putative oxidoreductase”, but not “NuoF”, 
           is found in the sequence. Sequences are classified as [FeFe] Group A3 if the domain “NuoF” 
           is found and [FeFe] Group A4 if the domain “HycB” is present.
           If none of the domains are found, the sequence is classified as A1. '''
        
        with open(out_file, "r") as f:
            reader=csv.reader(f, delimiter="\t")
            for row in reader:
                row[0]=row[0].split("@")[1]
                rows.append(row)
    def group_by_gene(manifest_folder, rows):
    
        groups=defaultdict(list)
        for file in os.listdir(manifest_folder):
            if not file.endswith(".txt"):
                continue
            with open(os.path.join(manifest_folder, file), "r") as f:
                records=f.read().split("\n")
            for record in records:
                if record=="":
                    continue
                recs=[row for row in rows if record in row[0]]
                groups[file].extend(recs)
        return groups
    def classify_group(group_rows):
        classification=""
        for row in group_rows:
            #follow the rules to classify the group
            #if family contains GltA, GltD, glutamate synthase small subunit, or putative oxidoreductase, but not NuoF, then it is Group A2
            if "GltA" in row[1] or "GltD" in row[1] or "glutamate synthase small subunit" in row[1] or "putative oxidoreductase" in row[1]:
                if "NuoF" not in row[1]:
                    classification="Group A2"
            #if family contains NuoF, then it is Group A3
            if "NuoF" in row[1]:
                classification="Group A3"
            #if family contains HycB, then it is Group A4
            if "HycB" in row[1]:
                classification="Group A4"
        #if none of the above, then it is Group A1
        if classification=="":
            classification="Group A1"
        #classification="Group A11"
        return classification
    groups=group_by_gene(os.path.join(tmp_folder, "manifests"), rows)
    
    classification_result={}
    for group in groups:
        classification_result[group]=classify_group(groups[group])
    #write the rows to the out_file
    with open(output_file, "w") as f:
        writer=csv.writer(f, delimiter="\t")
        for classification in classification_result:
            #remove the txt
            classification1=classification.replace(".txt", "")
            writer.writerow([classification1, classification_result[classification]])
    return output_file
    
                
    
def fix_hydrogenases_file(batch_folder, classification_result,tmp_prefix="tmp"):
    print("fix hydrogenases file", batch_folder, classification_result,tmp_prefix)
    #open the classification_result
    classification_dict={}
    with open(classification_result, "r") as f:
        reader=csv.reader(f, delimiter="\t")
        for row in reader:
            #if len of row is less than 2, continue
            if len(row)<2:
                continue
            file=row[0].split(".faa"    )[0]
            #add back the faa
            file=file+".faa"
            id=row[0].split(".faa")[1]
            #remove _ from the front of the id
            id=id[1:]
            if  file not in classification_dict:
                classification_dict[file]={}
            classification_dict[file][id]=row[1]
    #find the all_hydrogenase.txt in all_batch_folder
            
    for fold in [batch_folder]:
        print("fold in fix_hydrogenases_file", fold)
        #all batch folder is the parent of the fold
        all_batch_folder=os.path.dirname(batch_folder)
        fold=os.path.basename(batch_folder)
        
        
        #get the prodigal folder
        prodigal_folder=os.path.join(all_batch_folder, fold, tmp_prefix, "prodigal")
        print("prodigal_folder", prodigal_folder)
        #if prodigal folder does not exist, continue
        if not os.path.exists(prodigal_folder):
            print("ERROR", prodigal_folder, "does not exist")
            continue


        id_mappings=get_prodigal_mappings(prodigal_folder)
        #get  the hydrogenase file
        hydrogenase_file=os.path.join(all_batch_folder, fold, tmp_prefix, "all_hydrogenase.txt")
        #if file not exists, continue
        if not os.path.exists(hydrogenase_file):
            continue
        new_rows=[]
        #open the hydrogenase_file
        with open(hydrogenase_file, "r") as f:
            reader=csv.reader(f, delimiter="\t")
            for row in reader:
                file=row[0].split("_file_")[0]
                protein_id=row[0].split("_file_")[1]
                protein_id=protein_id.split(" ")[0]
                prodigal_id=id_mappings[file][protein_id]

                #find the entries in classification_dict with the file
                if file in classification_dict:
                    for id1 in classification_dict[file]:
                        if id1 ==prodigal_id:
                            orig_class= row[1].split("_-_")[-1]
                            new_class=classification_dict[file][id1]
                            #add [FeFe] to the new_class
                            new_class="[FeFe] "+new_class
                            #replace space with _
                            new_class=new_class.replace(" ", "_")
                            #replace the orig_class with the new_class
                            row[1]=row[1].replace(orig_class, new_class)
                            print("replaced", row[0],orig_class, "with", new_class)
                new_rows.append(row)
        #the new file  is the old file but with a _fixed suffix
        new_file=hydrogenase_file.replace(".txt", "_fixed.txt")
        print("NEW_FILE", new_file)
        #write the rows to the hydrogenase_file
        with open(new_file, "w") as f:
            writer=csv.writer(f, delimiter="\t")
            writer.writerows(new_rows)
        return new_file




                
    



batches_root_folder=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\batches"
def teset():

    prodigal_folder=r'C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\test_gene_organization\batch1\tmp\prodigal'
    hydrogenase_gene_id="GCF_000179555.1_ASM17955v1_genomic_prodigal.faa_file_NZ_AECZ01000043.1_7"
    tmp_folder=r'C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\test_gene_organization\batch1\tmp\gene_neighbourhood'
    neighbourhood_file=extract_gene_neighbourhood(hydrogenase_gene_id, prodigal_folder,tmp_folder=tmp_folder)
    diamond_exe_path=default_diamond_exe_path
    out_dir=tmp_folder
    diamond_output=run_diamond(diamond_exe_path, neighbourhood_file, homolog_folder, out_dir)
def teset2():
    root_dir=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\test_run\batches\batch_0"
    tmp_folder=os.path.join(root_dir, "neighbourhoods_tmp")
    compiled_hydrogenase_file=compile_batch_hydrogenase_file(root_dir, tmp_folder)
    
    cdd_search_script_path=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\cdd_search.py"
    result=batch_cdd_search(compiled_hydrogenase_file, tmp_folder, cdd_search_script_path,output_file=os.path.join(tmp_folder, "all_results.csv"))
    #fix the hydrogenases file
    fix_hydrogenases_file(root_dir, result)


def main(root_dir, args=None):
    #the cdd search script path is in the same parent dir as this script
    parent_dir=os.path.dirname(os.path.realpath(__file__))
    cdd_search_script_path=os.path.join(parent_dir, "cdd_search.py")
    print(f"extract_gene_neighbourhood.main( {root_dir},{args})")
    #if cdd search script path is None, use the default
    #it is in the same dir as this script
    this_dir=os.path.dirname(os.path.realpath(__file__))
    if cdd_search_script_path is None:
        
        cdd_search_script_path=os.path.join(this_dir, "cdd_search.py")
        print("cdd search script path", cdd_search_script_path)
    #compile the batch hydrogenase file
    tmp_folder=os.path.join(root_dir, "neighbourhoods_tmp")
    #if tmp folder does not exist, create it
    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)
    compiled_hydrogenase_file=compile_batch_hydrogenase_file(root_dir, tmpe_folder=tmp_folder, tmp_prefix=args.tmp_dir)
    print("compiled_hydrogenase_file", compiled_hydrogenase_file)
    
    output_file=os.path.join(tmp_folder, "all_results.csv")
    #run the batch cdd search
    result=batch_cdd_search(compiled_hydrogenase_file, tmp_folder, cdd_search_script_path,output_file=output_file)
    print("result", result)
    #fix the hydrogenases file
    result=fix_hydrogenases_file(root_dir, result,tmp_prefix=args.tmp_dir)
    print("EXTRACT GENE NEIGHBOURHOOD DONE",result)
    return result
import argparse
#teset2() 
#exit()
if __name__=="__main__":
    
    #make the argument parser
    args=argparse.ArgumentParser(description="Extract gene neighbourhood, run diamond blast, and classify the hydrogenases.")
    args.add_argument("--root_dir", help="The root directory of the batch folders.")
    #the default dir is the dir of this script 
    default_dir=os.path.dirname(os.path.realpath(__file__))
    default_cdd_search_script_path=os.path.join(default_dir, "cdd_search.py")
    args.add_argument("--cdd_search_script_path", help="The path to the cdd_search.py script.",default=default_cdd_search_script_path)
    args=args.parse_args()
    main(args.root_dir, args.cdd_search_script_path)




        
