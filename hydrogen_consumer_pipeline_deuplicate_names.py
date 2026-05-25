import os

from Bio import SeqIO




def deduplicate_names(folder):
    #rename headers that have duplicate names in other fasta files
    #get all the fasta filee
    unique_headers={}
    for file in os.listdir(folder):
        print("processing",file)
        if file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fna") or file.endswith(".fsa")\
        or file.endswith(".fas"):
            #get the headers
            replace_headers={}
            recs=list(SeqIO.parse(os.path.join(folder,file),'fasta'))
            for rec in recs:
                #print(rec.id)
                #if header is not in unique_headers, add it
                if rec.id not in unique_headers:
                    unique_headers[rec.id]=file
                #else, rename the header
                else:
                    #print("Duplicate header",rec.id,"in",file,"and",unique_headers[rec.id])
                    #get the new header, it is the old header witha  suffix 0,1,2,3,4,5,6,7,8,9...
                    suffix=2
                    while rec.id+str(suffix) in unique_headers:
                        suffix+=1
                    new_header=rec.id+str(suffix)
                    #add the new header to the replace_headers
                    replace_headers[rec.id]=new_header
                    #add the new header to the unique_headers
                    #print("Renaming",rec.id,"to",new_header,"in",file)
                    unique_headers[new_header]=file
                    rec.id=new_header
            SeqIO.write(recs,os.path.join(folder,file),'fasta')
    print("Done!")

if __name__ =="__main__":
    folder=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\genomes\isolates"
    deduplicate_names(folder)
    print("Done!")