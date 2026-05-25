import requests
import argparse
import time
from typing import Dict

def submit_search(filename: str, search_params: Dict) -> str:
    """Submit search to NCBI Batch CD-Search and return Request-ID (RID)."""
    url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
    files = {'queries': open(filename, 'rb')}
    response = requests.post(url, files=files, data=search_params)
    response.raise_for_status()
    
    rid_line = next((line for line in response.text.splitlines() if line.startswith("#cdsid")), None)
    if rid_line:
        return rid_line.split()[1]
    else:
        raise ValueError("Submitting the search failed, can't make sense of response.")

def check_search_status(rid: str) -> None:
    """Check the search status until it is completed."""
    url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
    params = {'cdsid': rid, 'tdata': 'hits'}
    while True:
        time.sleep(5)  # Wait for 5 seconds before checking the status again
        response = requests.post(url, data=params)
        response.raise_for_status()
        
        if "#status\t0" in response.text:
            print("Search has been completed.")
            break
        #this is the value of response '#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\tQM3-qcdsearch-19B1505E7575C03A-307C5DAA8D89455D\n#datatype\thitsConcise Results\n#status\t3\tmsg\tJob is still running\n'

        elif "#status\t3" in response.text:
            print("Search is still running, please wait.")
        else:
            print("Search failed with the following response:")
            #print(response.text)
            break
        # Add checks for other statuses if needed

def retrieve_results(rid: str, search_params: Dict) -> str:
    """Retrieve search results using the provided RID."""
    url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
    params = search_params.copy()
    params.update({'cdsid': rid})
    
    response = requests.post(url, data=params)
    response.raise_for_status()
    
    return response.text

def process_file(file, search_params: Dict, tmp_file,out_file,dummy_mode=False) :
    """Process a list of FASTA files and return a dictionary with their CD-Search results."""
    
    try:
            if dummy_mode:
                print("FILE:", file)
                out_file=make_dummy_results(file, out_file)
                return out_file
            filename = file
            print(f"Submitting {filename}")
            #rid = 'QM3-qcdsearch-3DA708068C359BA9-270413FC5AF7181D'
            rid=submit_search(filename, search_params)
            check_search_status(rid)
            result = retrieve_results(rid, search_params)
            out_file= parse_results(result, tmp_file,out_file)
            return out_file
            
    except Exception as e:
            raise e
    
import csv
import re
from collections import defaultdict
from Bio import SeqIO
import random
def make_dummy_results(fasta_file, out_file):
    possible_domains= ["GltA", "GltD", "glutamate synthase small subunit", "putative oxidoreductase", "NuoF", "HycB"]
    print(f"Generating dummy results for {fasta_file}")
    with open(fasta_file, "r") as f:
        records = list(SeqIO.parse(f, "fasta"))
    with open(out_file, "w") as f:
        for record in records:
            #randomly select a domain
            fam1=random.choice(possible_domains)
            fam2=random.choice(possible_domains)
            print("recid:", record.id, "fam1:", fam1, "fam2:", fam2)
            f.write(f"{record.id}\t{fam1}\t{fam2}\n")
    return out_file

def parse_results(result: str, tmp_file,out_file,max_domains=5) -> defaultdict:
    print("CDD parse results:",out_file)
    """Parse the CD-Search results and return a list of hits."""
    #save it as  a csv file
    #tmp_file="tmp.csv"
    with open(tmp_file, "w") as f:
        f.write(result)
    #open the file and read the content
    query_found=False
    res=defaultdict(list)
    with open(tmp_file, "r") as f:
        reader=csv.reader(f, delimiter="\t")
        for row in reader:
            #if len of row is 0, skip it
            if len(row)==0:
                continue
            #if row[0] is Query then we have found the start of the query
            if row[0]=="Query":
                query_found=True
                continue
            if query_found:
                #if len of row is less than 9 then we have reached the end of the query
                if len(row)<9:
                    query_found=False
                    continue
                query=row[0]
                #extract the gene{X} from the query
                #regex=re.compile(r"gene-\d+")
                gene=query.split("- >")[1].strip()
                
                #get the column I
                I=row[8]
                family=I
                print(gene, "is a member of family", family)
                res[gene].append(family)
    #write the results to a file
    print("WRITING THE OUTPUT TO", out_file)
    with open(out_file, "w") as f:
        for gene, family in res.items():
            fam_str="\t".join(family)[:max_domains]
            f.write(f"{gene}\t{fam_str}\n")
    return out_file




    #parse the content
    
        


def main():
    parser = argparse.ArgumentParser(description="Batch CD-Search with multiple FASTA files.")
    parser.add_argument('file', help='FASTA file to process.')
    #add the tmp file and out file
    parser.add_argument('tmp_file', help='Temporary file to store the results.')
    parser.add_argument('out_file', help='Output file to store the results.')
    #add the dummy mode
    parser.add_argument('--dummy', action='store_true', help='Use dummy mode to generate fake results.')

    args = parser.parse_args()
    
    # Define default search parameters
    search_params = {
        'dmode': 'rep',
        'clonly': 'false',
        'useid1': 'true',
        'maxhit': '250',
        'filter': 'true',
        'db': 'cdd',
        'evalue': '0.01',
        'tdata': 'hits',
    }
    file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\gene_organisation_dataset\cdd_tmp\batch0.faa"
    tmp_file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\gene_organisation_dataset\genomes\Alkaliphilus metalliredigens\neighbourhoods\[FeFe] Group A3_1_forward\tmp\cdd_tmp.csv"
    out_file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\gene_organisation_dataset\genomes\Alkaliphilus metalliredigens\neighbourhoods\[FeFe] Group A3_1_forward\tmp\cdd_search_results.csv"
    results = process_file(args.file, search_params, args.tmp_file,args.out_file,args.dummy)
    tmp_file1=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\gene_organisation_dataset\cdd_tmp\test_tmp1.csv"
    out_file1=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\gene_organisation_dataset\cdd_tmp\test_out1.csv"
    #results = process_file(file, search_params, tmp_file1,out_file1,False)
    

if __name__ == "__main__":
    #file=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\hydrogen_consumer_pipeline\gene_organisation_dataset\genomes\Alkaliphilus metalliredigens\neighbourhoods\[FeFe] Group A3_1_forward\tmp\combined.faa"
    main()