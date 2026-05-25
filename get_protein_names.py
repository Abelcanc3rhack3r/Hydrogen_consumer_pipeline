

import requests
import re
def fetch_protein_names(accession_number):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{accession_number}"
    response = requests.get(url)
    
    if response.status_code == 200:
        protein_data = response.json()
        #protein_data['protein']['submittedName'][0]['fullName']['value']
        #if recommendedx name in the protein data
        if 'recommendedName' in protein_data['protein']:
            protein_name = protein_data['protein']['recommendedName']['fullName']['value']
        elif 'submittedName' in protein_data['protein']:
            protein_name = protein_data['protein']['submittedName'][0]['fullName']['value']
        else:
            return "None"
            raise ValueError("No protein name found")
        return str(protein_name)
    else:
        return "None"
def tes():
    # Fetch protein names for each accession number in the 'protein match' column
    name1=fetch_protein_names('A0A3A8V0M6')
    name1
    print(name1)
    
    
def extract_name(row):
        #print the row columns
        #print("ROW COLUMNS:",row.keys())
        #print("ROW:",str(row))  
        protein_match=row["protein_match"]
        #if protein match is nan, return ""
        if pd.isnull(protein_match):
            return ""
        pattern = r"\{\{tr\|([A-Za-z0-9]+)\|[A-Za-z0-9_]+\}\}"

        # Use re.search to find the first occurrence of the pattern in the text
        match = re.search(pattern, protein_match)

        if match:
            accession_number = match.group(1)
            print("Accession Number:", accession_number)
            return fetch_protein_names(accession_number)
        else:
            return protein_match

import pandas as pd
def run():
    table_file=r"C:\Users\abel\Documents\mouse_gut\mouse_gut_catalogue\mgnify_,mouse_catalogue\new_results\hq_genomes3.csv"
    df = pd.read_csv(table_file, sep="\t")
    for index, row in df.iterrows():
        try:
            df.loc[index, "protein_name"] = extract_name(row)
        except Exception as e:
            print("Error fetching protein name for row", index)
    #save as csv
    df.to_csv(r"C:\Users\abel\Documents\mouse_gut\mouse_gut_catalogue\mgnify_,mouse_catalogue\new_results\hq_genomes_with_protein_names3.csv", index=False)
    df
def process(file):
    df = pd.read_csv(file, sep=",")
    #print the columns
    #print("COLUMNS NAMES:",df.columns)
    for index, row in df.iterrows():
        #try:
            df.loc[index, "protein_name"] = extract_name(row)
        #except Exception as e:
        #    print("Error fetching protein name for row", index,e)
    #save as csv
    df.to_csv(file.replace(".csv","_with_protein_names.csv"), index=False)
    return df