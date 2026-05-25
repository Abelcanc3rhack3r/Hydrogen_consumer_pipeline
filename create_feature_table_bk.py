'''create a feature table of
genome_name, is_hydrogen_consumer, phylum, class, order, family, genus, species, hydrogenase_groups, is_sulphate_reducer,is_acetogen,is_methanogen

'''
import argparse

hydrogen_consumer_pipeline_output_folder=r'C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\mouse_catalogues'
#hydrogen_classes_activity_file= r'C:\Users\abel\Documents\hydrogenases\hydrogenase_classes_activity.csv'
hydrogen_classes_activity_file= r"C:\Users\abel\Documents\hydrogen_consumers\hydrogenase\hydrogenase_classes_activity2.csv"
taxon_output=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\out.BAT.bin2classification_name.txt"

import csv
import re
import os
def combine_batches(hydrogen_consumer_pipeline_output_folder,out_file):
    #combine the tsvs in the batch folder into one tsv
    #get the batch folder
    batch_folder=hydrogen_consumer_pipeline_output_folder
    #get the files in the folder
    files=os.listdir(batch_folder)
    #get the tsv files
    tsv_files=[file for file in files if file.endswith('.tsv')]
    rows=[]
    header=None
    #open the first file
    for file in tsv_files:
        #open the file
        with open(os.path.join(batch_folder,file),'r') as f:
            reader=csv.reader(f,delimiter='\t')
            #get the header
            header=next(reader)
            #add the rows
            for row in reader:
                rows.append(row)
    #write the rows to the out file
    with open(out_file,'w') as f:
        writer=csv.writer(f,delimiter='\t')
        #write the header
        writer.writerow(header)
        #write the rows
        for row in rows:
            writer.writerow(row)

    return out_file
def make_classes_activity_dict(hydrogenase_classes_atcivity_file):
    #the data looks like this:
    '''
Class	Group	direction	oxygen tolerance	localization	Function
FeFe	A1	H2-evolution 	O2-labile 	Cytosolic 	Fermentative or photobiological evolution of H2�using reduced ferredoxin as an electron donor.
FeFe	A2	H2-uptake ?	O2-labile 	Cytosolic 	Unconfirmed role. It may couple H2�oxidation to NAD reduction, generating reductant for glutamate synthase.
FeFe	A3	Electron-bifurcation 	O2-labile 	Cytosolic 	Bifurcates electrons from H2�to NAD and ferredoxin. Reverse reaction can also occur resulting in reoxidation of redox carriers concomitant with fermentative H2�evolution.

    '''
    #make a dictionary of the (classes ,group)and their direction, localization, oxygen tolerance

    classes_activity_dict={}
    with open(hydrogenase_classes_atcivity_file,'r') as f:
        reader=csv.reader(f,delimiter=',')
        #get the header
        header=next(reader)
        #get the data
        for row in reader:
            class2=row[0].strip()
            group=row[1].strip()
            direction=row[2].strip()
            oxygen_tolerance=row[3].strip()
            localization=row[4].strip()
            function=row[5].strip()
            classes_activity_dict[(class2,group)]={'direction':direction,'oxygen_tolerance':oxygen_tolerance,'localization':localization,'function':function}
    return classes_activity_dict

def extract_hydrogenase_activity(cell, classes_activity_dict):
    #an example of a cell looks like this:
    #FeFe-WP_028027642.1 - Enterorhabdus mucosicola - [FeFe] Group A2
    #in this case the class is Fefe and the group is A2
    #get the class
    class_=cell.split('[')[1].split(']')[0]
    #if class is Fe, then return classes_activity_dict[(Fe,Hmd)]
    if class_=='Fe':
        return classes_activity_dict[('Fe','Hmd')]
    #get the group
    group=cell.split('Group ')[1]
    #refer to the hydrogenase classes activity file to get the direction
    direction=classes_activity_dict[(class_,group)]['direction'].strip()
    localization=classes_activity_dict[(class_,group)]['localization'].strip()
    oxygen_tolerance=classes_activity_dict[(class_,group)]['oxygen_tolerance'].strip()
    return direction,localization,oxygen_tolerance

def parse_combined_hyd_file(combined_hyd_tsv,hydrogen_classes_activity_file=hydrogen_classes_activity_file):
    #parse the combined tsv file into a dictionary
    #the data looks like this:
    '''
    filename	hydrogen_consumer	sulphate_reducer	methanogen	acetogen	consuming hydrogenase	bifurcating_or_bidirectional_hydrogenase	dsr_matches	mcr_matches	acetoscan_matches	hyDB_matches
MGBC000218.fna	False	False	False	False	False	False

MGBC000217.fna	True	False	False	False	True	False
										FeFe-WP_028027642.1 - Enterorhabdus mucosicola - [FeFe] Group A2
										NiFe-WP_028026422.1 - Enterorhabdus mucosicola - [NiFe] Group 1i
MGBC000219.fna	False	False	False	False	False	False

MGBC000344.fna	True	False	False	False	True	False
										NiFe-WP_013010514.1 - Denitrovibrio acetiphilus - [NiFe] Group 1c

    '''
    classes_activity_dict = make_classes_activity_dict(hydrogen_classes_activity_file)
    #open the file
    with open(combined_hyd_tsv,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        #get the header
        header=next(reader)
        #get the data
        current_genome=None
        data={}
        for row in reader:
            #if the first col is not empty, then it is a new genome
            #if len of row is 0,skip
            if len(row)==0:
                continue
            if row[0]!='':
                current_genome=row[0]
                data[current_genome]={}
                data[current_genome]['hydrogen_consumer']=row[1]
                data[current_genome]['sulphate_reducer']=row[2]
                data[current_genome]['methanogen']=row[3]
                data[current_genome]['acetogen']=row[4]
                data[current_genome]['consuming_hydrogenase']=row[5]
                data[current_genome]['bifurcating_or_bidirectional_hydrogenase']=row[6]
                data[current_genome]['dsr_matches']=[]
                data[current_genome]['mcr_matches']=[]
                data[current_genome]['acetoscan_matches']=[]
                data[current_genome]['hyDB_matches']=[]
            #if the first col is empty, then it is a hydrogenase
            else:
                #check the row[] to see ifthe dsr match is not empty
                if len(row) >7 and row[7]!='':
                    data[current_genome]['dsr_matches'].append(row[7])
                if len(row) >8 and row[8]!='':
                    data[current_genome]['mcr_matches'].append(row[8])
                if len(row) >9 and row[9]!='':
                    data[current_genome]['acetoscan_matches'].append(row[9])
                if len(row) >10 and row[10]!='':
                    #get the direction, localization, oxygen tolerance

                    direction,localization,oxygen_tolerance=extract_hydrogenase_activity(row[10],classes_activity_dict)
                    data[current_genome]['hyDB_matches'].append((row[10],direction,localization,oxygen_tolerance))
    return data

def make_row_data(genome_name,datadict,taxon_dict):

    datadict_row = []
    datadict1=datadict
    # for each data in the data dict
    for data in ['hydrogen_consumer', 'sulphate_reducer', 'methanogen', 'acetogen', 'consuming_hydrogenase',
                 'bifurcating_or_bidirectional_hydrogenase']:
        # add the data to the datadict row
        datadict_row.append(datadict1[data])

    # get the taxon dict1
    taxon_dict1 = taxon_dict[genome_name]
    taxon_rows = []
    # for each taxon level
    for taxon_level in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        # if the taxon level is in the taxon dict1
        if taxon_level in taxon_dict1:
            # add the taxon name to the taxon rows
            taxon_rows.append(taxon_dict1[taxon_level])
        # if the taxon level is not in the taxon dict1
        else:
            # add NA to the taxon rows
            taxon_rows.append('Unknown')
    return datadict_row,taxon_rows

def make_feature_table(combined_hyd_tsv, out_file, taxon_output=taxon_output, hydrogen_classes_activity_file=hydrogen_classes_activity_file):
    #open the taxon file and save the taxon to a dictionary
    #the data looks like this:
    '''
    # bin	classification	reason	lineage	lineage scores	full lineage names	SuperKingdom	Kingdom	SuperPhylum	Phylum	Class	Class	1	2	3
MGBC000001.fasta	taxid assigned	based on 5551/5683 ORFs	1;131567;2;1783272;1239;186801;186802;186803;186928;1898203	1.00;1.00;1.00;0.94;0.94;0.91;0.91;0.88;0.50;0.39	root (no rank): 1.00	cellular organisms (no rank): 1.00	Bacteria (superkingdom): 1.00	Terrabacteria group (clade): 0.94	Firmicutes (phylum): 0.94	Clostridia (class): 0.91	Clostridiales (order): 0.91	Lachnospiraceae (family): 0.88	unclassified Lachnospiraceae (no rank): 0.50	Lachnospiraceae bacterium (species): 0.39
MGBC000003.fasta	taxid assigned	based on 3663/3775 ORFs	1;131567;2;1783272;1239;186801;186802;186803	1.00;1.00;1.00;0.80;0.79;0.59;0.58;0.34	root (no rank): 1.00	cellular organisms (no rank): 1.00	Bacteria (superkingdom): 1.00	Terrabacteria group (clade): 0.80	Firmicutes (phylum): 0.79	Clostridia (class): 0.59	Clostridiales (order): 0.58	Lachnospiraceae (family): 0.34
MGBC000011.fasta	taxid assigned	based on 4388/4437 ORFs	1;131567;2;1783270;68336;976;200643;171549;815;816	1.00;0.99;0.99;0.87;0.87;0.87;0.87;0.87;0.74;0.73	root (no rank): 1.00	cellular organisms (no rank): 0.99	Bacteria (superkingdom): 0.99	FCB group (clade): 0.87	Bacteroidetes/Chlorobi group (clade): 0.87	Bacteroidetes (phylum): 0.87	Bacteroidia (class): 0.87	Bacteroidales (order): 0.87	Bacteroidaceae (family): 0.74	Bacteroides (genus): 0.73

    '''

    taxon_dict={}
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
                taxon_level=re.findall('\((.*?)\)',row[col])
                #if there is a taxon level
                if taxon_level:
                    #get the taxon name
                    taxon_name=row[col].split(':')[0]
                    #remove the words in the brackets using re
                    taxon_name=re.sub('\(.*?\)','',taxon_name)
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
    #parse the combined_hyd_tsv file
    parsed_dict=parse_combined_hyd_file(combined_hyd_tsv,hydrogen_classes_activity_file=hydrogen_classes_activity_file)
    #open the out file
    with open(out_file,'w') as f:
        #the result is in long format, with one row per dsr, mcr, acetoscan or hydrogenase protein
        #write the header
        #the header is :
        #genome_name, is_hydrogen_consumer, is_sulphate_reducer, is_acetogen, is_methanogen, consuming_hydrogenase, bifurcating or bidirectional hydrogenase,
        #kingdom, phylum, class, order, family, genus, species,
        #protein match, hyd_direction, hyd_localization, hyd_oxygen_tolerance (NA if protein match is not hydrogenase)
        writer=csv.writer(f,delimiter='\t')
        header=['genome_name','is_hydrogen_consumer','is_sulphate_reducer','is_acetogen','is_methanogen','consuming_hydrogenase','bifurcating_or_bidirectional_hydrogenase','kingdom','phylum','class','order','family','genus','species','protein_match','hyd_direction','hyd_localization','hyd_oxygen_tolerance']
        writer.writerow(header)
        #for each genome in the parsed dict
        for genome in parsed_dict:
            #get the data dict
            data_dict=parsed_dict[genome]
            #for each dsr in the genome
            for dsr in data_dict['dsr_matches']:
                hyd_direction='NA'
                hyd_localization='NA'
                hyd_oxygen_tolerance='NA'
                datadict_row, taxon_rows=make_row_data(genome_name=genome,datadict=data_dict,taxon_dict=taxon_dict)
                #add the genome name, datadict row, taxon rows, dsr, dsr direction, dsr localization, dsr oxygen tolerance
                writer.writerow([genome]+datadict_row+taxon_rows+[dsr,hyd_direction,hyd_localization,hyd_oxygen_tolerance])
            #for each mcr in the genome
            for mcr in data_dict['mcr_matches']:
                hyd_direction='NA'
                hyd_localization='NA'
                hyd_oxygen_tolerance='NA'
                datadict_row, taxon_rows=make_row_data(genome_name=genome,datadict=data_dict,taxon_dict=taxon_dict)
                #add the genome name, datadict row, taxon rows, mcr, mcr direction, mcr localization, mcr oxygen tolerance
                writer.writerow([genome]+datadict_row+taxon_rows+[mcr,hyd_direction,hyd_localization,hyd_oxygen_tolerance])
            #for each acetoscan in the genome
            for acetoscan in data_dict['acetoscan_matches']:
                hyd_direction='NA'
                hyd_localization='NA'
                hyd_oxygen_tolerance='NA'
                datadict_row, taxon_rows=make_row_data(genome_name=genome,datadict=data_dict,taxon_dict=taxon_dict)
                #add the genome name, datadict row, taxon rows, acetoscan, acetoscan direction, acetoscan localization, acetoscan oxygen tolerance
                writer.writerow([genome]+datadict_row+taxon_rows+[acetoscan,hyd_direction,hyd_localization,hyd_oxygen_tolerance])
            #for each hydrogenase in the genome
            for hydrogenase in data_dict['hyDB_matches']:
                #get the hydrogenase direction, localization and oxygen tolerance
                hyd_name=hydrogenase[0]
                hyd_direction=hydrogenase[1]
                hyd_localization=hydrogenase[2]
                hyd_oxygen_tolerance=hydrogenase[3]
                datadict_row, taxon_rows=make_row_data(genome_name=genome,datadict=data_dict,taxon_dict=taxon_dict)
                #add the genome name, datadict row, taxon rows, hydrogenase, hydrogenase direction, hydrogenase localization, hydrogenase oxygen tolerance
                writer.writerow([genome]+datadict_row+taxon_rows+[hyd_name,hyd_direction,hyd_localization,hyd_oxygen_tolerance])
    return out_file

def tes1():
    hydrogen_consumer_pipeline_output_folder=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\dummy_output"
    taxon_output=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\dummy_CAT_result.tsv"
    hydrogen_classes_activity_file= r"C:\Users\abel\Documents\hydrogen_consumers\hydrogenase\hydrogenase_classes_activity2.csv"
    out_file=combine_batches(hydrogen_consumer_pipeline_output_folder,
        os.path.join(hydrogen_consumer_pipeline_output_folder,"all_output.csv"))
    make_feature_table(combined_hyd_tsv=out_file,out_file=os.path.join(hydrogen_consumer_pipeline_output_folder,"feature_table_test.tsv"),
                       hydrogen_classes_activity_file=hydrogen_classes_activity_file,taxon_output=taxon_output)
    exit()
def arg_parser():

    args=argparse.ArgumentParser()
    #add the hydrogen consumer pipeline output folder
    args.add_argument("--hydrogen_consumer_pipeline_output_folder",help="the hydrogen consumer pipeline output folder",default=r'C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\mouse_catalogues')
    #add the hydrogen classes activity file
    args.add_argument("--hydrogen_classes_activity_file",help="the hydrogen classes activity file",default=r'C:\Users\abel\Documents\hydrogenases\hydrogenase_classes_activity.csv')
    #add the taxon output file
    args.add_argument("--taxon_output",help="the taxon output file",default=r"C:\Users\abel\Documents\mouse_gut\hydrogen_consumer_pipeline\out.BAT.bin2classification_name.txt")

    args.add_argument("--out_file",help="the output file",default=os.path.join(hydrogen_consumer_pipeline_output_folder,"feature_table_test.tsv"))
    return args.parse_args()
def main(hydrogen_consumer_pipeline_output_folder,hydrogen_classes_activity_file,taxon_output,out_file):
    #combine the files
    out_file1=combine_batches(hydrogen_consumer_pipeline_output_folder,
        os.path.join(hydrogen_consumer_pipeline_output_folder,"all_output.csv"))
    #make the feature table

    make_feature_table(combined_hyd_tsv=out_file1,out_file=out_file,hydrogen_classes_activity_file=hydrogen_classes_activity_file,taxon_output=taxon_output)
    return out_file
if __name__=="__main__":
    #tes1()
    args=arg_parser()
    #combine the files
    out_file = combine_batches(args.hydrogen_consumer_pipeline_output_folder,
        os.path.join(args.hydrogen_consumer_pipeline_output_folder, "all_output.csv"))
    #make the feature table
    make_feature_table(combined_hyd_tsv=out_file,out_file=args.out_file,hydrogen_classes_activity_file=args.hydrogen_classes_activity_file,taxon_output=args.taxon_output)

    #out_file=combine_batches(os.path.join(hydrogen_consumer_pipeline_output_folder,"all_output.csv"))
    #make_feature_table(combined_hyd_tsv=r"C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\results\all_output_new.csv"
    #                   ,out_file=os.path.join(hydrogen_consumer_pipeline_output_folder,"feature_table_test.tsv"))