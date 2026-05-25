import argparse

import pandas as pd
def is_acetogen(df,genome_name):
    #a genome is an acetogen if the protein match column contains FTFHS
    #filter the df for the genome_name
    df=df[df['genome_name']==genome_name]
    #check if the protein match column contains FTFHS
    for protein_match in df['protein_match']:
        if 'FTFHS' in str(protein_match):
            return True

def scan_acetogens(df):
    #unique genome names
    genome_names=df['genome_name'].unique()
    #set the dtype of is_acetogen column to boolean, create it if it does not exist
    if 'is_acetogen' not in df.columns:
        #create the is_acetogen column
        df['is_acetogen']=False
    df['is_acetogen']=df['is_acetogen'].astype('bool')
    #for each genome name
    for genome_name in genome_names:
        #check if it is an acetogen
        acetogen=is_acetogen(df,genome_name)
        #if it is an acetogen
        if acetogen:
            #set the is_acetogen column to True
            print(genome_name,'is an acetogen')
            df.loc[df['genome_name']==genome_name,'is_acetogen']=True
        else:
            #set the is_acetogen column to False
            df.loc[df['genome_name']==genome_name,'is_acetogen']=False

def is_sulphate_reducer(df,genome_name):
    #set the dtype of is_sulphate_reducer column to boolean
    df['is_sulphate_reducer']=df['is_sulphate_reducer'].astype('bool')
    #a genome is a sulphate reducer if the protein match column contains dissimilatory
    synonyms=['dissimilatory','dsr']
    #filter the df for the genome_name
    df=df[df['genome_name']==genome_name]
    #check if the protein match column contains DsrAB
    for protein_match in df['protein_match']:
        if any(synonym in str(protein_match) for synonym in synonyms):
            return True
    return False

def scan_sulphate_reducers(df):
    #set the dtype of is_sulphate_reducer column to boolean
    #create the is_sulphate_reducer column if it does not exist
    if 'is_sulphate_reducer' not in df.columns:
        #create the is_sulphate_reducer column
        df['is_sulphate_reducer']=False
    df['is_sulphate_reducer']=df['is_sulphate_reducer'].astype('bool')
    #unique genome names
    genome_names=df['genome_name'].unique()
    #for each genome name
    for genome_name in genome_names:
        #check if it is a sulphate reducer
        sulphate_reducer=is_sulphate_reducer(df,genome_name)
        #if it is a sulphate reducer
        if sulphate_reducer:
            #set the is_sulphate_reducer column to True
            print(genome_name,'is a sulphate reducer')
            df.loc[df['genome_name']==genome_name,'is_sulphate_reducer']=True
        else:
            #set the is_sulphate_reducer column to False
            df.loc[df['genome_name']==genome_name,'is_sulphate_reducer']=False
def is_methanogen(df,genome_name):
    #a genome is a methanogen if the protein match column contains Methyl-CoM reductase
    synonyms=['methyl coenzyme','methyl-coenzyme','coenzyme-B sulfo']
    #filter the df for the genome_name
    df=df[df['genome_name']==genome_name]
    #check if the protein match column contains Methyl-CoM reductase
    for protein_match in df['protein_match']:
        if any(synonym in str(protein_match) for synonym in synonyms):
            return True
    return False
def scan_methanogens(df):
    #set the dtype of is_methanogen column to boolean
    #create the is_methanogen column if it does not exist
    if 'is_methanogen' not in df.columns:
        #create the is_methanogen column
        df['is_methanogen']=False
    df['is_methanogen']=df['is_methanogen'].astype('bool')
    #unique genome names
    genome_names=df['genome_name'].unique()
    #for each genome name
    for genome_name in genome_names:
        #check if it is a methanogen
        methanogen=is_methanogen(df,genome_name)
        #if it is a methanogen
        if methanogen:
            #set the is_methanogen column to True
            print(genome_name,'is a methanogen')
            df.loc[df['genome_name']==genome_name,'is_methanogen']=True
        else:
            #set the is_methanogen column to False
            df.loc[df['genome_name']==genome_name,'is_methanogen']=False
def contains_uptake_hydrogenase(df,genome_name):
    #a genome contains a uptake hydrogenase if the hyd_direction column contains uptake
    #filter the df for the genome_name
    df=df[df['genome_name']==genome_name]
    #check if the hyd_direction column contains uptake
    for hyd_direction in df['hyd_direction']:
        if 'uptake' in str(hyd_direction):
            return True
    return False

def contains_evolution_hydrogenase(df,genome_name):
    #a genome contains a evolution hydrogenase if the hyd_direction column contains evolution
    #filter the df for the genome_name
    df=df[df['genome_name']==genome_name]
    #check if the hyd_direction column contains evolution
    for hyd_direction in df['hyd_direction']:
        if 'evolution' in str(hyd_direction):
            return True
    return False
def contains_bidirectional_or_bifurcating_hydrogenase(df,genome_name):
    #a genome contains a bidirectional or bifurcating hydrogenase if the hyd_direction column contains bidirectional or bifurcating
    #filter the df for the genome_name
    df=df[df['genome_name']==genome_name]
    #check if the hyd_direction column contains bidirectional or bifurcating
    for hyd_direction in df['hyd_direction']:
        if 'bidirectional' in str(hyd_direction) or 'bifurcating' in str(hyd_direction):
            return True
    return False
def contains_evolution_hydrogenase(df,genome_name):
    #a genome contains a evolution hydrogenase if the hyd_direction column contains evolution
    #filter the df for the genome_name
    df=df[df['genome_name']==genome_name]
    #check if the hyd_direction column contains evolution
    for hyd_direction in df['hyd_direction']:
        if 'evolution' in str(hyd_direction):
            return True
    return False
def scan_hydrogenases(df):
    #unique genome names
    genome_names=df['genome_name'].unique()
    #create the columns
    if 'contains_uptake_hydrogenase' not in df.columns:
        df['contains_uptake_hydrogenase']=False
    if 'contains_bidirectional_or_bifurcating_hydrogenase' not in df.columns:
        df['contains_bidirectional_or_bifurcating_hydrogenase']=False
    if 'contains_evolution_hydrogenase' not in df.columns:
        df['contains_evolution_hydrogenase']=False
    #set the dtype of the columns to boolean
    df['contains_uptake_hydrogenase']=df['contains_uptake_hydrogenase'].astype('bool')
    df['contains_bidirectional_or_bifurcating_hydrogenase']=df['contains_bidirectional_or_bifurcating_hydrogenase'].astype('bool')
    df['contains_evolution_hydrogenase']=df['contains_evolution_hydrogenase'].astype('bool')
    #for each genome name
    for genome_name in genome_names:
        #check if it contains a uptake hydrogenase
        uptake_hydrogenase=contains_uptake_hydrogenase(df,genome_name)
        #check if it contains a bidirectional or bifurcating hydrogenase
        bidirectional_or_bifurcating_hydrogenase=contains_bidirectional_or_bifurcating_hydrogenase(df,genome_name)
        #check if it contains a evolution hydrogenase
        evolution_hydrogenase=contains_evolution_hydrogenase(df,genome_name)
        #if it contains a uptake hydrogenase
        if uptake_hydrogenase:
            #set the consuming_hydrogenase column to uptake
            print(genome_name,'contains a uptake hydrogenase')
            df.loc[df['genome_name']==genome_name,'contains_uptake_hydrogenase']=True
        else:
            #set the consuming_hydrogenase column to False
            df.loc[df['genome_name']==genome_name,'contains_uptake_hydrogenase']=False
        #if it contains a bidirectional or bifurcating hydrogenase
        if bidirectional_or_bifurcating_hydrogenase:
            #set the consuming_hydrogenase column to bidirectional_or_bifurcating
            print(genome_name,'contains a bidirectional or bifurcating hydrogenase')
            df.loc[df['genome_name']==genome_name,'contains_bidirectional_or_bifurcating_hydrogenase']=True
        else:
            #set the consuming_hydrogenase column to False
            df.loc[df['genome_name']==genome_name,'contains_bidirectional_or_bifurcating_hydrogenase']=False
        #if it contains a evolution hydrogenase
        if evolution_hydrogenase:
            #set the consuming_hydrogenase column to evolution
            print(genome_name,'contains a evolution hydrogenase')
            df.loc[df['genome_name']==genome_name,'contains_evolution_hydrogenase']=True
        else:
            #set the consuming_hydrogenase column to False
            df.loc[df['genome_name']==genome_name,'contains_evolution_hydrogenase']=False

def is_hydrogen_consumer(df,genome_name):
    #a genome is a hydrogen consumer if
    #1) it has a uptake hydrogenase
    # or
    #2) it has a bidirectional hydrogenase or bifurcating hydrogenase, and it is a methanogen or acetogen, or it is a sulphate reducer
    #filter the df for the genome_name
    df=df[df['genome_name']==genome_name]
    #check if it has a uptake hydrogenase
    for contains_uptake_hydrogenase in df['contains_uptake_hydrogenase']:
        if contains_uptake_hydrogenase:
            return True
    #check if it has a bidirectional hydrogenase or bifurcating hydrogenase
    for contains_bidirectional_or_bifurcating_hydrogenase in df['contains_bidirectional_or_bifurcating_hydrogenase']:
        if contains_bidirectional_or_bifurcating_hydrogenase:
            #check if it is a methanogen or acetogen, or it is a sulphate reducer
            for is_methanogen in df['is_methanogen']:
                if is_methanogen:
                    return True
            for is_acetogen in df['is_acetogen']:
                if is_acetogen:
                    return True
            for is_sulphate_reducer in df['is_sulphate_reducer']:
                if is_sulphate_reducer:
                    return True
    return False
def scan_hydrogen_consumers(df):
    #set the dtype of is_hydrogen_consumer column to boolean
    #create the is_hydrogen_consumer column if it does not exist
    if 'is_hydrogen_consumer' not in df.columns:
        #create the is_hydrogen_consumer column
        df['is_hydrogen_consumer']=False
    df['is_hydrogen_consumer']=df['is_hydrogen_consumer'].astype('bool')
    #set the dtype of is_hydrogen_producer column to boolean
    #create the is_hydrogen_producer column if it does not exist
    if 'is_hydrogen_producer' not in df.columns:
        #create the is_hydrogen_producer column
        df['is_hydrogen_producer']=False
    #unique genome names
    genome_names=df['genome_name'].unique()
    #for each genome name
    for genome_name in genome_names:
        #check if it is a hydrogen consumer
        hydrogen_consumer=is_hydrogen_consumer(df,genome_name)
        #if it is a hydrogen consumer
        if hydrogen_consumer:
            #set the is_hydrogen_consumer column to True
            print(genome_name,'is a hydrogen consumer')
            df.loc[df['genome_name']==genome_name,'is_hydrogen_consumer']=True
        else:
            #set the is_hydrogen_consumer column to False
            df.loc[df['genome_name']==genome_name,'is_hydrogen_consumer']=False
        hydrogen_producer=contains_evolution_hydrogenase(df,genome_name)
        #if it is a hydrogen producer
        if hydrogen_producer:
            #set the is_hydrogen_producer column to True
            print(genome_name,'is a hydrogen producer')
            df.loc[df['genome_name']==genome_name,'is_hydrogen_producer']=True


def   arguments():
    args=argparse.ArgumentParser()
    args.add_argument("--input_file",help="the input file containing the feature table")
    args.add_argument("--out_file",help="the output file")
    args=args.parse_args()
    return args

def main(input_file,out_file):
        df=pd.read_csv(input_file, sep='\t')
        #drop the is_hydrogen_consumer , is_sulphate_reducer, is_methanogen, is_acetogen, consuming_hydrogenase, bifurcating_or_bidirectional_hydrogenase columns
        df=df.drop(['is_hydrogen_consumer','is_sulphate_reducer','is_methanogen','is_acetogen','consuming_hydrogenase','bifurcating_or_bidirectional_hydrogenase'],axis=1)
        scan_hydrogenases(df)
        scan_acetogens(df)
        scan_sulphate_reducers(df)
        scan_methanogens(df)
        scan_hydrogen_consumers(df)
        #df["is_hydrogen_producer"] = df["contains_evolution_hydrogenase"]
        #save df as csv
        df.to_csv(out_file, sep='\t',index=False)
        return out_file
if __name__ == '__main__':
    args=arguments()
    main(**args)
    exit()
    df=pd.read_csv(r'C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\results\feature_table_new1.csv', sep='\t')
    '''the df looks like this:
    genome_name	is_hydrogen_consumer	is_sulphate_reducer	is_methanogen	is_acetogen	consuming_hydrogenase	bifurcating_or_bidirectional_hydrogenase	kingdom	phylum	class	order	family	genus	species	protein_match	hyd_direction	hyd_localization	hyd_oxygen_tolerance
    single-wildmice_C1249	False	False	False	False	False	True	Bacteria	Firmicutes	Clostridia	Clostridiales	Unknown	Unknown	Unknown	NA	NA	NA	NA
    single-wildmice_C1245	False	False	False	False	False	True	Bacteria	Firmicutes	Clostridia	Clostridiales	Unknown	Unknown	Unknown	NA	NA	NA	NA
    single-wildmice_C1244	False	False	False	False	False	True	Bacteria	Bacteroidetes	Bacteroidia	Bacteroidales	Bacteroidaceae	Bacteroides	Unknown	NA	NA	NA	NA
    single-wildmice_C1240	False	False	False	False	False	True	Bacteria	Bacteroidetes	Bacteroidia	Bacteroidales	Bacteroidaceae	Bacteroides	Unknown	FeFe-WP_024986193.1 - Bacteroides acidifaciens - [FeFe] Group A3	Electron-bifurcation	Cytosolic	O2-labile
    single-wildmice_C1236	False	False	False	False	False	True	Bacteria	Firmicutes	Clostridia	Clostridiales	Unknown	Unknown	Unknown	NA	NA	NA	NA
    single-wildmice_C1235	False	False	False	False	False	True	Bacteria	Firmicutes	Clostridia	Clostridiales	Unknown	Unknown	Unknown	FeFe-WP_007862856.1 - Lachnoclostridium citroniae - [FeFe] Group C1	H2-sensing?	Cytosolic	O2-labile
    single-wildmice_C1235	False	False	False	False	False	True	Bacteria	Firmicutes	Clostridia	Clostridiales	Unknown	Unknown	Unknown	FeFe-WP_006778965.1 - Lachnoclostridium hathewayi - [FeFe] Group A3	Electron-bifurcation	Cytosolic	O2-labile
    single-wildmice_C1235	False	False	False	False	False	True	Bacteria	Firmicutes	Clostridia	Clostridiales	Unknown	Unknown	Unknown	FeFe-WP_002568888.1 - Lachnoclostridium bolteae - [FeFe] Group B	H2-evolution?	Cytosolic	O2-labile
    single-wildmice_C1235	False	False	False	False	False	True	Bacteria	Firmicutes	Clostridia	Clostridiales	Unknown	Unknown	Unknown	FeFe-WP_027640429.1 - Lachnoclostridium clostridioforme - [FeFe] Group C3	H2-sensing?	Cytosolic	O2-labile
    single-wildmice_C1233	False	False	False	False	False	True	Bacteria	Firmicutes	Clostridia	Clostridiales	Unknown	Unknown	Unknown	FeFe-WP_010251944.1 - Acetivibrio cellulolyticus - [FeFe] Group A3	Electron-bifurcation	Cytosolic	O2-labile
    
    
    '''

    #drop the is_hydrogen_consumer , is_sulphate_reducer, is_methanogen, is_acetogen, consuming_hydrogenase, bifurcating_or_bidirectional_hydrogenase columns
    df=df.drop(['is_hydrogen_consumer','is_sulphate_reducer','is_methanogen','is_acetogen','consuming_hydrogenase','bifurcating_or_bidirectional_hydrogenase'],axis=1)
    scan_acetogens(df)
    scan_sulphate_reducers(df)
    scan_methanogens(df)


    scan_hydrogen_consumers(df)
    df["is_hydrogen_producer"] = df["contains_evolution_hydrogenase"]
    #save df as csv
    df.to_csv(r'C:\Users\abel\Documents\hydrogen_consumers\hydrogen_consumer_pipeline_test\results\feature_table_fixed.csv', sep='\t',index=False)
