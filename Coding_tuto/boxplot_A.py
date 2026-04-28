import argparse
import sys
sys.path.append('../..')

parser = argparse.ArgumentParser(description='This script is used to plot the boxplots for the FPM and conservation score given in circpedia for our detected circRNA transcripts; and compare the set of selected vs. not-selected transcripts. Selected ones have more than x read counts and are in circpedia with RNaseR validation.')

parser.add_argument('--path_experiment',dest='path_experiment',default='/mnt/efs/home/npont/run_circexplorer2/output_circexplorer2/reproducibility',help='Provide the path of the experiment results.')

parser.add_argument('--path_db', dest='path_db', default='/mnt/efs/home/npont/references/circpedia/human_hg38_All_circRNA.csv', help='Provide the path of the database circpedia, which is a csv file downloadable at: http://yang-laboratory.com/circpedia/download')

args = parser.parse_args()
path_experiment = args.path_experiment
path_db = args.path_db



import os
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
import numpy as np
import seaborn as sns
from count_matrix_iso_positions import generate_count_matrix
import pandas as pd
from pandas import read_csv
from check_matches_circpedia_RNaseR import check_database_circpedia, check_database_circpedia_and_filter_seq_type, create_df_db_circpedia, intersection
import matplotlib.patches as mpatches

count_matrix=generate_count_matrix(path_experiment,False)

df_bed=create_df_db_circpedia(path_db)

sample1 = 'A1'
sample3 = 'A3'
#sample4 = 'B1'
#sample5 = 'B2'
#sample6 = 'B3'
#sample7 = 'C1'
#sample8 = 'C2'
#sample9 = 'C3'
sample1_ = 'A4'
sample2_ = 'A5'
sample3_ = 'A6'
#sample4_ = 'B4'
#sample5_ = 'B5'
#sample6_ = 'B6'
#sample7_ = 'C4'
#sample8_ = 'C5'
#sample9_ = 'C6'

def generate_list_genes(threshold):
    genes_sample1 = set(
    gene
    for gene in count_matrix[sample1][count_matrix[sample1] >= threshold].index)

    genes_sample3 = set(
    gene
    for gene in count_matrix[sample3][count_matrix[sample3] >= threshold].index)
    
    """
    genes_sample4 = set(
    gene
    for gene in count_matrix[sample4][count_matrix[sample4] >= threshold].index)

    genes_sample5 = set(
    gene
    for gene in count_matrix[sample5][count_matrix[sample5] >= threshold].index)

    genes_sample6 = set(
    gene
    for gene in count_matrix[sample6][count_matrix[sample6] >= threshold].index)

    genes_sample7 = set(
    gene
    for gene in count_matrix[sample7][count_matrix[sample7] >= threshold].index)

    genes_sample8 = set(
    gene
    for gene in count_matrix[sample8][count_matrix[sample8] >= threshold].index)

    genes_sample9 = set(
    gene
    for gene in count_matrix[sample9][count_matrix[sample9] >= threshold].index)
    """
    genes_sample1_ = set(
    gene
    for gene in count_matrix[sample1_][count_matrix[sample1_] >= threshold].index)

    genes_sample2_ = set(
    gene
    for gene in count_matrix[sample2_][count_matrix[sample2_] >= threshold].index)

    genes_sample3_ = set(
    gene
    for gene in count_matrix[sample3_][count_matrix[sample3_] >= threshold].index)
    """
    genes_sample4_ = set(
    gene
    for gene in count_matrix[sample4_][count_matrix[sample4_] >= threshold].index)

    genes_sample5_ = set(
    gene
    for gene in count_matrix[sample5_][count_matrix[sample5_] >= threshold].index)

    genes_sample6_ = set(
    gene
    for gene in count_matrix[sample6_][count_matrix[sample6_] >= threshold].index)

    genes_sample7_ = set(
    gene
    for gene in count_matrix[sample7_][count_matrix[sample7_] >= threshold].index)

    genes_sample8_ = set(
    gene
    for gene in count_matrix[sample8_][count_matrix[sample8_] >= threshold].index)

    genes_sample9_ = set(
    gene
    for gene in count_matrix[sample9_][count_matrix[sample9_] >= threshold].index)
    """

    # Convertion
    A1 = set(list(genes_sample1))
    if not A1:
        A1.add(0)
    A3 = set(list(genes_sample3))
    if not A3:
        A3.add(0)
    """
    B1 = set(list(genes_sample4))
    if not B1:
        B1.add(0)
    B2 = set(list(genes_sample5))
    if not B2:
        B2.add(0)
    B3 = set(list(genes_sample6))
    if not B3:
        B3.add(0)
    C1 = set(list(genes_sample7))
    if not C1:
        C1.add(0)
    C2 = set(list(genes_sample8))
    if not C2:
        C2.add(0)
    C3 = set(list(genes_sample9))
    if not C3:
        C3.add(0)
    """
    A4 = set(list(genes_sample1_))
    if not A4:
        A4.add(0)
    A5 = set(list(genes_sample2_))
    if not A5:
        A5.add(0)
    A6 = set(list(genes_sample3_))
    if not A6:
        A6.add(0)
    """
    B4 = set(list(genes_sample4_))
    if not B4:
        B4.add(0)
    B5 = set(list(genes_sample5_))
    if not B5:
        B5.add(0)
    B6 = set(list(genes_sample6_))
    if not B6:
        B6.add(0)
    C4 = set(list(genes_sample7_))
    if not C4:
        C4.add(0)
    C5 = set(list(genes_sample8_))
    if not C5:
        C5.add(0)
    C6 = set(list(genes_sample9_))
    if not C6:
        C6.add(0)
    """
    
    return A1,A3,A4,A5,A6

# Function to generate custom palette
def custom_palette(column_name):
    if 'not' in column_name:
        return 'blue'  # Color for not selected
    else:
        return 'green'   # Color for selected

thresholds=[10] #[1,2] #we will get the boxplots for these two thresholds candidates (based on plot_right.py and plot_left.py plots)

for threshold in thresholds:
    
    ###########################SELECTED SET##############################
    
    dict_list={'A1':[],'A3':[],'A4':[],'A5':[],'A6':[]}
               #,'B1':[],'B2':[],'B3':[],'B4':[],'B5':[],'B6':[],'C1':[],'C2':[],'C3':[],'C4':[],'C5':[],'C6':[]}
    
    A1_,A3_,A4_,A5_,A6_=generate_list_genes(threshold) #,B1_,B2_,B3_,B4_,B5_,B6_,C1_,C2_,C3_,C4_,C5_,C6_
    
    dict_list={'A1': A1_,'A3': A3_,'A4': A4_,'A5': A5_,'A6': A6_} #,'B1': B1_,'B2': B2_,'B3': B3_,'B4': B4_,'B5': B5_,'B6': B6_,'C1': C1_,'C2': C2_,'C3': C3_,'C4': C4_,'C5': C5_,'C6': C6_
    
    dict_list_intact={'A1': A1_,'A3': A3_,'A4': A4_,'A5': A5_,'A6': A6_} #,'B1': B1_,'B2': B2_,'B3': B3_,'B4': B4_,'B5': B5_,'B6': B6_,'C1': C1_,'C2': C2_,'C3': C3_,'C4': C4_,'C5': C5_,'C6': C6_
                                                                        
    dict_match_A={'A1': pd.DataFrame(), 'A3': pd.DataFrame(), 'A4': pd.DataFrame(), 'A5': pd.DataFrame(), 'A6': pd.DataFrame()}
    #dict_match_B={'B1': pd.DataFrame(), 'B2': pd.DataFrame(),'B3': pd.DataFrame(), 'B4': pd.DataFrame(), 'B5': pd.DataFrame(), 'B6': pd.DataFrame()}
    #dict_match_C={'C1': pd.DataFrame(), 'C2': pd.DataFrame(),'C3': pd.DataFrame(), 'C4': pd.DataFrame(), 'C5': pd.DataFrame(), 'C6': pd.DataFrame()}
    
    #store the whole database info of the detected transcripts found in circpedia with RNaseR, for each sample, in dict_match_A/B/C[sample_name]
    for sample,transcripts in dict_list.items():       
        if sample.startswith('A'):
            length,intersect,dict_match_A[f'{sample}']=check_database_circpedia_and_filter_seq_type(df_bed,transcripts,'null_id.txt','null_full.txt',True,False,True)
        """
        if sample.startswith('B'):
            length,intersect,dict_match_B[f'{sample}']=check_database_circpedia_and_filter_seq_type(df_bed,transcripts,'null_id.txt','null_full.txt',True,False,True)
        if sample.startswith('C'):
            length,intersect,dict_match_C[f'{sample}']=check_database_circpedia_and_filter_seq_type(df_bed,transcripts,'null_id.txt','null_full.txt',True,False,True)
        """    
        #replace the transcripts associated to a key i.e. a sample, with the list of those in circpedia&RNaseR+
        dict_list[f'{sample}']=intersect
        

    #filter the dataframes containing the info about the detected transcripts, to keep only conservation and FPM scores well formatted
    for sample, df_A in dict_match_A.items():
        dict_match_A[sample] = df_A[~df_A.index.duplicated(keep='first')]
        dict_match_A[sample] = dict_match_A[sample][['FPM', 'conservation']]
        dict_match_A[sample]['conservation'] = dict_match_A[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_match_A[sample]['conservation'] = dict_match_A[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_match_A[sample]['conservation'] = dict_match_A[sample]['conservation'].astype('int64')
        dict_match_A[sample].rename(columns={"FPM" : f'FPM_{sample}', "conservation" : f'conservation_{sample}'}, inplace=True)
    """
    # Update dict_match_B
    for sample, df_B in dict_match_B.items():
        dict_match_B[sample] = df_B[~df_B.index.duplicated(keep='first')]
        dict_match_B[sample] = dict_match_B[sample][['FPM', 'conservation']]
        dict_match_B[sample]['conservation'] = dict_match_B[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_match_B[sample]['conservation'] = dict_match_B[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_match_B[sample]['conservation'] = dict_match_B[sample]['conservation'].astype('int64')
        dict_match_B[sample].rename(columns={"FPM" : f'FPM_{sample}', "conservation" : f'conservation_{sample}'}, inplace=True)

    # Update dict_match_C
    for sample, df_C in dict_match_C.items():
        dict_match_C[sample] = df_C[~df_C.index.duplicated(keep='first')]
        dict_match_C[sample] = dict_match_C[sample][['FPM', 'conservation']]
        dict_match_C[sample]['conservation'] = dict_match_C[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_match_C[sample]['conservation'] = dict_match_C[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_match_C[sample]['conservation'] = dict_match_C[sample]['conservation'].astype('int64')
        dict_match_C[sample].rename(columns={"FPM" : f'FPM_{sample}', "conservation" : f'conservation_{sample}'}, inplace=True)
    """
    #Group FPM data to prepare the input for the boxplot    
    df_A_FPM=pd.concat([df[f'FPM_{sample}'] for sample,df in dict_match_A.items()],axis=1)
    print(f'df_A_FPM: {df_A_FPM.head(10)}')
    #df_B_FPM=pd.concat([df[f'FPM_{sample}'] for sample,df in dict_match_B.items()],axis=1)
    #df_C_FPM=pd.concat([df[f'FPM_{sample}'] for sample,df in dict_match_C.items()],axis=1)
    
    #Group conservation data to prepare for the boxplot
    df_A_conservation=pd.concat([df[f'conservation_{sample}'] for sample,df in dict_match_A.items()],axis=1)
    print(f'df_A_conservation: {df_A_conservation.head(10)}')    
    #df_B_conservation=pd.concat([df[f'conservation_{sample}'] for sample,df in dict_match_B.items()],axis=1)
    #df_C_conservation=pd.concat([df[f'conservation_{sample}'] for sample,df in dict_match_C.items()],axis=1)
    
    ############################NON-SELECTED SET############################
    
    dict_match_A_={'A1': pd.DataFrame(), 'A3': pd.DataFrame(), 'A4': pd.DataFrame(), 'A5': pd.DataFrame(), 'A6': pd.DataFrame()}
    #dict_match_B_={'B1': pd.DataFrame(), 'B2': pd.DataFrame(),'B3': pd.DataFrame(), 'B4': pd.DataFrame(), 'B5': pd.DataFrame(), 'B6': pd.DataFrame()}
    #dict_match_C_={'C1': pd.DataFrame(), 'C2': pd.DataFrame(),'C3': pd.DataFrame(), 'C4': pd.DataFrame(), 'C5': pd.DataFrame(), 'C6': pd.DataFrame()}
    
    dict_inferior_threshold={'A1':[],'A3':[],'A4':[],'A5':[],'A6':[]}
    #,'B1':[],'B2':[],'B3':[],'B4':[],'B5':[],'B6':[],'C1':[],'C2':[],'C3':[],'C4':[],'C5':[],'C6':[]}
    
    list_samples=['A1','A3','A4','A5','A6']
    #,'B1','B2','B3','B4','B5','B6','C1','C2','C3','C4','C5','C6']
    #Get the list of transcripts that did not pass the threshold
    for sample in list_samples:
        dict_inferior_threshold[f'{sample}']=count_matrix[~count_matrix[f'{sample}'].index.isin(dict_list_intact[f'{sample}'])].index
        print(f'List of transcripts that did not pass the threshold for sample {sample}: ', len(list(dict_inferior_threshold[f'{sample}'])))
    
    #Get the info of these transcripts from the db
    for sample,transcripts in dict_inferior_threshold.items():
        if sample.startswith('A'):
            print(f'check db on transcripts that did not pass the threshold for sample {sample}, store in dict_match_sample_')
            length,intersect,dict_match_A_[f'{sample}']=check_database_circpedia(df_bed,transcripts,'null.txt','null.txt',True,True)
        #if sample.startswith('B'):
            #print(f'check db on transcripts that did not pass the threshold for sample {sample}, store in dict_match_sample_')
            #length,intersect,dict_match_B_[f'{sample}']=check_database_circpedia(df_bed,transcripts,'null.txt','null.txt',True,True)
        #if sample.startswith('C'):
            #print(f'check db on transcripts that did not pass the threshold for sample {sample}, store in dict_match_sample_')
            #length,intersect,dict_match_C_[f'{sample}']=check_database_circpedia(df_bed,transcripts,'null.txt','null.txt',True,True)
            
    for sample, df_A_ in dict_match_A_.items():
        print(f'format correctly dict_match_sample_ columns for sample {sample}')
        dict_match_A_[sample] = df_A_[~df_A_.index.duplicated(keep='first')]
        dict_match_A_[sample] = dict_match_A_[sample][['FPM', 'conservation']]
        dict_match_A_[sample]['conservation'] = dict_match_A_[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_match_A_[sample]['conservation'] = dict_match_A_[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_match_A_[sample]['conservation'] = dict_match_A_[sample]['conservation'].astype('int64')
        #dict_match_A_[sample].rename(columns={"FPM" : f'FPM_{sample}', "conservation" : f'conservation_{sample}'}, inplace=True)
    
    """
    # Update dict_match_B_
    for sample, df_B_ in dict_match_B_.items():
        print(f'format correctly dict_match_sample_ columns for sample {sample}')
        dict_match_B_[sample] = df_B_[~df_B_.index.duplicated(keep='first')]
        dict_match_B_[sample] = dict_match_B_[sample][['FPM', 'conservation']]
        dict_match_B_[sample]['conservation'] = dict_match_B_[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_match_B_[sample]['conservation'] = dict_match_B_[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_match_B_[sample]['conservation'] = dict_match_B_[sample]['conservation'].astype('int64')
        #dict_match_B_[sample].rename(columns={"FPM" : f'FPM_{sample}', "conservation" : f'conservation_{sample}'}, inplace=True)

    # Update dict_match_C_
    for sample, df_C_ in dict_match_C_.items():
        print(f'format correctly dict_match_sample_ columns for sample {sample}')
        dict_match_C_[sample] = df_C_[~df_C_.index.duplicated(keep='first')]
        dict_match_C_[sample] = dict_match_C_[sample][['FPM', 'conservation']]
        dict_match_C_[sample]['conservation'] = dict_match_C_[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_match_C_[sample]['conservation'] = dict_match_C_[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_match_C_[sample]['conservation'] = dict_match_C_[sample]['conservation'].astype('int64')
        #dict_match_C_[sample].rename(columns={"FPM" : f'FPM_{sample}', "conservation" : f'conservation_{sample}'}, inplace=True)
    """
    print('dict_match_A_[A1] ie in circpedia but not with threshold, and filtered on fpm and conservation',dict_match_A_['A1'].head(10))    
    
    #Get the info of the transcripts that pass the thresholds but do not have a RNaseR validation
    dict_notRNaseR_A={'A1': pd.DataFrame(), 'A3': pd.DataFrame(), 'A4': pd.DataFrame(), 'A5': pd.DataFrame(), 'A6': pd.DataFrame()}
    """
    dict_notRNaseR_B={'B1': pd.DataFrame(), 'B2': pd.DataFrame(),'B3': pd.DataFrame(), 'B4': pd.DataFrame(), 'B5': pd.DataFrame(), 'B6': pd.DataFrame()}
    dict_notRNaseR_C={'C1': pd.DataFrame(), 'C2': pd.DataFrame(),'C3': pd.DataFrame(), 'C4': pd.DataFrame(), 'C5': pd.DataFrame(), 'C6': pd.DataFrame()}
    """
    
    #reminder: A1_,A3_,A4_ etc are the list of transcripts passing the threshold, without any check on circpedia done yet
    #first, just check the db on those transcripts (filter on those who do not have RNaseR treatment will be done in next steps)
    print(f'****get the info from db for transcripts passing the threshold, to further check RNaseR presence***')
    print(f'for sample A1')
    length,intersect,dict_notRNaseR_A['A1']=check_database_circpedia(df_bed,A1_,'null.txt','null.txt',True,True)
    print(f'for sample A3')
    length,intersect,dict_notRNaseR_A['A3']=check_database_circpedia(df_bed,A3_,'null.txt','null.txt',True,True)
    print(f'for sample A4')
    length,intersect,dict_notRNaseR_A['A4']=check_database_circpedia(df_bed,A4_,'null.txt','null.txt',True,True)
    print(f'for sample A5')
    length,intersect,dict_notRNaseR_A['A5']=check_database_circpedia(df_bed,A5_,'null.txt','null.txt',True,True)
    print(f'for sample A6')
    length,intersect,dict_notRNaseR_A['A6']=check_database_circpedia(df_bed,A6_,'null.txt','null.txt',True,True)
    
    """
    print('for B1')
    length,intersect,dict_notRNaseR_B['B1']=check_database_circpedia(df_bed,B1_,'null.txt','null.txt',True,True)
    print('for B2')
    length,intersect,dict_notRNaseR_B['B2']=check_database_circpedia(df_bed,B2_,'null.txt','null.txt',True,True)
    print('for B3')
    length,intersect,dict_notRNaseR_B['B3']=check_database_circpedia(df_bed,B3_,'null.txt','null.txt',True,True)
    print('for B4')
    length,intersect,dict_notRNaseR_B['B4']=check_database_circpedia(df_bed,B4_,'null.txt','null.txt',True,True)
    print('for B5')
    length,intersect,dict_notRNaseR_B['B5']=check_database_circpedia(df_bed,B5_,'null.txt','null.txt',True,True)
    print('for B6')
    length,intersect,dict_notRNaseR_B['B6']=check_database_circpedia(df_bed,B6_,'null.txt','null.txt',True,True)
    
    print('for C1')
    length,intersect,dict_notRNaseR_C['C1']=check_database_circpedia(df_bed,C1_,'null.txt','null.txt',True,True)
    print('for C2')
    length,intersect,dict_notRNaseR_C['C2']=check_database_circpedia(df_bed,C2_,'null.txt','null.txt',True,True)
    print('for C3')
    length,intersect,dict_notRNaseR_C['C3']=check_database_circpedia(df_bed,C3_,'null.txt','null.txt',True,True)
    print('for C4')
    length,intersect,dict_notRNaseR_C['C4']=check_database_circpedia(df_bed,C4_,'null.txt','null.txt',True,True)
    print('for C5')
    length,intersect,dict_notRNaseR_C['C5']=check_database_circpedia(df_bed,C5_,'null.txt','null.txt',True,True)
    print('for C6')
    length,intersect,dict_notRNaseR_C['C6']=check_database_circpedia(df_bed,C6_,'null.txt','null.txt',True,True)
    """
    
    #let's filter the df of each dict_notRNaseR, to keep only rows where no RNaseR is contained in seq_type column
    for sample,df in dict_notRNaseR_A.items():
        print(f'*****filter out transcripts not having a RNaseR treatment at all for sample {sample} *****')
        # Identify duplicated indices where at least one row contains 'RNaseR' in seq_type
        duplicated_indices = df[df['seq_type'] == 'RNaseR'].index.unique()
        print(f'index where one has RNaseR for sample {sample}')
        # Keep only rows where the index is not in duplicated_indices
        dict_notRNaseR_A[sample] = df[~df.index.isin(duplicated_indices)]
        dict_notRNaseR_A[sample] = df[~df.index.duplicated(keep='first')]
        dict_notRNaseR_A[sample] = dict_notRNaseR_A[sample][['FPM', 'conservation']]
        dict_notRNaseR_A[sample]['conservation'] = dict_notRNaseR_A[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_notRNaseR_A[sample]['conservation'] = dict_notRNaseR_A[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_notRNaseR_A[sample]['conservation'] = dict_notRNaseR_A[sample]['conservation'].astype('int64')
    
    """
    for sample,df in dict_notRNaseR_B.items():
        print(f'*****filter out transcripts not having a RNaseR treatment at all for sample {sample} *****')
        # Identify duplicated indices where at least one row contains 'RNaseR' in seq_type
        duplicated_indices = df[df['seq_type'] == 'RNaseR'].index.unique()
        # Keep only rows where the index is not in duplicated_indices
        dict_notRNaseR_B[sample] = df[~df.index.isin(duplicated_indices)]
        dict_notRNaseR_B[sample] = df[~df.index.duplicated(keep='first')]
        dict_notRNaseR_B[sample] = dict_notRNaseR_B[sample][['FPM', 'conservation']]
        dict_notRNaseR_B[sample]['conservation'] = dict_notRNaseR_B[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_notRNaseR_B[sample]['conservation'] = dict_notRNaseR_B[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_notRNaseR_B[sample]['conservation'] = dict_notRNaseR_B[sample]['conservation'].astype('int64')
        
    for sample,df in dict_notRNaseR_C.items():
        print(f'*****filter out transcripts not having a RNaseR treatment at all for sample {sample} *****')
        # Identify duplicated indices where at least one row contains 'RNaseR' in seq_type
        duplicated_indices = df[df['seq_type'] == 'RNaseR'].index.unique()
        # Keep only rows where the index is not in duplicated_indices
        dict_notRNaseR_C[sample] = df[~df.index.isin(duplicated_indices)]
        dict_notRNaseR_C[sample] = df[~df.index.duplicated(keep='first')]
        dict_notRNaseR_C[sample] = dict_notRNaseR_C[sample][['FPM', 'conservation']]
        dict_notRNaseR_C[sample]['conservation'] = dict_notRNaseR_C[sample]['conservation'].str.replace('MMU_CIRCpedia_', '', regex=False)
        dict_notRNaseR_C[sample]['conservation'] = dict_notRNaseR_C[sample]['conservation'].str.replace('species_specific', '0', regex=False)
        dict_notRNaseR_C[sample]['conservation'] = dict_notRNaseR_C[sample]['conservation'].astype('int64')
    """    
    print('dict_notRNaseR[A1] ie good threshold but not with RNaseR:', dict_notRNaseR_A['A1'].head(10))   
    
    ##concatenate over axis0 (rows) and then rename columns
    df_A1=pd.concat([dict_notRNaseR_A['A1'], dict_match_A_['A1']],axis=0)
    df_A1.rename(columns={"FPM" : 'FPM_A1', "conservation" : 'conservation_A1'}, inplace=True)
    df_A1 = df_A1.reset_index(drop=True)
    df_A3=pd.concat([dict_notRNaseR_A['A3'], dict_match_A_['A3']],axis=0)
    df_A3.rename(columns={"FPM" : 'FPM_A3', "conservation" : 'conservation_A3'}, inplace=True)
    df_A3 = df_A3.reset_index(drop=True)
    df_A4=pd.concat([dict_notRNaseR_A['A4'], dict_match_A_['A4']],axis=0)
    df_A4.rename(columns={"FPM" : 'FPM_A4', "conservation" : 'conservation_A4'}, inplace=True)
    df_A4 = df_A4.reset_index(drop=True)
    df_A5=pd.concat([dict_notRNaseR_A['A5'], dict_match_A_['A5']],axis=0)
    df_A5.rename(columns={"FPM" : 'FPM_A5', "conservation" : 'conservation_A5'}, inplace=True)
    df_A5 = df_A5.reset_index(drop=True)
    df_A6=pd.concat([dict_notRNaseR_A['A6'], dict_match_A_['A6']],axis=0)
    df_A6.rename(columns={"FPM" : 'FPM_A6', "conservation" : 'conservation_A6'}, inplace=True)
    df_A6 = df_A6.reset_index(drop=True)
    
    print(f'has created df_A1 till df_A6')
    
    print('df_A1 ie df from dict_notRNaseR and dict_match_A_, where they have been concat on rows, columns renamed and idx reset', df_A1.head(10))
    """
    df_B1=pd.concat([dict_notRNaseR_B['B1'], dict_match_B_['B1']],axis=0)
    df_B1.rename(columns={"FPM" : 'FPM_B1', "conservation" : 'conservation_B1'}, inplace=True)
    df_B2=pd.concat([dict_notRNaseR_B['B2'], dict_match_B_['B2']],axis=0)
    df_B2.rename(columns={"FPM" : 'FPM_B2', "conservation" : 'conservation_B2'}, inplace=True)
    df_B3=pd.concat([dict_notRNaseR_B['B3'], dict_match_B_['B3']],axis=0)
    df_B3.rename(columns={"FPM" : 'FPM_B3', "conservation" : 'conservation_B1'}, inplace=True)
    df_B4=pd.concat([dict_notRNaseR_B['B4'], dict_match_B_['B4']],axis=0)
    df_B4.rename(columns={"FPM" : 'FPM_B4', "conservation" : 'conservation_B4'}, inplace=True)
    df_B5=pd.concat([dict_notRNaseR_B['B5'], dict_match_B_['B5']],axis=0)
    df_B5.rename(columns={"FPM" : 'FPM_B5', "conservation" : 'conservation_B5'}, inplace=True)
    df_B6=pd.concat([dict_notRNaseR_B['B6'], dict_match_B_['B6']],axis=0)
    df_B6.rename(columns={"FPM" : 'FPM_B6', "conservation" : 'conservation_B6'}, inplace=True)
    df_B1 = df_B1.reset_index(drop=True)
    df_B2 = df_B2.reset_index(drop=True)
    df_B3 = df_B3.reset_index(drop=True)
    df_B4 = df_B4.reset_index(drop=True)
    df_B5 = df_B5.reset_index(drop=True)
    df_B6 = df_B6.reset_index(drop=True)

    print(f'has created df_B1 till df_B6')
    
    df_C1 = pd.concat([dict_notRNaseR_C['C1'], dict_match_C_['C1']],axis=0)
    df_C1.rename(columns={"FPM": 'FPM_C1', "conservation": 'conservation_C1'}, inplace=True)
    df_C2 = pd.concat([dict_notRNaseR_C['C2'], dict_match_C_['C2']],axis=0)
    df_C2.rename(columns={"FPM": 'FPM_C2', "conservation": 'conservation_C2'}, inplace=True)
    df_C3 = pd.concat([dict_notRNaseR_C['C3'], dict_match_C_['C3']],axis=0)
    df_C3.rename(columns={"FPM": 'FPM_C3', "conservation": 'conservation_C3'}, inplace=True)
    df_C4 = pd.concat([dict_notRNaseR_C['C4'], dict_match_C_['C4']],axis=0)
    df_C4.rename(columns={"FPM": 'FPM_C4', "conservation": 'conservation_C4'}, inplace=True)
    df_C5 = pd.concat([dict_notRNaseR_C['C5'], dict_match_C_['C5']],axis=0)
    df_C5.rename(columns={"FPM": 'FPM_C5', "conservation": 'conservation_C5'}, inplace=True)
    df_C6 = pd.concat([dict_notRNaseR_C['C6'], dict_match_C_['C6']],axis=0)
    df_C6.rename(columns={"FPM": 'FPM_C6', "conservation": 'conservation_C6'}, inplace=True)
    df_C1 = df_C1.reset_index(drop=True)
    df_C2 = df_C2.reset_index(drop=True)
    df_C3 = df_C3.reset_index(drop=True)
    df_C4 = df_C4.reset_index(drop=True)
    df_C5 = df_C5.reset_index(drop=True)
    df_C6 = df_C6.reset_index(drop=True)
    
    print(f'has created df_C1 till df_C6')
    """
    #Dataframes for the boxplot, ie grouping the info of replicates for FPM & for conservation separately
    df_not_selected_A_FPM=pd.concat([df_A1['FPM_A1'],df_A3['FPM_A3'],df_A4['FPM_A4'],df_A5['FPM_A5'],df_A6['FPM_A6']],axis=1)
    print(f'df_not_selected_A_FPM', df_not_selected_A_FPM.head(3))
    df_not_selected_A_FPM.rename(columns={"FPM_A1" : 'A1_not', "FPM_A3" : 'A3_not', "FPM_A4" : 'A4_not', "FPM_A5" : 'A5_not', "FPM_A6" : 'A6_not'}, inplace=True)
    #df_not_selected_B_FPM=pd.concat([df_B1['FPM_B1'],df_B2['FPM_B2'],df_B3['FPM_B3'],df_B4['FPM_B4'],df_B5['FPM_B5'],df_B6['FPM_B6']],axis=1)
    #df_not_selected_C_FPM=pd.concat([df_C1['FPM_C1'],df_C2['FPM_C2'],df_C3['FPM_C3'],df_C4['FPM_C4'],df_C5['FPM_C5'],df_C6['FPM_C6']],axis=1)
    
    
    df_not_selected_A_conservation = pd.concat([df_A1['conservation_A1'], df_A3['conservation_A3'], df_A4['conservation_A4'], df_A5['conservation_A5'], df_A6['conservation_A6']],axis=1)
    df_not_selected_A_conservation.rename(columns={"conservation_A1" : 'A1_not', "conservation_A3" : 'A3_not', "conservation_A4" : 'A4_not', "conservation_A5" : 'A5_not', "conservation_A6" : 'A6_not'}, inplace=True)
    
    #df_not_selected_B_conservation = pd.concat([df_B1['conservation_B1'], df_B2['conservation_B2'], df_B3['conservation_B3'], df_B4['conservation_B4'], df_B5['conservation_B5'], df_B6['conservation_B6']],axis=1)
    
    #df_not_selected_C_conservation = pd.concat([df_C1['conservation_C1'], df_C2['conservation_C2'], df_C3['conservation_C3'], df_C4['conservation_C4'], df_C5['conservation_C5'], df_C6['conservation_C6']],axis=1)
    
    ######try to group all boxplots######
    df_A_FPM = df_A_FPM.reset_index(drop=True)
    df_A_FPM.rename(columns={"FPM_A1" : 'A1', "FPM_A3" : 'A3', "FPM_A4" : 'A4', "FPM_A5" : 'A5', "FPM_A6" : 'A6'}, inplace=True)
    FPM_combined=pd.concat([df_not_selected_A_FPM,df_A_FPM],axis=1)
    
    df_A_conservation=df_A_conservation.reset_index(drop=True)
    df_A_conservation.rename(columns={"conservation_A1" : 'A1', "conservation_A3" : 'A3', "conservation_A4" : 'A4', "conservation_A5" : 'A5',"conservation_A6" : 'A6'},inplace=True)
    conservation_combined = pd.concat([df_A_conservation,df_not_selected_A_conservation],axis=1)
    
    #pairs = [("A1", "A3"),("A1","A4"),("A1","A5"),("A1","A6"),("A3","A4"),("A3","A5"),("A3","A6"),("A4","A5"),("A4","A6"),("A5","A6"),("A6_not","A6"),("A5_not","A5"),("A4_not","A4"),("A3_not","A3"),("A1_not","A1"),("A6_not","A1"),("A1_not", "A3_not"),("A1_not","A4_not"),("A1_not","A5_not"),("A1_not","A6_not"),("A3_not","A4_not"),("A3_not","A5_not"),("A3_not","A6_not"),("A4_not","A5_not"),("A4_not","A6_not"),("A5_not","A6_not")]
    pairs = [("A6_not","A5_not"),("A5_not","A4_not"),("A4_not","A3_not"),("A3_not","A1_not"),("A6_not","A6"),("A5_not","A5"),("A4_not","A4"),("A3_not","A3"),("A1_not","A1"),("A6_not","A1"),("A5_not","A3"),("A6","A5"),("A5","A4"),("A4","A3"),("A3","A1")]
    # Create custom palette mapping
    palette_map1 = {column: custom_palette(column) for column in FPM_combined.columns}

    # Plotting boxplot with custom palette
    fig1, ax1 = plt.subplots(figsize=(15, 10))
    box1 = sns.boxplot(data=FPM_combined, palette=palette_map1, ax=ax1, showfliers=False)
    ##statistical test##
    annot = Annotator(ax1, pairs=pairs, data=FPM_combined)
    annot.configure(comparisons_correction="Bonferroni",test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
    annot.apply_test()
    ax1, test_results = annot.annotate()
    ##legend##
    blue_patch = mpatches.Patch(color='blue', label='not selected set')
    green_patch = mpatches.Patch(color='green', label='selected set')
    plt.legend(handles=[blue_patch, green_patch])
    ax1.set_xlabel('Group of transcripts in the corresponding sample ("not" refers to transcripts not selected by the filtering criteria)')
    ax1.set_ylabel('FPM score')
    box1.set_title('FPM values for transcripts in the selected versus not-selected set. \n Plasma A. Outliers were discarded. \n The selected transcripts have more than 10 BSJ counts and are present \n in CIRCpedia with RNase R validation.')
    fig1.savefig(f"../output/A_group_boxplot_FPM_{threshold}BSJ.png")
    plt.close()
    
    # Create custom palette mapping
    palette_map2 = {column: custom_palette(column) for column in conservation_combined.columns}

    # Plotting boxplot with custom palette
    fig2, ax2 = plt.subplots(figsize=(15, 10))
    box2 = sns.boxplot(data=conservation_combined, palette=palette_map2, ax=ax2, showfliers=False)
    ##statistical test##
    annot = Annotator(ax2, pairs=pairs, data=conservation_combined)
    annot.configure(comparisons_correction="Bonferroni",test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
    annot.apply_test()
    ax2, test_results = annot.annotate()
    ##legend##
    blue_patch = mpatches.Patch(color='blue', label='not-selected set')
    green_patch = mpatches.Patch(color='green', label='selected set')
    plt.legend(handles=[blue_patch, green_patch])
    ax2.set_xlabel('Group of transcripts in the corresponding sample ("not" refers to transcripts not selected by the filtering criteria)')
    ax2.set_ylabel('Conservation score')
    box2.set_title('Conservation values for transcripts in the selected versus not-selected set. \n Plasma A. Outliers were discarded. \n The selected transcripts have more than 10 BSJ counts and are present \n in CIRCpedia with RNase R validation.')
    fig2.savefig(f"../output/A_group_boxplot_conservation_{threshold}BSJ.png")
    plt.close()
    
    #########################PLOT#####################################
    """
    ##PLASMA A##
    sns.set_style('white')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 10))
    
    box1 = sns.boxplot(data=df_A_FPM, palette='flare', ax=ax1,showfliers=False)
    box1.set_title('FPM values')
    box2 = sns.boxplot(data=df_A_conservation, palette='flare', ax=ax2, showfliers=False)
    box2.set_title('Conservation scores')
    fig.suptitle('FPM and conservation values for transcripts in the selected set, plasma A.')

    fig.savefig(f"../output/A_boxplot_select_bsj{threshold}.png")
    
    #############NON SELECTED#############
    ##PLASMA A##
    sns.set_style('white')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 10))
    
    box1 = sns.boxplot(data=df_not_selected_A_FPM, palette='flare', ax=ax1,showfliers=False)
    box1.set_title('FPM values')
    box2 = sns.boxplot(data=df_not_selected_A_conservation, palette='flare', ax=ax2, showfliers=False)
    box2.set_title('Conservation scores')
    fig.suptitle('FPM and conservation values for transcripts in the non-selected set, plasma A.')
    fig.savefig(f"../output/A_boxplot_not_select_bsj{threshold}.png")
    """
    