import argparse
import sys
sys.path.append('../..')
############################################################################################################
    # Script used to compute the correlation matrix based on the count matrix for circRNA #
        # Plot the correlation heatmap with hierarchical clustering #
        ## Only considering genes that are detected in all samples (A&B&C)  and in circpedia&RNaseR+##
############################################################################################################

parser = argparse.ArgumentParser(description='Plot the correlation heatmap, with pearson correlation coefficients computed on the count matrix. Clustering of the samples is also done according to the method entered at the command line. WE ONLY CONSIDER TRANSCRIPTS DETECTED ACROSS ALL SAMPLES COMMONLY AND IN CIRCPEDIA&RNaseR+.')

parser.add_argument('--path', dest='path', action='store', help='Parent directory of the circexplorer output of the experiment analyzed.')

parser.add_argument('--clustering', dest='clustering_method', default='single',
                    help='The clustering method can be chosen by the user. Possible methods: complete, single, average. Default: single.')

parser.add_argument('--level', dest='level', default = 'position', help='The user can choose to study at the gene level or at the isoform level. In the latter an isoform can be identified by its ID or chr:start-end. Possible choices: gene, id, position. Default: position.')

parser.add_argument('--path_db', dest='path_db', default='/mnt/efs/home/npont/references/circpedia/human_hg38_All_circRNA.csv', help='Provide the path of the database circpedia, which is a csv file downloadable at: http://yang-laboratory.com/circpedia/download')

args = parser.parse_args()
clustering_method=args.clustering_method
level=args.level
path_db=args.path_db

if args.path is not None:
    path_experiment = args.path
else:
    path_experiment='/mnt/efs/home/npont/run_circexplorer2/output_circexplorer2/reproducibility'

    
import seaborn as sb
#from count_matrix_genes import generate_count_matrix as generate_count_matrix_genes
from count_matrix_iso_positions import generate_count_matrix as generate_count_matrix_iso_position
from count_matrix_isoform import generate_count_matrix as generate_count_matrix_iso_id
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import numpy
from check_matches_circpedia_RNaseR import check_database_circpedia, check_database_circpedia_and_filter_seq_type, create_df_db_circpedia, intersection

"""
#count_matrix is generated out of all circularRNA_known files in count_matrix.py script
if level == 'gene':
    count_matrix=generate_count_matrix_genes(path_experiment,False)
elif level == 'id':
    count_matrix=generate_count_matrix_iso_id(path_experiment,False)
elif level == 'position':
    count_matrix=generate_count_matrix_iso_position(path_experiment,False)
"""

count_matrix=generate_count_matrix_iso_position(path_experiment,False)

# Exclude samples of poor quality
count_matrix=count_matrix.drop(columns=['A2'])

# Function to calculate Pearson correlation coefficient for non-zero pairs
def pearson_corr_nonzero(x, y):
    mask = (x != 0) | (y != 0)
    return x[mask].corr(y[mask])

print('before pre-filtering',count_matrix.head(20))
print(f'Total number of transcripts/genes: {len(count_matrix)}')

df_bed=create_df_db_circpedia(path_db)

#In oder to have only the genes/transcripts in the full intersection:
#need to pre-filter the count matrix to keep only rows where all cell values are different from zero
#which do have the min number of bsj counts ?
thresholds = [10]
for threshold in thresholds:
    list_idx=[]
    for idx, row in count_matrix.iterrows():
        #index of non-zero values
        result = row.to_numpy().nonzero()
        #go and search the value corresponding to that index
        values = row.iloc[result]
        #if all values in the row are different from zero we enter the following condition (if the gene of that row is detected in all samples)
        if len(values)==len(row):
            list_idx.append(idx)

    count_matrix=count_matrix.loc[list_idx]
    count_matrix=count_matrix[count_matrix>=threshold].dropna()
    print('count_matrix after filtering',count_matrix)
    print(f'size of count matrix before checking db: {len(count_matrix)}')
    print(f'Number of genes/transcripts commonly detected across all samples: {len(count_matrix)}')
    
    list_transcripts=set(count_matrix.index.values.tolist())
    
    #which are in the circpedia&RNaseR+ ?
    length,intersect=check_database_circpedia_and_filter_seq_type(df_bed,list_transcripts,'null.txt','null.txt',False,True,False)
    #keep only rows where the index is a transcript contained in intersect i.e. contained in circpedia&RNaseR+
    count_matrix=count_matrix.loc[list(intersect)]
    print('count matrix after checking db:', count_matrix.head(10))
    print(f'size count matrix after checking db: {len(count_matrix)}')
    
    # Get all combinations of column pairs
    column_combinations = list(itertools.combinations(count_matrix.columns, 2))

    # Create a dictionary to store correlation coefficients
    correlation_dict = {}

    # Iterate over column pairs and compute correlation coefficient
    for col1, col2 in column_combinations:
        corr_coefficient = pearson_corr_nonzero(count_matrix[col1], count_matrix[col2])
        correlation_dict[(col1, col2)] = corr_coefficient

    # Create a DataFrame from the correlation dictionary
    correlation_matrix = pd.DataFrame(index=count_matrix.columns, columns=count_matrix.columns)

    # Fill in the correlation matrix
    for (col1, col2), corr_coefficient in correlation_dict.items():
        correlation_matrix.at[col1, col2] = corr_coefficient
        correlation_matrix.at[col2, col1] = corr_coefficient

    # Fill diagonal with 1 since each column is perfectly correlated with itself
    corr = correlation_matrix.fillna(1.0)

    #print(corr)

    g = sb.clustermap(corr, cmap='coolwarm', figsize=(10, 8),
                   method=f'{clustering_method}', col_cluster=True, row_cluster=True,
                   annot=True, fmt=".2f", annot_kws={'size': 10})

    plt.savefig(f'fig/{clustering_method}_correlation_{level}_intersection_circpediaRNaseR_bsj{threshold}.png', dpi=300, bbox_inches='tight')
    plt.show()

