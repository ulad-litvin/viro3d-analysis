import argparse
import numpy as np
import pandas as pd

# get file paths from command line
parser = argparse.ArgumentParser(description='Create combined meta and stats files for mmseq foldseek clusters.')
parser.add_argument('--mmseq_stats', dest='file_mmseq_stats', required=True, help='Path to the mmseq_stats CSV file')
parser.add_argument('--mmseq_meta', dest='file_mmseq_meta', required=True, help='Path to the mmseq_meta CSV file')
parser.add_argument('--foldseek', dest='file_foldseek_results', required=True, help='Path to the foldseek easy-cluster TSV file')
parser.add_argument('--out', dest='out_file', required=True, help='Path to the output CSV files')

args = parser.parse_args()

# read mmseq cluster metadata csv file
file_mmseq_meta = args.file_mmseq_meta
mmseq_meta = pd.read_csv(file_mmseq_meta)

# read mmseq cluster stats data csv file
file_mmseq_stats = args.file_mmseq_stats
mmseq_stats = pd.read_csv(file_mmseq_stats)

# read foldseek easy-cluster tsv result file
file_foldseek_results = args.file_foldseek_results
col_names = ["cluster_representative", "cluster_member"]
foldseek_results = pd.read_csv(file_foldseek_results, sep="\t", names=col_names)

'''
STEP 1: Create a meta dataframe for combined mmseq+foldseek clusters
Meta dataframe contains one row for each cluster member
'''

# create a dictionary to store the combined metadata
combined_dict = {'cluster_representative':[],
                 'cluster_member':[],
                 'clustering_method':[]}

# iterate through foldseek cluster results
for index, row in foldseek_results.iterrows():
    
    # get the cluster representative and member
    cluster_foldseek_rep = row['cluster_representative']
    cluster_foldseek_member = row['cluster_member'][3:-8]

    # add the cluster representative and member to the combined metadata dictionary
    combined_dict['cluster_representative'].append(cluster_foldseek_rep)
    combined_dict['cluster_member'].append(cluster_foldseek_member)
    # add the foldseek clustering method to the combined metadata dictionary for these cluster members
    combined_dict['clustering_method'].append("foldseek")

    # for members that represent mmseq clusters
    # add mmseq cluster members to the combined metadata dictionary
    if cluster_foldseek_member in mmseq_stats['structure_representative_record_id'].unique():

        # use mmseq structure representative to get the mmseq cluster representative
        cluster_mmseq_rep = mmseq_stats[mmseq_stats['structure_representative_record_id'] == cluster_foldseek_member]['cluster_representative'].values[0]
        # use mmseq cluster representative to get all mmseq cluster members
        mmseq_meta_sub = mmseq_meta[mmseq_meta['cluster_representative'] == cluster_mmseq_rep]
        
        # add mmseq cluster members to the combined metadata dictionary 
        # add the mmseq clustering method
        for index2, row2 in mmseq_meta_sub.iterrows():
            
            cluster_mmseq_member = row2['cluster_member']
            
            if cluster_mmseq_member != cluster_foldseek_member:
                
                combined_dict['cluster_representative'].append(cluster_foldseek_rep)
                combined_dict['cluster_member'].append(cluster_mmseq_member)
                combined_dict['clustering_method'].append("mmseq")

# convert dictionary into a dataframe
foldseek_meta = pd.DataFrame(combined_dict)

# add additional columns from mmseq meeta
foldseek_meta = pd.merge(foldseek_meta, mmseq_meta.drop(['cluster_representative', 'cluster_id', 'pfam_desc'], axis=1), how='left', left_on='cluster_member', right_on='cluster_member')

# add cluster_member column
foldseek_meta['cluster_member'] = foldseek_meta['member_class'] + "-" + foldseek_meta['member_record_id'] + "_relaxed"

print("step 1 done")

'''
STEP 2: Create a stats dataframe from the meta dataframe
stats dataframe contains one row for each cluster
'''

# count the number of cluster_members for each cluster_representative, sort by count and make it a dataframe
foldseek_stats = foldseek_meta["cluster_representative"].value_counts().reset_index()
foldseek_stats.columns = ["cluster_representative", "cluster_size"]

# for each cluster_representative, get protlen_mean, plddt_mean, ptm_mean from df
cluster_means = foldseek_meta[['cluster_representative', 'protlen', 'plddt', 'ptm']].groupby("cluster_representative").mean().reset_index()

# merge cluster_counts and cluster_means dataframes
foldseek_stats = pd.merge(foldseek_stats, cluster_means, on="cluster_representative")

# round the values to 2 decimal places
foldseek_stats["protlen_mean"] = foldseek_stats["protlen"].round(2)
foldseek_stats["plddt_mean"] = foldseek_stats["plddt"].round(2)
foldseek_stats["ptm_mean"] = foldseek_stats["ptm"].round(2)
foldseek_stats = foldseek_stats.drop(columns=['protlen', 'plddt', 'ptm'])
foldseek_stats.columns = ["cluster_representative", "cluster_size", "protlen_mean", "plddt_mean", "ptm_mean"]

# add cluster_id column to the stats dataframe
foldseek_stats['cluster_id'] = foldseek_stats.index + 1

print("step 2 done")

'''
STEP 3: Find structure representative and most frequent genbank_name for each cluster
Structure representative is the structure with the highest pLDDT score from this cluster
'''

# find structure representative and most frequent genbank_name for each cluster
for index, row in foldseek_stats.iterrows():
    
    cluster_id = row["cluster_id"]
    cluster_rep = row["cluster_representative"]
    cluster_size = row["cluster_size"]
    
    # make a subset of foldseek_meta for this cluster
    df_sub = foldseek_meta[foldseek_meta["cluster_representative"] == cluster_rep].sort_values(by=["plddt", "ptm"], ascending=False)
    
    # find the structure representative (structure with the highest plddt score)
    structure_rep = df_sub.iloc[0]["cluster_member"]
    structure_rep_plddt = df_sub.iloc[0]["plddt"]
    
    # add the structure representative info to foldseek_stats
    foldseek_stats.loc[index, "structure_representative"] = structure_rep
    foldseek_stats.loc[index, "structure_representative_plddt"] = structure_rep_plddt

    # find the most frequent genbank_name in the cluster
    try:
        top_genbank_name = foldseek_meta[foldseek_meta['cluster_representative'] == cluster_rep]["genbank_name"].value_counts().index[0]
    except IndexError:
        top_genbank_name = np.nan
    
    # add the most frequent genbank_name to foldseek_stats
    foldseek_stats.loc[index, 'top_genbank_name'] = top_genbank_name

print("step 3 done")

'''
STEP 4: Save the meta and stats files
'''

# add cluster_id and structure_representative columns to the meta dataframe
foldseek_meta = pd.merge(foldseek_meta, foldseek_stats[['cluster_representative', 'cluster_id', 'structure_representative']], on='cluster_representative')
foldseek_meta.sort_values(by=['cluster_id'], ascending=True, inplace=True)
foldseek_meta.reset_index(drop=True, inplace=True)

# specify the column order for the meta dataframe
foldseek_meta = foldseek_meta[[
    "cluster_id",
    "cluster_representative",
    "structure_representative",
    "cluster_member",
    "clustering_method",
    "member_class",
    "member_record_id",
    "protlen",
    "plddt",
    "ptm",
    "ictv_sort",
    "ictv_species",
    "genbank_name"]]

# specify the column order for the stats dataframe
foldseek_stats = foldseek_stats[[
    "cluster_id",
    "cluster_size",
    "protlen_mean",
    "plddt_mean",
    "ptm_mean",
    "cluster_representative",
    "structure_representative",
    "structure_representative_plddt",
    "top_genbank_name"]]

# save the meta and stats dataframes to csv files
foldseek_meta.to_csv(args.out_file + "_meta.csv", index=False)
foldseek_stats.to_csv(args.out_file + "_stats.csv", index=False)

print("step 4 done")