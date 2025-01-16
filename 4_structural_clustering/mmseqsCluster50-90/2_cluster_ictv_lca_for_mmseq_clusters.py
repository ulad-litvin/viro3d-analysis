import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load ICTV taxonomy
file = "/home4/virosphere3D/ulad_analysis/tables/VMR_MSL38_v2.csv"
ictv_tax = pd.read_csv(file, index_col=0)
ictv_tax_sub = ictv_tax[['Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum',
                         'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', 'Family',
                         'Subfamily', 'Genus', 'Subgenus', 'Species']]
ictv_tax_sub.insert(0, 'Superkingdom', 'Viruses')

# Function to get lineage of a given taxonomy ID (mod to work with ICTV taxonomy)
"""
Created on Mon 23 Sep 11:24:38 BST 2019
Modified on Fri 4 Apr 13:14:24 BST 2024

@authors: sejmodha, ulad-litvin
"""

def get_rank(accession, ictv_sub_db=ictv_tax_sub):
    """A function to get the rank of a given ICTV taxon"""
    try:
        for col in ictv_sub_db.columns:
            if accession in ictv_sub_db[col].values:
                rank = col
    except ValueError:
        rank = None
        pass
    return rank

def get_ancestors(accession, ictv_sub_db=ictv_tax_sub):
    """A function to get lineage of a given taxonomy ID"""
    lineage = []
    try:
        for col in ictv_sub_db.columns:
            if accession in ictv_sub_db[col].values:
                lineage = ictv_sub_db[ictv_sub_db[col] == accession].iloc[0].dropna().to_list()
                index = lineage.index(accession)
                lineage = lineage[:index+1]
    except ValueError:
        lineage = None
        pass
    return lineage

def get_lca(taxon1,taxon2):
    """Function to get LCA between two taxonomy IDs"""
    if taxon1 is not None and taxon2 is not None:
        viralhits = 0
        ancestors_1 = get_ancestors(taxon1)[::-1]
        ancestors_2 = get_ancestors(taxon2)[::-1]
        for ancestor in ancestors_1:
            if ancestor in ancestors_2:
                return ancestor
    else:
        print(taxon1, taxon2)

def get_lca_list(taxalist):
    """A function to get LCA of a given list of taxonomy IDs"""
    taxon1 = taxalist.pop()
    while len(taxalist) > 0:
        taxon2 = taxalist.pop()
        lca = get_lca(taxon1, taxon2)
        taxon1 = lca
    return taxon1

# Get file paths from command line
parser = argparse.ArgumentParser(description='Process ICTV LCA for clusters.')
parser.add_argument('--stats', dest='file_rep', required=True, help='Path to the file_rep CSV file')
parser.add_argument('--meta', dest='file_mem', required=True, help='Path to the file_mem CSV file')
parser.add_argument('--out', dest='out_file', required=True, help='Path to the output CSV file')

args = parser.parse_args()

# load a file with info on cluster representatives
file_rep = args.file_rep
df_rep = pd.read_csv(file_rep)
df_rep

# load a file with taxonomy IDs of cluster members
file_mem = args.file_mem
df_mem = pd.read_csv(file_mem)
df_mem

# get LCA for each cluster
print('cluster_representative',
      'cluster_size',
      'ictv_lca',
      'ictv_lca_rank',
      sep='\t')

for index, row in df_rep.iterrows():

    clust_rep = row['cluster_representative']
    cluster_size = row['cluster_size']
    df_mem_sup = df_mem[df_mem['cluster_representative'] == clust_rep]

    taxalist = df_mem_sup['ictv_taxonomy'].tolist()

    if len(taxalist) > 1:
        df_rep.loc[index, 'ictv_lca'] = get_lca_list(taxalist)
    else:
        df_rep.loc[index, 'ictv_lca'] = taxalist[0]
    
    df_rep.loc[index, 'ictv_lca_rank'] = get_rank(df_rep.loc[index, 'ictv_lca'])
    
    print(clust_rep,
          cluster_size,
          df_rep.loc[index, 'ictv_lca'],
          df_rep.loc[index, 'ictv_lca_rank'],
          sep='\t')

# save the table with ICTV LCA for each cluster
out_file = args.out_file
df_rep.to_csv(out_file, index=False)

