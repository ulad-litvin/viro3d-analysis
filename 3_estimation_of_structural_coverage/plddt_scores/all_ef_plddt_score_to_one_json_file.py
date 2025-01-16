'''
Extract pLDDT scores from all ESMFold models and save them in one JSON file.
'''

import sys
import os
import json
from Bio import SeqIO
from biopandas.pdb import PandasPdb

def extract_plddt_scores_from_pdb_file(pdb_path: str):
    '''
    Extract the pLDDT scores from a PDB file
    and returns them as a list.
    '''
    
    # Load the PDB file
    ppdb = PandasPdb().read_pdb(pdb_path)

    # Extract the pLDDT scores from the CA (C-alpha) B-factor column
    plddt = ppdb.df['ATOM']['b_factor'].tolist()
    plddt =ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == 'CA']['b_factor'].tolist()
    
    return plddt

pdb_folder = "esmfold" # Path to the folder containing the ESMFold models in PDB format
pdb_file_names = [f for f in os.listdir(pdb_folder) if f.endswith('.pdb')]
pdb_file_paths = [os.path.join(pdb_folder, f) for f in pdb_file_names]

record_plddts = {}

n = 0

for pdb_file_path in pdb_file_paths:

    pdb_id = os.path.basename(pdb_file_path)[:-4]
    
    plddt = extract_plddt_scores_from_pdb_file(pdb_file_path)

    record_plddts[pdb_id] = plddt

    n += 1
    print(n, pdb_id, plddt)

output_file = "./all_ESMFold_plddt_scores.json"
with open(output_file, 'w') as f:
    json.dump(record_plddts, f)