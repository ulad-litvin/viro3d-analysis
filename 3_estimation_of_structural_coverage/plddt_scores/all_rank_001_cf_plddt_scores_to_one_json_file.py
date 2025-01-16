'''
Put pLDDT scores for all rank 001 ColabFold models into one json file.
'''

import os
import json

def plddt_scores_from_CF_json_file(json_file_path: str):
    '''
    Extract pLDDT scores from a ColabFold JSON file
    and returns them as a list.
    '''
    with open(json_file_path, 'r') as f:
        data = json.load(f)
    plddt = data["plddt"]
    return plddt

colabfold_archive = "colabfold" # Path to the folder containing the ColabFold pLDDT scores in JSON format
directories = [f for f in os.listdir(colabfold_archive) if f.startswith("mars")]

n = 0
record_plddts = {}

for directory in directories:

    json_files = [f for f in os.listdir(os.path.join(colabfold_archive, directory))
                  if f.endswith(".json") and "rank_001" in f]
    
    for json_file in json_files:
        
        n += 1
        json_file_path = os.path.join(colabfold_archive, directory, json_file)

        record_id = json_file.split("_scores_rank_001")[0]
        plddt = plddt_scores_from_CF_json_file(json_file_path)

        record_plddts[record_id] = plddt

        print(n, record_id, plddt)


output_file = "./all_rank_001_CF_plddt_scores.json"

with open(output_file, 'w') as f:
    json.dump(record_plddts, f)