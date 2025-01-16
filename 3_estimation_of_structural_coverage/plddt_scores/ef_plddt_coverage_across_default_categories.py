'''
Calculate the average pLDDT score and pLDDT coverage across default categories for all ESMFold/ ColabFold models.
'''

import pandas as pd
import json
import statistics

def plddt_mean(plddt: list, round_to: int = 1):
    '''
    Calculate the average pLDDT score.
    '''
    return round(statistics.mean(plddt), round_to)

def plddt_coverage_aa_number(plddt: list, lower_threshold: float,
                             upper_threshold: float, include_upper_threshold: bool = False):
    '''
    Calculate the coverage of a protein sequence in number of AA
    based on a pLDDT thresholds and return the number of AA.
    '''
    if include_upper_threshold:
        return len([i for i in plddt if i >= lower_threshold and i <= upper_threshold])
    else:
        return len([i for i in plddt if i >= lower_threshold and i < upper_threshold])

def plddt_coverage_across_default_categories(plddt: list):
    '''
    Calculate the coverage of a protein sequence in number of AA
    based on default pLDDT categories: [0, 50), [50, 70), [70, 90), [90, 100] -
    and return a list with the number of AA for each category.
    '''
    coverage = []
    coverage.append(plddt_coverage_aa_number(plddt, 0, 50))
    coverage.append(plddt_coverage_aa_number(plddt, 50, 70))
    coverage.append(plddt_coverage_aa_number(plddt, 70, 90))
    coverage.append(plddt_coverage_aa_number(plddt, 90, 100, include_upper_threshold=True))
    return coverage

file_meta = "../virus_proteomes/virus_records.csv"
df_meta = pd.read_csv(file_meta)

ef_json_file = "./all_ESMFold_plddt_scores.json"
with open(ef_json_file, 'r') as f:
    ef_plddt = json.load(f)

n = 0
print("index\trecord_id\tef_plddt_mean\tef_plddt_0_50\tef_plddt_50_70\tef_plddt_70_90\tef_plddt_90_100")
for record_id, plddt in ef_plddt.items():
    n += 1
    plddt_cov = plddt_coverage_across_default_categories(plddt)
    plddt_mean_value = plddt_mean(plddt)
    df_meta.loc[df_meta['header'] == record_id, 'ef_pLDDT_mean'] = plddt_mean_value
    df_meta.loc[df_meta['header'] == record_id, 'ef_pLDDT_0_50'] = plddt_cov[0]
    df_meta.loc[df_meta['header'] == record_id, 'ef_pLDDT_50_70'] = plddt_cov[1]
    df_meta.loc[df_meta['header'] == record_id, 'ef_pLDDT_70_90'] = plddt_cov[2]
    df_meta.loc[df_meta['header'] == record_id, 'ef_pLDDT_90_100'] = plddt_cov[3]
    print(n, record_id, plddt_mean_value, plddt_cov[0], plddt_cov[1], plddt_cov[2], plddt_cov[3], sep='\t')

df_meta.to_csv("../virus_proteomes/virus_records_ef_plddt.csv", index=False)