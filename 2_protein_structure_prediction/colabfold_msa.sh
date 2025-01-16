#!/bin/bash

for i in {0..8}
do
    time colabfold_search --mmseqs ~/anaconda3/envs/mmseqs/bin/mmseqs virosphere-fold-v1_set_${i}.fasta /db/colabfold/database/ virosphere-fold-v1_set_${i}_msas --threads 64
    cp -r virosphere-fold-v1_set_${i}_msas virosphere-fold-v1_set_${i}_msas_org
    python3 rename_a3m.py virosphere-fold-v1_set_${i}.fasta virosphere-fold-v1_set_${i}_msas_org
done