#!/bin/bash

for i in {0..8}
do
    time esm-fold -i virosphere-fold-v1_set_${i}.fasta -o virosphere-fold-v1_set_${i}_esmfold_fold
done