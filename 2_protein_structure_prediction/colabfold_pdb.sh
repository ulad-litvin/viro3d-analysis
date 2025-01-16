#!/bin/bash

for i in {0..8}
do
    time colabfold_batch virosphere-fold-v1_set_${i}_msas_org virosphere-fold-v1_set_${i}_pdbs
done