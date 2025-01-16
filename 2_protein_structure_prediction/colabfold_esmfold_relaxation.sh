#!/bin/bash

FOLDER="./cf_ef_rank1/"

for FILE in `ls $FOLDER`
do
    python relaxation_v4.py $FOLDER$FILE 3600
done
