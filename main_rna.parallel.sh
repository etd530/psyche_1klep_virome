#!/usr/bin/env bash

parallel --colsep '\t' -j procfile 'pre=$(basename {} .R1.fq.gz) && /usr/bin/time -o log/${pre}.time -v bash ./main_rna.sh {} > log/${pre}.main.stdout 2>log/${pre}.main.stderr' :::: input_rna.tsv
