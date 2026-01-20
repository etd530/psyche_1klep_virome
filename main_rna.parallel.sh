#!/usr/bin/env bash

parallel --colsep '\t' -j 15 'pre=$(basename {} .R1.fq.gz) && /usr/bin/time -o ${pre}.time -v bash ./main_rna.sh {} > ${pre}.main.stdout 2>${pre}.main.stderr' :::: input_rna.tsv
