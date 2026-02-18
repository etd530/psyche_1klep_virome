#!/usr/bin/env bash

parallel --colsep '\t' -j 10 'pre=$(basename {} .R1.fq.gz) && /usr/bin/time -o ${pre}.alt.time -v bash ./main_rna.alt.sh {} > ${pre}.main.alt.stdout 2>${pre}.main.alt.stderr' :::: input_rna.alt.tsv
