#!/usr/bin/env bash

parallel --colsep '\t' -j procfile 'pre=$(basename {} .fq.gz) && /usr/bin/time -o ${pre}.HIFI.time -v bash ./main.sh {} > ${pre}.main.HIFI.stdout 2>${pre}.main.HIFI.stderr' :::: input.tsv
