#!/usr/bin/env bash

parallel --colsep '\t' -j 1 "/usr/bin/time -o {1}.time -v bash ./mummer.sh {1} {2}" :::: mummer.parallel.tsv
