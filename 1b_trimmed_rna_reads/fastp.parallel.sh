#!/usr/bin/env bash
parallel -j 2 --colsep '\t' './fastp.sh {1} {2}' :::: fastp.parallel.tsv
