#!/usr/bin/env bash
parallel -j 2 --colsep '\t' './minimap2.trimmed_reads.sh {1} {2} 10' :::: minimap2.trimmed_reads.parallel.tsv
