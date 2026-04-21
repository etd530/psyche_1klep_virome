#!/usr/bin/env bash

#### VARS ####
reads=$1
threads=$2

reads_prefix=$(basename $reads .fq.gz)

#### MAIN ####
source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate viralflye

# Run metaFlye
flye --pacbio-hifi $reads --meta --threads $threads --out-dir ${reads_prefix}.metaFlye > ${reads_prefix}.metaFlye.out 2>${reads_prefix}.metaFlye.err && touch ${reads_prefix}.metaFlye.DONE || touch ${reads_prefix}.metaFlye.FAILED

viralFlye.py --dir ${reads_prefix}.metaFlye --hmm /data/databases/pfam-A/Pfam-A.hmm.gz --reads $reads --outdir ${reads_prefix}.viralFlye > ${reads_prefix}.viralFlye.out 2>${reads_prefix}.viralFlye.err && touch {reads_prefix}.viralFlye.DONE || touch {reads_prefix}.viralFlye.FAILED
