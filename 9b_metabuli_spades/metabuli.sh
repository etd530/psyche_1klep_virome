#!/usr/bin/env bash

source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate metabuli
set -euo pipefail
IFS=$'\n\t'

#### VARS ####
fasta=$1
dbdir=$2
threads=$3

prefix=$(basename $fasta .fa)
samplename=$(basename $fasta .fa | cut -f1,2 -d'.')

#### MAIN ####
metabuli classify $fasta $dbdir $samplename $prefix --threads $threads --max-ram 500 --seq-mode 3 --lineage 1 && touch ${prefix}.metabuli.DONE || touch ${prefix}.metabuli.FAILED
