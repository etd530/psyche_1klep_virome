#!/usr/bin/env bash

#### VARS ####
in_fasta=$1
out_fasta=$2
threads=$3

#### MAIN ####
source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate vclust

# Create a pre-alignment filter
vclust prefilter -i $in_fasta -o fltr.txt --min-ident 0.95 -t $threads

# Calculate ANI measures for genome pairs specified in the filter
vclust align -i $in_fasta -o ani.tsv --filter fltr.txt

# Cluster contigs into vOTUs using the MIUVIG thresholds (95% ANI and 85% MAF) and the Leiden algorithm
vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --algorithm leiden \
--metric ani --ani 0.95 --qcov 0.85 -out-repr
