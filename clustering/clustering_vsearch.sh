#!/usr/bin/env bash

#### VARS ####
in_fasta=$1
out_fasta=$2
threads=$3

#### MAIN ####
source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate ipyrad # this env has vsearch installed
# we use cluster_fast for initial length-based sorting, so centroids tend to be the longer (so probably more complete) contigs
vsearch --cluster_fast $in_fasta --id 0.95 --query_cov 0.85 --centroids $out_fasta --threads $threads --fasta_width 0 --clusters virus_clusters
