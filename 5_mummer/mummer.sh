#!/usr/bin/env bash

#### VARS ####
refgenome=$1
vircontigs=$2

prefix=$(basename $refgenome .sequences.fna)

#### MAIN ####
source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate mummer4
set -euo pipefail

# Run NUCMER
mkdir -p nucmer.${prefix}
nucmer ${refgenome} ${vircontigs} -t 15 --prefix ${prefix} && \
delta-filter -q ${prefix}.delta > ${prefix}.filter && \
# Get delta into coordinates
show-coords -rT ${prefix}.filter > ${prefix}.coords && \
tail -n +4 *.coords > ${prefix}.nucmer_alignment.tsv && \
mv ${prefix}* nucmer.${prefix}
touch nucmer.${prefix}.DONE || touch nucmer.${prefix}.FAILED
