#!/usr/bin/env bash

IFS=$'\n\t'
source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate spades
set -euo pipefail

#### VARS ####
reads1=$1
threads=$2

reads2=${reads1/.R1./.R2.}
prefix=$(basename $reads1 .non_human.R1.fq.gz)


#### MAIN ####
if [ -e ${prefix}.DONE ]; then    # check .DONE file exists, which means the run finished OK and we don't have to re-run it
    echo -e "SPAdes assembly for ${prefix} already done, not re-runnning!";
else
    rnaviralspades.py -1 ${reads1} -2 ${reads2} -o ${prefix} --threads ${threads} && touch ${prefix}.DONE || touch ${prefix}.FAILED
fi
