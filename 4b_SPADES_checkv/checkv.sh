#!/usr/bin/env bash

#### VARS ####
vs2_final_fa=$1
threads=$2
species=$(basename $vs2_final_fa .fa)

#### MAIN ####
source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate checkv
set -euo pipefail
checkv end_to_end $vs2_final_fa ${species}.checkv -t $threads -d /data/databases/checkv/checkv-db-v1.5 && touch ${species}.checkv.DONE || touch ${species}.checkv.FAILED
