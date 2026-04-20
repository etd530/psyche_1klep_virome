#!/bin/bash

source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate fastp

##### VARS ####
file1=$1
threads=$2

#### MAIN ####
mkdir -p fastp_reports
prefix1=`basename $file1 | sed -E 's/\.fq\.gz//g'`;

if [ -s ${prefix1}.trim.fq.gz ]; then
    echo "Trimmed reads file for $prefix1 already exists. Skipping";
    else
        fastp -i ${file1} -o ${prefix1}.trim.fq.gz --trim_poly_x --detect_adapter_for_pe --cut_front --cut_tail --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --html fastp_reports/${prefix1}.html --json fastp_reports/${prefix1}.json --thread ${threads};
fi
