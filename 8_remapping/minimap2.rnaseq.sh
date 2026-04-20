#!/usr/bin/env bash

IFS=$'\n\t'

source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate minimap2

set -euo pipefail

#### VARS ####
file1=$1                                               # ../1c_trimmed_rna_reads/acentria_ephemerella.ERR10123695.non_human.R1.fq.gz
ref_genome=$2                                          # ../4_checkv/acentria_ephemerella.ERR10123695.checkv_trinity.high_quality.viruses.fa
threads=$3

#### MAIN ####
sample=$(basename $file1 .non_human.R1.fq.gz);              # acentria_ephemerella.ERR10123695
file2=`echo $file1 | sed -E 's/\.R1\./\.R2\./g'`;      # ../1b_trimmed_rna_reads/acentria_ephemerella.ERR10123695.R2.trim.fq.gz
genome_prefix=$(basename $ref_genome .fa) # | cut -f3,4,5 -d'.'); # checkv_trinity.high_quality.viruses

(
# Map the RNA-Seq reads to the viral contigs
# We do not check mapping in proper pair bc it is poorly defined and also makes less sense in RNA-Seq since these can be spliced and would not map in proper pair anyway so programs do not always set that flag for RNA-Seq I think.
# We also do not filter by mapping quality to avoid discarding multimapping reads due to repeat regions, as these do belong to the contig despite not being able to be placed with confidence.
# Also note that, compared to standard DNA sequencing, we do NOT remove duplicate reads as in RNA-Seq there tend to be lots of duplicate reads due to the high expression levels of some genes, so you remove real biological variation

minimap2 -ax splice:sr -t ${threads} -R "@RG\tID:${sample}\tSM:${sample}" $ref_genome $file1 $file2 | samtools view -bh -F4 | samtools sort -@${threads} -o ${sample}.vs.${genome_prefix}.bam && samtools index ${sample}.vs.${genome_prefix}.bam && \
conda activate mosdepth && mosdepth -t ${threads} ${sample}.vs.${genome_prefix} ${sample}.vs.${genome_prefix}.bam && touch ${sample}.vs.${genome_prefix}.DONE || touch ${sample}.vs.${genome_prefix}.FAILED
) > ${sample}.vs.${genome_prefix}.stdout 2>${sample}.vs.${genome_prefix}.stderr

# Delete the initial rm -f ${sample}*.bam && rm -f ${sample}*.bai
