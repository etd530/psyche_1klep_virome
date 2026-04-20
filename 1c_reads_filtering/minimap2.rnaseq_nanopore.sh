#!/usr/bin/env bash

FS=$'\n\t'

source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate minimap2

set -euo pipefail

#### VARS ####
file1=$1                                               # ../1b_trimmed_rna_reads/acentria_ephemerella.ERR10123695.R1.trim.fq.gz
ref_genome=$2                                          # ../4_checkv/acentria_ephemerella.ERR10123695.checkv_trinity.high_quality.viruses.fa
threads=$3

sample=$(basename $file1 .trim.fq.gz)
genome_prefix=$(basename $ref_genome .fna.gz)

#### MAIN ####
if [ -e ${sample}.DONE ]; then
	echo "${sample}.DONE exists. Skipping run for ${sample}"
else
	# Map the RNA-Seq reads to the human genome, we lower the kmer value a bit to increase sensitivity for noisy direct RNA-seq
	minimap2 -ax splice -k14 -t $threads -R "@RG\tID:${sample}\tSM:${sample}" $ref_genome $file1 | samtools sort -@${threads} -o ${sample}.vs.${genome_prefix}.bam && samtools index ${sample}.vs.${genome_prefix}.bam && \
	mkdir -p ${sample}.star.tmp && mv ${sample}.vs.${genome_prefix}* ${sample}.star.tmp && \

	# Output only reads that are unmapped to the human genome
	samtools fastq -f 4 -@${threads} -o ${sample}.non_human.fq.gz -0 /dev/null -s /dev/null -n ${sample}.star.tmp/${sample}.vs.${genome_prefix}.bam && touch ${sample}.DONE || touch ${sample}.FAILED
fi
