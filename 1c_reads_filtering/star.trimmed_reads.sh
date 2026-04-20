#!/usr/bin/env bash

FS=$'\n\t'

source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate star2

set -euo pipefail

#### VARS ####
file1=$1                                               # ../1b_trimmed_rna_reads/acentria_ephemerella.ERR10123695.R1.trim.fq.gz
ref_genome=$2                                          # ../4_checkv/acentria_ephemerella.ERR10123695.checkv_trinity.high_quality.viruses.fa
threads=$3

sample=$(basename $file1 .R1.trim.fq.gz)
file2=${file1/.R1./.R2.}
genome_prefix=$(basename $ref_genome .fna.gz)

#### MAIN ####
if [ -e ${sample}.DONE ]; then
	echo "${sample}.DONE exists. Skipping run for ${sample}"
else
	# Map the RNA-Seq reads to the human genome
	STAR --runThreadN $threads --runMode alignReads --genomeDir /data/refgenomes/human_genomes/homo_sapiens.GCA_009914755.4_T2T-CHM13v2.0.genomefiles \
	--sjdbOverhang 150 --readFilesIn $file1 $file2 --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
	--outSAMmapqUnique 60 --outFileNamePrefix ${sample}.vs.${genome_prefix}. && \
	mv ${sample}.vs.${genome_prefix}.Aligned.sortedByCoord.out.bam ${sample}.vs.${genome_prefix}.bam && samtools index ${sample}.vs.${genome_prefix}.bam && \
	mkdir -p ${sample}.star.tmp && mv ${sample}.vs.${genome_prefix}* ${sample}.star.tmp && \

	# Output only reads that are unmapped and whose mate is also unmapped to the human genome
	samtools fastq -f 12 -@${threads} -1 ${sample}.non_human.R1.fq.gz -2 ${sample}.non_human.R2.fq.gz -0 /dev/null -s /dev/null -n ${sample}.star.tmp/${sample}.vs.${genome_prefix}.bam && touch ${sample}.DONE || touch ${sample}.FAILED
fi
