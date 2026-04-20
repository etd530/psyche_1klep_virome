#!/bin/bash

source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate minimap2

#### VARS ####
file1=$1
ref_genome=$2

#### MAIN ####
prefix1=`basename $file1 _1.trim.fq.gz`;
sample=`echo $prefix1 | sed -E 's/\.R[12]//g'`;
file2=`echo $file1 | sed -E 's/_1/_2/g'`;

(
# Get the non-Lepidoptera reads
minimap2 -ax splice:sr -t 20 -R "@RG\tID:${sample}\tSM:${sample}" $ref_genome $file1 $file2 | samtools sort -@20 -o ${sample}.bam && samtools index ${sample}.bam && \
samtools view -bh -f 12 -o ${sample}.UnUn.bam ${sample}.bam && `# both reads unmapped` \
samtools view -bh -f 4 -F 8 -o ${sample}.UnMa.bam ${sample}.bam && `# reads unmapped with mate mapped` \
samtools view -bh -f 8 -F 4 -o ${sample}.MaUn.bam ${sample}.bam && `# reads mapped with mate unmapped` \
samtools merge -o ${sample}.nonLep.bam ${sample}.UnUn.bam ${sample}.UnMa.bam ${sample}.MaUn.bam && \
samtools sort -@20 -n -o ${sample}.nonLep.sorted.bam ${sample}.nonLep.bam && \
samtools fastq -@20 ${sample}.nonLep.sorted.bam -1 ${sample}.nonLep.R1.fastq.gz -2 ${sample}.nonLep.R2.fastq.gz -n && \
touch ${sample}.minimap2.DONE || touch ${sample}.minimap2.FAILED

# Get the Lepidoptera reads and assess coverage with Mosdepth
samtools view -bh -f2 -q20 -F4 ${sample}.bam | samtools sort -@20 -o ${sample}.Lep_reads.bam && samtools index ${sample}.Lep_reads.bam && \
conda activate sambamba && \
sambamba markdup --overflow-list-size 600000 -t 20 ${sample}.Lep_reads.bam ${sample}.Lep_reads.MD.bam && \
conda activate minimap2 && samtools index ${sample}.Lep_reads.MD.bam && \
conda activate mosdepth && mosdepth -t 20 ${sample}.Lep_reads.MD ${sample}.Lep_reads.MD.bam && touch ${sample}.mosdepth.DONE || touch ${sample}.mosdepth.FAILED
) > ${sample}.minimap2.stdout 2>${sample}.minimap2.stderr

# rm -f ${sample}*.bam && rm -f ${sample}*.bai
