#!/usr/bin/env bash
set -euxo pipefail # can add -x for debugging
IFS=$'\n\t'

#### VARS ####
reads=$1             # 0c_rnaseq/acentria_ephemerella.ERR10123695.R1.fq.gz
threads=10

# Ensure reads and genome are set
if [ -z "$reads" ] || [ -z "$threads" ]; then
  echo "Usage: $0 <hifi_reads.fastq.gz> <threads=INT>"
  exit 1
fi

# Get base name of reads
reads_base=$(basename $reads)                # acentria_ephemerella.ERR10123695.R1.fq.gz

# Get reads prefix
reads_pre=$(basename $reads .R1.fq.gz)       # acentria_ephemerella.ERR10123695

# Get species name
species=$(echo $reads_pre | cut -f1 -d'.')   # acentria_ephemerella

# Get libary name
library=$(echo $reads_pre | cut -f2 -d'.')    # ERR10123695

# Get genome file from species name
genome=$(ls 0a_genomes/${species}*.fna.gz)   # 0a_genomes/acentria_ephemerella.GCA_943193645.1.sequences.renamed.fna.gz

# Get base name of genome
genome_pre=$(basename $genome)               # acentria_ephemerella.GCA_943193645.1.sequences.renamed.fna.gz

#### STEP 1 ####
# Trim the RNA-Seq reads
# cd 1b_trimmed_rna_reads
# /usr/bin/time -o ${reads_pre}.time -v bash ./fastp.sh ../$reads $threads # output: acentria_ephemerella.ERR10123695.R1.trim.fq.gz
# cd ..

#### STEP 1.5 ####
# Map the reads to the host genome to filter host reads
cd 1c_rna_host_filtering
# /usr/bin/time -o ${reads_pre}.time -v bash ./minimap2.trimmed_reads.sh ../1b_trimmed_rna_reads/${reads_pre}.R1.trim.fq.gz ../${genome} $threads
/usr/bin/time -o TEST.time -v bash ./minimap2.trimmed_reads.sh /data/rnaseq/psyche_lepidoptera/pararge_aegeria.ERR6363260.RENAMED.R1.fq.gz ../0a_genomes/pararge_aegeria.GCA_905163445.1.sequences.renamed.fna.gz 10
cd ..
