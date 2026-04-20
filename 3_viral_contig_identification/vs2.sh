#!/usr/bin/env bash

source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate vs2

#### VARS ####
mav=$1
threads=$2
species=`basename ${mav} .fasta.gz`
echo $species

#### MAIN ####
# we do not set --min-length since some viruses (e.g. Circoviridae) have very small genomes, only
# we use the --keep-original-seq flag since later we will run checkV
# we use rm-tmpdir all to remove all temporary files since there are many
# WARNING: we do not use --prep-for-dramv for now BUT it could be interesting in the future to annotate the viral genomes with DRAMv
# We include NCLDV (Nucleocytoviricota or nucleocytoplasmic large DNA viruses), which includes Ascoviridae which are known from Lepidoptera
# We include ssDNA (single-stranded DNA viruses) which is a general category and should catch some insect viruses?
# We include Lavidaviridae which are virophages that infect giant viruses and may also occur in insects? I found no examples for now

(
virsorter run -w ${species}.vs2.out -i ${mav} --include-groups "RNA,NCLDV,ssDNA,Lavidaviridae,Naldaviricetes" -j $threads --keep-original-seq --rm-tmpdir all && touch ${species}.vs2.DONE || touch ${species}.vs2.FAILED
) > ${species}.vs2.stdout 2>${species}.vs2.err
