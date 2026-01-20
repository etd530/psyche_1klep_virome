#!/usr/bin/env bash
set -euxo pipefail # can add -x for debugging
IFS=$'\n\t'

#### VARS ####
reads=$1                                                # 0b_hifi_reads/acleris_abietana.ERR14886516.fq.gz
threads=10

# Ensure reads and genome are set
if [ -z "$reads" ] || [ -z "$threads" ]; then
  echo "Usage: $0 <hifi_reads.fastq.gz> <threads=INT>"
  exit 1
fi

# Get base name of reads
reads_pre=$(basename $reads)                            # acleris_abietana.ERR14886516.fq.gz
reads_simple=$(basename ${reads} | cut -f1,2 -d'.')     # acleris_abietana.ERR14886516

# Get species name
species=$(echo $reads_pre | cut -f1 -d'.')              # acleris_abietana

# Get genome file from species name
genome=0a_genomes/${species}*.fna.gz                    # 0a_genomes/acleris_abietana.GCA_965276485.1.sequences.renamed.fna.gz

# Get base name of genome
genome_pre=$(basename $genome)                          # acleris_abietana.GCA_965276485.1.sequences.renamed.fna.gz

#### STEP 1 #### NOT USED ANYMORE BC REAL VIRUSES MAY MAP TO EVE REGIONS OR TEs
# Map the HiFi reads to the corresponding reference genome, 
# and split the Lepidoptera reads from those that most likely are not Lepidoptera
# cd 1a_mapped_hifi_reads
# /usr/bin/time -o ${reads_pre}.time -v bash ./minimap2.hifi.sh ../$reads ../$genome
# cd ..

#### STEP 2 ####
# Assemble the non-Lep reads using viralFlye
cd 2a_dna_virome_assembly
/usr/bin/time -o ${reads_pre}.time -v bash ./viralflye.sh ../$reads $threads
if [ -s ${reads_simple}.viralFlye/circulars_viralFlye.fasta ]; then
    mv ${reads_simple}.viralFlye/circulars_viralFlye.fasta ${reads_simple}.viralFlye/${reads_simple}.circulars_viralFlye.fasta
fi

if [ -s ${reads_simple}.viralFlye/linears_viralFlye.fasta ]; then
    mv ${reads_simple}.viralFlye/linears_viralFlye.fasta ${reads_simple}.viralFlye/${reads_simple}.linears_viralFlye.fasta
fi
cd ..

#### STEP 3 ####
# Identify viral contigs with VirSorter2 (only if something was assembled)
cd 3_viral_contig_identification
if [ -s ../2a_dna_virome_assembly/${reads_simple}.viralFlye/${reads_simple}.circulars_viralFlye.fasta ]; then
    /usr/bin/time -o ${reads_pre}.time -v bash ./vs2.sh ../2a_dna_virome_assembly/${reads_simple}.viralFlye/${reads_simple}.circulars_viralFlye.fasta $threads
    if [ -s ${reads_simple}.circulars_viralFlye.vs2.out/final-viral-combined.fa ]; then
        mv ${reads_simple}.circulars_viralFlye.vs2.out/final-viral-combined.fa ${reads_simple}.circulars_viralFlye.vs2.out/${reads_simple}.circulars.fa
        sed -Ei 's/\|\|/\./g' ${reads_simple}.circulars_viralFlye.vs2.out/${reads_simple}.circulars.fa
        sed -Ei 's/\|\|/\./g' ${reads_simple}.circulars_viralFlye.vs2.out/final-viral-score.tsv
    fi
else
    echo -e "WARNING: No circular contigs assembled for ${reads_simple}"
fi

if [ -s ../2a_dna_virome_assembly/${reads_simple}.viralFlye/${reads_simple}.linears_viralFlye.fasta ]; then
    /usr/bin/time -o ${reads_pre}.time -v bash ./vs2.sh ../2a_dna_virome_assembly/${reads_simple}.viralFlye/${reads_simple}.linears_viralFlye.fasta $threads
    if [ -s ${reads_simple}.linears_viralFlye.vs2.out/final-viral-combined.fa ]; then
        mv ${reads_simple}.linears_viralFlye.vs2.out/final-viral-combined.fa ${reads_simple}.linears_viralFlye.vs2.out/${reads_simple}.linears.fa
        sed -Ei 's/\|\|/\./g' ${reads_simple}.linears_viralFlye.vs2.out/${reads_simple}.linears.fa
        sed -Ei 's/\|\|/\./g' ${reads_simple}.linears_viralFlye.vs2.out/final-viral-score.tsv
    fi
else
    echo -e "WARNING: No linear contigs assembled for ${reads_simple}"
fi

cd ..

#### STEP 4 ####
# QC the viral contigs with checkV
cd 4_checkv
if [ -s ../3_viral_contig_identification/${reads_simple}.circulars_viralFlye.vs2.out/${reads_simple}.circulars.fa ]; then
    /usr/bin/time -o ${reads_pre}.time -v bash ./checkv.sh ../3_viral_contig_identification/${reads_simple}.circulars_viralFlye.vs2.out/${reads_simple}.circulars.fa $threads
    grep High-quality ${reads_simple}.circulars.checkv/quality_summary.tsv | cut -f1,3 | grep 'No' | cut -f1 > ${reads_simple}.checkv_circulars.high_quality.viruses.txt || true       # adding '|| true' prevents empty grep from causing an error
    grep High-quality ${reads_simple}.circulars.checkv/quality_summary.tsv | cut -f1,3 | grep 'Yes' | cut -f1 > ${reads_simple}.checkv_circulars.high_quality.proviruses.txt || true
    samtools faidx -r ${reads_simple}.checkv_circulars.high_quality.viruses.txt ${reads_simple}.circulars.checkv/viruses.fna > ${reads_simple}.checkv_circulars.high_quality.viruses.fa || true # empty faidx also returns error
    samtools faidx -r ${reads_simple}.checkv_circulars.high_quality.proviruses.txt ${reads_simple}.circulars.checkv/proviruses.fna > ${reads_simple}.checkv_circulars.high_quality.proviruses.fa || true
else
    echo -e "WARNING: No circular contigs identified as viral by VirSorter2 for ${reads_simple}"
fi

if [ -s ../3_viral_contig_identification/${reads_simple}.linears_viralFlye.vs2.out/${reads_simple}.linears.fa ]; then
    /usr/bin/time -o ${reads_pre}.time -v bash ./checkv.sh ../3_viral_contig_identification/${reads_simple}.linears_viralFlye.vs2.out/${reads_simple}.linears.fa $threads
    grep High-quality ${reads_simple}.linears.checkv/quality_summary.tsv | cut -f1,3 | grep 'No' | cut -f1 > ${reads_simple}.checkv_linears.high_quality.viruses.txt || true       # adding '|| true' prevents empty grep from causing an error
    grep High-quality ${reads_simple}.linears.checkv/quality_summary.tsv | cut -f1,3 | grep 'Yes' | cut -f1 > ${reads_simple}.checkv_linears.high_quality.proviruses.txt || true
    samtools faidx -r ${reads_simple}.checkv_linears.high_quality.viruses.txt ${reads_simple}.linears.checkv/viruses.fna > ${reads_simple}.checkv_linears.high_quality.viruses.fa || true # empty faidx also returns error
    samtools faidx -r ${reads_simple}.checkv_linears.high_quality.proviruses.txt ${reads_simple}.linears.checkv/proviruses.fna > ${reads_simple}.checkv_linears.high_quality.proviruses.fa || true
else
    echo -e "WARNING: No linear contigs identified as viral by VirSorter2 for ${reads_simple}"
fi

cd ..

# #### STEP 5 ####
# # Align the viral contigs to the host's reference genome using MUMMER to ensure they are not EVEs
# cd 5_mummer

# # # Decompress genome
# ucomp_genome=$(basename ${genome} .gz)   # acentria_ephemerella*.fna
# gunzip -c ../${genome} > ${ucomp_genome}

# # Run for circular viruses if any exist
# if [ -s ../4_checkv/${reads_simple}.checkv_circulars.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./mummer.sh ../4_checkv/${reads_simple}.checkv_circulars.high_quality.viruses.fa ${ucomp_genome}
# else
#     echo -e "WARNING: No high-quality circular viral contigs for ${reads_simple} to screen against host genome"
# fi

# # Run for circular proviruses if any exist
# if [ -s ../4_checkv/${reads_simple}.checkv_circulars.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./mummer.sh ../4_checkv/${reads_simple}.checkv_circulars.high_quality.proviruses.fa ${ucomp_genome}
# else
#     echo -e "WARNING: No high-quality circular proviral contigs for ${reads_simple} to screen against host genome"
# fi

# # Run for linear viruses if any exist
# if [ -s ../4_checkv/${reads_simple}.checkv_linears.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./mummer.sh ../4_checkv/${reads_simple}.checkv_linears.high_quality.viruses.fa ${ucomp_genome}
# else
#     echo -e "WARNING: No high-quality linear viral contigs for ${reads_simple} to screen against host genome"
# fi

# # Run for linear proviruses if any exist
# if [ -s ../4_checkv/${reads_simple}.checkv_linears.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./mummer.sh ../4_checkv/${reads_simple}.checkv_linears.high_quality.proviruses.fa ${ucomp_genome}
# else
#     echo -e "WARNING: No high-quality linear proviral contigs for ${reads_simple} to screen against host genome"
# fi

# # Remove decompressed genome to save space
# rm ${ucomp_genome}
# cd ..

#### STEP 6 ####
# # BLAST the viral contigs against the NT database
# set +euxo pipefail      # disable strict error handling for conda activation
# source /home/software/miniforge3/etc/profile.d/conda.sh && conda activate paleofinder
# set -euxo pipefail

# cd 6b_nucleotide_blast
# # Run for linear viruses if any exist
# if [ -s ../4_checkv/${reads_simple}.checkv_linears.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_simple}.time -v bash ./blastn.sh ../4_checkv/${reads_simple}.checkv_linears.high_quality.viruses.fa /data/databases/nt/nt $threads
#     # Summarise BLAST results
#     blast_get_mrca.py ${reads_simple}.checkv_linears.high_quality.viruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239
# else
#     echo -e "WARNING: No high-quality linear free viruses to BLAST for ${reads_simple}"
# fi

# # Run for linear proviruses if any exist
# if [ -s ../4_checkv/${reads_simple}.checkv_linears.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_simple}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_linears.high_quality.proviruses.fa /data/databases/nt/nt $threads
#     # Summarise BLAST results
#     blast_get_mrca.py ${reads_simple}.checkv_linears.high_quality.proviruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239
# else
#     echo -e "WARNING: No high-quality linear proviruses to BLAST for ${reads_pre}"
# fi

# # Run for circular viruses if any exist
# if [ -s ../4_checkv/${reads_simple}.checkv_circulars.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_simple}.time -v bash ./blastn.sh ../4_checkv/${reads_simple}.checkv_circulars.high_quality.viruses.fa /data/databases/nt/nt $threads
#     # Summarise BLAST results
#     blast_get_mrca.py ${reads_simple}.checkv_circulars.high_quality.viruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239
# else
#     echo -e "WARNING: No high-quality circular free viruses to BLAST for ${reads_simple}"
# fi

# # Run for circular proviruses if any exist
# if [ -s ../4_checkv/${reads_simple}.checkv_circulars.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_simple}.time -v bash ./blastn.sh ../4_checkv/${reads_simple}.checkv_circulars.high_quality.proviruses.fa /data/databases/nt/nt $threads
#     # Summarise BLAST results
#     blast_get_mrca.py ${reads_simple}.checkv_circulars.high_quality.proviruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239
# else
#     echo -e "WARNING: No high-quality circular proviruses to BLAST for ${reads_simple}"
# fi

# cd ..


### STEP 6b ####
# BLAST the viral contigs against a nucleotide database of viral genomes
# cd 6c_nucleotide_blast_viruses_db
# # Run for linear viruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_linears.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_linears.high_quality.viruses.fa /data/databases/nt_viruses/nt_10239_viruses $threads
# else
#     echo -e "WARNING: No high-quality linear free viruses to BLAST for ${reads_pre}"
# fi
# # Summarise BLAST results
# blast_get_mrca.py ${reads_pre}.checkv_linears.high_quality.viruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239

# # Run for linear proviruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_linears.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_linears.high_quality.proviruses.fa /data/databases/nt_viruses/nt_10239_viruses $threads
# else
#     echo -e "WARNING: No high-quality linear proviruses to BLAST for ${reads_pre}"
# fi
# # Summarise BLAST results
# blast_get_mrca.py ${reads_pre}.checkv_trinity.high_quality.proviruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239

# # RUn for circular viruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_circulars.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_circulars.high_quality.viruses.fa /data/databases/nt_viruses/nt_10239_viruses $threads
# else
#     echo -e "WARNING: No high-quality circular free viruses to BLAST for ${reads_pre}"
# fi
# # Summarise BLAST results
# blast_get_mrca.py ${reads_pre}.checkv_circulars.high_quality.viruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239

# # Run for circular proviruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_circulars.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_circulars.high_quality.proviruses.fa /data/databases/nt_viruses/nt_10239_viruses $threads
# else
#     echo -e "WARNING: No high-quality circular proviruses to BLAST for ${reads_pre}"
# fi
# # Summarise BLAST results
# blast_get_mrca.py ${reads_pre}.checkv_circulars.high_quality.proviruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239
# cd ..

#### STEP 8 ####
# Classify the viral contigs with Metabuli
cd 9_metabuli

# Run for linear viruses if any exist
if [ -s ../4_checkv/${reads_simple}.checkv_linears.high_quality.viruses.fa ]; then
    /usr/bin/time -o ${reads_simple}.time -v bash ./metabuli.sh ../4_checkv/${reads_simple}.checkv_linears.high_quality.viruses.fa /data/databases/metabuli_refseq_viruses/refseq_virus/ $threads
else
    echo -e "WARNING: No high-quality linear free viruses to classify for ${reads_simple}"
fi

# Run for linear proviruses if any exist
if [ -s ../4_checkv/${reads_simple}.checkv_linears.high_quality.proviruses.fa ]; then
    /usr/bin/time -o ${reads_simple}.time -v bash ./metabuli.sh ../4_checkv/${reads_simple}.checkv_linears.high_quality.proviruses.fa /data/databases/metabuli_refseq_viruses/refseq_virus/ $threads
else
    echo -e "WARNING: No high-quality linear proviruses to classify for ${reads_simple}"
fi

# Run for circular viruses if any exist
if [ -s ../4_checkv/${reads_simple}.checkv_circulars.high_quality.viruses.fa ]; then
    /usr/bin/time -o ${reads_simple}.time -v bash ./metabuli.sh ../4_checkv/${reads_simple}.checkv_circulars.high_quality.viruses.fa /data/databases/metabuli_refseq_viruses/refseq_virus/ $threads
else
    echo -e "WARNING: No high-quality circular free viruses to classify for ${reads_simple}"
fi

# Run for circular proviruses if any exist
if [ -s ../4_checkv/${reads_simple}.checkv_circulars.high_quality.proviruses.fa ]; then
    /usr/bin/time -o ${reads_simple}.time -v bash ./metabuli.sh ../4_checkv/${reads_simple}.checkv_circulars.high_quality.proviruses.fa /data/databases/metabuli_refseq_viruses/refseq_virus/ $threads
else
    echo -e "WARNING: No high-quality circular proviruses to classify for ${reads_simple}"
fi
