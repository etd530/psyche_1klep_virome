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

#### STEP 2 #####
# # Assemble the trimmed reads using Trinity
# cd 2c_rna_assembly
# /usr/bin/time -o ${reads_pre}.time -v bash ./trinity.sh ../1b_trimmed_rna_reads/${reads_pre}.R1.trim.fq.gz $threads
# rm -r ${reads_pre}.trinity/ # remove the Trinity assembly folder to save space, we keep only the final FASTA file
# seqtk rename ${species}.${library}.trinity.Trinity.fasta ${species}.${library}.Trinity_ > ${species}.${library}.Trinity.tmp && mv ${species}.${library}.Trinity.tmp ${species}.${library}.Trinity.fasta
# cd ..
# #### output: acentria_ephemerella.ERR10123695.trinity
# ####         acentria_ephemerella.ERR10123695.trinity.Trinity.fasta

#### STEP 2b ####
# Assemble with rnaviralSPAdes
# cd 2d_rna_assembly_rnaviralSPAdes
# /usr/bin/time -o ${reads_pre}.time -v bash ./rnaviralSPAdes.sh ../1b_trimmed_rna_reads/${reads_pre}.R1.trim.fq.gz $threads
# cd ..

#### STEP 3 ####
# # Identify viral contigs with VirSorter2 (only if something was assembled)
# cd 3_viral_contig_identification
# if [ -s ../2c_rna_assembly/${reads_pre}.Trinity.fasta ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./vs2.sh ../2c_rna_assembly/${reads_pre}.Trinity.fasta $threads
#     if [ -s ../3_viral_contig_identification/${reads_pre}.Trinity.vs2.out/final-viral-combined.fa ]; then
#         mv ${reads_pre}.Trinity.vs2.out/final-viral-combined.fa ../3_viral_contig_identification/${reads_pre}.Trinity.vs2.out/${reads_pre}.Trinity.fa
#         sed -Ei 's/\|\|/\./g' ${reads_pre}.Trinity.vs2.out/${reads_pre}.Trinity.fa
#         sed -Ei 's/\|\|/\./g' ${reads_pre}.Trinity.vs2.out/final-viral-score.tsv
#     fi
# else
#     echo -e "WARNING: No contigs assembled for ${reads_pre}"
# fi
# cd ..

#### output: acentria_ephemerella.ERR10123695.Trinity.vs2.out/final-viral-combined.fa

#### STEP 3b ####
# # Identify viral contigs with VirSorter2 from SPADES output (only if something was assembled)
# cd 3b_SPADES_viral_contig_identification
# if [ -s ../2d_rna_assembly_rnaviralSPAdes/${reads_pre}/${reads_pre}.rnaviralspades.fasta ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./vs2.sh ../2d_rna_assembly_rnaviralSPAdes/${reads_pre}/${reads_pre}.rnaviralspades.fasta $threads
# else
#     echo -e "WARNING: No contigs assembled for ${reads_pre}"
# fi
# cd ..

#### STEP 4 ####
# QC the viral contigs with checkV
# cd 4_checkv
# if [ -s ../3_viral_contig_identification/${reads_pre}.Trinity.vs2.out/${reads_pre}.Trinity.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./checkv.sh ../3_viral_contig_identification/${reads_pre}.Trinity.vs2.out/${reads_pre}.Trinity.fa $threads
#     if [ -s ${reads_pre}.Trinity.checkv/proviruses.fna ]; then
#         sed -Ei 's/_1 .+$//g' ${reads_pre}.Trinity.checkv/proviruses.fna;
#     fi
#     # Keep contigs classified as high-quality, with at least one viral gene, and without warnings
#     awk -F"\t" '{if (($8 == "High-quality" || $9 == "High-quality") && $14 == "" && $6 >= 1 && $3 == "No") print $1}' ${reads_pre}.Trinity.checkv/quality_summary.tsv > ${reads_pre}.checkv_trinity.high_quality.viruses.txt
#     awk -F"\t" '{if (($8 == "High-quality" || $9 == "High-quality") && $14 == "" && $6 >= 1 && $3 == "Yes") print $1}' ${reads_pre}.Trinity.checkv/quality_summary.tsv > ${reads_pre}.checkv_trinity.high_quality.proviruses.txt
#     samtools faidx -r ${reads_pre}.checkv_trinity.high_quality.viruses.txt ${reads_pre}.Trinity.checkv/viruses.fna > ${reads_pre}.checkv_trinity.high_quality.viruses.fa || true # empty faidx also returns error
#     samtools faidx -r ${reads_pre}.checkv_trinity.high_quality.proviruses.txt ${reads_pre}.Trinity.checkv/proviruses.fna > ${reads_pre}.checkv_trinity.high_quality.proviruses.fa || true
# else
#     echo -e "WARNING: No contigs identified as viral by VirSorter2 for ${reads_pre}"
# fi
# cd ..

#### STEP 4b ####
# QC the viral contigs from SPAdes with checkV
# cd 4b_SPADES_checkv
# if [ -s ../3b_SPADES_viral_contig_identification/${reads_pre}.rnaviralspades.vs2.out/${reads_pre}.rnaviralspades.vs2.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./checkv.sh ../3b_SPADES_viral_contig_identification/${reads_pre}.rnaviralspades.vs2.out/${reads_pre}.rnaviralspades.vs2.fa $threads
#     if [ -s ${reads_pre}.rnaviralspades.vs2.checkv/proviruses.fna ]; then
#         sed -Ei 's/_1 .+$//g' ${reads_pre}.rnaviralspades.vs2.checkv/proviruses.fna;
#     fi
#     # Keep contigs classified as high-quality, with at least one viral gene, and without warnings
#     awk -F"\t" '{if (($8 == "High-quality" || $9 == "High-quality") && $14 == "" && $6 >= 1 && $3 == "No") print $1}' ${reads_pre}.rnaviralspades.vs2.checkv/quality_summary.tsv > ${reads_pre}.checkv_spades.high_quality.viruses.txt
#     awk -F"\t" '{if (($8 == "High-quality" || $9 == "High-quality") && $14 == "" && $6 >= 1 && $3 == "Yes") print $1}' ${reads_pre}.rnaviralspades.vs2.checkv/quality_summary.tsv > ${reads_pre}.checkv_spades.high_quality.proviruses.txt
#     samtools faidx -r ${reads_pre}.checkv_spades.high_quality.viruses.txt ${reads_pre}.rnaviralspades.vs2.checkv/viruses.fna > ${reads_pre}.checkv_spades.high_quality.viruses.fa || true # empty faidx also returns error
#     samtools faidx -r ${reads_pre}.checkv_spades.high_quality.proviruses.txt ${reads_pre}.rnaviralspades.vs2.checkv/proviruses.fna > ${reads_pre}.checkv_spades.high_quality.proviruses.fa || true
# else
#     echo -e "WARNING: No contigs identified as viral by VirSorter2 for ${reads_pre}"
# fi
# cd ..

#### STEP 5 ####
# Remap the RNA-Seq reads to the viral contigs using minimap2
# cd 8_remapping
# # Remap to free viruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa ]; then
#     contig_list=$(grep '^>' ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa | sed 's/^>//g')
#     for contig in $contig_list; do
#         samtools faidx ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa ${contig} > ${contig}.fa
#         /usr/bin/time -o ${contig}.time -v bash ./minimap2.rnaseq.sh ../1b_trimmed_rna_reads/${reads_pre}.R1.trim.fq.gz ${contig}.fa $threads
#         rm ${contig}.fa
#         breadth=$(tail -n2 ${reads_pre}.vs.${contig}.mapq20.mosdepth.global.dist.txt | head -n1 | cut -f3)
#         if (( $(echo "$breadth == 1" | bc -l) )); then
#             echo $contig >> ${reads_pre}.checkv_trinity.high_quality.viruses.breadth100.txt
#         fi
#         rm ${reads_pre}.vs.${contig}*.bam
#         rm ${reads_pre}.vs.${contig}*.bai
#         samtools faidx -r ${reads_pre}.checkv_trinity.high_quality.viruses.breadth100.txt ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa > ${reads_pre}.checkv_trinity.high_quality.viruses.breadth100.fa || true
#     done
    

# else
#     echo -e "WARNING: No high-quality free viruses to remap for ${reads_pre}"
# fi 

# # Remap to proviruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa ]; then
#     contig_list=$(grep '^>' ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa | sed 's/^>//g')
#     for contig in $contig_list; do
#         samtools faidx ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa ${contig} > ${contig}.fa
#         /usr/bin/time -o ${contig}.time -v bash ./minimap2.rnaseq.sh ../1b_trimmed_rna_reads/${reads_pre}.R1.trim.fq.gz ${contig}.fa $threads
#         rm ${contig}.fa
#         breadth=$(tail -n2 ${reads_pre}.vs.${contig}.mapq20.mosdepth.global.dist.txt | head -n1 | cut -f3)
#         if (( $(echo "$breadth == 1" | bc -l ) )); then
#             echo $contig >> ${reads_pre}.checkv_trinity.high_quality.proviruses.breadth100.txt
#         fi
#         rm ${reads_pre}.vs.${contig}.*.bam
#         rm ${reads_pre}.vs.${contig}.*.bai
#         samtools faidx -r ${reads_pre}.checkv_trinity.high_quality.proviruses.breadth100.txt ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa > ${reads_pre}.checkv_trinity.high_quality.proviruses.breadth100.fa || true
#     done
# else
#     echo -e "WARNING: No high-quality proviruses to remap for ${reads_pre}"
# fi
# cd ..

#### STEP 6 ####
# Cluster the viral contigs to remove the ones that are likely the same virus
# cd clustering
# # Make file with all viruses to cluster
# cat ../8_remapping/*fa > viruses_all_pre-clustering.fa
# ./clustering_vsearch.sh viruses_all_pre-clustering.fa viruses_cluster.greedy_097.fasta $threads

# #### STEP 7 ####
# # Align the de-clustered viral contigs to the host's reference genome using MUMMER to ensure they are not EVEs
# cd 5_mummer
# for contig in $(grep '^>' ../clustering/viruses_cluster.greedy_097.fasta | sed 's/^>//g'); do
#     samtools faidx ../clustering/viruses_cluster.greedy_097.fasta ${contig} > ${contig}.fa
#     sp=$(echo $contig | cut -f1-d'.')
#     genome=../0a_genomes/${sp}*.fna.gz
#     # Decompress genome
#     ucomp_genome=$(basename ${genome} .gz)   # acentria_ephemerella*.fna
#     gunzip -c ../${genome} > ${ucomp_genome}
#     # Run MUMMER
#     /usr/bin/time -o ${reads_pre}.time -v bash ./mummer.sh ${contig}.fa ${ucomp_genome}
# done
# # Remove decompressed genome to save space
# rm *fna
# cd ..

#### STEP 8 ####
# # Classify the viral contigs with Metabuli
# cd 9_metabuli
# # Run for free viruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./metabuli.sh ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa /data/databases/metabuli_refseq_viruses/refseq_virus/ $threads
# else
#     echo -e "WARNING: No high-quality free viruses to classify for ${reads_pre}"
# fi

# # Run for proviruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./metabuli.sh ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa /data/databases/metabuli_refseq_viruses/refseq_virus/ $threads
# else
#     echo -e "WARNING: No high-quality proviruses to classify for ${reads_pre}"
# fi
# cd ..

#### STEP 8b ####
# # Classify the viral contigs from SPADES with Metabuli
cd 9b_metabuli_spades
# Run for free viruses if any exist
if [ -s ../4b_SPADES_checkv/${reads_pre}.checkv_spades.high_quality.viruses.fa ]; then
    /usr/bin/time -o ${reads_pre}.time -v bash ./metabuli.sh ../4b_SPADES_checkv/${reads_pre}.checkv_spades.high_quality.viruses.fa /data/databases/metabuli_refseq_viruses/refseq_virus/ $threads
else
    echo -e "WARNING: No high-quality free viruses to classify for ${reads_pre}"
fi

# Run for proviruses if any exist
if [ -s ../4b_SPADES_checkv/${reads_pre}.checkv_spades.high_quality.proviruses.fa ]; then
    /usr/bin/time -o ${reads_pre}.time -v bash ./metabuli.sh ../4b_SPADES_checkv/${reads_pre}.checkv_spades.high_quality.proviruses.fa /data/databases/metabuli_refseq_viruses/refseq_virus/ $threads
else
    echo -e "WARNING: No high-quality proviruses to classify for ${reads_pre}"
fi
cd ..

#### STEP 9 ####
# Extract the proteins predicted by checkV for the contigs and do a tblastn against a viral NT
# WARNING: BLASTP IS SUPER SLOW, SHOULD CHANGE TO DIAMOND OR WILL TAKE FOREVER
# cd 6a_tblastn
# mkdir -p ${species}

# # Run for free viruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa ]; then
#     virus_list=$(grep '^>' ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa | sed 's/^>//g')
#     for virus in $virus_list; do
#         grep "${virus}" ../4_checkv/${reads_pre}.Trinity.checkv/tmp/proteins.faa | sed -E 's/^>//g' | cut -f1 -d' ' > ${species}/${species}.${virus}.protein_ids.txt
#         samtools faidx -r ${species}/${species}.${virus}.protein_ids.txt ../4_checkv/${reads_pre}.Trinity.checkv/tmp/proteins.faa > ${species}/${species}.${virus}.proteins.faa
#         # /usr/bin/time -o ${species}.${virus}.time -v bash ./tblastn.sh ${species}/${species}.${virus}.proteins.faa /data/databases/nt_viruses/nt.10239_viruses $threads
#         /usr/bin/time -o ${species}.${virus}.blastp.time -v bash ./blastp.sh ${species}/${species}.${virus}.proteins.faa /data/databases/nr_viruses/nr.10239_viruses $threads
#     done
# else
#     echo -e "WARNING: No high-quality free viruses to tblastn for ${reads_pre}"
# fi

# # Run for proviruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa ]; then
#     provirus_list=$(grep '^>' ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa | sed 's/^>//g')
#     for provirus in $provirus_list; do
#         grep "${provirus}" ../4_checkv/${reads_pre}.Trinity.checkv/tmp/proteins.faa | sed -E 's/^>//g' | cut -f1 -d' ' > ${species}/${species}.${provirus}.protein_ids.txt
#         samtools faidx -r ${species}/${species}.${provirus}.protein_ids.txt ../4_checkv/${reads_pre}.Trinity.checkv/tmp/proteins.faa > ${species}/${species}.${provirus}.proteins.faa
#         # /usr/bin/time -o ${species}.${provirus}.time -v bash ./tblastn.sh ${species}/${species}.${provirus}.proteins.faa /data/databases/nt_viruses/nt.10239_viruses $threads
#         /usr/bin/time -o ${species}.${provirus}.blastp.time -v bash ./blastp.sh ${species}/${species}.${provirus}.proteins.faa /data/databases/nr_viruses/nr.10239_viruses $threads

#     done
# else
#     echo -e "WARNING: No high-quality proviruses to tblastn for ${reads_pre}"
# fi
# cd ..

#### STEP 10 ####
# # Annotate domains of the viral proteins with interproscan
# cd 11_interproscan
# for fasta in ../6a_tblastn/${species}/*.proteins.faa; do
#     cp $fasta .
#     protein_base=$(basename $fasta)
#     sed -Ei 's/\*//g' $protein_base
#     protein_pre=$(basename $fasta .proteins.faa)
#     /usr/bin/time -o ${protein_pre}.time -v bash ./interproscan.sh $protein_base $threads
# done








################################################################################################################

# #### STEP 6 #### PROBABLY NOT NEEDE SINCE WHOLE NT IS BAD FOR ID OF VIRUSES, ESPECIALLY PROVIRUSES
# # BLAST the viral contigs against the NT database
# cd 6b_nucleotide_blast
# # Run for free viruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa /data/databases/nt/nt $threads
# else
#     echo -e "WARNING: No high-quality free viruses to BLAST for ${reads_pre}"
# fi

# # Summarise BLAST results
# blast_get_mrca.py ${reads_pre}.checkv_trinity.high_quality.viruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239

# # Run for proviruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa /data/databases/nt/nt $threads
# else
#     echo -e "WARNING: No high-quality proviruses to BLAST for ${reads_pre}"
# fi
# # Summarise BLAST results
# blast_get_mrca.py ${reads_pre}.checkv_trinity.high_quality.proviruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239
# cd ..

#### STEP 6b #### MAYBE NOT NEEDED EITHER SINCE WE HAVE METABULI FOR IDENTIFICATION
# BLAST the viral contigs against a nucleotide database of viral genomes
# cd 6c_nucleotide_blast_viruses_db
# # Run for free viruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_trinity.high_quality.viruses.fa /data/databases/nt_viruses/nt_10239_viruses $threads
# else
#     echo -e "WARNING: No high-quality free viruses to BLAST for ${reads_pre}"
# fi

# # Summarise BLAST results
# blast_get_mrca.py ${reads_pre}.checkv_trinity.high_quality.viruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239

# # Run for proviruses if any exist
# if [ -s ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa ]; then
#     /usr/bin/time -o ${reads_pre}.time -v bash ./blastn.sh ../4_checkv/${reads_pre}.checkv_trinity.high_quality.proviruses.fa /data/databases/nt_viruses/nt_10239_viruses $threads
# else
#     echo -e "WARNING: No high-quality proviruses to BLAST for ${reads_pre}"
# fi
# # Summarise BLAST results
# blast_get_mrca.py ${reads_pre}.checkv_trinity.high_quality.proviruses.vs.nt.out /data/users/etode/paleofinder/PaleoFinder --colnum 15 --taxid 10239
# cd ..
