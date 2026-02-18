#!/usr/bin/env bash
set -euo pipefail # can add -x for debugging
IFS=$'\n\t'

# Initialize output file with header
# file=""
# for f in 3_viral_contig_identification/*/final-viral-score.tsv; do
#   if [ -e "$f" ]; then
#     file="$f"
#     break
#   fi
# done
# vs2_header=$(head -n1 $file)

file=""
for f in 4b_SPADES_checkv/*/quality_summary.tsv; do
  if [ -e "$f" ]; then
    file="$f"
    break
  fi
done
checkv_header=$(head -n1 $file)

echo -e "Species\tLibrary\tSource\tContig_ID\tseqname\tRNA\tNCLDV\tssDNA\tLavidaviridae\tNaldaviricetes\tmax_score\tmax_score_group\tlength\thallmark\tviral\tcellular\t${checkv_header}\tblastn_viral_mrca\tmetabuli_score\tmetabuli_classification" > psyche_viruses_summary.SPADES.tsv

# Iterate over contigs that were classified as high quality by checkV
for file in 4b_SPADES_checkv/*fa; do
	# Check file has some sequences in it
	if [ -s $file ]; then
		# Get the species and library name
		species=$(basename $file | cut -f1 -d'.')
		library=$(basename $file | cut -f2 -d'.')
		echo -e "${file}"
		# Determine if RNA or DNA based on filename
		if [[ $file == *spades* ]]; then
			source="RNA"
		else
			source="DNA"
		fi
		echo $source
		# List the contig IDs in the file and iterate over them
		contigs_list=$(grep '^>' $file | sed 's/>//g')
		for contig in $contigs_list; do
			vs2_line=$(grep "$contig\b" 3b_SPADES_viral_contig_identification/${species}.${library}*/final-viral-score.tsv | cut -f2 -d':')
			checkv_line=$(grep "$contig\b" 4b_SPADES_checkv/${species}.${library}*/quality_summary.tsv | cut -f2- -d':')
			if grep -q "$contig\b" 6b_nucleotide_blast/${species}.${library}*.vs.nt.out.mrcas.tsv; then
				blastn_mrca=$(grep "$contig\b" 6b_nucleotide_blast/${species}.${library}*.vs.nt.out.mrcas.tsv | cut -f2)
			else
				blastn_mrca="No blastn hits"
			fi
			if grep -q "$contig\b" 9b_metabuli_spades/${species}.${library}/${species}.${library}*viruses_classifications.tsv; then
				metabuli_score=$(grep "$contig\b" 9b_metabuli_spades/${species}.${library}/${species}.${library}*viruses_classifications.tsv | cut -f5)
				metabuli_class=$(grep "$contig\b" 9b_metabuli_spades/${species}.${library}/${species}.${library}*viruses_classifications.tsv | cut -f7)
			else
				metabuli_score="No metabuli score"
				metabuli_class="No metabuli classification"
			fi
			echo -e "${species}\t${library}\t${source}\t${contig}\t${vs2_line}\t${checkv_line}\t${blastn_mrca}\t${metabuli_score}\t${metabuli_class}" >> psyche_viruses_summary.SPADES.tsv
		done
	fi
done

# Verify that the number of contigs in the summary matches the number of high-quality contigs
summary_contigs=$(($(wc -l psyche_viruses_summary.SPADES.tsv | cut -f1 -d' ') - 1))
checkv_contigs=$(cat 4b_SPADES_checkv/*.fa | grep -c '^>')
