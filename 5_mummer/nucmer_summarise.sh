#!/usr/bin/env bash
set -euo pipefail # can add -x for debugging
IFS=$'\n\t'

for folder in nucmer*fa; do
	species=$(basename $folder | cut -f2 -d'.')
	echo $species
	library=$(basename $folder | cut -f3 -d'.')
	contigs_list=$(cut -f8 ${folder}/*nucmer_alignment.tsv | grep -v '[TAGS]' | wc -l)
	for contig in $contigs_list; do
		echo $contig
	done
done
