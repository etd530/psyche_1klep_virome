#!/usr/bin/env bash

#### VARS ####
query=$1
db=$2
threads=$3

query_prefix=$(basename $query .fa)
db_prefix=$(basename $db)

#### MAIN ####
# Do diamond search of all 'good-quality' or unknown quality viruses obtained from CheckV against the Orthornavirae database of Katy
diamond blastx --very-sensitive -e 1e-5 --threads $threads --outfmt 6 qseqid qlen sallseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames --db $db --query $query --out ${query}.vs.${db}
