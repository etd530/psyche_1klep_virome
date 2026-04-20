#!/usr/bin/env bash

zcat SRR525187_1.trim.fq.gz | sed -E "s/^@.+\.([0-9]+) .+$/@SRR525187:105:HWKHMDRXX:1:2101:\1:1031 1:N:0:NNNNNNNNNN+NNNNNNNNNN/g" > SRR525187_1.trim.renamed.fq
zcat SRR525187_2.trim.fq.gz | sed -E "s/^@.+\.([0-9]+) .+$/@SRR525187:105:HWKHMDRXX:1:2101:\1:1031 1:N:0:NNNNNNNNNN+NNNNNNNNNN/g" > SRR525187_2.trim.renamed.fq
gzip SRR525187_1.trim.renamed.fq
gzip SRR525187_2.trim.renamed.fq
