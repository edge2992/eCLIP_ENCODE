#!/bin/bash
# blastpを動かすスクリプトのサンプル
#
INPUT_FILE=/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed.fasta
DB_NAME=/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/blast/db/RbpHomoDB
OUTPUT_FILE=/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/blast/reviewed.tsv

blastp \
  -num_threads 6 \
  -query $INPUT_FILE \
  -db $DB_NAME \
  -out $OUTPUT_FILE \
  -outfmt "6 qseqid sseqid pident evalue bitscore score" \
  -evalue 1000
