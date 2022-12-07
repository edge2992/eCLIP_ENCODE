#!/bin/bash
# blastpを動かすスクリプトのサンプル
#
INPUT_FILE=/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed.fasta
DB_DIR=/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/blast/db

pushd $DB_DIR

makeblastdb -in $INPUT_FILE -dbtype prot \
  -out RbpHomoDB -parse_seqids

popd
