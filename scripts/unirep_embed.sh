#!/bin/bash

INPUT_FASTA="/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed.fasta"
OUTPUT_FILE="/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_unirep.npz"

pipenv run \
  tape-embed unirep \
  ${INPUT_FASTA} \
  ${OUTPUT_FILE} \
  babbler-1900 \
  --batch_size 2 \
  --tokenizer unirep
