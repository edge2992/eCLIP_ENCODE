#!/bin/bash

INPUT_FASTA="/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed.fasta"
OUTPUT_FILE="/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_tape.npz"

pipenv run \
  tape-embed transformer \
  ${INPUT_FASTA} \
  ${OUTPUT_FILE} \
  bert-base \
  --batch_size 2 \
  --tokenizer iupac
