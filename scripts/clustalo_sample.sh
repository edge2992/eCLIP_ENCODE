#!/bin/bash
# clustaloを動かすスクリプトのサンプル
INPUTFILE=./data/uniprot/reviewed.fasta
LABEL=reviewed_output

echo $INPUTFILE

clustalo \
    --seqtype protein \
    --MAC-RAM 32000 \
    --distmat-out  ${LABEL}.distmat \
    --infile ${INPUTFILE} --threads 10 \
    --full \
    --verbose \
    --outfmt clustal \
    --resno \
    --output-order input-order \
    --outfile ${LABEL}.clustal_num \
