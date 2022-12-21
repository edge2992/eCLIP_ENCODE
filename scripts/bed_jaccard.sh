#!/bin/bash

# sortされている必要がある
# デフォルトのファイルはソートされていない
DATADIR="/home/edge2992/H/MYWORK/eCLIP_ENCODE/data/eCLIP"
TARGET1="${DATADIR}/AARS/K562/ENCFF684FVK.bed"
TARGET2="${DATADIR}/AATF/K562/ENCFF349HTU.bed"

SORTEDDIR="/home/edge2992/H/MYWORK/eCLIP_ENCODE/scripts/tmp/"

SORTED_TARGET1="${SORTEDDIR}$(basename $TARGET1 .bed).sorted.bed"
SORTED_TARGET2="${SORTEDDIR}$(basename $TARGET2 .bed).sorted.bed"
echo $SORTED_TARGET1
sort -k1,1 -k2,2n $TARGET1 > $SORTED_TARGET1
sort -k1,1 -k2,2n $TARGET2 > $SORTED_TARGET2


bedtools jaccard \
  -a $SORTED_TARGET1 \
  -b $SORTED_TARGET2 



# カウントする場合
# intersectBed \
#   -a /home/edge2992/H/MYWORK/eCLIP_ENCODE/data/eCLIP/AARS/K562/ENCFF684FVK.bed \
#   -b /home/edge2992/Resource/gencode.v24.annotation.gene.gtf \
#   -wa -c -s -loj > tmp/count_ENCFF684FVK_s.bed


