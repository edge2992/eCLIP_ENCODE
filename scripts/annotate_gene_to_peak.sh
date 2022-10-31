#!/bin/bash

intersectBed \
  -a /home/edge2992/H/MYWORK/eCLIP_ENCODE/data/eCLIP/AARS/K562/ENCFF684FVK.bed \
  -b /home/edge2992/Resource/gencode.v24.annotation.gene.gtf \
  -wa -s -loj > tmp/annotated_ENCFF684FVK_s.bed

# カウントする場合
# intersectBed \
#   -a /home/edge2992/H/MYWORK/eCLIP_ENCODE/data/eCLIP/AARS/K562/ENCFF684FVK.bed \
#   -b /home/edge2992/Resource/gencode.v24.annotation.gene.gtf \
#   -wa -c -s -loj > tmp/count_ENCFF684FVK_s.bed


