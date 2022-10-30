#!/bin/bash
#
#

/home/edge2992/tools/subread-2.0.3-Linux-x86_64/bin/featureCounts \
  -T 5 \
  -a /home/edge2992/Resource/gencode.v24.chr_patch_hapl_scaff.basic.annotation.gtf \
  -o counts.txt \
  /home/edge2992/H/MYWORK/eCLIP_ENCODE/data/eCLIP/AARS/K562/ENCFF128BDG.bed

  # -a /home/edge2992/Resource/gencode.v24.long_noncoding_RNAs.gtf \
