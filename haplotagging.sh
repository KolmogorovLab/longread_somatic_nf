#!/usr/bin/env bash

# To install WhatsHap using Conda

 conda create -n whatshap-env whatshap 
 conda activate whatshap-env
 
# Input and output variables
REF=/data/KolmogorovLab/references/grch38_chr.fasta
PHASED_VCF=/data/KolmogorovLab/julieZelaya/clair3_output.vcf
NORMAL_BAM=
TUMOR_BAM=

# To haplotag Normal and Tumor

whatshap haplotag --reference $REF $PHASED_VCF normal.bam -o normal.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=4

whatshap haplotag --reference $REF $PHASED_VCF tumor.bam -o tumor.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=4
