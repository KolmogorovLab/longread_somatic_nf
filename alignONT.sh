#!/usr/bin/env bash
# Command to submit job: sbatch minimap2.sh


# Batch job
set -e

# Load the needed modules
module load minimap2
module load samtools

# Input and output variables
REF=/data/KolmogorovLab/references/grch38_chr.fasta
READS=/data/KolmogorovLab/CellLinesNanoQ20/nano_raw_bams/1945.bam
ALIGNED_BAM_OUTPUT=/data/KolmogorovLab/julieZelaya/1945_aligned.bam

# Align ONT long reads with methylation annotation
samtools fastq -TMm,Ml,MM,ML $READS | minimap2 -ax map-ont -k 17 -t 30 -y --eqx $REF - | samtools sort -@4 -m 4G $ALIGNED_BAM_OUTPUT





