#!/usr/bin/env bash

# Allocate interactive job in biowulf
sinteractive --gres=gpu:p100:1,lscratch:10 --mem=16g -c4 --time=8:00:00

# Load the singularity module 
module load singularity

# Pull Clair3 Docker and convert to singularity
singularity pull docker://google/deepvariant:latest
docker pull kishwars/pepper_deepvariant:latest

# Input and output variables
REF=/data/KolmogorovLab/references/grch38_chr.fasta
INPUT_DIR=/data/KolmogorovLab/julieZelaya/1945_aligned.bam
OUTPUT_DIR=/data/KolmogorovLab/julieZelaya/clair3_output
THREADS=4
MODEL_NAME=ont

# Run Clair3
# Generates phased VCF Using Clair3

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR} \    
  --ref_fn=$REF \       
  --threads=${THREADS} \               
  --platform="ont" \                   
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR}.vcf \
  --enable_phasing \
  --longphase_for_phasing


# Exit the interactive session
exit
