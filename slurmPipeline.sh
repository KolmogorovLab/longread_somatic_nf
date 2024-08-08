#!/bin/bash
#SBATCH --job-name=nextflow_pipeline
#SBATCH --output=nextflow_pipeline.%j.out
#SBATCH --error=nextflow_pipeline.%j.err
#SBATCH --time=24:00:00 
#SBATCH --cpus-per-task=40
#SBATCH --mem=70G 
#SBATCH --partition=standard 

#Load modules 
module load java 
module load nextflow 

# Run Nextflow pipeline
nextflow run masterPipeline.nf
