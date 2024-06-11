#! /bin/bash
#SBATCH --job-name=nf_longReadCancerPipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --gres=lscratch:200
#SBATCH --time=24:00:00

module load nextflow
export NXF_SINGULARITY_CACHEDIR=/data/$USER/nxf_singularity_cache;
export SINGULARITY_CACHEDIR=/data/$USER/.singularity;
export TMPDIR=/lscratch/$SLURM_JOB_ID
export NXF_JVM_ARGS="-Xms2g -Xmx4g"

nextflow run nf-core/rnaseq -r 3.13.2 \
-profile biowulf \
