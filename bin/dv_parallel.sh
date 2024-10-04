#!/bin/bash
set -eu

###Script for running DeepVariant tumor-only on CPU
###Inteded to run inside the DeepSomatic Docker/Singularity container: google/deepvariant:1.6.1
###e.g. singularity run -B /data deepvariant.1.6.1.sif ./dv_paralel.sh ARGS
###Should be submitted to exlusive Slurm instance (by adding --exclusive to the submission command)

if [ "$#" -ne 4 ]; then
    echo "Usage dv_parallel.sh reads reference_fasta out_dir sample"
    exit 1
fi

#Command line arguments
READS="$1"
REF="$2"
OUT_DIR="$3"
SAMPLE="$4"

#remove training slash if any, postprocess_variants does not like it
OUT_DIR=${OUT_DIR%/}

#Change with caution. Optimized for a single exlusive node with 56 CPUs
#EXAMPLES_THREADS should be divisible by TF_JOBS
EXAMPLES_THREADS=56
TF_JOBS=4

INTERMEDIATE=${OUT_DIR}/intermediate
EXAMPLES=${INTERMEDIATE}/examples
VAR_OUT=${INTERMEDIATE}/variants
VCF_OUT=${INTERMEDIATE}/vcf
TMP_DIR=${INTERMEDIATE}/tmp
mkdir -p $EXAMPLES $VAR_OUT $VCF_OUT $TMP_DIR

#may be important to set TMP_DIR in case there are many files generated
export TMPDIR=${TMP_DIR}
#by default, biowulf only allwys 1024 file descriptors, which will not be enough, increase.
ulimit -u 10240 -n 16384

###DeepVariant stages (you can get these command lines by running run_deepvariant with --dry_run

#make examples
time seq 0 $((EXAMPLES_THREADS - 1)) | parallel -j ${EXAMPLES_THREADS} -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \
    --mode calling --ref "${REF}" --reads "${READS}" \
    --examples "${EXAMPLES}/make_examples.tfrecord@${EXAMPLES_THREADS}.gz" \
    --add_hp_channel --alt_aligned_pileup "diff_channels" --max_reads_per_partition "600" \
    --min_mapping_quality "5" --parse_sam_aux_fields --partition_size "25000" --phase_reads \
    --pileup_image_width "199" --norealign_reads --sort_by_haplotypes --track_ref_reads \
    --vsc_min_fraction_indels "0.12" --vsc_min_fraction_snps "0.08" --task {}

#call variants
parallel --link -j ${TF_JOBS} /opt/deepvariant/bin/call_variants \
    --outfile "${VAR_OUT}/call_variants_output_{2}.tfrecord.gz" \
    --examples {1} --checkpoint "/opt/models/ont_r104" \
    ::: `find "${EXAMPLES}" -name "make_examples*.gz"` ::: `seq -w 0 $((EXAMPLES_THREADS - 1))`

#postprocess variants
parallel -j 1 /opt/deepvariant/bin/postprocess_variants \
    --ref "$REF" --infile "${VAR_OUT}/call_variants_output_{}.tfrecord.gz" \
    --outfile "${VCF_OUT}/postprocess_variants_{}.vcf.gz" --cpus "1" --sample_name "${SAMPLE}" \
    ::: `seq -w 0 $((EXAMPLES_THREADS - 1))`

#concatenate vcf
bcftools concat -a ${VCF_OUT}/*.vcf.gz | bcftools sort | bgzip > "${OUT_DIR}/dv.merged.vcf.gz"
tabix "${OUT_DIR}/dv.merged.vcf.gz"
