#!/bin/bash
set -eu

###Script for running DeepSomatic tumor-only on CPU
###Inteded to run inside the DeepSomatic Docker/Singularity container: google/deepsomatic:1.7.0
###e.g. singularity run -B /data deepsomatic_1.7.0.sif ./ds_paralel_tumor_normal.sh ARGS
###Should be submitted to exlusive Slurm instance (by adding --exclusive to the submission command)

if [ "$#" -ne 6 ]; then
    echo "Usage ds_parallel_tumor_only.sh reads_tumor reads_normal reference_fasta out_dir tumor_sample normal_sample"
    exit 1
fi

#Command line arguments
READS_TUMOR="$1"
READS_NORMAL="$2"
REF="$3"
OUT_DIR="$4"
TUMOR_SAMPLE="$5"
NORMAL_SAMPLE="$6"

#remove trailing slash if any, postprocess_variants does not like it
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

###DeepSomatic stages (you can get these command lines by running run_deepsomatic with --dry_run

#make examples
time seq 0 $((EXAMPLES_THREADS - 1)) | parallel -j ${EXAMPLES_THREADS} -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples_somatic \
    --mode calling --ref "${REF}" --reads_tumor "${READS_TUMOR}" --reads_normal "${READS_NORMAL}" \
    --examples "${EXAMPLES}/make_examples_somatic.tfrecord@${EXAMPLES_THREADS}.gz" \
    --checkpoint "/opt/models/deepsomatic/ont" --alt_aligned_pileup "diff_channels" \
    --min_mapping_quality "5" --parse_sam_aux_fields --partition_size "25000" --phase_reads \
    --pileup_image_width "99" --norealign_reads --sample_name_normal "${NORMAL_SAMPLE}" \
    --sample_name_tumor "${TUMOR_SAMPLE}" --sort_by_haplotypes --track_ref_reads \
    --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample "0.6" \
    --vsc_max_fraction_snps_for_non_target_sample "0.6" --vsc_min_fraction_indels "0.1" --vsc_min_fraction_snps "0.05" --task {}

#call variants
parallel --plus -j ${TF_JOBS} /opt/deepvariant/bin/call_variants \
    --outfile "${VAR_OUT}/call_variants_output_{0#}.tfrecord.gz" \
    --examples {} --checkpoint "/opt/models/deepsomatic/ont" \
    ::: `find "${EXAMPLES}" -name "make_examples*.gz"`

#postprocess variants
parallel --plus -j 1 /opt/deepvariant/bin/postprocess_variants \
    --ref "$REF" --infile "${VAR_OUT}/call_variants_output_{0#}.tfrecord.gz" \
    --outfile "${VCF_OUT}/postprocess_variants_{0#}.vcf.gz" --process_somatic=true \
    ::: `find "${VAR_OUT}" -name "*.gz"`


#concatenate vcf
bcftools concat -a ${VCF_OUT}/*.vcf.gz | bcftools sort | bgzip > "${OUT_DIR}/ds.merged.vcf.gz"
tabix "${OUT_DIR}/ds.merged.vcf.gz"
