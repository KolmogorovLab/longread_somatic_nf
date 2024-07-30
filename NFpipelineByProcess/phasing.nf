#!/usr/bin/env nextflow

/* 
 *Define the workflow for haplotagging/phasing
 */
workflow {
    ref = Channel.fromPath("/data/KolmogorovLab/references/grch38_chr.fasta")
    phased_vcf = Channel.fromPath("/data/KolmogorovLab/julieZelaya/LongReadCancerPipeline/clair3_output.vcf")
    normal_bam = Channel.fromPath("/data/KolmogorovLab/julieZelaya/LongReadCancerPipeline/1954BL.haplotagged.alignment.bam)
    tumor_bam = Channel.fromPath("/data/KolmogorovLab/julieZelaya/LongReadCancerPipeline/1954.haplotagged.alignment.bam)

    phasing(ref, phased_vcf, normal_bam, tumor_bam)
}

/*
 * Process to run Whatshap
 */
process phasing {
    label 'phasing'

    input:
    path ref
    path phased_vcf
    path normal_bam
    path tumor_bam

    output:
    path 'haplotagged'

    script:
    """
    whatshap haplotag --reference $ref $phased_vcf $normal_bam -o 1954BL.haplotagged.bam --ignore-read-groups --tag-suppplementary --skip-missing-contigs --output-threads=4
    """
}
