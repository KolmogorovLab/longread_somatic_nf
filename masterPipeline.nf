#!/usr/bin/env nextflow

/*
 * Define the workflow
 */
workflow {
    reads = Channel.fromPath("/data/KolmogorovLab/julieZelaya/LongReadCancerPipeline/HCC1954BL.chr1.bam")
    ref = Channel.fromPath("/data/KolmogorovLab/references/grch38_chr.fasta")
    phased_vcf = Channel.fromPath("/data/KolmogorovLab/julieZelaya/LongReadCancerPipeline/clair3_output.vcf")
    
    sortedBam = align(ref, reads)
    variantCalling(ref, sortedBam)
}

/*
 * Align reads using minimap2 and sort BAM using samtools
 */
process align {
    label 'alignment'

    input:
    path ref
    path reads

    output:
    path 'aligned.bam'

    script:
    """
    samtools fastq -TMm,M1,MM,ML ${reads} |
    minimap2 -ax map-ont -k 17 -t 30 -y --eqx ${ref} - |
    samtools sort -@30 -m 4G > aligned.bam
    """
}

/*
 * Process to run Clair3
 */
process variantCalling {
    tag 'clair3'

    input:
    path ref
    path bam

    output:
    path 'clair3_output'

    script:
    """
    clair3 \
    --bam_fn=${bam} \
    --ref_fn=${ref} \
    --threads=20 \
    --platform="ont" \
    --model_path="/opt/models/ont" \
    --output=clair3_output \
    --enable_phasing \
    --longphase_for_phasing
    """
}

/*
 * Process to run Whatshap
 */

