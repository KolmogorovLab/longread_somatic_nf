#!/usr/bin/env nextflow

/*
 * Define the workflow
 */
workflow {
    reads = Channel.fromPath(params.reads)
    ref = Channel.fromPath(params.ref)

    // Capture both the aligned BAM and BAI files from the align process
    (sortedBam, indexedBai) = align(ref, reads)

    // Pass the BAM file to the variantCalling process
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
    path 'aligned.bam.bai'

    script:
    """
    samtools fastq -TMm,M1,MM,ML ${reads} |
    minimap2 -ax map-ont -k 17 -t 30 -y --eqx ${ref} - |
    samtools sort -@30 -m 4G > aligned.bam
    samtools index aligned.bam
    """
}

/*
 * Process to run Clair3
 */
process variantCalling {
    label 'variantCalling'

    input:
    path ref
    path sortedBam 

    output:
    path 'clair3_output'

    script:
    """
    singularity exec \
    ${baseDir}/clair3_latest.sif \
    clair3 \
    --bam_fn=${sortedBam} \
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
    whatshap haplotag --reference $ref $phased_vcf $normal_bam
    """
}
