#!/usr/bin/env nextflow

/*
 * Define the workflow
 */
workflow {
    reads = Channel.fromPath(params.reads)
    ref = Channel.fromPath(params.ref)

    // Capture both the aligned BAM and BAI files from the align process
    (alignedBam, indexedBai, indexedRef) = align(ref, reads)

    // Pass the BAM file to the variantCalling process
    vcf = variantCalling(ref, alignedBam)

    // Pass the VCF and BAM files to the phasing process
    phasing(ref, alignedBam, vcf)
}

/*
 * Align reads using minimap2, sort BAM using samtools, and create BAM index
 */
process align {
    label 'alignment'

    input:
    path ref
    path reads

    output:
    tuple path('aligned.bam'), path('aligned.bam.bai'), path('ref.bai')


    script:
    """
    samtools fastq -TMm,M1,MM,ML ${reads} |
    minimap2 -ax map-ont -k 17 -t 30 -y --eqx ${ref} - |
    samtools sort -@30 -m 4G -o aligned.bam -
    samtools index aligned.bam -
    samtools faidx ${ref}
    """
}

/*
 * Process to run Clair3
 */
process variantCalling {
    label 'variantCalling'

    input:
    path ref
    path alignedBam 

    output:
    path 'clair3_output.vcf'

    script:
    """
    run_clair3.sh \
    ${baseDir}/clair3_latest.sif \
    --bam_fn=${alignedBam} \
    --ref_fn=${ref} \
    --threads=20 \
    --platform="ont" \
    --model_path="/opt/models/ont" \
    --output=clair3_output.vcf \
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
    path alignedBam
    path 'clair3_output.vcf'

    output:
    path 'haplotagged.vcf'

    script:
    """
    whatshap haplotag --reference ${ref} clair3_output.vcf aligned.bam -o haplotagged.vcf

    """
}
