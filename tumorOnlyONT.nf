#!/usr/bin/env nextflow

/*
 * Define the workflow
 */
workflow tumorOnlyWorkflow {
    take:
        reads
        ref

    main:
        (alignedBam, indexedBai, indexedRef) = align(ref, reads)
        (phasedVcf, vcfIndex) = variantCalling(ref, alignedBam, indexedBai, indexedRef)
        haplotaggedVcf = haplotag_bam(ref, indexedRef, phasedVcf, alignedBam, indexedBai, vcfIndex)

    emit:
        alignedBam
        phasedVcf
        haplotaggedVcf
}

workflow {
    tumorOnlyWorkflow(channel.from(params.reads), channel.from(params.ref))
}

/*
 * Align reads using minimap2, sort BAM using samtools, and create BAM index
 */
process align {
    label 'alignment'
    container 'docker://quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1'
    cpus 30
    memory '120 GB'
    time '24.h'

    input:
        path ref
        path reads

    output:
        path 'aligned.bam', emit: alignment
        path 'aligned.bam.bai', emit: alignment_index
        path "${ref}.fai", emit: reference_index
          
    script:
        """  
        samtools fastq -TMm,Ml,MM,ML ${reads} | minimap2 -ax map-ont -k 17 -t 30 -y --eqx ${ref} - | samtools sort -@4 -m 4G > aligned.bam
        samtools index -@8 aligned.bam
        samtools faidx ${ref}
        """
}   
    
/*
 * Process to run Clair3 
 */ 
process variantCalling {
    label 'variantCalling'
    container 'docker://hkubal/clair3:v1.0.10'
    cpus 30
    memory '120 G'
    time '24.h'

    input:
        path ref
        path alignedBam
        path indexedBai
        path indexedRef 

    output:
        path 'clair3_output/pileup.vcf.gz'
        path 'clair3_output/pileup.vcf.gz.tbi'

    
    script:
        """
        /opt/bin/run_clair3.sh \
            --bam_fn=${alignedBam} \
            --ref_fn=${indexedRef} \
            --threads=20 \
            --platform="ont" \
            --model_path="/opt/models/ont" \
            --output="clair3_output" \
            --enable_phasing \
            --longphase_for_phasing
        """
}

/*
 * Process to run Whatshap
 */
process haplotag_bam {
    label 'haplotag_bam'
    container 'docker://quay.io/biocontainers/whatshap:2.3--py310h84f13bb_2'
    cpus 10
    memory '64 G'
    time '10.h'

    input:
        path ref
        path indexedRef
        path vcf
        path alignedBam
        path indexedBai
        path vcfIndex

    output:
        path 'haplotagged.bam'

    script:
        """
        whatshap haplotag --reference ${ref} ${vcf} ${alignedBam} -o 'haplotagged.bam' --ignore-read-groups \
            --tag-supplementary --skip-missing-contigs --output-threads 4
        """
}

