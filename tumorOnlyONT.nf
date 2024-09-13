#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

/*
 * Main workflow
 */
workflow tumorOnlyWorkflow {
    take:
        reads
        reference

    main:
        alignMinimap2(reference, reads)
        callClair3(reference, alignMinimap2.out.bam, alignMinimap2.out.bam_idx, alignMinimap2.out.ref_idx)
        haplotagWhatshap(reference, alignMinimap2.out.ref_idx, callClair3.out.vcf, alignMinimap2.out.bam, 
                         alignMinimap2.out.bam_idx)

    emit:
        phasedVcf = callClair3.out.vcf
        haplotaggedBam = haplotagWhatshap.out.bam

    publish:
        phasedVcf >> "phased_vcf"
        haplotaggedBam >> "haplotagged_bam"
}

/*
 * Entry point
 */
workflow {
    if (!params.reads || !params.reference || !params.outdir) {
        error """
              ERROR: Some required arguments are not defined.
              Usage: tumorOnlyONT.nf --reads PATH --reference PATH --outdir PATH
              """.stripIndent()
    }

    tumorOnlyWorkflow(channel.from(params.reads), channel.from(params.reference))
}

output {
    directory params.outdir
    mode "copy"
}

//--- Individual processes description ---

/*
 * Align reads using minimap2, sort BAM using samtools, and create BAM index
 */
process alignMinimap2 {
    def threads = 28

    container 'docker://quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1'
    cpus threads
    memory '128 GB'
    time '24.h'

    input:
        path ref
        path reads

    output:
        path 'aligned.bam', emit: bam
        path 'aligned.bam.bai', emit: bam_idx
        path "${ref}.fai", emit: ref_idx
          
    script:
        """  
        samtools fastq -TMm,Ml,MM,ML ${reads} | minimap2 -ax map-ont -k 17 -t ${threads} -y --eqx ${ref} - | samtools sort -@4 -m 4G > aligned.bam
        samtools index -@8 aligned.bam
        samtools faidx ${ref}
        """
}   
    
/*
 * Process to run Clair3 
 */ 
process callClair3 {
    def threads = 28

    container 'docker://hkubal/clair3:v1.0.10'
    cpus threads
    memory '128 G'
    time '24.h'

    input:
        path ref
        path alignedBam
        path indexedBai
        path indexedRef 

    output:
        path 'clair3_output/merge_output.vcf.gz', emit: vcf

    
    script:
        """
        /opt/bin/run_clair3.sh \
            --bam_fn=${alignedBam} \
            --ref_fn=${indexedRef} \
            --threads=${threads} \
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
process haplotagWhatshap {
    container 'docker://quay.io/biocontainers/whatshap:2.3--py310h84f13bb_2'
    cpus 8
    memory '64 G'
    time '10.h'

    input:
        path ref
        path indexedRef
        path vcf
        path alignedBam
        path indexedBai

    output:
        path 'haplotagged.bam', emit: bam

    script:
        """
        tabix -@4 ${vcf}
        whatshap haplotag --reference ${ref} ${vcf} ${alignedBam} -o 'haplotagged.bam' --ignore-read-groups \
            --tag-supplementary --skip-missing-contigs --output-threads 4
        """
}

