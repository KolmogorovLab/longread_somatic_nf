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
        vntrAnnotation
        svPanelOfNormals
        clair3Model

    main:
        alignMinimap2(reference, reads)
        callClair3(alignMinimap2.out.bam, alignMinimap2.out.bam_idx, reference, alignMinimap2.out.ref_idx, clair3Model)
        phaseLongphase(alignMinimap2.out.bam, alignMinimap2.out.bam_idx, reference, 
                       alignMinimap2.out.ref_idx, callClair3.out.vcf)
        haplotagWhatshap(reference, alignMinimap2.out.ref_idx, phaseLongphase.out.phasedVcf, alignMinimap2.out.bam, 
                         alignMinimap2.out.bam_idx)
        severusTumorOnly(haplotagWhatshap.out.bam, haplotagWhatshap.out.bam_idx, callClair3.out.vcf, 
                         vntrAnnotation, svPanelOfNormals)

    emit:
        phasedVcf = callClair3.out.vcf
        haplotaggedBam = haplotagWhatshap.out.bam
        severusSomaticVcf = severusTumorOnly.out.severusSomaticVcf

    publish:
        phasedVcf >> "phased_vcf"
        haplotaggedBam >> "haplotagged_bam"
        severusSomaticVcf >> "severus"
}

/*
 * Entry point
 */
workflow {
    if (!params.reads || !params.reference || !params.outdir || 
        !params.vntr || !params.sv_pon || !params.clair3_model) {
        error """
              ERROR: Some required arguments are not defined.
              Usage: tumorOnlyONT.nf --reads PATH --reference PATH --outdir PATH 
                                     --vntr PATH --sv_pon PATH --clair3_model PATH
              """.stripIndent()
    }

    tumorOnlyWorkflow(Channel.fromPath(params.reads), Channel.fromPath(params.reference), 
                      Channel.fromPath(params.vntr), Channel.fromPath(params.sv_pon),
                      Channel.fromPath(params.clair3_model))
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
    
process callClair3 {
    def threads = 28

    container 'docker://hkubal/clair3:v1.0.10'
    cpus threads
    memory '128 G'
    time '24.h'

    input:
        path alignedBam
        path indexedBai
        path reference
        path referenceIdx
        path modelPath

    output:
        path 'clair3_output/merge_output.vcf.gz', emit: vcf

    
    script:
        """
        /opt/bin/run_clair3.sh \
            --bam_fn=${alignedBam} \
            --ref_fn=${reference} \
            --threads=${threads} \
            --platform="ont" \
            --model_path=${modelPath} \
            --output="clair3_output" \
        """
}

process phaseLongphase {
    def threads = 10

    container 'docker://mkolmogo/longphase:1.7.3'
    cpus threads
    memory '64 G'
    time '4.h'

    input:
        path alignedBam
        path indexedBai
        path reference
        path referenceIdx
        path vcf

    output:
        path 'longphase.vcf.gz', emit: phasedVcf

    
    script:
        """
        longphase phase -s ${vcf} -b ${alignedBam} -r ${reference} -t ${threads} -o longphase --ont
        bgzip longphase.vcf
        """
}

process haplotagWhatshap {
    container 'docker://mkolmogo/whatshap:2.3'
    cpus 8
    memory '64 G'
    time '10.h'

    input:
        path reference
        path referenceIdx
        path phasedVcf
        path alignedBam
        path indexedBai

    output:
        path 'haplotagged.bam', emit: bam
        path 'haplotagged.bam.bai', emit: bam_idx

    script:
        """
        tabix ${phasedVcf}
        whatshap haplotag --reference ${reference} ${phasedVcf} ${alignedBam} -o 'haplotagged.bam' --ignore-read-groups \
            --tag-supplementary --skip-missing-contigs --output-threads 4
        samtools index -@8 haplotagged.bam
        """
}

process severusTumorOnly {
    def threads = 28

    container 'docker://mkolmogo/severus:dev1.2'
    cpus threads
    memory '128 G'
    time '8.h'

    input:
        path tumorBam
        path tumorBamIdx
        path phasedVcf
        path vntrBed
        path panelOfNormals

    output:
        path 'severus_out/somatic_SVs/severus_somatic.vcf', emit: severusSomaticVcf

    script:
        """
        severus --target-bam ${tumorBam} --out-dir severus_out -t ${threads} --phasing-vcf ${phasedVcf} \
            --vntr-bed ${vntrBed} --PON ${panelOfNormals}
        """
}
