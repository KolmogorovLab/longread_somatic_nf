#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

include { alignMinimap2; callClair3; phaseLongphase; 
          haplotagWhatshap; severusTumorOnly } from "./processes/processes.nf"

/*
 * Main workflow
 */
workflow tumorOnlyOntWorkflow {
    take:
        reads
        reference
        vntrAnnotation
        svPanelOfNormals
        clair3Model

    main:
        alignMinimap2(reference, reads.collect())
        callClair3(alignMinimap2.out.bam, alignMinimap2.out.bam_idx, reference, alignMinimap2.out.ref_idx, clair3Model)
        phaseLongphase(alignMinimap2.out.bam, alignMinimap2.out.bam_idx, reference, 
                       alignMinimap2.out.ref_idx, callClair3.out.vcf)
        haplotagWhatshap(reference, alignMinimap2.out.ref_idx, phaseLongphase.out.phasedVcf, alignMinimap2.out.bam, 
                         alignMinimap2.out.bam_idx)
        severusTumorOnly(haplotagWhatshap.out.bam, haplotagWhatshap.out.bam_idx, phaseLongphase.out.phasedVcf, 
                         vntrAnnotation, svPanelOfNormals)

    emit:
        phasedVcf = callClair3.out.vcf
        haplotaggedBam = haplotagWhatshap.out.bam
        severusFullOutput = severusTumorOnly.out.severusFullOutput

    publish:
        phasedVcf >> "phased_vcf"
        haplotaggedBam >> "haplotagged_bam"
        severusFullOutput >> "severus"
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

    readsChannel = Channel.fromPath(params.reads.split(" ").toList(), checkIfExists: true)
    readsChannel.collect().view{it -> "Input read files: $it"}

    tumorOnlyOntWorkflow(readsChannel, 
                         Channel.fromPath(params.reference, checkIfExists: true), 
                         Channel.fromPath(params.vntr, checkIfExists: true), 
                         Channel.fromPath(params.sv_pon, checkIfExists: true),
                         Channel.fromPath(params.clair3_model, checkIfExists: true))
}

output {
    directory params.outdir
    mode "copy"
}
