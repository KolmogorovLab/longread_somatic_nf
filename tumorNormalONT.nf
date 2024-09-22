#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

include { callClair3; phaseLongphase; severusTumorNormal } from "./processes/processes.nf"
include { alignMinimap2 as alignTumor } from "./processes/processes.nf"
include { alignMinimap2 as alignNormal } from "./processes/processes.nf"
include { haplotagWhatshap as haplotagNormal } from "./processes/processes.nf"
include { haplotagWhatshap as haplotagTumor } from "./processes/processes.nf"

/*
 * Main workflow
 */
workflow tumorNormalOntWorkflow {
    take:
        readsTumor
        readsNormal
        reference
        vntrAnnotation
        clair3Model

    main:
        alignTumor(reference, readsTumor)
        alignNormal(reference, readsNormal)
        callClair3(alignNormal.out.bam, alignNormal.out.bam_idx, reference, alignNormal.out.ref_idx, clair3Model)
        phaseLongphase(alignNormal.out.bam, alignNormal.out.bam_idx, reference, 
                       alignNormal.out.ref_idx, callClair3.out.vcf)
        haplotagNormal(reference, alignNormal.out.ref_idx, phaseLongphase.out.phasedVcf, alignNormal.out.bam, 
                       alignNormal.out.bam_idx)
        haplotagTumor(reference, alignTumor.out.ref_idx, phaseLongphase.out.phasedVcf, alignTumor.out.bam, 
                      alignTumor.out.bam_idx)
        severusTumorNormal(haplotagTumor.out.bam, haplotagTumor.out.bam_idx, 
                           haplotagNormal.out.bam, haplotagNormal.out.bam_idx,
                           phaseLongphase.out.phasedVcf, vntrAnnotation)

    emit:
        phasedVcf = callClair3.out.vcf
        haplotaggedTumor = haplotagTumor.out.bam
        haplotaggedNormal = haplotagNormal.out.bam
        severusSomaticVcf = severusTumorNormal.out.severusSomaticVcf

    publish:
        phasedVcf >> "phased_vcf"
        haplotaggedTumor >> "haplotagged_bam_tumor"
        haplotaggedNormal >> "haplotagged_bam_normal"
        severusSomaticVcf >> "severus"
}

/*
 * Entry point
 */
workflow {
    if (!params.reads_tumor || !params.reads_normal || !params.reference || !params.outdir || 
        !params.vntr || !params.clair3_model) {
        error """
              ERROR: Some required arguments are not defined.
              Usage: tumorNormalONT.nf --reads_tumor PATH --reads_normal PATH --reference PATH --outdir PATH 
                                       --vntr PATH --clair3_model PATH
              """.stripIndent()
    }

    tumorNormalOntWorkflow(Channel.fromPath(params.reads_tumor), Channel.fromPath(params.reads_normal),
                         Channel.fromPath(params.reference), Channel.fromPath(params.vntr),
                         Channel.fromPath(params.clair3_model))
}

output {
    directory params.outdir
    mode "copy"
}
