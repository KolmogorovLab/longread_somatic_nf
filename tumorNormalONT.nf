#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

include { callClair3; phaseLongphase; severusTumorNormal; wakhanTumorNormal; deepsomaticTumorNormal } from "./processes/processes.nf"
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
        alignTumor(reference, readsTumor.collect())
        alignNormal(reference, readsNormal.collect())
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
        wakhanTumorNormal(haplotagTumor.out.bam, haplotagTumor.out.bam_idx, reference, phaseLongphase.out.phasedVcf,
                          severusTumorNormal.out.severusSomaticVcf)
        deepsomaticTumorNormal(alignTumor.out.bam, alignTumor.out.bam_idx, alignNormal.out.bam, alignNormal.out.bam_idx,
                               reference, alignTumor.out.ref_idx)

    emit:
        phasedVcf = phaseLongphase.out.phasedVcf
        haplotaggedTumor = haplotagTumor.out.bam
        haplotaggedTumorIdx = haplotagTumor.out.bam_idx
        haplotaggedNormal = haplotagNormal.out.bam
        haplotaggedNormalIdx = haplotagNormal.out.bam_idx
        severusFullOutput = severusTumorNormal.out.severusFullOutput
        wakhanFullOutput = wakhanTumorNormal.out.wakhanOutput
        deepsomaticOutput = deepsomaticTumorNormal.out.deepsomaticOutput

    publish:
        phasedVcf >> "phased_vcf"
        haplotaggedTumor >> "haplotagged_bam_tumor"
        haplotaggedTumorIdx >> "haplotagged_bam_tumor"
        haplotaggedNormal >> "haplotagged_bam_normal"
        haplotaggedNormalIdx >> "haplotagged_bam_normal"
        severusFullOutput >> "severus"
        wakhanFullOutput >> "wakhan"
        deepsomaticOutput >> "deepsomatic"
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

    tumorChannel = Channel.fromPath(params.reads_tumor.split(" ").toList(), checkIfExists: true)
    tumorChannel.view{it -> "Tumor reads: $it"}

    normalChannel = Channel.fromPath(params.reads_normal.split(" ").toList(), checkIfExists: true)
    normalChannel.view{it -> "Normal reads: $it"}

    tumorNormalOntWorkflow(tumorChannel, normalChannel,
                           Channel.fromPath(params.reference, checkIfExists: true), 
                           Channel.fromPath(params.vntr, checkIfExists: true),
                           Channel.fromPath(params.clair3_model, checkIfExists: true))
}

output {
    directory params.outdir
    mode "copy"
}
