#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

include { alignMinimap2; callClair3; phaseLongphase; deepsomaticTumorOnly;
          wakhanTumorOnly } from "./processes/processes.nf"
          
include {haplotagWhatshap as haplotagWhatshap1} from "./processes/processes.nf"
include {haplotagWhatshap as haplotagWhatshap2} from "./processes/processes.nf"
include {severusTumorOnly as severusTumorOnly1} from "./processes/processes.nf"
include {severusTumorOnly as severusTumorOnly2} from "./processes/processes.nf"

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
        haplotagWhatshap1(reference, alignMinimap2.out.ref_idx, phaseLongphase.out.phasedVcf, alignMinimap2.out.bam, 
                         alignMinimap2.out.bam_idx)
        severusTumorOnly1(haplotagWhatshap1.out.bam, haplotagWhatshap1.out.bam_idx, phaseLongphase.out.phasedVcf, 
                         vntrAnnotation, svPanelOfNormals)
        wakhanTumorOnly(haplotagWhatshap1.out.bam, haplotagWhatshap1.out.bam_idx, reference, phaseLongphase.out.phasedVcf,
                        severusTumorOnly1.out.severusSomaticVcf)
        haplotagWhatshap2(reference, alignMinimap2.out.ref_idx, wakhanTumorOnly.out.rephasedVcf, alignMinimap2.out.bam, 
                         alignMinimap2.out.bam_idx)
        severusTumorOnly2(haplotagWhatshap2.out.bam, haplotagWhatshap2.out.bam_idx, wakhanTumorOnly.out.rephasedVcf, 
                         vntrAnnotation, svPanelOfNormals)
        deepsomaticTumorOnly(alignMinimap2.out.bam, alignMinimap2.out.bam_idx, reference, alignMinimap2.out.ref_idx)

    emit:
        phasedVcf = phaseLongphase.out.phasedVcf
        haplotaggedBam = haplotagWhatshap2.out.bam
        haplotaggedBamidx = haplotagWhatshap2.out.bam_idx
        severusFullOutput = severusTumorOnly2.out.severusFullOutput
        wakhanFullOutput = wakhanTumorOnly.out.wakhanOutput
        deepsomaticOutput = deepsomaticTumorOnly.out.deepsomaticOutput

    publish:
        phasedVcf >> "phased_vcf"
        haplotaggedBam >> "haplotagged_bam"
        haplotaggedBamidx >> "haplotagged_bam"
        severusFullOutput >> "severus"
        wakhanFullOutput >> "wakhan"
        deepsomaticOutput >> "deepsomatic"
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
    readsChannel.view{it -> "Input reads: $it"}

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
