
/*
 * Define the workflow
 */
workflow {
    ref = Channel.fromPath("/data/KolmogorovLab/references/grch38_chr.fasta")
    bam = Channel.fromPath("/data/KolmogorovLab/julieZelaya/LongReadCancerPipeline/results/aligned.bam")

    variantCalling(ref, bam)
}

/*
 * Process to run Clair3
 */
process variantCalling {
    label 'variantCalling'

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

