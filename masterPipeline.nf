#!/usr/bin/env nextflow

/*
 * Define the workflow
 */
workflow {
    // Create channels for reads and reference
    reads = Channel.fromPath(params.reads)
    ref = Channel.fromPath(params.ref)

    // Alignment
    (alignedBam, indexedBai, indexedRef) = align(ref, reads)

    // Pass the BAM file, BAI file, and reference index to the variantCalling process
    (vcf, vcfIndex) = variantCalling(ref, alignedBam, indexedBai, indexedRef)

    // Haplotagging
    haplotaggedVcf = phasing(ref, indexedRef, vcf, alignedBam, indexedBai, vcfIndex)

    // Run Severus on the phased data

    runSeverus(haplotaggedVcf, alignedBam, indexedBai)
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
    path 'aligned.bam'
    path 'aligned.bam.bai'
    path "${ref}.fai"
          
    script:
    """  
    samtools fastq -TMm,M1,MM,ML ${reads} |
    minimap2 -ax map-ont -k 17 -t 30 -y --eqx ${ref} - |
    samtools sort -@30 -m 4G -o aligned.bam -
    samtools index aligned.bam
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
    path indexedBai
    path indexedRef	

    output:
    path 'clair3_output/pileup.vcf.gz'
    path 'clair3_output/pileup.vcf.gz.tbi'

    script:
    """
    run_clair3.sh \
    ${baseDir}/clair3_latest.sif \
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
process phasing {
    label 'phasing'

    input:
    path ref
    path indexedRef
    path vcf
    path alignedBam
    path indexedBai
    path vcfIndex

    output:
    path 'haplotagged.vcf'

    script:
    """
    whatshap haplotag -o haplotagged.vcf --reference ${ref} --ignore-read-groups 'pileup.vcf.gz' 'aligned.bam'
    """
}

/*

* Process to run Severus

*/

process runSeverus {
    label 'severus'
 
    input:
    path haplotaggedVcf
    path alignedBam
    path indexedBai

    output:
    path 'severus_out/*'

    script:
    """
    # Install Severus if not already installed
    if ! conda list -n severus_env | grep -q "severus"; then
    conda create -n severus_env 

    # Activate the environment and run Severus
    conda activate severus_env

    # Single Sample SV Calling
    severus --target-bam ${alignedBam} --out-dir severus_out -t 16 --phasing-vcf ${haplotaggedVcf} \
    --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
    """

}
