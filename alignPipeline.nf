#!/usr/bin/env nextflow

/*
 * Define the workflow
 */
workflow {
    reads = Channel.fromPath("/data/KolmogorovLab/julieZelaya/LongReadCancerPipeline/HCC1954BL.chr1.bam")
    ref = Channel.fromPath("/data/KolmogorovLab/references/grch38_chr.fasta")

    convertedFastq = convertToFastq(reads)
    alignedSam = alignReads(convertedFastq, ref)
    sortedBam = sortBam(alignedSam)
}

/*
 * Convert BAM to FASTQ using samtools
 */
process convertToFastq {
    label 'alignment'    

    input:
    path reads

    output:
    path 'reads.fastq'

    script:
    """
    samtools fastq -TMm,M1,MM,ML ${reads} > reads.fastq
    """
}

/*
 * Align reads using minimap2
 */
process alignReads {
    label 'alignment'    

    input:
    path fastq
    path ref

    output:
    path 'aligned.sam'

    script:
    """
    minimap2 -ax map-ont -k 17 -t 30 -y --eqx ${ref} ${fastq} > aligned.sam
    """
}

/*
 * Sort BAM using samtools
 */
process sortBam {
    label 'alignment'
    publishDir '/results', mode: 'copy'    

    input:
    path sam

    output:
    path 'aligned.bam'

    script:
    """
    samtools view -Sb ${sam} > aligned.unsorted.bam
    samtools sort -@ 30 -m 4G -o aligned.bam aligned.unsorted.bam
    """
}

