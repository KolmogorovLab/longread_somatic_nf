#!/usr/bin/env nextflow

params.ref = '/data/KolmogorovLab/references/grch38_chr.fasta'
params.reads = '/data/KolmogorovLab/CellLinesNanoQ20/nano_raw_bams/1945.bam'
params.output_dir = '/data/KolmogorovLab/julieZelaya/clair3_output'
params.threads = 4
params.model_name = 'ont'

process align {
    tag 'alignment'

    input:
    path ref
    path reads

    output:
    path 'alignment.bam'

    script:
    """
    module load minimap2
    module load samtools
    
    samtools fastq -TMm,Ml,MM,ML ${reads} | 
    minimap2 -ax map-ont -k 17 -t 30 -y --eqx ${ref} - | 
    samtools sort -@4 -m 4G -o alignment.bam
    """
}

process run_clair3 {
    tag 'clair3'

    input:
    path ref
    path bam

    output:
    path 'clair3_output'

    script:
    """
    module load singularity

    singularity pull docker://kishwars/pepper_deepvariant:latest
    
    mkdir -p clair3_output

    singularity exec docker://kishwars/pepper_deepvariant:latest \
    /opt/bin/run_clair3.sh \
    --bam_fn=${bam} \
    --ref_fn=${ref} \
    --threads=${params.threads} \
    --platform="ont" \
    --model_path="/opt/models/${params.model_name}" \
    --output=clair3_output \
    --enable_phasing \
    --longphase_for_phasing
    """
}

process haplotag {
    tag 'haplotag'

    input:
    path ref
    path phased_vcf
    path normal_bam
    path tumor_bam

    output:
    path 'normal.haplotagged.bam'
    path 'tumor.haplotagged.bam'

    script:
    """
    conda create -n whatshap-env whatshap
    source activate whatshap-env

    whatshap haplotag --reference ${ref} ${phased_vcf} ${normal_bam} -o normal.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=${params.threads}
    
    whatshap haplotag --reference ${ref} ${phased_vcf} ${tumor_bam} -o tumor.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=${params.threads}
    """
}

workflow {
    ref = file(params.ref)
    reads = file(params.reads)

    alignment = align(ref: ref, reads: reads)

    clair3_output = run_clair3(ref: ref, bam: alignment)

    haplotagged_bams = haplotag(ref: ref, phased_vcf: clair3_output, normal_bam: alignment, tumor_bam: alignment)
}

