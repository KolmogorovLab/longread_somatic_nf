# Nextflow Workflow: Long-Read Sequencing Analysis
Nextflow processes contain a collection of scripts that perform alignment, small-Variant calling, and haplotagging.
*Current version only runs tumor cells. Input of normal and tumor cells under development*
*Extended pipeline with use of Severus and Wakhan under development*


## Dependencies
- Unix operating system
- Bash 3.2
- Java 11 or later


## Install Nextflow
1. Install Nextflow by copying and pasting the following snippet in your shell terminal: 

```
curl -fsSL get.nextflow.io | bash
```

It will download the 'nextflow' application launcher in your working directory.


2. Make Nextflow executable:

```
chmod +x nextflow
```  

4. Move Nextflow into an executable path:

```
sudo mv nextflow /usr/local/bin
```  
  
6. Confirm that Nextflow is installed correctly:

```
nextflow info 
```

## Inputs and Parameters
The current version of the pipeline takes tumor-only long-read cell lines in bam format. The `nextflow.config` file specifies the path of the `reads` as well as the path of the reference fasta genome, labeled as `ref`. These two parameters are **required** and must be updated with the proper paths of the input long-reads and reference file. The current configuration file shows the following:

```
params {
    reads = "/data/KolmogorovLab/julieZelaya/LongReadCancerPipeline/HCC1954.chr22.bam"
    ref = "/data/KolmogorovLab/references/grch38_chr.fasta"
}
```

## Running Nextflow Pipeline
The following files located in this repository are **required** to run the Long-Read Analysis Pipeline: 
1. Either method below to run the pipeline
    - `masterPipeline.nf` to run the pipeline directly through nextflow
    - `slurmPipeline.sh` to run the nextflow pipeline through slurm
2. Configuration file with established parameters, environment and container configurations
   - `nextflow.config` 


### Citations
P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:10.1038/nbt.3820

