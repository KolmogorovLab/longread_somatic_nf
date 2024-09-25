# Nextflow Workflow: Long-Read Somatic Variant Calling
Nextflow processes contain a collection of scripts that perform alignment, small-Variant calling, and haplotagging.

Currently pipeline includes the following:
* Alignment with minimap2
* Small variant calling with Clair3
* Phasing with longphase
* Somatic SV calling with Severus
* CNA calling with Wakhan
  
## Running on Biowulf/Slurm

Note that the commands below need to be run either in interactive session, or using `sbatch`, but not on the head node. The nextflow engine automatically submits jobs via Slurm scheduler.

### Tumor-normal version

```
module load nextflow
nextflow run tumorNormalONT.nf --reads_tumor TUMOR_BAM --reads_normal NORMAL_BAM --reference REF_FASTA --outdir NF_OUT --vntr VNTR_BED --clair3_model CLAIR3_MODEL
```

### Tumor-only version

```
module load nextflow
nextflow run ~/projects/LongReadCancerPipeline/tumorOnlyONT.nf --reads READS_BAM --reference REF_FASTA --outdir NF_OUT --vntr VNTR_BED --sv_pon SEVERUS_PON --clair3_model CLAIR3_MODEL
```

### Some caveats

* Input reads are in (unmapped) bam format. You can specify multiple input files, but the whole argument needs to be wrapped in *single quotes* (`'`): `--reads_tumor 'BAM1 BAM2'`, or `reads_tumor 'BAM_DIR/*bam'`
* If you are running multiple instances of nextflow (e.g. processing multiple samples), each instance needs to be run from *separate* working directory. This is because they create `work` folder for intermediate files and `.nextflow` logs.
* If pipeline failed for some reason it can be restarted by adding `-resume` argument (note single dash). You can modify the pipeline and restart, and it will attemnt to reuse the results if possible. Restart should happen inside the same working directory.

## Resources to learn about Nextflow

* Highly suggest to go over this introduction if this is your first time with Nextflow: https://training.nextflow.io/

## Tips for developing/debugging

* To run the pipeline locally, comment out `process.executor = 'slurm'` in the `nextflow.config` file
* Nextflow will report running/completed jobs as follows `[11/4852c3] tumorOnlyOntWorkflow:alignMinimap2 (1) [  0%] 0 of 1`. You can find the associated files for this job at the (unique) working directory: `work/11/4852c3...`.
* There, you'll see staged input files and everything that job output. There are also hidden files (start with `.`, use `ls -la` to show).
* In the working dir: `.command.sh` shows the executed command line (run inside of singularity container). `.command.out` - stdout, `.command.err` - stderr.
* Each process is run inside a singularity/docker container (referenced in the process definition). Build scripts for some of the custom docker containers are in `docker` folder.
