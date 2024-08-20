# Nextflow Workflow
Nextflow processes contain a collection of scripts that perform alignment, small-Variant calling, and haplotagging.
*Current version only runs tumor cells. Input of normal and tumor cells under development*
*Extended pipeline with use of Severus and Wakhan under development*

## Dependencies
- Unix operating system
- Bash 3.2
- Java 11 or later

## Install Nextflow
1. Install Nextflow by copying and pasting the following snippet in your shell terminal: 

    curl -fsSL get.nextflow.io | bash

It will download the 'nextflow' application launcher in your working directory.


2. Make Nextflow executable:

    chmod +x nextflow
   

4. Move Nextflow into an executable path:

    sudo mv nextflow /usr/local/bin
   
  
6. Confirm that Nextflow is installed correctly:

    nextflow info 
   
# Running Nextflow Pipeline

