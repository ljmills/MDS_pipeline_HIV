# MDS_pipeline_HIV
User-friendly pipeline for mutational analysis of HIV using ultra-accurate maximum-depth sequencing


# User-friendly Pipeline for Analysis of Ultra-Accurate Maximum-Depth Sequencing 
This repository contains the pipline described in Development of a user-friendly pipeline for mutational anal-yses of HIV using ultra-accurate maximum-depth sequencing 

There are two versions of the alignment and mutation calling pipline, the first is a BASH script designed to be used with a Linux computing cluster with job submission. The second is a Galaxy pipeline that can be uploaded to any Galaxy instance and once the tools are installed used.

The included R scripts take the output from bcftools mpileup from either the Galaxy or BASH scripts and performs the hotspot identification and generates the data for the mutational profile analysis seen in the publication. 

# Reference Genome FASTA 
You will need a FASTA formatted version of the reference genome that you are working with. This could be the sequence from a specific plasmid or a reference genome downloaded from a repository such as ENSEMBL or the UCSC genome browser. In the case of this paper we also hard masked (replaced some bp with Ns) across highly similar regions so we didn't get cross-mapping to regions we were not targeting. We used bedtools maskfasta and a bed file containing the regions we wanted masked i.e. 3' UTR region. 

# Galaxy Workflow 

# BASH scripts (least user-friendly) 
## Software Dependencies 
- *[Java 8](https://java.com/en/download/)* (aka Java 1.8) or later 
- *[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)* verion 0.33 or later 
- *[fgbio](https://github.com/fulcrumgenomics/fgbio)*
- samtools 1.9 or later
- bwa
- picard tools
- bcftools 1.9 or later
- bedtools 2.29 or later

## Script overview
All scripts will need to be edited to be used with your compute enviroment. Edits could include paths to software, paths to specific files needed and edits to how we extract sample names from the FASTQ files. They also include SLURM submission information ( sbatch ) which will change based on if you are running the SLURM scheduler or not. 

## indexNewGenome.sh
Helper script to create bwa, picard tools and samtools indices from FASTA reference genomes. 

## by_Sample
