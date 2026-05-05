# ARMS-MBON CO1 Processing Pipeline
This repository contains the scripts and workflow used to process ARMS-MBON COI metabarcoding data from 2022-2025. The processed results are used to generate summary tables, visualizations, taxonomic assignments and supporting material for a data paper.

The workflow is designed to run on the Dardel HPC system using SLURM and an Apptainer sandbox. The main pipeline combines R scripts, command-line bioinformatics tools, and Python-based BOLDigger3 taxonomic assignment. 

## 1. Project overview
The aim of this project is to process a large ARMS-MBON dataset using an existing R-script pipeline. The pipeline starts from paired-end COI FASTQ files and produces curated MOTU tables and taxonomy tables that can be used for summaries, graphs and a final data paper.

The workflow does the following:

Raw paired-end COI FASTQ files

-> Primer removal with cutadapt

-> Quality filtering and ASV inference with DADA2

-> Chimera and singleton removal

-> Removal of mitochondrial pseudogenes / nuMTs using MACSE

-> Negative control / blank correction

-> ASV dereplication and abundance-labelled FASTA creation

-> Clustering ASVs into MOTUs using swarm

-> MOTU table creation

-> LULU curation

-> Taxonomic assignment with BOLDigger3

-> Final count and taxonomy tables

-> Graphs and data paper

## 2. Software and tools
The pipeline uses the following tools:

HPC / container system:
- Dardel HPC
- SLURM
- Apptainer
- PDC module environment

R packages:
- dada2
- Shortread
- Biostrings
- ggplot2
- readxl
- data.table
- tidyverse
- lulu

Command-line tools:
- cutadapt
- MACSE
- swarm
- BLAST+
- BOLDigger3

## Getting repository on Dardel
