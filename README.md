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

## 3. Getting repository on Dardel
To get repository on Dardel -> Clone this GitHub

`cd /cfs/klemming/projects/...`

`git clone <GITHUB_REPOSITORY_URL> ARMS-scripts`

Another option is to use Rsync to copy existing directory from home to project storage

## 4. Input data organisation
The first R script expects FASTQ files to be organised by batch under:

`novaseq/COI/fastq_files/`

The data files can be ordered by different batches inside the /fastq-files directory

For example: `noveseq/COI/fastq_files/Batch_1`

To get the FASTQ files in the correct place, use the Rscript `new_symlink.R`

## 5. Metadata
The metadata file used for the current 2022-2025 ARMS-MBON processing is:

`metadata/generated_meta/COI_batch.3.4.5.csv`

This metadata is used later in the pipeline, especially for blank correction.

The metadata also contains sample identifiers that can be used for downstream plotting, for example to extract deployment and retrieval dates from MaterialSampleID.

## 6. Pipeline steps
### Step 1: Primer removal, filtering, trimming and ASV inference
Script: `novaseq/scripts/project_scripts/loessErrfun_mod4_sol.R`

Run command: `Rscript novaseq/scripts/project_scripts/loessErrfun_mod4_sol.R -d /envs/git_env/bin/cutadapt`

This script performs:
- batch detection
- paired FASTQ detection
- primer checking
- primer removal with cutadapt
- quality filtering and trimming with DADA2
- ASV inference
- ASV sequence table creation

Main outputs include:

`novaseq/COI/seqtab_Batch..._mod4.rds`

`novaseq/COI/track_Batch..._mod4.txt`

`novaseq/COI/quality_forward.jpg`

`novaseq/COI/quality_reverse.jpg`

### Step 2: Chimera and singleton removal
Script: `novaseq/scripts/original_scripts/COI_chimera.R`

Run command: `Rscript novaseq/scripts/original_scripts/COI_chimera.R`

This script merges batch sequence tables, filters expected COI sequence lengths, removes chimeras, and removes singleton ASVs

Main outputs include:

`novaseq/COI/COI_nochim_nosingle_ASVs.fa`

`novaseq/COI/COI_nochim_nosingle_ASVs.fa`

### Step 3: nuMT / pseudogene removal
Script: `novaseq/scripts/original_scripts/MACSE_align_pseudo.R`

Run command: `Rscript novaseq/scripts/original_scripts/MACSE_align_pseudo.R -d /envs/git_env/share/macse-2.07-0/macse_v2.07.jar`

This step uses MACSE to detect likely mitochondrial pseudogenes / nuMTs by checking whether ASVs align as valid coding COI sequences.

Main outputs include:

`novaseq/COI/pseudo/pseudo.combined.names.txt`

`novaseq/COI/pseudo/nonpseudo.combined.names.txt`

`novaseq/COI/pseudo/COI_nochim_nosingle_nopseudo.fa`

### Step 4: Negative control / blank correction
Script: `novaseq/scripts/original_scripts/COI_blank_corr.R`

Run command: `Rscript novaseq/scripts/original_scripts/COI_blank_corr.R`

This step uses metadata to identify negative controls and removes ASVs that are likely contaminants.

Main outputs include:

`novaseq/COI/blank_corr/asv_contaminants_COI.txt`

`novaseq/COI/blank_corr/asv_no_contaminants_COI.txt`

`novaseq/COI/blank_corr/no_contam_headers_COI.txt`

`novaseq/COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa`  

### Step 5: Dereplication headers
Script: `novaseq/scripts/original_scripts/dereplication_headers.R`

Run command: `Rscript novaseq/scripts/original_scripts/dereplication_headers.R`

This step creates FASTA headers containing ASV abundance information, which is needed for swarm clustering.

Main outputs include:

`novaseq/COI/ASV_dereplicated.txt`

`novaseq/COI/COI_dereplicated_ASVs.fa`


### Step 6: MOTU clustering with swarm
Command: `/envs/git_env/bin/swarm \
  -d 13 
  -i swarm/internal.txt 
  -o swarm/output.txt 
  -s swarm/statistics.txt 
  -u swarm/uclust.txt 
  -w swarm/COI_cluster_reps.fa 
  COI_dereplicated_ASVs.fa`

This clusters ASVs into MOTUs using swarm with distance parameter -d 13.

Main outputs include:

`novaseq/COI/swarm/output.txt`

`novaseq/COI/swarm/uclust.txt`

`novaseq/COI/swarm/COI_cluster_reps.fa`

### Step 7: MOTU table creation
Script: `novaseq/scripts/original_scripts/MOTU_tables.R`

Run command: `Rscript novaseq/scripts/original_scripts/MOTU_tables.R`

This script converts ASV-level count information and swarm clustering output into a MOTU table.

Main outputs include:

`novaseq/COI/MOTU/motu_table_COI.txt`

### Step 8: BLAST match list for LULU



























