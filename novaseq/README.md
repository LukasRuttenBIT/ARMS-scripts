

## Steps of pipeline:

### Filter and trim using cutadapt 

`Rscript novaseq/scripts/original_scripts/COI_trimming_filtering.R -d CUTADAPT_DIR`

### Remove chimeras and singletons

`Rscript scripts/original_scripts/COI_chimera.R`

### Remove nuMTs using MACSE

`Rscript scripts/original_scripts/MACSE_align_pseudo.R -d PATH/TO/MACSE.jar`

(macOS) `grep -w -A 1 -f novaseq/COI/pseudo/nonpseudo.combined.names.txt novaseq/COI/COI_nochim_nosingle_ASVs.fa > novaseq/COI/pseudo/COI_nochim_nosingle_nopseudo.fa`

(Windows/GNU) `grep -w -A 1 -f COI/pseudo/nonpseudo.combined.names.txt COI/COI_nochim_nosingle_ASVs.fa > COI/pseudo/COI_nochim_nosingle_nopseudo.fa  --no-group-separator`

### Blank correction 

`Rscript novaseq/scripts/original_scripts/COI_blank_corr.R`

`grep -w -A 1 -f novaseq/COI/blank_corr/no_contam_headers_COI.txt novaseq/COI/pseudo/COI_nochim_nosingle_nopseudo.fa > novaseq/COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa`

### Clustering ASVs into MOTUs

`Rscript novaseq/scripts/original_scripts/dereplication_headers.R`

`grep -v "^--" novaseq/COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa | awk 'NR%2==0' | paste -d'\n' novaseq/COI/ASV_dereplicated.txt - > novaseq/COI/COI_dereplicated_ASVs.fa`

`cd novaseq/COI/`

`PATH/TO/SWARM -d 13 -i swarm/internal.txt -o swarm/output.txt -s swarm/statistics.txt -u swarm/uclust.txt -w swarm/COI_cluster_reps.fa COI_dereplicated_ASVs.fa`

`/opt/miniconda3/envs/main_env/bin/swarm -d 13 -i swarm/internal.txt -o swarm/output.txt -s swarm/statistics.txt -u swarm/uclust.txt -w swarm/COI_cluster_reps.fa COI_dereplicated_ASVs.fa`

`cd ../..`

### LULU curation

`Rscript novaseq/scripts/original_scripts/MOTU_tables.R`

`awk -F'_' '{print $1}' novaseq/COI/swarm/COI_cluster_reps.fa > novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa`

`makeblastdb -in novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa -parse_seqids -dbtype nucl`

`blastn -db novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa -outfmt '6 qseqid sseqid pident' -out novaseq/COI/MOTU/match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa`

`Rscript novaseq/scripts/original_scripts/LULU_curation.R`

(macOS) `grep -w -A 1 -f novaseq/COI/MOTU/lulu_curated_headers.txt novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa | grep -v "^--" > novaseq/COI/MOTU/COI_cluster_reps_lulu_curated.fa`

(Windows/GNU) `grep -w -A 1 -f COI/MOTU/lulu_curated_headers.txt COI/MOTU/COI_cluster_reps_lulu_ready.fa --no-group-separator > COI/MOTU/COI_cluster_reps_lulu_curated.fa`

### Taxonomic assignment using BOLDigger3

`pip install boldigger3`

`pip install lxml_html_clean`

`boldigger3 identify novaseq/COI/MOTU/COI_cluster_reps_lulu_curated.fa --db 1 --mode 3`

`Rscript scripts/project_scripts/BOLDigger_tax_table.R`
