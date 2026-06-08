#!/usr/bin/env Rscript

library(dada2)

# specify input directory 
input_dir <- file.path("novaseq", "COI")

# specify output directory 
output_dir <- file.path("novaseq", "COI")

# load sequencing batch numbers
batch_dir <- file.path("novaseq", "COI", "fastq_files")

# Only use real batch folders, e.g. Batch_1, Batch_2, ..., Batch_11
batch_list <- list.files(batch_dir, pattern = "^Batch_[0-9]+$")

# Sort batches numerically
batch_list <- batch_list[order(as.numeric(gsub("^Batch_", "", batch_list)))]

# Extract batch numbers
coi_run <- gsub("^Batch_", "", batch_list)

message("Batches used for chimera removal: ", paste(batch_list, collapse = ", "))
message("Batch numbers: ", paste(coi_run, collapse = ", "))

# Load the sequence tables of the different sequence runs

rds_list <- list()
for (run in coi_run) {
  run_file <- file.path(input_dir, paste("seqtab_Batch", run, "_mod4.rds", sep = ""))
  rds_run <- readRDS(run_file)
  rds_list[[as.character(run)]] <- rds_run
}

# Merge sequence tables

if (length(rds_list) > 1) {
  merged <- do.call(
    mergeSequenceTables,
    c(unname(rds_list), list(repeats = "sum"))
  )

  # Keep sequence reads with a length of 310, 313 or 316 bp only.
  seqtab.filtered <- merged[, nchar(colnames(merged)) %in% c(310, 313, 316)]

} else { 
  merged <- as.matrix(rds_list[[1]])

  # Keep sequence reads with a length of 310, 313 or 316 bp only.
  seqtab.filtered <- merged[, nchar(colnames(merged)) %in% c(310, 313, 316)]
}

saveRDS(merged, file = file.path(output_dir, "merged_seqtab_coi.rds"))

saveRDS(seqtab.filtered, file = file.path(output_dir, "seqtab_filtered_coi.rds"))

# Remove chimeras #

seqtab.nochim <- removeBimeraDenovo(seqtab.filtered, multithread=T, verbose=TRUE)

# Save sequence table with the non-chimeric sequences as RDS file:

saveRDS(seqtab.nochim, file = file.path(output_dir, "seqtab_nochim_coi.rds"))

# It is possible that a large fraction of the total number of UNIQUE SEQUENCES will be chimeras.
# However, this is usually not the case for the majority of the READS.
# Calculate percentage of the reads that were non-chimeric.

non_chimeric <- (sum(seqtab.nochim)/sum(merged))*100
message(non_chimeric, "% of sequences were non-chimeric")

# Remove singletons from the non-chimeric ASVs

#Transform counts to numeric (as they will most likely be integers)
mode(seqtab.nochim) = "numeric"

# Subset columns with counts of > 1 and save to file
seqtab.nochim.nosingle <- seqtab.nochim[,colSums(seqtab.nochim) > 1]
saveRDS(seqtab.nochim.nosingle, file = file.path(output_dir, "seqtab_nochim_nosingle_coi.rds"))

# Write a fasta file of the final, non-chimeric , non-singleton sequences with short >ASV... type headers

asv_seqs <- colnames(seqtab.nochim.nosingle)
asv_headers <- vector(dim(seqtab.nochim.nosingle)[2], mode = "character")
for (i in 1:dim(seqtab.nochim.nosingle)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = file.path(output_dir, "COI_nochim_nosingle_ASVs.fa"))

# Write an ASV count table of the final, non-chimeric, non-singleton sequences with short >ASV... type names

colnames(seqtab.nochim.nosingle) <- paste0("ASV", seq(ncol(seqtab.nochim.nosingle)))

ASV_counts<-t(seqtab.nochim.nosingle) # transposing table

write.table(ASV_counts, file = file.path(output_dir, "COI_ASV_counts_nosingle.txt"), sep = "\t", quote = F, col.names = NA)

## Track reads through the entire dada2 pipeline ##

# Track reads through the chimera and singleton removal step.

# Read non-chimeric and non-singleton table again, as it has been modified
seqtab.nochim.nosingle <- readRDS(file = file.path(output_dir, "seqtab_nochim_nosingle_coi.rds"))

track_nochim_nosingle <- cbind(rowSums(seqtab.filtered), rowSums(seqtab.nochim), rowSums(seqtab.nochim.nosingle))
colnames(track_nochim_nosingle) <- c("length_filt", "nonchim", "nosingle")

# Read tracking tables of the single runs and combine

track_list <- list()
for (run in coi_run) {
  track_file <- file.path(input_dir, paste("track_Batch", run, "_mod4.txt", sep = ""))
  txt_track <- read.table(track_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  txt_track$sample_id <- rownames(txt_track)
  track_list[[as.character(run)]] <- txt_track
}

tracks <- do.call(rbind, unname(track_list))

# Sum duplicate samples across batches, matching repeats = "sum" above
tracks <- aggregate(. ~ sample_id, data = tracks, FUN = sum)

rownames(tracks) <- tracks$sample_id
tracks$sample_id <- NULL

# Combine all tracking tables

tracks <- tracks[order(match(rownames(tracks), rownames(track_nochim_nosingle))),]
track_coi <- cbind(tracks, track_nochim_nosingle)

# Calculate percentages for each step compared to input

track_coi$cutadapt_perc <- (track_coi$cutadapt / track_coi$input)*100
track_coi$filtered_perc <- (track_coi$filtered / track_coi$input)*100
track_coi$denoisedF_perc <- (track_coi$denoisedF / track_coi$input)*100
track_coi$denoisedR_perc <- (track_coi$denoisedR / track_coi$input)*100
track_coi$merged_perc <- (track_coi$merged / track_coi$input)*100
track_coi$length_filt_perc <- (track_coi$length_filt / track_coi$input)*100
track_coi$nonchim_perc <- (track_coi$nonchim / track_coi$input)*100
track_coi$nosingle_perc <- (track_coi$nosingle / track_coi$input)*100

# Save final read tracking table to file

write.table(track_coi, file = file.path(output_dir, "track_COI.txt"), sep = "\t", col.names = NA)
