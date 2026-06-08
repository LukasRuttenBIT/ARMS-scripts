#!/usr/bin/env Rscript

# =============================================================================
# MACSE_align_pseudo.R — nuMT removal via pseudo-alignment
#
# Tries to align each ASV against reference COI sequences across five
# genetic codes (gc5 → gc4 → gc2 → gc9 → gc13). ASVs that fail alignment
# under all five codes are classified as nuMTs and excluded.
#
# Usage:
#   Rscript run/pipeline/MACSE_align_pseudo.R \
#       --datadir novaseq \
#       -d /path/to/macse_v2.07.jar
#
# Input:  {datadir}/COI/COI_nochim_nosingle_ASVs.fa
# Output: {datadir}/COI/pseudo/nonpseudo.combined.names.txt
#         {datadir}/COI/pseudo/pseudo.combined.names.txt
#         {datadir}/COI/pseudo/gc{2,4,5,9,13}/  (intermediate MACSE files)
# =============================================================================

library(Biostrings)
library(hiReadsProcessor)
library(dplyr)
library(seqinr)
library(parallel)
library(argparse)

parser <- ArgumentParser(description = "nuMT removal using MACSE")
parser$add_argument("-d", "--directory", required = TRUE,
                    help = "Path to macse_v2.07.jar")
parser$add_argument("--datadir", default = "novaseq",
                    help = "Data directory containing COI/ subdir (default: novaseq)")
parser$add_argument("-c", "--cores", type = "integer", default = 20,
                    help = "Cores for parallel MACSE split processing (default: 20)")
args <- parser$parse_args()

macse_jar <- args$directory
datadir   <- args$datadir
n_cores   <- args$cores

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

path      <- file.path(datadir, "COI", "pseudo")
path_gc2  <- file.path(path, "gc2")
path_gc4  <- file.path(path, "gc4")
path_gc5  <- file.path(path, "gc5")
path_gc9  <- file.path(path, "gc9")
path_gc13 <- file.path(path, "gc13")

for (p in c(path, path_gc2, path_gc4, path_gc5, path_gc9, path_gc13))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)

# Reference alignments (from Daraghmeh 2024, in metadata/)
macse_gc2  <- "metadata/macse_taxonomic_alignment/MACSE_BOLD_gc2_marine_taxa_aligned_NT.fasta"
macse_gc4  <- "metadata/macse_taxonomic_alignment/MACSE_BOLD_gc4_marine_taxa_aligned_NT.fasta"
macse_gc5  <- "metadata/macse_taxonomic_alignment/MACSE_BOLD_gc5_marine_taxa_aligned_NT.fasta"
macse_gc9  <- "metadata/macse_taxonomic_alignment/MACSE_BOLD_gc9_marine_taxa_aligned_NT.fasta"
macse_gc13 <- "metadata/macse_taxonomic_alignment/MACSE_BOLD_Tunicata_gc13_aligned_NT.fasta"

nochim_fa <- file.path(datadir, "COI", "COI_nochim_nosingle_ASVs.fa")
if (!file.exists(nochim_fa)) stop("Input fasta not found: ", nochim_fa)

fastafile <- read.fasta(nochim_fa, seqtype = "DNA", as.string = TRUE)
message("Loaded ", length(fastafile), " ASVs from ", nochim_fa)

# ---------------------------------------------------------------------------
# Helper: run MACSE enrichAlignment on split files, return non-pseudo names
# ---------------------------------------------------------------------------

run_macse_gc <- function(input_fa, gc_code, gc_ref, path_split, n_splits, cores = n_splits) {
  if (!dir.exists(path_split)) dir.create(path_split, recursive = TRUE)
  x <- readDNAStringSet(input_fa)
  if (length(x) == 0) return(list(nonpseudo = data.frame(), pseudo = data.frame()))

  splitSeqsToFiles(x, n_splits, "fasta", "splitseqs", path_split)
  split_files <- sort(list.files(path_split, pattern = "\\.fasta$", full.names = TRUE))

  get_name  <- function(f) strsplit(basename(f), "\\.fasta")[[1]][1]
  snames    <- unname(sapply(split_files, get_name))
  outputAA  <- file.path(path_split, paste0(snames, "_AA_gc", gc_code, ".fa"))
  outputNT  <- file.path(path_split, paste0(snames, "_NT_gc", gc_code, ".fa"))
  outputST  <- file.path(path_split, paste0(snames, "_stats_gc", gc_code, ".csv"))

  message("  Running ", length(split_files), " splits in parallel (cores=",
          min(length(split_files), cores), ")...")
  mclapply(seq_along(split_files), function(i) {
    system2("java", args = c(
      "-jar", macse_jar,
      "-prog", "enrichAlignment",
      "-align", gc_ref,
      "-seq",   split_files[i],
      paste0("-gc_def ", gc_code,
             " -maxSTOP_inSeq 0 -output_only_added_seq_ON =TRUE",
             " -fixed_alignment_ON =TRUE -maxDEL_inSeq 3",
             " -maxFS_inSeq 0 -maxINS_inSeq 0"),
      " -out_AA ", outputAA[i],
      " -out_NT ", outputNT[i],
      " -out_tested_seq_info ", outputST[i]
    ), stdout = FALSE, stderr = FALSE)
  }, mc.cores = min(length(split_files), cores))

  stat_files <- sort(list.files(path_split,
                     pattern = paste0("_stats_gc", gc_code, "\\.csv$"),
                     full.names = TRUE))

  nonpseudo <- data.frame()
  pseudo    <- data.frame()
  for (f in stat_files) {
    res <- read.table(f, h = TRUE, sep = ";")
    nonpseudo <- rbind(nonpseudo, res[res$added == "yes", ])
    pseudo    <- rbind(pseudo,    res[res$added == "no",  ])
  }
  list(nonpseudo = nonpseudo, pseudo = pseudo)
}

# ---------------------------------------------------------------------------
# gc5 — first pass on full fasta
# ---------------------------------------------------------------------------

message("\n--- Genetic code 5 ---")
path_gc5_split <- file.path(path_gc5, "splitseqs")
res5 <- run_macse_gc(nochim_fa, 5, macse_gc5, path_gc5_split, 20, cores = n_cores)

# Remaining putative nuMTs passed to gc4
pseudo5_names <- as.character(res5$pseudo$name)
if (length(pseudo5_names) > 0) {
  pseudo5_fa <- file.path(path_gc5_split, "gc5_pseudo.fasta")
  write.fasta(fastafile[names(fastafile) %in% pseudo5_names],
              names = pseudo5_names, file.out = pseudo5_fa)
} else {
  pseudo5_fa <- NULL
}

# ---------------------------------------------------------------------------
# gc4
# ---------------------------------------------------------------------------

message("\n--- Genetic code 4 ---")
nonpseudo4 <- pseudo4 <- data.frame()
if (!is.null(pseudo5_fa) && file.exists(pseudo5_fa)) {
  path_gc4_split <- file.path(path_gc4, "splitseqs")
  res4 <- run_macse_gc(pseudo5_fa, 4, macse_gc4, path_gc4_split, 20, cores = n_cores)
  nonpseudo4 <- res4$nonpseudo
  pseudo4    <- res4$pseudo

  pseudo4_names <- as.character(pseudo4$name)
  if (length(pseudo4_names) > 0) {
    pseudo4_fa <- file.path(path_gc4_split, "gc4_pseudo.fasta")
    write.fasta(fastafile[names(fastafile) %in% pseudo4_names],
                names = pseudo4_names, file.out = pseudo4_fa)
  } else { pseudo4_fa <- NULL }
} else { pseudo4_fa <- NULL }

# ---------------------------------------------------------------------------
# gc2
# ---------------------------------------------------------------------------

message("\n--- Genetic code 2 ---")
nonpseudo2 <- pseudo2 <- data.frame()
if (!is.null(pseudo4_fa) && file.exists(pseudo4_fa)) {
  path_gc2_split <- file.path(path_gc2, "splitseqs")
  res2 <- run_macse_gc(pseudo4_fa, 2, macse_gc2, path_gc2_split, 20, cores = n_cores)
  nonpseudo2 <- res2$nonpseudo
  pseudo2    <- res2$pseudo

  pseudo2_names <- as.character(pseudo2$name)
  if (length(pseudo2_names) > 0) {
    pseudo2_fa <- file.path(path_gc2_split, "gc2_pseudo.fasta")
    write.fasta(fastafile[names(fastafile) %in% pseudo2_names],
                names = pseudo2_names, file.out = pseudo2_fa)
  } else { pseudo2_fa <- NULL }
} else { pseudo2_fa <- NULL }

# ---------------------------------------------------------------------------
# gc9
# ---------------------------------------------------------------------------

message("\n--- Genetic code 9 ---")
nonpseudo9 <- pseudo9 <- data.frame()
if (!is.null(pseudo2_fa) && file.exists(pseudo2_fa)) {
  path_gc9_split <- file.path(path_gc9, "splitseqs")
  res9 <- run_macse_gc(pseudo2_fa, 9, macse_gc9, path_gc9_split, 20, cores = n_cores)
  nonpseudo9 <- res9$nonpseudo
  pseudo9    <- res9$pseudo

  pseudo9_names <- as.character(pseudo9$name)
  if (length(pseudo9_names) > 0) {
    pseudo9_fa <- file.path(path_gc9_split, "gc9_pseudo.fasta")
    write.fasta(fastafile[names(fastafile) %in% pseudo9_names],
                names = pseudo9_names, file.out = pseudo9_fa)
  } else { pseudo9_fa <- NULL }
} else { pseudo9_fa <- NULL }

# ---------------------------------------------------------------------------
# gc13
# ---------------------------------------------------------------------------

message("\n--- Genetic code 13 ---")
nonpseudo13 <- pseudo13 <- data.frame()
if (!is.null(pseudo9_fa) && file.exists(pseudo9_fa)) {
  path_gc13_split <- file.path(path_gc13, "splitseqs")
  res13 <- run_macse_gc(pseudo9_fa, 13, macse_gc13, path_gc13_split, 20, cores = n_cores)
  nonpseudo13 <- res13$nonpseudo
  pseudo13    <- res13$pseudo
}

# ---------------------------------------------------------------------------
# Combine results across all genetic codes
# ---------------------------------------------------------------------------

pseudo.combined    <- rbind(res5$pseudo, pseudo4, pseudo2, pseudo9, pseudo13)
nonpseudo.combined <- rbind(res5$nonpseudo, nonpseudo4, nonpseudo2, nonpseudo9, nonpseudo13)

pseudo.names    <- paste0(">", as.character(unique(pseudo.combined$name)))
nonpseudo.names <- paste0(">", as.character(unique(nonpseudo.combined$name)))

write.table(pseudo.names,    file.path(path, "pseudo.combined.names.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(nonpseudo.names, file.path(path, "nonpseudo.combined.names.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

message("\n===== MACSE_align_pseudo.R complete =====")
message("  nuMTs removed  : ", length(unique(pseudo.combined$name)))
message("  ASVs retained  : ", length(unique(nonpseudo.combined$name)))
message("  nonpseudo list : ", file.path(path, "nonpseudo.combined.names.txt"))
