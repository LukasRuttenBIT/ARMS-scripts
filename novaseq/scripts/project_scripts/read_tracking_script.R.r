# =============================================================================
# Phase 1 QC — Read tracking visualization
# Requires: ggplot2, dplyr, tidyr, scales, ggrepel
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)

setwd("/cfs/klemming/projects/supr/naiss2025-23-46/Lukas/ARMS-scripts")

# ---- SETTINGS ----------------------------------------------------------------
MISEQ_DIR   <- "miseq"    # folder with track_Run*.txt files
NOVASEQ_DIR <- "novaseq"  # folder with track_Batch*.txt files
OUTPUT_PDF  <- "phase1_qc.pdf"
MERGE_WARN  <- 20    # flag samples below this merge %
DENOISE_WARN <- 30   # flag samples below this denoised %

# ---- LOAD DATA ---------------------------------------------------------------
read_tracking <- function(dir, instrument) {
  files <- list.files(dir, pattern = "^track_.*\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  lapply(files, function(f) {
    label <- sub("^track_(.+)\\.txt$", "\\1", basename(f))
    df    <- read.table(f, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, check.names = FALSE)
    df$batch      <- label
    df$instrument <- instrument
    df
  }) |> bind_rows()
}

track <- bind_rows(
  read_tracking(NOVASEQ_DIR, "NovaSeq")
)

# Sort batches sensibly (Run1 < Run2 … < Batch1 < Batch2 …)
batch_order <- c(
  paste0("Batch", sort(as.integer(sub("Batch", "", grep("^Batch", unique(track$batch), value = TRUE)))))
)


track <- track |>
  filter(!is.na(input), input > 0) |>
  mutate(
    sample_type = ifelse(grepl("^CEB_", sample), "control", "real"),
    merge_pct   = merged   / input * 100,
    denoise_pct = denoised / input * 100,
    batch       = factor(batch, levels = intersect(batch_order, unique(batch)))
  )

real <- filter(track, sample_type == "real")

# Colour palette — one colour per batch
n_batches  <- nlevels(track$batch)
batch_cols <- setNames(
  colorRampPalette(c("#2166ac","#4dac26","#d01c8b","#f1a340","#018571","#7b3294"))(n_batches),
  levels(track$batch)
)

base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        strip.background  = element_rect(fill = "grey92"),
        legend.position   = "right")

# =============================================================================
# PLOT 1 — Total reads processed per batch
# =============================================================================
totals <- real |>
  group_by(instrument, batch) |>
  summarise(across(c(input, cutadapt, merged, chimera, denoised), 
                 \(x) sum(x, na.rm = TRUE)),
          n_samples = n(), .groups = "drop") |>
  pivot_longer(cols = c(input, cutadapt, merged, chimera, denoised),
               names_to = "step", values_to = "reads") |>
  mutate(step = factor(step,
                       levels = c("input","cutadapt","merged","chimera","denoised"),
                       labels = c("Input","Cutadapt","Merged","Chimera-free","Denoised")))

p1 <- ggplot(totals, aes(x = step, y = reads / n_samples / 1e6, fill = batch)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~instrument, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = batch_cols) +
  scale_y_continuous(labels = label_comma(suffix = " M")) +
  labs(title  = "Mean reads per sample at each processing step",
       subtitle = "Real samples only",
       x = NULL, y = "Million reads (mean per sample)", fill = "Batch / Run") +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# =============================================================================
# PLOT 2 — Mean retention profile (% of input) per batch
# =============================================================================
retention <- real |>
  filter(!is.na(merged)) |>
  group_by(instrument, batch) |>
  summarise(
    Cutadapt    = mean(cutadapt  / input * 100, na.rm = TRUE),
    Merged      = mean(merged    / input * 100, na.rm = TRUE),
    `Chim-free` = mean(chimera   / input * 100, na.rm = TRUE),
    Denoised    = mean(denoised  / input * 100, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(Cutadapt, Merged, `Chim-free`, Denoised),
               names_to = "step", values_to = "pct") |>
  mutate(step = factor(step, levels = c("Cutadapt","Merged","Chim-free","Denoised")))

p2 <- ggplot(retention, aes(x = step, y = pct, colour = batch, group = batch)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~instrument, ncol = 2) +
  scale_colour_manual(values = batch_cols) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(title    = "Mean read retention at each step (% of input)",
       subtitle = "Real samples only — one line per batch",
       x = NULL, y = "% of input reads", colour = "Batch / Run") +
  base_theme

# =============================================================================
# PLOT 3 — Merge rate distribution per batch
# =============================================================================
p3 <- real |>
  filter(!is.na(merge_pct)) |>
  ggplot(aes(x = batch, y = merge_pct, fill = batch)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6, colour = "grey30") +
  geom_hline(yintercept = MERGE_WARN, linetype = "dashed", colour = "red", linewidth = 0.6) +
  facet_wrap(~instrument, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = batch_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(title    = "Merge rate distribution per batch",
       subtitle = paste0("Dashed line = ", MERGE_WARN, "% warning threshold"),
       x = NULL, y = "Merge rate (% of input reads)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# =============================================================================
# PLOT 4 — Denoised rate distribution per batch
# =============================================================================
p4 <- real |>
  filter(!is.na(denoise_pct)) |>
  ggplot(aes(x = batch, y = denoise_pct, fill = batch)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6, colour = "grey30") +
  geom_hline(yintercept = DENOISE_WARN, linetype = "dashed", colour = "red", linewidth = 0.6) +
  facet_wrap(~instrument, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = batch_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(title    = "Denoised rate distribution per batch",
       subtitle = paste0("Dashed line = ", DENOISE_WARN, "% warning threshold"),
       x = NULL, y = "Denoised rate (% of input reads)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# =============================================================================
# PLOT 5 — ESV richness (n_esv) per batch
# =============================================================================
p5 <- real |>
  filter(!is.na(n_esv), n_esv > 0) |>
  ggplot(aes(x = batch, y = n_esv, fill = batch)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6, colour = "grey30") +
  facet_wrap(~instrument, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = batch_cols, guide = "none") +
  scale_y_log10(labels = label_comma()) +
  labs(title    = "ESV richness per sample per batch",
       subtitle = "Log scale — proxy for local diversity before MOTU clustering",
       x = NULL, y = "Number of ESVs (log scale)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# =============================================================================
# SAVE PDF
# =============================================================================
pdf(OUTPUT_PDF, width = 12, height = 7)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()

message("Saved: ", OUTPUT_PDF)
