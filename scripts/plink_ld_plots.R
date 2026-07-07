#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[grep("^--file=", full_args)][1]
script_name <- ifelse(is.na(script_arg), "this_script.R", basename(sub("^--file=", "", script_arg)))

if (length(args) < 2 || length(args) > 3) {
  stop(sprintf("Usage: Rscript %s <ld_decay.ld> <output_dir> [reference.fai]", script_name))
}

input_ld   <- args[[1]]
output_dir <- args[[2]]
input_fai  <- if (length(args) == 3) args[[3]] else NULL

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

ld <- read.table(input_ld, header = TRUE, stringsAsFactors = FALSE)
# Standardise column names (plink outputs CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2)
colnames(ld) <- toupper(colnames(ld))

if (!all(c("CHR_A", "BP_A", "CHR_B", "BP_B", "R2") %in% colnames(ld))) {
  stop("Unexpected columns in .ld file: ", paste(colnames(ld), collapse = ", "))
}

ld <- ld |>
  mutate(R2 = as.numeric(R2)) |>
  filter(!is.na(R2))

chromosomes <- sort(unique(c(ld$CHR_A, ld$CHR_B)))

# ── Chromosome lengths from .fai (cols: name, length, ...) ───────────────────
if (!is.null(input_fai)) {
  fai <- read.table(input_fai, header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("CHR", "length", "offset", "linebases", "linewidth"))
  chr_limits <- fai |>
    filter(CHR %in% chromosomes) |>
    select(CHR, length)
} else {
  chr_limits <- data.frame(
    CHR    = chromosomes,
    length = tapply(c(ld$BP_A, ld$BP_B), c(ld$CHR_A, ld$CHR_B), max)[chromosomes]
  )
}

# blank data to anchor each facet from 1 to chromosome length
chr_blanks <- bind_rows(
  chr_limits |> mutate(pos_x = 1L,      pos_y = 1L),
  chr_limits |> mutate(pos_x = length,  pos_y = length)
)

# ── Heatmap: pairwise r² per chromosome ──────────────────────────────────────
# Mirror A↔B so the heatmap is symmetric
ld_sym <- bind_rows(
  ld |> select(CHR = CHR_A, pos_x = BP_A, pos_y = BP_B, R2),
  ld |> select(CHR = CHR_B, pos_x = BP_B, pos_y = BP_A, R2)
) |>
  filter(CHR %in% chromosomes)

# Bin genomic positions: geom_tile defaults to width/height = 1 bp, which is
# invisible at chromosomal scales. Round to ~1/200 of each chromosome's span.
n_bins    <- 200L
bin_sizes <- ld_sym |>
  group_by(CHR) |>
  summarise(bin_size = diff(range(c(pos_x, pos_y))) / n_bins, .groups = "drop")
ld_binned <- ld_sym |>
  left_join(bin_sizes, by = "CHR") |>
  mutate(
    pos_x = round(pos_x / bin_size) * bin_size,
    pos_y = round(pos_y / bin_size) * bin_size
  ) |>
  group_by(CHR, pos_x, pos_y, bin_size) |>
  summarise(R2 = mean(R2, na.rm = TRUE), .groups = "drop")

p_heat <- ggplot(ld_binned, aes(x = pos_x, y = pos_y, fill = R2)) +
  geom_blank(data = chr_blanks, aes(x = pos_x, y = pos_y), inherit.aes = FALSE) +
  geom_tile(aes(width = bin_size, height = bin_size)) +
  facet_wrap(~CHR, scales = "free", ncol = 2) +
  scale_fill_gradientn(
    colours = c("white", "#4DBBD5", "#E64B35"),
    limits  = c(0, 1),
    name    = expression(r^2)
  ) +
  labs(
    x = "Position (bp)",
    y = "Position (bp)",
    title = "Pairwise LD heatmap"
  ) +
  theme_dark(base_size = 10) +
  theme(
    strip.text      = element_text(size = 8),
    panel.spacing   = unit(0.4, "lines"),
    axis.text       = element_text(size = 7)
  )

# ── Decay plot: mean r² vs pairwise distance ──────────────────────────────────
decay <- ld |>
  mutate(distance = abs(BP_B - BP_A)) |>
  group_by(CHR_A, distance) |>
  summarise(mean_r2 = mean(R2, na.rm = TRUE), .groups = "drop") |>
  rename(CHR = CHR_A)

p_decay <- ggplot(decay, aes(x = distance, y = mean_r2)) +
  geom_point(size = 0.6, alpha = 0.5, colour = "#4DBBD5") +
  geom_smooth(method = "loess", span = 0.1, se = FALSE,
              colour = "grey40", linewidth = 0.7, linetype = "dashed") +
  facet_wrap(~CHR, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Pairwise distance (bp)",
    y = expression(Mean~r^2),
    title = "LD decay"
  ) +
  theme_light(base_size = 11) +
  theme(strip.text = element_text(size = 8))

# ── Save ──────────────────────────────────────────────────────────────────────
n_chr    <- length(chromosomes)
heat_h   <- max(4, ceiling(n_chr / 2) * 4)
decay_h  <- max(4, ceiling(n_chr / 2) * 3)

ggsave(file.path(output_dir, "ld_heatmap.pdf"),   p_heat,  width = 10, height = heat_h,  limitsize = FALSE)
ggsave(file.path(output_dir, "ld_decay_plot.pdf"), p_decay, width = 10, height = decay_h, limitsize = FALSE)

message("Plots written to ", output_dir)
