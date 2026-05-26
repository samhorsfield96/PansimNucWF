library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggsci)
library(vcfR)

# Usage:
#   Rscript plot_sfs_nuc.R <vcf_file> [gff_annotation] [output_prefix]
#
# Reads a multi-sample VCF where sample names follow the convention
# pop_<pop>_gen_<gen>_genome_<id>.  For each (population, generation) group,
# computes per-site allele frequencies and plots the site frequency spectrum
# (SFS) as a histogram.
# If a GFF annotation file is provided, sites are grouped by the feature_type
# of the overlapping element; otherwise all sites are labelled "all".
#
# args[2] is interpreted as a GFF file if it is an existing file path, and as
# the output prefix otherwise (allowing the GFF to be omitted).

args       <- commandArgs(trailingOnly = TRUE)
args       <- args[!grepl("^--", args)]
if (length(args) < 1L)
  stop("Usage: Rscript plot_sfs_nuc.R <vcf_file> [gff_annotation] [output_prefix]")

vcf_file   <- args[1L]
gff_file   <- if (length(args) >= 2L && nchar(args[2L]) > 0L && file.exists(args[2L])) args[2L] else NULL
out_prefix <- if (!is.null(gff_file)) {
  if (length(args) >= 3L) args[3L] else "sfs_nuc"
} else {
  if (length(args) >= 2L && nchar(args[2L]) > 0L) args[2L] else "sfs_nuc"
}

# ── Attribute parser ──────────────────────────────────────────────────────────
parse_attrs <- function(attr_str) {
  pairs <- strsplit(attr_str, ";", fixed = TRUE)[[1L]]
  kv    <- strsplit(pairs, "=", fixed = TRUE)
  keys  <- vapply(kv, `[[`, character(1L), 1L)
  vals  <- vapply(kv, function(x) if (length(x) >= 2L) x[[2L]] else NA_character_,
                  character(1L))
  setNames(vals, keys)
}

# ── GFF annotation reader ─────────────────────────────────────────────────────
# Returns a data.frame with columns: chrom, start, end, feature_type.
read_annotation_gff <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nchar(lines) > 0L & !startsWith(lines, "#")]
  if (length(lines) == 0L) return(NULL)
  rows <- lapply(lines, function(line) {
    f  <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(f) < 9L) return(NULL)
    a  <- parse_attrs(f[9L])
    ft <- a[["feature_type"]]
    if (is.na(ft)) return(NULL)
    data.frame(
      chrom        = f[1L],
      start        = as.integer(f[4L]),
      end          = as.integer(f[5L]),
      feature_type = ft,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(Filter(Negate(is.null), rows))
}

# ── Read VCF ──────────────────────────────────────────────────────────────────
message("Reading VCF: ", vcf_file)
vcf    <- read.vcfR(vcf_file, verbose = FALSE)
gt_mat <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

if (ncol(gt_mat) == 0L) stop("No samples found in VCF.")

variants <- data.frame(
  chrom = getCHROM(vcf),
  pos   = as.integer(getPOS(vcf)),
  stringsAsFactors = FALSE
)
message("Loaded ", nrow(variants), " variant site(s) across ",
        ncol(gt_mat), " sample(s).")

# ── Parse sample names ────────────────────────────────────────────────────────
# Expected format: pop_<pop>_gen_<gen>_genome_<id>
sample_meta <- data.frame(sample_name = colnames(gt_mat),
                          stringsAsFactors = FALSE)
sample_meta$pop_id    <- suppressWarnings(
  as.integer(sub(".*pop_(\\d+).*",    "\\1", sample_meta$sample_name)))
sample_meta$gen_id    <- suppressWarnings(
  as.integer(sub(".*gen_(\\d+).*",    "\\1", sample_meta$sample_name)))
sample_meta$genome_id <- suppressWarnings(
  as.integer(sub(".*genome_(\\d+).*", "\\1", sample_meta$sample_name)))

bad <- is.na(sample_meta$pop_id) | is.na(sample_meta$gen_id)
if (any(bad)) {
  warning("Excluding ", sum(bad), " sample(s) whose names do not match ",
          "pop_<N>_gen_<N>_genome_<N>: ",
          paste(sample_meta$sample_name[bad], collapse = ", "))
  sample_meta <- sample_meta[!bad, ]
}
if (nrow(sample_meta) == 0L) stop("No valid samples after name parsing.")

# ── Compute per-site allele frequencies ──────────────────────────────────────
message("Computing per-site allele frequencies...")

groups <- distinct(sample_meta, pop_id, gen_id)

sfs_list <- lapply(seq_len(nrow(groups)), function(i) {
  pop <- groups$pop_id[i]
  gen <- groups$gen_id[i]

  group_samples <- sample_meta$sample_name[
    sample_meta$pop_id == pop & sample_meta$gen_id == gen]
  gt_sub <- gt_mat[, group_samples, drop = FALSE]

  # For each variant site, split GT strings into alleles and compute frequencies
  site_df <- bind_rows(lapply(seq_len(nrow(gt_sub)), function(j) {
    gt_vec  <- gt_sub[j, , drop = TRUE]
    alleles <- unlist(strsplit(
      gt_vec[!is.na(gt_vec) & gt_vec != "./." & gt_vec != "."],
      "[/|]"
    ))
    alleles <- alleles[alleles != "."]
    if (length(alleles) < 2L)
      return(data.frame(major_freq = NA_real_, minor_freq = NA_real_))
    freq <- sort(as.numeric(table(alleles)) / length(alleles), decreasing = TRUE)
    data.frame(
      major_freq = freq[[1L]],
      minor_freq = if (length(freq) > 1L) freq[[2L]] else 0.0
    )
  }))

  cbind(
    data.frame(pop_id = pop, gen_id = gen,
               chrom  = variants$chrom, pos = variants$pos,
               stringsAsFactors = FALSE),
    site_df
  )
})

all_sites <- bind_rows(sfs_list) |>
  filter(!is.na(major_freq), minor_freq > 0)

if (nrow(all_sites) == 0L) stop("No variable sites found in VCF.")

# ── Annotate variants with feature types ─────────────────────────────────────
if (!is.null(gff_file)) {
  message("Annotating variants with features from: ", gff_file)
  gff <- read_annotation_gff(gff_file)

  all_sites$feature_type <- mapply(function(chr, pos) {
    hits <- gff$feature_type[
      gff$chrom == chr & gff$start <= pos & gff$end >= pos]
    if (length(hits) == 0L) "intergenic" else hits[[1L]]
  }, all_sites$chrom, all_sites$pos)
} else {
  all_sites$feature_type <- "all"
}

sfs_data <- all_sites

# ── Plot ──────────────────────────────────────────────────────────────────────
n_genomes_per_group <- sample_meta |>
  group_by(pop_id, gen_id) |>
  summarise(n = n(), .groups = "drop")

n_bins <- max(n_genomes_per_group$n)

stacked_sfs_data <- melt(sfs_data, measure.vars = c("major_freq", "minor_freq"))
stacked_sfs_data$variable <- as.character(stacked_sfs_data$variable)
stacked_sfs_data$variable[stacked_sfs_data$variable == "major_freq"] <- "Major allele"
stacked_sfs_data$variable[stacked_sfs_data$variable == "minor_freq"] <- "Minor allele"

p_minor_density <- ggplot(sfs_data, aes(x = minor_freq, fill = feature_type)) +
  geom_histogram(
    bins     = n_bins,
    aes(y = ..density..),
    boundary = 0,
    colour   = NA,
    position = "identity",
    alpha    = 0.8
  ) +
  facet_grid(
    rows = vars(interaction(pop_id, gen_id,
                             sep = " / gen=", lex.order = TRUE)),
    cols = vars(feature_type),
    labeller = labeller(
      .rows = function(x) paste0("pop=", sub(" / gen=", "  gen=", x))
    ),
    scales = "free_y"
  ) +
  scale_x_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x     = "Minor allele frequency",
    y     = "Density"
  ) +
  scale_fill_npg() +
  theme_light(base_size = 11) +
  theme(legend.position = "none",
        strip.text      = element_text(size = 9))

p_both_density <- ggplot(stacked_sfs_data, aes(x = value, fill = variable)) +
  geom_histogram(
    bins     = n_bins,
    aes(y = ..density..),
    boundary = 0,
    colour   = NA,
    position = "identity",
    alpha    = 0.8
  ) +
  facet_grid(
    rows = vars(interaction(pop_id, gen_id,
                            sep = " / gen=", lex.order = TRUE)),
    cols = vars(feature_type),
    labeller = labeller(
      .rows = function(x) paste0("pop=", sub(" / gen=", "  gen=", x))
    ),
    scales = "free_y"
  ) +
  scale_x_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x     = "Allele frequency",
    y     = "Density",
    fill  = "Allele type"
  ) +
  scale_fill_npg() +
  theme_light(base_size = 11) +
  theme(strip.text = element_text(size = 9))

# ── Save ──────────────────────────────────────────────────────────────────────
n_pops  <- n_distinct(sfs_data$pop_id)
n_gens  <- n_distinct(sfs_data$gen_id)
n_types <- n_distinct(sfs_data$feature_type)

pdf_w <- max(9, 4.5 * n_types)
pdf_h <- max(6, 3 * n_pops * n_gens)

out_pdf <- paste0(out_prefix, "_density_minor_alleles.pdf")
ggsave(out_pdf, plot = p_minor_density, width = pdf_w, height = pdf_h)
message("Saved: ", out_pdf)

out_pdf <- paste0(out_prefix, "_density_all_alleles.pdf")
ggsave(out_pdf, plot = p_both_density, width = pdf_w, height = pdf_h)
message("Saved: ", out_pdf)

out_csv <- paste0(out_prefix, "_sfs.csv")
write.csv(sfs_data, out_csv, row.names = FALSE)
message("Saved: ", out_csv)
