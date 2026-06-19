#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[grep("^--file=", full_args)][1]
script_name <- ifelse(is.na(script_arg), "this_script.R", basename(sub("^--file=", "", script_arg)))

if (length(args) < 3L) {
  stop(sprintf("Usage: Rscript %s <input.vcf.gz> <output.tsv> <output.pdf> [recombination_threshold]", script_name))
}

input_vcf               <- args[[1L]]
output_tsv              <- args[[2L]]
output_pdf              <- args[[3L]]
recombination_threshold <- if (length(args) >= 4L) as.numeric(args[[4L]]) else 0.9

suppressPackageStartupMessages({
  library(vcfR)
  library(ape)
  library(pegas)
})

vcf <- read.vcfR(input_vcf, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

if (is.null(dim(gt)) || ncol(gt) == 0 || nrow(gt) == 0) {
  stop("No genotype data found in VCF.")
}

# Replace missing/diploid GTs; genomes are always haploid
gt[gt == "." | gt == "./." | gt == ".|."] <- NA
gt <- sub("[/|].*$", "", gt)   # "0/0" -> "0", already-haploid "0" stays "0"

message("GT unique values: ", paste(sort(unique(as.vector(gt))), collapse = ", "))
message("Variants: ", nrow(gt), "  Samples: ", ncol(gt))

# ── Filter to last generation ─────────────────────────────────────────────────
# Sample names are expected to follow pop_<N>_gen_<N>_genome_<N>.
# If parseable, retain only samples from the final generation.
sample_names <- colnames(gt)
gen_vals <- suppressWarnings(
  as.integer(sub(".*_gen_(\\d+)_genome_.*", "\\1", sample_names))
)
if (!all(is.na(gen_vals))) {
  last_gen   <- max(gen_vals, na.rm = TRUE)
  keep       <- !is.na(gen_vals) & gen_vals == last_gen
  gt         <- gt[, keep, drop = FALSE]
  message(sprintf("Filtering to last generation (%d). Samples retained: %d / %d.",
                  last_gen, sum(keep), length(keep)))
} else {
  message("Sample names do not contain generation info; using all samples.")
}

# gt is variants x samples; transpose to samples x variants
allele_mat <- t(gt)

# In a variant-only VCF, absence of a call means the sample carries the
# reference allele ("0") — fill NA accordingly rather than discarding samples.
n_missing <- rowSums(is.na(allele_mat))
message("Samples with missing calls (treated as REF): ", sum(n_missing > 0), " / ", nrow(allele_mat))
message("Samples retained: ", nrow(allele_mat))

# Encode allele indices as distinct DNA bases for pegas
encode <- function(x) {
  alleles <- sort(unique(na.omit(as.vector(x))))
  bases   <- c("a", "c", "g", "t")[seq_along(alleles)]
  out     <- x
  for (i in seq_along(alleles)) out[x == alleles[i]] <- bases[i]
  out
}
dna_mat <- encode(allele_mat)
message("DNA base unique values: ", paste(sort(unique(as.vector(dna_mat))), collapse = ", "))

# Build DNAbin as a named list of character vectors (reliable path)
dna_list <- lapply(seq_len(nrow(dna_mat)), function(i) dna_mat[i, ])
names(dna_list) <- rownames(dna_mat)
dna_bin <- as.matrix(as.DNAbin(dna_list))
message("DNAbin dimensions: ", nrow(dna_bin), " x ", ncol(dna_bin))

# ── Compute haplotypes manually ───────────────────────────────────────────────
# Collapse each sample's alleles to a string and count unique patterns
hap_strings  <- apply(dna_mat, 1, paste, collapse = "")
unique_haps  <- sort(unique(hap_strings))
freq_table   <- table(factor(hap_strings, levels = unique_haps))
freq         <- as.integer(freq_table)
n_haps       <- length(unique_haps)
total_haplotypes <- sum(freq)
message("Unique haplotypes: ", n_haps, "  Total individuals: ", total_haplotypes)

if (total_haplotypes <= 0 || n_haps == 0) {
  stop("No haplotypes were inferred from the provided VCF genotypes.")
}

hap_labels <- paste0("H", seq_len(n_haps))

summary_df <- data.frame(
  haplotype = hap_labels,
  count     = freq,
  frequency = freq / total_haplotypes
)

# ── Haplotype network ─────────────────────────────────────────────────────────
# Build a haplotype-class DNAbin object manually so haploNet() accepts it
unique_row_idx <- match(unique_haps, hap_strings)
hap_dnabin     <- dna_bin[unique_row_idx, , drop = FALSE]
rownames(hap_dnabin) <- hap_labels
index          <- lapply(unique_haps, function(h) which(hap_strings == h))
class(hap_dnabin) <- c("haplotype", "DNAbin")
attr(hap_dnabin, "freq")  <- freq
attr(hap_dnabin, "index") <- index

net <- haploNet(hap_dnabin)

# ── Recombination detection (profile-based) ───────────────────────────────────
# For each haplotype build a mutation profile — a set of "pos:allele" tokens
# for every site that carries a non-reference allele (GT != "0").  A haplotype
# is flagged as recombinant when two others (potential parents) can be
# identified: each parent must have >= recombination_threshold of its mutations
# present in the candidate, and the two parents must each carry at least one
# mutation the other lacks.

# Haplotype allele matrix (rows = unique haplotypes, cols = sites)
hap_alleles <- allele_mat[unique_row_idx, , drop = FALSE]

# Named list: haplotype label -> character vector of "pos:allele" mutation tokens
# allele_mat contains raw GT allele indices ("0" = REF, "1", "2", ... = ALT)
hap_profiles <- setNames(
  lapply(seq_len(n_haps), function(i) {
    hv  <- hap_alleles[i, ]
    idx <- which(hv != "0" & !is.na(hv))
    if (length(idx) == 0L) return(character(0L))
    paste0(idx, ":", hv[idx])
  }),
  hap_labels
)

find_recombinant_parents <- function(profile_c, all_profiles_named, threshold = 0.9) {
  if (threshold <= 0 || length(all_profiles_named) < 2L) return(NULL)
  candidate_names <- Filter(
    function(nm) {
      a <- all_profiles_named[[nm]]
      length(a) > 0L && (sum(a %in% profile_c) / length(a)) >= threshold
    },
    names(all_profiles_named)
  )
  if (length(candidate_names) < 2L) return(NULL)
  n <- length(candidate_names)
  for (i in seq_len(n - 1L)) {
    a <- all_profiles_named[[candidate_names[i]]]
    for (j in seq(i + 1L, n)) {
      b <- all_profiles_named[[candidate_names[j]]]
      if (!any(!a %in% b) || !any(!b %in% a)) next
      return(c(candidate_names[i], candidate_names[j]))
    }
  }
  NULL
}

# For each haplotype check all others as potential parents
recomb_node_set      <- integer(0L)
recomb_arcs          <- list()
recomb_parent_labels <- vector("list", n_haps)   # index -> c(p1_label, p2_label)

for (i in seq_len(n_haps)) {
  others  <- hap_profiles[hap_labels[-i]]
  parents <- find_recombinant_parents(hap_profiles[[i]], others,
                                      threshold = recombination_threshold)
  if (!is.null(parents)) {
    recomb_node_set           <- union(recomb_node_set, i)
    recomb_parent_labels[[i]] <- parents
    p1 <- match(parents[1L], hap_labels)
    p2 <- match(parents[2L], hap_labels)
    recomb_arcs <- c(recomb_arcs, list(c(i, p1)), list(c(i, p2)))
  }
}

n_recombinants <- length(recomb_node_set)
message("Recombinant haplotypes detected: ", n_recombinants)

# Add parents column to summary and write TSV
summary_df$parents <- vapply(seq_len(n_haps), function(i) {
  p <- recomb_parent_labels[[i]]
  if (is.null(p)) NA_character_ else paste(p, collapse = ",")
}, character(1L))
write.table(summary_df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# ── Plots (single PDF, one page per figure) ───────────────────────────────────

node_bg                  <- rep("white", n_haps)
node_bg[recomb_node_set] <- "tomato"

# Quadratic bezier arc — draws a curved dashed line between two points
draw_arc <- function(x0, y0, x1, y1, col = "red2", lwd = 1.5, lty = 2,
                     n = 60L, curve = 0.25) {
  dx <- x1 - x0; dy <- y1 - y0
  mx <- (x0 + x1) / 2 - dy * curve
  my <- (y0 + y1) / 2 + dx * curve
  t  <- seq(0, 1, length.out = n)
  lines((1-t)^2*x0 + 2*(1-t)*t*mx + t^2*x1,
        (1-t)^2*y0 + 2*(1-t)*t*my + t^2*y1,
        col = col, lwd = lwd, lty = lty)
}

pdf(output_pdf, width = 8, height = 7)

# Page 1: Haplotype network — coloured nodes + dashed arcs for recombination
# plot.haploNet returns node coordinates invisibly (rows 1:n_haps = real haplotypes)
if (n_haps <= 1L) {
  plot.new()
  title(main = "Haplotype network")
  text(0.5, 0.5,
       sprintf("Only %d unique haplotype — network cannot be drawn.", n_haps),
       cex = 1.2)
  coords <- NULL
} else {
  coords <- plot(
    net,
    size          = sqrt(freq),
    bg            = node_bg,
    labels        = TRUE,
    show.mutation = 1,
    main          = paste0("Haplotype network  (recombinants = ", n_recombinants, ")")
  )
}
if (n_recombinants > 0L && !is.null(coords)) {
  coord_labels <- rownames(coords)
  for (pair in recomb_arcs) {
    r1 <- match(hap_labels[pair[1L]], coord_labels)
    r2 <- match(hap_labels[pair[2L]], coord_labels)
    if (is.na(r1) || is.na(r2)) next
    draw_arc(coords[r1, 1L], coords[r1, 2L],
             coords[r2, 1L], coords[r2, 2L])
  }
}

# Page 2: Recombinant haplotype summary
plot.new()
title(main = paste0("Recombinant haplotypes  (threshold = ", recombination_threshold, ")"))
if (n_recombinants == 0L) {
  text(0.5, 0.5, "No recombinant haplotypes detected", cex = 1.2)
} else {
  lines_out <- vapply(recomb_node_set, function(i) {
    parents <- recomb_parent_labels[[i]]
    sprintf("%s  <-  %s + %s", hap_labels[i], parents[1L], parents[2L])
  }, character(1L))
  text(0.5, seq(0.9, 0.1, length.out = length(lines_out)),
       lines_out, cex = 1.0)
}
invisible(dev.off())