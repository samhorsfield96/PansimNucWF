#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[grep("^--file=", full_args)][1]
script_name <- ifelse(is.na(script_arg), "this_script.R", basename(sub("^--file=", "", script_arg)))

if (length(args) != 3) {
  stop(sprintf("Usage: Rscript %s <input.vcf.gz> <output.tsv> <output.pdf>", script_name))
}

input_vcf  <- args[[1]]
output_tsv <- args[[2]]
output_pdf <- args[[3]]

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

# gt is variants x samples; transpose to samples x variants
allele_mat <- t(gt)

# Drop samples with any missing allele
n_missing <- rowSums(is.na(allele_mat))
message("Samples with missing calls: ", sum(n_missing > 0), " / ", nrow(allele_mat))
complete   <- n_missing == 0
if (sum(complete) == 0) stop("All samples have missing genotypes at one or more sites.")
allele_mat <- allele_mat[complete, , drop = FALSE]
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

write.table(summary_df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

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

# ── Recombination detection (four-gamete test + Hudson-Kaplan Rmin) ───────────
# Four-gamete test: a pair of biallelic sites is incompatible when all four
# two-locus gametes (00, 01, 10, 11) are observed, implying recombination
# (or recurrent mutation) between them.
four_gamete_test <- function(mat) {
  n_sites <- ncol(mat)
  ig <- matrix(FALSE, nrow = n_sites, ncol = n_sites)
  for (i in seq_len(n_sites - 1L)) {
    for (j in seq(i + 1L, n_sites)) {
      if (length(unique(paste0(mat[, i], mat[, j]))) == 4L)
        ig[i, j] <- ig[j, i] <- TRUE
    }
  }
  ig
}

# Rmin (Hudson & Kaplan 1985): minimum recombination events via greedy
# interval point cover — each incompatible pair [i,j] is an interval that
# must be "hit" by at least one recombination event.
rmin_estimate <- function(ig) {
  pairs <- which(ig & upper.tri(ig), arr.ind = TRUE)
  if (nrow(pairs) == 0L) return(0L)
  pairs <- pairs[order(pairs[, 2L]), , drop = FALSE]
  n_events <- 0L; last_point <- -1L
  for (k in seq_len(nrow(pairs))) {
    if (last_point < pairs[k, 1L]) {
      n_events   <- n_events + 1L
      last_point <- pairs[k, 2L]
    }
  }
  n_events
}

ig             <- four_gamete_test(allele_mat)
rmin           <- rmin_estimate(ig)
n_incompatible <- sum(ig[upper.tri(ig)])

message("Minimum recombination events (Rmin): ", rmin)
message("Incompatible site pairs (four-gamete test): ", n_incompatible)

# ── Plots (single PDF, one page per figure) ───────────────────────────────────

# Haplotype allele matrix (rows = unique haplotypes, cols = sites)
hap_alleles <- allele_mat[unique_row_idx, , drop = FALSE]

# Build recombination annotations from the incompatibility matrix:
#   recomb_node_set — haplotype indices carrying the minority (putative
#                     recombinant) gamete at any incompatible site pair.
#   recomb_arcs     — pairs of indices to connect with dashed arcs; each
#                     recombinant node is linked to its closest neighbour
#                     outside its gamete class.
recomb_node_set <- integer(0)
recomb_arcs     <- list()
drawn_conn      <- matrix(FALSE, n_haps, n_haps)

if (n_incompatible > 0L) {
  inc_pairs <- which(ig & upper.tri(ig), arr.ind = TRUE)
  for (k in seq_len(nrow(inc_pairs))) {
    si <- inc_pairs[k, 1L]; sj <- inc_pairs[k, 2L]
    gametes  <- paste0(hap_alleles[, si], hap_alleles[, sj])
    gtab     <- sort(table(gametes))
    rare_g   <- names(gtab)[1L]          # putative recombinant gamete
    rare_idx <- which(gametes == rare_g)
    recomb_node_set <- union(recomb_node_set, rare_idx)

    # Parental gametes: each shares exactly one allele with the rare gamete.
    # e.g. rare = "01"  →  parent A has "00" (gave allele at si),
    #                        parent B has "11" (gave allele at sj).
    g_split <- strsplit(rare_g, "")[[1L]]
    a_si    <- g_split[1L]
    a_sj    <- g_split[2L]
    other_si <- setdiff(unique(substr(gametes, 1L, 1L)), a_si)[1L]
    other_sj <- setdiff(unique(substr(gametes, 2L, 2L)), a_sj)[1L]
    parent_A_g <- paste0(a_si,    other_sj)   # same allele at si
    parent_B_g <- paste0(other_si, a_sj)      # same allele at sj

    pA <- which(gametes == parent_A_g)
    pB <- which(gametes == parent_B_g)
    if (length(pA) == 0L || length(pB) == 0L) next

    # Pick the most-frequent representative from each parental group
    h1 <- pA[which.max(freq[pA])]
    h2 <- pB[which.max(freq[pB])]
    if (!drawn_conn[h1, h2]) {
      drawn_conn[h1, h2] <- drawn_conn[h2, h1] <- TRUE
      recomb_arcs <- c(recomb_arcs, list(c(h1, h2)))
    }
  }
}

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
coords <- plot(
  net,
  size          = sqrt(freq),
  bg            = node_bg,
  labels        = TRUE,
  show.mutation = 1,
  main          = paste0("Haplotype network  (Rmin = ", rmin, ")")
)
if (n_incompatible > 0L && isTRUE(nrow(coords) >= n_haps)) {
  for (pair in recomb_arcs) {
    draw_arc(coords[pair[1L], 1L], coords[pair[1L], 2L],
             coords[pair[2L], 1L], coords[pair[2L], 2L])
  }
}
legend(
  "bottomright",
  legend = c(paste0(hap_labels, " (n=", freq, ")"),
             if (length(recomb_node_set) > 0L) "Putative recombinant" else NULL,
             if (length(recomb_arcs)     > 0L) "Parental source pair" else NULL),
  pch    = c(rep(21L, n_haps),
             if (length(recomb_node_set) > 0L) 21L else NULL,
             if (length(recomb_arcs)     > 0L) NA  else NULL),
  lty    = c(rep(NA,  n_haps),
             if (length(recomb_node_set) > 0L) NA else NULL,
             if (length(recomb_arcs)     > 0L) 2L else NULL),
  pt.bg  = c(node_bg,
             if (length(recomb_node_set) > 0L) "tomato" else NULL,
             if (length(recomb_arcs)     > 0L) NA       else NULL),
  col    = c(rep("black", n_haps),
             if (length(recomb_node_set) > 0L) "black" else NULL,
             if (length(recomb_arcs)     > 0L) "red2"  else NULL),
  bty    = "n",
  cex    = 0.7
)

# Page 2: Four-gamete incompatibility matrix
if (n_incompatible > 0L) {
  image(
    seq_len(ncol(ig)), seq_len(nrow(ig)), ig,
    col  = c("white", "steelblue"),
    xlab = "Site index", ylab = "Site index",
    main = paste0("Four-gamete incompatibility  (Rmin = ", rmin, ")")
  )
} else {
  plot.new()
  title(main = paste0("Four-gamete incompatibility  (Rmin = ", rmin, ")"))
  text(0.5, 0.5, "No incompatible site pairs detected", cex = 1.2)
}
invisible(dev.off())