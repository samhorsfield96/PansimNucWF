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

pdf(output_pdf, width = 8, height = 7)
plot(
  net,
  size          = sqrt(freq),
  labels        = TRUE,
  show.mutation = 1,
  main          = "Haplotype network"
)
legend(
  "bottomright",
  legend = paste0(hap_labels, " (n=", freq, ")"),
  bty    = "n",
  cex    = 0.7
)
invisible(dev.off())
