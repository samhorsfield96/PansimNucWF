#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[grep("^--file=", full_args)][1]
script_name <- ifelse(is.na(script_arg), "script", basename(sub("^--file=", "", script_arg)))

if (length(args) != 2) {
  stop(sprintf("Usage: Rscript %s <input.vcf.gz> <output.tsv>", script_name))
}

input_vcf <- args[[1]]
output_tsv <- args[[2]]

suppressPackageStartupMessages({
  library(vcfR)
  library(pegas)
})

vcf <- read.vcfR(input_vcf, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

if (is.null(dim(gt)) || ncol(gt) == 0 || nrow(gt) == 0) {
  stop("No genotype data found in VCF.")
}

gt <- gsub("\\|", "/", gt)
gt_df <- as.data.frame(t(gt))
loci_data <- as.loci(gt_df, allele.sep = "/")
haps <- haplotype(loci_data)

freq <- attr(haps, "freq")
total_count <- sum(freq)

if (total_count <= 0) {
  stop("No haplotypes were inferred from the provided VCF genotypes.")
}

summary_df <- data.frame(
  haplotype = labels(haps),
  count = as.integer(freq),
  frequency = as.numeric(freq) / total_count
)

write.table(summary_df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
