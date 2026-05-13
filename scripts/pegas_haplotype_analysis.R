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
  library(pegas)
})

vcf <- read.vcfR(input_vcf, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

if (is.null(dim(gt)) || ncol(gt) == 0 || nrow(gt) == 0) {
  stop("No genotype data found in VCF.")
}

# Replace missing genotypes with NA
gt[gt == "." | gt == "./."] <- NA
gt_df <- as.data.frame(t(gt))
# Haploid GT values are single alleles ("0", "1"); no allele separator needed
loci_data <- as.loci(gt_df)
haps <- haplotype(loci_data)

freq <- attr(haps, "freq")
total_haplotypes <- sum(freq)

if (total_haplotypes <= 0) {
  stop("No haplotypes were inferred from the provided VCF genotypes.")
}

summary_df <- data.frame(
  haplotype = labels(haps),
  count = as.integer(freq),
  frequency = as.numeric(freq) / total_haplotypes
)

write.table(summary_df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# ── Haplotype network ─────────────────────────────────────────────────────────
net <- haploNet(haps)

pdf(output_pdf, width = 8, height = 7)
plot(
  net,
  size       = sqrt(freq),        # node area proportional to count
  pie        = freq / total_haplotypes,
  labels     = TRUE,
  show.mutation = 1,
  main       = "Haplotype network"
)
legend(
  "bottomright",
  legend = paste0(labels(haps), " (n=", as.integer(freq), ")"),
  bty    = "n",
  cex    = 0.7
)
invisible(dev.off())
