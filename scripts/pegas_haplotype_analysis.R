args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript <script> <input.vcf.gz> <output.tsv>")
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
