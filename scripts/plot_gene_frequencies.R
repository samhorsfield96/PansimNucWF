#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggsci)
})

# Usage:
#   Rscript plot_gene_frequencies.R <gff_dir> <output_dir> [options]
#
# Options:
#   --gene-type TYPE       GFF feature type used for genes (default: gene)
#   --final-generation      restrict analysis to the final generation

args <- commandArgs(trailingOnly = TRUE)

take_flag <- function(flag, args, default) {
  index <- match(flag, args)
  if (is.na(index)) return(list(value = default, args = args))
  if (index == length(args)) stop(flag, " requires a value.")
  list(value = args[[index + 1L]], args = args[-c(index, index + 1L)])
}

take_bool <- function(flag, args, default = FALSE) {
  index <- match(flag, args)
  if (is.na(index)) return(list(value = default, args = args))
  next_arg <- if (index < length(args)) tolower(args[[index + 1L]]) else ""
  if (next_arg %in% c("true", "false")) {
    return(list(value = identical(next_arg, "true"),
                args = args[-c(index, index + 1L)]))
  }
  list(value = TRUE, args = args[-index])
}

if (length(args) < 2L) {
  stop("Usage: Rscript plot_gene_frequencies.R <gff_dir> <output_dir> ",
       "[--gene-type TYPE] [--final-generation]")
}

gff_dir <- args[[1L]]
output_dir <- args[[2L]]
args <- args[-c(1L, 2L)]

flag <- take_flag("--gene-type", args, "gene")
gene_type <- flag$value
args <- flag$args
flag <- take_bool("--final-generation", args)
final_generation_only <- flag$value
args <- flag$args

if (length(args) > 0L) {
  stop("Unknown argument: ", paste(args, collapse = " "))
}
if (!dir.exists(gff_dir)) stop("GFF directory does not exist: ", gff_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gff_files <- list.files(
  gff_dir,
  pattern="^pop_\\d+_gen_\\d+_genome_\\d+\\.gff",
  full.names=TRUE
)

if (final_generation_only) {
  gens <- as.integer(sub(".*_gen_(\\d+)_genome_.*", "\\1",
                         basename(gff_files)))
  last_gen <- max(gens, na.rm=TRUE)
  gff_files <- gff_files[gens == last_gen]
  message("Restricting to generation ", last_gen)
}

if (length(gff_files) == 0L) stop("No .gff or .gff.gz files found in ", gff_dir)

read_gff <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    dt <- fread(cmd = paste("gzip -cd --", shQuote(path)), sep = "\t",
                header = FALSE, comment.char = "#", showProgress = FALSE)
  } else {
    dt <- fread(path, sep = "\t", header = FALSE, comment.char = "#",
                showProgress = FALSE)
  }
  if (nrow(dt) == 0L) return(NULL)
  if (ncol(dt) < 5L) stop("GFF has fewer than five columns: ", path)
  setnames(dt, c("seqid", "source", "feature_type", "start", "end",
                 "score", "strand", "phase", "attributes")[seq_len(ncol(dt))])
  dt[, start := as.numeric(start)]
  dt[, end := as.numeric(end)]
  dt[, attribute_feature_type := sub(
    ".*(?:^|;)feature_type=([^;]+).*", "\\1", attributes, perl = TRUE
  )]
  dt[, attribute_feature_id := sub(
    ".*(?:^|;)feature_id=([^;]+).*", "\\1", attributes, perl = TRUE
  )]
  dt[, attribute_feature_broken := sub(
    ".*(?:^|;)feature_broken=([^;]+).*", "\\1", attributes, perl = TRUE
  )]
  dt[is.na(attributes) | !grepl("(?:^|;)feature_type=", attributes, perl = TRUE),
     attribute_feature_type := NA_character_]
  dt[is.na(attributes) | !grepl("(?:^|;)feature_id=", attributes, perl = TRUE),
     attribute_feature_id := NA_character_]
  dt[is.na(attributes) | !grepl("(?:^|;)feature_broken=", attributes, perl = TRUE),
     attribute_feature_broken := NA_character_]
  dt[, feature_broken := tolower(attribute_feature_broken) == "true"]
  dt[is.finite(start) & is.finite(end) & start <= end]
}

build_genes <- function(features) {
  # In PansimNuc annotations, exon/intron rows belonging to one gene share
  # the same attributes:feature_id value. Collapse those rows into one span.
  grouped_rows <- features[
    feature_type %in% c("exon", "intron") & !is.na(attribute_feature_type)
  ]
  if (nrow(grouped_rows) > 0L) {
    genes <- grouped_rows[
      , .(
        start = min(start),
        end = max(end),
        # a gene is broken if any of its exon/intron rows are broken
        broken = any(feature_broken, na.rm = TRUE)
      ),
      by = .(seqid, gene_id = attribute_feature_id)
    ]
  } else {
    genes <- features[feature_type == gene_type,
                      .(seqid, start, end, gene_id = attribute_feature_id,
                        broken = feature_broken)]
  }
  genes <- genes[broken != TRUE]
  genes[, broken := NULL]
  genes
}

ids_from_path <- function(path) {
  ids <- regmatches(
    basename(path),
    regexec("pop_(\\d+)_gen_(\\d+)_genome_(\\d+)", basename(path))
  )[[1]]
  list(
    pop_id = as.integer(ids[2]),
    generation = as.integer(ids[3]),
    genome_id = as.integer(ids[4])
  )
}

parse_gff <- function(path) {
  features <- read_gff(path)
  if (is.null(features)) {
    warning("Skipping empty GFF: ", path)
    return(NULL)
  }

  genes <- build_genes(features)
  genes <- genes[!is.na(gene_id)]
  if (nrow(genes) == 0L) {
    warning("Skipping GFF with no grouped exon/intron or ", gene_type,
            " features: ", path)
    return(NULL)
  }

  ids <- ids_from_path(path)
  unique(genes[, .(
    pop_id = ids$pop_id,
    generation = ids$generation,
    genome_id = ids$genome_id,
    gene_id
  )])
}

message("Parsing GFFs...")

all_data <- rbindlist(
  lapply(gff_files, parse_gff),
  use.names=TRUE,
  fill=TRUE
)

if (nrow(all_data) == 0)
  stop("No '", gene_type, "' features found.")

# number of genomes sampled per population/generation
genome_counts <- unique(
  all_data[, .(pop_id, generation, genome_id)]
)[
  ,
  .(n_genomes=.N),
  by=.(pop_id, generation)
]

gene_counts <- all_data[
  ,
  .(n_present=.N),
  by=.(pop_id, generation, gene_id)
]

gene_frequencies <- merge(
  gene_counts,
  genome_counts,
  by=c("pop_id", "generation")
)

gene_frequencies[, frequency := n_present / n_genomes]

freq_dist <- gene_frequencies[
  ,
  .(n_genes=.N),
  by=.(pop_id,
       generation,
       frequency)
]

fwrite(gene_frequencies,
       file.path(output_dir, "gene_frequencies_per_gene.csv"))

fwrite(freq_dist,
       file.path(output_dir, "gene_frequencies_distribution.csv"))

# ── Plots ─────────────────────────────────────────────────────────────────────

p_dist <- ggplot(gene_frequencies,
                 aes(x = frequency)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = pal_npg()(1)) +
  facet_wrap(generation ~ pop_id, labeller = label_both, scales = "free_y") +
  labs(
    x = "Gene frequency", y = "Number of genes") +
  theme_light()

ggsave(file.path(output_dir, "gene_frequencies_frequency_dist.pdf"),
       p_dist, width = 10, height = 6)
