#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Usage:
#   Rscript plot_gene_feature_distances.R <gff_dir> <output_dir> [options]
#
# Options:
#   --gene-type TYPE       GFF feature type used for genes (default: gene)
#   --te-types TYPE,...     comma-separated TE feature types
#                           (default: TE-CUT,TE-COPY)
#   --bins N                number of heatmap bins per axis (default: 50)
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
  stop("Usage: Rscript plot_gene_feature_distances.R <gff_dir> <output_dir> ",
       "[--gene-type TYPE] [--te-types TYPE,...] [--bins N]")
}

gff_dir <- args[[1L]]
output_dir <- args[[2L]]
args <- args[-c(1L, 2L)]

flag <- take_flag("--gene-type", args, "gene")
gene_type <- flag$value
args <- flag$args
flag <- take_flag("--te-types", args, "TE-CUT,TE-COPY")
te_types <- strsplit(flag$value, ",", fixed = TRUE)[[1L]]
args <- flag$args
flag <- take_flag("--bins", args, "50")
n_bins <- as.integer(flag$value)
args <- flag$args
flag <- take_bool("--final-generation", args)
final_generation_only <- flag$value
args <- flag$args

if (length(args) > 0L || !is.finite(n_bins) || n_bins < 2L) {
  stop("Unknown argument or invalid --bins value.")
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

extract_feature_distance <- function(genes, features, label) {
  if (nrow(features) == 0L) {
    genes[[paste0("upstream_", label)]] <- NA_real_
    genes[[paste0("downstream_", label)]] <- NA_real_
    return(genes)
  }

  setkey(features, seqid, start, end)
  upstream <- features[genes, on = .(seqid, end < start), mult = "last",
                       .(gene_row = i.gene_row, feature_end = x.end), nomatch = NA]
  downstream <- features[genes, on = .(seqid, start > end), mult = "first",
                         .(gene_row = i.gene_row, feature_start = x.start), nomatch = NA]

  genes[[paste0("upstream_", label)]] <- genes$start - upstream$feature_end
  genes[[paste0("downstream_", label)]] <- downstream$feature_start - genes$end
  genes[[paste0("upstream_", label)]][genes[[paste0("upstream_", label)]] < 0] <- 0
  genes[[paste0("downstream_", label)]][genes[[paste0("downstream_", label)]] < 0] <- 0
  genes
}

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
  # the same attributes:feature_type value. Collapse those rows into one span.
  grouped_rows <- features[
    feature_type %in% c("exon", "intron") & !is.na(attribute_feature_type)
  ]
  if (nrow(grouped_rows) > 0L) {
    genes <- grouped_rows[
      , .(
        start = min(start),
        end = max(end),
        gene_feature_type = first(attribute_feature_id)
      ),
      by = .(seqid, gene_group = attribute_feature_id)
    ]
  } else {
    genes <- features[feature_type == gene_type,
                      .(seqid, start, end, gene_feature_id = feature_id)]
  }
  genes[, gene_row := .I]
  genes
}

make_distance_plot <- function(data, distance_label, file_label, median_x, median_y) {
  data <- data[is.finite(log_upstream) & is.finite(log_downstream)]
  if (nrow(data) == 0L) {
    return(ggplot() + theme_void() +
             labs(title = paste("No", distance_label, "distances available for", file_label)))
  }

  x_limits <- range(data$log_downstream, finite = TRUE)
  y_limits <- range(data$log_upstream, finite = TRUE)
  x_pad <- max(diff(x_limits) * 0.04, 0.1)
  y_pad <- max(diff(y_limits) * 0.04, 0.1)
  x_limits <- x_limits + c(-x_pad, x_pad)
  y_limits <- y_limits + c(-y_pad, y_pad)

  ggplot(data, aes(x = log_downstream, y = log_upstream)) +
    geom_bin2d(bins = n_bins) +
    geom_vline(xintercept = median_x, linetype = "dashed", colour = "black") +
    geom_hline(yintercept = median_y, linetype = "dashed", colour = "black") +
    annotate("text", x = x_limits[1] + x_pad, y = y_limits[2] - y_pad,
         label = "QSL", hjust = 0, vjust = 1,
         colour = "black", fontface = "bold") +
    annotate("text", x = x_limits[2] - x_pad, y = y_limits[2] - y_pad,
         label = "QLL", hjust = 1, vjust = 1,
         colour = "black", fontface = "bold") +
    annotate("text", x = x_limits[1] + x_pad, y = y_limits[1] + y_pad,
         label = "QSS", hjust = 0, vjust = 0,
         colour = "black", fontface = "bold") +
    annotate("text", x = x_limits[2] - x_pad, y = y_limits[1] + y_pad,
         label = "QLS", hjust = 1, vjust = 0,
         colour = "black", fontface = "bold") +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits) +
    scale_fill_viridis_c(name = "Genes (N)", trans = "sqrt") +
    labs(
      title = paste(distance_label, "distance from each gene"),
      subtitle = sprintf("n = %d genes; dashed lines are median log10 distances", nrow(data)),
      x = "log10(downstream distance, bp)",
      y = "log10(upstream distance, bp)"
    ) +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
}

process_gff <- function(path) {
  features <- read_gff(path)
  if (is.null(features)) {
    warning("Skipping empty GFF: ", path)
    return(NULL)
  }

  genes <- build_genes(features)
  if (nrow(genes) == 0L) {
    warning("Skipping GFF with no grouped exon/intron or ", gene_type,
            " features: ", path)
    return(NULL)
  }
  gene_features <- genes[, .(seqid, start, end)]
  te_features <- features[feature_type %in% te_types, .(seqid, start, end)]

  genes <- extract_feature_distance(genes, gene_features, "gene")
  genes <- extract_feature_distance(genes, te_features, "te")
  genes[, `:=`(
    log_upstream_gene = log10(upstream_gene + 1),
    log_downstream_gene = log10(downstream_gene + 1),
    log_upstream_te = log10(upstream_te + 1),
    log_downstream_te = log10(downstream_te + 1)
  )]

  assign_quadrants <- function(upstream, downstream) {
    upstream_median <- median(upstream, na.rm = TRUE)
    downstream_median <- median(downstream, na.rm = TRUE)
    result <- rep(NA_character_, length(upstream))
    result[downstream < downstream_median & upstream > upstream_median] <- "Q1"
    result[downstream > downstream_median & upstream > upstream_median] <- "Q2"
    result[downstream > downstream_median & upstream < upstream_median] <- "Q3"
    result[downstream < downstream_median & upstream < upstream_median] <- "Q4"
    result
  }
  genes[, `:=`(
    quadrant_gene = assign_quadrants(log_upstream_gene, log_downstream_gene),
    quadrant_te = assign_quadrants(log_upstream_te, log_downstream_te)
  )]

  plot_data <- rbind(
    genes[, .(gene_row, seqid, start, end, distance_type = "next gene",
              log_upstream = log_upstream_gene, log_downstream = log_downstream_gene)],
    genes[, .(gene_row, seqid, start, end, distance_type = "next TE",
              log_upstream = log_upstream_te, log_downstream = log_downstream_te)]
  )
  plot_data <- plot_data[is.finite(log_upstream) & is.finite(log_downstream)]
  if (nrow(plot_data) == 0L) stop("No complete gene distance pairs found in ", path)

  base_name <- sub("\\.gff(\\.gz)?$", "", basename(path), ignore.case = TRUE)
  fwrite(genes, file.path(output_dir, paste0(base_name, "_gene_distances.csv")))

  gene_data <- plot_data[distance_type == "next gene"]
  te_data <- plot_data[distance_type == "next TE"]
  pdf(file.path(output_dir, paste0(base_name, "_gene_distance_quadrants.pdf")),
      width = 13, height = 6.5, onefile = TRUE)
  print(make_distance_plot(gene_data, "Next-gene", base_name,
                 median(gene_data$log_downstream), median(gene_data$log_upstream)))
  print(make_distance_plot(te_data, "Next-TE", base_name,
                 median(te_data$log_downstream), median(te_data$log_upstream)))
  dev.off()
  message("Wrote plots for ", basename(path))
  genes[, genome_file := basename(path)]
  genes
}

all_genes <- Filter(Negate(is.null), lapply(gff_files, process_gff))
if (length(all_genes) > 0L) {
  combined <- rbindlist(all_genes, use.names = TRUE, fill = TRUE)
  fwrite(combined, file.path(output_dir, "all_genomes_gene_distances.csv"))

  combined_plot_data <- rbind(
    combined[, .(genome_file, gene_row, seqid, start, end,
                 distance_type = "next gene",
                 log_upstream = log_upstream_gene,
                 log_downstream = log_downstream_gene)],
    combined[, .(genome_file, gene_row, seqid, start, end,
                 distance_type = "next TE",
                 log_upstream = log_upstream_te,
                 log_downstream = log_downstream_te)]
  )
  combined_plot_data <- combined_plot_data[
    is.finite(log_upstream) & is.finite(log_downstream)
  ]

  if (nrow(combined_plot_data) > 0L) {
    combined_gene_data <- combined_plot_data[distance_type == "next gene"]
    combined_te_data <- combined_plot_data[distance_type == "next TE"]
    pdf(file.path(output_dir, "all_genomes_gene_distance_quadrants.pdf"),
        width = 13, height = 6.5, onefile = TRUE)
    print(make_distance_plot(
      combined_gene_data, "Next-gene", "all genomes",
      median(combined_gene_data$log_downstream),
      median(combined_gene_data$log_upstream)
    ))
    print(make_distance_plot(
      combined_te_data, "Next-TE", "all genomes",
      median(combined_te_data$log_downstream),
      median(combined_te_data$log_upstream)
    ))
    dev.off()
    message("Wrote combined plots for ", length(all_genes), " genome(s)")
  }
}

