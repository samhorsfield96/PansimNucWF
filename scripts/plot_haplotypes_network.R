library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(ggpattern)
library(vcfR)
library(tools)
library(networkD3)
library(htmlwidgets)

# Usage:
#   Rscript plot_haplotypes_new.R [vcf_file] [gff_dir] [output_prefix] [top_n] [recombination_threshold]
#
# Arguments:
#   vcf_file                – path to VCF/VCF.gz file
#   gff_dir                 – directory containing GFF files
#   output_prefix           – prefix for output files       (default: haplotypes)
#   top_n                   - Number of top-changing haplotypes to return. 0 keeps all.
#   recombination_threshold - proportion of mutations in a haplotype that must be present to
#                             identify a new haplotype as being a recombinant. Default 0.9 (90%).
#                             Set to 0 to disable recombinant detection and classification.
#
# Sample names must follow: pop_<N>_gen_<N>_genome_<N>
# If they do not, all samples are assigned population_id=0 and generation=0.
#
# Output: an igraph RDS and an interactive HTML allele-sharing network.

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[!grepl("^--", args)]
vcf_file <- if (length(args) >= 1) args[1] else stop("VCF file path required as first argument")
gff_dir <- if (length(args) >= 2) args[2] else "."
outpref  <- if (length(args) >= 3) args[3] else "haplotypes"
top_n    <- if (length(args) >= 4) as.integer(args[4]) else 5L
recombination_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 1.0

vcf_file <- "/Users/samhorsfield/Software/PansimNuc_analysis/pansimnuc_runs/neutral_w_recombination_mig_0.01_N50_Gen50/PansimNucWF/variants/filtered_variants_refgt.vcf.gz"
gff_dir <-  "/Users/samhorsfield/Software/PansimNuc_analysis/pansimnuc_runs/neutral_w_recombination_mig_0.01_N50_Gen50"

# ── Parse sample names → population_id / generation / genome_id ──────────────

parse_sample_name <- function(nm) {
  m <- regmatches(nm, regexpr("^pop_(\\d+)_gen_(\\d+)_genome_(\\d+)$", nm, perl = TRUE))
  if (length(m) == 0L) return(list(population_id = NA, generation = NA, genome_id = NA))
  parts <- as.integer(strsplit(sub("^pop_", "", m), "_gen_|_genome_")[[1L]])
  list(population_id = parts[1L], generation = parts[2L], genome_id = parts[3L])
}

parse_attrs <- function(attr_str) {
  pairs <- strsplit(attr_str, ";", fixed = TRUE)[[1L]]
  kv    <- strsplit(pairs, "=", fixed = TRUE)
  keys  <- vapply(kv, `[[`, character(1L), 1L)
  vals  <- vapply(kv, function(x) if (length(x) >= 2L) x[[2L]] else NA_character_,
                  character(1L))
  setNames(vals, keys)
}

# read the start of gff files to get selection coefficients for mutations, if available
get_gff_data <- function(gff_dir) {
  gff_files <- list.files(gff_dir, pattern = "\\.gff(\\.gz)?$", full.names = TRUE)
  sel_coeffs <- data.frame(chrom = character(), pos = integer(), sel_coeff = numeric(), stringsAsFactors = FALSE)
  rows <- lapply(gff_files, function(gff_file) {
    # remove all extensions including compression
    base <- file_path_sans_ext(file_path_sans_ext(basename(gff_file)))
    parsed <- parse_sample_name(base)
    con = file(gff_file, "r")
    while(TRUE) {
      line = readLines(con, n = 1)
      if (!startsWith(line, "#")) {
        break
      }
    }
    close(con) 
    
    if (length(line) > 0) {
      f <- strsplit(line, "\t")[[1]]
      a <- parse_attrs(f[9L])
      contig_name  <- f[1L]
      contig_index <- suppressWarnings(
        as.integer(sub("contig_", "", contig_name)) - 1L
      )
      data.frame(
        sample_name = base,
        population = parsed$population_id,
        generation = parsed$generation,
        genome = parsed$genome_id,
        contig_index = contig_index,
        start        = as.integer(f[4L]),   # GFF is 1-based
        end          = as.integer(f[5L]),
        strand       = f[7L],
        element_id   = suppressWarnings(as.integer(a[["element_id"]])),
        feature_type = a[["feature_type"]],
        log_sel_coeff = suppressWarnings(as.numeric(a[["log_genome_selection_coefficient"]])),
        log_sel_prob = suppressWarnings(as.numeric(a[["log_genome_selection_probability"]])),
        stringsAsFactors = FALSE
      )
    }
  })
  bind_rows(Filter(Negate(is.null), rows))
}

# ── Read VCF and build per-genome mutation profiles ───────────────────────────

# make a note of the reference log sel coefficient
reference_log_sel_coeff <- 0.0

genome_profiles_rds <- paste0(outpref, "_genome_profiles.rds")
if (!file.exists(genome_profiles_rds)) {
  message("Reading VCF: ", vcf_file)
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  gff_data <- get_gff_data(gff_dir)
  
  gt_mat <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  # Guard against single-variant VCFs which drop the matrix dimension
  if (!is.matrix(gt_mat)) {
    gt_mat <- matrix(gt_mat, nrow = 1L, dimnames = list(NULL, names(gt_mat)))
  }
  
  chrom      <- getCHROM(vcf)
  pos        <- getPOS(vcf)
  alt        <- getALT(vcf)
  mut_tokens <- paste(chrom, pos, alt, sep = ":")
  
  # Any non-ref, non-missing call is treated as carrying the alt allele
  ref_gts <- c("0", "0/0", "0|0")
  has_alt  <- !is.na(gt_mat) & !gt_mat %in% ref_gts
  
  sample_names <- colnames(gt_mat)
  
  n_unparsed <- sum(vapply(sample_names,
                           function(nm) is.na(parse_sample_name(nm)[1]), logical(1L)))
  if (n_unparsed > 0L) {
    message(n_unparsed, " sample name(s) could not be parsed as pop_N_gen_N_genome_N; ",
            "assigning population_id=0, generation=0.")
  }
  
  sample_meta <- do.call(rbind, lapply(seq_along(sample_names), function(i) {
    nm      <- sample_names[i]
    parsed  <- parse_sample_name(nm)
    toks    <- mut_tokens[has_alt[, i]]
    mut_sig <- if (length(toks) == 0L) NA_character_ else paste(sort(toks), collapse = ";")
    if (is.na(parsed[1])) {
      #determine if reference, take this as the reference sel coefficient
      if (is.na(mut_sig))
      {
        log_sel_coeff <- gff_data$log_sel_coeff[gff_data$sample_name == nm]
        population_id = -1
        generation = -1
        reference_log_sel_coeff <<- if (length(log_sel_coeff) > 0) mean(log_sel_coeff) else 0.0
      } else
      {
        log_sel_coeff = 0.0
        population_id = 0L
        generation = 0L
      }
      
      # write output
      data.frame(sample = nm, population_id = 0L, generation = 0L,
                 genome_id = i, mut_sig = mut_sig, log_sel_coeff = if (length(log_sel_coeff) > 0) log_sel_coeff else 0.0,
                 stringsAsFactors = FALSE)
    } else {
      log_sel_coeff <- gff_data$log_sel_coeff[gff_data$sample_name == nm]
      data.frame(sample = nm, population_id = parsed$population_id,
                 generation = parsed$generation, genome_id = parsed$genome_id,
                 mut_sig = mut_sig, log_sel_coeff = if (length(log_sel_coeff) > 0) log_sel_coeff else 0.0,
                 stringsAsFactors = FALSE)
    }
  }))
  
  saveRDS(sample_meta, genome_profiles_rds)
} else {
  sample_meta <- readRDS(genome_profiles_rds)
  reference_coeffs <- sample_meta$log_sel_coeff[sample_meta$population_id == -1]
  if (length(reference_coeffs) > 0L) reference_log_sel_coeff <- mean(reference_coeffs)
}

if (nrow(sample_meta) == 0L) stop("No samples found in VCF.")
message(sprintf("Loaded profiles for %d sample(s) across %d generation(s).",
                nrow(sample_meta), length(unique(sample_meta$generation))))

# The reference genome's log selection coefficient is the fixed baseline.
sample_meta$sel_coeff <- exp(sample_meta$log_sel_coeff - reference_log_sel_coeff)

# Edges describe shared observed alternate alleles, not confirmed parentage.
build_allele_network <- function(sample_meta) {
  profiles <- vapply(sample_meta$mut_sig, function(sig) {
    if (is.na(sig) || !nzchar(sig)) return("")
    paste(sort(unique(strsplit(sig, ";", fixed = TRUE)[[1L]])), collapse = ";")
  }, character(1L))

  nodes <- sample_meta %>%
    mutate(profile = profiles) %>%
    group_by(population_id, generation, profile) %>%
    summarise(
      count = n(),
      #samples = paste(sort(sample), collapse = ","),
      genome_ids = paste(sort(genome_id), collapse = ","),
      .groups = "drop"
    ) %>%
    arrange(population_id, generation, profile) %>%
    mutate(name = paste0("node_", row_number())) %>%
    select(name, everything())
  mutations <- lapply(nodes$profile, function(sig) {
    if (!nzchar(sig)) character(0L) else strsplit(sig, ";", fixed = TRUE)[[1L]]
  })

  edge_list <- list()
  for (pop in unique(nodes$population_id)) {
    pop_nodes <- which(nodes$population_id == pop)
    for (gen in unique(nodes$generation[pop_nodes])) {
      previous <- pop_nodes[nodes$generation[pop_nodes] == gen - 1L]
      current <- pop_nodes[nodes$generation[pop_nodes] == gen]
      if (length(previous) == 0L) next

      previous_snps <- unlist(mutations[previous], use.names = FALSE)
      previous_index <- split(rep(previous, lengths(mutations[previous])), previous_snps)
      for (child in current) {
        child_snps <- mutations[[child]]
        candidate_indices <- unique(unlist(previous_index[mutations[[child]]], use.names = FALSE))
        for (parent in candidate_indices) {
          shared <- intersect(mutations[[parent]], mutations[[child]])
          edge_list[[length(edge_list) + 1L]] <- data.frame(
            from = nodes$name[parent],
            to = nodes$name[child],
            shared_alleles = paste(shared, collapse = ";"),
            n_shared_alleles = length(shared),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  edges <- if (length(edge_list) > 0L) {
    bind_rows(edge_list)
  } else {
    data.frame(from = character(), to = character(), shared_alleles = character(),
               n_shared_alleles = integer(), stringsAsFactors = FALSE)
  }
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  parent_counts <- igraph::degree(graph, mode = "in")
  first_generations <- ave(nodes$generation, nodes$population_id, FUN = min)
  igraph::V(graph)$type <- ifelse(
    parent_counts >= 2L, "recombinant",
    ifelse(parent_counts == 1L, "mutant",
           ifelse(nodes$generation == first_generations, "founder", "unresolved"))
  )
  graph
}

allele_network <- build_allele_network(sample_meta)
saveRDS(allele_network, paste0(outpref, "_allele_network.rds"))

network_nodes <- igraph::as_data_frame(allele_network, what = "vertices")
network_nodes$plot_row <- ave(seq_len(nrow(network_nodes)), network_nodes$generation,
                              FUN = seq_along)
network_layout <- cbind(network_nodes$generation, network_nodes$plot_row)
type_colours <- c(founder = "grey70", mutant = "steelblue",
                  recombinant = "tomato", unresolved = "orange") 
 
plot(
  allele_network, layout = network_layout,
  vertex.color = type_colours[network_nodes$type],
  vertex.size = 10 + 3 * sqrt(network_nodes$count),
  vertex.label = paste0("P", network_nodes$population_id, " n=", network_nodes$count),
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.3,
  main = "Shared alternate alleles between consecutive generations"
)

build_interactive_allele_network <- function(graph) {
  nodes <- igraph::as_data_frame(graph, what = "vertices")
  edges <- igraph::as_data_frame(graph, what = "edges")
  nodes$label <- paste0("P", nodes$population_id, " gen ", nodes$generation,
                        " n=", nodes$count)
  links <- data.frame(
    source = match(edges$from, nodes$name) - 1L,
    target = match(edges$to, nodes$name) - 1L,
    shared_count = edges$n_shared_alleles
  )
  if (anyNA(links$source) || anyNA(links$target)) {
    stop("Network edges refer to unknown nodes.")
  }

  # networkD3 warns when node zero has no edges, even though indices are valid.
  widget <- withCallingHandlers(
    networkD3::forceNetwork(
      Links = links, Nodes = nodes, Source = "source", Target = "target",
      Value = "shared_count", NodeID = "label", Nodesize = "count", Group = "type",
      linkColour = if (nrow(links) == 0L) character() else "#666666",
      arrows = TRUE, zoom = TRUE, legend = TRUE, opacity = 0.8
    ),
    warning = function(w) {
      if (identical(conditionMessage(w),
                    "It looks like Source/Target is not zero-indexed. This is required in JavaScript and so your plot may not render.")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  widget$x$nodes$details <- paste0(
    nodes$label, "\nType: ", nodes$type, "\nGenomes: ", nodes$genome_ids,
    "\nAlleles: ", ifelse(nzchar(nodes$profile), nodes$profile, "(none)")
  )
  widget$x$links$details <- if (nrow(edges) > 0L) {
    paste0("Shared alleles: ", edges$shared_alleles)
  } else {
    character()
  }
  htmlwidgets::onRender(widget, "
    function(el) {
      d3.select(el).selectAll('.node').append('title')
        .text(function(d) { return d.details; });
      d3.select(el).selectAll('.link').append('title')
        .text(function(d) { return d.details; });
    }
  ")
}

interactive_network <- build_interactive_allele_network(allele_network)
htmlwidgets::saveWidget(interactive_network, file = paste0(outpref, "_allele_network.html"),
                        selfcontained = FALSE)
