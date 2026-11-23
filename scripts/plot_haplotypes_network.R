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
#   Rscript plot_haplotypes_network.R [vcf_file] [gff_dir] [output_prefix] [top_n] [recombination_threshold]
#
# Arguments:
#   vcf_file                – path to VCF/VCF.gz file
#   gff_dir                 – directory containing GFF files
#   output_prefix           – prefix for output files       (default: haplotypes)
#   top_n                   - Number of top-changing haplotypes to return. 0 keeps all.
#   recombination_threshold - accepted for compatibility; network classification uses
#                             shared alleles rather than this threshold.
#
# Sample names must follow: pop_<N>_gen_<N>_genome_<N>
# If they do not, all samples are assigned population_id=0 and generation=0.
#
# Output: a haplotype summary CSV, four PDFs, an igraph RDS and an interactive HTML network.

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[!grepl("^--", args)]
vcf_file <- if (length(args) >= 1) args[1] else stop("VCF file path required as first argument")
gff_dir <- if (length(args) >= 2) args[2] else "."
outpref  <- if (length(args) >= 3) args[3] else "haplotypes"
top_n    <- if (length(args) >= 4) as.integer(args[4]) else 5L
recombination_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 1.0

vcf_file <- "/Users/samhorsfield/Software/PansimNuc_analysis/laplace_selection_mean0_scale0.1_w_recombination_mig_0.1_N50_Gen50_mu1-e6/PansimNucWF/variants/filtered_variants_refgt.vcf.gz"
gff_dir <-  "/Users/samhorsfield/Software/PansimNuc_analysis/laplace_selection_mean0_scale0.1_w_recombination_mig_0.1_N50_Gen50_mu1-e6"

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
  reference_coeffs <- sample_meta$log_sel_coeff[sample_meta$generation == 0]
  if (length(reference_coeffs) > 0L) reference_log_sel_coeff <- mean(reference_coeffs)
}

if (nrow(sample_meta) == 0L) stop("No samples found in VCF.")
message(sprintf("Loaded profiles for %d sample(s) across %d generation(s).",
                nrow(sample_meta), length(unique(sample_meta$generation))))

# The reference genome's log selection coefficient is the fixed baseline.
sample_meta$sel_coeff <- exp(sample_meta$log_sel_coeff - reference_log_sel_coeff)

# Edges describe shared observed alternate alleles, not confirmed parentage.
# Exact copies retain their predecessor's type. A single-parent mutant can add
# alleles absent from the previous generation in the same population. Other
# types count incoming candidate predecessors, not confirmed parents.
normalise_profile <- function(sig) {
  if (is.na(sig) || !nzchar(sig)) return("")
  paste(sort(unique(strsplit(sig, ";", fixed = TRUE)[[1L]])), collapse = ";")
}

#TODO - issue here is that I am trying to reinvent the wheel. Perhaps there is some tool e.g. plink that can help with identifying haplotype networks, from this I can infer ancestry and recombination?

build_allele_network <- function(sample_meta) {
  profiles <- vapply(sample_meta$mut_sig, normalise_profile, character(1L))

  nodes <- sample_meta %>%
    mutate(profile = profiles) %>%
    group_by(population_id, generation, profile) %>%
    summarise(
      count = n(),
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

          # ensure strict gain of mutations, meaning child must be same length or more than parent
          if (length(child_snps) < length(mutations[[parent]])) {
            next
          }

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

  # Profile origins and edges are inferred from observations, not a pedigree.
  origins <- nodes %>%
    arrange(generation, population_id) %>%
    group_by(profile) %>%
    slice_head(n = 1L) %>%
    ungroup()
  origin_index <- match(nodes$profile, origins$profile)
  source_populations <- origins$population_id[origin_index]
  first_generations <- ave(nodes$generation, nodes$population_id, FUN = min)
  types <- rep("unresolved", nrow(nodes))
  migrant_lineage <- rep(FALSE, nrow(nodes))
  parent_nodes <- rep("", nrow(nodes))
  source_nodes <- rep("", nrow(nodes))
  generation_key <- paste(nodes$population_id, nodes$generation, sep = ":")
  generation_nodes <- split(seq_len(nrow(nodes)), generation_key)
  generation_alleles <- lapply(generation_nodes, function(indices) {
    unique(unlist(mutations[indices], use.names = FALSE))
  })
  for (i in order(nodes$generation, nodes$population_id, nodes$name)) {
    current_haplotype <- nodes$profile[i]
    parents <- as.integer(igraph::neighbors(graph, i, mode = "in"))
    previous_key <- paste(nodes$population_id[i], nodes$generation[i] - 1L, sep = ":")
    previous <- generation_nodes[[previous_key]]
    identical_parent <- previous[nodes$profile[previous] == current_haplotype]
    mutation_parent <- integer()
    if (length(identical_parent) == 0L && length(parents) > 0L) {
      previous_alleles <- generation_alleles[[previous_key]]
      inherited <- intersect(mutations[[i]], previous_alleles)
      novel <- setdiff(mutations[[i]], previous_alleles)
      if (length(novel) > 0L) {
        mutation_parent <- parents[vapply(parents, function(parent) {
          setequal(mutations[[parent]], inherited)
        }, logical(1L))]
      }
    }
    selected_parents <- if (length(identical_parent) > 0L) {
      identical_parent
    } else if (length(mutation_parent) == 1L) {
      mutation_parent
    } else {
      parents
    }
    parent_lineages <- migrant_lineage[selected_parents]
    is_migrant <- nodes$population_id[i] != source_populations[i]
    migrant_lineage[i] <- is_migrant || any(parent_lineages)
    if (is_migrant) {
      source_nodes[i] <- origins$name[origin_index[i]]
    } else if (length(identical_parent) > 0L) {
      parent_nodes[i] <- parent_nodes[identical_parent[1L]]
      source_nodes[i] <- source_nodes[identical_parent[1L]]
    } else {
      parent_nodes[i] <- paste(nodes$name[selected_parents], collapse = ";")
      inherited_sources <- unlist(strsplit(source_nodes[selected_parents], ";", fixed = TRUE),
                                  use.names = FALSE)
      source_nodes[i] <- paste(unique(inherited_sources[nzchar(inherited_sources)]),
                               collapse = ";")
    }

    types[i] <- if (is_migrant) {
      "migrant"
    } else if (length(identical_parent) > 0L) {
      types[identical_parent[1L]]
    } else if (length(selected_parents) >= 2L) {
      if (any(parent_lineages) && any(!parent_lineages)) "migrant-recombinant" else "recombinant"
    } else if (length(selected_parents) == 1L) {
      if (parent_lineages[1L]) "migrant-de-novo" else "mutant"
    } else if (nodes$generation[i] == first_generations[i]) {
      "founder"
    } else if (length(selected_parents) == 0L) {
      "mutant"
    }
    else {
      "unresolved"
    }
  }
  igraph::V(graph)$type <- types
  igraph::V(graph)$origin_population_id <- source_populations
  igraph::V(graph)$migrant_lineage <- migrant_lineage
  igraph::V(graph)$parent_nodes <- parent_nodes
  igraph::V(graph)$source_nodes <- source_nodes
  graph
}

allele_network <- build_allele_network(sample_meta)
saveRDS(allele_network, paste0(outpref, "_allele_network.rds"))

network_haplotypes <- igraph::as_data_frame(allele_network, what = "vertices")
global_profiles <- network_haplotypes %>%
  arrange(generation, population_id, profile) %>%
  distinct(profile, .keep_all = TRUE) %>%
  mutate(
    prefix = case_when(
      profile == "" ~ "REF",
      grepl("recombinant", type, fixed = TRUE) ~ "R",
      TRUE ~ "M"
    ),
    number = ave(seq_along(prefix), prefix, FUN = seq_along),
    base_haplotype_id = ifelse(prefix == "REF", prefix, paste0(prefix, number))
  ) %>%
  select(profile, base_haplotype_id)
network_haplotypes <- network_haplotypes %>%
  left_join(global_profiles, by = "profile") %>%
  mutate(haplotype_id = paste0("P", population_id, base_haplotype_id))

node_labels <- setNames(network_haplotypes$haplotype_id, network_haplotypes$name)
node_populations <- setNames(as.character(network_haplotypes$population_id), network_haplotypes$name)
resolve_node_refs <- function(refs, values) {
  if (!nzchar(refs)) return(NA_character_)
  ids <- strsplit(refs, ";", fixed = TRUE)[[1L]]
  result <- unname(values[ids])
  if (anyNA(result)) stop("Haplotype metadata refers to unknown network nodes.")
  paste(result, collapse = ",")
}

selection_by_profile <- sample_meta %>%
  mutate(profile = vapply(mut_sig, normalise_profile, character(1L))) %>%
  group_by(population_id, generation, profile) %>%
  summarise(sel_coeff = mean(sel_coeff, na.rm = TRUE), .groups = "drop")

hap_data <- network_haplotypes %>%
  left_join(selection_by_profile, by = c("population_id", "generation", "profile")) %>%
  group_by(population_id, generation) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup() %>%
  mutate(
    profile_str = ifelse(profile == "", "NA", profile),
    parents = vapply(parent_nodes, resolve_node_refs, character(1L), values = node_labels),
    source_population_id = vapply(source_nodes, resolve_node_refs, character(1L),
                                  values = node_populations),
    source_haplotype_id = vapply(source_nodes, resolve_node_refs, character(1L),
                                 values = node_labels)
  ) %>%
  select(population_id, generation, haplotype_id, base_haplotype_id, profile_str,
         freq, type, parents, sel_coeff, source_population_id, source_haplotype_id)

hap_summary <- hap_data %>%
  group_by(population_id, haplotype_id, type, profile_str) %>%
  summarise(
    first_generation = min(generation),
    peak_freq = max(freq),
    mean_sel_coeff = mean(sel_coeff, na.rm = TRUE),
    source_population_id = { v <- source_population_id[!is.na(source_population_id)]; if (length(v)) v[1L] else NA_character_ },
    source_haplotype_id = { v <- source_haplotype_id[!is.na(source_haplotype_id)]; if (length(v)) v[1L] else NA_character_ },
    parents = { v <- parents[!is.na(parents)]; if (length(v)) v[1L] else NA_character_ },
    .groups = "drop"
  ) %>%
  arrange(population_id, first_generation)
write.csv(hap_summary, paste0(outpref, "_haplotype_summary.csv"), row.names = FALSE)

hap_data <- hap_data %>%
  group_by(population_id) %>%
  complete(generation = unique(generation), haplotype_id = unique(haplotype_id),
           fill = list(freq = 0)) %>%
  group_by(population_id, haplotype_id) %>%
  fill(type, base_haplotype_id, profile_str, sel_coeff, .direction = "downup") %>%
  ungroup()

if (top_n > 0L) {
  top_haps <- hap_data %>%
    arrange(population_id, haplotype_id, generation) %>%
    group_by(population_id, haplotype_id, type) %>%
    summarise(max_val = max(freq), .groups = "drop") %>%
    group_by(population_id, type) %>%
    slice_max(max_val, n = top_n, with_ties = FALSE) %>%
    ungroup()
  hap_data <- hap_data %>% semi_join(top_haps, by = c("population_id", "haplotype_id"))
}

type_stack_order <- c("reference", "founder", "mutant", "recombinant", "migrant",
                      "migrant-recombinant", "migrant-de-novo", "migrant-recombinant-de-novo",
                      "unresolved")
hap_id_levels <- hap_data %>%
  distinct(haplotype_id, type) %>%
  mutate(type = factor(type, levels = type_stack_order)) %>%
  arrange(type, haplotype_id) %>%
  distinct(haplotype_id) %>%
  pull(haplotype_id)
hap_data <- hap_data %>% mutate(haplotype_id = factor(haplotype_id, levels = hap_id_levels))

type_colour_values <- c(
  reference = "#3C5488FF", founder = "#4DBBD5", mutant = "#E64B35",
  recombinant = "#00A087", migrant = "#F39B7FFF",
  "migrant-recombinant" = "#7E6148FF", "migrant-de-novo" = "#8491B4FF",
  "migrant-recombinant-de-novo" = "#B09C85FF", unresolved = "grey70"
)
type_pattern_values <- c(
  reference = "none", founder = "none", mutant = "stripe", recombinant = "crosshatch",
  migrant = "circle", "migrant-recombinant" = "weave", "migrant-de-novo" = "wave",
  "migrant-recombinant-de-novo" = "stripe", unresolved = "none"
)
add_facets <- function(p) {
  if (length(unique(hap_data$population_id)) > 1L) {
    p + facet_grid(~ population_id, labeller = as_labeller(function(x) paste0("population_id: ", x)))
  } else {
    p
  }
}
base_theme <- theme_light(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

p_lines <- ggplot(hap_data, aes(x = generation, y = freq, colour = type,
                                group = interaction(haplotype_id, population_id))) +
  geom_line(alpha = 0.8) +
  labs(x = "Generation", y = "Haplotype frequency", colour = "Haplotype") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = type_colour_values, name = "Haplotype") + base_theme
ggsave(paste0(outpref, "_haplotype_freq.pdf"), plot = add_facets(p_lines), width = 8, height = 6)

p_area <- ggplot(hap_data, aes(x = generation, y = freq, fill = type, group = haplotype_id)) +
  geom_area(position = "stack", colour = NA, alpha = 0.8) +
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Haplotype") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = type_colour_values, name = "Haplotype") + base_theme
ggsave(paste0(outpref, "_haplotype_composition.pdf"), plot = add_facets(p_area), width = 8, height = 6)

p_hap <- ggplot(hap_data, aes(x = generation, y = freq, fill = base_haplotype_id,
                              group = haplotype_id)) +
  geom_area(position = "stack", colour = NA, alpha = 0.8) +
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Haplotype ID") +
  scale_y_continuous(limits = c(0, 1)) + base_theme
ggsave(paste0(outpref, "_per_haplotype_composition.pdf"), plot = add_facets(p_hap), width = 8, height = 6)

p_sel <- ggplot(hap_data, aes(x = generation, y = freq, fill = sel_coeff, colour = sel_coeff,
                              pattern = type, pattern_colour = type, group = haplotype_id)) +
  geom_area_pattern(position = "stack", colour = NA, alpha = 0.8, pattern_density = 0.35,
                    pattern_spacing = 0.025, pattern_fill = NA) +
  scale_pattern_manual(values = type_pattern_values, name = "Haplotype type") +
  scale_pattern_colour_manual(values = setNames(rep("grey30", length(type_pattern_values)),
                                                names(type_pattern_values)),
                              name = "Haplotype type") +
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Selection coefficient") +
  scale_y_continuous(limits = c(0, 1)) + base_theme
ggsave(paste0(outpref, "_sel_coeff_composition.pdf"), plot = add_facets(p_sel), width = 8, height = 6)

#TESTING
network_nodes <- igraph::as_data_frame(allele_network, what = "vertices")
network_nodes$plot_row <- ave(seq_len(nrow(network_nodes)), network_nodes$generation,
                              FUN = seq_along)
network_layout <- cbind(network_nodes$generation, network_nodes$plot_row)
type_colours <- c(founder = "grey70", mutant = "steelblue",
                  recombinant = "tomato", unresolved = "orange",
                  migrant = "purple", `migrant-de-novo` = "darkgreen",
                  `migrant-recombinant` = "brown")
 
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
    "\nOrigin population: ", nodes$origin_population_id,
    "\nMigrant lineage: ", nodes$migrant_lineage,
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
