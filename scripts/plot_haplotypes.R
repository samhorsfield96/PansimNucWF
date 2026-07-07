library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(ggpattern)
library(vcfR)
library(tools)

# Usage:
#   Rscript plot_haplotypes.R [vcf_file] [output_prefix] [top_n] [recombination_threshold]
#
# Arguments:
#   vcf_file                – path to VCF/VCF.gz file
#   output_prefix           – prefix for output files       (default: haplotypes)
#   top_n                   - Number of top-changing haplotypes to return. 0 keeps all.
#   recombination_threshold - proportion of mutations in a haplotype that must be present to
#                             identify a new haplotype as being a recombinant. Default 0.9 (90%).
#                             Set to 0 to disable recombinant detection and classification.
#
# Sample names must follow: pop_<N>_gen_<N>_genome_<N>
# If they do not, all samples are assigned population_id=0 and generation=0.
#
# Output: multiple PDFs and a CSV describing haplotype dynamics across generations.

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[!grepl("^--", args)]
vcf_file <- if (length(args) >= 1) args[1] else stop("VCF file path required as first argument")
gff_dir <- if (length(args) >= 2) args[2] else "."
outpref  <- if (length(args) >= 3) args[3] else "haplotypes"
top_n    <- if (length(args) >= 4) as.integer(args[4]) else 5L
recombination_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 0.9

# ── Parse sample names → population_id / generation / genome_id ──────────────

parse_sample_name <- function(nm) {
  m <- regmatches(nm, regexpr("^pop_(\\d+)_gen_(\\d+)_genome_(\\d+)$", nm, perl = TRUE))
  if (length(m) == 0L) return(NULL)
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
  gff_files <- list.files(gff_dir, pattern = "\\.gff$", full.names = TRUE)
  sel_coeffs <- data.frame(chrom = character(), pos = integer(), sel_coeff = numeric(), stringsAsFactors = FALSE)
  rows <- lapply(gff_files, function(gff_file) {
    base <- file_path_sans_ext(basename(gff_file))
    parsed <- parse_sample_name(base)
    if (is.null(parsed)) {
      # skip files with unparseable names, but warn if any found
      return(NULL)
    }
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
        stringsAsFactors = FALSE
      )
    }
  })
  bind_rows(Filter(Negate(is.null), rows))
}

# ── Read VCF and build per-genome mutation profiles ───────────────────────────

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
                           function(nm) is.null(parse_sample_name(nm)), logical(1L)))
  if (n_unparsed > 0L) {
    message(n_unparsed, " sample name(s) could not be parsed as pop_N_gen_N_genome_N; ",
            "assigning population_id=0, generation=0.")
  }

  sample_meta <- do.call(rbind, lapply(seq_along(sample_names), function(i) {
    nm      <- sample_names[i]
    parsed  <- parse_sample_name(nm)
    toks    <- mut_tokens[has_alt[, i]]
    mut_sig <- if (length(toks) == 0L) NA_character_ else paste(sort(toks), collapse = ";")
    if (is.null(parsed)) {
      data.frame(sample = nm, population_id = 0L, generation = 0L,
                 genome_id = i, mut_sig = mut_sig, log_sel_coeff = 0,
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
}

if (nrow(sample_meta) == 0L) stop("No samples found in VCF.")
message(sprintf("Loaded profiles for %d sample(s) across %d generation(s).",
                nrow(sample_meta), length(unique(sample_meta$generation))))

# ── Helper functions ──────────────────────────────────────────────────────────
# checked, all good
parse_sig <- function(s) if (nchar(s) == 0) character(0) else strsplit(s, ";")[[1]]

# ── Whole-genome haplotype functions ─────────────────────────────────────────
# profile_c           : parsed character vector of mutation tokens for one genome
# all_profiles_named  : named list (profile_str -> parsed vec) of all known profiles
# threshold           : proportion of mutations in a haplotype that must be present to 
#                       identify a new haplotype as being a recombinant. Default 0.9 (90%). 
#                       Set to 0 to disable recombinant detection and classification.
# Returns a length-2 character vector of the two parent profile strings when C
# is a recombinant (union of their mutation sets equals C's, each contributes
# at least one exclusive mutation), or NULL otherwise.
# checked, all good
find_recombinant_parents <- function(profile_c, all_profiles_named, threshold = 0.9) {
  if (threshold <= 0 || length(all_profiles_named) < 2) return(NULL)

  profile_c_len <- length(profile_c)
  candidate_names <- Filter(
    function(nm) {
      a <- all_profiles_named[[nm]]
      length(a) > 0 &&
        (sum(a %in% profile_c) / length(a)) >= threshold
    },
    names(all_profiles_named)
  )
  if (length(candidate_names) < 2) return(NULL)

  n <- length(candidate_names)
  for (i in seq_len(n - 1)) {
    a_name <- candidate_names[[i]]
    a      <- all_profiles_named[[a_name]]
    for (j in seq(i + 1, n)) {
      b_name <- candidate_names[[j]]
      b      <- all_profiles_named[[b_name]]
      if (!any(!a %in% b) || !any(!b %in% a)) next
      # A and B each contribute >= threshold proportion of their mutations to C
      return(c(a_name, b_name))
    }
  }
  NULL
}

# Assign a globally-consistent base label (e.g. "M1", "R2", "REF") to every
# mutation profile, processing all populations together in generation order.
# The same mutation profile gets the same number across all populations, so
# migration is immediately visible: P0M3 and P1M3 carry identical mutations.
build_global_hap_labels <- function(sample_meta, recombination_threshold) {
  all_genome_profiles <- sample_meta %>%
    group_by(population_id, genome_id, generation) %>%
    summarise(
      profile_str = {
        sigs <- mut_sig[!is.na(mut_sig)]
        if (length(sigs) == 0L) "NA" else paste(sigs, collapse = ";")
      },
      .groups = "drop"
    )

  all_generations <- sort(unique(all_genome_profiles$generation))
  first_gen       <- all_generations[1L]
  adjustment      <- if (length(all_generations) > 1L) 1L else 0L

  known_profiles <- list()
  hap_labels     <- list()
  counters       <- list()

  new_label <- function(prefix) {
    counters[[prefix]] <<- if (is.null(counters[[prefix]])) 1L else counters[[prefix]] + 1L
    paste0(prefix, counters[[prefix]])
  }

  for (gen in all_generations) {
    prev_profiles      <- unique(all_genome_profiles$profile_str[
      all_genome_profiles$generation == (gen - adjustment)])
    prev_profiles      <- prev_profiles[!is.na(prev_profiles)]
    all_profiles_named <- setNames(lapply(prev_profiles, parse_sig), prev_profiles)

    curr_profiles <- unique(all_genome_profiles$profile_str[
      all_genome_profiles$generation == gen])

    for (prof_str in curr_profiles) {
      if (prof_str %in% names(known_profiles)) next
      prof_vec <- parse_sig(prof_str)
      if (prof_str == "NA") {
        lbl <- "REF"
      } else if (gen == first_gen) {
        lbl <- new_label("M")
      } else {
        parents <- find_recombinant_parents(prof_vec, all_profiles_named, recombination_threshold)
        lbl <- if (!is.null(parents)) new_label("R") else new_label("M")
      }
      known_profiles[[prof_str]] <- prof_vec
      hap_labels[[prof_str]]     <- lbl
    }
  }
  hap_labels
}

# Classify whole-genome haplotypes for one population.
# pop_df must have mut_sig (from assign_element_sigs) and log_sel_coeff columns.
# Returns: generation, haplotype_id, profile_str, sequence (concatenated), freq, type, sel_coeff
# checked, all good
classify_genome_haplotypes <- function(pop_df, population_id = 0L, global_labels) {
  generations <- sort(unique(pop_df$generation))

  genome_profiles <- pop_df %>%
    group_by(genome_id, generation) %>%
    summarise(
      profile_str = {
        sigs <- mut_sig[!is.na(mut_sig)]
        if (length(sigs) == 0L) "NA" else paste(sigs, collapse = ";")
      },
      sel_coeff   = sum(log_sel_coeff, na.rm = TRUE),
      .groups     = "drop"
    )
  
  known_profiles <- list()   # profile_str -> parsed vec  (profiles seen so far)
  known_types    <- list()   # profile_str -> haplotype type string
  known_parents  <- list()   # profile_str -> "P1,P2" or NA
  hap_labels     <- list()   # profile_str -> full label (with population prefix)
  
  # determine how many generations present, adjust which generation to look for recombinants
  if (length(generations) > 1)
  {
    adjustment = 1
  } else {
    adjustment = 0
  }

  rows <- list()
  for (gen in generations) {
    gen_data <- genome_profiles[genome_profiles$generation == gen, ]
    n_total  <- nrow(gen_data)
    if (n_total == 0) next

    # look for recombinants in the context of profiles only in prior generation
    all_profile_names  <- names(table(genome_profiles[genome_profiles$generation == (gen - adjustment), ]$profile_str))
    all_profiles_named <- setNames(
      lapply(all_profile_names, parse_sig),
      all_profile_names
    )

    prof_tbl    <- table(gen_data$profile_str)
    sel_by_prof <- split(gen_data$sel_coeff, gen_data$profile_str)

    for (prof_str in names(prof_tbl)) {
      freq      <- prof_tbl[[prof_str]] / n_total
      sel_coeff <- mean(unlist(sel_by_prof[[prof_str]]), na.rm = TRUE)

      if (!prof_str %in% names(known_profiles)) {
        prof_vec       <- parse_sig(prof_str)
        parent_str     <- NA_character_
        if (prof_str == "NA") {
          htype <- "reference"
        } else {
          parent_profiles <- find_recombinant_parents(prof_vec, all_profiles_named, threshold = recombination_threshold)
          if (!is.null(parent_profiles)) {
            parent_str <- paste(parent_profiles, collapse = "||")
            htype <- if (gen == generations[1L]) "founder" else "recombinant"
          } else {
            htype <- if (gen == generations[1L]) "founder" else "mutant"
          }
        }
        base_label <- global_labels[[prof_str]]
        if (is.null(base_label)) base_label <- paste0("?", substr(prof_str, 1L, 8L))
        known_profiles[[prof_str]] <- prof_vec
        known_types[[prof_str]]    <- htype
        known_parents[[prof_str]]  <- parent_str
        hap_labels[[prof_str]]     <- paste0("P", population_id, base_label)
      }

      rows[[length(rows) + 1]] <- data.frame(
        generation        = gen,
        haplotype_id      = hap_labels[[prof_str]],
        base_haplotype_id = global_labels[[prof_str]],
        profile_str       = prof_str,
        freq              = freq,
        type              = known_types[[prof_str]],
        parents           = known_parents[[prof_str]],
        sel_coeff         = sel_coeff,
        stringsAsFactors  = FALSE
      )
    }
  }
  result <- bind_rows(rows)

  # Resolve parent profile strings to haplotype labels now that all labels are assigned.
  result$parents <- vapply(result$parents, function(ps) {
    if (is.na(ps) || nchar(ps) == 0) return(NA_character_)

    profs  <- strsplit(ps, "\\|\\|")[[1]]
    labels <- vapply(profs, function(p) {
      lbl <- hap_labels[[p]]
      if (is.null(lbl)) NA_character_ else lbl
    }, character(1L))
    paste(labels[!is.na(labels)], collapse = ",")
  }, character(1L))

  result
}

# ── Classify whole-genome haplotypes ─────────────────────────────────────────
message("Classifying whole-genome haplotypes per population...")

hap_data_rds_file <- paste0(outpref, "_hap_data.rds")
if (!file.exists(hap_data_rds_file)) { 
  global_labels <- build_global_hap_labels(sample_meta, recombination_threshold)
  hap_data <- sample_meta %>%
    group_by(population_id) %>%
    group_modify(~ classify_genome_haplotypes(.x, .y$population_id, global_labels)) %>%
    ungroup()
  
  # Fill in zero-frequency rows so that every haplotype appears in every
  # generation (needed for correct stacked areas and unbroken lines).
  hap_data <- hap_data %>%
    group_by(population_id) %>%
    complete(
      generation   = unique(generation),
      haplotype_id = unique(haplotype_id),
      fill         = list(freq = 0)
    ) %>%
    group_by(population_id, haplotype_id) %>%
    fill(type, base_haplotype_id, profile_str, sel_coeff, .direction = "downup") %>%
    ungroup()
  
  saveRDS(hap_data, hap_data_rds_file)
} else {
  hap_data <- readRDS(hap_data_rds_file)
}

# ── Ensure migration columns exist ────────────────────────────────────────────
if (!"source_population_id" %in% colnames(hap_data)) hap_data$source_population_id <- NA_integer_
if (!"source_haplotype_id"  %in% colnames(hap_data)) hap_data$source_haplotype_id  <- NA_character_

# ── Detect cross-population haplotype migration ───────────────────────────────
# A haplotype is classified as a migrant when the same mutation profile
# (profile_str) is found in another population at a strictly later generation.
# The source population is the one where the profile first appeared.
if (length(unique(hap_data$population_id)) > 1) {
  message("Detecting cross-population haplotype migration...")

  # First generation each non-reference profile appears (freq > 0) per population
  profile_origins <- hap_data %>%
    filter(freq > 0, !is.na(profile_str), profile_str != "NA") %>%
    group_by(population_id, profile_str) %>%
    summarise(
      first_gen    = min(generation),
      haplotype_id = haplotype_id[which.min(generation)],
      .groups      = "drop"
    )

  # Profiles shared across more than one population
  multi_pop_profiles <- profile_origins %>%
    group_by(profile_str) %>%
    filter(n() > 1) %>%
    ungroup()

  if (nrow(multi_pop_profiles) > 0) {
    # Source = earliest first appearance; ties broken by lowest population_id
    source_origins <- multi_pop_profiles %>%
      group_by(profile_str) %>%
      arrange(first_gen, population_id) %>%
      slice(1L) %>%
      ungroup() %>%
      rename(source_population_id = population_id,
             source_haplotype_id  = haplotype_id,
             source_first_gen     = first_gen)

    # Founding generation per population (earliest generation present)
    founding_gens <- hap_data %>%
      group_by(population_id) %>%
      summarise(founding_gen = min(generation), .groups = "drop")

    # Migrants: same profile appeared strictly later in a different population,
    # BUT only if it did not arrive in the founding generation of that population
    # (individuals present at founding are classified as founders, not migrants).
    migrant_entries <- multi_pop_profiles %>%
      inner_join(source_origins, by = "profile_str") %>%
      inner_join(founding_gens, by = "population_id") %>%
      filter(population_id != source_population_id,
             first_gen > source_first_gen,
             first_gen > founding_gen) %>%
      select(-founding_gen)

    if (nrow(migrant_entries) > 0) {
      for (i in seq_len(nrow(migrant_entries))) {
        prof    <- migrant_entries$profile_str[i]
        pop     <- migrant_entries$population_id[i]
        src_pop <- migrant_entries$source_population_id[i]
        src_hap <- migrant_entries$source_haplotype_id[i]
        mask    <- hap_data$profile_str == prof & hap_data$population_id == pop
        # haplotype_id is kept as-is: the shared global label number already
        # links the migrant to its source (e.g. P0M3 migrated to P1M3).
        hap_data$type[mask]                 <- "migrant"
        hap_data$source_population_id[mask] <- src_pop
        hap_data$source_haplotype_id[mask]  <- src_hap
      }
      message(sprintf(
        "Identified %d migrant haplotype(s) from %d source population(s).",
        nrow(migrant_entries),
        length(unique(migrant_entries$source_population_id))
      ))
    } else {
      message("No migration events detected (shared profiles appeared simultaneously).")
    }
  } else {
    message("No haplotype profiles shared across populations.")
  }
}

# ── Summary table ─────────────────────────────────────────────────────────────
hap_summary <- hap_data %>%
  group_by(population_id, haplotype_id, type, profile_str) %>%
  summarise(
    first_generation     = min(generation[freq > 0]),
    peak_freq            = max(freq),
    mean_sel_coeff       = mean(sel_coeff, na.rm = TRUE),
    source_population_id = { v <- source_population_id[!is.na(source_population_id)]; if (length(v)) v[1L] else NA_integer_ },
    source_haplotype_id  = { v <- source_haplotype_id[!is.na(source_haplotype_id)];  if (length(v)) v[1L] else NA_character_ },
    .groups              = "drop"
  ) %>%
  arrange(population_id, first_generation)

write.csv(hap_summary,
          file      = paste0(outpref, "_haplotype_summary.csv"),
          row.names = FALSE)

# ── Filter to top N haplotypes per type (0 = keep all) ───────────────────────
if (top_n > 0L) {
  message(sprintf("Retaining top %d haplotype(s) per type by peak frequency.", top_n))
  top_haps <- hap_data %>%
    arrange(population_id, haplotype_id, generation) %>%
    group_by(population_id, haplotype_id, type) %>%
    summarise(max_val = max(freq), .groups = "drop") %>%
    group_by(population_id, type) %>%
    slice_max(max_val, n = top_n, with_ties = FALSE) %>%
    ungroup()

  hap_data <- hap_data %>%
    semi_join(top_haps, by = c("population_id", "haplotype_id"))
}

# Order haplotype_id factor levels by type so same-coloured areas stack
# contiguously (avoids interleaved colour breaks in stacked area charts).
# Population-prefixed labels (e.g. P0M3, P1M3) are unique, so distinct() is
# sufficient — no deduplication needed.
type_stack_order <- c("reference", "founder", "mutant", "recombinant", "migrant")
hap_id_levels <- hap_data %>%
  distinct(haplotype_id, type) %>%
  mutate(type = factor(type, levels = type_stack_order)) %>%
  arrange(type, haplotype_id) %>%
  pull(haplotype_id)
hap_data <- hap_data %>%
  mutate(haplotype_id = factor(haplotype_id, levels = hap_id_levels))

# ── Plotting helpers ──────────────────────────────────────────────────────────
n_pops        <- length(unique(hap_data$population_id))
has_multi_pop <- n_pops > 1

type_colour_values <- c(
  reference   = "#3C5488FF",
  founder     = "#4DBBD5",
  mutant      = "#E64B35",
  recombinant = "#00A087",
  migrant     = "#F39B7FFF"
)

type_colour_scale <- scale_colour_manual(values = type_colour_values, name = "Haplotype")
type_fill_scale   <- scale_fill_manual(values = type_colour_values,   name = "Haplotype")

add_facets <- function(p) {
  if (has_multi_pop) {
    p + facet_grid(
      ~ population_id,
      labeller = as_labeller(function(x) paste0("population_id: ", x))
    )
  } else {
    p
  }
}

base_theme <- theme_light(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

# ── Plot 1: haplotype frequency lines, coloured by type ──────────────────────
message("Plotting haplotype frequency lines...")

p_lines <- ggplot(
  hap_data,
  aes(
    x      = generation,
    y      = freq,
    colour = type,
    group  = interaction(haplotype_id, population_id)
  )
) +
  geom_line(alpha = 0.8) +
  labs(x = "Generation", y = "Haplotype frequency", colour = "Haplotype") +
  scale_y_continuous(limits = c(0, 1)) +
  type_colour_scale +
  base_theme

p_lines <- add_facets(p_lines)
ggsave(paste0(outpref, "_haplotype_freq.pdf"), plot = p_lines, width = 8, height = 6)

# ── Plot 2: stacked area chart of haplotype composition ──────────────────────
message("Plotting stacked haplotype composition areas...")

p_area <- ggplot(
  hap_data,
  aes(
    x    = generation,
    y    = freq,
    fill = type,
    group = haplotype_id
  )
) +
  geom_area(position = "stack", colour = NA, alpha = 0.8) +
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Haplotype") +
  scale_y_continuous(limits = c(0, 1)) +
  type_fill_scale +
  base_theme

p_area <- add_facets(p_area)
ggsave(paste0(outpref, "_haplotype_composition.pdf"), plot = p_area, width = 8, height = 6)

# ── Plot 3: stacked area chart of top changing haplotype composition ──────────────────────
message("Plotting changing haplotype composition stacked areas...")

p_area <- ggplot(
  hap_data,
  aes(
    x    = generation,
    y    = freq,
    fill = base_haplotype_id,
    group = haplotype_id
  )
) +
  geom_area(position = "stack", colour = NA, alpha = 0.8) +
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Haplotype ID") +
  scale_y_continuous(limits = c(0, 1)) +
  base_theme

p_area <- add_facets(p_area)
ggsave(paste0(outpref, "_per_haplotype_composition.pdf"), plot = p_area, width = 8, height = 6)

# ── Plot 4: top hits with selection coefficients + haplotype-type hatching ───
type_pattern_values <- c(
  reference   = "none",
  founder     = "none",
  mutant      = "stripe",
  recombinant = "crosshatch",
  migrant     = "circle"
)

p_sel <- ggplot(
  hap_data,
  aes(
    x              = generation,
    y              = freq,
    fill           = sel_coeff,
    colour           = sel_coeff,
    pattern        = type,
    pattern_colour = type,
    group          = haplotype_id
  )
) +
  geom_area_pattern(
    position        = "stack",
    colour          = NA,
    alpha           = 0.8,
    pattern_density = 0.35,
    pattern_spacing = 0.025,
    pattern_fill    = NA
  ) +
  scale_pattern_manual(
    values = type_pattern_values,
    name   = "Haplotype type"
  ) +
  scale_pattern_colour_manual(
    values = c(
      reference   = "grey30",
      founder     = "grey30",
      mutant      = "grey30",
      recombinant = "grey30",
      migrant     = "grey30"
    ),
    name = "Haplotype type"
  ) +
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Selection coefficient") +
  scale_y_continuous(limits = c(0, 1)) +
  base_theme

p_sel <- add_facets(p_sel)
ggsave(paste0(outpref, "_sel_coeff_composition.pdf"), plot = p_sel, width = 8, height = 6)

message(sprintf(
  "Done. %d haplotypes tracked (%d founder, %d mutant, %d recombinant, %d migrant).",
  nrow(hap_summary),
  sum(hap_summary$type == "founder"),
  sum(hap_summary$type == "mutant"),
  sum(hap_summary$type == "recombinant"),
  sum(hap_summary$type == "migrant")
))

