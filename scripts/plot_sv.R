#!/usr/bin/env Rscript
# sv_plot.R — Progressive Mauve-style SV plot for PansimNuc
#
# Each genome is drawn as a single horizontal track. Contigs within each
# genome are concatenated in numeric order with a small gap between them,
# mirroring the Progressive Mauve display style. Homologous elements
# (shared element_id) are connected by synteny ribbons across genomes.
#
# Usage:
#   Rscript sv_plot.R <root.gff> <sim.directory> [options]
#
# Options:
#   --out FILE       output file (default: sv_plot.pdf)
#   --width N        plot width in inches (default: 16)
#   --height N       plot height in inches (default: n_genomes x 2.5)
#   --types T,...    comma-separated feature types to display (default: all)
#   --link-types T,. comma-separated feature types to draw links for (default: all)
#   --no-links       suppress synteny ribbons
#   --gap N          bp gap inserted between contigs (default: 500)

suppressPackageStartupMessages({
  library(gggenomes)
  library(ggnewscale)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(tools)
  library(stringi)
  library(data.table)
})

# ── helpers ───────────────────────────────────────────────────────────────────

extract_attr <- function(x, key) {
  out <- sub(
    paste0(".*(?:^|;)", key, "=([^;]+).*"),
    "\\1",
    x,
    perl=TRUE
  )
  out[out == x] <- NA_character_
  out
}

parse_gff <- function(path, bin_id) {
  ids <- regmatches(
    basename(path),
    regexec("pop_(\\d+)_gen_(\\d+)_genome_(\\d+)", basename(path))
  )[[1]]
  
  pop_id <- as.integer(ids[2])
  generation <- as.integer(ids[3])
  genome_id <- as.integer(ids[4])
  
  dt <- fread(
    path,
    sep="\t",
    header=FALSE,
    comment.char="#",
    showProgress=FALSE
  )
  
  if (nrow(dt) == 0)
    return(NULL)
  
  setnames(
    dt,
    c("contig_name","source","feature_type","start","end",
      "score","strand","phase","attributes")
  )
  
  dt[, element_id :=
       as.integer(extract_attr(attributes, "element_id"))]
  
  dt[, log_genome_selection_coefficient :=
       as.numeric(extract_attr(
         attributes,
         "log_genome_selection_coefficient"))]
  
  dt[, log_element_selection_coefficient :=
       as.numeric(extract_attr(
         attributes,
         "log_element_selection_coefficient"))]
  
  dt[, multiplier :=
       as.numeric(extract_attr(attributes, "multiplier"))]
  
  dt[, feature_id :=
       as.numeric(extract_attr(attributes, "feature_id"))]
  
  dt[, `:=`(
    pop_id = pop_id,
    generation = generation,
    genome_id = genome_id,
    bin_id = bin_id
  )]
  
  dt[, .(
    bin_id,
    pop_id,
    generation,
    genome_id,
    contig_name,
    start,
    end,
    strand,
    feature_type,
    element_id,
    feature_id,
    multiplier,
    log_element_selection_coefficient,
    log_genome_selection_coefficient
  )]
}

# ── CLI argument parsing ──────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  cat(
    "Usage: Rscript sv_plot.R <root.gff> <sim0.gff> [sim1.gff ...]\n",
    "  [--out sv_plot.pdf] [--width 16] [--height auto]\n",
    "  [--types exon,intron,intergenic,TE-CUT,TE-COPY]\n",
    "  [--link-types exon,intron] [--no-links] [--gap 500] [--max-alignments 20]\n"
  )
  quit(status = 1L)
}

take_flag <- function(flag, args, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(list(val = default, args = args))
  val  <- args[i + 1L]
  args <- args[-c(i, i + 1L)]
  list(val = val, args = args)
}
take_switch <- function(flag, args) {
  i <- match(flag, args)
  if (is.na(i)) return(list(val = FALSE, args = args))
  list(val = TRUE, args = args[-i])
}

r <- take_flag("--out",        args, "sv_plot.pdf"); out_file   <- r$val; args <- r$args
r <- take_flag("--width",      args, NULL);          p_width    <- if (is.null(r$val)) NULL else as.numeric(r$val); args <- r$args
r <- take_flag("--height",     args, NULL);          p_height   <- if (is.null(r$val)) NULL else as.numeric(r$val); args <- r$args
r <- take_flag("--types",      args, NULL);          keep_types <- if (is.null(r$val)) NULL else strsplit(r$val, ",")[[1L]]; args <- r$args
r <- take_flag("--link-types", args, NULL);          link_types <- if (is.null(r$val)) NULL else strsplit(r$val, ",")[[1L]]; args <- r$args
r <- take_flag("--gap",        args, "500");         contig_gap <- as.integer(r$val); args <- r$args
r <- take_switch("--no-links", args);                no_links   <- r$val; args <- r$args
r <- take_flag("--max-alignments", args, "20");   max_alignments <- as.integer(r$val); args <- r$args
r <- take_switch("--final-generation", args);      final_generation_only <- r$val; args <- r$args

root_path <- args[1L]
sim_directory <- args[-1L]

# ── read data ─────────────────────────────────────────────────────────────────
sim_paths <- list.files(sim_directory, pattern = "*.gff", full.names = TRUE)
sim_paths <- sim_paths[sim_paths != root_path]

if (final_generation_only) {
  generations <- suppressWarnings(
    as.integer(sub(".*_gen_(\\d+)_genome_.*", "\\1", basename(sim_paths)))
  )
  if (!all(is.na(generations))) {
    last_generation <- max(generations, na.rm = TRUE)
    sim_paths <- sim_paths[!is.na(generations) & generations == last_generation]
    message("Restricting SV plot input to generation ", last_generation)
  }
}

# randomly downsample to max_alignments if there are more than max_alignments
if (length(sim_paths) > max_alignments) {
  sim_paths <- sample(root_path, max_alignments)
}

sim_paths <- c(root_path, sim_paths)

message("Reading all GFFs...")
all_feats <- rbindlist(
  mapply(
    parse_gff,
    path = sim_paths,
    bin_id = tools::file_path_sans_ext(basename(sim_paths)),
    SIMPLIFY = FALSE
  ),
  use.names = TRUE,
  fill = TRUE
)

if (is.null(all_feats) || nrow(all_feats) == 0L) {
  stop("No features loaded. Check that the GFF files are valid PansimNuc output.")
}

if (!is.null(keep_types)) {
  all_feats <- filter(all_feats, feature_type %in% keep_types)
}

# ── linearize contigs (Progressive Mauve style) ──────────────────────────────
# Within each genome, sort contigs by their numeric suffix (contig_1 < contig_2
# < ...), concatenate them end-to-end with `contig_gap` bp between each pair,
# and shift all feature coordinates into that linearized space. Each genome
# becomes a single seq_id so gggenomes draws it as one track.

# Extract contig number once
all_feats[, contig_num :=
            as.integer(stringi::stri_extract_first_regex(contig_name, "\\d+"))]


contig_info <- all_feats[
  ,
  .(contig_len = max(end)),
  by = .(bin_id, contig_name, contig_num)
]

contig_gap <- max(contig_info$contig_len) * 0.1

# Calculate offsets
contig_info[
  ,
  `:=`(
    offset = cumsum(data.table::shift(contig_len + contig_gap, fill = 0)),
    contig_rank = seq_len(.N)
  ),
  by = bin_id
]


# ─────────────────────────────────────────────────────────────────────────────
# Contig blocks
# ─────────────────────────────────────────────────────────────────────────────

contig_blocks <- contig_info[
  ,
  .(
    bin_id,
    seq_id = bin_id,
    start = offset + 1L,
    end = offset + contig_len,
    contig_name,
    contig_shade = factor((contig_rank - 1L) %% 2L)
  )
]


# ─────────────────────────────────────────────────────────────────────────────
# Shift feature coordinates
# ─────────────────────────────────────────────────────────────────────────────

setkey(contig_info, bin_id, contig_name)

all_feats[
  contig_info,
  offset := i.offset,
  on = .(bin_id, contig_name)
]

all_feats[
  ,
  `:=`(
    seq_id = bin_id,
    start = start + offset,
    end = end + offset
  )
]


# ─────────────────────────────────────────────────────────────────────────────
# Sequence table
# ─────────────────────────────────────────────────────────────────────────────

seqs <- all_feats[
  ,
  .(
    seq_id = first(bin_id),
    length = max(end)
  ),
  by = bin_id
]


seqs[
  ,
  sort_key := fifelse(
    bin_id == "root",
    -1L,
    as.integer(stringi::stri_extract_first_regex(bin_id, "\\d+"))
  )
]

setorder(seqs, sort_key)
seqs[, sort_key := NULL]


# ─────────────────────────────────────────────────────────────────────────────
# Links
# ─────────────────────────────────────────────────────────────────────────────

links <- NULL


if (!no_links) {
  
  ordered_genomes <- seqs$seq_id
  
  # Split once instead of repeatedly filtering
  feat_split <- split(all_feats, by = "seq_id")
  
  
  links_list <- vector(
    "list",
    length(ordered_genomes) - 1L
  )
  
  
  for (i in seq_len(length(ordered_genomes) - 1L)) {
    
    upper <- feat_split[[ordered_genomes[i]]][
      ,
      .(
        seq_id,
        start,
        end,
        strand1 = strand,
        element_id,
        feature_type
      )
    ]
    
    
    if (!is.null(link_types)) {
      upper <- upper[
        feature_type %chin% link_types
      ]
    }
    
    
    lower <- feat_split[[ordered_genomes[i + 1L]]][
      ,
      .(
        seq_id2 = seq_id,
        start2 = start,
        end2 = end,
        strand2 = strand,
        element_id,
        feature_type2 = feature_type
      )
    ]
    
    
    links_list[[i]] <- merge(
      upper,
      lower,
      by = "element_id",
      allow.cartesian = TRUE
    )[
      ,
      .(
        element_id,
        seq_id,
        start,
        end,
        strand1,
        seq_id2,
        start2,
        end2,
        strand2,
        feature_type,
        feature_type2
      )
    ]
  }
  
  
  links <- rbindlist(
    links_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  
  if (nrow(links) > 0) {
    links[
      ,
      strand := fifelse(
        strand1 == strand2,
        "+",
        "-"
      )
    ]
  } else {
    links <- NULL
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# Colours
# ─────────────────────────────────────────────────────────────────────────────

feature_colors <- c(
  exon       = "#4DAF4A",
  intron     = "#984EA3",
  intergenic = "#999999",
  "TE-CUT"   = "#E41A1C",
  "TE-COPY"  = "#FF7F00"
)


extra_types <- setdiff(
  all_feats[, unique(feature_type)],
  names(feature_colors)
)


if (length(extra_types) > 0) {
  feature_colors <- c(
    feature_colors,
    setNames(
      hue_pal()(length(extra_types)),
      extra_types
    )
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# Plot dimensions
# ─────────────────────────────────────────────────────────────────────────────

n_bins <- nrow(seqs)

max_seq_len <- max(seqs$length)


if (is.null(p_width))
  p_width <- min(max(max_seq_len / 2.5e5, 7), 49.9)

if (is.null(p_height))
  p_height <- min(max(4.0, n_bins * 0.3), 49.9)


# ─────────────────────────────────────────────────────────────────────────────
# Build gggenomes plot
# ─────────────────────────────────────────────────────────────────────────────

if (!no_links &&
    !is.null(links) &&
    nrow(links) > 0) {
  
  p <- gggenomes(
    seqs = as.data.frame(seqs),
    genes = as.data.frame(all_feats),
    links = as.data.frame(links),
    feats = list(
      contigs = as.data.frame(contig_blocks)
    )
  )
  
} else {
  
  p <- gggenomes(
    seqs = as.data.frame(seqs),
    genes = as.data.frame(all_feats),
    feats = list(
      contigs = as.data.frame(contig_blocks)
    )
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# Plot layers
# ─────────────────────────────────────────────────────────────────────────────

p <- p +
  geom_feat(
    data = feats("contigs"),
    aes(colour = contig_shade),
    alpha = 0.12,
    linewidth = NA,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "0" = "#AAAAAA",
      "1" = "#555555"
    ),
    guide = "none"
  ) +
  geom_seq() +
  geom_seq_label()


if (!no_links &&
    !is.null(links) &&
    nrow(links) > 0) {
  
  p <- p +
    geom_link(
      aes(fill = feature_type),
      alpha = 0.28,
      colour = NA
    ) +
    scale_fill_manual(
      values = feature_colors,
      na.value = "grey60",
      name = "Feature type"
    ) +
    new_scale_fill()
}


p <- p +
  geom_gene(
    aes(fill = feature_type)
  ) +
  scale_fill_manual(
    values = feature_colors,
    na.value = "grey60",
    name = "Feature type"
  ) +
  theme_gggenomes_clean()


# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────

message(
  sprintf(
    "Writing %s (%.1f x %.1f in)",
    out_file,
    p_width,
    p_height
  )
)

ggsave(
  out_file,
  p,
  width = p_width,
  height = p_height
)

message("Done.")
