library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
args <- args[!grepl("^--", args)]

input_dir <- if (length(args) >= 1) args[1] else "."
outpref <- if (length(args) >= 2) args[2] else "final_DFE_plot"

DFE_files <- list.files(
    input_dir,
    pattern="^pop_\\d+_final_dfe.csv",
    full.names=TRUE
)

if (length(DFE_files) == 0)
  stop("No final DFE files found.")

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

parse_file <- function(path) {
  ids <- regmatches(
    basename(path),
    regexec("pop_(\\d+)_final_dfe.csv", basename(path))
  )[[1]]
  
  pop_id <- as.integer(ids[2])
  pop_id <- paste0("Population ", pop_id)
  
  dt <- fread(
    path,
    sep=",",
    header=TRUE,
    comment.char="#",
    showProgress=FALSE
  )
  
  if (nrow(dt) == 0)
    return(NULL)
  
  dt$selection_coefficient <- dt$selection_coefficient + 1.0
  
  dt[, `:=`(
    pop_id = pop_id
  )]
}

all_data <- rbindlist(
  lapply(DFE_files, parse_file),
  use.names=TRUE,
  fill=TRUE
)

if (nrow(all_data) == 0)
  stop("No features found.")  

p.split <- ggplot(all_data, aes(x = selection_coefficient, after_stat(ncount), fill = feature_type)) +
  geom_histogram() +
  facet_grid(pop_id ~ feature_type) +
  scale_fill_npg() +
  scale_x_continuous(limits = c(0, NA)) +
  geom_vline(xintercept = 1.0, colour = "black", linetype="dotted") +
  labs(
    x = "Selection Coefficient", y = "Scaled Count", fill = "Feature Type") +
  theme_light()

ggsave(file.path(paste0(outpref, "_", "split.png")),
       p.split, width = 10, height = 6)

p.total <- ggplot(all_data, aes(x = selection_coefficient, after_stat(ncount))) +
  geom_histogram() +
  facet_grid(pop_id ~ .) +
  scale_fill_npg() +
  geom_vline(xintercept = 1.0, colour = "black", linetype="dotted") +
  labs(
    x = "Selection Coefficient", y = "Scaled Count", fill = "Feature Type") +
  theme_light()

ggsave(file.path(paste0(outpref, "_", "total.png")),
       p.total, width = 5, height = 6)
