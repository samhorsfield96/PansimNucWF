library(data.table)
library(ggplot2)
library(ggsci)

args <- commandArgs(trailingOnly = TRUE)
args <- args[!grepl("^--", args)]
final_generation_only <- any(commandArgs(trailingOnly = TRUE) == "--final-generation")

output_dir <- if (length(args) >= 1) args[1] else "."
outpref <- if (length(args) >= 2) args[2] else "te_copy_numbers"

if (!dir.exists(output_dir))
    stop("Output directory does not exist: ", output_dir)

gff_files <- list.files(
    output_dir,
    pattern="^pop_\\d+_gen_\\d+_genome_\\d+\\.gff$",
    full.names=TRUE
)

if (length(gff_files) == 0)
    stop("No GFF files found.")

if (final_generation_only) {
    gens <- as.integer(sub(".*_gen_(\\d+)_genome_.*", "\\1",
                           basename(gff_files)))
    last_gen <- max(gens, na.rm=TRUE)
    gff_files <- gff_files[gens == last_gen]
    message("Restricting to generation ", last_gen)
}

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

parse_gff <- function(path) {

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
        c("seqid","source","feature_type","start","end",
          "score","strand","phase","attributes")
    )

    dt <- dt[
        feature_type %chin% c("TE-CUT","TE-COPY")
    ]

    if (nrow(dt) == 0)
        return(NULL)

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

    dt[, `:=`(
        pop_id = pop_id,
        generation = generation,
        genome_id = genome_id
    )]

    dt[, .(
        pop_id,
        generation,
        genome_id,
        feature_type,
        element_id,
        log_genome_selection_coefficient,
        log_element_selection_coefficient,
        multiplier
    )]
}

message("Parsing GFFs...")

all_data <- rbindlist(
    lapply(gff_files, parse_gff),
    use.names=TRUE,
    fill=TRUE
)

if (nrow(all_data) == 0)
    stop("No TE features found.")

copy_counts <- all_data[
    ,
    .(
        copies=.N,
        log_genome_selection_coefficient=
            mean(log_genome_selection_coefficient),
        log_element_selection_coefficient=
            mean(log_element_selection_coefficient),
        multiplier=mean(multiplier)
    ),
    by=.(pop_id,
         generation,
         genome_id,
         feature_type,
         element_id)
]

copy_dist <- copy_counts[
    ,
    .(frequency=.N),
    by=.(pop_id,
         generation,
         feature_type,
         copies)
]

mean_copies <- copy_counts[
    ,
    .(
        mean_copies=mean(copies),
        sd_copies=sd(copies),
        n_genomes=.N
    ),
    by=.(pop_id,
         generation,
         feature_type,
         element_id)
]

total_load <- copy_counts[
    ,
    .(total_copies=sum(copies)),
    by=.(pop_id,
         generation,
         genome_id,
         feature_type)
]

total_load_summary <- total_load[
    ,
    .(
        mean_load=mean(total_copies),
        sd_load=sd(total_copies),
        median_load=median(total_copies)
    ),
    by=.(pop_id,
         generation,
         feature_type)
]

fwrite(copy_counts,
       paste0(outpref, "_per_genome.csv"))

fwrite(copy_dist,
       paste0(outpref, "_distribution.csv"))

fwrite(mean_copies,
       paste0(outpref, "_mean_per_element.csv"))

fwrite(total_load_summary,
       paste0(outpref, "_total_load.csv"))