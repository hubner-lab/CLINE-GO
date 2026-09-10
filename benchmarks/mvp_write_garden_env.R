#!/usr/bin/env Rscript
# mvp_write_garden_env.R -- write ALL garden environment tables for a seed up front.
#
# The per-garden driver (mvp_garden_run_par.sh) rewrote a single
# {PROJECT}_env_garden.tsv before each Snakemake invocation, so a seed's 42
# gardens could only be run one Snakemake call at a time. With the {scenario}
# dimension in place a garden is just one more scenario, so every garden gets
# its own table here and one invocation covers all of them.
#
# A garden is a UNIFORM transplant environment: every source site is evaluated
# at the garden's (opt1, opt2), which is why each output table repeats one pair
# of values down all rows -- the site column carries the source identity.
#
# Scenario name == file basename, so the tables are named after the garden id
# and the pipeline's scenario wildcard reads back as `deme_001`, `novel_0.25`, ...
#
# Usage:
#   Rscript benchmarks/mvp_write_garden_env.R --seed=1232218 [--outdir=DIR] [--gardens=TSV]
#
# Writes: {INDIR}/gardens/{garden_id}.tsv  (INDIR = data/mvp/MVP{seed})

suppressPackageStartupMessages(library(data.table))

ROOT <- Sys.getenv("PIPELINE_ROOT", "/mnt/data/eugene/CLINE-GO")

args <- commandArgs(trailingOnly = TRUE)
kv <- function(k, d = NULL) {
    hit <- grep(paste0("^--", k, "="), args, value = TRUE)
    if (!length(hit)) return(d)
    sub(paste0("^--", k, "="), "", hit[1])
}

SEED <- kv("seed")
if (is.null(SEED)) stop("--seed=<seed> is required")

PROJ    <- paste0("MVP", SEED)
INDIR   <- file.path(ROOT, "data/mvp", PROJ)
GARDENS <- kv("gardens", file.path(ROOT, "benchmarks/mvp_eval/offset09",
                                   paste0("gardens_", SEED, ".tsv")))
OUTDIR  <- kv("outdir", file.path(INDIR, "gardens"))
PRESENT <- file.path(INDIR, paste0(PROJ, "_env_present.tsv"))

for (f in c(GARDENS, PRESENT)) {
    if (!file.exists(f)) stop("missing input: ", f)
}
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

gard <- fread(GARDENS, colClasses = c("garden_id" = "character"))
pres <- fread(PRESENT, colClasses = c("site" = "character"))

# Predictor columns are whatever the present table carries besides the key.
predictors <- setdiff(names(pres), "site")
if (length(predictors) != 2) {
    stop("expected 2 predictor columns in ", PRESENT, "; found: ",
         paste(predictors, collapse = ", "))
}
if (!all(c("opt1", "opt2") %in% names(gard))) {
    stop("gardens table must carry opt1/opt2 columns: ", GARDENS)
}

message(sprintf("== %s: %d gardens x %d sites ==", PROJ, nrow(gard), nrow(pres)))

for (i in seq_len(nrow(gard))) {
    g   <- gard$garden_id[i]
    out <- data.table(site = pres$site)
    out[[predictors[1]]] <- gard$opt1[i]
    out[[predictors[2]]] <- gard$opt2[i]
    fwrite(out, file.path(OUTDIR, paste0(g, ".tsv")), sep = "\t")
}

message(sprintf("wrote %d tables to %s", nrow(gard), OUTDIR))
message("scenario names: ", paste(head(gard$garden_id, 3), collapse = ", "), ", ...")
