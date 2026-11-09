#!/usr/bin/env Rscript
# pregea_hit_overlap.R: do the rungs of the LFMM and EMMAX ladders call the same SNPs?
#
# The ladders hold counts, so two rungs with equal counts may call different SNPs. This counts
# Bonferroni hits per predictor as the production run calls them: LFMM divides each predictor's
# chi-squared by that predictor's own lambda, EMMAX applies no recalibration, and each predictor is
# its own Bonferroni family. A hit is one SNP on one predictor.
#
#   Rscript manuscript/arabidopsis/pregea_hit_overlap.R <project> [--lfmm-k FROM:TO]
#
# Reads <project>_results/_intermediate/pregea/{lfmm,emmax}/<rung>/pvalues.tsv and the two ladder
# tables. The rungs are every rung a ladder table lists; --lfmm-k keeps the LFMM rungs K FROM to TO,
# which must all be in it. Stops unless each rung's count per predictor equals the ladder's hits_bonf.
# Writes ANALYSIS/<project>/pregea/rung_hit_overlap.tsv, one row per rung:
#   engine rung reference n_hits n_snps n_shared_with_reference jaccard_with_reference
#
# The reference is K = k_best for LFMM and n_pcs = 0 for EMMAX. These are pruned markers, so the
# counts rank the rungs and are not a forecast of the whole-genome run.
#
# Run from the repo root, in the image:
#   docker run --rm --name arab-hit-overlap --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp \
#     -e OPENBLAS_NUM_THREADS=1 -v $PWD:/data -w /data cline-go:6cdedcc \
#     Rscript manuscript/arabidopsis/pregea_hit_overlap.R <project> [--lfmm-k FROM:TO]

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(yaml)
})

usage <- "usage: pregea_hit_overlap.R <project> [--lfmm-k FROM:TO]"
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 3L) || (length(args) == 3L && args[2] != "--lfmm-k")) stop(usage)
project <- args[1]
lfmm_k <- NULL
if (length(args) == 3L) {
  bounds <- suppressWarnings(as.integer(strsplit(args[3], ":", fixed = TRUE)[[1]]))
  if (length(bounds) != 2L || anyNA(bounds) || bounds[1] > bounds[2]) stop(usage)
  lfmm_k <- seq(bounds[1], bounds[2])
}
bonferroni_alpha <- 0.05

results <- here(paste0(project, "_results"))
config <- read_yaml(here(paste0("config_", project, ".yaml")))
out_file <- here("ANALYSIS", project, "pregea", "rung_hit_overlap.tsv")
key_columns <- c("SNPID", "chr", "pos")

recalibrate <- function(p) {
  chisq <- qchisq(p, df = 1, lower.tail = FALSE)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
  pchisq(chisq / lambda, df = 1, lower.tail = FALSE)
}

rung_hits <- function(engine, rung, apply_gc, ladder) {
  pvalues <- fread(file.path(results, "_intermediate", "pregea", engine, rung, "pvalues.tsv"))
  hits <- rbindlist(lapply(setdiff(names(pvalues), key_columns), function(predictor) {
    p <- pvalues[[predictor]]
    if (apply_gc) p <- recalibrate(p)
    called <- pvalues$SNPID[which(p <= bonferroni_alpha / sum(!is.na(p)))]
    expected <- ladder[rung_value == sub("^(K|npc)", "", rung) & trait == predictor, hits_bonf]
    if (length(called) != expected) stop(engine, " ", rung, " ", predictor, ": ", length(called),
                                         " hits, the ladder has ", expected)
    data.table(snp = called, predictor = predictor)
  }))
  hits[, rung := rung]
}

read_ladder <- function(engine) {
  fread(file.path(results, "PreGEA", "tables", engine, paste0(engine, "_ladder.tsv")),
        colClasses = c(rung_value = "character"))
}

# Every rung the ladder lists, in numeric order, or only `values` (all of which it must list).
ladder_rungs <- function(ladder, prefix, values = NULL) {
  listed <- sort(unique(ladder$rung_value_num))
  if (is.null(values)) values <- listed
  absent <- setdiff(values, listed)
  if (length(absent) > 0L) stop("the ", ladder$engine[1], " ladder has no rung ",
                                paste0(prefix, absent, collapse = ", "))
  paste0(prefix, values)
}

compare_rungs <- function(engine, rungs, reference, apply_gc, ladder) {
  if (!reference %in% rungs) stop(engine, ": the reference rung ", reference, " is not among ",
                                  paste(rungs, collapse = ", "))
  hits <- rbindlist(lapply(rungs, function(r) rung_hits(engine, r, apply_gc, ladder)))
  hits[, key := paste(snp, predictor)]
  reference_keys <- hits[rung == reference, key]
  rbindlist(lapply(rungs, function(r) {
    keys <- hits[rung == r, key]
    shared <- length(intersect(keys, reference_keys))
    data.table(engine = engine, rung = r, reference = reference, n_hits = length(keys),
               n_snps = uniqueN(hits[rung == r, snp]), n_shared_with_reference = shared,
               jaccard_with_reference = signif(shared / length(union(keys, reference_keys)), 3))
  }))
}

lfmm_ladder <- read_ladder("lfmm")
emmax_ladder <- read_ladder("emmax")
result <- rbind(
  compare_rungs("lfmm", ladder_rungs(lfmm_ladder, "K", lfmm_k), paste0("K", config$sNMF$k_best),
                apply_gc = TRUE, lfmm_ladder),
  compare_rungs("emmax", ladder_rungs(emmax_ladder, "npc"), "npc0", apply_gc = FALSE, emmax_ladder)
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
fwrite(result, out_file, sep = "\t")
message("Every rung reproduces the ladder's hits_bonf. Wrote ", out_file)
