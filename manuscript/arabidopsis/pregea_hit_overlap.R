#!/usr/bin/env Rscript
# pregea_hit_overlap.R: do the rungs of the LFMM and EMMAX ladders call the same SNPs?
#
# The ladders hold counts, so two rungs with equal counts may call different SNPs. This counts
# Bonferroni hits per predictor as the production run calls them: LFMM divides each predictor's
# chi-squared by that predictor's own lambda, EMMAX applies no recalibration, and each predictor is
# its own Bonferroni family. A hit is one SNP on one predictor.
#
# Reads _intermediate/pregea/{lfmm,emmax}/<rung>/pvalues.tsv and the two ladder tables. Stops unless
# each rung's count per predictor equals the ladder's hits_bonf. Writes
# ANALYSIS/Arabidopsis1135/pregea/rung_hit_overlap.tsv, one row per rung:
#   engine rung reference n_hits n_snps n_shared_with_reference jaccard_with_reference
#
# The reference is K = k_best for LFMM and n_pcs = 0 for EMMAX. These are pruned markers, so the
# counts rank the rungs and are not a forecast of the whole-genome run.
#
# Run from the repo root, in the image:
#   docker run --rm --name arab-hit-overlap --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp \
#     -e OPENBLAS_NUM_THREADS=1 -v $PWD:/data -w /data cline-go:6cdedcc \
#     Rscript manuscript/arabidopsis/pregea_hit_overlap.R

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(yaml)
})

project <- "Arabidopsis1135"
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

compare_rungs <- function(engine, rungs, reference, apply_gc) {
  ladder <- fread(file.path(results, "PreGEA", "tables", engine, paste0(engine, "_ladder.tsv")),
                  colClasses = c(rung_value = "character"))
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

lfmm_rungs <- paste0("K", 4:8)
emmax_rungs <- paste0("npc", 0:10)
result <- rbind(
  compare_rungs("lfmm", lfmm_rungs, paste0("K", config$sNMF$k_best), apply_gc = TRUE),
  compare_rungs("emmax", emmax_rungs, "npc0", apply_gc = FALSE)
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
fwrite(result, out_file, sep = "\t")
message("Every rung reproduces the ladder's hits_bonf. Wrote ", out_file)
