#!/usr/bin/env Rscript
# mvp_write_hd_configs.R -- configs for the LD-window plan's Phase 4 (simulation sets, option B).
#
# Two projects per seed, both from config_MVP<seed>_c1.yaml (the onefit cell: EMMAX n_pcs = LFMM
# K = RDA condition_pcs = k_best, 99 permutations):
#
#   MVP<seed>hd    LFMM K = k_best, EMMAX n_pcs = k_best, RDA uncorrected (condition_pcs 0), and the
#                  six declared sets of mvp_build_snp_sets.R / Phase 3 (toppct 0.25 % per predictor,
#                  RDA 0.5 % by override; Union, Cross-method >=2, >=3) at auto_genome_wide under the
#                  critical rule. Regions take the same distance (snp_clumping_distance).
#   MVP<seed>hdrc  RDA conditioned on k_best PCs only (Figure S1), no sets. Cloned from the hd tree
#                  after mode=structure (benchmarks/mvp_clone_project.sh <seed>hd rc).
#
# The suffix follows the seed (MVP<seed>hd, not MVPhd<seed>): every MVP tool builds
# paste0("MVP", seed, suffix), mvp_clone_project.sh included.
#
# Changed keys against _c1, both projects: project_name, cpu: 2, GEA.configs (as above),
# GEA.clumping_ld_rule: critical, GEA.snp_clumping_distance: auto_genome_wide,
# LDdecay.min_samples: 99999 (only the genome-wide "All" group, as the corrected estimator's out600
# run; decided with the user 2026-10-27), LDdecay.scope: genome_wide (as _c1). hd also gets
# GEA.snp_sets. PreGEA is kept in the file but never run.
#
# Usage: Rscript benchmarks/mvp_write_hd_configs.R <seeds_csv|file with one seed per line>
#   run from the repo root (inside the image: /pipeline)

suppressPackageStartupMessages(library(yaml))

a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 1)
seeds <- if (file.exists(a[1])) trimws(readLines(a[1])) else strsplit(a[1], ",")[[1]]
seeds <- seeds[nzchar(seeds)]

rda_pct <- list(RDA = list(type = "toppct", value = 0.5))
pct_set <- function(name, methods = NULL, strategy = "Union") {
    s <- list(name = name, threshold_type = "toppct", threshold_value = 0.25)
    if (!is.null(methods)) s$methods <- methods
    if (is.null(methods) || "RDA" %in% methods) s$overrides <- rda_pct
    s$strategy <- strategy
    s$clumping_distance <- "auto_genome_wide"
    s
}
snp_sets <- list(
    pct_set("solo_lfmm",  methods = list("LFMM")),
    pct_set("solo_emmax", methods = list("EMMAX")),
    pct_set("solo_rda",   methods = list("RDA")),
    pct_set("union"),
    pct_set("best",       strategy = "Cross-method >=2"),
    pct_set("intersect3", strategy = "Cross-method >=3")
)
method_entry <- function(m, params) list(method = m, adjust = "bonf", threshold = "0.05", params = params)

for (s in seeds) {
    src <- sprintf("config_MVP%s_c1.yaml", s)
    if (!file.exists(src)) stop("missing base config ", src)
    base <- yaml::read_yaml(src)
    stopifnot(identical(as.character(base$project_name), paste0("MVP", s)))
    k <- as.integer(base$sNMF$k_best)
    perm <- Filter(function(e) e$method == "RDA", base$GEA$configs)[[1]]$params$permutations
    stopifnot(length(k) == 1, !is.na(k), identical(as.integer(perm), 99L))
    for (sfx in c("hd", "hdrc")) {
        cfg <- base
        cfg$project_name <- paste0("MVP", s, sfx)
        # cpu 2: the old onefit / rdaunc runs passed CPU = 2 to rda.R (Snakemake metadata of
        # MVP<seed>_results); unset, cpu now defaults to host cores - 2 (190 on lab)
        cfg$cpu <- 2L
        cfg$LDdecay$min_samples <- 99999L
        cfg$LDdecay$scope <- "genome_wide"
        cfg$GEA$clumping_ld_rule <- "critical"
        cfg$GEA$snp_clumping_distance <- "auto_genome_wide"
        if (sfx == "hd") {
            cfg$GEA$configs <- list(method_entry("EMMAX", list(n_pcs = k)),
                                    method_entry("LFMM",  list(K = k)),
                                    method_entry("RDA",   list(condition_pcs = 0L, permutations = 99L)))
            cfg$GEA$snp_sets <- snp_sets
        } else {
            cfg$GEA$configs <- list(method_entry("RDA", list(condition_pcs = k, permutations = 99L)))
            cfg$GEA$snp_sets <- NULL
        }
        dst <- sprintf("config_MVP%s%s.yaml", s, sfx)
        hdr <- c(sprintf("# Written by benchmarks/mvp_write_hd_configs.R from %s -- do not hand-edit.", src),
                 sprintf("# LD-window plan Phase 4: %s (k_best %d).", cfg$project_name, k))
        writeLines(c(hdr, yaml::as.yaml(cfg)), dst)
    }
}
message("wrote hd + hdrc configs for ", length(seeds), " seed(s)")
