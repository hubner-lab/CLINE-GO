#!/usr/bin/env Rscript
# mvp_write_p3_configs.R -- configs for the LD-window plan's Phase 3 check: the framework's
# declared SNP sets (GEA.snp_sets) must equal the sets benchmarks/mvp_build_snp_sets.R built
# for offset13, on the same old p-value trio.
#
# Each config is config_MVP<seed>_c1.yaml with only these keys changed:
#   project_name                 MVP<seed><suffix>
#   GEA.clumping_ld_rule         critical (explicit, never the default)
#   GEA.snp_clumping_distance    auto_genome_wide
#   LDdecay.min_samples          99999 -- only the genome-wide "All" group, as the corrected
#                                estimator's out600 run (work/ld_phase3_2026-10-26); the
#                                per-deme groups are never read
#   GEA.snp_sets                 the six sets of mvp_build_snp_sets.R at 5 kb, plus best_ld and
#                                intersect3_ld at auto_genome_wide (the seed's critical distance)
#
# The set rules mirror mvp_build_snp_sets.R:129-139/:226-268: top 0.25 % per predictor for LFMM
# and EMMAX, 0.5 % for RDA's single multivariate column (0.25 % x 2 predictors), agreement by
# methods within the clumping distance. RDA's rule goes in through `overrides`, because a
# set-level threshold pins a non-univariate method to its registry rule (common.smk, snp-set
# parser).
#
# Usage: Rscript benchmarks/mvp_write_p3_configs.R <seeds_csv> [suffix=p3]
#   run from the repo root (inside the image: /pipeline)

suppressPackageStartupMessages(library(yaml))

args   <- commandArgs(trailingOnly = TRUE)
seeds  <- strsplit(args[1], ",")[[1]]
suffix <- if (length(args) >= 2) args[2] else "p3"
stopifnot(length(seeds) > 0, nzchar(suffix))

rda_pct  <- list(RDA = list(type = "toppct", value = 0.5))
pct_set  <- function(name, methods = NULL, strategy = "Union", cd = 5000L) {
    s <- list(name = name, threshold_type = "toppct", threshold_value = 0.25)
    if (!is.null(methods)) s$methods <- methods
    if (is.null(methods) || "RDA" %in% methods) s$overrides <- rda_pct
    s$strategy <- strategy
    s$clumping_distance <- cd
    s
}
snp_sets <- list(
    pct_set("solo_lfmm",  methods = list("LFMM")),
    pct_set("solo_emmax", methods = list("EMMAX")),
    pct_set("solo_rda",   methods = list("RDA")),
    pct_set("union"),
    pct_set("best",       strategy = "Cross-method >=2"),
    pct_set("intersect3", strategy = "Cross-method >=3"),
    pct_set("best_ld",       strategy = "Cross-method >=2", cd = "auto_genome_wide"),
    pct_set("intersect3_ld", strategy = "Cross-method >=3", cd = "auto_genome_wide")
)

for (s in seeds) {
    src <- sprintf("config_MVP%s_c1.yaml", s)
    dst <- sprintf("config_MVP%s%s.yaml", s, suffix)
    if (!file.exists(src)) stop("missing base config ", src)
    cfg <- yaml::read_yaml(src)
    stopifnot(identical(as.character(cfg$project_name), s) ||
              identical(as.character(cfg$project_name), paste0("MVP", s)))
    cfg$project_name <- paste0("MVP", s, suffix)
    cfg$GEA$clumping_ld_rule <- "critical"
    cfg$GEA$snp_clumping_distance <- "auto_genome_wide"
    cfg$LDdecay$min_samples <- 99999L
    cfg$GEA$snp_sets <- snp_sets
    hdr <- c(sprintf("# Written by benchmarks/mvp_write_p3_configs.R from %s -- do not hand-edit.", src),
             "# LD-window plan Phase 3: declared sets vs mvp_build_snp_sets.R on the old p-value trio.")
    writeLines(c(hdr, yaml::as.yaml(cfg)), dst)
    message("wrote ", dst)
}
