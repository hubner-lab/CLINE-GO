#!/usr/bin/env Rscript
# =============================================================================
# mvp_block_report.R -- the step-12 report for one SS-Clines demography block.
#
# WHY THIS EXISTS. The SS-Clines arm is a COMPLETE BLOCK: 3 genic levels x 4
# architecture sub-levels x 10 replicates = 120, with no threshold applied to any
# realized statistic. The legacy arm was band-selected on r2_pc1_temp, which meant
# selecting on an OUTCOME. A complete block has no selection rule inside it to
# attack -- so confounding stops being a criterion and becomes a reported
# covariate. This script is what reports it.
#
# Three products, per benchmarks/SSCLINES_RUNBOOK.md step 12:
#   1. per-cell medians          -- genic level x architecture sub-level
#   2. realized covariate table  -- what the paper reports INSTEAD of selecting on
#   3. accuracy vs final_LA      -- the continuous re-plot that replaced the
#                                   two-cloud equal-S/unequal-S figure
# Plus a block-completeness table, because an underfilled marker panel is a
# RESULT (n = 0 SNPs, recall = 0), not missing data, and the row arithmetic has
# to be shown rather than assumed.
#
# THE SUB-LEVEL IS NOT `demog_level_sub`. That column holds the DEMOGRAPHY
# (N-variable_m-variable), which is constant within a block -- grouping by it
# gives one cell. The four architecture sub-levels are the suffix of the
# `architecture` column: {pleiotropy, no-pleiotropy} x {equal-S, unequal-S}.
#
# Reads only tables already on disk. Fits nothing, runs no pipeline mode.
#
# Usage:
#   OFFSET_DIR=offset12_ssclines_b1 MVP_ARM=primary_ssclines \
#   MVP_ADDED=ssclines_nvar_mvar MVP_N_EXPECT=120 \
#   FIG_OUT=/pipeline/benchmarks/mvp_eval/figures_ssclines_b1 \
#     Rscript /pipeline/benchmarks/mvp_block_report.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset12_ssclines_b1")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_ssclines_b1"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))

rd <- function(...) {
    f <- file.path(...)
    if (!file.exists(f)) stop("MISSING: ", f)
    fread(f, colClasses = c("seed" = "character"))
}
emit <- function(dt, stem) {
    fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
    message(sprintf("  wrote %s.tsv (%d rows)", stem, nrow(dt)))
    invisible(dt)
}
save_fig <- function(p, stem, w, h) {
    clinego_save_both(file.path(OUT, stem), p, w = w, h = h); message("  OK ", stem)
}

ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")   # corpus typo
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
# Same exclusion rule as mvp_oracle_stats.R, and for the same reason: RDA-corrected
# predicts backwards on a few percent of gardens. It is kept in the per-cell table
# (where the reader can see it) and dropped from the final_LA regression.
WORKING <- c("GFoffset", "LFMM2offset", "RDA-uncorrected")
ORACLE_SET <- "adaptive"          # the true causal loci
BEST_SET   <- "gea_best"          # the >=2-of-3 recommended panel

# ------------------------------------------------------------------- manifest
seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
PRIM  <- mvp_prim(seeds)
PRIM[, arch_lab := factor(arch_level, levels = ARCH_LEVELS, labels = ARCH_LABELS)]
stopifnot(nrow(PRIM) == mvp_n_expect(), !any(is.na(PRIM$arch_lab)))

# Architecture sub-level = `architecture` with the "<arch_level>_" prefix removed.
PRIM[, sub_level := sub("^[^_]+_", "", architecture)]
message("architecture sub-levels: ", paste(sort(unique(PRIM$sub_level)), collapse = " | "))
# The complete-block design asserted, not assumed. If a block is not 3 x 4 x 10 it
# is not complete, and every "no selection rule inside the block" claim built on
# this report is false.
CELLS <- PRIM[, .N, by = .(arch_lab, sub_level)]
setorder(CELLS, arch_lab, sub_level)
if (uniqueN(CELLS$N) != 1L || nrow(CELLS) != 12L)
    message("!! block is NOT the complete 3 x 4 x 10 design -- cells: ",
            paste(sprintf("%s/%s=%d", CELLS$arch_lab, CELLS$sub_level, CELLS$N),
                  collapse = ", "))
stopifnot(nrow(CELLS) == 12L)
message(sprintf("design: %d cells x %d replicates = %d",
                nrow(CELLS), CELLS$N[1], sum(CELLS$N)))

# ------------------------------------------------------------------ scored data
G <- rd(EVAL, OFF, "garden_performance.tsv")
G <- G[seed %in% PRIM$seed]
stopifnot(uniqueN(G$seed) == nrow(PRIM))
message(sprintf("garden_performance: %d rows, %d seeds, %d marker sets, %d methods",
                nrow(G), uniqueN(G$seed), uniqueN(G$marker_set), uniqueN(G$method_label)))

# Cross-check against the scoring log, which medians over the landscape gardens
# directly. If this table disagrees with logs_score_<cohort>.log the arm/cohort
# filter is selecting the wrong seeds and every number below is wrong.
LAND <- G[garden_type == "landscape" & control == FALSE]
CHK  <- dcast(LAND[, .(tau = median(tau)), by = .(marker_set, method_label)],
              marker_set ~ method_label, value.var = "tau")
emit(CHK, "block_marker_set_medians")
message("\n=== median tau per marker set x offset method (landscape, non-control) ===")
message("    cross-check this against logs_score_<cohort>.log")
print(CHK)

# =============================================================================
# 1. per-cell medians -- genic level x architecture sub-level
# =============================================================================
# Aggregation stated once: median over the landscape gardens within
# (seed, set, method), then median over replicates within the cell. A replicate
# is one unit however many gardens it carries.
SEED_MED <- LAND[, .(tau = median(tau)),
                 by = .(seed, marker_set, method_label)]
SEED_MED[, accuracy := -tau]        # NEVER abs(): an anti-predicting model must
                                    # stay negative, not become accurate
SEED_MED <- merge(SEED_MED, PRIM[, .(seed, arch_lab, sub_level, final_LA)], by = "seed")

CELL_MED <- SEED_MED[, .(n_replicates = uniqueN(seed),
                         accuracy = median(accuracy),
                         q25 = quantile(accuracy, .25),
                         q75 = quantile(accuracy, .75)),
                     by = .(arch_lab, sub_level, marker_set, method_label)]
setorder(CELL_MED, arch_lab, sub_level, marker_set, method_label)
emit(CELL_MED, "block_cell_medians")

# The headline cut of the same table: the recommended panel against the oracle,
# per cell, in paired delta units.
PAIR <- dcast(SEED_MED[method_label %in% WORKING &
                       marker_set %in% c(ORACLE_SET, BEST_SET)],
              seed + arch_lab + sub_level + method_label ~ marker_set,
              value.var = "accuracy")
PAIR[, delta := get(BEST_SET) - get(ORACLE_SET)]
PAIR <- PAIR[is.finite(delta)]
CELL_DELTA <- PAIR[, .(n = .N,
                       median_delta = median(delta),
                       q25 = quantile(delta, .25), q75 = quantile(delta, .75),
                       reach_pct = 100 * mean(delta >= 0)),
                   by = .(arch_lab, sub_level)]
setorder(CELL_DELTA, arch_lab, sub_level)
emit(CELL_DELTA, "block_cell_delta_best_vs_oracle")
message("\n=== >=2-of-3 panel minus causal loci, per cell (>=0 means it reached) ===")
print(CELL_DELTA)

CELL_HEAT <- ggplot(CELL_DELTA, aes(sub_level, arch_lab, fill = median_delta)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%+.3f\n%.0f%% reach", median_delta, reach_pct)),
              size = 3, colour = "grey15") +
    scale_fill_gradient2(low = CLINEGO_REMOVED, mid = "grey92", high = CLINEGO_RETAINED,
                         midpoint = 0, name = "median delta") +
    theme_clinego() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(x = "architecture sub-level", y = NULL,
         title = "Recommended panel minus causal loci, per cell",
         subtitle = sprintf("%s -- %d replicates per cell, 3 working offset methods",
                            mvp_arm_label(), CELLS$N[1]))
save_fig(CELL_HEAT, "block_cell_delta_heatmap", w = 9, h = 4.5)

# =============================================================================
# 2. realized covariates -- reported, never selected on
# =============================================================================
COVARS <- c("r2_pc1_temp", "r2_pc1_sal", "meanFst", "K_authors",
            "final_LA", "n_snps", "n_causal_maf01")
stopifnot(all(COVARS %in% names(PRIM)))
COV_LONG <- melt(PRIM[, c("seed", "arch_lab", "sub_level", COVARS), with = FALSE],
                 id.vars = c("seed", "arch_lab", "sub_level"),
                 variable.name = "covariate", value.name = "value")
COV_TAB <- COV_LONG[, .(n = .N, min = min(value), q25 = quantile(value, .25),
                        median = median(value), q75 = quantile(value, .75),
                        max = max(value)),
                    by = .(covariate, arch_lab)]
COV_POOL <- COV_LONG[, .(n = .N, min = min(value), q25 = quantile(value, .25),
                         median = median(value), q75 = quantile(value, .75),
                         max = max(value)),
                     by = .(covariate)][, arch_lab := "pooled"]
emit(rbind(COV_TAB, COV_POOL, use.names = TRUE), "block_realized_covariates")
message("\n=== realized covariates, pooled over the block ===")
print(COV_POOL[, .(covariate, min, median, max)])

# The confounding one gets its own panel, because it is the statistic the legacy
# arm SELECTED on and this arm does not. Showing the realized spread is the whole
# argument that the block has no selection rule inside it.
CONF <- ggplot(PRIM, aes(r2_pc1_temp)) +
    geom_histogram(bins = 30, fill = CLINEGO_NEUTRAL, colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = c(0.20, 0.75), linetype = "dashed",
               colour = CLINEGO_THRESHOLD, linewidth = 0.5) +
    facet_wrap(~ arch_lab) +
    theme_clinego() +
    labs(x = expression(tau^2 * "(PC1, temperature)"), y = "replicates",
         title = "Realized structure-environment confounding across the block",
         subtitle = sprintf("%s -- dashed lines mark the legacy arm's selection band, applied here to NOTHING",
                            mvp_arm_label()))
save_fig(CONF, "block_confounding_realized", w = 9, h = 3.6)

# =============================================================================
# 3. accuracy vs final_LA -- the continuous re-plot
# =============================================================================
# Replaces the equal-S / unequal-S two-cloud figure. Symmetry is a two-level knob
# (SIGMA_K_2 in {0.5, 4.0}) acting THROUGH final_LA, so plotting the continuous
# covariate turns two clouds into one relationship and keeps the selection regime
# out of the analysis as a factor.
LA <- SEED_MED[method_label %in% WORKING & marker_set %in% c(ORACLE_SET, BEST_SET)]
LA[, panel := factor(marker_set, levels = c(ORACLE_SET, BEST_SET),
                     labels = c("causal loci", ">=2-of-3 panel"))]
LA_SEED <- LA[, .(accuracy = median(accuracy)),
              by = .(seed, arch_lab, final_LA, panel, method_label)]
emit(LA_SEED, "block_accuracy_vs_LA")

RHO <- LA_SEED[, .(rho = suppressWarnings(cor(final_LA, accuracy, method = "spearman")),
                   n = .N),
               by = .(panel, method_label, arch_lab)]
emit(RHO, "block_accuracy_vs_LA_rho")
message("\n=== Spearman rho(final_LA, accuracy), pooled over architecture ===")
print(dcast(LA_SEED[, .(rho = suppressWarnings(cor(final_LA, accuracy, method = "spearman"))),
                    by = .(panel, method_label)],
            panel ~ method_label, value.var = "rho"))

# geom_smooth, not geom_point-only: 120 replicates x 3 methods is dense enough
# that the cloud hides the relationship (CLAUDE.md code-structure rule 6).
P_LA <- ggplot(LA_SEED, aes(final_LA, accuracy, colour = panel, fill = panel)) +
    geom_point(alpha = 0.25, size = 1) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.8, alpha = 0.15) +
    facet_grid(method_label ~ arch_lab) +
    scale_color_clinego() + scale_fill_clinego() +
    theme_clinego() + theme(legend.position = "bottom") +
    labs(x = "final local adaptation (final_LA, deposit-provided)",
         y = expression("accuracy = " * -tau),
         colour = NULL, fill = NULL,
         title = "Offset accuracy against realized local adaptation",
         subtitle = sprintf("%s -- %d replicates; LA is a covariate, not a selection criterion",
                            mvp_arm_label(), nrow(PRIM)))
save_fig(P_LA, "block_accuracy_vs_LA", w = 10, h = 7.5)

# =============================================================================
# 4. block completeness -- an empty panel is a result, not missing data
# =============================================================================
N_GARDENS <- uniqueN(G$garden_id)
N_METHODS <- uniqueN(G$method_label)
N_SETS    <- uniqueN(G$marker_set)
FULL      <- nrow(PRIM) * N_GARDENS * N_SETS * N_METHODS
PRESENT   <- G[, .N, by = .(seed, marker_set)][, .N]
EXPECTED  <- nrow(PRIM) * N_SETS
DROPPED   <- EXPECTED - PRESENT

COMPLETE <- data.table(
    quantity = c("replicates", "gardens", "marker sets", "offset methods",
                 "rows if every (seed, set) ran", "rows actually scored",
                 "(seed, set) pairs expected", "(seed, set) pairs present",
                 "(seed, set) pairs dropped", "rows per dropped pair",
                 "rows accounted for by drops"),
    value = c(nrow(PRIM), N_GARDENS, N_SETS, N_METHODS, FULL, nrow(G),
              EXPECTED, PRESENT, DROPPED, N_GARDENS * N_METHODS,
              DROPPED * N_GARDENS * N_METHODS))
emit(COMPLETE, "block_completeness")
message("\n=== block completeness ===")
print(COMPLETE)
if (FULL - nrow(G) != DROPPED * N_GARDENS * N_METHODS)
    message("!! row arithmetic does NOT close -- ", FULL - nrow(G), " missing but ",
            DROPPED * N_GARDENS * N_METHODS, " explained by dropped panels")

# Which (seed, panel) pairs are missing, so the list can be diffed against the
# step-6 pre-screen (panel_underfilled_<cohort>.tsv) rather than trusted.
GRID <- CJ(seed = PRIM$seed, marker_set = unique(G$marker_set), unique = TRUE)
MISS <- GRID[!G[, .(seed, marker_set)], on = .(seed, marker_set)]
setorder(MISS, marker_set, seed)
emit(MISS, "block_dropped_panels")
message(sprintf("dropped (seed, panel) pairs: %d -- %s", nrow(MISS),
                paste(sprintf("%s=%d", MISS[, .N, by = marker_set]$marker_set,
                              MISS[, .N, by = marker_set]$N), collapse = ", ")))

message("\nBlock report written to ", OUT)
