#!/usr/bin/env Rscript
# =============================================================================
# mvp_figures.R -- publication-ready figures and tables for the CLINE-GO
# genetic-offset benchmark.
#
# TWO ARGUMENTS, KEPT SEPARATE ON PURPOSE:
#
#   A. PIPELINE PERFORMANCE (figures 1-3). What the scenario-aware maladaptation
#      arm costs versus the per-garden driver it replaced. Every number is read
#      from benchmarks/mvp_perf_measured.tsv, which records only measurements
#      taken from run logs -- no estimates, and each row carries its source.
#
#   B. SCIENTIFIC RESULT (figures 4-6). What the benchmark says about genomic
#      offset: which marker panels work, how that depends on genetic
#      architecture, and what predicts offset accuracy. Read from the scoring
#      tables in benchmarks/mvp_eval/offset09/.
#
# Conflating the two would be the easy mistake: "our pipeline is fast" and "GEA
# panels recover the oracle" are different claims with different evidence, and a
# reviewer should be able to check them independently.
#
# Every figure ships as PNG + SVG, and every figure's underlying numbers ship as
# a TSV of the same stem, so nothing in a panel is unverifiable from the repo.
#
# Usage:
#   Rscript benchmarks/mvp_figures.R [--outdir=DIR] [--evaldir=DIR]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/theme_clinego.R"))

args <- commandArgs(trailingOnly = TRUE)
kv <- function(k, d) {
    hit <- grep(paste0("^--", k, "="), args, value = TRUE)
    if (!length(hit)) return(d)
    sub(paste0("^--", k, "="), "", hit[1])
}
EVAL <- kv("evaldir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset09"))
OUT  <- kv("outdir",  file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/figures"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Save PNG + SVG + the numbers behind the panel. A figure whose table is missing
# is a figure nobody can check.
save_fig <- function(p, stem, data, w = 8, h = 5) {
    ggsave(file.path(OUT, paste0(stem, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
    ggsave(file.path(OUT, paste0(stem, ".svg")), p, width = w, height = h, bg = "white")
    fwrite(data, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
    message(sprintf("  %-34s %d rows", stem, nrow(data)))
}

message("== A. pipeline performance ==")
P <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_perf_measured.tsv"))
getv <- function(m, s, v) P[metric == m & stage == s & variant == v]$value[1]

# ---- Fig 1: compute per seed, by rule, before/after -------------------------
# Stacked bars rather than one total, so it is visible WHERE the cost went and
# which method still dominates. geometric_offset is deliberately left tall: its
# LFMM2 fit cannot be hoisted out of the scenario loop (LEA::genetic.gap() fits
# and projects in one call), and hiding that would overstate the result.
f1 <- rbindlist(list(
    data.table(variant = "old driver\n(per-garden)", rule = "all rules (median 917 s x 42 gardens)",
               seconds = getv("compute_per_seed", "offset_sweep", "old_driver")),
    data.table(variant = "scenario pipeline", rule = "gradient_forest_offset",
               seconds = getv("rule_total", "gradient_forest_offset", "new_postbatch")),
    data.table(variant = "scenario pipeline", rule = "geometric_offset",
               seconds = getv("rule_total", "geometric_offset", "new_prebatch")),
    data.table(variant = "scenario pipeline", rule = "rda_offset (x2 methods)",
               seconds = getv("rule_total", "rda_offset", "new_prebatch")),
    data.table(variant = "scenario pipeline", rule = "stage_custom_climate_future",
               seconds = getv("rule_total", "stage_custom_climate_future", "new_prebatch"))
))
tot_old <- f1[variant %like% "old", sum(seconds)]
tot_new <- f1[variant %like% "scenario", sum(seconds)]

p1 <- ggplot(f1, aes(variant, seconds / 3600, fill = rule)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(CLINEGO_CATEGORICAL[1:4], CLINEGO_NEUTRAL)) +
    labs(title = "Compute per replicate for one garden sweep",
         subtitle = sprintf("%.1f h -> %.1f h  (%.1fx less compute), 42 gardens x 7 panels x 4 methods",
                            tot_old / 3600, tot_new / 3600, tot_old / tot_new),
         x = NULL, y = "CPU-hours per replicate", fill = NULL) +
    theme_clinego() +
    guides(fill = guide_legend(nrow = 2))
save_fig(p1, "fig1_compute_per_seed", f1)

# ---- Fig 2: what batching predict() did -------------------------------------
# The whole 42-scenario sweep was gated on ONE job: Gradient Forest projecting an
# ~850 MB forest 42 times. Each scenario is only 100 raster cells, so the cost was
# per-call overhead, not data volume -- which is why stacking the scenarios into a
# single predict() call collapses it without changing a single output value.
f2 <- rbindlist(list(
    data.table(stat = "total (7 panels)", phase = "per-scenario calls",
               seconds = getv("rule_total", "gradient_forest_offset", "new_prebatch")),
    data.table(stat = "total (7 panels)", phase = "one batched call",
               seconds = getv("rule_total", "gradient_forest_offset", "new_postbatch")),
    data.table(stat = "slowest job (neutral_all)", phase = "per-scenario calls",
               seconds = getv("rule_max", "gradient_forest_offset", "new_prebatch")),
    data.table(stat = "slowest job (neutral_all)", phase = "one batched call",
               seconds = getv("rule_max", "gradient_forest_offset", "new_postbatch")),
    data.table(stat = "median job", phase = "per-scenario calls",
               seconds = getv("rule_median", "gradient_forest_offset", "new_prebatch")),
    data.table(stat = "median job", phase = "one batched call",
               seconds = getv("rule_median", "gradient_forest_offset", "new_postbatch"))
))
f2[, phase := factor(phase, levels = c("per-scenario calls", "one batched call"))]
f2[, stat := factor(stat, levels = c("total (7 panels)", "slowest job (neutral_all)", "median job"))]

p2 <- ggplot(f2, aes(stat, seconds, fill = phase)) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    geom_text(aes(label = paste0(seconds, " s")), position = position_dodge(0.7),
              vjust = -0.35, size = 3, colour = CLINEGO_COL$fg) +
    scale_fill_manual(values = c(CLINEGO_REMOVED, CLINEGO_RETAINED)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Gradient Forest projection: 42 separate predict() calls vs one",
         subtitle = "Outputs byte-identical before and after; each scenario is only 100 raster cells, so this was call overhead",
         x = NULL, y = "seconds", fill = NULL) +
    theme_clinego()
save_fig(p2, "fig2_gf_batched_predict", f2)

# ---- Fig 3: orchestration shape ---------------------------------------------
# Log scale: the counts differ by orders of magnitude, and a linear axis would
# render the new bars invisible rather than small.
f3 <- rbindlist(list(
    data.table(quantity = "Snakemake invocations", variant = "old driver",
               count = getv("invocations_per_seed", "offset_sweep", "old_driver")),
    data.table(quantity = "Snakemake invocations", variant = "scenario pipeline",
               count = getv("invocations_per_seed", "offset_sweep", "new")),
    data.table(quantity = "offset jobs", variant = "old driver",
               count = getv("offset_jobs_per_seed", "offset_sweep", "old_driver")),
    data.table(quantity = "offset jobs", variant = "scenario pipeline",
               count = getv("offset_jobs_per_seed", "offset_sweep", "new")),
    data.table(quantity = "clone lanes", variant = "old driver",
               count = getv("clone_lanes_per_seed", "offset_sweep", "old_driver")),
    data.table(quantity = "clone lanes", variant = "scenario pipeline",
               count = getv("clone_lanes_per_seed", "offset_sweep", "new"))
))
f3[, variant := factor(variant, levels = c("old driver", "scenario pipeline"))]
p3 <- ggplot(f3, aes(quantity, pmax(count, 0.5), fill = variant)) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    geom_text(aes(label = count), position = position_dodge(0.7), vjust = -0.35,
              size = 3, colour = CLINEGO_COL$fg) +
    scale_fill_manual(values = c(CLINEGO_REMOVED, CLINEGO_RETAINED)) +
    scale_y_log10(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Orchestration per replicate",
         subtitle = "log scale; clone lanes plotted at 0.5 where the true value is 0",
         x = NULL, y = "count per replicate (log10)", fill = NULL) +
    theme_clinego()
save_fig(p3, "fig3_orchestration", f3)

message("== B. scientific result ==")

read_if <- function(f) {
    p <- file.path(EVAL, f)
    if (!file.exists(p)) { message("  SKIP (missing): ", f); return(NULL) }
    fread(p)
}

# ---- Fig 4: median tau by method x marker panel ------------------------------
# Kendall's tau between predicted offset and realised common-garden fitness.
# MORE NEGATIVE IS BETTER: a good offset predicts LOW fitness. The sign is kept
# raw rather than flipped, so the numbers match the tables and the source paper.
mk <- read_if("marker_set_comparison_by_method.tsv")
gp <- read_if("garden_performance.tsv")
if (!is.null(gp) && nrow(gp)) {
    setDT(gp)
    landscape <- if ("type" %in% names(gp)) gp[type == "landscape"] else gp
    if ("is_control" %in% names(landscape)) landscape <- landscape[is_control == FALSE]
    f4 <- landscape[, .(median_tau = median(tau, na.rm = TRUE), n = .N),
                    by = .(method_label, marker_set)]
    p4 <- ggplot(f4, aes(reorder(marker_set, median_tau), median_tau, fill = method_label)) +
        geom_col(position = position_dodge(0.8), width = 0.72) +
        coord_flip() +
        scale_fill_manual(values = CLINEGO_CATEGORICAL[1:4]) +
        labs(title = "Offset accuracy by marker panel and method",
             subtitle = "Kendall's tau vs common-garden fitness; more negative = better prediction",
             x = NULL, y = expression(median~Kendall*"'"*s~tau), fill = NULL) +
        theme_clinego()
    save_fig(p4, "fig4_tau_by_panel_method", f4, w = 9, h = 5.5)

    # ---- Fig 5: the architecture claim --------------------------------------
    # The headline: whether curating markers matters at all depends on genetic
    # architecture. Oligogenic replicates gain from curation; highly-polygenic
    # ones barely do. Averaging over architecture is what hides this.
    man <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
                 colClasses = c(seed = "character"))
    if ("seed" %in% names(landscape)) {
        # Seeds are numeric-looking IDs, so fread infers integer on whichever side
        # lacks an explicit colClasses -- and an integer/character join throws.
        # Force both to character before merging (CLAUDE.md, R script conventions).
        landscape[, seed := as.character(seed)]
        man[, seed := as.character(seed)]
        L <- merge(landscape, man[, .(seed, arch_level, arm)], by = "seed", all.x = TRUE)
        L <- L[arm == "primary" & !is.na(arch_level)]
        f5 <- L[, .(median_tau = median(tau, na.rm = TRUE), n_models = .N),
                by = .(arch_level, marker_set, method_label)]
        p5 <- ggplot(f5, aes(marker_set, median_tau, fill = method_label)) +
            geom_col(position = position_dodge(0.8), width = 0.72) +
            facet_wrap(~ arch_level, ncol = 1) +
            coord_flip() +
            scale_fill_manual(values = CLINEGO_CATEGORICAL[1:4]) +
            labs(title = "Marker curation pays off only for some architectures",
                 subtitle = "SS-Mtn replicates only; more negative = better",
                 x = NULL, y = expression(median~Kendall*"'"*s~tau), fill = NULL) +
            theme_clinego()
        save_fig(p5, "fig5_tau_by_architecture", f5, w = 9, h = 9)
    }
}

# ---- Fig 6: what predicts offset accuracy -----------------------------------
# Local adaptation predicts offset accuracy; F_ST does not. One point per
# replicate (n <= 32), so a scatter is the honest geometry here -- this is not
# WGS-scale data where points would collapse into noise.
tl <- read_if("tau_vs_LA.tsv")
if (!is.null(tl) && nrow(tl)) {
    setDT(tl)
    ycol <- intersect(c("tau", "median_tau"), names(tl))[1]
    xcands <- intersect(c("final_LA", "LA", "meanFst"), names(tl))
    if (!is.na(ycol) && length(xcands)) {
        f6 <- melt(tl, id.vars = intersect(c("seed", "method_label", ycol), names(tl)),
                   measure.vars = xcands, variable.name = "property", value.name = "x")
        p6 <- ggplot(f6, aes(x, get(ycol))) +
            geom_point(colour = CLINEGO_NEUTRAL, alpha = 0.75, size = 1.8) +
            geom_smooth(method = "lm", se = TRUE, colour = CLINEGO_THRESHOLD, linewidth = 0.6) +
            facet_wrap(~ property, scales = "free_x") +
            labs(title = "Local adaptation predicts offset accuracy; F_ST does not",
                 subtitle = "one point per replicate x method; more negative tau = better",
                 x = NULL, y = expression(Kendall*"'"*s~tau)) +
            theme_clinego()
        save_fig(p6, "fig6_tau_vs_properties", f6, w = 9, h = 4.5)
    }
}

message(sprintf("\nfigures + tables written to %s", OUT))
