#!/usr/bin/env Rscript
# pregea_recommend.R — reads every preGEA ladder table and emits ONE tidy
# recommendation TSV. Every row records the RULE that produced it, so the
# GEA method picker's pre-fill (Shiny) is never a black box.
#
# Deliberately NOT wired as a @pregea sentinel in workflow/methods/gea.py's
# _resolve_param_default(): that function runs at Snakemake config-PARSE
# time, before any DAG exists — coupling it to this file would make every
# mode's config parse depend on pregea_recommendations.tsv existing on disk.
# Instead: preGEA writes this artifact -> Shiny's GEA method editor reads it
# and pre-fills fields (with evidence_metric/evidence_value as help text) ->
# the user confirms/edits -> the LITERAL value is written into
# config_<PROJECT>.yaml. The pipeline stays dumb; the config stays a
# complete, reproducible record.

suppressPackageStartupMessages(library(data.table))
source("/pipeline/scripts/R/lib/pregea_calibration.R")  # lambda_in_band()

args <- commandArgs(trailingOnly = TRUE)
################################################################################
LFMM_LADDER_TSV  <- args[1]
EMMAX_LADDER_TSV <- args[2]
RDA_LADDER_TSV   <- args[3]
RDA_COLLIN_TSV   <- args[4]
RDA_AXIS_TSV     <- args[5]                  # per-axis-per-rung detail table (condition_pcs, axis, axis_eig, axis_p)
K_BEST           <- args[6]
LAMBDA_TOL       <- as.numeric(args[7])
DEFLATION_FLOOR  <- as.numeric(args[8])
VIF_MAX          <- as.numeric(args[9])      # real config value now (was hardcoded 10.0 below)
AXIS_ALPHA       <- as.numeric(args[10])     # real config value now (was hardcoded 0.05 below)
OUT_TSV          <- args[11]
################################################################################
# VARPART_FRACTIONS (old arg 6) dropped — varpart moved to the climate module
# in the split; the recommender must not depend on a climate-mode artifact.

dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)

read_or_null <- function(path) if (identical(path, "NULL") || !file.exists(path)) NULL else fread(path, sep = "\t", header = TRUE)

rows <- list()
add_rec <- function(method, param, value, value_type, rule, evidence_metric, evidence_value,
                    alternatives = "", ladder_table = "", ladder_plot = "", note = "") {
    rows[[length(rows) + 1]] <<- data.table(
        method = method, param = param, recommended_value = as.character(value), value_type = value_type,
        rule = rule, evidence_metric = evidence_metric, evidence_value = as.character(evidence_value),
        alternatives = alternatives, ladder_table = ladder_table, ladder_plot = ladder_plot, note = note)
}

################################################################################
# LFMM.K — min hist_flatness_ks among rungs whose hits_qval is >= 50% of the
# ladder max (not yet collapsed). NEVER uses lambda (C.0 — genomic.control
# recalibration forces lambda toward 1 regardless of K correctness).
################################################################################
lfmm <- read_or_null(LFMM_LADDER_TSV)
if (!is.null(lfmm)) {
    pooled <- lfmm[trait == "__pooled__"]
    if (nrow(pooled) > 0) {
        max_hits <- suppressWarnings(max(pooled$hits_qval, na.rm = TRUE))
        eligible <- pooled[!is.na(hits_qval) & (is.infinite(max_hits) | hits_qval >= 0.5 * max_hits)]
        if (nrow(eligible) == 0) eligible <- pooled
        eligible <- eligible[order(hist_flatness_ks, rung_value_num)]
        best <- eligible[1]
        add_rec("LFMM", "K", best$rung_value_num, "int", "min hist_flatness_ks among non-collapsed rungs",
               "hist_flatness_ks", best$hist_flatness_ks,
               alternatives = paste(pooled$rung_value_num, collapse = ","),
               ladder_table = "PreGEA/tables/lfmm/lfmm_ladder.tsv",
               ladder_plot = "PreGEA/plots/lfmm/lfmm_pvalue_histogram_grid.png",
               note = "lambda_GC was deliberately NOT used for K selection — genomic.control recalibration forces it toward 1 regardless of whether K is correct (docs/rda_research.md C.0).")
    }
}

################################################################################
# EMMAX.n_pcs — smallest n_pcs with lambda in [deflation_floor, 1+lambda_tol];
# fallback: min |lambda-1|.
################################################################################
emmax <- read_or_null(EMMAX_LADDER_TSV)
if (!is.null(emmax)) {
    pooled <- emmax[trait == "__pooled__"]
    if (nrow(pooled) > 0) {
        in_range <- pooled[lambda_in_band(lambda_gc, DEFLATION_FLOOR, LAMBDA_TOL)]
        if (nrow(in_range) > 0) {
            in_range <- in_range[order(rung_value_num)]
            best <- in_range[1]
            add_rec("EMMAX", "n_pcs", best$rung_value_num, "int",
                   sprintf("smallest n_pcs with lambda in [%.2f, %.2f]", DEFLATION_FLOOR, 1 + LAMBDA_TOL),
                   "lambda_gc", best$lambda_gc,
                   alternatives = paste(pooled$rung_value_num, collapse = ","),
                   ladder_table = "PreGEA/tables/emmax/emmax_ladder.tsv",
                   ladder_plot = "PreGEA/plots/emmax/emmax_lambda_vs_n_pcs.png")
        } else {
            pooled2 <- pooled[!is.na(lambda_gc)]
            if (nrow(pooled2) > 0) {
                pooled2 <- pooled2[order(abs(lambda_gc - 1))]
                best <- pooled2[1]
                add_rec("EMMAX", "n_pcs", best$rung_value_num, "int", "min |lambda-1| (no rung reached calibration band)",
                       "lambda_gc", best$lambda_gc,
                       alternatives = paste(pooled$rung_value_num, collapse = ","),
                       ladder_table = "PreGEA/tables/emmax/emmax_ladder.tsv",
                       ladder_plot = "PreGEA/plots/emmax/emmax_lambda_vs_n_pcs.png",
                       note = "no rung reached calibration")
            }
        }
    }
}

################################################################################
# RDA.condition_pcs / RDA.axes — smallest rung with max_vif<vif_max,
# n_aliased==0, anova_full_p<axis_alpha, R2adj within 1 SE of ladder max.
################################################################################
rda_ladder <- read_or_null(RDA_LADDER_TSV)
rda_axis   <- read_or_null(RDA_AXIS_TSV)
if (!is.null(rda_ladder)) {
    ok_rows <- rda_ladder[status == "ok"]
    if (nrow(ok_rows) > 0) {
        r2_sd <- stats::sd(ok_rows$r2_adj, na.rm = TRUE)
        r2_max <- max(ok_rows$r2_adj, na.rm = TRUE)
        se <- if (is.na(r2_sd) || nrow(ok_rows) < 2) 0 else r2_sd / sqrt(nrow(ok_rows))
        qualifying <- ok_rows[max_vif < VIF_MAX & n_aliased == 0 &
                              !is.na(anova_full_p) & anova_full_p < AXIS_ALPHA &
                              r2_adj >= r2_max - se]
        # The fallback rung passed none of the rule's tests, so it gets its own rule text: the
        # rule string used to claim it passed significance (SIMDATA's pc2, anova p 0.42; 3602c5).
        rda_rule <- "smallest rung passing VIF/aliasing/significance, R2adj within 1 SE of ladder max"
        rda_note <- ""
        if (nrow(qualifying) == 0) {
            qualifying <- ok_rows[order(-r2_adj)][1]
            rda_rule <- "max R2adj (no rung passed VIF/aliasing/significance within 1 SE of ladder max)"
            rda_note <- sprintf(paste0("no rung passed the rule; smallest full-model anova p %.3g ",
                                       "(axis_alpha %g), smallest max VIF %.3g (vif_max %g)"),
                                suppressWarnings(min(ok_rows$anova_full_p, na.rm = TRUE)), AXIS_ALPHA,
                                suppressWarnings(min(ok_rows$max_vif, na.rm = TRUE)), VIF_MAX)
        }
        qualifying <- qualifying[order(condition_pcs)]
        best <- qualifying[1]
        add_rec("RDA", "condition_pcs", best$condition_pcs, "int", rda_rule,
               "r2_adj", best$r2_adj,
               alternatives = paste(ok_rows$condition_pcs, collapse = ","),
               ladder_table = "PreGEA/tables/rda/rda_condition_ladder.tsv",
               ladder_plot = "PreGEA/plots/rda/rda_model_comparison.png",
               note = rda_note)
        # axis_p_str: per-axis p-values at the recommended rung, from the
        # detailed axis table (RDA_AXIS_TSV) — richer evidence than just the
        # n_axes_sig count already carried by the ladder table itself.
        axis_p_str <- ""
        if (!is.null(rda_axis)) {
            rows_here <- rda_axis[condition_pcs == best$condition_pcs][order(axis)]
            if (nrow(rows_here) > 0) axis_p_str <- paste(sprintf("%.4g", rows_here$axis_p), collapse = ",")
        }
        add_rec("RDA", "axes", max(2, best$n_axes_sig), "int", "n_axes_sig at the recommended condition_pcs rung, floored at 2",
               "n_axes_sig", best$n_axes_sig,
               alternatives = axis_p_str,
               ladder_table = sprintf("PreGEA/tables/rda/models/pc%d/axis_anova.tsv", best$condition_pcs),
               ladder_plot = sprintf("PreGEA/plots/rda/models/pc%d/axis_screeplot.png", best$condition_pcs),
               note = "floored at 2 -- covRob's pairwiseGK estimator needs >= 2 loading columns (verified empirically in rda.R). alternatives column holds the per-axis p-values at this rung.")
    }
}

collin <- read_or_null(RDA_COLLIN_TSV)
if (!is.null(collin)) {
    kept <- collin[action == "kept"]
    if (nrow(kept) > 0) {
        dropped <- collin[action == "dropped"]
        add_rec("RDA", "predictor_set", paste(kept$predictor, collapse = ","), "str",
               "pre-screened + VIF-passing predictor set from the collinearity table",
               "n_kept", nrow(kept),
               note = if (nrow(dropped) > 0)
                   sprintf("Dropped for collinearity: %s. Never auto-applied — this is a suggestion, GEA.configs[RDA].params.predictor_set stays 'auto' unless the user opts in.",
                          paste(dropped$predictor, collapse = ", "))
                   else "No predictors dropped for collinearity.")
    }
}

################################################################################
# Write
################################################################################
out <- if (length(rows) > 0) rbindlist(rows) else data.table(
    method = character(), param = character(), recommended_value = character(),
    value_type = character(), rule = character(), evidence_metric = character(),
    evidence_value = character(), alternatives = character(), ladder_table = character(),
    ladder_plot = character(), note = character())
fwrite(out, OUT_TSV, sep = "\t", quote = FALSE)
message("INFO: preGEA recommendations: ", nrow(out), " row(s) written to ", OUT_TSV)
