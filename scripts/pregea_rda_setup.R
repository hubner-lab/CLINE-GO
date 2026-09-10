#!/usr/bin/env Rscript
# pregea_rda_setup.R — RDA exploratory setup on the LD-pruned set (preGEA
# Block 3). Diagnostics only (CLAUDE.md Rule 5) — does NOT call
# find_sig_snps.R, write a p-value table, or touch any GEA output.
#
# Deliberately NOT a wildcard fan-out: every Condition()-PC rung refits the
# SAME Y (the LD-pruned imputed matrix) with a different Condition() PC
# count — a fan-out would re-read the genotype matrix N times for no
# caching benefit, and the deliverable is a single comparison table.
#
# Four analyses:
#   A. Predictor collinearity pre-screen (|r| < collinearity_r) + post-fit
#      VIF (A16, a REAL gate now — a rung whose max_vif >= vif_max is
#      flagged status="vif_exceeded" and excluded from "ok" rows used
#      downstream, not just cosmetically rescaled into a plot line).
#      NEVER drops below MIN_PREDICTORS — keeps the least-correlated pair
#      and records the deviation.
#   B. Condition()-PC ladder: for each condition_pcs in [cmin, cmax], fit a
#      partial RDA, run anova.cca (full/by-axis/by-term), call candidates via
#      the SAME scripts/R/lib/rdadapt.R the GEA-mode RDA scan uses. EVERY
#      rung also gets its own biplot + axis screeplot + axis anova table
#      (models/pc{n}/...) so the Shiny RDA tab can let the user inspect any
#      swept value, not just a single "best rung" (module-split follow-up —
#      the old single best-rung-only rda_biplot_best.png/rda_axis_screeplot.png
#      are retired in favor of this).
#   C. Plain-vs-partial candidate-count comparison (one uncorrected fit vs
#      every partial rung).
#   D. ordiR2step forward-selection path — cumulative adjusted R2 per step,
#      not just the final variable set (A17).

library(data.table)
library(dplyr)
library(vegan)
library(robust)
library(qvalue)
library(stringr)
library(ggplot2)
library(ggrepel)
library(scattermore)

source("/pipeline/scripts/R/utils/theme_clinego.R")
source("/pipeline/scripts/R/utils/emmax_core.R")   # load_pca_covariates()
source("/pipeline/scripts/R/lib/rdadapt.R")         # rdadapt() — shared with rda.R

args <- commandArgs(trailingOnly = TRUE)
################################################################################
LFMM_PRUNED     <- args[1]
VCFSNP_PRUNED   <- args[2]
CLIMATE         <- args[3]
PCA             <- args[4]
SAMPLES_ORDER   <- args[5]
CLIMATE_VALID   <- args[6]
PREDICTORS      <- args[7]
COND_PCS_MIN    <- as.integer(args[8])
COND_PCS_MAX    <- as.integer(args[9])
COLLIN_R        <- as.numeric(args[10])
VIF_MAX         <- as.numeric(args[11])
MIN_PREDICTORS  <- as.integer(args[12])
AXIS_ALPHA      <- as.numeric(args[13])
PERMUTATIONS    <- as.numeric(args[14])
ORDIR2STEP_PIN  <- as.numeric(args[15])
R2_PERMUTATIONS <- as.numeric(args[16])
SEED            <- as.numeric(args[17])
CPU             <- as.numeric(args[18])
OUT_COLLIN_TSV  <- args[19]
OUT_LADDER_TSV  <- args[20]
OUT_AXIS_TSV    <- args[21]
OUT_FWD_TSV     <- args[22]
PLOT_DIR        <- args[23]
TABLE_DIR       <- args[24]
INTER_DIR       <- args[25]
################################################################################

# Per-model artifact paths — one set per condition_pcs value, mirroring
# common.smk's rda_model_dir()/rda_model_biplot()/rda_model_screeplot()/
# rda_model_axis_anova() Python helpers exactly (same PLOT_DIR/TABLE_DIR
# convention: models/pc{n}/...). Any drift here breaks Snakemake's declared
# rule outputs, not just a cosmetic mismatch.
rda_model_plot_dir  <- function(n) file.path(PLOT_DIR, "models", paste0("pc", n))
rda_model_table_dir <- function(n) file.path(TABLE_DIR, "models", paste0("pc", n))

for (d in c(dirname(OUT_COLLIN_TSV), dirname(OUT_LADDER_TSV), PLOT_DIR, TABLE_DIR, INTER_DIR))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("INFO: preGEA RDA setup — condition_pcs range [", COND_PCS_MIN, ",", COND_PCS_MAX, "]")

predictor_names <- str_split(PREDICTORS, ',')[[1]] %>% str_trim()

################################################################################
# 1. Load data — same shape as rda.R
################################################################################
climate_dt <- fread(CLIMATE, sep = '\t', header = TRUE)
missing_preds <- setdiff(predictor_names, colnames(climate_dt))
if (length(missing_preds) > 0) stop("preGEA RDA setup: predictor(s) not found: ", paste(missing_preds, collapse = ", "))
Env <- as.data.frame(climate_dt[, ..predictor_names])  # data.frame, not data.table:
                                                        # Env[, active, drop=FALSE] below needs
                                                        # base-R `[` semantics for a variable
                                                        # column-name vector (data.table's `[`
                                                        # requires the `..active` NSE prefix instead)
n_samples <- nrow(Env)

Y_full <- fread(LFMM_PRUNED, sep = ' ', header = FALSE) %>% as.matrix()
if (nrow(Y_full) != n_samples) stop("preGEA RDA setup: genotype/climate row mismatch (", nrow(Y_full), " vs ", n_samples, ")")
col_sds <- apply(Y_full, 2, sd, na.rm = TRUE)
kept_idx <- which(col_sds > 0 & !is.na(col_sds))
Y <- Y_full[, kept_idx, drop = FALSE]
m <- ncol(Y)
message("INFO: m=", m, " markers (LD-pruned, ", length(col_sds) - m, " zero-variance dropped), n=", n_samples, " samples")

climate_valid_ids <- fread(CLIMATE_VALID, header = FALSE, colClasses = "character")$V2

################################################################################
# 2. Predictor collinearity pre-screen (A16) — never below MIN_PREDICTORS
#
# SITE-WEIGHTED, not sample-weighted. Env has one row per SAMPLE and a site's
# predictor values are identical for every sample sequenced there, so cor() over
# Env rows weights each site by its sample size — a site with 30 individuals
# counts 30x against a site with 1, even though both contribute exactly one
# environment. That is pseudoreplication, and unlike the z-scoring in
# download_climate_present.R (a per-column affine transform, which cannot change
# any linear fit) this correlation drives a DISCRETE decision: which predictors
# enter every RDA rung below. Simulation, 8 predictors x 11 sites, 200 draws of
# the loop below: with balanced n the two weightings agree 200/200 and the
# correlation matrices are identical; at n 1-9 per site the median max
# |r_sample - r_site| is 0.25 and they keep a DIFFERENT predictor set in 78% of
# draws; at SIMDATA-like imbalance (n 1-30) 93%.
# unique() IS the site-level table here — predictor values are constant within a
# site, so distinct rows are distinct environments and no metadata join is needed.
################################################################################
Env_site <- unique(as.data.frame(Env))
if (nrow(Env_site) < 3) {
    # Fewer than 3 distinct environments: a site-level correlation is degenerate
    # (2 points make every |r| exactly 1). Fall back to the per-sample matrix so
    # the screen still runs, and say so — the RDA that follows is barely
    # interpretable at this design size either way.
    message("WARNING: only ", nrow(Env_site), " distinct site environment(s) across ",
            n_samples, " samples — site-level correlation is degenerate, ",
            "falling back to the per-sample correlation matrix for the collinearity screen.")
    cor_mat <- cor(as.data.frame(Env), use = "pairwise.complete.obs")
} else {
    cor_mat <- cor(Env_site, use = "pairwise.complete.obs")
}
active <- predictor_names
collin_rows <- list()
repeat {
    if (length(active) <= MIN_PREDICTORS) break
    sub <- cor_mat[active, active, drop = FALSE]
    diag(sub) <- 0
    worst <- which(abs(sub) == max(abs(sub)), arr.ind = TRUE)[1, ]
    r_worst <- sub[worst[1], worst[2]]
    if (abs(r_worst) < COLLIN_R) break
    p1 <- rownames(sub)[worst[1]]; p2 <- colnames(sub)[worst[2]]
    # Drop whichever of the pair is, on average, more correlated with the
    # REST of the active set (keeps the more "independent" predictor).
    mean_r <- function(p) mean(abs(sub[p, setdiff(active, p)]))
    drop_p <- if (mean_r(p1) >= mean_r(p2)) p1 else p2
    collin_rows[[length(collin_rows) + 1]] <- data.table(
        predictor = drop_p, action = "dropped",
        reason = sprintf("|r|=%.3f with %s (>= collinearity_r=%.2f)", r_worst,
                         setdiff(c(p1, p2), drop_p), COLLIN_R))
    active <- setdiff(active, drop_p)
}
if (length(active) < MIN_PREDICTORS) {
    # Should not happen given the MIN_PREDICTORS floor above, but never let a
    # collinearity screen silently produce an unfittable RDA (_validate_rda_semantics
    # already enforces >=2 at the GEA-config level; preGEA must not contradict it).
    stop("preGEA RDA setup: collinearity screen left ", length(active),
        " predictor(s), below MIN_PREDICTORS=", MIN_PREDICTORS,
        " — widen collinearity_r or reduce min_predictors.")
}
for (p in active) collin_rows[[length(collin_rows) + 1]] <- data.table(predictor = p, action = "kept", reason = "")
collin_dt <- rbindlist(collin_rows)
message("INFO: collinearity screen (",
       if (nrow(Env_site) < 3) "per-sample fallback, " else "site-weighted, ",
       nrow(Env_site), " distinct environments / ", n_samples, " samples) kept ",
       length(active), "/", length(predictor_names), " predictors: ",
       paste(active, collapse = ","))

################################################################################
# 3. Structure-correction PCs (shared helper — k<=0 = no Condition())
################################################################################
get_covariates <- function(k) load_pca_covariates(PCA, k, SAMPLES_ORDER, climate_valid_ids)

fit_partial_rda <- function(cond_pcs) {
    covariates <- if (cond_pcs > 0) get_covariates(cond_pcs) else NULL
    Variables <- as.data.frame(Env[, active, drop = FALSE])
    rhs_terms <- active
    if (!is.null(covariates)) {
        Variables <- cbind(Variables, as.data.frame(covariates))
        rhs_terms <- c(rhs_terms, paste0("Condition(", paste(colnames(covariates), collapse = "+"), ")"))
    }
    f <- as.formula(paste("Y ~", paste(rhs_terms, collapse = " + ")))
    list(fit = rda(f, data = Variables, scale = TRUE), variables = Variables)
}

################################################################################
# 4. Condition()-PC ladder (Block 3B) — refit once per rung
################################################################################
axis_rows   <- list()
ladder_rows <- list()

# Diagnostic candidate rule for THIS exploratory table only (NOT the applied
# GEA rule — that lives in GEA.configs and is applied by rda.R itself).
# Fixed Bonferroni-on-m, matching rda.R's own sensitivity-ladder convention.
diagnostic_bonf <- function(pvals, mm) sum(pvals < 0.05 / mm, na.rm = TRUE)

# Per-rung model artifacts (biplot + axis screeplot + axis anova table) —
# written for EVERY swept condition_pcs value, not just a single "best rung"
# (module-split follow-up: the Shiny RDA tab now lets the user pick which
# model to inspect). `fit`/`axis_pvals`/`axis_eigs` are NULL/empty for rungs
# that failed to fit or reach rank>=2 — the empty-state placeholders explain
# why rather than leaving a 0-byte file.
write_model_artifacts <- function(cond_pcs, fit = NULL, axis_pvals = NULL, axis_eigs = NULL, reason = NULL) {
    p_dir <- rda_model_plot_dir(cond_pcs); t_dir <- rda_model_table_dir(cond_pcs)
    dir.create(p_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(t_dir, recursive = TRUE, showWarnings = FALSE)

    # Axis anova table
    if (!is.null(axis_pvals) && length(axis_pvals) > 0) {
        axis_here <- data.table(condition_pcs = cond_pcs, axis = seq_along(axis_pvals),
                                axis_eig = axis_eigs, axis_p = axis_pvals)
    } else {
        axis_here <- data.table(condition_pcs = cond_pcs, axis = NA_integer_,
                                axis_eig = NA_real_, axis_p = NA_real_)
    }
    fwrite(axis_here, file.path(t_dir, "axis_anova.tsv"), sep = "\t", quote = FALSE)

    # Axis screeplot
    if (!is.null(axis_pvals) && length(axis_pvals) > 0) {
        g_scree <- ggplot(axis_here, aes(x = factor(axis), y = axis_eig, fill = axis_p < AXIS_ALPHA)) +
            geom_col() +
            scale_fill_manual(values = c(`TRUE` = CLINEGO_RETAINED, `FALSE` = CLINEGO_THRESHOLD), guide = "none") +
            labs(x = "Constrained axis", y = "Eigenvalue",
                title = "RDA constrained-axis screeplot",
                subtitle = sprintf("condition_pcs=%d", cond_pcs)) +
            theme_clinego()
    } else {
        g_scree <- clinego_empty_plot(reason %||% "No constrained axes available")
    }
    clinego_save_both(file.path(p_dir, "axis_screeplot"), g_scree, w = 7, h = 5)

    # Site-scores biplot — plain scatter is correct here (site/sample scores,
    # N=samples not SNPs; Rule 6 does not apply, see docs/rda_research.md C.2).
    if (!is.null(fit) && fit$CCA$rank >= 2) {
        site_scores <- as.data.frame(scores(fit, display = "sites", choices = 1:2))
        colnames(site_scores) <- c("RDA1", "RDA2")
        biplot_scores <- as.data.frame(fit$CCA$biplot[, 1:2, drop = FALSE])
        colnames(biplot_scores) <- c("RDA1", "RDA2")
        biplot_scores$predictor <- rownames(fit$CCA$biplot)
        arrow_scale <- 0.8 * max(abs(site_scores$RDA1), abs(site_scores$RDA2)) /
            max(abs(biplot_scores$RDA1), abs(biplot_scores$RDA2), 1e-6)
        g_biplot <- ggplot(site_scores, aes(x = RDA1, y = RDA2)) +
            geom_point(color = CLINEGO_NEUTRAL, size = 2) +
            geom_segment(data = biplot_scores,
                        aes(x = 0, y = 0, xend = RDA1 * arrow_scale, yend = RDA2 * arrow_scale),
                        inherit.aes = FALSE, color = CLINEGO_THRESHOLD, arrow = arrow(length = unit(0.2, "cm"))) +
            ggrepel::geom_text_repel(data = biplot_scores,
                                     aes(x = RDA1 * arrow_scale, y = RDA2 * arrow_scale, label = predictor),
                                     inherit.aes = FALSE, color = CLINEGO_THRESHOLD, size = 3) +
            labs(x = "RDA1", y = "RDA2", title = "RDA site scores biplot",
                subtitle = sprintf("condition_pcs=%d", cond_pcs)) +
            theme_clinego()
    } else {
        g_biplot <- clinego_empty_plot(reason %||% "rank<2: no biplot available")
    }
    clinego_save_both(file.path(p_dir, "biplot"), g_biplot, w = 7, h = 6)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

for (cond_pcs in COND_PCS_MIN:COND_PCS_MAX) {
    message("INFO: rung condition_pcs=", cond_pcs)
    fitted <- tryCatch(fit_partial_rda(cond_pcs), error = function(e) {
        message("WARNING: rung condition_pcs=", cond_pcs, " fit failed: ", e$message); NULL
    })
    if (is.null(fitted)) {
        ladder_rows[[length(ladder_rows) + 1]] <- data.table(
            condition_pcs = cond_pcs, status = "fit_failed", note = "rda() fit raised an error")
        write_model_artifacts(cond_pcs, reason = "rda() fit failed at this condition_pcs value")
        next
    }
    fit <- fitted$fit
    aliased <- fit$CCA$alias
    n_aliased <- if (is.null(aliased)) 0L else length(aliased)
    K_max <- fit$CCA$rank
    r2adj <- suppressWarnings(RsquareAdj(fit)$adj.r.squared)
    vif_vals <- tryCatch(vif.cca(fit), error = function(e) NA_real_)
    max_vif  <- suppressWarnings(max(vif_vals, na.rm = TRUE))
    mean_vif <- suppressWarnings(mean(vif_vals, na.rm = TRUE))

    if (K_max < 2) {
        ladder_rows[[length(ladder_rows) + 1]] <- data.table(
            condition_pcs = cond_pcs, n_predictors = length(active),
            predictors = paste(active, collapse = ","), r2 = NA_real_, r2_adj = r2adj,
            max_vif = max_vif, mean_vif = mean_vif, n_aliased = n_aliased, rank = K_max,
            n_axes_sig = NA_integer_, anova_full_F = NA_real_, anova_full_p = NA_real_,
            term_min_p = NA_real_, term_max_p = NA_real_,
            n_candidates_plain = NA_integer_, n_candidates_partial = NA_integer_,
            gif_lambda = NA_real_, status = "rank_lt_2",
            note = "rank<2: covRob requires >=2 loading columns (verified empirically in rda.R)")
        write_model_artifacts(cond_pcs, fit = fit, reason = "rank<2: no constrained axes to show")
        next
    }

    set.seed(SEED)
    anova_full <- tryCatch(anova.cca(fit, permutations = PERMUTATIONS, parallel = max(1, CPU)),
                           error = function(e) NULL)
    set.seed(SEED)
    anova_axis <- tryCatch(anova.cca(fit, by = "axis", permutations = PERMUTATIONS), error = function(e) NULL)
    set.seed(SEED)
    anova_term <- tryCatch(anova.cca(fit, by = "terms", permutations = PERMUTATIONS), error = function(e) NULL)

    anova_full_F <- if (!is.null(anova_full)) as.data.frame(anova_full)[["F"]][1] else NA_real_
    anova_full_p <- if (!is.null(anova_full)) as.data.frame(anova_full)[["Pr(>F)"]][1] else NA_real_

    axis_pvals <- if (!is.null(anova_axis)) as.data.frame(anova_axis)[["Pr(>F)"]][seq_len(K_max)] else rep(NA_real_, K_max)
    axis_eigs  <- if (!is.null(anova_axis)) as.data.frame(anova_axis)[["Variance"]][seq_len(K_max)] else rep(NA_real_, K_max)
    for (a in seq_len(K_max)) {
        axis_rows[[length(axis_rows) + 1]] <- data.table(
            condition_pcs = cond_pcs, axis = a, axis_eig = axis_eigs[a], axis_p = axis_pvals[a])
    }
    n_axes_sig <- sum(axis_pvals < AXIS_ALPHA, na.rm = TRUE)

    term_p <- if (!is.null(anova_term)) {
        tp <- as.data.frame(anova_term)[["Pr(>F)"]]
        tp[is.finite(tp)]
    } else numeric(0)
    term_min_p <- if (length(term_p) > 0) min(term_p) else NA_real_
    term_max_p <- if (length(term_p) > 0) max(term_p) else NA_real_

    # K for candidate calling — same rule rda.R uses: floor at 2, cap at K_max.
    K_sel <- max(2, min(n_axes_sig, K_max))
    rd <- tryCatch(rdadapt(fit, K_sel), error = function(e) {
        message("WARNING: rdadapt failed for condition_pcs=", cond_pcs, ": ", e$message); NULL
    })
    n_cand_partial <- if (!is.null(rd)) diagnostic_bonf(rd$p.values, m) else NA_integer_
    gif_lambda     <- if (!is.null(rd)) rd$gif_lambda else NA_real_

    # VIF as a REAL gate (Phase 4 fix): a rung whose max_vif reaches vif_max
    # is flagged and excluded from "ok" rows used by the recommender and the
    # best-rung selection below — previously vif_max only rescaled a plot
    # line and never filtered anything.
    vif_ok <- is.na(max_vif) || max_vif < VIF_MAX
    status <- if (vif_ok) "ok" else "vif_exceeded"
    note   <- if (vif_ok) "" else sprintf("max_vif=%.2f exceeds vif_max=%.2f", max_vif, VIF_MAX)

    ladder_rows[[length(ladder_rows) + 1]] <- data.table(
        condition_pcs = cond_pcs, n_predictors = length(active),
        predictors = paste(active, collapse = ","), r2 = suppressWarnings(RsquareAdj(fit)$r.squared),
        r2_adj = r2adj, max_vif = max_vif, mean_vif = mean_vif, n_aliased = n_aliased, rank = K_max,
        n_axes_sig = n_axes_sig, anova_full_F = anova_full_F, anova_full_p = anova_full_p,
        term_min_p = term_min_p, term_max_p = term_max_p,
        n_candidates_plain = NA_integer_,  # filled in after the plain (cond_pcs=0) reference fit below
        n_candidates_partial = n_cand_partial, gif_lambda = gif_lambda, status = status, note = note)
    write_model_artifacts(cond_pcs, fit = fit, axis_pvals = axis_pvals, axis_eigs = axis_eigs)
}

ladder_dt <- rbindlist(ladder_rows, fill = TRUE)
axis_dt   <- rbindlist(axis_rows,  fill = TRUE)

# Plain (uncorrected, condition_pcs=0) reference candidate count — computed
# ONCE, not per rung, and back-filled into every ladder row for the
# plain-vs-partial comparison column (Block 3C).
plain_fitted <- tryCatch(fit_partial_rda(0), error = function(e) NULL)
n_cand_plain <- NA_integer_
if (!is.null(plain_fitted) && plain_fitted$fit$CCA$rank >= 2) {
    K_plain <- min(2, plain_fitted$fit$CCA$rank)
    rd_plain <- tryCatch(rdadapt(plain_fitted$fit, K_plain), error = function(e) NULL)
    if (!is.null(rd_plain)) n_cand_plain <- diagnostic_bonf(rd_plain$p.values, m)
}
if (nrow(ladder_dt) > 0) ladder_dt[, n_candidates_plain := n_cand_plain]

fwrite(collin_dt, OUT_COLLIN_TSV, sep = "\t", quote = FALSE)
fwrite(ladder_dt, OUT_LADDER_TSV, sep = "\t", quote = FALSE)
fwrite(axis_dt,   OUT_AXIS_TSV,   sep = "\t", quote = FALSE)
message("INFO: Wrote collinearity (", nrow(collin_dt), " rows), ladder (", nrow(ladder_dt),
       " rows), axis (", nrow(axis_dt), " rows) tables")

################################################################################
# 5. ordiR2step forward-selection path (A17) — cumulative R2adj per step,
#    not just the final variable set. Unconditioned (no Condition() term),
#    matching Capblancq & Forester 2021's reference code (A.3).
################################################################################
fwd_dt <- data.table(step = integer(), variable = character(),
                     r2_adj_cumulative = numeric(), r2_adj_full_ceiling = numeric(),
                     df = numeric(), F = numeric(), p = numeric())
Variables_full <- as.data.frame(Env[, active, drop = FALSE])
mod0 <- tryCatch(rda(Y ~ 1, data = Variables_full), error = function(e) NULL)
modfull <- tryCatch(rda(as.formula(paste("Y ~", paste(active, collapse = " + "))),
                        data = Variables_full, scale = TRUE), error = function(e) NULL)
full_ceiling <- if (!is.null(modfull)) suppressWarnings(RsquareAdj(modfull)$adj.r.squared) else NA_real_

if (!is.null(mod0) && !is.null(modfull)) {
    step_res <- tryCatch(
        vegan::ordiR2step(mod0, scope = formula(modfull), Pin = ORDIR2STEP_PIN,
                          R2scope = TRUE, R2permutations = R2_PERMUTATIONS, trace = FALSE),
        error = function(e) { message("WARNING: ordiR2step failed: ", e$message); NULL }
    )
    if (!is.null(step_res)) {
        selected_terms <- attr(terms(step_res), "term.labels")
        # Re-fit incrementally so cumulative R2adj is computed directly from
        # the model, not parsed out of ordiR2step's internal $anova table
        # (whose exact column shape varies across vegan versions).
        cum_terms <- character(0)
        for (i in seq_along(selected_terms)) {
            cum_terms <- c(cum_terms, selected_terms[i])
            step_fit <- tryCatch(
                rda(as.formula(paste("Y ~", paste(cum_terms, collapse = " + "))),
                   data = Variables_full, scale = TRUE),
                error = function(e) NULL)
            if (is.null(step_fit)) next
            r2adj_cum <- suppressWarnings(RsquareAdj(step_fit)$adj.r.squared)
            av <- tryCatch(anova.cca(step_fit, permutations = R2_PERMUTATIONS), error = function(e) NULL)
            fwd_dt <- rbind(fwd_dt, data.table(
                step = i, variable = selected_terms[i], r2_adj_cumulative = r2adj_cum,
                r2_adj_full_ceiling = full_ceiling,
                df = if (!is.null(av)) as.data.frame(av)$Df[1] else NA_real_,
                F  = if (!is.null(av)) as.data.frame(av)[["F"]][1] else NA_real_,
                p  = if (!is.null(av)) as.data.frame(av)[["Pr(>F)"]][1] else NA_real_))
        }
    }
}
fwrite(fwd_dt, OUT_FWD_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote ordiR2step forward-selection path (", nrow(fwd_dt), " steps)")

################################################################################
# 6. Plots
################################################################################
# NOTE: the old top-level "best rung only" collinearity heatmap, axis
# screeplot, and biplot are RETIRED — every rung already got its own
# biplot/screeplot/axis-anova table inside the main loop above (module split:
# the Shiny RDA tab now selects between them). The predictor collinearity
# heatmap itself is redundant with climate/plots/correlation_heatmap.png —
# which predictors got dropped for collinearity is still fully recorded in
# OUT_COLLIN_TSV (predictor, action, reason).
ok_rows <- ladder_dt[status == "ok"]
best_rung <- if (nrow(ok_rows) > 0) ok_rows[which.max(r2_adj)]$condition_pcs else NA_integer_

# -- Model comparison panel: r2_adj / max_vif / n_axes_sig / anova_full_p
#    across the whole condition_pcs ladder, each on its own facet with its
#    own y-scale + threshold reference line. Replaces the old
#    rda_condition_ladder.png dual-axis-by-rescaling hack (a VIF line with no
#    readable scale).
if (nrow(ladder_dt) > 0) {
    cmp_long <- rbindlist(list(
        ladder_dt[, .(condition_pcs, metric = "r2_adj", value = r2_adj, threshold = NA_real_)],
        ladder_dt[, .(condition_pcs, metric = "max_vif", value = max_vif, threshold = VIF_MAX)],
        ladder_dt[, .(condition_pcs, metric = "n_axes_sig", value = as.numeric(n_axes_sig), threshold = NA_real_)],
        ladder_dt[, .(condition_pcs, metric = "anova_full_p", value = anova_full_p, threshold = AXIS_ALPHA)]
    ))
    cmp_long[, metric := factor(metric, levels = c("r2_adj", "max_vif", "n_axes_sig", "anova_full_p"))]
    g_cmp <- ggplot(cmp_long, aes(x = condition_pcs, y = value)) +
        geom_line(color = CLINEGO_NEUTRAL) + geom_point(color = CLINEGO_NEUTRAL, size = 2) +
        geom_hline(aes(yintercept = threshold), color = CLINEGO_THRESHOLD, linetype = "dashed", na.rm = TRUE) +
        facet_wrap(~ metric, scales = "free_y",
                  labeller = as_labeller(c(r2_adj = "Adjusted R2", max_vif = "Max VIF",
                                          n_axes_sig = "# significant axes", anova_full_p = "Full-model p"))) +
        labs(x = "Condition() PC count", y = NULL, title = "RDA model comparison across the Condition()-PC ladder",
            subtitle = sprintf("dashed lines: vif_max=%.1f, axis_alpha=%.2g", VIF_MAX, AXIS_ALPHA)) +
        theme_clinego_grid()
    clinego_save_both(file.path(PLOT_DIR, "rda_model_comparison"), g_cmp, w = 9, h = 6)
} else {
    clinego_save_both(file.path(PLOT_DIR, "rda_model_comparison"),
                    clinego_empty_plot("No condition_pcs rung fit successfully"), w = 9, h = 6)
}

# -- ordiR2step path: cumulative R2adj per step + full-model ceiling ------
if (nrow(fwd_dt) > 0) {
    g_fwd <- ggplot(fwd_dt, aes(x = step, y = r2_adj_cumulative)) +
        geom_hline(aes(yintercept = r2_adj_full_ceiling), color = CLINEGO_THRESHOLD, linetype = "dashed") +
        geom_line(color = CLINEGO_NEUTRAL) + geom_point(color = CLINEGO_NEUTRAL, size = 2) +
        ggrepel::geom_text_repel(aes(label = variable), size = 3, color = CLINEGO_COL$fg) +
        labs(x = "Forward-selection step", y = "Cumulative adjusted R2",
            title = "ordiR2step forward-selection path",
            subtitle = sprintf("Pin=%.3g | dashed = full-model ceiling (%.4f)", ORDIR2STEP_PIN, full_ceiling)) +
        theme_clinego()
    clinego_save_both(file.path(PLOT_DIR, "rda_ordir2step_path"), g_fwd, w = 7, h = 5)
} else {
    clinego_save_both(file.path(PLOT_DIR, "rda_ordir2step_path"),
                    clinego_empty_plot(sprintf("No variable passed Pin=%.3g (full-model R2adj=%.4f)",
                                             ORDIR2STEP_PIN, if (is.na(full_ceiling)) 0 else full_ceiling)),
                    w = 7, h = 5)
}

message("INFO: preGEA RDA setup complete. Best rung by R2adj: condition_pcs=", best_rung)
