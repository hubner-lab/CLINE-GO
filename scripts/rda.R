#!/usr/bin/env Rscript
# RDA genome scan — multivariate GEA method (Capblancq & Forester 2021 "Swiss
# Army Knife for landscape genomics" / Capblancq et al. 2018 rdadapt).
#
# Fits a PARTIAL RDA across ALL configured climate predictors jointly
# (Condition(PC1..PCk) from the LD-pruned LEA PCA for structure correction),
# then calls candidate SNPs via a robust-Mahalanobis distance of SNP loadings
# over the retained constrained axes (rdadapt), GIF-rescaled, converted to a
# chi-squared p-value. B6 (docs/rda_research.md, 2026-07-29 session): a SECOND,
# UNCONSTRAINED fit (same predictors, no Condition()) is run whenever
# condition_pcs > 0, K pinned to the partial fit's K, and the two rdadapt()
# p-value vectors are combined via pmax() before being written as
# "climate_multivariate" — a SNP only stays significant if BOTH fits call it,
# the continuous-p-value analogue of Capblancq & Forester 2021's
# Reduce(intersect, ...) tutorial pattern. Unlike EMMAX/LFMM this is genuinely
# multivariate — RDA emits ONE p-value per SNP (the "climate_multivariate"
# pseudo-trait), not one column per predictor. See docs/rda_research.md Part A
# for the full decision table and literature citations behind every design
# choice below (referenced inline as A1, A6, A11, etc.), and the B6 row in
# Part B for the candidate-intersection rationale (scoped there to the offset;
# applied here to the GEA scan itself, this session's deliberate broadening).
#
# Full spec: docs/rda_research.md Part A. Literature default (A6, verified
# from Capblancq & Forester 2021's own code comments) is to fit the ordination
# on the FULL marker set — they explicitly chose NOT to pre-prune, handling LD
# entirely post-hoc via clumping (done downstream by create_regions.R, not by
# this script). Pruned-set fitting here is an explicit, loudly-labeled
# ENGINEERING FALLBACK for when the full-matrix fit is infeasible, never a
# silent substitution.

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
source("/pipeline/scripts/R/utils/emmax_core.R")  # load_pca_covariates()
source("/pipeline/scripts/R/utils/pval_threshold.R")  # compute_pval_threshold()
source("/pipeline/scripts/R/lib/rdadapt.R")  # rdadapt() — shared with preGEA's RDA-setup block

args <- commandArgs(trailingOnly = TRUE)
################################################################################
LFMM_FULL      <- args[1]                    # full marker set, imputed, coord-valid (space-sep, no header)
LFMM_PRUNED    <- args[2]                    # LD-pruned equivalent (fallback fit only)
CLIMATE        <- args[3]                    # climate_present_site_scaled.tsv (tab, header, col1=sample)
VCFSNP_FULL    <- args[4]                    # *.vcfsnp for the FULL marker set (space-sep, no header, V1=chr V2=pos)
VCFSNP_PRUNED  <- args[5]                    # *.vcfsnp for the LD-pruned marker set (fallback fit only —
                                              # a pruned fit has fewer SNPs than the full set, so it needs
                                              # its OWN vcfsnp; using the full one would silently mismatch
                                              # row counts between the genotype matrix and the output table)
PCA            <- args[6]                    # LEA .projections (space-sep, no header), LD-pruned
SAMPLES_ORDER  <- args[7]                    # full-cohort sample order (one per line)
CLIMATE_VALID  <- args[8]                    # plink --keep list (FID IID) of climate-valid samples
PREDICTORS     <- args[9]                    # comma-separated predictor names, e.g. "bio_1,bio_2,bio_3"
N_COND_PCS     <- as.numeric(args[10])       # 0 disables Condition()
AXES           <- args[11]                   # "auto" or an integer >= 2
AXIS_ALPHA     <- as.numeric(args[12])
PERMUTATIONS   <- as.numeric(args[13])
FIT_MODE       <- args[14]                   # "auto" | "full" | "pruned"
MAX_SNPS_FULL  <- as.numeric(args[15])
MAX_GB         <- as.numeric(args[16])
SEED           <- as.numeric(args[17])
CPU            <- as.numeric(args[18])
K_BEST         <- args[19]                   # sNMF K — filenames/titles only, NOT the RDA axis count
OUT_PVALUES    <- args[20]
OUT_CANDIDATES <- args[21]
OUT_DIAGNOSTICS<- args[22]
OUT_ANOVA      <- args[23]
PLOT_DIR       <- args[24]
# Appended after the outputs (positions 25-26) deliberately: rda.R takes
# positional args with no schema check, so inserting mid-list would silently
# shift OUT_PVALUES..PLOT_DIR. The out-of-group position is the cost of that
# safety. Single caller: workflow/rules/gea.smk. These carry the candidate
# significance rule (GEA.configs[RDA].adjust/threshold) — linking the
# candidates table to the same rule find_sig_snps.R applies downstream.
CAND_ADJUST    <- args[25]                   # "bonf" | "qval" | "top" | "custom"
CAND_VALUE     <- suppressWarnings(as.numeric(args[26]))
################################################################################

if (!CAND_ADJUST %in% c("bonf", "qval", "top", "custom")) {
    stop("RDA: unknown adjust '", CAND_ADJUST, "' — expected bonf|qval|top|custom. ",
         "This must match GEA.configs[RDA].adjust; see workflow/methods/gea.py.")
}
if (is.na(CAND_VALUE) || CAND_VALUE <= 0) {
    stop("RDA: candidate threshold value '", args[26], "' is not a positive number.")
}

# Candidate-table cap — not a CLI arg (an internal safety valve, not a
# scientific parameter). At WGS density with a heavy-tailed null the raw
# candidate list can run into the millions of rows; the wide p-value table
# (OUT_PVALUES) is always complete regardless of this cap.
MAX_CANDIDATES_WRITTEN <- 100000

dir.create(dirname(OUT_PVALUES), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_CANDIDATES), recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# diag: accumulates key/value rows written to OUT_DIAGNOSTICS at the end —
# every guard/decision below records here so the diagnostics file is the
# single, always-written record of what actually happened in this run (the
# anti-silent-substitution mechanism; see docs/rda_research.md A4/A6).
diag <- list()
add_diag <- function(key, value) diag[[key]] <<- value

message("INFO: Starting RDA genome scan (multivariate)")
message(paste0("INFO: Predictors: ", PREDICTORS))
message(paste0("INFO: condition_pcs=", N_COND_PCS, "  axes=", AXES,
               "  fit_mode=", FIT_MODE))

predictor_names <- str_split(PREDICTORS, ',')[[1]] %>% str_trim()
add_diag("predictors", paste(predictor_names, collapse = ","))
add_diag("n_predictors", length(predictor_names))
add_diag("condition_pcs", N_COND_PCS)

# Provenance rows — every CLI parameter that isn't otherwise recorded below,
# so the diagnostics table alone answers "what was this run actually run
# with" without cross-referencing the Snakemake shell command or config YAML.
add_diag("axes_requested",     AXES)
add_diag("axis_alpha",         AXIS_ALPHA)
add_diag("permutations",       PERMUTATIONS)
add_diag("seed",               SEED)
add_diag("fit_mode_requested", FIT_MODE)
add_diag("full_fit_max_snps",  MAX_SNPS_FULL)
add_diag("full_fit_max_gb",    MAX_GB)
add_diag("k_best",             K_BEST)
add_diag("cpu",                CPU)

################################################################################
# 1. Climate predictors
################################################################################
message("INFO: Reading climate predictors")
climate_dt <- fread(CLIMATE, sep = '\t', header = TRUE)
missing_preds <- setdiff(predictor_names, colnames(climate_dt))
if (length(missing_preds) > 0) {
    stop("RDA: requested predictor(s) not found in climate table: ",
         paste(missing_preds, collapse = ", "))
}
Env <- climate_dt[, ..predictor_names]
n_samples <- nrow(Env)
add_diag("n_samples", n_samples)

################################################################################
# 2. Structure-correction PCs (Condition()) — A10
################################################################################
# Reuses the exact same helper EMMAX uses, including its k<=0 early-return fix
# (scripts/R/utils/emmax_core.R) — so condition_pcs=0 is a first-class "no
# structure correction" mode, not a crash.
climate_valid_ids <- fread(CLIMATE_VALID, header = FALSE, colClasses = "character")$V2
covariates <- load_pca_covariates(PCA, N_COND_PCS, SAMPLES_ORDER, climate_valid_ids)
if (!is.null(covariates)) {
    add_diag("structure_proxy", "LD-pruned LEA PCA (not intergenic-restricted) — A10 deviation")
    if (nrow(covariates) != n_samples) {
        stop("RDA: PCA covariate row count (", nrow(covariates),
             ") does not match climate row count (", n_samples, ")")
    }
} else {
    add_diag("structure_proxy", "none (condition_pcs=0 — uncorrected scan)")
}

################################################################################
# 3. Size gate — choose full vs LD-pruned marker set (A6)
################################################################################
message("INFO: Determining marker set (fit_mode=", FIT_MODE, ")")
count_matrix_dims <- function(path) {
    # First line only — column count via whitespace split, row count via wc -l.
    first_line <- readLines(path, n = 1)
    ncols <- length(strsplit(trimws(first_line), "\\s+")[[1]])
    nrows <- as.numeric(system(paste("wc -l <", shQuote(path)), intern = TRUE))
    c(n = nrows, m = ncols)
}
dims_full <- count_matrix_dims(LFMM_FULL)
m_full <- dims_full[["m"]]
est_bytes_full <- 8 * dims_full[["n"]] * m_full * 4
add_diag("n_markers_available", m_full)
add_diag("est_bytes_full", est_bytes_full)

FALLBACK_WARNING <- paste0(
    "================================================================\n",
    "WARNING: RDA ORDINATION FIT FELL BACK TO THE LD-PRUNED MARKER SET\n",
    "================================================================\n",
    "This is an ENGINEERING FALLBACK, not a scientific alternative.\n",
    "Capblancq & Forester (2021, MEE 12:2298-2309) state in their own code\n",
    "comments they \"preferred not to [prune]\" and handled LD entirely\n",
    "post-hoc via clumping. Pre-pruning is a documented DEVIATION.\n",
    "Trigger: m=", m_full, " > full_fit_max_snps=", MAX_SNPS_FULL,
    " (or) est. size ", round(est_bytes_full / 1e9, 2), "GB > ", MAX_GB, "GB\n",
    "Candidate p-values below are NOT comparable to a full-set fit.\n",
    "================================================================"
)

use_pruned <- FALSE
if (FIT_MODE == "pruned") {
    use_pruned <- TRUE
} else if (FIT_MODE == "full") {
    use_pruned <- FALSE
} else { # auto
    use_pruned <- (m_full > MAX_SNPS_FULL) || (est_bytes_full > MAX_GB * 1e9)
}

if (use_pruned) {
    message(FALLBACK_WARNING)
    LFMM_PATH <- LFMM_PRUNED
    VCFSNP    <- VCFSNP_PRUNED
    add_diag("fallback_warning", FALLBACK_WARNING)
} else {
    LFMM_PATH <- LFMM_FULL
    VCFSNP    <- VCFSNP_FULL
    add_diag("fallback_warning", "")
}
add_diag("fit_mode", if (use_pruned) "pruned" else "full")

################################################################################
# 4. Read genotype matrix
################################################################################
message("INFO: Reading genotype matrix from ", LFMM_PATH)
Y_full <- fread(LFMM_PATH, sep = ' ', header = FALSE) %>% as.matrix()
if (nrow(Y_full) != n_samples) {
    stop("RDA: genotype matrix has ", nrow(Y_full), " rows but climate table has ",
         n_samples, " rows — row alignment is broken.")
}
if (max(Y_full, na.rm = TRUE) > 2) {
    stop("RDA: genotype matrix contains values > 2 — this looks like an ",
         "UNIMPUTED LFMM file (LEA writes 9 for missing). RDA requires an ",
         "imputed matrix (*_imp*). Check the Snakemake input wiring.")
}
m_fitted <- ncol(Y_full)
add_diag("n_markers_fitted_before_variance_filter", m_fitted)

################################################################################
# 5. Zero-variance guard
################################################################################
col_sds <- apply(Y_full, 2, sd, na.rm = TRUE)
zero_var_idx <- which(col_sds == 0 | is.na(col_sds))
kept_idx <- setdiff(seq_len(m_fitted), zero_var_idx)
add_diag("n_zero_variance_dropped", length(zero_var_idx))
message(paste0("INFO: Dropped ", length(zero_var_idx),
               " zero-variance SNP(s) from the fit (climate-valid subset)"))
# Zero-variance SNPs are correctly excluded from m for Bonferroni purposes —
# a SNP untested in the climate-valid subset must not inflate the denominator.
# See docs/rda_research.md A2/A3; do not "fix" this back to the nominal m.
Y <- Y_full[, kept_idx, drop = FALSE]
m <- ncol(Y)
add_diag("n_markers_fitted", m)

################################################################################
# 6. Sample-count guard
################################################################################
min_n <- length(predictor_names) + N_COND_PCS + 2
if (n_samples <= min_n) {
    stop("RDA: only ", n_samples, " samples but need > ", min_n,
         " (n_predictors + condition_pcs + 2). Reduce condition_pcs or add samples.")
}

################################################################################
# 7. Fit the partial RDA (A8, A9, A12)
################################################################################
message("INFO: Fitting RDA")
Variables <- as.data.frame(Env)
rhs_terms <- predictor_names
if (!is.null(covariates)) {
    Variables <- cbind(Variables, as.data.frame(covariates))
    rhs_terms <- c(rhs_terms, paste0("Condition(", paste(colnames(covariates), collapse = "+"), ")"))
}
rda_formula <- as.formula(paste("Y ~", paste(rhs_terms, collapse = " + ")))
message(paste0("INFO: Formula: ", deparse(rda_formula)))
RDA_env <- rda(rda_formula, data = Variables, scale = TRUE)

################################################################################
# 8. Aliasing guard (empirically required — SIMDATA aliases bio_3 with
#    Condition(PC1+PC2+PC3) on 3 sites)
################################################################################
aliased_terms <- RDA_env$CCA$alias
if (!is.null(aliased_terms) && length(aliased_terms) > 0) {
    message(paste0("WARNING: RDA predictor(s) aliased (linearly dependent given ",
                   "the Condition() terms): ", paste(aliased_terms, collapse = ", ")))
    add_diag("aliased_terms", paste(aliased_terms, collapse = ","))
} else {
    add_diag("aliased_terms", "")
}

################################################################################
# 9. Constrained-axis count / hard floor (covRob needs >= 2 loading columns —
#    verified empirically: K=1 crashes with "'data' must have at least two columns")
################################################################################
K_max <- RDA_env$CCA$rank
add_diag("rda_axes_max", K_max)
if (K_max < 2) {
    stop("RDA: only ", K_max, " constrained axis/axes retained after Condition() — ",
         "rdadapt's covRob step needs >= 2 loading columns. Reduce condition_pcs ",
         "(currently ", N_COND_PCS, ") or add more predictors.")
}

################################################################################
# 10. RsquareAdj + VIF
################################################################################
r2adj <- RsquareAdj(RDA_env)$adj.r.squared
add_diag("adj_r_squared", r2adj)
vif_vals <- tryCatch(vif.cca(RDA_env), error = function(e) {
    message("WARNING: vif.cca() failed: ", e$message)
    setNames(rep(NA_real_, length(predictor_names)), predictor_names)
})
add_diag("max_vif", suppressWarnings(max(vif_vals, na.rm = TRUE)))
vif_flagged <- names(vif_vals)[!is.na(vif_vals) & vif_vals > 10]
add_diag("vif_flagged_predictors", paste(vif_flagged, collapse = ","))

################################################################################
# 11. ANOVA — full, by-axis, by-margin. All tryCatch-wrapped; by="margin" is
#     empirically degenerate on small/structured datasets (SIMDATA: Df=0,
#     F=Inf for every predictor) — never trust blindly, status column marks it.
################################################################################
set.seed(SEED)
anova_full <- tryCatch(
    anova.cca(RDA_env, permutations = PERMUTATIONS, parallel = max(1, CPU)),
    error = function(e) { message("WARNING: full-model anova.cca failed: ", e$message); NULL }
)
set.seed(SEED)
anova_axis <- tryCatch(
    anova.cca(RDA_env, by = "axis", permutations = PERMUTATIONS),
    error = function(e) { message("WARNING: by-axis anova.cca failed: ", e$message); NULL }
)
set.seed(SEED)
anova_margin <- tryCatch(
    anova.cca(RDA_env, by = "margin", permutations = PERMUTATIONS),
    error = function(e) { message("WARNING: by-margin anova.cca failed: ", e$message); NULL }
)

anova_row <- function(test_name, tbl) {
    if (is.null(tbl)) return(data.table())
    dt <- as.data.table(as.data.frame(tbl), keep.rownames = "term")
    dt[, test := test_name]
    # Standardize column names across the three anova.cca outputs (Df/Variance/
    # F/Pr(>F) -> df/variance/F/p). Only rename columns that are actually present.
    old_names <- c("Df", "Variance", "F", "Pr(>F)")
    new_names <- c("df", "variance", "F", "p")
    present <- old_names %in% colnames(dt)
    setnames(dt, old_names[present], new_names[present])
    has_df <- "df" %in% colnames(dt)
    dt[, status := {
        f_finite <- is.finite(F)
        ok <- if (has_df) (f_finite & df != 0) else f_finite
        ifelse(ok, "ok", "degenerate")
    }]
    dt
}
anova_dt <- rbindlist(list(
    anova_row("full",   anova_full),
    anova_row("axis",   anova_axis),
    anova_row("margin", anova_margin)
), fill = TRUE)
if (nrow(anova_dt) > 0) anova_dt[, model := "partial"]
n_margin_degenerate <- if (nrow(anova_dt) > 0) sum(anova_dt$test == "margin" & anova_dt$status == "degenerate") else 0L
add_diag("n_margin_rows_degenerate", n_margin_degenerate)

################################################################################
# 12. Choose K (A11)
################################################################################
axis_pvals <- if (!is.null(anova_axis) && "Pr(>F)" %in% colnames(as.data.frame(anova_axis))) {
    as.data.frame(anova_axis)[["Pr(>F)"]][seq_len(K_max)]
} else {
    rep(NA_real_, K_max)
}
K_floored <- FALSE
if (AXES == "auto") {
    K_sel <- sum(axis_pvals < AXIS_ALPHA, na.rm = TRUE)
    if (K_sel < 2) { K_floored <- TRUE; K_sel <- 2 }
    K_sel <- min(K_sel, K_max)
    K_selection <- "auto"
} else {
    K_sel <- as.integer(AXES)
    if (K_sel < 2) stop("RDA: axes must be >= 2 (parse-time validation should have caught this).")
    K_sel <- min(K_sel, K_max)
    K_selection <- "fixed"
}
message(paste0("INFO: Selected K=", K_sel, " constrained axes (K_max=", K_max,
               ", selection=", K_selection, ", floored=", K_floored, ")"))
add_diag("rda_axes", K_sel)
add_diag("K_selection", K_selection)
add_diag("K_floored", K_floored)

################################################################################
# 13. rdadapt on the PARTIAL (structure-corrected) fit — extracted to
#     scripts/R/lib/rdadapt.R so preGEA's RDA-setup block (Block 3) can share
#     the exact same candidate-calling math instead of duplicating it.
#     Verbatim port from Capblancq/RDA-genome-scan/Script_RDA.R (see that
#     file's header for the drop=FALSE deviation note).
################################################################################
message("INFO: Running rdadapt on partial fit (robust Mahalanobis distance over retained axes)")
rdadapt_partial <- rdadapt(RDA_env, K_sel)
K_partial <- K_sel
add_diag("gif_lambda", rdadapt_partial$gif_lambda)          # headline lambda = partial fit's (A8/A9)
add_diag("gif_lambda_partial", rdadapt_partial$gif_lambda)
add_diag("k_partial", K_partial)

################################################################################
# 14. B6 — second, UNCONSTRAINED rda() fit (same predictors, no Condition())
#     + intersection-style candidate combination. Capblancq & Forester 2021's
#     tutorial calls Reduce(intersect, ...) on two independently-thresholded
#     candidate lists (partial + unconstrained); the continuous-p-value
#     equivalent is an elementwise max — a SNP only stays maximally
#     significant if BOTH fits call it significant. K is PINNED to K_partial
#     rather than letting the unconstrained fit choose its own axis count via
#     anova.cca(by="axis"): two fits at different K sit on different
#     chi-squared df and are not comparable (this session's design decision —
#     see docs/rda_research.md B6). Skipped when condition_pcs=0: the partial
#     fit IS already unconstrained, so a second identical fit is dead work
#     (Code Structure Rule 3).
################################################################################
unc_error_msg <- NULL
K_unconstrained <- NA_integer_
k_pin_capped <- FALSE
rdadapt_unconstrained <- NULL
anova_dt_unc <- data.table()

if (N_COND_PCS == 0) {
    unconstrained_fit_status <- "skipped_condition_pcs_zero (partial fit already unconstrained)"
} else {
    message("INFO: Fitting unconstrained RDA (B6 — no Condition() term)")
    RDA_unc <- tryCatch(
        rda(as.formula(paste("Y ~", paste(predictor_names, collapse = " + "))),
            data = Variables, scale = TRUE),
        error = function(e) { unc_error_msg <<- conditionMessage(e); NULL }
    )
    if (is.null(RDA_unc)) {
        unconstrained_fit_status <- paste0("failed: ", unc_error_msg)
    } else {
        K_max_unc <- RDA_unc$CCA$rank
        K_unconstrained <- min(K_partial, K_max_unc)
        k_pin_capped <- K_unconstrained < K_partial
        if (K_unconstrained < 2) {
            unconstrained_fit_status <- paste0("failed: only ", K_max_unc,
                                               " constrained axis/axes retained (< 2)")
            K_unconstrained <- NA_integer_
        } else {
            message("INFO: Running rdadapt on unconstrained fit (K=", K_unconstrained,
                   ", pinned from partial K=", K_partial, ")")
            rdadapt_unconstrained <- rdadapt(RDA_unc, K_unconstrained)
            unconstrained_fit_status <- "fitted"

            # Same 3-way anova.cca as section 11, on the unconstrained fit —
            # tagged model="unconstrained" below so RDA_anova.tsv carries both
            # models' diagnostics without a second output file.
            set.seed(SEED)
            anova_full_unc <- tryCatch(
                anova.cca(RDA_unc, permutations = PERMUTATIONS, parallel = max(1, CPU)),
                error = function(e) { message("WARNING: unconstrained full-model anova.cca failed: ", e$message); NULL }
            )
            set.seed(SEED)
            anova_axis_unc <- tryCatch(
                anova.cca(RDA_unc, by = "axis", permutations = PERMUTATIONS),
                error = function(e) { message("WARNING: unconstrained by-axis anova.cca failed: ", e$message); NULL }
            )
            set.seed(SEED)
            anova_margin_unc <- tryCatch(
                anova.cca(RDA_unc, by = "margin", permutations = PERMUTATIONS),
                error = function(e) { message("WARNING: unconstrained by-margin anova.cca failed: ", e$message); NULL }
            )
            anova_dt_unc <- rbindlist(list(
                anova_row("full",   anova_full_unc),
                anova_row("axis",   anova_axis_unc),
                anova_row("margin", anova_margin_unc)
            ), fill = TRUE)
            if (nrow(anova_dt_unc) > 0) anova_dt_unc[, model := "unconstrained"]
        }
    }
}
add_diag("unconstrained_fit_status", unconstrained_fit_status)
add_diag("k_unconstrained", K_unconstrained)
add_diag("k_pin_capped", k_pin_capped)
add_diag("gif_lambda_unconstrained",
         if (!is.null(rdadapt_unconstrained)) rdadapt_unconstrained$gif_lambda else NA_real_)

if (!is.null(rdadapt_unconstrained)) {
    p_combined <- pmax(rdadapt_partial$p.values, rdadapt_unconstrained$p.values)
    # Report whichever fit's distance is "binding" (drove the higher, less-
    # significant p-value per SNP) — the pmax winner, not always the partial fit.
    driving_unc <- rdadapt_unconstrained$p.values >= rdadapt_partial$p.values
    mahalanobis_combined <- ifelse(driving_unc, rdadapt_unconstrained$distance, rdadapt_partial$distance)
} else {
    p_combined <- rdadapt_partial$p.values
    mahalanobis_combined <- rdadapt_partial$distance
}
q_combined_res <- qvalue_with_fallback(p_combined)
add_diag("qvalue_method", q_combined_res$method)

# rdadapt_res is now the COMBINED result — every downstream section (15-21)
# that reads p.values/q.values/distance operates on the B6 combination, not a
# single fit. p.values_partial/p.values_unconstrained/distance_partial/
# distance_unconstrained carry the per-fit values through to the candidates
# table (section 17) so a reader can see which fit drove each SNP's call.
rdadapt_res <- list(
    p.values = p_combined, q.values = q_combined_res$qvalues, distance = mahalanobis_combined,
    gif_lambda = rdadapt_partial$gif_lambda, qvalue_method = q_combined_res$method,
    p.values_partial = rdadapt_partial$p.values,
    p.values_unconstrained = if (!is.null(rdadapt_unconstrained)) rdadapt_unconstrained$p.values else rep(NA_real_, m),
    distance_partial = rdadapt_partial$distance,
    distance_unconstrained = if (!is.null(rdadapt_unconstrained)) rdadapt_unconstrained$distance else rep(NA_real_, m)
)

################################################################################
# 15. Marker-envelope statement (A7, A15) — always written, regardless of m.
################################################################################
envelope_status <- if (m >= 1e4 && m <= 1e5) "within_validated" else if (m < 1e4) "below_validated" else "ABOVE_VALIDATED"
envelope_note <- paste0(
    "This run uses m = ", m, " markers. All published RDA/rdadapt validation is ",
    "at 10^4-10^5 markers (Forester et al. 2018: 10,000 simulated / 42,587 wolf ",
    "SNPs; Capblancq & Forester 2021: tens of thousands post-MAF). No RDA GEA ",
    "validation exists at WGS density in any mating system. The chi-squared ",
    "calibration and Bonferroni threshold are extrapolations beyond the ",
    "validated envelope when status != 'within_validated'."
)
add_diag("marker_envelope_status", envelope_status)
add_diag("marker_envelope_note", envelope_note)
message(paste0("INFO: marker_envelope_status=", envelope_status))

################################################################################
# 16-17. Candidate calling + per-predictor correlation / assigned predictor (A18)
################################################################################
bonf_thresh_001 <- 0.01 / m
bonf_thresh_005 <- 0.05 / m
n_bonf_001 <- sum(rdadapt_res$p.values < bonf_thresh_001, na.rm = TRUE)
n_bonf_005 <- sum(rdadapt_res$p.values < bonf_thresh_005, na.rm = TRUE)
n_qval_005 <- sum(rdadapt_res$q.values < 0.05, na.rm = TRUE)
add_diag("n_candidates_bonf_0.01", n_bonf_001)
add_diag("n_candidates_bonf_0.05", n_bonf_005)
add_diag("n_candidates_q_0.05", n_qval_005)

# p-value histogram flatness (A4 diagnostic — reported, never enforced): chi-sq
# goodness-of-fit of a 100-bin histogram against uniform, EXCLUDING the first
# bin (a spike there is the desired signal, not a departure from flat).
hist_counts <- hist(rdadapt_res$p.values, breaks = seq(0, 1, length.out = 101), plot = FALSE)$counts
expected_per_bin <- length(rdadapt_res$p.values) / 100
flatness_chisq <- sum((hist_counts[-1] - expected_per_bin)^2 / expected_per_bin)
add_diag("pval_hist_flatness_chisq", flatness_chisq)

################################################################################
# 16b. Sensitivity ladder (DIAGNOSTIC ONLY — not the applied rule).
#      docs/rda_research.md A.4(1) records the literature as unresolved:
#      Capblancq et al. 2018 use q<0.1; Capblancq & Forester 2021 use 0.01/m.
#      These four counts say "what each published rule would give on this run";
#      the rule actually applied is candidate_rule below.
################################################################################
n_qval_010 <- sum(rdadapt_res$q.values < 0.10, na.rm = TRUE)
add_diag("n_candidates_q_0.10", n_qval_010)
add_diag("sensitivity_ladder_note",
         paste0("n_candidates_bonf_* / n_candidates_q_* are DIAGNOSTIC counts under ",
                "published alternative rules, NOT the applied rule. Applied rule = ",
                "candidate_rule. Literature is unresolved (2018: q<0.1; 2021: ",
                "Bonferroni 0.01/m) — docs/rda_research.md A.4(1)."))

################################################################################
# 16c. Applied candidate rule — SHARED with find_sig_snps.R.
#      compute_pval_threshold() is the same function find_sig_snps.R calls on
#      OUT_PVALUES, so the candidates table and *_sig_snps_{adjust}.tsv are the
#      same SNP set by construction. That equality is preserved by B6 (section
#      14): p_combined is computed ONCE, above, and both this section and
#      section 18's OUT_PVALUES write read the same rdadapt_res$p.values
#      (=p_combined) — there is no second, independent combination step. The
#      equality further rests on the pre-existing invariant: section 18 below
#      re-expands variance-dropped SNPs as NA, and compute_pval_threshold()
#      NA-drops before counting, so its n_tested equals m here. Never write 0
#      instead of NA in section 18.
################################################################################
cand_res <- compute_pval_threshold(rdadapt_res$p.values, CAND_ADJUST, CAND_VALUE)

add_diag("candidate_adjust",           CAND_ADJUST)
add_diag("candidate_threshold_value",  CAND_VALUE)
add_diag("candidate_rule",             paste0(CAND_ADJUST, "_", CAND_VALUE))
add_diag("candidate_threshold",        cand_res$threshold)
add_diag("candidate_threshold_status", cand_res$status)
add_diag("candidate_n_tested",         cand_res$n_tested)
add_diag("candidate_n_na_dropped",     cand_res$n_na_dropped)

# The q-value COLUMN in the candidates table comes from section 14's
# qvalue_with_fallback(p_combined) call (storey -> lambda=0 -> BH), applied to
# the COMBINED p-values, NOT from pval_threshold.R's bare qvalue()
# (max_pvalue_fdr, no fallback). Under adjust="qval" the applied cutoff is
# computed via the shared function (so it agrees with find_sig_snps.R), but it
# can be a DIFFERENT estimate than the q_value column whenever the fallback
# chain fires — record the divergence rather than silently reconciling it.
# See docs/rda_research.md A.4.
add_diag("candidate_qvalue_engine",
         if (identical(CAND_ADJUST, "qval"))
             paste0("pval_threshold::max_pvalue_fdr — bare qvalue(), NO fallback chain. ",
                    "The q_value COLUMN in the candidates table is qvalue_with_fallback's ",
                    "(method=", rdadapt_res$qvalue_method, ", applied to the B6-combined ",
                    "p-values) and may use a different estimator.")
         else "n/a (adjust != qval)")

if (!identical(cand_res$status, "ok")) {
    message("WARNING: RDA candidate threshold unavailable (status=", cand_res$status,
            ") for rule ", CAND_ADJUST, " ", CAND_VALUE,
            " on ", cand_res$n_tested, " tests. Writing header-only candidates table. ",
            "find_sig_snps.R will reach the same conclusion on the same p-values.")
    is_candidate <- rep(FALSE, length(rdadapt_res$p.values))
} else {
    # Inclusive '<=' and NA->FALSE, byte-identical to sig_snps.R's mask. The
    # cutoff compute_pval_threshold() returns is a member of the call set under
    # 'qval'/'top' (see that function's header), so `<` would drop the boundary
    # SNP and its ties here and break the set equality asserted in 16c.
    is_candidate <- rdadapt_res$p.values <= cand_res$threshold
    is_candidate[is.na(is_candidate)] <- FALSE
}
n_candidates_total <- sum(is_candidate)
add_diag("n_candidates_total", n_candidates_total)

# SNP identifiers, in the SAME order as Y's original columns (before the
# zero-variance drop) — read once, positionally matched to genotype columns,
# exactly as lfmm.R does.
vcfsnp_dt <- fread(VCFSNP, sep = ' ', header = FALSE) %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chr = as.character(chr), SNPID = paste0(chr, ':', pos)) %>%
    dplyr::select(SNPID, chr, pos)
if (nrow(vcfsnp_dt) != m_fitted) {
    stop("RDA: vcfsnp row count (", nrow(vcfsnp_dt), ") does not match genotype ",
         "column count (", m_fitted, ") before variance filtering.")
}
vcfsnp_kept <- vcfsnp_dt[kept_idx, ]

candidate_idx <- which(is_candidate)
if (length(candidate_idx) > MAX_CANDIDATES_WRITTEN) {
    message(paste0("WARNING: ", length(candidate_idx), " candidates found; capping ",
                   "the written side table at ", MAX_CANDIDATES_WRITTEN,
                   " by Mahalanobis distance (see rank column)."))
}
add_diag("n_candidates_written", min(length(candidate_idx), MAX_CANDIDATES_WRITTEN))

if (length(candidate_idx) > 0) {
    cand_loadings <- RDA_env$CCA$v[candidate_idx, 1:K_sel, drop = FALSE]
    colnames(cand_loadings) <- paste0("RDA", seq_len(K_sel))

    # Per-predictor correlation, candidates only (A18) — includes aliased
    # predictors (aliasing affects the ordination fit, not marginal correlation;
    # verified fine on SIMDATA's aliased bio_3).
    cor_mat <- sapply(predictor_names, function(p) {
        cor(Y[, candidate_idx, drop = FALSE], Variables[[p]], use = "pairwise.complete.obs")
    })
    if (is.null(dim(cor_mat))) cor_mat <- matrix(cor_mat, ncol = length(predictor_names),
                                                 dimnames = list(NULL, predictor_names))
    colnames(cor_mat) <- paste0("cor_", predictor_names)
    assigned_predictor <- predictor_names[apply(abs(cor_mat), 1, which.max)]
    assigned_cor <- cor_mat[cbind(seq_len(nrow(cor_mat)), apply(abs(cor_mat), 1, which.max))]

    candidates_dt <- data.table(vcfsnp_kept[candidate_idx, ],
                                mahalanobis    = rdadapt_res$distance[candidate_idx],
                                p_value        = rdadapt_res$p.values[candidate_idx],
                                q_value        = rdadapt_res$q.values[candidate_idx],
                                pval_threshold = cand_res$threshold,
                                # B6 split columns — which fit drove the combined call for this SNP.
                                mahalanobis_partial       = rdadapt_res$distance_partial[candidate_idx],
                                mahalanobis_unconstrained = rdadapt_res$distance_unconstrained[candidate_idx],
                                p_value_partial           = rdadapt_res$p.values_partial[candidate_idx],
                                p_value_unconstrained     = rdadapt_res$p.values_unconstrained[candidate_idx])
    candidates_dt <- cbind(candidates_dt, as.data.table(cand_loadings), as.data.table(cor_mat))
    candidates_dt[, assigned_predictor := assigned_predictor]
    candidates_dt[, assigned_cor := assigned_cor]
    candidates_dt <- candidates_dt[order(-mahalanobis)]
    candidates_dt[, rank := .I]
    if (nrow(candidates_dt) > MAX_CANDIDATES_WRITTEN) {
        candidates_dt <- candidates_dt[seq_len(MAX_CANDIDATES_WRITTEN)]
    }
} else {
    message("INFO: No candidate SNPs at ", CAND_ADJUST, " ", CAND_VALUE,
            " (raw p cutoff ", format(cand_res$threshold, digits = 3),
            ", status=", cand_res$status, ") — writing header-only candidates table")
    candidates_dt <- data.table(
        SNPID = character(), chr = character(), pos = integer(),
        mahalanobis = numeric(), p_value = numeric(), q_value = numeric(),
        pval_threshold = numeric(),
        mahalanobis_partial = numeric(), mahalanobis_unconstrained = numeric(),
        p_value_partial = numeric(), p_value_unconstrained = numeric()
    )
    for (i in seq_len(K_sel)) candidates_dt[[paste0("RDA", i)]] <- numeric()
    for (p in predictor_names) candidates_dt[[paste0("cor_", p)]] <- numeric()
    candidates_dt[, assigned_predictor := character()]
    candidates_dt[, assigned_cor := numeric()]
    candidates_dt[, rank := integer()]
}

################################################################################
# 18. Write wide p-value table — ALL m_fitted rows (dropped SNPs re-expanded as NA)
################################################################################
message("INFO: Writing p-value table: ", OUT_PVALUES)
full_pvals <- rep(NA_real_, m_fitted)
full_pvals[kept_idx] <- rdadapt_res$p.values
pval_out <- vcfsnp_dt
pval_out[, climate_multivariate := full_pvals]
setorder(pval_out, chr, pos)
fwrite(pval_out, OUT_PVALUES, sep = '\t', quote = FALSE)

################################################################################
# 19. Write candidates side table
################################################################################
message("INFO: Writing candidates table: ", OUT_CANDIDATES)
fwrite(candidates_dt, OUT_CANDIDATES, sep = '\t', quote = FALSE)

################################################################################
# 20. Write diagnostics (long key/value) + anova (long) tables
################################################################################
message("INFO: Writing diagnostics table: ", OUT_DIAGNOSTICS)
# NOTE: data.table()'s `key` argument is reserved (sets the table's key, not a
# column name) — build positionally and rename to avoid a silent misparse.
diag_dt <- data.table(names(diag), sapply(diag, function(x) paste(x, collapse = ",")))
setnames(diag_dt, c("key", "value"))
fwrite(diag_dt, OUT_DIAGNOSTICS, sep = '\t', quote = FALSE)

message("INFO: Writing anova table: ", OUT_ANOVA)
anova_dt <- rbindlist(list(anova_dt, anova_dt_unc), fill = TRUE)  # B6 — partial + unconstrained rows
fwrite(anova_dt, OUT_ANOVA, sep = '\t', quote = FALSE)

################################################################################
# 21. Plots (inline — the fitted cca object is needed by every panel; RDA has
#     no persisted model artifact to serialize/re-read, unlike Gradient Forest)
################################################################################
message("INFO: Generating plots in ", PLOT_DIR)

# --- Screeplot: constrained eigenvalues, retained vs dropped axes ---
eig_dt <- data.table(axis = seq_along(RDA_env$CCA$eig), eigenvalue = RDA_env$CCA$eig)
eig_dt[, retained := axis <= K_sel]
axis_p_dt <- data.table(axis = seq_len(K_max), p = axis_pvals[seq_len(K_max)])
eig_dt <- merge(eig_dt, axis_p_dt, by = "axis", all.x = TRUE)
g_scree <- ggplot(eig_dt, aes(x = factor(axis), y = eigenvalue, fill = retained)) +
    geom_col() +
    geom_text(aes(label = ifelse(!is.na(p), sprintf("p=%.3g", p), "")),
              vjust = -0.3, size = 3, color = CLINEGO_COL$fg) +
    scale_fill_manual(values = c(`TRUE` = CLINEGO_RETAINED, `FALSE` = CLINEGO_REMOVED), guide = "none") +
    labs(x = "Constrained axis", y = "Eigenvalue",
         title = "RDA constrained eigenvalues",
         subtitle = paste0("K=", K_sel, " (", K_selection, "), K_max=", K_max)) +
    theme_clinego()
ggsave(O_scree <- file.path(PLOT_DIR, paste0("rda_screeplot_K", K_BEST, ".png")), g_scree,
       width = 7, height = 5, dpi = 300)
ggsave(sub("\\.png$", ".svg", O_scree), g_scree, width = 7, height = 5,
       device = svglite::svglite, bg = "transparent")

# --- p-value histogram (binned by construction — Rule 6 compliant at any m) ---
g_hist <- ggplot(data.table(p = rdadapt_res$p.values), aes(x = p)) +
    geom_histogram(breaks = seq(0, 1, length.out = 101), fill = CLINEGO_NEUTRAL, color = "black") +
    geom_hline(yintercept = expected_per_bin, linetype = "dashed", color = CLINEGO_THRESHOLD) +
    labs(x = "p-value", y = "Count",
         title = "RDA candidate p-value histogram",
         subtitle = sprintf("flatness chi-sq=%.1f (excl. first bin) — expect flat + spike at 0 (A4)",
                            flatness_chisq)) +
    theme_clinego()
ggsave(O_hist <- file.path(PLOT_DIR, paste0("rda_pvalue_histogram_K", K_BEST, ".png")), g_hist,
       width = 7, height = 5, dpi = 300)
ggsave(sub("\\.png$", ".svg", O_hist), g_hist, width = 7, height = 5,
       device = svglite::svglite, bg = "transparent")

# --- SNP loadings biplot (RDA1 x RDA2): geom_bin2d, NOT geom_hex (hexbin not
#     installed) and NOT geom_point for the SNP cloud (Code Structure Rule 6) ---
if (K_sel >= 2) {
    loadings_all <- as.data.table(RDA_env$CCA$v[, 1:2, drop = FALSE])
    setnames(loadings_all, c("RDA1", "RDA2"))
    biplot_scores <- as.data.frame(RDA_env$CCA$biplot[, 1:2, drop = FALSE])
    colnames(biplot_scores) <- c("RDA1", "RDA2")
    biplot_scores$predictor <- rownames(RDA_env$CCA$biplot)
    arrow_scale <- 0.8 * max(abs(loadings_all$RDA1), abs(loadings_all$RDA2)) /
        max(abs(biplot_scores$RDA1), abs(biplot_scores$RDA2))

    g_biplot <- ggplot(loadings_all, aes(x = RDA1, y = RDA2)) +
        geom_bin2d(bins = 80) +
        scale_fill_viridis_c(trans = "log10", name = "SNP density") +
        geom_segment(data = biplot_scores,
                     aes(x = 0, y = 0, xend = RDA1 * arrow_scale, yend = RDA2 * arrow_scale),
                     inherit.aes = FALSE, color = CLINEGO_THRESHOLD,
                     arrow = arrow(length = unit(0.2, "cm"))) +
        ggrepel::geom_text_repel(data = biplot_scores,
                                 aes(x = RDA1 * arrow_scale, y = RDA2 * arrow_scale, label = predictor),
                                 inherit.aes = FALSE, color = CLINEGO_THRESHOLD, size = 3) +
        labs(x = "RDA1", y = "RDA2", title = "RDA SNP loadings biplot") +
        theme_clinego()
    if (nrow(candidates_dt) > 0) {
        cand_plot_dt <- candidates_dt[, .(RDA1 = RDA1, RDA2 = RDA2, assigned_predictor)]
        g_biplot <- g_biplot +
            scattermore::geom_scattermore(data = cand_plot_dt, aes(x = RDA1, y = RDA2, color = assigned_predictor),
                                          inherit.aes = FALSE, pointsize = 4) +
            scale_color_clinego(name = "Assigned predictor")
    }
} else {
    # Unreachable given the K>=2 hard floor above (step 9) — kept as a defensive
    # empty-state placeholder, matching plot_qc_processing.R's convention.
    g_biplot <- ggplot() +
        annotate("text", x = 0, y = 0, label = "K < 2 — biplot unavailable") +
        theme_void()
}
ggsave(O_biplot <- file.path(PLOT_DIR, paste0("rda_loadings_biplot_K", K_BEST, ".png")), g_biplot,
       width = 7, height = 6, dpi = 300)
ggsave(sub("\\.png$", ".svg", O_biplot), g_biplot, width = 7, height = 6,
       device = svglite::svglite, bg = "transparent")

message("INFO: RDA genome scan complete.")
message(paste0("INFO: m=", m, " K=", K_sel, " (max=", K_max, ") gif=",
               round(rdadapt_res$gif_lambda, 4),
               " | B6 unconstrained_fit=", unconstrained_fit_status,
               " | APPLIED ", CAND_ADJUST, "_", CAND_VALUE, " -> ", n_candidates_total,
               " candidates | ladder: bonf0.01=", n_bonf_001, " bonf0.05=", n_bonf_005,
               " q0.05=", n_qval_005, " q0.10=", n_qval_010))
