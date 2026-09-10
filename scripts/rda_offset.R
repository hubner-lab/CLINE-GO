#!/usr/bin/env Rscript
# RDA genetic offset — Redundancy Analysis genomic offset (Capblancq & Forester
# 2021, "Swiss Army Knife for landscape genomics").
#
# NOT a reuse of the GEA-scan RDA object (scripts/rda.R). "Compute once, reuse
# downstream" does not apply between the GEA scan and the offset — the two
# fits differ on TWO axes simultaneously (docs/rda_research.md Part B.0):
#   GEA scan  : ALL SNPs,                partial (Condition(PCs))
#   Offset    : CANDIDATE SNPs only,     no Condition() by default
# A second rda() fit is scientifically required, not an optimization to skip.
# Row IDs below (B1, B2, ...) reference docs/rda_research.md Part B's decision
# table — read that table before touching the weighting/projection math.

library(data.table)
library(dplyr)
library(vegan)
library(terra)
library(stringr)
library(ggplot2)

source("/pipeline/scripts/R/utils/theme_clinego.R")
source("/pipeline/scripts/R/utils/emmax_core.R")  # load_pca_covariates() — condition_pcs > 0 only

args <- commandArgs(trailingOnly = TRUE)
################################################################################
# Arg order (26 positional):
LFMM_IMP_FULL  <- args[1]   # full imputed genotype matrix (.lfmm, rows=samples cols=SNPs)
VCFSNP         <- args[2]   # *.vcfsnp for the full marker set (space-sep, no header, V1=chr V2=pos)
REMOVED        <- args[3]   # removed SNPs table (no header)
CANDIDATE_SNPS <- args[4]   # curated snp_set_file(run_label): selected_snps.tsv, SNPID column
ENV_SITE_PRES  <- args[5]   # climate_present_site.tsv (per-sample rows, unscaled — B11 needs raw values)
ENV_SITE_FUT   <- args[6] %>% str_split(',') %>% unlist   # climate_future_..._site.tsv, ONE PER SCENARIO
ENV_ALL_PRES   <- args[7]   # climate_present_all.tsv (per raster-cell rows, has ID column)
ENV_ALL_FUT    <- args[8] %>% str_split(',') %>% unlist   # climate_future_..._all.tsv, ONE PER SCENARIO
PRES_RASTER    <- args[9]   # present climate raster (.tif, spatial template)
SAMPLES        <- args[10]  # metadata_climate_valid.tsv (site, sample, latitude, longitude, ...)
PREDICTORS     <- args[11] %>% str_split(',') %>% unlist
AXES           <- args[12]  # "auto" or an integer >= 2
AXIS_ALPHA     <- as.numeric(args[13])
PERMUTATIONS   <- as.numeric(args[14])
CONDITION_PCS  <- as.integer(args[15])   # 0 = canonical (B7); >0 = documented deviation
PCA            <- args[16]  # LEA .projections (LD-pruned) — only read if CONDITION_PCS > 0
SAMPLES_ORDER  <- args[17]  # full-cohort sample order — only read if CONDITION_PCS > 0
CLIMATE_VALID  <- args[18]  # plink --keep list (FID IID) — only read if CONDITION_PCS > 0
SEED           <- as.numeric(args[19])
CPU            <- as.numeric(args[20])
OUT_SITE_TSV   <- args[21] %>% str_split(',') %>% unlist   # ONE PER SCENARIO
OUT_MAP_TSV    <- args[22] %>% str_split(',') %>% unlist   # ONE PER SCENARIO
OUT_RASTER_TIF <- args[23] %>% str_split(',') %>% unlist   # ONE PER SCENARIO
OUT_IMPORTANCE <- args[24]  # scenario-free — describes the fit, not the projection
OUT_DIAGNOSTICS<- args[25]  # scenario-free
PLOT_DIR       <- args[26]  # scenario-free
SCENARIOS      <- if (length(args) >= 27 && nzchar(args[27])) {
    args[27] %>% str_split(',') %>% unlist
} else {
    paste0('scenario_', seq_along(ENV_SITE_FUT))
}
################################################################################
# The RDA is fitted on PRESENT climate only; a future scenario enters solely at
# the projection step (section 9). So the fit, its axis selection, its weights
# and its diagnostics are computed ONCE and every scenario is projected from
# them — which is the whole point of accepting scenario lists here.
N_SCEN <- length(ENV_SITE_FUT)
stopifnot(length(ENV_ALL_FUT)    == N_SCEN,
          length(OUT_SITE_TSV)   == N_SCEN,
          length(OUT_MAP_TSV)    == N_SCEN,
          length(OUT_RASTER_TIF) == N_SCEN,
          length(SCENARIOS)      == N_SCEN)

for (p in unique(dirname(c(OUT_SITE_TSV, OUT_MAP_TSV, OUT_RASTER_TIF, OUT_DIAGNOSTICS)))) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(SEED)
message('INFO: RDA genetic offset (second, candidate-only RDA fit)')
message(paste0('INFO: axes=', AXES, '  condition_pcs=', CONDITION_PCS, '  snp_set=', CANDIDATE_SNPS))

# diag: accumulates key/value rows written to OUT_DIAGNOSTICS — same
# anti-silent-substitution discipline as scripts/rda.R.
diag <- list()
add_diag <- function(key, value) diag[[key]] <<- value
add_diag("snp_set_name", basename(dirname(CANDIDATE_SNPS)))
add_diag("predictors", paste(PREDICTORS, collapse = ","))
add_diag("axes_requested", AXES)
add_diag("axis_alpha", AXIS_ALPHA)
add_diag("permutations", PERMUTATIONS)
add_diag("condition_pcs", CONDITION_PCS)
add_diag("seed", SEED)
add_diag("b6_deviation_note",
         paste0("Candidate set = curated run_label SNP set (snp_set_file), NOT the ",
                "canonical B6 intersection of partial + unconstrained rdadapt outliers. ",
                "Justified by B5 + Lind et al. 2025 (marker-set choice ~0.95% of forecast ",
                "variance) — see docs/rda_research.md Part B."))
add_diag("interpretation_note",
         paste0("Report as a RELATIVE SITE RANKING only (B20); unvalidated absent ",
                "common-garden data (B21). Not comparable in absolute magnitude across ",
                "runs/datasets."))

################################################################################
# 1. Load full imputed genotype matrix + SNP identifiers
################################################################################
message('INFO: Loading full imputed LFMM matrix')
lfmm_imp <- as.matrix(fread(LFMM_IMP_FULL))
message(paste0('INFO: LFMM matrix: ', nrow(lfmm_imp), ' samples x ', ncol(lfmm_imp), ' SNPs'))
if (max(lfmm_imp, na.rm = TRUE) > 2) {
    stop("RDA offset: genotype matrix contains values > 2 — this looks like an ",
         "UNIMPUTED LFMM file (LEA writes 9 for missing). Requires an imputed matrix ",
         "(*_imp*). Check the Snakemake input wiring.")
}

vcfsnp_dt <- fread(VCFSNP, header = FALSE) %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chr = as.character(chr), chrpos = paste0(chr, ':', pos))
vcfsnp_full <- vcfsnp_dt$chrpos

if (ncol(lfmm_imp) != nrow(vcfsnp_dt)) {
    stop(sprintf(
        paste0('FATAL: Column count mismatch - lfmm_imp has %d columns but vcfsnp has %d rows. ',
               'Imputation must preserve SNP order.'),
        ncol(lfmm_imp), nrow(vcfsnp_dt)
    ))
}
if (!is.null(colnames(lfmm_imp)) && colnames(lfmm_imp)[1] != '' &&
        !startsWith(colnames(lfmm_imp)[1], 'V')) {
    expected <- vcfsnp_full[1]
    actual   <- colnames(lfmm_imp)[1]
    if (!grepl(sub('.*:', '', expected), actual, fixed = TRUE)) {
        stop(sprintf(
            'FATAL: First column name "%s" does not match first SNP "%s". Column order mismatch.',
            actual, expected
        ))
    }
}
message('INFO: Alignment assertion passed')

removed_dt <- fread(REMOVED, header = FALSE)
removed <- if (!is.null(removed_dt) && nrow(removed_dt) > 0) {
    removed_dt %>% dplyr::select(V1, V2) %>% setNames(c('chr', 'pos')) %>%
        dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>% .$chrpos
} else character(0)

################################################################################
# 2. Candidate SNPs — subset the matrix (B5, B6: unlike geometric_offset.R,
#    which must never subset because LEA needs the full matrix for latent-
#    factor estimation; RDA has no such requirement).
################################################################################
sig_dt   <- fread(CANDIDATE_SNPS, header = TRUE, colClasses = c('SNPID' = 'character'))
sig_snps <- sig_dt$SNPID
n_candidates_requested <- length(sig_snps)
add_diag("n_candidates_requested", n_candidates_requested)

candidate_idx <- which((!vcfsnp_full %in% removed) & (vcfsnp_full %in% sig_snps))
n_candidates_matched <- length(candidate_idx)
add_diag("n_candidates_matched", n_candidates_matched)
message(paste0('INFO: Candidate loci: ', n_candidates_matched, ' of ', n_candidates_requested,
               ' requested (', ncol(lfmm_imp), ' total SNPs)'))
if (n_candidates_matched == 0) {
    stop('FATAL: No candidate loci found. Check SNPID format in ', CANDIDATE_SNPS,
         ' (expected chr:pos) vs vcfsnp file.')
}

Y_cand <- lfmm_imp[, candidate_idx, drop = FALSE]

# Zero-variance guard on the candidate subset (rda.R precedent, section 5)
col_sds <- apply(Y_cand, 2, sd, na.rm = TRUE)
zero_var_idx <- which(col_sds == 0 | is.na(col_sds))
add_diag("n_zero_variance_dropped", length(zero_var_idx))
if (length(zero_var_idx) > 0) {
    Y_cand <- Y_cand[, setdiff(seq_len(ncol(Y_cand)), zero_var_idx), drop = FALSE]
}
message(paste0('INFO: ', ncol(Y_cand), ' candidate SNP(s) after zero-variance drop'))
if (ncol(Y_cand) < 3) {
    stop('FATAL: Only ', ncol(Y_cand), ' candidate SNP(s) survive zero-variance filtering — ',
         'too few for a stable RDA fit. Curate a larger SNP set (', CANDIDATE_SNPS, ') from ',
         'the GEA tab before running rda_offset.')
}

################################################################################
# 3. Climate predictors — B11: standardize ONCE from the present site matrix,
#    capture the constants, apply identically to present AND future.
################################################################################
env_site_pres_raw <- fread(ENV_SITE_PRES) %>% dplyr::select(!!PREDICTORS)
env_all_pres  <- fread(ENV_ALL_PRES) %>% dplyr::select(ID, !!PREDICTORS)

n_samples <- nrow(env_site_pres_raw)
add_diag("n_samples", n_samples)
add_diag("n_scenarios", N_SCEN)
if (nrow(Y_cand) != n_samples) {
    stop(paste0('FATAL: Sample-count mismatch -- env_site_pres has ', n_samples,
                ' rows, genotype matrix has ', nrow(Y_cand), ' rows. These must match ',
                '(both subset to the same climate-valid sample set). Re-run ',
                'subset_lfmm_climate / filter_climate_valid_samples and check for a ',
                'stale intermediate file.'))
}

Env_scaled <- scale(as.matrix(env_site_pres_raw), center = TRUE, scale = TRUE)
center_env <- attr(Env_scaled, 'scaled:center')
scale_env  <- attr(Env_scaled, 'scaled:scale')
Variables  <- as.data.frame(Env_scaled)
colnames(Variables) <- PREDICTORS

################################################################################
# 4. Structure-correction PCs (B7 — off by default; condition_pcs > 0 is a
#    documented deviation from the canonical uncorrected-offset construction)
################################################################################
covariates <- NULL
if (CONDITION_PCS > 0) {
    climate_valid_ids <- fread(CLIMATE_VALID, header = FALSE, colClasses = "character")$V2
    covariates <- load_pca_covariates(PCA, CONDITION_PCS, SAMPLES_ORDER, climate_valid_ids)
    if (!is.null(covariates) && nrow(covariates) != n_samples) {
        stop("RDA offset: PCA covariate row count (", nrow(covariates),
             ") does not match sample count (", n_samples, ")")
    }
    add_diag("condition_deviation",
             paste0("condition_pcs=", CONDITION_PCS, " > 0 — DEVIATION from B7 (canonical ",
                    "offset RDA carries NO Condition()). Conditioning removes exactly the ",
                    "climate-correlated structure variance the adaptive index needs to ",
                    "project (docs/rda_research.md B7)."))
}

################################################################################
# 5. Sample-count guard (mirrors rda.R section 6)
################################################################################
min_n <- length(PREDICTORS) + CONDITION_PCS + 2
if (n_samples <= min_n) {
    stop("RDA offset: only ", n_samples, " samples but need > ", min_n,
         " (n_predictors + condition_pcs + 2). Reduce condition_pcs or add samples.")
}

################################################################################
# 6. Fit the offset RDA — CANDIDATE SNPs, no Condition() unless CONDITION_PCS>0 (B7, B8)
################################################################################
rhs_terms <- PREDICTORS
if (!is.null(covariates)) {
    Variables <- cbind(Variables, as.data.frame(covariates))
    rhs_terms <- c(rhs_terms, paste0("Condition(", paste(colnames(covariates), collapse = "+"), ")"))
}
rda_formula <- as.formula(paste("Y_cand ~", paste(rhs_terms, collapse = " + ")))
message(paste0("INFO: Formula: ", deparse(rda_formula)))
RDA_off <- rda(rda_formula, data = Variables, scale = TRUE)

aliased_terms <- RDA_off$CCA$alias
if (!is.null(aliased_terms) && length(aliased_terms) > 0) {
    message(paste0("WARNING: RDA offset predictor(s) aliased: ", paste(aliased_terms, collapse = ", ")))
    add_diag("aliased_terms", paste(aliased_terms, collapse = ","))
} else {
    add_diag("aliased_terms", "")
}

r2adj <- tryCatch(RsquareAdj(RDA_off)$adj.r.squared, error = function(e) NA_real_)
add_diag("adj_r_squared", r2adj)

################################################################################
# 7. Constrained-axis count + hard floor. A 1-axis offset collapses to a
#    single scaled climate distance — exactly the "genomics added nothing"
#    degenerate case B.5 warns against. Floor at 2 (rda.R:361 precedent);
#    if the fit can't support 2 axes, stop rather than silently degrade.
################################################################################
K_max <- RDA_off$CCA$rank
add_diag("rda_axes_max", K_max)
if (K_max < 2) {
    stop("RDA offset: only ", K_max, " constrained axis/axes retained — need >= 2 for a ",
         "meaningful multi-axis offset (see docs/rda_research.md B.5). Curate more candidate ",
         "SNPs or predictors, or reduce condition_pcs (currently ", CONDITION_PCS, ").")
}

set.seed(SEED)
anova_axis <- tryCatch(
    anova.cca(RDA_off, by = "axis", permutations = PERMUTATIONS),
    error = function(e) { message("WARNING: by-axis anova.cca failed: ", e$message); NULL }
)
axis_pvals <- if (!is.null(anova_axis) && "Pr(>F)" %in% colnames(as.data.frame(anova_axis))) {
    as.data.frame(anova_axis)[["Pr(>F)"]][seq_len(K_max)]
} else rep(NA_real_, K_max)

K_floored <- FALSE
if (AXES == "auto") {
    K_sel <- sum(axis_pvals < AXIS_ALPHA, na.rm = TRUE)
    if (K_sel < 2) { K_floored <- TRUE; K_sel <- 2 }
    K_selection <- "auto"
} else {
    K_sel <- as.integer(AXES)
    if (K_sel < 2) stop("RDA offset: axes must be >= 2.")
    K_selection <- "fixed"
}
K_sel <- min(K_sel, K_max)

# B22: cap at 2-3 when site count is small.
n_sites <- length(unique(fread(SAMPLES, colClasses = "character")$site))
if (n_sites < 10 && K_sel > 3) {
    message(paste0("INFO: n_sites=", n_sites, " (small) — capping K_sel at 3 per B22"))
    K_sel <- 3
}
message(paste0("INFO: Selected K=", K_sel, " constrained axes (K_max=", K_max,
               ", selection=", K_selection, ", floored=", K_floored, ")"))
add_diag("rda_axes_used", K_sel)
add_diag("K_selection", K_selection)
add_diag("K_floored", K_floored)
add_diag("n_sites", n_sites)

################################################################################
# 8. Weights — B2/B3: denominator sums over ALL constrained axes, then take
#    the first K_sel. Weight the SCORES before dist(), so they enter squared.
################################################################################
w <- (RDA_off$CCA$eig / sum(RDA_off$CCA$eig))[seq_len(K_sel)]
add_diag("weights", paste(round(w, 4), collapse = ","))

################################################################################
# 9. Project present/future site + raster scores via loadings (B4/G1/G5 —
#    NEVER predict(type="lc"); CCA$biplot are biplot scores, not regression
#    coefficients). Port of rda_offset_sites() from docs/rda_research.md:238-248.
################################################################################
v <- rownames(RDA_off$CCA$biplot)
B <- RDA_off$CCA$biplot[v, seq_len(K_sel), drop = FALSE]

project_scores <- function(env_df, center, scale) {
    # as.matrix() FIRST, then subset by name — env_df can be a data.table (from
    # fread/dplyr::select), whose `[.data.table` treats a bare symbol like `v`
    # as a column-name literal (NSE), not "the columns named by variable v".
    # Matrix indexing has no such NSE and takes a character vector directly.
    z <- scale(as.matrix(env_df)[, v, drop = FALSE], center[v], scale[v])
    z %*% B
}

# Present-side scores and the raster/sample scaffolding are scenario-invariant.
sp_site <- project_scores(env_site_pres_raw, center_env, scale_env)
clim_pres <- rast(PRES_RASTER)
samples <- fread(SAMPLES,
                 colClasses = c('site' = 'character', 'sample' = 'character',
                                'latitude' = 'numeric', 'longitude' = 'numeric'))
site_ids <- samples[, c('site', 'sample')] %>% unique()

################################################################################
# 10-11. Per scenario: project the future, write site TSV + raster + map TSV
################################################################################
for (i in seq_len(N_SCEN)) {
    message(paste0('INFO: [', i, '/', N_SCEN, '] scenario ', SCENARIOS[i]))

    env_site_fut_raw <- fread(ENV_SITE_FUT[i]) %>% dplyr::select(!!PREDICTORS)
    env_all_fut      <- fread(ENV_ALL_FUT[i])  %>% dplyr::select(ID, !!PREDICTORS)

    if (nrow(env_site_fut_raw) != n_samples) {
        stop(paste0('FATAL: Sample-count mismatch -- env_site_pres has ', n_samples,
                    ' rows but ', ENV_SITE_FUT[i], ' has ', nrow(env_site_fut_raw), '.'))
    }

    # Sites (present/future) — B11 constants applied to BOTH.
    sf_site <- project_scores(env_site_fut_raw, center_env, scale_env)
    offset_site <- sqrt(rowSums(sweep((sp_site - sf_site), 2, w, `*`)^2))

    # Raster cells — B14: identical NA mask across present + future.
    ok <- complete.cases(env_all_pres[, -1]) & complete.cases(env_all_fut[, -1])
    env_all_pres_ok <- env_all_pres[ok, ]
    env_all_fut_ok  <- env_all_fut[ok, ]
    message(paste0('INFO: Non-NA raster cells: ', sum(ok)))

    sp_land <- project_scores(env_all_pres_ok[, -1], center_env, scale_env)
    sf_land <- project_scores(env_all_fut_ok[,  -1], center_env, scale_env)
    offset_landscape <- sqrt(rowSums(sweep((sp_land - sf_land), 2, w, `*`)^2))

    # ---- Write site TSV (byte-compatible with GF/GeoOff output) ----
    go_site <- site_ids %>%
        cbind(genetic_offset = offset_site) %>%
        dplyr::arrange(desc(genetic_offset))
    fwrite(go_site, OUT_SITE_TSV[i], sep = '\t')
    message(paste0('INFO: Saved site GO values: ', OUT_SITE_TSV[i]))
    message(paste0('INFO: Site offset range: ', round(min(offset_site), 4), ' - ',
                   round(max(offset_site), 4)))

    # ---- Write raster + map TSV (mirrors geometric_offset.R) ----
    # rast_offset is re-derived from clim_pres each iteration: terra SpatRasters
    # wrap external pointers, so reusing a hoisted template risks mutating it.
    rast_offset <- clim_pres[[PREDICTORS[1]]]
    rast_offset[] <- NA
    rast_offset[env_all_pres_ok$ID] <- offset_landscape
    writeRaster(rast_offset, filename = OUT_RASTER_TIF[i], overwrite = TRUE,
                gdal = c('INTERLEAVE=BAND', 'COMPRESS=LZW'))
    message(paste0('INFO: Saved offset raster: ', OUT_RASTER_TIF[i]))

    GO_matrix <- matrix(values(rast_offset), ncol = ncol(rast_offset), byrow = TRUE)
    GO_matrix %>% as.data.table() %>% fwrite(OUT_MAP_TSV[i], sep = '\t', col.names = FALSE)
    message(paste0('INFO: Saved map GO values: ', OUT_MAP_TSV[i]))
}

################################################################################
# 12. Importance plot — per-predictor contribution: sum of |biplot loading|
#     weighted by axis weight (parallels geometric_offset.R's eigenvalue-
#     weighted-loading contribution, section 11).
################################################################################
tryCatch({
    contrib <- as.numeric(abs(B) %*% w)
    names(contrib) <- v
    png(OUT_IMPORTANCE, width = 800, height = 500)
    par(mar = c(8, 4, 3, 1))
    barplot(sort(contrib, decreasing = TRUE), las = 2, col = CLINEGO_RETAINED,
            main = 'RDA Offset: Predictor Contribution',
            ylab = 'Weighted |biplot loading|', cex.names = 0.85)
    dev.off()
    message(paste0('INFO: Saved importance plot: ', OUT_IMPORTANCE))
}, error = function(e) {
    message(paste0('WARNING: Importance plot failed (', conditionMessage(e), '). Writing placeholder.'))
    png(OUT_IMPORTANCE, width = 400, height = 200)
    par(mar = c(1, 1, 2, 1))
    plot.new()
    text(0.5, 0.5, 'Importance plot unavailable\n(see log for details)', cex = 1.2, col = 'grey40')
    dev.off()
})

################################################################################
# 13. Screeplot (mirrors rda.R:588-605)
################################################################################
eig_dt <- data.table(axis = seq_along(RDA_off$CCA$eig), eigenvalue = RDA_off$CCA$eig)
eig_dt[, retained := axis <= K_sel]
axis_p_dt <- data.table(axis = seq_len(K_max), p = axis_pvals[seq_len(K_max)])
eig_dt <- merge(eig_dt, axis_p_dt, by = "axis", all.x = TRUE)
g_scree <- ggplot(eig_dt, aes(x = factor(axis), y = eigenvalue, fill = retained)) +
    geom_col() +
    geom_text(aes(label = ifelse(!is.na(p), sprintf("p=%.3g", p), "")),
              vjust = -0.3, size = 3, color = CLINEGO_COL$fg) +
    scale_fill_manual(values = c(`TRUE` = CLINEGO_RETAINED, `FALSE` = CLINEGO_REMOVED), guide = "none") +
    labs(x = "Constrained axis", y = "Eigenvalue", title = "RDA offset constrained eigenvalues",
         subtitle = paste0("K=", K_sel, " (", K_selection, "), K_max=", K_max)) +
    theme_clinego()
ggsave(file.path(PLOT_DIR, "rda_axis_screeplot.png"), g_scree, width = 7, height = 5, dpi = 300)
ggsave(file.path(PLOT_DIR, "rda_axis_screeplot.svg"), g_scree, width = 7, height = 5,
       device = svglite::svglite, bg = "transparent")

################################################################################
# 14. Write diagnostics (long key/value)
################################################################################
diag_dt <- data.table(names(diag), sapply(diag, function(x) paste(x, collapse = ",")))
setnames(diag_dt, c("key", "value"))
fwrite(diag_dt, OUT_DIAGNOSTICS, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved diagnostics: ', OUT_DIAGNOSTICS))

message('INFO: RDA genetic offset complete.')
message(paste0('INFO: n_candidates=', ncol(Y_cand), ' K=', K_sel, ' (max=', K_max,
               ') adj_r2=', round(r2adj, 4), ' | offset range ',
               round(min(offset_site), 4), '-', round(max(offset_site), 4)))
