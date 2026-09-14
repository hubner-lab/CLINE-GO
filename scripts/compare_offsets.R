#!/usr/bin/env Rscript
# compare_offsets.R — Extrapolation-focused comparison of two (or more) genetic offset models.
#
# Called on-demand by the Shiny app when the user opens a comparison tab.
# Writes a stats.json + PNG/TSV cache to:
#   {PROJECT}_results/_intermediate/model_compare/{comparison_key}/
# Novelty (ExDet NT1/NT2) is model-independent → cached under:
#   {PROJECT}_results/_intermediate/novelty/{predictors_hash}_{scenario}/
#
# Args (positional):
#   1  CACHE_DIR      — absolute path to the per-comparison cache directory
#   2  NOVELTY_DIR    — absolute path to the per-scenario novelty cache directory
#   3  MODEL_A_SITE   — genetic_offset_site.tsv for Model A
#   4  MODEL_A_MAP    — genetic_offset_map.tsv   for Model A (matrix, no header)
#   5  MODEL_A_LABEL  — display label (e.g. "GF / bonf05_union_spatial")
#   6  MODEL_B_SITE   — genetic_offset_site.tsv for Model B
#   7  MODEL_B_MAP    — genetic_offset_map.tsv   for Model B
#   8  MODEL_B_LABEL  — display label
#   9  ENV_SITE_PRES  — present climate at sampling sites
#  10  ENV_ALL_PRES   — present climate for all raster cells (reference)
#  11  ENV_ALL_FUT    — future  climate for all raster cells
#  12  PREDICTORS     — comma-separated predictor names
#  13  RASTER_TIF     — present raster for spatial template
#  14  TOP_K          — integer: how many top-risk sites to use for Jaccard (default 10)
#  15  N_MODELS_EXTRA — count of extra model paths following (space-joined, for N-way)
#  16  METADATA_VALID — metadata_climate_valid.tsv (site, sample, latitude, longitude, ...):
#                       maps the sample-keyed ENV_SITE_PRES rows to sites (one reference row
#                       per site for ExDet) and supplies site coordinates for the Dutilleul test
#  17+ (optional)     — extra model site TSVs (N-way only)

library(data.table)
library(dplyr)
library(terra)
library(jsonlite)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

CACHE_DIR     <- args[1]
NOVELTY_DIR   <- args[2]
MODEL_A_SITE  <- args[3]
MODEL_A_MAP   <- args[4]
MODEL_A_LABEL <- args[5]
MODEL_B_SITE  <- args[6]
MODEL_B_MAP   <- args[7]
MODEL_B_LABEL <- args[8]
ENV_SITE_PRES <- args[9]
ENV_ALL_PRES  <- args[10]
ENV_ALL_FUT   <- args[11]
PREDICTORS    <- str_split(args[12], ',')[[1]]
RASTER_TIF    <- args[13]
TOP_K         <- as.integer(args[14])
N_EXTRA       <- as.integer(args[15])
METADATA_VALID <- args[16]
extra_sites   <- if (N_EXTRA > 0 && length(args) > 16) args[17:(16 + N_EXTRA)] else character(0)
if (is.na(METADATA_VALID) || !file.exists(METADATA_VALID))
    stop(paste0('arg 16 (metadata_climate_valid.tsv) missing or not found: ', METADATA_VALID))

if (is.na(TOP_K) || TOP_K < 1) TOP_K <- 10L

dir.create(CACHE_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(NOVELTY_DIR, recursive = TRUE, showWarnings = FALSE)

message('INFO: compare_offsets.R')
message(paste0('INFO: Model A: ', MODEL_A_LABEL))
message(paste0('INFO: Model B: ', MODEL_B_LABEL))

# ═════════════════════════════════════════════════════════════════════════════
# 1. Load per-site offset tables
# ═════════════════════════════════════════════════════════════════════════════
load_site <- function(path, label) {
    dt <- fread(path, colClasses = c(site = 'character', sample = 'character'))
    # normalise column name (GF uses 'genetic_offset'; future proofing)
    off_col <- setdiff(names(dt), c('site', 'sample', 'latitude', 'longitude'))
    if (length(off_col) == 0) stop(paste0('No offset column in ', path))
    setnames(dt, off_col[1], 'offset')
    # Aggregate to one row per site (mean over samples) — prevents Cartesian join
    # when merging two model tables that both have multiple samples per site.
    dt[, .(offset = mean(offset, na.rm = TRUE), model = label), by = site]
}

site_a <- load_site(MODEL_A_SITE, MODEL_A_LABEL)
site_b <- load_site(MODEL_B_SITE, MODEL_B_LABEL)

# Sample -> site map and per-site coordinates. Every site-level quantity below (the ExDet
# reference, the Dutilleul coordinates) is built from this, never from the sample-keyed
# climate table: download_climate_present.R forbids distance/covariance consumers on the
# per-sample table because N samples per site replicate each site N times.
meta_valid <- fread(METADATA_VALID, colClasses = c(site = 'character', sample = 'character'))
if (!all(c('site', 'sample') %in% names(meta_valid)))
    stop('metadata_climate_valid.tsv must have site and sample columns')
site_coords <- if (all(c('latitude', 'longitude') %in% names(meta_valid))) {
    meta_valid[, .(latitude = mean(latitude, na.rm = TRUE),
                   longitude = mean(longitude, na.rm = TRUE)), by = site]
} else NULL

# Join on site
combined_site <- merge(
    site_a[, .(site, offset_a = offset)],
    site_b[, .(site, offset_b = offset)],
    by = 'site'
)
message(paste0('INFO: Matched sites: ', nrow(combined_site)))

# ═════════════════════════════════════════════════════════════════════════════
# 2. Compute ExDet NT1 + NT2 novelty (inline — ~30 lines, no dsmextra dep)
#    Reference set = present climate at SAMPLED SITES (training space)
# ═════════════════════════════════════════════════════════════════════════════
novelty_raster_path <- file.path(NOVELTY_DIR, 'exdet_novelty.tif')
novelty_table_path  <- file.path(NOVELTY_DIR, 'exdet_novelty.tsv')

# Cache hit requires the post-2026-09-13 schema (pct_cells_type1/2): a raster cached by the
# pre-fix code carries the inverted novelty class and must be recomputed, not reused.
.novelty_cache_current <- function() {
    if (!file.exists(novelty_raster_path) || !file.exists(novelty_table_path)) return(FALSE)
    all(c('pct_cells_type1', 'pct_cells_type2') %in% names(fread(novelty_table_path, nrows = 0)))
}

if (!.novelty_cache_current()) {
    message('INFO: Computing ExDet novelty (NT1 + NT2)')

    # Reference = present climate at the sampled SITES, one row per site. ENV_SITE_PRES is
    # sample-keyed (46 rows on SIMDATA for 9 sites), so its rows are mapped to sites through
    # the metadata and averaged (climate is constant within a site, so the mean is the value).
    env_site_raw <- fread(ENV_SITE_PRES, colClasses = c(sample = 'character'))
    if (!'sample' %in% names(env_site_raw)) stop('ENV_SITE_PRES has no sample column')
    env_site_raw[, site := meta_valid$site[match(sample, meta_valid$sample)]]
    if (anyNA(env_site_raw$site))
        stop(paste0(sum(is.na(env_site_raw$site)), ' sample(s) in ', basename(ENV_SITE_PRES),
                    ' have no row in ', basename(METADATA_VALID)))
    env_site <- env_site_raw[, lapply(.SD, mean, na.rm = TRUE), by = site, .SDcols = PREDICTORS]
    env_all  <- fread(ENV_ALL_PRES)[, c('ID', PREDICTORS), with = FALSE]
    env_fut  <- fread(ENV_ALL_FUT)[,  c('ID', PREDICTORS), with = FALSE]
    ok       <- complete.cases(env_all[, -1]) & complete.cases(env_fut[, -1])
    env_fut_ok <- env_fut[ok]

    ref_mat  <- as.matrix(env_site[, ..PREDICTORS])   # one row per site
    n_ref    <- nrow(ref_mat)
    if (n_ref < length(PREDICTORS) + 2)
        message(paste0('WARNING: only ', n_ref, ' reference sites for ', length(PREDICTORS),
                       ' predictors — the NT2 covariance is rank-deficient or near it'))
    ref_min  <- apply(ref_mat, 2, min, na.rm = TRUE)
    ref_max  <- apply(ref_mat, 2, max, na.rm = TRUE)
    ref_cov  <- cov(ref_mat)
    ref_mean <- colMeans(ref_mat)

    fut_mat  <- as.matrix(env_fut_ok[, -1])  # exclude ID column

    # ExDet (Mesgaran, Cousens & Webber 2014, Divers Distrib 20:1147).
    # NT1 = sum over predictors of the univariate departure UD_j, each <= 0: 0 inside the
    # training range, -(distance beyond the range)/(range) outside. NT1 < 0 <=> type-1
    # (univariate) novelty. Sum, not min: two predictors each 10 % out are more novel than one.
    ref_range <- ref_max - ref_min
    ref_range[ref_range == 0] <- NA_real_   # invariant predictor: no univariate departure defined
    nt1_vals <- apply(fut_mat, 1, function(x) {
        below <- pmin((x - ref_min) / ref_range, 0)  # negative if below min
        above <- pmin((ref_max - x) / ref_range, 0)  # negative if above max
        sum(below + above, na.rm = TRUE)
    })

    # NT2 = squared Mahalanobis distance to the reference centroid, divided by the LARGEST
    # squared Mahalanobis distance among the reference sites themselves, so that NT2 > 1
    # means "outside the reference hull in covariance space" (type-2, combinatorial
    # novelty) on the same scale for every dataset. Before 2026-09-13 the raw D^2 was
    # written and no type-2 threshold existed (audit SC2).
    inv_cov  <- tryCatch(solve(ref_cov), error = function(e) MASS::ginv(ref_cov))
    d2_of    <- function(m) {
        delta <- sweep(m, 2, ref_mean, '-')
        apply(delta, 1, function(d) as.numeric(d %*% inv_cov %*% d))
    }
    ref_d2_max <- max(d2_of(ref_mat))
    nt2_vals   <- d2_of(fut_mat) / ref_d2_max

    # Combined ExDet band: NT1 where the cell is outside the univariate range (< 0), else the
    # NT2 ratio. Analog <=> 0 <= ExDet <= 1; novel <=> ExDet < 0 (type 1) or ExDet > 1 (type 2).
    exdet    <- ifelse(nt1_vals < 0, nt1_vals, nt2_vals)
    novel_t1 <- nt1_vals < 0
    novel_t2 <- !novel_t1 & nt2_vals > 1

    # Write to raster — use [[1]] (band index) so this is robust regardless of
    # which band names the template raster carries. The template is used only for
    # its spatial geometry (CRS, extent, resolution), never for its values.
    clim_pres  <- rast(RASTER_TIF)
    nt1_rast   <- clim_pres[[1]]; nt1_rast[] <- NA
    nt2_rast   <- clim_pres[[1]]; nt2_rast[] <- NA
    exdet_rast <- clim_pres[[1]]; exdet_rast[] <- NA

    nt1_rast[env_fut_ok$ID]   <- nt1_vals
    nt2_rast[env_fut_ok$ID]   <- nt2_vals
    exdet_rast[env_fut_ok$ID] <- exdet

    novelty_stack <- c(nt1_rast, nt2_rast, exdet_rast)
    names(novelty_stack) <- c('NT1', 'NT2', 'ExDet')

    writeRaster(novelty_stack, novelty_raster_path, overwrite = TRUE,
                gdal = c('INTERLEAVE=BAND', 'COMPRESS=LZW'))

    # Before 2026-09-13 this counted exdet >= 0 as novel, i.e. every IN-range cell (NT2 is
    # never negative) — the complement of the ExDet definition (audit SC2).
    pct_novel <- round(100 * mean(novel_t1 | novel_t2, na.rm = TRUE), 1)
    pct_t1    <- round(100 * mean(novel_t1, na.rm = TRUE), 1)
    pct_t2    <- round(100 * mean(novel_t2, na.rm = TRUE), 1)
    max_nt2   <- round(max(nt2_vals, na.rm = TRUE), 2)
    fwrite(data.table(
        pct_cells_novel  = pct_novel,
        pct_cells_type1  = pct_t1,
        pct_cells_type2  = pct_t2,
        max_nt2          = max_nt2,
        n_cells_valid    = sum(ok),
        n_reference_sites = n_ref
    ), novelty_table_path, sep = '\t')

    message(paste0('INFO: Novelty computed. Novel cells: ', pct_novel, '% (type 1 ', pct_t1,
                   '%, type 2 ', pct_t2, '%), max NT2 ratio: ', max_nt2,
                   ', reference sites: ', n_ref))
} else {
    message('INFO: Novelty cache hit (schema current)')
}

# ═════════════════════════════════════════════════════════════════════════════
# 3. Per-site rank concordance (Tab 3)
#    Spearman + Kendall + Dutilleul modified.ttest (autocorr-corrected p)
# ═════════════════════════════════════════════════════════════════════════════
message('INFO: Computing per-site rank concordance')

rk_a  <- rank(combined_site$offset_a)
rk_b  <- rank(combined_site$offset_b)
spear <- cor.test(rk_a, rk_b, method = 'spearman', exact = FALSE)
kend  <- cor.test(combined_site$offset_a, combined_site$offset_b, method = 'kendall', exact = FALSE)

# Site coordinates for the spatial-autocorrelation-corrected correlation. Until 2026-09-13
# this merged `site` against the climate table's `sample` column and looked for lat/lon
# columns it does not have, so the test never ran (audit B21).
dutilleul_p <- NA_real_
dutilleul_n_eff <- NA_real_
if (!is.null(site_coords)) {
    site_env <- merge(combined_site, site_coords, by = 'site', all.x = FALSE)
    tryCatch({
        if (!requireNamespace('SpatialPack', quietly = TRUE)) stop('SpatialPack not installed')
        if (nrow(site_env) < 4) stop('fewer than 4 sites with coordinates')
        coords <- as.matrix(site_env[, .(longitude, latitude)])
        mt <- SpatialPack::modified.ttest(site_env$offset_a,
                                          site_env$offset_b,
                                          coords, nclass = 13)
        dutilleul_p   <- mt$p.value
        dutilleul_n_eff <- mt$dof
    }, error = function(e) {
        message(paste0('WARNING: Dutilleul modified.ttest failed: ', conditionMessage(e)))
    })
}

# ═════════════════════════════════════════════════════════════════════════════
# 4. Top-K high-risk site stability — Jaccard (Tab 4)
# ═════════════════════════════════════════════════════════════════════════════
k_use <- min(TOP_K, nrow(combined_site))
top_a <- combined_site[order(-offset_a)]$site[seq_len(k_use)]
top_b <- combined_site[order(-offset_b)]$site[seq_len(k_use)]
# With k >= n_sites both top-k sets are "every site" and the Jaccard is 1 by construction
# (audit B21): report NA rather than a meaningless perfect agreement.
jaccard <- if (k_use >= nrow(combined_site)) NA_real_ else
    length(intersect(top_a, top_b)) / length(union(top_a, top_b))

rank_stability_dt <- combined_site[, .(
    site,
    rank_a = rank(-offset_a),
    rank_b = rank(-offset_b)
)]
rank_stability_dt[, rank_diff := abs(rank_a - rank_b)]
fwrite(rank_stability_dt, file.path(CACHE_DIR, 'rank_stability.tsv'), sep = '\t')

# ═════════════════════════════════════════════════════════════════════════════
# 5. Per-cell disagreement × novelty (Tab 2 — HEADLINE)
#    Rank-transform each model's raster; |Δrank| disagrement map
#    Binned by ExDet novelty class + Spearman(|Δrank|, NT2) with subsampling
# ═════════════════════════════════════════════════════════════════════════════
message('INFO: Computing per-cell disagreement map')

map_a_raw <- as.matrix(fread(MODEL_A_MAP, header = FALSE))
map_b_raw <- as.matrix(fread(MODEL_B_MAP, header = FALSE))

vals_a <- as.vector(t(map_a_raw))
vals_b <- as.vector(t(map_b_raw))

# Rank-transform (ties.method = average, NA stays NA)
valid  <- !is.na(vals_a) & !is.na(vals_b)
rk_a_c <- vals_a; rk_b_c <- vals_b
rk_a_c[valid] <- rank(vals_a[valid], ties.method = 'average')
rk_b_c[valid] <- rank(vals_b[valid], ties.method = 'average')
delta_rank <- abs(rk_a_c - rk_b_c)

# Write disagreement raster — [[1]] for geometry only (band-name-agnostic)
clim_pres     <- rast(RASTER_TIF)
disagree_rast <- clim_pres[[1]]
disagree_rast[] <- NA
disagree_rast[seq_along(delta_rank)] <- delta_rank
writeRaster(disagree_rast, file.path(CACHE_DIR, 'disagree_rank.tif'),
            overwrite = TRUE, gdal = c('INTERLEAVE=BAND', 'COMPRESS=LZW'))

# Spearman(|Δrank|, NT2) with spatial subsampling (cap at 50k for speed)
novelty_stack <- rast(novelty_raster_path)
nt2_vals      <- values(novelty_stack[['NT2']])
ok_cells      <- !is.na(delta_rank) & !is.na(nt2_vals)
n_ok          <- sum(ok_cells)
SUB_CAP       <- 50000L
if (n_ok > SUB_CAP) {
    set.seed(42)
    idx <- sample(which(ok_cells), SUB_CAP)
} else {
    idx <- which(ok_cells)
}
disagree_nt2_cor <- tryCatch(
    cor(delta_rank[idx], nt2_vals[idx], method = 'spearman', use = 'complete.obs'),
    error = function(e) NA_real_
)

# Bin by novelty class: analog <=> 0 <= ExDet <= 1; non-analog <=> ExDet < 0 (outside the
# univariate range) or ExDet > 1 (outside the reference Mahalanobis hull)
exdet_vals <- values(novelty_stack[['ExDet']])
novelty_class <- ifelse(exdet_vals < 0 | exdet_vals > 1, 'non-analog', 'analog')
binned_dt <- data.table(
    delta_rank    = delta_rank[ok_cells],
    novelty_class = novelty_class[ok_cells],
    nt2           = nt2_vals[ok_cells]
)
binned_summary <- binned_dt[, .(
    mean_delta  = mean(delta_rank),
    median_delta = median(delta_rank),
    n           = .N
), by = novelty_class]
fwrite(binned_summary, file.path(CACHE_DIR, 'disagree_by_novelty.tsv'), sep = '\t')
fwrite(combined_site,  file.path(CACHE_DIR, 'site_offsets.tsv'), sep = '\t')

# ═════════════════════════════════════════════════════════════════════════════
# 6. N-way concordance: Kendall's W + ensemble offset (Tab 5)
# ═════════════════════════════════════════════════════════════════════════════
kendall_w    <- NA_real_
kendall_w_p  <- NA_real_
nway_n_models <- 2L + length(extra_sites)

if (nway_n_models >= 2) {
    all_site_offsets <- list(
        setNames(combined_site$offset_a, combined_site$site),
        setNames(combined_site$offset_b, combined_site$site)
    )
    for (ep in extra_sites) {
        tryCatch({
            ext <- load_site(ep, 'extra')
            merged_ext <- merge(combined_site[, .(site)], ext, by = 'site')
            all_site_offsets <- c(all_site_offsets,
                list(setNames(merged_ext$offset, merged_ext$site)))
        }, error = function(e) {
            message(paste0('WARNING: extra model site table ', ep, ' skipped: ', conditionMessage(e)))
        })
    }
    nway_n_models <- length(all_site_offsets)   # models actually loaded, not requested
    # Build rank matrix (sites × models)
    sites_common <- Reduce(intersect, lapply(all_site_offsets, names))
    if (length(sites_common) >= 2 && length(all_site_offsets) >= 2) {
        rank_mat <- sapply(all_site_offsets, function(v) rank(-v[sites_common]))
        # Kendall's W: vegan::kendall.global(Y) takes OBJECTS in rows and JUDGES in columns —
        # here sites x models, i.e. rank_mat as built. It returns a single element
        # `Concordance_analysis` (matrix; rows W, F, Prob.F, Chi2, Prob.perm). Until
        # 2026-09-13 the matrix was passed transposed (9 sites judging 2 models: W = 0.00
        # instead of 0.82 on SIMDATA) and two non-existent fields were read, whose NULL then
        # crashed `if (!is.na(NULL))` below before stats.json was ever written (audit B2/ST5).
        if (requireNamespace('vegan', quietly = TRUE)) {
            tryCatch({
                set.seed(42)
                kw <- vegan::kendall.global(rank_mat, nperm = 999)
                ca <- kw$Concordance_analysis
                kendall_w   <- as.numeric(ca['W', 1])
                kendall_w_p <- as.numeric(ca['Prob.perm', 1])
            }, error = function(e) {
                message(paste0('WARNING: Kendall W via vegan failed: ', conditionMessage(e)))
            })
        }
        # Ensemble: mean rank across models → normalised consensus offset
        mean_rank <- rowMeans(rank_mat)
        nway_dt   <- data.table(site = sites_common, mean_rank = mean_rank)
        fwrite(nway_dt, file.path(CACHE_DIR, 'nway_mean_rank.tsv'), sep = '\t')
    }
}

# ═════════════════════════════════════════════════════════════════════════════
# 7. Persist stats.json
# ═════════════════════════════════════════════════════════════════════════════
# NULL-safe: `if (!is.na(NULL))` is `if (logical(0))` and errors, so every optional statistic
# goes through this instead of a bare is.na() test.
num_or_na <- function(x, f, digits) {
    if (is.null(x) || length(x) == 0 || is.na(x)) NA else f(x, digits)
}
stats <- list(
    model_a_label       = MODEL_A_LABEL,
    model_b_label       = MODEL_B_LABEL,
    n_sites             = nrow(combined_site),
    spearman_rho        = round(spear$estimate, 4),
    spearman_p          = signif(spear$p.value, 4),
    kendall_tau         = round(kend$estimate, 4),
    kendall_p           = signif(kend$p.value, 4),
    dutilleul_p         = num_or_na(dutilleul_p, signif, 4),
    dutilleul_n_eff     = num_or_na(dutilleul_n_eff, round, 1),
    jaccard_top_k       = num_or_na(jaccard, round, 4),
    top_k               = k_use,
    disagree_nt2_spear  = num_or_na(disagree_nt2_cor, round, 4),
    kendall_w           = num_or_na(kendall_w, round, 4),
    kendall_w_p         = num_or_na(kendall_w_p, signif, 4),
    n_models_nway       = nway_n_models,
    novelty_cache       = NOVELTY_DIR
)
writeLines(toJSON(stats, auto_unbox = TRUE, pretty = TRUE),
           file.path(CACHE_DIR, 'stats.json'))
message(paste0('INFO: Done. Stats written to ', file.path(CACHE_DIR, 'stats.json')))
