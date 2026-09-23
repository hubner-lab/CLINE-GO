library(dplyr)
library(data.table)
library(vegan)
library(stringr)
library(geosphere)
library(ggplot2)
library(qs)
# Shared plot theme + palette (CLAUDE.md rule 9) — theme_clinego(),
# scale_fill_clinego(), CLINEGO_THRESHOLD.
source("/pipeline/scripts/R/utils/theme_clinego.R")

args = commandArgs(trailingOnly=TRUE)
# Permutation p-values must reproduce run to run. 42 = the pipeline-wide seed
# (hardcoded here like most scripts until it becomes one config parameter).
set.seed(42)
####################################
SAMPLES = args[1]
CLUSTERS = args[2]
ENV = args[3]
PREDICTORS_SELECTED = args[4] %>% str_split(',') %>% unlist
PLOT_DIR = args[5]
INTER_DIR = args[6]
TABLES_DIR = args[7]
####################################

#################################### Functions

analyze_mantel <- function(geo, env, clust) {
  # Calculate distance matrices
  env_dist <- vegdist(env,
                      method = "euclidean",
                      binary = FALSE,
                      diag = FALSE,
                      upper = FALSE,
                      na.rm = FALSE)

  clust_dist <- vegdist(clust,
                        method = 'hellinger')

  geo_dist <- distm(geo, fun = distVincentyEllipsoid)

  # Simple Mantel tests
  ibd_full <- mantel(clust_dist,
                     geo_dist,
                     method = 'pearson',
                     permutations = 999)

  ibe_full <- mantel(clust_dist,
                     env_dist,
                     method = 'pearson',
                     permutations = 999)

  # Partial Mantel tests
  ibd_partial <- mantel.partial(clust_dist,
                                geo_dist,
                                env_dist,
                                method = 'pearson',
                                permutations = 999)

  ibe_partial <- mantel.partial(clust_dist,
                                env_dist,
                                geo_dist,
                                method = 'pearson',
                                permutations = 999)

  # NO VARIANCE PARTITION HERE. This used to square the four r values into a
  # "Geography Only / Environment Only / Geography x Environment / Unexplained"
  # pie, clamp each at 0, derive "shared" by Venn subtraction and renormalise to
  # 1 (twice: here and again in the plot). None of that is a partition: a Mantel
  # r correlates two vectors of PAIRWISE DISTANCES, so r^2 is not a fraction of
  # genomic variance; a partial Mantel r is a partial correlation (a share of the
  # residual), not a semipartial, so the Venn algebra does not hold; squaring
  # drops the sign, so a NEGATIVE IBD r was drawn as a positive geography share;
  # and r_geo = r_env = 0.8 with both partials 0 rendered "100% shared, 0%
  # unexplained". The genuine partition of the same question (dbMEM + partial RDA
  # on adjusted R2) is mode=climate's climate/tables/varpart/variance_partition.tsv.
  # The four statistics are reported as what they are.
  stats_dt <- data.table::data.table(
    test         = c('IBD (geography)', 'IBE (environment)',
                     'IBD | environment', 'IBE | geography'),
    type         = c('mantel', 'mantel', 'partial_mantel', 'partial_mantel'),
    mantel_r     = c(ibd_full$statistic, ibe_full$statistic,
                     ibd_partial$statistic, ibe_partial$statistic),
    p_value      = c(ibd_full$signif, ibe_full$signif,
                     ibd_partial$signif, ibe_partial$signif),
    permutations = c(ibd_full$permutations, ibe_full$permutations,
                     ibd_partial$permutations, ibe_partial$permutations)
  )

  for (i in seq_len(nrow(stats_dt))) {
    message(sprintf('INFO: %-20s r = %+.3f, p = %.3f',
                    stats_dt$test[i], stats_dt$mantel_r[i], stats_dt$p_value[i]))
  }

  results <- list(
    mantel_tests = list(
      ibd_full = ibd_full,
      ibe_full = ibe_full,
      ibd_partial = ibd_partial,
      ibe_partial = ibe_partial
    ),
    statistics = stats_dt
  )

  return(results)
}

# Bar chart of the four Mantel statistics (replaces the variance pie — see the
# comment in analyze_mantel()). Value labels are data, not commentary (rule 8).
create_mantel_plot <- function(results) {
  # as.data.frame, not the data.table: $<- on a data.table emits the
  # shallow-copy warning, and the test suite's warning baseline is 0.
  st <- as.data.frame(results$statistics)
  st$test <- factor(st$test, levels = rev(st$test))
  st$label <- sprintf('r = %+.3f   p = %.3f', st$mantel_r, st$p_value)
  st$fill_grp <- ifelse(st$type == 'partial_mantel', 'Partial Mantel', 'Mantel')

  lo <- min(0, min(st$mantel_r, na.rm = TRUE)) - 0.38
  hi <- max(0, max(st$mantel_r, na.rm = TRUE)) + 0.38

  ggplot(st, aes(x = test, y = mantel_r, fill = fill_grp)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, colour = CLINEGO_THRESHOLD) +
    geom_text(aes(label = label, hjust = ifelse(mantel_r >= 0, -0.08, 1.08)),
              size = 3.4) +
    coord_flip() +
    scale_fill_clinego() +
    scale_y_continuous(limits = c(lo, hi)) +
    labs(x = NULL, y = 'Mantel correlation (r)', fill = NULL) +
    theme_clinego(base_size = 12) +
    theme(legend.position = 'bottom')
}

########################################## Main

# Load + process
geo_raw <- fread(SAMPLES, colClasses = c("site" = "character", "sample" = "character"))
geo_cols <- colnames(geo_raw)

lat_col <- geo_cols[grepl("^lat", geo_cols, ignore.case = TRUE)]
lon_col <- geo_cols[grepl("^lon", geo_cols, ignore.case = TRUE)]

if (length(lat_col) == 0) stop("No latitude column found (expected column starting with 'lat')")
if (length(lon_col) == 0) stop("No longitude column found (expected column starting with 'lon')")

if (length(lat_col) > 1) { message('WARNING: Multiple latitude columns, using: ', lat_col[1]); lat_col <- lat_col[1] }
if (length(lon_col) > 1) { message('WARNING: Multiple longitude columns, using: ', lon_col[1]); lon_col <- lon_col[1] }

#=============================================================================
# SINGLE-SITE GUARD
#=============================================================================
# Every Mantel test here correlates a genetic distance matrix against a
# geographic or environmental one. With one sampling site the geographic matrix
# is all zeros and every per-sample climate column is constant, so both the
# geographic and the environmental term are undefined -- the statistics come back
# NaN and propagate into the plot and the statistics table.
#
# This has to run BEFORE the zero-variance predictor check further down: with one
# site that check drops every predictor and stop()s first, reporting a predictor
# problem for what is really a sampling-design one. Write the declared plot and
# exit 0 -- population stats are opt-in and must not abort structure mode.
n_sites <- length(unique(geo_raw$site))
if (n_sites < 2) {
  message(sprintf(paste0('WARNING: Mantel test skipped -- %d sampling site(s), need >= 2. ',
                         'Geographic and environmental distances are both undefined ',
                         'within a single site.'), n_sites))
  skip_plot <- ggplot() +
    annotate('text', x = 0, y = 0, size = 5, colour = 'grey50',
             label = sprintf('Mantel test unavailable\n(%d sampling site, need >= 2)', n_sites)) +
    theme_void()
  ggsave(paste0(PLOT_DIR, 'mantel_test.png'), skip_plot, width = 10, height = 8, dpi = 150, bg = 'white')
  ggsave(paste0(PLOT_DIR, 'mantel_test.svg'), skip_plot, width = 10, height = 8,
         device = svglite::svglite, bg = 'white', fix_text_size = FALSE)
  qsave(list(status = 'skipped_single_site', n_sites = n_sites),
        paste0(INTER_DIR, 'mantel_test.qs'))
  # mantel_statistics.tsv is a declared rule output — write the empty shape so the
  # opt-in skip stays a skip instead of failing the rule on a missing file.
  dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(test = character(0), type = character(0), mantel_r = numeric(0),
                    p_value = numeric(0), permutations = numeric(0), n_sites = integer(0)),
         paste0(TABLES_DIR, 'mantel_statistics.tsv'), sep = '\t')
  quit(status = 0)
}

geo <- geo_raw %>%
  dplyr::select(longitude = all_of(lon_col), latitude = all_of(lat_col))

if (any(is.na(geo))) {
  message('WARNING: Missing values in geographic coordinates, removing affected rows')
  geo <- na.omit(geo)
}

message(sprintf('INFO: Geographic data - Longitude: [%.2f, %.2f], Latitude: [%.2f, %.2f]',
                min(geo$longitude), max(geo$longitude), min(geo$latitude), max(geo$latitude)))

env_raw <- fread(ENV) %>%
  dplyr::select(all_of(PREDICTORS_SELECTED))

# Align the ancestry matrix to the metadata BY SAMPLE ID, not by position.
# clusters_K{k}.tsv covers every sample that reached sNMF, while SAMPLES here is
# metadata_climate_valid.tsv — the climate-valid subset — so the two tables differ
# in length whenever any sample lacks coordinates or climate values (46 vs 47 on
# the shipped test dataset). The previous positional read subset `clust` with a
# logical vector one element shorter than its own row count, which R recycles, so
# every ancestry row after the dropped sample was paired with the wrong site's
# geography and climate — silently, since nothing compared the lengths.
clust_raw <- fread(CLUSTERS, colClasses = c("sample" = "character", "site" = "character"))
clust_idx <- match(geo_raw$sample, clust_raw$sample)
if (anyNA(clust_idx)) {
  stop('ERROR: ', sum(is.na(clust_idx)), ' sample(s) in ', basename(SAMPLES),
       ' have no row in ', basename(CLUSTERS),
       ' — cannot align the ancestry matrix to the climate/geography tables.')
}
if (nrow(clust_raw) != length(clust_idx)) {
  message(sprintf('INFO: ancestry table has %d samples, %d are climate-valid — subsetting by sample ID',
                  nrow(clust_raw), length(clust_idx)))
}
clust <- clust_raw[clust_idx, ] %>%
  dplyr::select(-sample, -site)

# Remove rows with NA climate values (affects geo, env, clust equally)
complete_rows <- complete.cases(env_raw)
if (any(!complete_rows)) {
  n_removed <- sum(!complete_rows)
  message(sprintf('WARNING: Removing %d samples with missing climate values (%d remain)',
                  n_removed, sum(complete_rows)))
  env_raw <- env_raw[complete_rows, ]
  geo <- geo[complete_rows, ]
  clust <- clust[complete_rows, ]
}

# SITE-LEVEL AGGREGATION. Every one of the three matrices below is per SAMPLE,
# but IBD/IBE is a question about SITES: geographic coordinates and climate values
# are identical for every sample sequenced at a site, so a 30-sample site
# contributes 435 within-site pairs at distance ~0 to both the geographic and the
# environmental matrix, and Mantel's permutation test counts those as independent
# units — the p-value is inflated by sampling effort, not by geography. The same
# imbalance also weights scale()'s centre and sd by sample size, and unlike the
# other consumers of the climate table this one turns the scaled columns into a
# EUCLIDEAN DISTANCE, where column scale is meaningful and the weighting is
# therefore NOT absorbed (cf. the note in download_climate_present.R).
# So: collapse to one row per site first, then scale. Ancestry is averaged per
# site (Q rows sum to 1, so the mean is still a valid ancestry composition —
# the standard population-level unit for a Hellinger distance); coordinates are
# averaged too, in case a site's samples carry slightly jittered coordinates.
site_vec <- geo_raw$site[complete_rows]
stopifnot(length(site_vec) == nrow(env_raw), nrow(geo) == nrow(env_raw),
          nrow(clust) == nrow(env_raw))

site_levels <- unique(site_vec)
n_sites     <- length(site_levels)
agg_mean <- function(df) {
    as.data.frame(do.call(rbind, lapply(site_levels, function(s)
        colMeans(as.matrix(df[site_vec == s, , drop = FALSE]), na.rm = TRUE))))
}
message(sprintf('INFO: collapsing %d samples to %d sites for the Mantel tests',
                length(site_vec), n_sites))
env_raw <- agg_mean(env_raw)
geo     <- agg_mean(geo)
clust   <- agg_mean(clust)

if (n_sites < 4) {
    stop('ERROR: Mantel test needs at least 4 sites (', n_sites, ' present) — ',
         'a distance matrix over 3 sites has 3 pairs and no permutation test is ',
         'meaningful. Disable Population.calc_stats or add sites.')
}
if (n_sites < 8) {
    message(sprintf(paste0('WARNING: only %d sites — %d pairwise distances, and the ',
                           'permutation test cannot resolve p below 1/%d. Read the ',
                           'Mantel statistics as descriptive, not as a test.'),
                    n_sites, n_sites * (n_sites - 1) / 2, factorial(n_sites)))
}

env <- scale(env_raw)

# Drop zero-variance predictors (scale() produces NaN when sd=0)
zero_var <- apply(env, 2, function(x) all(is.nan(x)))
if (any(zero_var)) {
  dropped <- colnames(env)[zero_var]
  message(sprintf('WARNING: Dropping %d zero-variance predictor(s): %s',
                  length(dropped), paste(dropped, collapse = ', ')))
  env <- env[, !zero_var, drop = FALSE]
}

if (ncol(env) == 0) {
  stop('ERROR: All predictors have zero variance — cannot compute environmental distances. ',
       'Check that PREDICTORS_SELECTED contains variables with variation across sites.')
}

# Run Mantel analysis
results <- analyze_mantel(geo, env, clust)

# Plot the four Mantel statistics
gMantel <- create_mantel_plot(results)

# Save. The statistics table is a DECLARED OUTPUT: CLAUDE.md rule 1 forbids
# reading plots, and before this the r values and permutation p-values existed
# only inside the PNG and the .qs.
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(cbind(results$statistics, n_sites = n_sites),
       paste0(TABLES_DIR, 'mantel_statistics.tsv'), sep = '\t')
ggsave(paste0(PLOT_DIR, 'mantel_test.png'), gMantel, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(paste0(PLOT_DIR, 'mantel_test.svg'), gMantel,
       device = svglite::svglite, bg = 'white', fix_text_size = FALSE)
qsave(results, paste0(INTER_DIR, 'mantel_test.qs'))

message('INFO: Mantel test complete')
