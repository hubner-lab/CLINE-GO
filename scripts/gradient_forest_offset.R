library(qs)
library(data.table)
library(dplyr)
library(terra)
library(stringr)
library(gradientForest)
args = commandArgs(trailingOnly=TRUE)
###############
# FUTURE_ALL and the three OUTPUT_* args are COMMA-SEPARATED LISTS, one entry per
# future scenario, in the same order. A single-scenario call passes one path in
# each and behaves exactly as before scenarios existed.
#
# Everything that depends only on the model and the PRESENT climate -- the model
# load, the present-climate prediction, the template raster, the sample coords --
# is computed once and reused across scenarios. That setup dominates runtime for
# small SNP panels, which is why a 42-garden sweep projects from one job instead
# of re-loading the same model 42 times.
GF_PATH = args[1]
PREDICTORS_SELECTED = args[2] %>% str_split(',') %>% unlist
FUTURE_ALL = args[3] %>% str_split(',') %>% unlist
PRESENT_ALL = args[4]
PRESENT_RASTER = args[5]
SAMPLES = args[6]
OUTPUT_RASTER = args[7] %>% str_split(',') %>% unlist
OUTPUT_MAP_VALUES = args[8] %>% str_split(',') %>% unlist
OUTPUT_SITE_VALUES = args[9] %>% str_split(',') %>% unlist
SCENARIOS = if (length(args) >= 10 && nzchar(args[10]) && args[10] != 'NULL') {
  args[10] %>% str_split(',') %>% unlist
} else {
  paste0('scenario_', seq_along(FUTURE_ALL))
}
# args[11]: predict.gradientForest(extrap=): TRUE = linear continuation of each turnover
# function beyond the training range (package default), FALSE = plateau at the range limit.
# Config key Maladaptation.methods.gradient_forest.extrap; see common.smk for the rationale.
EXTRAP = if (length(args) >= 11) toupper(args[11]) == 'TRUE' else TRUE
# args[12]: scenario-free diagnostics table (quantity, value): the extrap setting, the
# training envelope per predictor and the share of raster cells outside it.
OUTPUT_DIAGNOSTICS = if (length(args) >= 12) args[12] else 'NULL'
###############

n_scen <- length(FUTURE_ALL)
stopifnot(length(OUTPUT_RASTER) == n_scen,
          length(OUTPUT_MAP_VALUES) == n_scen,
          length(OUTPUT_SITE_VALUES) == n_scen,
          length(SCENARIOS) == n_scen)

message(paste0('INFO: Calculating Genetic Offset for ', n_scen, ' scenario(s)'))

# ---- loop-invariant setup: model, present climate, template raster, coords ----

# Load GF model
gf <- qread(GF_PATH)
if (!inherits(gf, 'gradientForest')) {
  stop(paste0('GF model at ', GF_PATH, ' is not a fitted gradientForest object',
              if (is.list(gf) && !is.null(gf$reason)) paste0(' (', gf$status, ': ', gf$reason, ')') else ''))
}

# Load climate data (select only predictors + ID)
present_climate_all <- fread(PRESENT_ALL) %>% dplyr::select(ID, !!PREDICTORS_SELECTED)
message(paste0('INFO: Present climate cells: ', nrow(present_climate_all)))

message(paste0('INFO: extrapolation beyond the training range: ',
               if (EXTRAP) 'linear (extrap = TRUE, package default)' else 'plateau (extrap = FALSE)'))

# Training envelope per predictor, from the fitted object's own predictor table. Cells
# outside it are where extrap= decides the offset, so their share is reported.
train_X   <- as.data.frame(gf$X)[, PREDICTORS_SELECTED, drop = FALSE]
train_min <- vapply(train_X, min, numeric(1), na.rm = TRUE)
train_max <- vapply(train_X, max, numeric(1), na.rm = TRUE)
outside_range <- function(tab) {
  # per-predictor logical matrix: TRUE where the cell lies outside [train_min, train_max]
  m <- as.matrix(tab[, ..PREDICTORS_SELECTED])
  sweep(m, 2, train_min, '<') | sweep(m, 2, train_max, '>')
}
diag_rows <- list()
put <- function(k, v) diag_rows[[length(diag_rows) + 1]] <<- data.table(quantity = k, value = as.character(v))
put('extrap', EXTRAP)
put('n_training_rows', nrow(train_X))   # rows of the fitted predictor table (samples, or sites for a frequency fit)
for (p in PREDICTORS_SELECTED) {
  put(paste0('train_min__', p), signif(train_min[[p]], 6))
  put(paste0('train_max__', p), signif(train_max[[p]], 6))
}
pres_out <- outside_range(present_climate_all)
put('pct_present_cells_outside_range', round(100 * mean(rowSums(pres_out) > 0), 1))
for (p in PREDICTORS_SELECTED) put(paste0('pct_present_cells_outside_range__', p), round(100 * mean(pres_out[, p]), 1))

# Predict cumulative importance for present -- identical for every scenario
pred <- predict(gf, present_climate_all[, -1], extrap = EXTRAP)

clim_present <- rast(PRESENT_RASTER)

samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))
coords <- vect(data.frame(long = samples$longitude, lat = samples$latitude),
               geom = c("long", "lat"), crs = "EPSG:4326")
site_ids <- samples[, c('site', 'sample')] %>% unique

# ---- batched future prediction ----
#
# predict.gradientForest applies the fitted splines row by row, so rows are
# independent and stacking every scenario into ONE call gives identical values
# to calling it per scenario. It is far cheaper: on a large panel the per-call
# overhead against an ~850 MB forest dominated the whole sweep (one job spent
# 2300 s on 42 separate calls).
#
# Scenario tables are all derived from the same present raster grid, so they
# normally share a row count -- but that is not enforced anywhere upstream, so
# offsets are tracked explicitly and each scenario is sliced by its own length.
future_tabs <- lapply(FUTURE_ALL, function(f) {
  fread(f) %>% dplyr::select(ID, !!PREDICTORS_SELECTED)
})
future_n <- vapply(future_tabs, nrow, integer(1))
message(paste0('INFO: Future climate cells per scenario: ',
               paste(range(future_n), collapse = '-'),
               ' (total ', sum(future_n), ')'))

for (i in seq_along(future_tabs)) {
  fut_out <- outside_range(future_tabs[[i]])
  put(paste0('pct_future_cells_outside_range__', SCENARIOS[i]), round(100 * mean(rowSums(fut_out) > 0), 1))
  for (p in PREDICTORS_SELECTED)
    put(paste0('pct_future_cells_outside_range__', SCENARIOS[i], '__', p), round(100 * mean(fut_out[, p]), 1))
}
if (OUTPUT_DIAGNOSTICS != 'NULL') {
  fwrite(rbindlist(diag_rows), OUTPUT_DIAGNOSTICS, sep = '\t')
  message(paste0('INFO: Saved projection diagnostics: ', OUTPUT_DIAGNOSTICS))
}

pred_future_all <- predict(gf, rbindlist(lapply(future_tabs, function(x) x[, -1])), extrap = EXTRAP)
future_end   <- cumsum(future_n)
future_start <- future_end - future_n + 1L

# ---- per-scenario output ----

for (i in seq_len(n_scen)) {
  message(paste0('INFO: [', i, '/', n_scen, '] scenario ', SCENARIOS[i]))

  future_climate_all <- future_tabs[[i]]
  pred.future <- pred_future_all[future_start[i]:future_end[i], , drop = FALSE]

  # Calculate genetic offset: GO = sqrt(sum((future - present)^2))
  genetic_offset <- sapply(1:ncol(pred.future), function(j) {
    (pred.future[, j] - pred[, j])^2
  }) %>%
    rowSums(na.rm = TRUE) %>%
    sqrt()

  message(paste0('INFO: Genetic offset range: ', round(min(genetic_offset, na.rm = TRUE), 4),
                 ' - ', round(max(genetic_offset, na.rm = TRUE), 4)))

  # Create genetic offset raster. Re-derived from clim_present each iteration
  # rather than copied from a hoisted template: terra SpatRasters wrap external
  # pointers, so `x <- template; x[i] <- v` risks mutating the template too.
  rast_offset <- clim_present[[PREDICTORS_SELECTED[1]]]
  rast_offset[present_climate_all$ID] <- genetic_offset

  # Save offset raster
  writeRaster(rast_offset,
              filename = OUTPUT_RASTER[i],
              overwrite = TRUE,
              gdal = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
  message(paste0('INFO: Saved offset raster: ', OUTPUT_RASTER[i]))

  # Save GO values for whole map (matrix format)
  GO_values <- matrix(values(rast_offset), ncol = ncol(rast_offset), byrow = TRUE)
  GO_values %>%
    as.data.table %>%
    fwrite(OUTPUT_MAP_VALUES[i], sep = '\t', col.names = FALSE)
  message(paste0('INFO: Saved map GO values: ', OUTPUT_MAP_VALUES[i]))

  # Save GO values per sampling site
  go_site <- terra::extract(rast_offset, coords)[, -1] %>%
    as.data.frame %>%
    cbind(site_ids, .) %>%
    setNames(c('site', 'sample', 'genetic_offset')) %>%
    dplyr::arrange(desc(genetic_offset))

  go_site %>%
    fwrite(OUTPUT_SITE_VALUES[i], sep = '\t')
  message(paste0('INFO: Saved site GO values: ', OUTPUT_SITE_VALUES[i]))
}

message('INFO: Genetic offset calculation complete')
