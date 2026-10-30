library(terra)
library(data.table)
library(dplyr)
library(stringr)

source("/pipeline/scripts/R/lib/climate_cells.R")  # climate_cell_table()

args = commandArgs(trailingOnly=TRUE)
#################################
SAMPLES = args[1]
MODEL_RASTERS = args[2]  # comma-separated list of per-model .tif files
N_MODELS = args[3] %>% as.numeric
PRESENT_RASTER = args[4]  # present climate raster for cell count validation
PRESENT_ALL = args[5]     # present climate all values for cell count validation
OUTPUT_RASTER = args[6]
OUTPUT_ALL = args[7]
OUTPUT_SITE = args[8]
OUTPUT_NA_EXCLUDED = args[9]
#################################

# Approximate distance (km) from an NA point to the nearest non-NA cell -- see
# download_climate_present.R for the identical helper (kept duplicated, not shared, matching
# this script's existing standalone-Rscript convention).
FUN_nearest_valid_distance_km <- function(raster_layer, lon, lat, window = 0.1, max_window = 3.2) {
  pt <- vect(data.frame(long = lon, lat = lat), geom = c("long", "lat"), crs = "EPSG:4326")
  w <- window
  while (w <= max_window) {
    win_rast <- tryCatch(
      terra::crop(raster_layer, ext(lon - w, lon + w, lat - w, lat + w)),
      error = function(e) NULL
    )
    if (!is.null(win_rast)) {
      valid_pts <- terra::as.points(win_rast, na.rm = TRUE)
      if (length(valid_pts) > 0) {
        return(round(min(terra::distance(pt, valid_pts)) / 1000, 2))  # meters -> km
      }
    }
    w <- w * 4
  }
  return(NA_real_)
}

message('INFO: Merging future climate rasters from individual models')

# Load samples for site extraction
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# Load all per-model rasters
model_files <- MODEL_RASTERS %>% str_split(',') %>% unlist
message(paste0('INFO: Loading ', length(model_files), ' model rasters'))

clim_list <- lapply(model_files, function(f) {
  message(paste0('INFO:   Loading: ', f))
  rast(f)
})

# Get layer names from first model
layer_names <- names(clim_list[[1]])
n_layers <- length(layer_names)
message(paste0('INFO: ', n_layers, ' bioclimatic layers per model'))

# Average across models
if (length(clim_list) == 1) {
  clim_future <- clim_list[[1]]
} else {
  # Stack all models and average per layer
  all_layers <- rast(clim_list)
  clim_future <- tapp(all_layers, rep(1:n_layers, N_MODELS), fun = "mean", na.rm = TRUE)
  names(clim_future) <- layer_names
}

message('INFO: Model averaging complete')

# Every land pixel with its cell ID, for raster remapping
clim_future_all <- climate_cell_table(clim_future)

# Validate cell count matches present climate
message('INFO: Validating cell count matches present climate...')
present_all <- fread(PRESENT_ALL)

future_cells <- nrow(clim_future_all)
present_cells <- nrow(present_all)

message(paste0('INFO: Present climate cells: ', present_cells))
message(paste0('INFO: Future climate cells: ', future_cells))

if (future_cells != present_cells) {
  stop(paste0('ERROR: Cell count mismatch - present: ', present_cells, ', future: ', future_cells))
}
message('INFO: Cell count validation passed')

# Save merged raster
writeRaster(clim_future,
            filename = OUTPUT_RASTER,
            overwrite = TRUE,
            gdal = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
message(paste0('INFO: Saved merged raster: ', OUTPUT_RASTER))

# Save all values
clim_future_all %>%
  fwrite(OUTPUT_ALL, sep = '\t')
message(paste0('INFO: Saved all values: ', OUTPUT_ALL))

# Extract site values
coords <- vect(data.frame(long = samples$longitude, lat = samples$latitude),
               geom = c("long", "lat"), crs = "EPSG:4326")

clim_future_site <- terra::extract(clim_future, coords)[, -1] %>% as.data.frame

# Validate: check for samples with NA future climate values. SAMPLES is already narrowed to
# the climate-valid tier upstream (metadata_climate_valid), so this is normally a no-op defense
# check -- it only fires on a residual future-only NA (e.g. a different NoData mask in a CMIP6
# model raster than in the present-climate raster). Always excluded, never a hard stop -- same
# treatment as download_climate_present.R.
na_rows <- which(rowSums(is.na(clim_future_site)) > 0)

excluded_dt <- data.table(sample = character(), site = character(), latitude = numeric(),
                          longitude = numeric(), reason = character(), distance_km = numeric())

if (length(na_rows) > 0) {
    bad_samples <- samples[na_rows, ]
    reasons <- character(nrow(bad_samples))
    distances <- rep(NA_real_, nrow(bad_samples))
    message("WARNING: Future climate extraction returned NA for the following samples:")
    for (i in seq_len(nrow(bad_samples))) {
        is_placeholder <- bad_samples$latitude[i] == 0 && bad_samples$longitude[i] == 0
        reasons[i] <- if (is_placeholder) "likely placeholder (0,0)" else "falls on ocean/NoData pixel"
        dist_msg <- ""
        if (!is_placeholder) {
            distances[i] <- FUN_nearest_valid_distance_km(clim_future[[1]],
                                                            bad_samples$longitude[i], bad_samples$latitude[i])
            if (!is.na(distances[i])) dist_msg <- paste0(" (~", distances[i], " km to nearest valid pixel)")
        }
        message(paste0("  - ", bad_samples$sample[i], " (", bad_samples$site[i],
                       ") at (", bad_samples$latitude[i], ", ", bad_samples$longitude[i], "): ", reasons[i], dist_msg))
    }

    warning(paste0("Future climate extraction returned NA for ", length(na_rows), " sample(s) -- ",
                   "excluding from climate-VALUE-dependent steps. See the future ",
                   "climate_na_excluded table. Fix coordinates in input metadata if unexpected."))

    excluded_dt <- data.table(sample = bad_samples$sample, site = bad_samples$site,
                              latitude = bad_samples$latitude, longitude = bad_samples$longitude,
                              reason = reasons, distance_km = distances)

    clim_future_site <- clim_future_site[-na_rows, , drop = FALSE]
    samples <- samples[-na_rows, ]
}

fwrite(excluded_dt, OUTPUT_NA_EXCLUDED, sep = '\t', quote = FALSE)

# Save site values (with sample column for traceability)
cbind(sample = samples$sample, clim_future_site) %>%
  fwrite(OUTPUT_SITE, sep = '\t')
message(paste0('INFO: Saved site values: ', OUTPUT_SITE))
message('INFO: Future climate merge complete')
