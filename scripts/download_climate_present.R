library(terra)
library(data.table)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
#################################
SAMPLES = args[1]
CLIMATE_EXTENT = args[2]
GAP = args[3] %>% as.numeric
DATA_DIR = args[4]      # directory for storing downloaded climate data
RESOLUTION = args[5] %>% as.numeric
RASTER_DIR = args[6]    # for raster output
TABLES_DIR = args[7]    # for tsv outputs
#################################

#################################### Functions
FUN_define_crop_region <- function(samples, gap) {
  min_lat <- min(samples$latitude) - gap
  max_lat <- max(samples$latitude) + gap
  min_long <- min(samples$longitude) - gap
  max_long <- max(samples$longitude) + gap
  return(c(min_long, max_long, min_lat, max_lat))
}

# Download WorldClim data directly from geodata.ucdavis.edu
FUN_download_worldclim <- function(data_dir, resolution = 0.5) {
  # Convert resolution to WorldClim format
  resolution_format <- switch(as.character(resolution),
                              '0.5' = "30s",
                              '2.5' = "2.5m",
                              '5'  = "5m",
                              '10' = "10m")

  tif_dir <- file.path(data_dir, paste0("wc2.1_", resolution_format))
  zip_file <- file.path(data_dir, paste0("wc2.1_", resolution_format, "_bio.zip"))

  # Check if tif files already exist
  expected_tifs <- file.path(tif_dir, paste0("wc2.1_", resolution_format, "_bio_", 1:19, ".tif"))

  if (all(file.exists(expected_tifs))) {
    message("INFO: All bioclimatic tif files already exist, skipping download")
    return(tif_dir)
  }

  # Create directory if needed
  if (!dir.exists(tif_dir)) {
    dir.create(tif_dir, recursive = TRUE)
  }

  # Download zip if not exists
  if (!file.exists(zip_file)) {
    url <- paste0("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_",
                  resolution_format, "_bio.zip")
    message(paste0("INFO: Downloading WorldClim bioclimatic data from: ", url))
    old_timeout <- getOption("timeout")
    options(timeout = Inf)
    on.exit(options(timeout = old_timeout), add = TRUE)
    download.file(url, zip_file, mode = "wb", quiet = FALSE)
  } else {
    message("INFO: Zip file already exists, skipping download")
  }

  # Unzip using system unzip (R's unzip can fail on large files)
  message("INFO: Extracting tif files...")
  exit_code <- system2("unzip", args = c("-o", "-j", shQuote(zip_file), "-d", shQuote(tif_dir)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(exit_code, "status")) && attr(exit_code, "status") != 0) {
    # Fallback to R unzip
    warning("System unzip failed, trying R unzip...")
    unzip(zip_file, exdir = tif_dir)
  }

  return(tif_dir)
}

# Approximate distance (km) from an NA point to the nearest non-NA cell, by cropping a
# small window around just that point (independent of the overall CLIMATE_EXTENT crop) and
# searching outward. Cheap: only ever operates on a handful of NA points, never the full raster.
# Returns NA if no valid cell is found within a few widenings (degrades gracefully, never errors).
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

FUN_present_climate <- function(samples,
                                data_dir,
                                crop_region,
                                resolution = 0.5) {

  # Download/check WorldClim data
  tif_dir <- FUN_download_worldclim(data_dir, resolution)

  # Construct folder path based on resolution
  resolution_format <- switch(as.character(resolution),
                              '0.5' = "30s",
                              '2.5' = "2.5m",
                              '5'  = "5m",
                              '10' = "10m")

  clim.list <- dir(tif_dir,
                   full.names = TRUE,
                   pattern = paste0("^wc2\\.1_", resolution_format, "_bio_[0-9]+\\.tif$"))
  message(clim.list %>% str)

  if (length(clim.list) == 0) {
    stop(paste0("ERROR: No bioclimatic tif files found in ", tif_dir))
  }

  # Load as SpatRaster
  clim.layer <- rast(clim.list)
  message('INFO: Present climate SpatRaster loaded')
  message(samples %>% str)

  # Create SpatVector for site coordinates
  coords <- vect(data.frame(long = samples$longitude, lat = samples$latitude),
                 geom = c("long", "lat"), crs = "EPSG:4326")

  # Crop
  if (length(crop_region) == 1 && crop_region == 'world') {
    clim.present <- clim.layer
  } else {
    clim.present <- terra::crop(clim.layer, ext(crop_region))
  }

  names(clim.present) <- gsub(paste0("wc2\\.1_", resolution_format, "_"), "", names(clim.present))

  # Extract data for coordinates
  clim.present.coords <- terra::extract(clim.present, coords)[, -1] %>%
    as.data.frame

  # Extract all pixel values with cell IDs for raster remapping
  # terra::extract with cell indices doesn't return ID column (unlike raster package)
  clim.land <- terra::extract(clim.present, 1:ncell(clim.present))
  clim.land$ID <- 1:ncell(clim.present)
  clim.land <- na.omit(clim.land) %>%
    dplyr::select(ID, everything())

  return(list(RasterStack = clim.present,
              SiteValues = clim.present.coords,
              AllValues = clim.land))
}

################################### MAIN
# Load Sample info
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# Define crop regions
if (CLIMATE_EXTENT == "auto") {
  crop_region <- FUN_define_crop_region(samples, gap = GAP)
  message(paste0('INFO: Extract climate data in region minLongitude=', crop_region[1],
                 '; maxLongitude=', crop_region[2],
                 '; minLatitude=', crop_region[3],
                 '; maxLatitude=', crop_region[4]))
} else if (CLIMATE_EXTENT == "world") {
  crop_region <- 'world'
} else {
  crop_region <- CLIMATE_EXTENT %>% str_split(',') %>% unlist %>% as.numeric
}

message(crop_region %>% str)

# Download and process bioclimatic raster and variables
clim_present <- FUN_present_climate(samples,
                                    DATA_DIR,
                                    crop_region,
                                    RESOLUTION)

message('INFO: Loading climate present data complete, saving to the disk')
message(clim_present %>% str)

# Validate: check for samples with NA climate values (e.g. a coordinate landing on an
# ocean/NoData raster pixel). Always excluded from climate-VALUE-dependent steps only
# (GEA/gradient_forest/geometric_offset/Mantel) -- see filter_climate_valid_samples -- while
# staying in GWAS/phenotype/structure/PCA/sNMF and coordinate-only plotting (piemaps).
# Never a hard stop: the user can fix coordinates or ignore it (see climate_na_excluded.tsv).
site_values <- clim_present$SiteValues
na_rows <- which(rowSums(is.na(site_values)) > 0)

excluded_dt <- data.table(sample = character(), site = character(), latitude = numeric(),
                          longitude = numeric(), reason = character(), distance_km = numeric())

if (length(na_rows) > 0) {
    bad_samples <- samples[na_rows, ]
    reasons <- character(nrow(bad_samples))
    distances <- rep(NA_real_, nrow(bad_samples))
    message("WARNING: Climate extraction returned NA for the following samples:")
    for (i in seq_len(nrow(bad_samples))) {
        is_placeholder <- bad_samples$latitude[i] == 0 && bad_samples$longitude[i] == 0
        reasons[i] <- if (is_placeholder) "likely placeholder (0,0)" else "falls on ocean/NoData pixel"
        dist_msg <- ""
        if (!is_placeholder) {
            # Distance to coast is meaningless for a (0,0) data-entry placeholder -- skip it.
            distances[i] <- FUN_nearest_valid_distance_km(clim_present$RasterStack[[1]],
                                                            bad_samples$longitude[i], bad_samples$latitude[i])
            if (!is.na(distances[i])) dist_msg <- paste0(" (~", distances[i], " km to nearest valid pixel)")
        }
        message(paste0("  - ", bad_samples$sample[i], " (", bad_samples$site[i],
                       ") at (", bad_samples$latitude[i], ", ", bad_samples$longitude[i], "): ", reasons[i], dist_msg))
    }

    warning(paste0("Climate extraction returned NA for ", length(na_rows), " sample(s) -- ",
                   "excluding from climate-VALUE-dependent steps. See climate_na_excluded.tsv. ",
                   "Fix coordinates in input metadata if this is unexpected."))
    message('WARNING: MAF/missingness/biallelic filtering was computed on the full processing-module ',
            'cohort BEFORE this exclusion. Removing samples here can leave those filters inconsistent ',
            'with the actual GEA sample set -- a SNP that passed MAF/missingness genome-wide can become ',
            'monomorphic or zero-variance once these samples are excluded, and is NOT re-filtered. This ',
            'can produce degenerate model fits (e.g. exact p=0/-log10(p)=Inf) for a small number of SNPs, ',
            'especially for methods sensitive to zero-variance genotypes (e.g. EMMAX). Check GEA results ',
            'if they look anomalous for a subset of SNPs.')

    excluded_dt <- data.table(sample = bad_samples$sample, site = bad_samples$site,
                              latitude = bad_samples$latitude, longitude = bad_samples$longitude,
                              reason = reasons, distance_km = distances)

    # Drop the NA rows from site_values AND samples in lockstep, before any table is written,
    # so climate_present_site(.scaled).tsv stay self-consistent (fewer rows, same relative
    # order) with what filter_climate_valid_samples will later produce for the genotype side.
    site_values <- site_values[-na_rows, , drop = FALSE]
    samples <- samples[-na_rows, ]
}

fwrite(excluded_dt, paste0(TABLES_DIR, 'climate_na_excluded.tsv'), sep = '\t', quote = FALSE)

# Save raster
writeRaster(clim_present$RasterStack,
            filename = paste0(RASTER_DIR, 'climate_present_rasterstack.tif'),
            overwrite = TRUE,
            gdal = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

# Save all values
clim_present$AllValues %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_all.tsv'), sep = '\t')

# Save site values (with sample column for traceability)
cbind(sample = samples$sample, site_values) %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_site.tsv'), sep = '\t')

# Save site values scaled (with sample column)
#
# NOTE on the weighting: these rows are per SAMPLE, and a site's climate values are
# repeated for every sample sequenced there, so scale()'s centre and sd are weighted
# by sampling effort rather than by the environment. That is INERT for every current
# consumer of this file, because per-column centring and scaling is an affine
# transform and every consumer is invariant to one: EMMAX/GAPIT regress one trait
# column at a time; lfmm.R re-scales and passes a single column; rda.R and
# pregea_varpart.R (including compute_Px(), which re-scales internally and
# correlates) use the block only through its column span; pregea_rda_setup.R
# correlates it site-wise, and correlation is affine-invariant too. Verified on a
# deliberately unbalanced design (9 sites, n 30/2/2/15/1/1/20/1/1): sample- vs
# site-weighted z-scores give identical lm p-values, projection matrices equal to
# 2e-15, and an identical constrained variance fraction. Re-scaling site-wise would
# rewrite this file and force a re-run of every GEA rule for zero numeric change.
#
# It stops being inert for any consumer that reads predictor SCALE as meaningful —
# a distance between environments, a penalised/ridge fit, a PCA of the predictors.
# Same test: the sample- and site-weighted environmental distance matrices are only
# Spearman 0.95 apart. Do not add such a consumer without switching to a site-level
# centre/sd first. The one distance consumer in the pipeline, mantel_test.R, reads
# the UNSCALED climate_present_site.tsv and does its own scaling — which is exactly
# why it had to be moved to site-level rows there (see that script's own note). (The screen in pregea_rda_setup.R had the related problem via
# repeated ROWS rather than scale — correlation is unaffected by the z-scoring but
# not by the duplication — and is now computed over distinct sites there, as is the
# climate block of plot_correlation_heatmap.R.)
cbind(sample = samples$sample, site_values %>% dplyr::mutate_all(scale)) %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_site_scaled.tsv'), sep = '\t')
