#!/usr/bin/env Rscript
# Unified PieMap plotting script
# Handles both climate-based piemaps and genetic offset piemaps
library(qs)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(terra)
library(ggspatial)
library(ggnewscale)
library(svglite)
library(viridis)

source("/pipeline/scripts/R/lib/piemap_geometry.R")

args = commandArgs(trailingOnly=TRUE)
#################################
RASTER_PATH = args[1]       # path to .tif file
RASTER_LAYER = args[2]      # layer name (e.g., "bio_1") or "1" for single-layer rasters
RASTER_LEGEND = args[3]     # legend title for raster (e.g., "bio_1" or "Genetic Offset")
SAMPLES = args[4]           # metadata file
CLUSTERS = args[5]          # clusters file
SIZE_TRAIT = args[6]        # trait file for pie sizing, or "NULL" for uniform
SIZE_TRAIT_NAME = args[7]   # name for size legend, or "NULL"
PIE_ALPHA = args[8] %>% as.numeric
POP_LABEL = args[9]
POP_LABEL_SIZE = args[10] %>% as.numeric
PLOT_DIR = args[11]
INTER_DIR = args[12]
OUTPUT_PREFIX = args[13]
REGIONMAP_EXTENT = if (length(args) >= 14) args[14] else "NULL"  # Optional: "xmin,xmax,ymin,ymax" or "NULL"
PIE_SCALE = if (length(args) >= 15) as.numeric(args[15]) else 1.0  # Pie size scale factor (1.0 = 100%)
PIE_FILL_LABEL = if (length(args) >= 16) args[16] else "Clusters"  # Legend label for pie slices
USE_POINTS = if (length(args) >= 17) toupper(substr(args[17], 1, 1)) == 'T' else FALSE  # Use points instead of pie charts
#################################

set.seed(42)

message(paste0('INFO: Creating PieMap: ', OUTPUT_PREFIX))

# Color palette for pie charts
# Colorblind-safe (Wong 2011 + Paul Tol), optimized for visibility on yellow/orange climate rasters
my_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
               "#332288", "#882255", "#44AA99", "#AA4499", "#999933",
               "#E69F00", "#661100", "#88CCEE", "#CC6677", "#117733",
               "#DDDDDD", "#855C75", "#D9AF6B", "#736F4C", "#526A83", "#625377")

############################# Functions

#' Calculate autoscaled pie size range based on map extent
#' @param raster_layer SpatRaster layer to get extent from
#' @return List with min_radius, max_radius, mid_radius (in map units)
calculate_pie_size_range <- function(raster_layer) {
  map_extent <- ext(raster_layer)
  x_range <- xmax(map_extent) - xmin(map_extent)
  y_range <- ymax(map_extent) - ymin(map_extent)

  # Use the smaller dimension to ensure pies fit well
  map_size <- min(c(x_range, y_range))

  # Define min/max as proportion of map size
  min_radius <- map_size * 0.02   # 2% of map size
  max_radius <- map_size * 0.06   # 6% of map size
  mid_radius <- (min_radius + max_radius) / 2

  message(paste0('INFO: Map extent - X: ', round(x_range, 2), ', Y: ', round(y_range, 2)))
  message(paste0('INFO: Autoscaled pie sizes - min: ', round(min_radius, 4),
                 ', mid: ', round(mid_radius, 4), ', max: ', round(max_radius, 4)))

  return(list(
    min_radius = min_radius,
    max_radius = max_radius,
    mid_radius = mid_radius
  ))
}

#' Scale trait values to pie radius using min-max normalization
#' @param trait_values Numeric vector of trait values
#' @param min_radius Minimum pie radius
#' @param max_radius Maximum pie radius
#' @return Numeric vector of scaled radii
scale_trait_to_radius <- function(trait_values, min_radius, max_radius) {
  trait_min <- min(trait_values, na.rm = TRUE)
  trait_max <- max(trait_values, na.rm = TRUE)

  # Handle case where all values are the same
  if (trait_max == trait_min) {
    return(rep((min_radius + max_radius) / 2, length(trait_values)))
  }

  # Min-max scaling to [min_radius, max_radius]
  scaled <- min_radius +
    (trait_values - trait_min) / (trait_max - trait_min) *
    (max_radius - min_radius)

  return(scaled)
}

#' Create a PieMap plot
#' @param samples Sample metadata with site, sample, latitude, longitude
#' @param clusters Cluster assignments per sample
#' @param raster_layer SpatRaster layer for background
#' @param trait Optional data.frame with site and trait value columns for pie sizing
#' @param trait_name Name for size legend (e.g., "Tajima's D")
#' @param pie_size_range List with min_radius, max_radius, mid_radius
#' @param pie_alpha Alpha for pie charts
#' @param pop_labels Show population labels
#' @param pop_label_size Size of population labels
#' @param legend_map_title Title for raster legend
#' @return ggplot object
create_piemap <- function(samples,
                          clusters,
                          raster_layer,
                          trait = NULL,
                          trait_name = NULL,
                          pie_size_range = NULL,
                          pie_alpha = 0.6,
                          pop_labels = FALSE,
                          pop_label_size = 5,
                          legend_map_title = 'Value',
                          pie_fill_label = 'Clusters',
                          use_points = FALSE) {

  # Calculate pie size range if not provided
  if (is.null(pie_size_range)) {
    pie_size_range <- calculate_pie_size_range(raster_layer)
  }

  # Get cluster/haplotype columns (everything except sample and site)
  cluster_cols <- setdiff(colnames(clusters), c('sample', 'site'))
  n_clust <- length(cluster_cols)

  # Process cluster info by population: the site-mean slices, drawn once per site at the
  # centroid of its samples (cdb1a1), so jittered per-sample coordinates no longer give a site
  # one pie per distinct coordinate. A site with no coordinate is not drawn.
  clusters_by_pop <- clusters %>%
    dplyr::select(-sample) %>%
    dplyr::group_by(site) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::inner_join(site_centroids(samples), by = 'site')
  no_coord <- setdiff(unique(as.character(clusters$site)), clusters_by_pop$site)
  if (length(no_coord) > 0)
    message('INFO: ', length(no_coord), ' site(s) without a coordinate not drawn: ',
            paste(no_coord, collapse = ', '))

  # Apply pie sizing based on trait or use midpoint for uniform size
  if (is.null(trait) || trait_name == 'NULL' || is.null(trait_name)) {
    # No trait: use midpoint size for all pies
    clusters_by_pop$pie_radius <- pie_size_range$mid_radius
    message('INFO: Using uniform midpoint pie size (no trait provided)')
    trait_name <- NULL
  } else {
    # Trait provided: min-max scale trait values to pie radius
    clusters_by_pop <- clusters_by_pop %>%
      dplyr::left_join(trait, by = 'site') %>%
      dplyr::rename(trait_value = colnames(trait)[2])

    clusters_by_pop$pie_radius <- scale_trait_to_radius(
      clusters_by_pop$trait_value,
      pie_size_range$min_radius,
      pie_size_range$max_radius
    )
    message(paste0('INFO: Scaled ', trait_name, ' values to pie radius [',
                   round(pie_size_range$min_radius, 4), ', ',
                   round(pie_size_range$max_radius, 4), ']'))
  }

  # Build plot with viridis plasma palette (colorblind-friendly)
  gPlot <- ggplot() +
    layer_spatial(raster_layer, aes(fill = after_stat(band1))) +
    scale_fill_viridis_c(option = "plasma", na.value = NA) +
    labs(fill = legend_map_title)

  if (use_points) {
    # Point mode ("clear map"): tiny, barely-visible dots marking sample
    # locations only, so the raster and geography read clearly under dense
    # sampling where pies would overlap. NOTE: this is ~30-400 sample sites,
    # not per-SNP data, so the WGS-density plotting rule (no geom_point at
    # SNP scale) does not apply here.
    # Size = PIE_SCALE * base (1.2pt), floor at 0.5pt to stay clickable.
    point_size <- max(0.5, 1.2 * PIE_SCALE)

    gPlot <- gPlot +
      geom_point(data = clusters_by_pop,
                 aes(x = longitude, y = latitude),
                 shape = 21, fill = 'grey15', color = 'white',
                 size = point_size, stroke = 0.15, alpha = 0.85) +
      theme_bw()

    message(paste0('INFO: Using point mode (size=', round(point_size, 2), ', PIE_SCALE=', PIE_SCALE, ')'))
  } else {
    # Round pies (f3fcce): coord_sf stretches latitude by 1/cos(mid-latitude) on a lon/lat map,
    # so each pie's longitude radius is stretched by the same factor (piemap_geometry.R). The
    # panel's latitude range is the raster's plus the pies'.
    pies <- data.frame(x = clusters_by_pop$longitude, y = clusters_by_pop$latitude,
                       r = clusters_by_pop$pie_radius, group = clusters_by_pop$site)
    pies <- cbind(pies, as.data.frame(clusters_by_pop)[, cluster_cols, drop = FALSE])
    y_range <- range(c(as.vector(terra::ext(raster_layer))[3:4],
                       pies$y - pies$r, pies$y + pies$r), na.rm = TRUE)
    x_stretch <- pie_x_stretch(y_range, terra::is.lonlat(raster_layer, warn = FALSE))
    wedges <- pie_wedges(pies, cluster_cols, x_stretch)
    message('INFO: Pie longitude radius x', round(x_stretch, 4),
            ' (mid-latitude ', round(mean(y_range), 3), ')')

    gPlot <- gPlot +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = wedges,
                   aes(x = x, y = y, group = wedge, fill = slice),
                   color = 'black', linewidth = 0.4,
                   alpha = pie_alpha) +
      scale_fill_manual(values = setNames(my_colors[seq_len(n_clust)], cluster_cols),
                        breaks = cluster_cols, drop = FALSE) +
      theme_bw() +
      labs(fill = pie_fill_label, x = NULL, y = NULL)   # the polygons' x/y, as scatterpie drew none
  }

  # Add population labels
  if (pop_labels) {
    label_data <- clusters_by_pop %>%
      dplyr::distinct(site, .keep_all = TRUE)
    gPlot <- gPlot +
      geom_text_repel(data = label_data,
                      aes(x = longitude, y = latitude, label = site),
                      alpha = 0.7,
                      color = 'black',
                      size = pop_label_size,
                      box.padding = 0.5,
                      point.padding = 0.3,
                      segment.color = 'grey50',
                      segment.alpha = 0.5,
                      max.overlaps = Inf)
  }

  # Add size legend showing actual trait values
  if (!is.null(trait_name)) {
    trait_vals <- clusters_by_pop$trait_value
    trait_min <- min(trait_vals, na.rm = TRUE)
    trait_max <- max(trait_vals, na.rm = TRUE)
    trait_mid <- (trait_min + trait_max) / 2

    # Create break values and corresponding point sizes (in mm for geom_point)
    size_min <- 2
    size_max <- 6
    size_mid <- (size_min + size_max) / 2

    legend_data <- data.frame(
      x = rep(mean(clusters_by_pop$longitude), 3),
      y = rep(mean(clusters_by_pop$latitude), 3),
      trait_value = c(trait_min, trait_mid, trait_max),
      point_size = c(size_min, size_mid, size_max)
    )

    # Format legend values based on trait type
    # Deduplicate when min == max (all values identical)
    break_vals <- unique(c(trait_min, trait_mid, trait_max))
    if (grepl("Pi|pi|diversity", trait_name, ignore.case = TRUE)) {
      break_labels <- formatC(break_vals, format = "e", digits = 2)
    } else {
      break_labels <- formatC(break_vals, format = "f", digits = 2)
    }

    gPlot <- gPlot +
      geom_point(data = legend_data,
                 aes(x = x, y = y, size = trait_value),
                 alpha = 0) +
      scale_size_continuous(
        name = trait_name,
        range = c(size_min, size_max),
        breaks = break_vals,
        labels = break_labels
      ) +
      guides(size = guide_legend(override.aes = list(alpha = 1, color = 'grey30')))

    message(paste0('INFO: Added size legend for ', trait_name,
                   ' [', round(trait_min, 3), ', ', round(trait_max, 3), ']'))
  }

  return(gPlot)
}

################################# MAIN

# Load raster
raster_data <- rast(RASTER_PATH)

# Extract layer
if (nlyr(raster_data) > 1) {
  # Multi-layer raster (e.g., climate stack)
  if (!RASTER_LAYER %in% names(raster_data)) {
    stop(paste0("ERROR: Layer '", RASTER_LAYER, "' not found in raster. Available: ",
                paste(names(raster_data), collapse = ", ")))
  }
  rlayer <- raster_data[[RASTER_LAYER]]
  message(paste0('INFO: Using layer "', RASTER_LAYER, '" from multi-layer raster'))
} else {
  # Single-layer raster (e.g., genetic offset)
  rlayer <- raster_data
  message('INFO: Using single-layer raster')
}

# Load data
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# sample stays character too (2d4533): the zoomed maps filter clusters by sample, and an id such
# as 001 read as the number 1 matched none, so the zoom drew no pies.
cl_classes <- c(site = "character")
if ("sample" %in% names(fread(CLUSTERS, nrows = 0))) cl_classes["sample"] <- "character"
clusters <- fread(CLUSTERS, colClasses = cl_classes)

# Load size trait if provided
if (SIZE_TRAIT != "NULL" && file.exists(SIZE_TRAIT)) {
  size_trait <- fread(SIZE_TRAIT, colClasses = c("site" = "character"))
  message(paste0('INFO: Loaded size trait from ', SIZE_TRAIT))
} else {
  size_trait <- NULL
  message('INFO: No size trait provided, using uniform pie size')
}

# Parse population label setting
pop_labels <- toupper(substr(POP_LABEL, 1, 1)) == 'T'

# Detect near-zero spatial variance of the climate variable at sampling sites
site_coords <- unique(samples[, .(site, longitude, latitude)])
site_extracted <- terra::extract(rlayer, cbind(site_coords$longitude, site_coords$latitude))
site_climate_vals <- site_extracted[, 1]
site_climate_vals <- site_climate_vals[!is.na(site_climate_vals)]
climate_range <- if (length(site_climate_vals) >= 2) diff(range(site_climate_vals)) else 0
climate_mean  <- if (length(site_climate_vals) >= 1) mean(site_climate_vals) else 0
has_spatial_variance <- climate_range > max(0.01, abs(climate_mean) * 0.01)

if (!has_spatial_variance) {
  flag_file <- file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, "_no_spatial_variance.flag"))
  writeLines(c(
    paste0("variable=", RASTER_LEGEND),
    paste0("range=",    round(climate_range, 6)),
    paste0("mean=",     round(climate_mean,  6))
  ), flag_file)
  message(paste0("WARNING: Climate variable '", RASTER_LEGEND,
                 "' shows near-zero spatial variance across sampling sites ",
                 "(range=", round(climate_range, 6), "). ",
                 "Piemap generated but may not be informative. ",
                 "Flag written: ", flag_file))
}

# Calculate pie size range based on map extent
pie_size_range <- calculate_pie_size_range(rlayer)

# Apply user-defined scale factor to pie sizes
if (PIE_SCALE != 1.0) {
    message(paste0('INFO: Applying pie scale factor: ', PIE_SCALE, ' (', PIE_SCALE * 100, '%)'))
    pie_size_range$min_radius <- pie_size_range$min_radius * PIE_SCALE
    pie_size_range$max_radius <- pie_size_range$max_radius * PIE_SCALE
    pie_size_range$mid_radius <- pie_size_range$mid_radius * PIE_SCALE
    message(paste0('INFO: Scaled pie sizes - min: ', round(pie_size_range$min_radius, 4),
                   ', mid: ', round(pie_size_range$mid_radius, 4),
                   ', max: ', round(pie_size_range$max_radius, 4)))
}

# Create plot
gPlot <- create_piemap(
  samples = samples,
  clusters = clusters,
  raster_layer = rlayer,
  trait = size_trait,
  trait_name = if (SIZE_TRAIT_NAME != "NULL") SIZE_TRAIT_NAME else NULL,
  pie_size_range = pie_size_range,
  pie_alpha = PIE_ALPHA,
  pop_labels = pop_labels,
  pop_label_size = POP_LABEL_SIZE,
  legend_map_title = RASTER_LEGEND,
  pie_fill_label = PIE_FILL_LABEL,
  use_points = USE_POINTS
)

# Save outputs (full extent)
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, '.png')), gPlot)
ggsave(file.path(PLOT_DIR, paste0(OUTPUT_PREFIX, '.svg')), gPlot,
       device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
qsave(gPlot, file.path(INTER_DIR, paste0(OUTPUT_PREFIX, '.qs')))

message(paste0('INFO: PieMap complete: ', OUTPUT_PREFIX))

# Generate zoomed plot if REGIONMAP_EXTENT is provided
if (REGIONMAP_EXTENT != "NULL" && REGIONMAP_EXTENT != "") {
    message("INFO: Generating zoomed regionmap plot...")

    # Parse extent coordinates: "xmin,xmax,ymin,ymax"
    coords <- as.numeric(strsplit(REGIONMAP_EXTENT, ",")[[1]])
    if (length(coords) != 4) {
        stop("ERROR: REGIONMAP_EXTENT must be 'xmin,xmax,ymin,ymax', got: ", REGIONMAP_EXTENT)
    }

    xmin_val <- coords[1]
    xmax_val <- coords[2]
    ymin_val <- coords[3]
    ymax_val <- coords[4]

    # Create extent object for cropping
    zoom_extent <- ext(xmin_val, xmax_val, ymin_val, ymax_val)

    # Crop raster to zoom extent
    rlayer_zoom <- crop(rlayer, zoom_extent)

    # Filter samples to those within zoom extent
    samples_zoom <- samples %>%
        filter(longitude >= xmin_val & longitude <= xmax_val) %>%
        filter(latitude >= ymin_val & latitude <= ymax_val)

    if (nrow(samples_zoom) == 0) {
        message("WARNING: No samples within zoom extent, skipping zoomed plot")
    } else {
        message(paste0("INFO: ", nrow(samples_zoom), " samples within zoom extent"))

        # Filter clusters to match zoomed samples
        clusters_zoom <- clusters %>% filter(sample %in% samples_zoom$sample)

        # Filter size trait if provided
        size_trait_zoom <- NULL
        if (!is.null(size_trait)) {
            size_trait_zoom <- size_trait %>% filter(site %in% samples_zoom$site)
        }

        # Recalculate pie size range for zoomed extent (maintains visual proportions)
        pie_size_range_zoom <- calculate_pie_size_range(rlayer_zoom)

        # Apply user-defined scale factor to zoomed pie sizes
        if (PIE_SCALE != 1.0) {
            pie_size_range_zoom$min_radius <- pie_size_range_zoom$min_radius * PIE_SCALE
            pie_size_range_zoom$max_radius <- pie_size_range_zoom$max_radius * PIE_SCALE
            pie_size_range_zoom$mid_radius <- pie_size_range_zoom$mid_radius * PIE_SCALE
            message(paste0('INFO: Applied pie scale factor to zoomed plot: ', PIE_SCALE, ' (', PIE_SCALE * 100, '%)'))
        }

        # Create zoomed plot
        gPlot_zoom <- create_piemap(
            samples = samples_zoom,
            clusters = clusters_zoom,
            raster_layer = rlayer_zoom,
            trait = size_trait_zoom,
            trait_name = if (SIZE_TRAIT_NAME != "NULL") SIZE_TRAIT_NAME else NULL,
            pie_size_range = pie_size_range_zoom,
            pie_alpha = PIE_ALPHA,
            pop_labels = pop_labels,  # Keep same label setting
            pop_label_size = POP_LABEL_SIZE,
            legend_map_title = paste0(RASTER_LEGEND, " (Zoomed)"),
            pie_fill_label = PIE_FILL_LABEL,
            use_points = USE_POINTS
        )

        # Create coordinate subdirectory for zoomed plots
        coords_str <- gsub("\\.", "p", REGIONMAP_EXTENT)  # Replace . with p
        coords_str <- gsub(",", "_", coords_str)          # Replace , with _

        # Save zoomed outputs in coordinate subdirectory
        zoom_dir <- file.path(PLOT_DIR, "zoom", coords_str)
        dir.create(zoom_dir, recursive = TRUE, showWarnings = FALSE)
        ggsave(file.path(zoom_dir, paste0(OUTPUT_PREFIX, '.png')), gPlot_zoom)
        ggsave(file.path(zoom_dir, paste0(OUTPUT_PREFIX, '.svg')), gPlot_zoom,
               device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
        qsave(gPlot_zoom, file.path(INTER_DIR, paste0(OUTPUT_PREFIX, '_zoom_', coords_str, '.qs')))

        message(paste0('INFO: Zoomed PieMap complete: ', file.path(zoom_dir, OUTPUT_PREFIX)))
    }
}
