#!/usr/bin/env Rscript
# run_haplotype_viz.R
# Generate crosshap_viz, phenotype boxplots, haplotype piemaps at selected epsilon
#
# Args:
#   1  selected_regions  - Regions TSV
#   2  metadata_file     - Aligned metadata (site, sample, lat, lon, traits...)
#   3  epsilon_selected  - Selected epsilon value
#   4  metadata_type     - 'site' or 'cluster_K{N}'
#   5  tag               - Combined tag for filenames
#   6  intermediate_dir  - Dir with .qs HapObjects
#   7  plots_dir         - Output dir for plots
#   8  tables_dir        - Output dir for tables
#   9  raster_path       - Climate raster .grd for piemap background
#  10  raster_layer      - Layer name (e.g., "bio_1")
#  11  piemap_args       - Comma-separated: alpha,pop_label,pop_label_size,pie_scale
#  12  regionmap_extent  - "xmin,xmax,ymin,ymax" or "NULL"
#  13  clusters_file     - Cluster assignments TSV (for ancestry piemaps) or 'NULL'

args <- commandArgs(trailingOnly = TRUE)
selected_regions <- args[1]
metadata_file    <- args[2]
epsilon_selected <- as.numeric(args[3])
metadata_type    <- args[4]
tag              <- args[5]
intermediate_dir <- args[6]
plots_dir        <- args[7]
tables_dir       <- args[8]
raster_path      <- args[9]
raster_layer_name <- args[10]
piemap_args      <- strsplit(args[11], ",")[[1]]
regionmap_extent <- if (length(args) >= 12) args[12] else "NULL"
clusters_file    <- if (length(args) >= 13) args[13] else "NULL"

pie_alpha     <- as.numeric(piemap_args[1])
pop_label     <- piemap_args[2]
pop_label_size <- as.numeric(piemap_args[3])
pie_scale     <- if (length(piemap_args) >= 4) as.numeric(piemap_args[4]) else 1.0
use_points    <- if (length(piemap_args) >= 5) piemap_args[5] else "F"

library(data.table)
library(crosshap)
library(qs)
library(ggplot2)
library(ggdist)
library(ggpubr)
library(svglite)

message("=== Haplotype Visualization ===")
message("Epsilon: ", epsilon_selected)
message("Tag: ", tag)

# Read regions and metadata
regions <- fread(selected_regions)
meta <- fread(metadata_file, colClasses = c("site" = "character", "sample" = "character"))
sample_col <- names(meta)[2]
site_col   <- names(meta)[1]
trait_cols  <- names(meta)[-(1:4)]   # after site/sample/latitude/longitude (1d4be2); 5:ncol named longitude a trait when there was none

message("Regions: ", nrow(regions))
message("Traits: ", paste(trait_cols, collapse = ", "))

# Process each region
all_assignments <- list()

for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  rid <- reg$region_id

  hap_file <- file.path(intermediate_dir, paste0("Region_", rid, "_HapObject.qs"))
  if (!file.exists(hap_file)) {
    message("WARNING: No HapObject for region ", rid, ", skipping")
    next
  }

  message("\n--- Region ", rid, " ---")
  hap_obj <- qread(hap_file)

  # === Extract haplotype assignments ===
  # HapObject is a nested list: hap_obj$Haplotypes_MGmin{N}_E{eps}$Indfile
  # Find the element matching the selected epsilon
  indfile <- NULL

  tryCatch({
    eps_pattern <- paste0("_E", epsilon_selected, "$")
    eps_names <- grep(eps_pattern, names(hap_obj), value = TRUE)
    if (length(eps_names) > 0) {
      eps_data <- hap_obj[[eps_names[1]]]
      message("Accessing HapObject element: ", eps_names[1])
      if (!is.null(eps_data$Indfile) && is.data.frame(eps_data$Indfile)) {
        indfile <- as.data.table(eps_data$Indfile)
        message("Extracted Indfile: ", nrow(indfile), " individuals, columns: ",
                paste(names(indfile), collapse = ", "))
      }
    } else {
      message("WARNING: No results for epsilon=", epsilon_selected,
              " in HapObject. Available: ", paste(names(hap_obj), collapse = ", "))
    }
  }, error = function(e) {
    message("WARNING: Could not extract Indfile for region ", rid, ": ", e$message)
  })

  if (is.null(indfile) || nrow(indfile) == 0) {
    message("WARNING: No haplotype assignments for region ", rid, ", skipping downstream")
    next
  }

  # Standardize Indfile columns
  # crosshap Indfile typically has: Ind, hapgroup (or similar)
  ind_col <- grep("^Ind$|^ind$|^Individual$", names(indfile), value = TRUE)
  hap_col <- grep("^hap|^Hap|^group|^Group|^cluster", names(indfile), value = TRUE)

  if (length(ind_col) == 0) ind_col <- names(indfile)[1]
  if (length(hap_col) == 0) hap_col <- names(indfile)[ncol(indfile)]

  # crosshap uses periods in sample names (replaces dashes), revert for matching
  indfile[[ind_col[1]]] <- gsub("\\.", "-", indfile[[ind_col[1]]])

  assignments <- data.table(
    sample = indfile[[ind_col[1]]],
    haplotype = paste0("Hap_", indfile[[hap_col[1]]]),
    region_id = rid
  )

  # Join with metadata
  assignments <- merge(assignments, meta, by.x = "sample", by.y = sample_col, all.x = TRUE)

  all_assignments[[as.character(rid)]] <- assignments

  # === crosshap_viz + boxplots per trait ===
  for (trait in trait_cols) {
    # crosshap_viz: update Pheno column in HapObject for this trait, then render
    tryCatch({
      message("  Generating crosshap_viz for trait '", trait, "'...")
      hap_obj_trait <- hap_obj
      for (eps_name in names(hap_obj_trait)) {
        if (!is.null(hap_obj_trait[[eps_name]]$Indfile)) {
          # crosshap uses periods in sample names; match against assignments
          ind_names <- gsub("\\.", "-", hap_obj_trait[[eps_name]]$Indfile$Ind)
          trait_vals <- meta[[trait]][match(ind_names, meta[[sample_col]])]
          hap_obj_trait[[eps_name]]$Indfile$Pheno <- trait_vals
        }
      }
      p_viz <- crosshap_viz(hap_obj_trait, epsilon = epsilon_selected)
      p_viz <- p_viz & ggplot2::labs(colour = trait, fill = trait)
      viz_file <- file.path(plots_dir, paste0("Region_", rid, "_crosshap_viz_", trait))
      qsave(p_viz, paste0(viz_file, ".qs"))
      ggsave(paste0(viz_file, ".png"), p_viz, width = 14, height = 8, dpi = 150)
      ggsave(paste0(viz_file, ".svg"), p_viz, width = 14, height = 8)
      message("  Saved crosshap_viz: ", viz_file)
    }, error = function(e) {
      message("  WARNING: crosshap_viz failed for trait '", trait, "': ", e$message)
    })

    # Boxplot
    tryCatch({
      plot_data <- assignments[!is.na(get(trait))]
      if (nrow(plot_data) < 3) {
        message("  Skipping boxplot for trait '", trait, "' (too few samples)")
        next
      }

      n_haps <- length(unique(plot_data$haplotype))
      if (n_haps < 2) {
        message("  Skipping boxplot for trait '", trait, "' (only 1 haplotype group)")
        next
      }

      p_box <- ggplot(plot_data, aes(x = haplotype, y = .data[[trait]], fill = haplotype)) +
        ggdist::stat_halfeye(
          adjust = 0.5, width = 0.6, .width = 0, justification = -0.2,
          point_colour = NA, alpha = 0.6
        ) +
        geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.05, alpha = 0.5, size = 1.5) +
        labs(
          title = paste0("Region ", rid, " - ", trait),
          x = "Haplotype Group", y = trait,
          fill = "Haplotype"
        ) +
        theme_minimal() +
        theme(legend.position = "none")

      if (n_haps == 2) {
        p_box <- p_box + stat_compare_means(method = "wilcox.test", label = "p.format")
      } else if (n_haps > 2) {
        p_box <- p_box + stat_compare_means(method = "kruskal.test", label = "p.format")
      }

      box_file <- file.path(plots_dir, paste0("Region_", rid, "_boxplot_", trait))
      qsave(p_box, paste0(box_file, ".qs"))
      ggsave(paste0(box_file, ".png"), p_box, width = 8, height = 6, dpi = 150)
      ggsave(paste0(box_file, ".svg"), p_box, width = 8, height = 6)
    }, error = function(e) {
      message("  WARNING: Boxplot failed for trait '", trait, "': ", e$message)
    })
  }

  # === Haplotype piemaps (one per trait) ===
  # Create clusters-format file for plot_piemap.R reuse:
  # columns: sample, site, C1, C2, ... (one-hot encoded haplotype groups)
  tryCatch({
    hap_levels <- sort(unique(assignments$haplotype))
    n_haps <- length(hap_levels)

    if (n_haps >= 2) {
      # Create one-hot encoding
      hap_clusters <- data.table(
        sample = assignments$sample,
        site = assignments[[site_col]]
      )
      for (h in seq_along(hap_levels)) {
        hap_clusters[, (hap_levels[h]) := as.integer(assignments$haplotype == hap_levels[h])]
      }

      # Write clusters file for piemap
      clusters_out <- file.path(tables_dir, paste0("Region_", rid, "_haplotype_clusters.tsv"))
      fwrite(hap_clusters, clusters_out, sep = "\t", quote = FALSE)

      # Generate piemaps per trait using plot_piemap.R
      for (trait in trait_cols) {
        # Create site-level trait means file
        trait_data <- assignments[!is.na(get(trait)), .(value = mean(get(trait), na.rm = TRUE)),
                                  by = site_col]
        setnames(trait_data, c("site", trait))

        trait_file <- file.path(tables_dir, paste0("Region_", rid, "_", trait, "_site_means.tsv"))
        fwrite(trait_data, trait_file, sep = "\t", quote = FALSE)

        # Call plot_piemap.R
        piemap_prefix <- paste0("haplotype_piemap_region_", rid, "_", trait)
        cmd <- paste0(
          "Rscript /pipeline/scripts/plot_piemap.R ",
          raster_path, " ", raster_layer_name, " ", raster_layer_name, " ",
          metadata_file, " ", clusters_out, " ",
          trait_file, " \"", trait, "\" ",
          pie_alpha, " ", pop_label, " ", pop_label_size, " ",
          plots_dir, " ", intermediate_dir, " ",
          piemap_prefix, " ", regionmap_extent, " ", pie_scale, " Haplotypes ", use_points
        )
        message("  Generating piemap: ", piemap_prefix)
        system(cmd)
      }
    }
  }, error = function(e) {
    message("WARNING: Piemap generation failed for region ", rid, ": ", e$message)
  })
}

# === Write combined haplotype assignment table ===
if (length(all_assignments) > 0) {
  combined <- rbindlist(all_assignments, fill = TRUE)
  assign_file <- file.path(tables_dir, "Haplotype_assignments.tsv")
  fwrite(combined, assign_file, sep = "\t", quote = FALSE)
  message("\nSaved haplotype assignments: ", assign_file, " (", nrow(combined), " rows)")

  # Write haplotype frequency per site per region
  freq_list <- list()
  for (rid in names(all_assignments)) {
    a <- all_assignments[[rid]]
    if (nrow(a) == 0) next
    freq <- a[, .N, by = .(get(site_col), haplotype, region_id)]
    setnames(freq, "get", "site")
    total <- a[, .N, by = .(get(site_col), region_id)]
    setnames(total, c("site", "region_id", "total"))
    freq <- merge(freq, total, by = c("site", "region_id"))
    freq[, frequency := N / total]
    freq_list[[rid]] <- freq
  }
  if (length(freq_list) > 0) {
    freq_all <- rbindlist(freq_list)
    freq_file <- file.path(tables_dir, "Haplotype_frequencies.tsv")
    fwrite(freq_all, freq_file, sep = "\t", quote = FALSE)
    message("Saved haplotype frequencies: ", freq_file)
  }
}

message("\n=== Haplotype visualization complete ===")
