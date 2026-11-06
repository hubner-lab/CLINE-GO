#!/usr/bin/env Rscript
# k_evidence_plot.R: the figures behind the choice of sNMF.k_best.
#
# Reads ANALYSIS/Arabidopsis1135/k_evidence/, written by k_evidence.py, and the pipeline's
# cross-entropy table. Writes into the same folder:
#   k_evidence_curve.png  cross-entropy of every run, its drop from K-1 against the run-to-run
#                         spread, and how many of the 10 runs reach the best run's solution
#   k_evidence_map.png    accessions by main cluster at each K of the command line
#
# The pipeline's cross-entropy plot draws the best run per K and its structure plots keep metadata
# order, so neither shows spread, reproducibility or geography.
#
# Run from the repo root, in the image:
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp -v $PWD:/data -w /data \
#     cline-go:6cdedcc Rscript manuscript/arabidopsis/k_evidence_plot.R 6 11

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(terra)
  library(here)
})
source("/pipeline/scripts/R/utils/theme_clinego.R")

map_k <- as.integer(commandArgs(trailingOnly = TRUE))
if (length(map_k) == 0 || anyNA(map_k)) stop("k_evidence_plot.R: name the K values to map")

evidence_dir <- here("ANALYSIS", "Arabidopsis1135", "k_evidence")
read_evidence <- function(name, ...) fread(file.path(evidence_dir, name), ...)

curve <- read_evidence("curve.tsv")
solutions <- read_evidence("solutions.tsv")
clusters <- read_evidence("clusters.tsv")
assignments <- read_evidence("assignments.tsv", colClasses = c(sample = "character"))
runs <- fread(here("Arabidopsis1135_results", "PreStructure", "tables", "cross_entropy_K2-12.tsv"))

k_axis <- scale_x_continuous(breaks = curve$K, name = "K")

plot_cross_entropy <- ggplot(mapping = aes(K)) +
  geom_line(data = curve, aes(y = mean), color = CLINEGO_NEUTRAL) +
  geom_point(data = runs, aes(y = cross_entropy), color = CLINEGO_NEUTRAL, size = 1.2) +
  k_axis +
  labs(y = "Cross-entropy, all 10 runs") +
  theme_clinego()

pooled_sd <- curve[, sqrt(mean(sd^2))]
plot_step <- ggplot(curve[!is.na(step)], aes(K, step)) +
  geom_col(fill = CLINEGO_NEUTRAL, width = 0.7) +
  geom_hline(yintercept = pooled_sd, color = CLINEGO_THRESHOLD, linetype = "dashed") +
  k_axis +
  labs(y = "Drop in mean cross-entropy from K-1") +
  theme_clinego()

plot_solution <- ggplot(solutions[holds_best_run == TRUE], aes(K, n_runs)) +
  geom_col(fill = CLINEGO_NEUTRAL, width = 0.7) +
  k_axis +
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = c(0, 0.02))) +
  labs(y = "Runs reaching the best run's solution") +
  theme_clinego()

ggsave(file.path(evidence_dir, "k_evidence_curve.png"),
       plot_cross_entropy + plot_step + plot_solution + plot_annotation(tag_levels = "a"),
       width = 13, height = 4.2, dpi = 200)

# Land from the BIO1 layer, sampled to 0.25 degrees: the image has no vector coastline.
xlim <- c(-30, 145)
ylim <- c(10, 70)
bio1 <- rast(here("data", "wc2.1_30s", "wc2.1_30s_bio_1.tif"))
template <- rast(ext(xlim[1], xlim[2], ylim[1], ylim[2]), resolution = 0.25, crs = crs(bio1))
land <- as.data.frame(resample(bio1, template, method = "near"), xy = TRUE, na.rm = TRUE)
setDT(land)

two_countries <- function(countries) {
  vapply(strsplit(countries, ", ", fixed = TRUE),
         function(parts) paste(head(parts, 2), collapse = ", "), character(1))
}

plot_k_map <- function(k) {
  labelled <- clusters[K == k][order(-n_main)]
  labelled[, label := sprintf("%s  n = %d: %s", cluster, n_main, two_countries(countries))]
  points <- assignments[, .(longitude, latitude, cluster = get(paste0("K", k)))]
  points[, cluster := factor(cluster, levels = labelled$cluster)]
  ggplot() +
    geom_raster(data = land, aes(x, y), fill = "grey90") +
    geom_point(data = points, aes(longitude, latitude, color = cluster), size = 1) +
    scale_color_manual(values = setNames(clinego_cluster_palette(k), labelled$cluster),
                       labels = setNames(labelled$label, labelled$cluster), name = NULL) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.6))) +
    coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = sprintf("K = %d", k), x = NULL, y = NULL) +
    theme_clinego() +
    theme(legend.position = "right", legend.key.height = unit(0.5, "cm"))
}

ggsave(file.path(evidence_dir, "k_evidence_map.png"),
       wrap_plots(lapply(map_k, plot_k_map), ncol = 1),
       width = 14, height = 5.6 * length(map_k), dpi = 200)
message("k_evidence_plot.R: figures in ", evidence_dir)
