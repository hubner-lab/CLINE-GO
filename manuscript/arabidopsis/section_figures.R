#!/usr/bin/env Rscript
# section_figures.R: over-complete figure drafts for the Arabidopsis section, for the user to prune
# (plan previously-we-made-the-zazzy-wadler, Phase 4). Every candidate panel is drawn, and one assembled
# main-figure draft holds them all; the final figures are rendered after the user's choices.
#
#   Rscript manuscript/arabidopsis/section_figures.R [project]      (default Arabidopsis1135it3)
#
# Reads ANALYSIS/Arabidopsis1135/offset/<project>/{rho_scale,garden_means}.tsv (rho_scale.R,
# offset_eval.R), offset/<project>/survival/garden_survival.tsv (survival_eval.R),
# data/arabidopsis1135/garden_sites.tsv, ANALYSIS/Arabidopsis1135/time_memory/time_memory.tsv
# (collect_resources.R), and the PNGs of the K map, the Miami plots and the pipeline's LD decay plot.
# Writes ANALYSIS/<project>/section/figures/drafts/<stem>.{png,svg}:
#   rho_scale            mean rho (bars: IQR) of every predictor: GrENE-net's as its text states them, its
#                        replicate prediction from Table S12, our climate distance, 18 runs and the draws
#   level3               the union against its 20 random draws per method: (a) mean rho, (b) per garden,
#                        the union minus the draws' mean
#   garden_gf_climate    the 30 garden means: Gradient Forest union against climate distance
#   survival             garden offset (and climate distance) against the share of fitness plots alive in
#                        year 3, 29 gardens
#   cost_compact         wall time and peak memory per mode, iteration 3 and the stress test only
#   design_map           accessions by sNMF cluster at the project's k_best, the GrENE-net founders outlined,
#                        and the 30 gardens
#   main_assembled       one draft of the main figure with every candidate panel (PNG only)
# and the final figures after the user's review (2026-11-15) into ANALYSIS/<project>/section/figures/:
#   main_arabidopsis     a design map with Europe enlarged beside it; b the Miami redrawn on the pipeline's
#                        own background cloud and overlay coordinates (miami_with_rda.sh <project> rda-second),
#                        with the classic genes the scans found marked; c mean rho of GrENE-net's genomic
#                        offset beside CLINE-GO's three offset methods on the union set; d per garden,
#                        CLINE-GO's Gradient Forest accuracy against GrENE-net's genomic offset accuracy
#   main_panel_a to main_panel_d   the panels alone
#   s_rho_scale_all      mean rho of every set per method beside GrENE-net's predictors
#   s_agreement          per garden, each method's accuracy against GrENE-net's offset (union set), with
#                        Spearman's rho and p; every set's agreement is a table, not a panel
#   s_cost (drawn here), s_ld_decay (a copy of the pipeline's LD decay plot); no survival figure (user, 2026-11-15)
# No random-draw comparison and no climate distance in the final figures (user, 2026-11-15).
# Colours: theme_clinego.R's constants; Okabe-Ito for categories, the map gradient for garden BIO1.
# Run from the repo root, in the image.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(here)
  library(terra)
})
source(here("scripts", "R", "utils", "theme_clinego.R"))

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1L) args[1] else "Arabidopsis1135it3"
offset_dir <- here("ANALYSIS", "Arabidopsis1135", "offset", project)
draft_dir <- here("ANALYSIS", project, "section", "figures", "drafts")
dir.create(draft_dir, recursive = TRUE, showWarnings = FALSE)

method_label <- c(gradient_forest = "Gradient Forest", geometric_offset = "Geometric offset",
                  rda_offset = "RDA offset")
set_order <- c("gea.union", "gea.emmax", "gea.lfmm", "gea.rda", "gea.two_of_three", "gea.three_of_three")
set_label <- c(gea.union = "union set", gea.emmax = "EMMAX set", gea.lfmm = "LFMM set", gea.rda = "RDA set",
               gea.two_of_three = "2-of-3 set", gea.three_of_three = "3-of-3 set")
draw_pattern <- "^random\\.(.+)\\.([0-9]+)$"
experiment_colors <- c("Trained on GrENE-net's evolved populations" = CLINEGO_CATEGORICAL[1],
                       "No GrENE-net data" = CLINEGO_CATEGORICAL[5])
out_of_range <- c(26L, 60L)

save_draft <- function(stem, plot, w, h) {
  clinego_save_both(file.path(draft_dir, stem), plot, w = w, h = h, dpi = 200)
  message("wrote ", stem)
}

rho_scale <- fread(file.path(offset_dir, "rho_scale.tsv"))
garden_means <- fread(file.path(offset_dir, "garden_means.tsv"))
gardens <- fread(here("data", "arabidopsis1135", "garden_sites.tsv"), select = c("garden", "bio_1"))
setnames(gardens, "bio_1", "bio_garden")
climate_rho <- rho_scale[predictor == "ours.climate_distance", mean]

# Rho scale ----------------------------------------------------------------------------------------------

singles <- rho_scale[source != "CLINE-GO random draw"]
singles[, group := fcase(startsWith(source, "GrENE-net"), "GrENE-net, leave-one-out",
                         source == "CLINE-GO baseline", "CLINE-GO, climate",
                         default = paste("CLINE-GO,", method_label[method]))]
singles[, row_label := fcase(
  predictor == "grene.replicate_prediction", "Replicate prediction (Table S12)",
  predictor == "grene.genomic_offset_gaussian", "Genomic offset + Gaussian parameters (text)",
  predictor == "grene.genomic_offset", "Genomic offset, BIO1 (text)",
  predictor == "grene.climate_distance", "Climate distance, BIO1, ERA5 (text)",
  source == "CLINE-GO baseline", "Climate distance, BIO1, WorldClim",
  default = set_label[set])]
singles[, order := fcase(predictor == "grene.replicate_prediction", 1L,
                         predictor == "grene.genomic_offset_gaussian", 2L,
                         predictor == "grene.genomic_offset", 3L,
                         predictor == "grene.climate_distance", 4L,
                         default = 10L + match(set, set_order, nomatch = 0L))]
draws <- rho_scale[source == "CLINE-GO random draw"]
draws[, group := paste("CLINE-GO,", method_label[method])]
draws[, row_label := "20 random draws, union size"]
draws[, order := 20L]
draw_range <- draws[, .(lo = min(mean), hi = max(mean), mean = mean(mean)), by = .(group, row_label, order,
                                                                                  uses_experiment_data)]
group_levels <- c("GrENE-net, leave-one-out", "CLINE-GO, climate", paste("CLINE-GO,", method_label))
group_strip <- setNames(c("GrENE-net\nleave-one-out", "CLINE-GO climate", "CLINE-GO\nGradient Forest",
                          "CLINE-GO\ngeometric offset", "CLINE-GO\nRDA offset"), group_levels)
level_table <- unique(rbind(singles[, .(group, row_label, order)], draws[, .(group, row_label, order)]))
setorder(level_table, order)
row_levels <- rev(unique(level_table$row_label))
for (dt in list(singles, draws, draw_range)) {
  dt[, group := factor(group, levels = group_levels)]
  dt[, row_label := factor(row_label, levels = row_levels)]
  dt[, experiment := factor(fifelse(uses_experiment_data == "yes", names(experiment_colors)[1],
                                    names(experiment_colors)[2]), levels = names(experiment_colors))]
}

rho_plot <- function(s, d, r, base_size = 10) {
  ggplot() +
    geom_vline(xintercept = 0, colour = CLINEGO_COL$muted, linewidth = 0.3) +
    geom_vline(xintercept = climate_rho, colour = CLINEGO_THRESHOLD, linewidth = 0.4, linetype = "dashed") +
    geom_errorbar(data = s, aes(y = row_label, xmin = q1, xmax = q3, colour = experiment),
                  width = 0.25, linewidth = 0.4, orientation = "y") +
    geom_point(data = s, aes(y = row_label, x = mean, colour = experiment, shape = source), size = 2.4) +
    geom_segment(data = r, aes(y = row_label, yend = row_label, x = lo, xend = hi, colour = experiment),
                 linewidth = 2.2, alpha = 0.35) +
    geom_point(data = d, aes(y = row_label, x = mean, colour = experiment), size = 1, alpha = 0.7,
               position = position_jitter(height = 0.12, width = 0, seed = 1)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(group_strip)) +
    scale_colour_manual(values = experiment_colors, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c("GrENE-net text" = 18, "GrENE-net Table S12" = 5, "CLINE-GO run" = 16,
                                  "CLINE-GO baseline" = 15), guide = "none") +
    labs(x = "Spearman's rho with first-year fitness, mean of 325 plots (bars: IQR)", y = NULL) +
    theme_clinego_grid(base_size = base_size) +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
}
p_rho <- rho_plot(singles, draws, draw_range)
save_draft("rho_scale", p_rho, w = 9, h = 8)

# Level 3: union against its draws -------------------------------------------------------------------------

union_runs <- rho_scale[source == "CLINE-GO run" & set == "gea.union"]
union_runs[, method_name := factor(method_label[method], levels = rev(method_label))]
draws[, method_name := factor(method_label[method], levels = rev(method_label))]
p_l3a <- ggplot() +
  geom_vline(xintercept = climate_rho, colour = CLINEGO_THRESHOLD, linewidth = 0.4, linetype = "dashed") +
  geom_point(data = draws, aes(x = mean, y = method_name), colour = CLINEGO_COL$secondary, size = 1.4,
             alpha = 0.6, position = position_jitter(height = 0.15, width = 0, seed = 2)) +
  geom_point(data = union_runs, aes(x = mean, y = method_name), colour = CLINEGO_CATEGORICAL[6],
             shape = 18, size = 4.5) +
  labs(x = "Mean rho over 325 plots (diamond: union set; dots: its 20 MAF-matched random draws)", y = NULL,
       tag = "a") +
  theme_clinego(base_size = 10)

garden_draws <- garden_means[grepl(draw_pattern, set) & sub(draw_pattern, "\\1", set) == "gea.union",
                             .(draw_mean = mean(rho_offset)), by = .(method, garden)]
garden_union <- garden_means[set == "gea.union", .(method, garden, union = rho_offset)]
garden_diff <- garden_union[garden_draws, on = .(method, garden)]
garden_diff[, difference := union - draw_mean]
garden_diff[, method_name := factor(method_label[method], levels = rev(method_label))]
p_l3b <- ggplot(garden_diff, aes(x = difference, y = method_name)) +
  geom_vline(xintercept = 0, colour = CLINEGO_THRESHOLD, linewidth = 0.4, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = CLINEGO_COL$secondary, fill = NA) +
  geom_point(colour = CLINEGO_CATEGORICAL[6], size = 1.3, alpha = 0.7,
             position = position_jitter(height = 0.12, width = 0, seed = 3)) +
  labs(x = "Per garden: union's mean rho minus the mean of its 20 draws (30 gardens)", y = NULL, tag = "b") +
  theme_clinego(base_size = 10)
p_level3 <- p_l3a / p_l3b
save_draft("level3", p_level3, w = 8, h = 6)

# Per garden: Gradient Forest union against climate distance -------------------------------------------

gf_garden <- garden_means[method == "gradient_forest" & set == "gea.union"][gardens, on = "garden"]
gf_garden[, flag := fifelse(garden %in% out_of_range, "outside the fitted range", "inside")]
p_garden <- ggplot(gf_garden, aes(x = rho_climate, y = rho_offset)) +
  geom_abline(slope = 1, intercept = 0, colour = CLINEGO_THRESHOLD, linewidth = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = CLINEGO_COL$muted, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = CLINEGO_COL$muted, linewidth = 0.3) +
  geom_point(aes(fill = bio_garden, shape = flag), size = 2.8, colour = "black", stroke = 0.3) +
  geom_text_repel(aes(label = garden), size = 2.6, colour = CLINEGO_COL$secondary, seed = 4,
                  min.segment.length = 0.2) +
  scale_fill_gradientn(colours = CLINEGO_MAP_GRADIENT, name = "Garden BIO1 (°C)") +
  scale_shape_manual(values = c("inside" = 21, "outside the fitted range" = 24), name = NULL) +
  labs(x = "Climate distance: garden mean rho", y = "Gradient Forest, union set: garden mean rho") +
  theme_clinego(base_size = 10) +
  theme(legend.position = "right")
save_draft("garden_gf_climate", p_garden, w = 7.5, h = 5.2)

# Survival ----------------------------------------------------------------------------------------------

survival <- fread(file.path(offset_dir, "survival", "garden_survival.tsv"))
survival_long <- melt(survival[, .(garden, bio_garden, share_3, offset_union_gradient_forest,
                                   offset_union_geometric_offset, offset_union_rda_offset, climate_distance)],
                      id.vars = c("garden", "bio_garden", "share_3"), variable.name = "predictor",
                      value.name = "x")
survival_long[, predictor := factor(predictor,
  levels = c("offset_union_gradient_forest", "offset_union_geometric_offset", "offset_union_rda_offset",
             "climate_distance"),
  labels = c("Gradient Forest, union set", "Geometric offset, union set", "RDA offset, union set",
             "Climate distance, BIO1"))]
survival_long[, flag := fifelse(garden %in% out_of_range, "outside the fitted range", "inside")]
p_survival <- ggplot(survival_long[!is.na(share_3)], aes(x = x, y = share_3)) +
  geom_point(aes(fill = bio_garden, shape = flag), size = 2.4, colour = "black", stroke = 0.3) +
  geom_text_repel(aes(label = garden), size = 2.1, colour = CLINEGO_COL$secondary, seed = 5,
                  min.segment.length = 0.2, max.overlaps = 12) +
  facet_wrap(~ predictor, scales = "free_x", nrow = 1) +
  scale_fill_gradientn(colours = CLINEGO_MAP_GRADIENT, name = "Garden BIO1 (°C)") +
  scale_shape_manual(values = c("inside" = 21, "outside the fitted range" = 24), name = NULL) +
  labs(x = "Founder-population value at the garden (p0-weighted mean over 218 founders)",
       y = "Share of fitness plots alive in year 3") +
  theme_clinego_grid(base_size = 10) +
  theme(legend.position = "right")
save_draft("survival", p_survival, w = 13, h = 4.2)

# Cost, compact ---------------------------------------------------------------------------------------

time_memory <- fread(here("ANALYSIS", "Arabidopsis1135", "time_memory", "time_memory.tsv"))
cost <- time_memory[project %in% c("Arabidopsis1135it3", "Arabidopsis1135stress")]
cost[, project_label := sprintf("%s: %s samples x %s SNPs",
                                fifelse(project == "Arabidopsis1135it3", "iteration 3", "stress test"),
                                formatC(n_samples, format = "d", big.mark = ","),
                                formatC(n_snps, format = "d", big.mark = ","))]
mode_levels <- rev(unique(cost[project == "Arabidopsis1135it3", label]))
cost[, y := match(label, mode_levels) + fifelse(project == "Arabidopsis1135it3", 0.19, -0.19)]
cost[, hours := wall_seconds / 3600]
cost[, time_text := fifelse(wall_seconds < 3600, sprintf("%.0f min", wall_seconds / 60), sprintf("%.1f h", hours))]
cost[status != "completed", time_text := sprintf("%s, %s at %s", time_text, status, failed_rule)]
cost[, memory_text := fifelse(is.na(peak_anon_gib), "not sampled", sprintf("%.0f GiB", peak_anon_gib))]
cost_colors <- setNames(CLINEGO_CATEGORICAL[c(4, 5)], unique(cost$project_label))
mode_axis <- scale_y_continuous(breaks = seq_along(mode_levels), labels = mode_levels, expand = expansion(add = 0.5))
p_time <- ggplot(cost) +
  geom_rect(aes(xmin = 0, xmax = hours, ymin = y - 0.17, ymax = y + 0.17, fill = project_label)) +
  geom_text(aes(x = hours, y = y, label = time_text), hjust = -0.08, size = 2.5, colour = CLINEGO_COL$secondary) +
  mode_axis +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  scale_fill_manual(values = cost_colors, name = NULL) +
  labs(x = "Wall time, hours (120 cores)", y = NULL, tag = "a") +
  theme_clinego(base_size = 10)
p_memory <- ggplot(cost) +
  geom_rect(aes(xmin = 0, xmax = fcoalesce(peak_anon_gib, 0), ymin = y - 0.17, ymax = y + 0.17,
                fill = project_label)) +
  geom_text(aes(x = fcoalesce(peak_anon_gib, 0), y = y, label = memory_text), hjust = -0.08, size = 2.5,
            colour = CLINEGO_COL$secondary) +
  mode_axis +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = cost_colors, name = NULL) +
  labs(x = "Peak anonymous memory, GiB (sampled every 15 s)", y = NULL, tag = "b") +
  theme_clinego(base_size = 10) +
  theme(axis.text.y = element_blank())
p_cost <- (p_time | p_memory) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_draft("cost_compact", p_cost, w = 11, h = 5.5)

# Design map: accessions by cluster, founders and gardens -------------------------------------------

k_best <- as.integer(yaml::read_yaml(here(paste0("config_", project, ".yaml")))$sNMF$k_best)
assignments <- fread(here("ANALYSIS", project, "k_evidence", "assignments.tsv"), colClasses = c(sample = "character"))
clusters <- fread(here("ANALYSIS", project, "k_evidence", "clusters.tsv"))[K == k_best][order(-n_main)]
founders <- fread(here("data", "arabidopsis1135", "founder_sites.tsv"), select = c("sample", "in_metadata"))
garden_xy <- fread(here("data", "arabidopsis1135", "garden_sites.tsv"), select = c("garden", "longitude", "latitude"))
points <- assignments[, .(sample, longitude, latitude, cluster = get(paste0("K", k_best)))]
points[, cluster := factor(cluster, levels = clusters$cluster)]
points[, founder := sample %in% founders[in_metadata == TRUE, sample]]
xlim <- c(-125, 145)
ylim <- c(10, 70)
bio1 <- rast(here("data", "wc2.1_30s", "wc2.1_30s_bio_1.tif"))
template <- rast(ext(xlim[1], xlim[2], ylim[1], ylim[2]), resolution = 0.5, crs = crs(bio1))
land <- as.data.frame(resample(bio1, template, method = "near"), xy = TRUE, na.rm = TRUE)
p_map <- ggplot() +
  geom_raster(data = land, aes(x, y), fill = "grey90") +
  geom_point(data = points[founder == FALSE], aes(longitude, latitude, colour = cluster), size = 0.8) +
  geom_point(data = points[founder == TRUE], aes(longitude, latitude, fill = cluster), shape = 21,
             colour = "black", stroke = 0.35, size = 1.7) +
  geom_point(data = garden_xy, aes(longitude, latitude), shape = 24, fill = "white", colour = "black",
             size = 2.2, stroke = 0.6) +
  geom_text_repel(data = garden_xy, aes(longitude, latitude, label = garden), size = 2.2, seed = 6,
                  colour = CLINEGO_COL$fg, min.segment.length = 0.1, max.overlaps = 40) +
  scale_colour_manual(values = setNames(clinego_cluster_palette(k_best), clusters$cluster), name = NULL) +
  scale_fill_manual(values = setNames(clinego_cluster_palette(k_best), clusters$cluster), guide = "none") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.4))) +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("K = %d; outlined: GrENE-net founders (%d); triangles: the 30 gardens", k_best,
                          points[founder == TRUE, .N])) +
  theme_clinego(base_size = 10)
save_draft("design_map", p_map, w = 13, h = 4.6)

# Main figure, assembled draft -----------------------------------------------------------------------

raster_panel <- function(path, tag) {
  if (!file.exists(path)) {
    return(wrap_elements(full = clinego_empty_plot(paste("missing:", basename(path)))) + labs(tag = tag))
  }
  wrap_elements(full = grid::rasterGrob(png::readPNG(path), interpolate = TRUE)) + labs(tag = tag)
}
k_map <- here("ANALYSIS", project, "report", "figures", "k_map_it3.png")
miami <- here("ANALYSIS", project, "figures", "miami", "miami_combined_K11.png")
main <- (p_map + labs(tag = "a") | raster_panel(miami, "b")) /
  (p_rho + labs(tag = "c") | ((p_l3a + labs(tag = "d")) / (p_l3b + labs(tag = "e")))) /
  (p_garden + labs(tag = "f") | p_survival + labs(tag = "g")) +
  plot_layout(heights = c(0.55, 1.1, 0.75))
ggsave(file.path(draft_dir, "main_assembled.png"), main, width = 24, height = 20, dpi = 120)
message("wrote main_assembled")

# Final figures, after the user's review of 2026-11-15 ------------------------------------------------
# The idea to carry: offsets from public data only predict garden fitness at a level close to GrENE-net's own
# offset and follow it garden by garden. No random draws and no climate distance here.

final_dir <- here("ANALYSIS", project, "section", "figures")
save_final <- function(stem, plot, w, h) {
  clinego_save_both(file.path(final_dir, stem), plot, w = w, h = h, dpi = 300)
  message("wrote final ", stem)
}
experiment_short <- c("Trained on GrENE-net's evolved populations", "No GrENE-net data")
level_label <- c(
  grene.replicate_prediction = "Other plots of the same garden",
  grene.genomic_offset_gaussian = "Local adaptation model: each accession's fitness curve",
  grene.genomic_offset = "Genomic offset",
  ours.gradient_forest.gea.union = "Gradient Forest",
  ours.geometric_offset.gea.union = "Geometric offset",
  ours.rda_offset.gea.union = "RDA offset")

# a: the design map, with Europe enlarged beside it.
europe <- list(x = c(-11, 31), y = c(35.5, 62))
garden_xy[, in_europe := longitude >= europe$x[1] & longitude <= europe$x[2] &
                         latitude >= europe$y[1] & latitude <= europe$y[2]]
land_eu <- as.data.frame(resample(bio1, rast(ext(europe$x[1], europe$x[2], europe$y[1], europe$y[2]),
                                             resolution = 0.1, crs = crs(bio1)), method = "near"),
                         xy = TRUE, na.rm = TRUE)
cluster_colours <- setNames(clinego_cluster_palette(k_best), clusters$cluster)
map_layers <- function(land_df, label_gardens, label_size) {
  list(
    geom_raster(data = land_df, aes(x, y), fill = "grey90"),
    geom_point(data = points[founder == FALSE], aes(longitude, latitude, colour = cluster), size = 0.6),
    geom_point(data = points[founder == TRUE], aes(longitude, latitude, fill = cluster), shape = 21,
               colour = "black", stroke = 0.3, size = 1.4),
    geom_point(data = garden_xy, aes(longitude, latitude), shape = 24, fill = "white", colour = "black",
               size = 2, stroke = 0.55),
    geom_text_repel(data = label_gardens, aes(longitude, latitude, label = garden), size = label_size, seed = 6,
                    colour = CLINEGO_COL$fg, min.segment.length = 0.1, max.overlaps = 50, box.padding = 0.3),
    scale_fill_manual(values = cluster_colours, guide = "none")
  )
}
p_world <- ggplot() +
  map_layers(land, garden_xy[in_europe == FALSE], 2.3) +
  annotate("rect", xmin = europe$x[1], xmax = europe$x[2], ymin = europe$y[1], ymax = europe$y[2],
           fill = NA, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
  scale_colour_manual(values = cluster_colours, name = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.2))) +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(x = NULL, y = NULL, tag = "a",
       subtitle = sprintf("sNMF K = %d; outlined: GrENE-net founders (%d); triangles: the 30 gardens", k_best,
                          points[founder == TRUE, .N])) +
  theme_clinego(base_size = 9)
p_europe <- ggplot() +
  map_layers(land_eu, garden_xy[in_europe == TRUE], 2.5) +
  scale_colour_manual(values = cluster_colours, guide = "none") +
  coord_quickmap(xlim = europe$x, ylim = europe$y, expand = FALSE) +
  labs(x = NULL, y = NULL, subtitle = "Europe, enlarged") +
  theme_clinego(base_size = 9) +
  theme(panel.border = element_rect(colour = CLINEGO_THRESHOLD, fill = NA, linewidth = 0.6))
row_map <- (p_world | p_europe) + plot_layout(widths = c(2.3, 1))

# b: the Miami on the pipeline's own background cloud (the PNG and overlay coordinates the app uses), with the
# significant calls on top and the classic genes the scans found marked by dashed lines.
miami_dir <- here("ANALYSIS", project, "figures", "miami_rda_second")
coords <- jsonlite::fromJSON(file.path(miami_dir, sprintf("miami_combined_K%d_coords.json", k_best)))
background <- png::readPNG(file.path(miami_dir, sprintf("miami_combined_K%d_background.png", k_best)))
chr_offsets <- unlist(coords$chr_offsets)
results_dir <- here(paste0(project, "_results"))
read_sig <- function(module, method) {
  path <- file.path(results_dir, module, "tables", "methods", method,
                    sprintf("%s_pvalues_K%d_sig_snps_bonf_0.05.tsv", method, k_best))
  fread(path, select = c("chr", "pos", "pvalue", "method", "trait"), colClasses = c(chr = "character"))[
    , module := module][]
}
sig <- rbind(read_sig("GEA", "EMMAX"), read_sig("GEA", "LFMM"), read_sig("GEA", "RDA"),
             read_sig("GWAS", "EMMAX"), read_sig("GWAS", "BLINK"))
sig[, x := pos + chr_offsets[chr]]
sig[, y := fifelse(module == "GEA", -log10(pvalue), log10(pvalue))]
if (anyNA(sig$x)) stop("a significant SNP lies on a chromosome the Miami coordinates do not list")
traits <- c(coords$gea_traits, coords$gwas_traits)
method_shapes <- c(EMMAX = 16, LFMM = 17, RDA = 15, BLINK = 18)   # get_method_shapes()'s order (manhattan_utils.R)
sig[, trait := factor(trait, levels = traits)]
sig[, method := factor(method, levels = names(method_shapes))]
found_genes <- c("FT", "FRI", "FLC", "DOG1")
classic <- fread(here("data", "arabidopsis1135", "lists", "classic_genes.tsv"), colClasses = c(chr = "character"))
genes <- classic[label %in% found_genes][, x := (start + end) / 2 + chr_offsets[chr]]
if (nrow(genes) != length(found_genes)) stop("a classic gene is missing from classic_genes.tsv")
chr_mid <- chr_offsets + unlist(coords$chr_lengths)[names(chr_offsets)] / 2
p_miami <- ggplot() +
  annotation_raster(background, xmin = coords$x_range[1], xmax = coords$x_range[2],
                    ymin = coords$y_range[1], ymax = coords$y_range[2], interpolate = TRUE) +
  geom_hline(yintercept = c(coords$gea_threshold_y, coords$gwas_threshold_y), colour = CLINEGO_SIGNIFICANCE,
             linetype = "dashed", linewidth = 0.3) +
  geom_vline(data = genes, aes(xintercept = x), linetype = "22", colour = "grey15", linewidth = 0.4) +
  geom_point(data = sig, aes(x, y, colour = trait, shape = method), size = 1.4) +
  geom_label(data = genes, aes(x = x, y = coords$y_range[2] * 0.9, label = label), fontface = "italic",
             size = 3, label.size = 0, fill = "white", alpha = 0.9, label.padding = unit(0.12, "lines")) +
  annotate("text", x = coords$x_range[1] + diff(coords$x_range) * 0.008, y = coords$y_range[2] * 0.98,
           label = "GEA", hjust = 0, vjust = 1, fontface = "bold", size = 3, colour = CLINEGO_COL$fg) +
  annotate("text", x = coords$x_range[1] + diff(coords$x_range) * 0.008, y = coords$y_range[1] * 0.98,
           label = "GWAS", hjust = 0, vjust = 0, fontface = "bold", size = 3, colour = CLINEGO_COL$fg) +
  scale_x_continuous(breaks = chr_mid, labels = names(chr_offsets), limits = coords$x_range, expand = c(0, 0)) +
  scale_y_continuous(limits = coords$y_range, labels = function(v) abs(v), expand = c(0, 0)) +
  scale_colour_manual(values = clinego_trait_colors(traits), limits = traits, name = "Predictor or trait") +
  scale_shape_manual(values = method_shapes, limits = names(method_shapes), name = "Method") +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 1.8), order = 2),
         shape = guide_legend(nrow = 2, override.aes = list(size = 1.8), order = 1)) +
  labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +
  theme_clinego(base_size = 9) +
  # Wider than the background's own 4200 x 2400 px (aspect 0.571), so the panel spans the joint figure; the
  # cloud is stretched horizontally by about 1.36, which the dense background hides.
  theme(aspect.ratio = 0.42, legend.position = "bottom", legend.box = "horizontal",
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6.5))

# c: how well each offset ranks founders by first-year fitness, mean over the plots.
lollipop <- function(rows, base_size = 9) {
  ggplot(rows, aes(y = row_label)) +
    geom_vline(xintercept = 0, colour = CLINEGO_COL$muted, linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = mean, yend = row_label, colour = experiment), linewidth = 0.9) +
    geom_point(aes(x = mean, colour = experiment), size = 2.6) +
    geom_text(aes(x = mean, label = sprintf("%.3f", mean)), hjust = -0.35, size = 2.9, colour = CLINEGO_COL$fg) +
    scale_colour_manual(values = setNames(CLINEGO_CATEGORICAL[c(1, 5)], experiment_short), name = NULL,
                        drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Mean Spearman's rho with first-year fitness, 325 plots", y = NULL) +
    theme_clinego_grid(base_size = base_size) +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
}
level_rows <- function(predictors, run_label) {
  rows <- rho_scale[predictor %in% predictors]
  rows[, group := fifelse(startsWith(source, "GrENE-net"), "GrENE-net",
                          fifelse(predictor %in% names(level_label), "CLINE-GO",
                                  paste("CLINE-GO\n", method_label[method])))]
  rows[, row_label := run_label(predictor, set)]
  rows[, experiment := factor(experiment_short[fifelse(uses_experiment_data == "yes", 1L, 2L)],
                              levels = experiment_short)]
  rows[, row_label := factor(row_label, levels = rev(unique(row_label)))]
  rows
}
# GrENE-net's local adaptation model (0.415) is left out of panel c: it learns each accession's own fitness
# curve from the other gardens, so the text gives it as the study's upper reference (user, 2026-11-15).
main_predictors <- c("grene.genomic_offset", "ours.gradient_forest.gea.union",
                     "ours.geometric_offset.gea.union", "ours.rda_offset.gea.union")
c_rows <- level_rows(main_predictors, function(predictor, set) level_label[predictor])
c_rows <- c_rows[match(main_predictors, predictor)]
c_rows[, row_label := factor(as.character(row_label), levels = rev(level_label[main_predictors]))]
c_rows[, group := factor(group, levels = c("GrENE-net", "CLINE-GO"))]
p_level <- lollipop(c_rows)

# d: per garden, Gradient Forest's accuracy against GrENE-net's genomic offset accuracy (their per-plot results,
# averaged per garden: the only per-garden source of their offset), with Spearman's rho and p over the gardens
# from garden_agreement.tsv (rho_scale.R).
agreement_table <- fread(file.path(offset_dir, "garden_agreement.tsv"))[source == "CLINE-GO run"]
format_p <- function(p) ifelse(p < 0.001, "italic(p) < 0.001",
                               sprintf("italic(p) == '%s'", ifelse(p < 0.1, sprintf("%.3f", p), sprintf("%.2f", p))))
agreement_label <- function(rows) {
  rows[, label := sprintf("rho == '%.2f' * ',' ~~ %s", gardens_spearman, format_p(gardens_p))]
}
agreement_rows <- garden_means[set == "gea.union" & !grepl(draw_pattern, set)][gardens, on = "garden"]
agree_scatter <- function(rows, base_size = 9) {
  ggplot(rows, aes(s12_offset, rho_offset)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = CLINEGO_THRESHOLD, linewidth = 0.4,
                linetype = "dashed") +
    geom_point(aes(fill = bio_garden), shape = 21, colour = "black", stroke = 0.3, size = 2.3) +
    geom_text_repel(aes(label = garden), size = 2.1, colour = CLINEGO_COL$secondary, seed = 7,
                    min.segment.length = 0.2, max.overlaps = 30) +
    scale_fill_gradientn(colours = CLINEGO_MAP_GRADIENT, name = "Garden\nBIO1 (°C)") +
    labs(x = "GrENE-net genomic offset: mean rho per garden", y = "CLINE-GO: mean rho per garden") +
    theme_clinego(base_size = base_size) +
    theme(legend.position = "right")
}
gf_label <- agreement_label(agreement_table[method == "gradient_forest" & set == "gea.union"])
p_agree <- agree_scatter(agreement_rows[method == "gradient_forest"]) +
  geom_text(data = gf_label, aes(x = -Inf, y = Inf, label = label), parse = TRUE, hjust = -0.08, vjust = 1.4,
            size = 3.2, colour = CLINEGO_COL$fg, inherit.aes = FALSE) +
  labs(y = "CLINE-GO Gradient Forest:\nmean rho per garden")

row_offsets <- (p_level + labs(tag = "c") | p_agree + labs(tag = "d")) + plot_layout(widths = c(1.25, 1))
main_final <- row_map / (p_miami + labs(tag = "b")) / row_offsets + plot_layout(heights = c(0.8, 1.35, 0.95))
save_final("main_arabidopsis", main_final, w = 10, h = 11)
save_final("main_panel_a", row_map, w = 10, h = 3.8)
save_final("main_panel_b", p_miami, w = 10, h = 5.2)
save_final("main_panel_c", p_level, w = 6, h = 3.2)
save_final("main_panel_d", p_agree, w = 5.2, h = 3.6)

# Supplement a: every set per method beside GrENE-net's predictors.
all_predictors <- c("grene.replicate_prediction", "grene.genomic_offset_gaussian", "grene.genomic_offset",
                    unlist(lapply(names(method_label), function(m) paste0("ours.", m, ".", set_order))))
s_rows <- level_rows(all_predictors, function(predictor, set) {
  fifelse(predictor %in% names(level_label) & startsWith(predictor, "grene."), level_label[predictor],
          set_label[set])
})
s_rows <- s_rows[match(all_predictors, predictor)]
s_rows[, group := factor(fifelse(startsWith(predictor, "grene."), "GrENE-net",
                                 paste0("CLINE-GO\n", method_label[method])),
                         levels = c("GrENE-net", paste0("CLINE-GO\n", method_label)))]
s_rows[, row_label := factor(paste(group, row_label, sep = "|"), levels = rev(paste(group, row_label, sep = "|")))]
p_level_all <- lollipop(s_rows) + scale_y_discrete(labels = function(v) sub(".*\\|", "", v))
save_final("s_rho_scale_all", p_level_all, w = 8, h = 7)

# Supplement b: agreement with GrENE-net's offset for each method, and for every set.
# The agreement for every set is a table on the review page (garden_agreement.tsv), not a panel.
agreement_rows[, method_name := factor(method_label[method], levels = method_label)]
method_labels <- agreement_label(agreement_table[set == "gea.union"])
method_labels[, method_name := factor(method_label[method], levels = method_label)]
p_agree_methods <- agree_scatter(agreement_rows) + facet_wrap(~ method_name, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  geom_text(data = method_labels, aes(x = -Inf, y = Inf, label = label), parse = TRUE, hjust = -0.08,
            vjust = 1.4, size = 3, colour = CLINEGO_COL$fg, inherit.aes = FALSE) +
  theme_clinego_grid(base_size = 9) + theme(legend.position = "right")
save_final("s_agreement", p_agree_methods, w = 11, h = 4)

# The garden survival test (survival_eval.R) has no final figure: the user dropped it on 2026-11-15. The text
# says that CLINE-GO's offset agrees with GrENE-net's, whose link to survival their Fig. S66 shows.

# Supplement c, d: LD decay (the pipeline's plot) and cost.
save_final("s_cost", p_cost, w = 11, h = 5.5)
ld_stem <- here(paste0(project, "_results"), "Structure", "plots", "ld_decay", "ld_decay_genome_wide")
copied <- c(file.copy(paste0(ld_stem, ".png"), file.path(final_dir, "s_ld_decay.png"), overwrite = TRUE),
            file.copy(paste0(ld_stem, ".svg"), file.path(final_dir, "s_ld_decay.svg"), overwrite = TRUE))
if (!all(copied)) stop("could not copy the LD decay plot")
message("copied s_ld_decay")
