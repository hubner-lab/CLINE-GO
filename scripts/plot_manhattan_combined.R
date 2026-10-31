#!/usr/bin/env Rscript
# Combined Manhattan plot for all traits and methods
# Shows all GWAS results on a single plot with:
# - Color = trait (Okabe-Ito colorblind-safe palette)
# - Shape = method (up to 10 distinct shapes)
# - Every (SNP, trait, method) detection shown individually (no "Multi" collapse)
# - Regions marked with rectangles (from Regions_climate_combined.tsv)
# - Always uses scattermore for background SNPs

library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(qvalue)
library(scales)
library(scattermore)
library(ggnewscale)
library(jsonlite)

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/utils/manhattan_utils.R")
source("/pipeline/scripts/R/utils/io_pvalues.R")
source("/pipeline/scripts/R/utils/scattermore_threads.R")

args = commandArgs(trailingOnly=TRUE)
################################
ASSOC_FILES_STR = args[1]   # Comma-separated "METHOD:ADJUST:FILEPATH"
PREDICTORS = args[2]        # Comma-separated trait names
Kbest = args[3] %>% as.numeric
PLOT_DIR = args[4]          # Output directory
REGIME    = if (length(args) >= 5) tolower(args[5]) else "snp"  # "snp" or "wza"
################################

is_wza   <- REGIME == "wza"
name_pfx <- if (is_wza) "wza_" else ""

message('INFO: Combined Manhattan plot for all traits and methods')
message(paste0('INFO: K = ', Kbest, ' [regime=', REGIME, ']'))

# Parse ASSOC_FILES_STR into named list
assoc_info <- parse_assoc_files_str(ASSOC_FILES_STR)

if (length(assoc_info) == 0) {
    stop("No valid association files provided")
}

message(paste0('INFO: Found ', length(assoc_info), ' methods: ', paste(names(assoc_info), collapse = ', ')))

# Parse predictors
traits <- str_split(PREDICTORS, ',')[[1]] %>% str_trim()
message(paste0('INFO: Traits: ', paste(traits, collapse = ', ')))

################################ Main

# Load and threshold all association data
message('INFO: Loading association data for all methods and traits')
assoc_result    <- load_assoc_data(assoc_info, traits)
plot_data_all   <- assoc_result$data
method_thresholds <- assoc_result$thresholds

# Filter to methods+traits where threshold computation succeeded (status="ok").
# Methods returning NA (too_few_tests) are plotted but excluded from the threshold line.
ok_thresholds <- Filter(function(x) !is.na(x), method_thresholds)
if (length(ok_thresholds) == 0) {
    message('WARNING: No methods produced a valid threshold — threshold line omitted')
    min_threshold <- NA_real_
} else {
    min_threshold <- min(unlist(ok_thresholds))
    message(paste0('INFO: Minimum threshold (most stringent): -log10(p) = ', round(min_threshold, 2)))
}

# Prepare manhattan data (using first trait for chromosome info)
# ORDERING ASSUMPTION: this requires traits[1] to exist in names(assoc_info)[1]'s
# table. Holds for every planned config (a multivariate method's pseudo-trait,
# e.g. RDA's climate_multivariate, is appended AFTER real predictors by
# combine_predictors() in common.smk, and RDA is never the sole/first configured
# method) — if that ever changes, reference_data comes back empty and
# prepare_manhattan_data() breaks on it.
reference_data <- plot_data_all %>%
    dplyr::filter(trait == traits[1], method == names(assoc_info)[1]) %>%
    dplyr::select(SNPID, chr, pos, pvalue)

manhattan_prep <- prepare_manhattan_data(reference_data, "chr", "pos", "pvalue")
chr_info <- manhattan_prep$chr_info

# Add cumulative positions to all data
plot_data_all <- plot_data_all %>%
    left_join(chr_info %>%
              dplyr::mutate(chr = as.character(chr_f)) %>%
              dplyr::select(chr, tot),
              by = "chr") %>%
    dplyr::mutate(pos_cum = pos + tot)

# Get palettes
trait_colors  <- get_trait_colors(traits)
method_shapes <- get_method_shapes(names(assoc_info))

# Calculate y-axis limit — extend to include threshold line when it exceeds data cloud
y_max <- max(plot_data_all$log10p, na.rm = TRUE) * 1.1
if (!is.na(min_threshold)) y_max <- max(y_max, min_threshold * 1.1)

message('INFO: Generating combined Manhattan plot (simple version)')

# Background SNPs: ALL SNPs, ALL methods (threshold-independent), full data cloud.
# Sig markers overlay on top; toggling a marker off leaves its faint background dot.
# Every SNP (significant or not, every method) is baked in — background = complete static Manhattan.
# Drop non-finite log10p (NA WZA windows; +Inf from p==0) — scattermore aborts on them
# (geom_point in the sig overlay tolerates non-finite, so df_significant is left intact).
df_background <- plot_data_all %>% dplyr::filter(is.finite(log10p))

# Significant SNPs: every (SNP, trait, method) combination shown individually
df_significant <- plot_data_all %>%
    dplyr::filter(is_significant)

n_sig_points <- nrow(df_significant)
n_sig_unique <- length(unique(df_significant$SNPID))
message(paste0('INFO: Total significant points (SNP x trait x method): ', n_sig_points))
message(paste0('INFO: Unique significant SNPs: ', n_sig_unique))

# Simple plot — per-trait background scattermore (chr gaps for separation)
p_simple <- ggplot()

# Per-trait background SNPs (scattermore, first method only)
for (t in traits) {
    df_t <- df_background %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_simple <- p_simple +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = log10p, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 1500), interpolate = FALSE)
}

p_simple <- p_simple +
    scale_color_identity() +
    # 3. Significant SNPs - color by trait, shape by method
    ggnewscale::new_scale_color() +
    geom_point(data = df_significant,
               aes(x = pos_cum, y = log10p, color = trait, shape = method),
               alpha = 0.9, size = 2.5, stroke = 0.5) +
    scale_color_manual(values = trait_colors, name = "Trait") +
    scale_shape_manual(values = method_shapes, name = "Method") +
    scale_x_continuous(
        labels = chr_info$chr_f,
        breaks = chr_info$center,
        expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.05)),
        limits = c(0, y_max)
    ) +
    labs(
        title    = "Combined Association Analysis (All Traits & Methods)",
        subtitle = paste0("K = ", Kbest, " | ", n_sig_unique, " significant SNPs, ",
                         n_sig_points, " detections across ",
                         length(names(assoc_info)), " methods"),
        x = "Chromosome",
        y = expression(-log[10](p-value))
    ) +
    theme_manhattan() +
    theme(
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.5, "cm")
    ) +
    guides(
        color = guide_legend(order = 1),
        shape = guide_legend(order = 2)
    )

# Add threshold line to static export only (Shiny draws it dynamically from bonferroni_y)
if (!is.na(min_threshold)) {
    p_simple <- p_simple +
        geom_hline(yintercept = min_threshold, linetype = "dashed",
                   color = CLINEGO_SIGNIFICANCE, linewidth = 0.5, alpha = 0.5)   # rule 9
}

# Save simple plot
simple_base <- paste0("manhattan_", name_pfx, "combined_K", Kbest)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".png")), p_simple,
       width = 14, height = 5, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".svg")), p_simple,
       width = 14, height = 5)
message(paste0('INFO: Saved combined simple plot: ', simple_base))

################################ Background PNG + Coords JSON (for Shiny plotly overlay)

message('INFO: Generating combined background Manhattan for Shiny overlay')

x_min_data <- min(plot_data_all$pos_cum)
x_max_data <- max(plot_data_all$pos_cum)
x_span     <- x_max_data - x_min_data
bg_x_lo    <- x_min_data - 0.01 * x_span
bg_x_hi    <- x_max_data + 0.01 * x_span
bg_y_lo    <- 0
bg_y_hi    <- y_max * 1.05

p_bg_comb <- ggplot()
for (t in traits) {
    df_t_bg <- df_background %>% dplyr::filter(trait == t)
    df_t_bg$point_color <- unname(trait_colors[t])
    p_bg_comb <- p_bg_comb +
        geom_scattermore(data = df_t_bg,
                         aes(x = pos_cum, y = log10p, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 1500), interpolate = FALSE)
}

p_bg_comb <- p_bg_comb +
    scale_color_identity() +
    scale_x_continuous(limits = c(bg_x_lo, bg_x_hi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(bg_y_lo, bg_y_hi), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

bg_comb_base <- paste0("manhattan_", name_pfx, "combined_K", Kbest, "_background")
ggsave(file.path(PLOT_DIR, paste0(bg_comb_base, ".png")), p_bg_comb,
       width = 14, height = 5, dpi = 300)
message(paste0('INFO: Saved combined background Manhattan: ', bg_comb_base, '.png'))

coords_list_comb <- list(
    chr_offsets    = setNames(as.list(chr_info$tot),     as.character(chr_info$chr_f)),
    chr_lengths    = setNames(as.list(chr_info$chr_len), as.character(chr_info$chr_f)),
    gap_fraction   = 0.02,
    x_range        = c(bg_x_lo, bg_x_hi),
    y_range        = c(bg_y_lo, bg_y_hi),
    bonferroni_y   = if (!is.na(min_threshold)) min_threshold else jsonlite::unbox(NA),
    plot_width_px  = 4200L,
    plot_height_px = 1500L
)
coords_file_comb <- file.path(PLOT_DIR, paste0("manhattan_", name_pfx, "combined_K", Kbest, "_coords.json"))
jsonlite::write_json(coords_list_comb, coords_file_comb, auto_unbox = TRUE, digits = 6)
message(paste0('INFO: Saved combined coords JSON: ', basename(coords_file_comb)))

################################ QQ Plot (Combined)

message('INFO: Generating combined QQ plot')

# For each (trait, method) combo: sort p-values, compute expected uniform quantiles
qq_data <- plot_data_all %>%
    group_by(trait, method) %>%
    arrange(pvalue) %>%
    mutate(
        expected = -log10(ppoints(n())),
        observed = -log10(pvalue)
    ) %>%
    ungroup()

# 95% confidence interval band (under uniform null)
n_snps_ci <- max(table(paste(qq_data$trait, qq_data$method)))
ci_data <- data.frame(
    expected = -log10(ppoints(n_snps_ci)),
    ci_lower = -log10(qbeta(0.975, 1:n_snps_ci, n_snps_ci:1)),
    ci_upper = -log10(qbeta(0.025, 1:n_snps_ci, n_snps_ci:1))
)

p_qq <- ggplot() +
    geom_ribbon(data = ci_data,
                aes(x = expected, ymin = ci_lower, ymax = ci_upper),
                fill = "grey80", alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    # Non-significant points
    geom_point(data = qq_data %>% filter(!is_significant),
               aes(x = expected, y = observed, color = trait, shape = method),
               alpha = 0.5, size = 1.5) +
    # Significant points (larger)
    geom_point(data = qq_data %>% filter(is_significant),
               aes(x = expected, y = observed, color = trait, shape = method),
               size = 2.5) +
    scale_color_manual(values = trait_colors, name = "Trait") +
    scale_shape_manual(values = method_shapes, name = "Method") +
    labs(
        title = "Combined QQ Plot (All Traits & Methods)",
        subtitle = paste0("K = ", Kbest, " | ", length(names(assoc_info)), " methods, ",
                         length(traits), " traits"),
        x = expression(Expected ~ -log[10](p-value)),
        y = expression(Observed ~ -log[10](p-value))
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        plot.margin   = margin(4, 8, 4, 4)
    ) +
    guides(
        color = guide_legend(order = 1),
        shape = guide_legend(order = 2)
    )

qq_base <- paste0("qq_", name_pfx, "combined_K", Kbest)
ggsave(file.path(PLOT_DIR, paste0(qq_base, ".png")), p_qq,
       width = 7, height = 7, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(qq_base, ".svg")), p_qq,
       width = 7, height = 7)
message(paste0('INFO: Saved combined QQ plot: ', qq_base))

message('INFO: Combined Manhattan plot complete')
