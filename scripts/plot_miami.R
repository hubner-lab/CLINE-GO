#!/usr/bin/env Rscript
# Static Miami plot combining GEA (top) and GWAS (bottom)
# GEA data plotted with positive -log10(p), GWAS with negative -log10(p)
# Color = trait (Okabe-Ito), Shape = method, scattermore mandatory for background

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
GEA_FILES_STR  = args[1]       # Comma-sep "METHOD:ADJUST:FILEPATH"
GWAS_FILES_STR = args[2]       # Comma-sep "METHOD:ADJUST:FILEPATH"
PREDICTORS_GEA = args[3]       # Comma-sep trait names (bio_1,bio_12)
PREDICTORS_GWAS = args[4]      # Comma-sep trait names (height,flowering_time)
Kbest          = args[5] %>% as.numeric
PLOT_DIR       = args[6]
REGIME         = if (length(args) >= 7) tolower(args[7]) else "snp"  # "snp" or "wza"
################################

is_wza   <- REGIME == "wza"
name_pfx <- if (is_wza) "wza_" else ""

message('INFO: Static Miami plot (GEA top, GWAS bottom)')
message(paste0('INFO: K = ', Kbest, ' [regime=', REGIME, ']'))

################################ Parse input files

gea_info  <- parse_assoc_files_str(GEA_FILES_STR)
gwas_info <- parse_assoc_files_str(GWAS_FILES_STR)

gea_traits  <- str_split(PREDICTORS_GEA,  ',')[[1]] %>% str_trim()
gwas_traits <- str_split(PREDICTORS_GWAS, ',')[[1]] %>% str_trim()

message(paste0('INFO: GEA methods: ',  paste(names(gea_info),  collapse = ', ')))
message(paste0('INFO: GWAS methods: ', paste(names(gwas_info), collapse = ', ')))
message(paste0('INFO: GEA traits: ',   paste(gea_traits,  collapse = ', ')))
message(paste0('INFO: GWAS traits: ',  paste(gwas_traits, collapse = ', ')))

################################ Load data

gea_result  <- load_assoc_data(gea_info,  gea_traits,  "GEA")
gwas_result <- load_assoc_data(gwas_info, gwas_traits, "GWAS")

gea_data  <- gea_result$data
gwas_data <- gwas_result$data

# Calculate thresholds (min per panel)
gea_threshold  <- min(unlist(gea_result$thresholds))
gwas_threshold <- min(unlist(gwas_result$thresholds))
message(paste0('INFO: GEA threshold: -log10(p) = ',  round(gea_threshold,  2)))
message(paste0('INFO: GWAS threshold: -log10(p) = ', round(gwas_threshold, 2)))

################################ Prepare positions

# Use GEA data for chromosome coordinate system
reference_data <- gea_data %>%
    dplyr::filter(trait == gea_traits[1], method == names(gea_info)[1]) %>%
    dplyr::select(SNPID, chr, pos, pvalue)

manhattan_prep <- prepare_manhattan_data(reference_data, "chr", "pos", "pvalue")
chr_info <- manhattan_prep$chr_info

# Add cumulative positions to both panels
add_cum_pos <- function(df) {
    df %>%
        left_join(chr_info %>%
                  dplyr::mutate(chr = as.character(chr_f)) %>%
                  dplyr::select(chr, tot),
                  by = "chr") %>%
        dplyr::mutate(pos_cum = pos + tot)
}

gea_data  <- add_cum_pos(gea_data)
gwas_data <- add_cum_pos(gwas_data)

# For Miami: GEA positive, GWAS negative
gea_data  <- gea_data  %>% dplyr::mutate(miami_y =  log10p)
gwas_data <- gwas_data %>% dplyr::mutate(miami_y = -log10p)

################################ Palettes

# All traits (union of GEA + GWAS) for consistent colors
all_traits  <- unique(c(gea_traits, gwas_traits))
trait_colors <- get_trait_colors(all_traits)

# All methods (union) for consistent shapes
all_methods  <- unique(c(names(gea_info), names(gwas_info)))
method_shapes <- get_method_shapes(all_methods)

################################ Y-axis limits

y_max_gea  <- max(gea_data$log10p,  na.rm = TRUE) * 1.1
y_max_gwas <- max(gwas_data$log10p, na.rm = TRUE) * 1.1
y_limit    <- max(y_max_gea, y_max_gwas)

################################ Build simple Miami plot

message('INFO: Generating simple Miami plot')

# Background SNPs: ALL SNPs, ALL methods (threshold-independent), full cloud per panel.
# Sig markers overlay on top; threshold lines drawn live in Shiny.
# Every SNP is baked in — background = complete static Manhattan per panel.
gea_bg  <- gea_data
gwas_bg <- gwas_data

# Sig SNPs
gea_sig  <- gea_data  %>% dplyr::filter(is_significant)
gwas_sig <- gwas_data %>% dplyr::filter(is_significant)

n_gea_sig  <- length(unique(gea_sig$SNPID))
n_gwas_sig <- length(unique(gwas_sig$SNPID))
message(paste0('INFO: GEA significant: ',  n_gea_sig,  ' unique SNPs, ', nrow(gea_sig),  ' detections'))
message(paste0('INFO: GWAS significant: ', n_gwas_sig, ' unique SNPs, ', nrow(gwas_sig), ' detections'))

# Simple Miami — per-trait background scattermore (chr gaps for separation)
p_miami <- ggplot()

# GEA background per trait
for (t in gea_traits) {
    df_t <- gea_bg %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_miami <- p_miami +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = miami_y, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 2400), interpolate = FALSE)
}

# GWAS background per trait
for (t in gwas_traits) {
    df_t <- gwas_bg %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_miami <- p_miami +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = miami_y, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 2400), interpolate = FALSE)
}

p_miami <- p_miami +
    scale_color_identity() +
    # GEA sig SNPs
    ggnewscale::new_scale_color() +
    geom_point(data = gea_sig,
               aes(x = pos_cum, y = miami_y, color = trait, shape = method),
               alpha = 0.9, size = 2.5, stroke = 0.5) +
    # GWAS sig SNPs
    geom_point(data = gwas_sig,
               aes(x = pos_cum, y = miami_y, color = trait, shape = method),
               alpha = 0.9, size = 2.5, stroke = 0.5) +
    scale_color_manual(values = trait_colors, name = "Trait") +
    scale_shape_manual(values = method_shapes, name = "Method") +
    # Threshold lines
    geom_hline(yintercept =  gea_threshold,  linetype = "dashed",
               color = "red", linewidth = 0.4, alpha = 0.5) +
    geom_hline(yintercept = -gwas_threshold, linetype = "dashed",
               color = "red", linewidth = 0.4, alpha = 0.5) +
    # Bold zero line separating GEA and GWAS
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
    # Axes
    scale_x_continuous(
        labels = chr_info$chr_f,
        breaks = chr_info$center,
        expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
        limits = c(-y_limit, y_limit),
        labels = function(x) abs(x),
        expand = expansion(mult = c(0.02, 0.02))
    ) +
    # Labels
    labs(
        title = "Miami Plot: GEA (top) vs GWAS (bottom)",
        subtitle = paste0("K = ", Kbest, " | GEA: ", n_gea_sig, " sig SNPs | GWAS: ", n_gwas_sig, " sig SNPs"),
        x = "Chromosome",
        y = expression(-log[10](p-value))
    ) +
    # Annotate panel labels
    annotate("text", x = -Inf, y = y_limit * 0.95, label = "GEA",
             hjust = -0.3, fontface = "bold", size = 4, color = "grey30") +
    annotate("text", x = -Inf, y = -y_limit * 0.95, label = "GWAS",
             hjust = -0.3, fontface = "bold", size = 4, color = "grey30") +
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

# Save simple Miami
simple_base <- paste0("miami_", name_pfx, "combined_K", Kbest)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".png")), p_miami,
       width = 14, height = 8, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".svg")), p_miami,
       width = 14, height = 8)
message(paste0('INFO: Saved simple Miami plot: ', simple_base))

################################ Background PNG + Coords JSON (for Shiny plotly overlay)

message('INFO: Generating Miami background for Shiny overlay')

all_pos_cum <- c(gea_data$pos_cum, gwas_data$pos_cum)
x_min_data  <- min(all_pos_cum)
x_max_data  <- max(all_pos_cum)
x_span      <- x_max_data - x_min_data
bg_x_lo     <- x_min_data - 0.01 * x_span
bg_x_hi     <- x_max_data + 0.01 * x_span

# y: limits = c(-y_limit, y_limit), expand = mult c(0.02, 0.02)
bg_y_lo <- -y_limit * 1.04
bg_y_hi <-  y_limit * 1.04

p_bg_miami <- ggplot()

for (t in gea_traits) {
    df_t <- gea_bg %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_bg_miami <- p_bg_miami +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = miami_y, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 2400), interpolate = FALSE)
}

for (t in gwas_traits) {
    df_t <- gwas_bg %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_bg_miami <- p_bg_miami +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = miami_y, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 2400), interpolate = FALSE)
}

p_bg_miami <- p_bg_miami +
    scale_color_identity() +
    # Zero line baked in (structural); threshold lines drawn live in Shiny from coords JSON.
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
    scale_x_continuous(limits = c(bg_x_lo, bg_x_hi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(bg_y_lo, bg_y_hi), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

bg_miami_base <- paste0("miami_", name_pfx, "combined_K", Kbest, "_background")
ggsave(file.path(PLOT_DIR, paste0(bg_miami_base, ".png")), p_bg_miami,
       width = 14, height = 8, dpi = 300)
message(paste0('INFO: Saved Miami background: ', bg_miami_base, '.png'))

coords_miami <- list(
    chr_offsets      = setNames(as.list(chr_info$tot),     as.character(chr_info$chr_f)),
    chr_lengths      = setNames(as.list(chr_info$chr_len), as.character(chr_info$chr_f)),
    gap_fraction     = 0.02,
    x_range          = c(bg_x_lo, bg_x_hi),
    y_range          = c(bg_y_lo, bg_y_hi),
    gea_threshold_y  = gea_threshold,
    gwas_threshold_y = -gwas_threshold,
    gea_traits       = as.list(gea_traits),
    gwas_traits      = as.list(gwas_traits),
    is_miami         = TRUE,
    plot_width_px    = 4200L,
    plot_height_px   = 2400L
)
coords_miami_file <- file.path(PLOT_DIR, paste0("miami_", name_pfx, "combined_K", Kbest, "_coords.json"))
jsonlite::write_json(coords_miami, coords_miami_file, auto_unbox = TRUE, digits = 6)
message(paste0('INFO: Saved Miami coords JSON: ', basename(coords_miami_file)))

message('INFO: Miami plot complete')
