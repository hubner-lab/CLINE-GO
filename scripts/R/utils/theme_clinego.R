# theme_clinego.R — shared ggplot2 theme + palette matched to the Shiny
# app's "Precision Genomics" bslib theme (scripts/clinego.app/R/app_theme.R).
# Usage: source("/pipeline/scripts/R/utils/theme_clinego.R")
# Requires: ggplot2 loaded by the sourcing script.
#
# Rolled out to the Processing module first (plot_qc_processing.R). Intent is
# to extend this to the rest of the pipeline's plot scripts over time — see
# CLAUDE.md "Code Structure Rules". Colors only for now: no base_family is set
# (Docker installs no fonts; setting one would silently fall back to default
# sans while implying a match that isn't there — same trap as the existing
# dead base_family='Helvetica' settings elsewhere).

# Core palette — mirrors app_theme.R's bs_theme() colors exactly.
CLINEGO_COL <- list(
    primary   = "#1B7A6E",  # deep teal-green
    secondary = "#4A5568",  # slate
    success   = "#38A169",
    info      = "#3182CE",
    warning   = "#D69E2E",  # amber
    danger    = "#E53E3E",  # red
    fg        = "#2D3748",
    muted     = "#718096"
)

# Semantic aliases — use these in plot code, not raw hex, so the meaning of a
# color stays consistent across every processing plot (and future ones).
# Data-point/fill colors below are the last 3 of wesanderson::wes_palette
# ("Rushmore1") = c("#E1BD6D","#EABE94","#0B775E","#35274A","#F2300F") — user
# preference, chosen over the app-matched teal/slate/amber. Note the middle
# value is a dark plum/purple, not literally blue.
CLINEGO_RETAINED  <- "#0B775E"  # green — sample/SNP passed a filter
CLINEGO_REMOVED   <- "#F2300F"  # red   — sample/SNP dropped by a filter
CLINEGO_THRESHOLD <- "#35274A"  # dark plum/purple — threshold/cutoff lines
CLINEGO_NEUTRAL   <- CLINEGO_RETAINED  # plain/no-flag data points — same green as
                                   # retained by request (no more grey dots)

# Okabe-Ito colorblind-safe categorical palette (also used by
# scripts/R/utils/manhattan_utils.R — kept identical here to converge the
# pipeline's two independent categorical palettes onto one going forward).
CLINEGO_CATEGORICAL <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                        "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Extended categorical palette for many-level grouping (K-clusters, sites) where
# 8 colors isn't enough — converges the 21-color `my.colors` duplicated across
# plot_structure.R / plot_pca_structure.R onto one shared palette. Colorblind-
# oriented (Wong 2011 + Paul Tol); no grey (CLINEGO_NEUTRAL already owns "no flag").
CLINEGO_CLUSTERS <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
                    "#332288", "#882255", "#44AA99", "#AA4499", "#999933",
                    "#E69F00", "#661100", "#88CCEE", "#CC6677", "#117733",
                    "#6A3D9A", "#855C75", "#D9AF6B", "#736F4C", "#526A83", "#625377")

#' n distinct colors from CLINEGO_CLUSTERS; interpolates when n exceeds the palette.
#' @noRd
clinego_cluster_palette <- function(n) {
    if (n <= length(CLINEGO_CLUSTERS)) CLINEGO_CLUSTERS[seq_len(n)]
    else grDevices::colorRampPalette(CLINEGO_CLUSTERS)(n)
}

#' Shared ggplot2 theme — "Publication Classic"
#'
#' theme_classic() base (axis lines, no gridlines) + Precision Genomics sizing.
#' No base_family: device-default sans until the Docker image ships fonts.
#' @noRd
theme_clinego <- function(base_size = 11) {
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
        plot.title      = ggplot2::element_text(size = base_size + 2, face = "bold", color = CLINEGO_COL$fg),
        plot.subtitle   = ggplot2::element_text(size = base_size - 2, color = CLINEGO_COL$muted),
        axis.title      = ggplot2::element_text(size = base_size, face = "bold", color = CLINEGO_COL$fg),
        axis.text       = ggplot2::element_text(color = CLINEGO_COL$secondary),
        axis.line       = ggplot2::element_line(color = CLINEGO_COL$secondary, linewidth = 0.4),
        axis.ticks      = ggplot2::element_line(color = CLINEGO_COL$secondary, linewidth = 0.4),
        legend.position = "bottom",
        legend.title    = ggplot2::element_text(size = base_size - 1, color = CLINEGO_COL$fg),
        legend.text     = ggplot2::element_text(size = base_size - 2, color = CLINEGO_COL$secondary)
    )
}

#' Shared ggplot2 theme for FACETED GRIDS (PreGEA ladder histogram/QQ grids
#' and similar multi-panel diagnostics). theme_clinego() alone leaves
#' facet panels borderless on a plain white background — hard to tell where
#' one rung/trait panel ends and the next begins. This adds a visible panel
#' border, a shaded+bold strip label (reads as a header, not stray text), and
#' breathing room between panels. Canvas size for these grids is already
#' computed by the caller (e.g. plot_pregea_ladder.R's 2.2*n_traits+2 width),
#' so bumping base_size from the old 9 to the default 11 here is safe.
#' @noRd
theme_clinego_grid <- function(base_size = 11) {
    theme_clinego(base_size = base_size) +
    ggplot2::theme(
        panel.border     = ggplot2::element_rect(colour = CLINEGO_COL$secondary, fill = NA, linewidth = 0.5),
        axis.line        = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "#EDF2F7", colour = CLINEGO_COL$secondary),
        strip.text       = ggplot2::element_text(face = "bold", size = base_size, color = CLINEGO_COL$fg),
        panel.spacing    = grid::unit(0.6, "lines")
    )
}

#' Categorical color scale using the Okabe-Ito palette above 8 levels of grouping.
#' @noRd
scale_color_clinego <- function(...) {
    ggplot2::scale_color_manual(values = CLINEGO_CATEGORICAL, ...)
}

#' @noRd
scale_fill_clinego <- function(...) {
    ggplot2::scale_fill_manual(values = CLINEGO_CATEGORICAL, ...)
}

#' Save a ggplot as both PNG and SVG with one call — the pair every
#' pipeline plot script writes. Hoisted out of pregea_varpart.R (the only
#' script that had it) so pregea_rda_setup.R's 5 duplicated ggsave() pairs
#' can use the same helper (CLAUDE.md: avoid redundancy aggressively).
#' @noRd
clinego_save_both <- function(path_stem, plot, w = 7, h = 5, dpi = 300) {
    ggplot2::ggsave(paste0(path_stem, ".png"), plot, width = w, height = h, dpi = dpi)
    ggplot2::ggsave(paste0(path_stem, ".svg"), plot, width = w, height = h,
                    device = svglite::svglite, bg = "transparent")
}

#' Placeholder plot for an empty/skipped analysis state — explains WHY there
#' is nothing to show instead of leaving a blank or (worse) a 0-byte file.
#' Bare file.create() stubs (the old pregea_rda_setup.R pattern) produce
#' 0-byte PNGs that render as broken images, not empty-but-informative ones.
#' @noRd
clinego_empty_plot <- function(msg) {
    ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = msg, color = CLINEGO_COL$muted, size = 4) +
        ggplot2::theme_void()
}
