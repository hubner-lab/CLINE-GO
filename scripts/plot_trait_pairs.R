#!/usr/bin/env Rscript
# Pairwise scatter/density/correlation grid for phenotypic traits (mode=traits).
#
# Traits only — never climate (climate factors are characterized by their own
# correlogram in mode=climate; mixing them here would make a k x k grid that
# nobody can read).
#
# Hand-rolled with ggplot2 + ggpubr::ggarrange rather than GGally::ggpairs:
# GGally is not in the image and adding it needs a Dockerfile change + rebuild.
# TODO (logged as a follow-up): swap to GGally::ggpairs once the image carries it.
#
# Layout per cell, standard pairs-plot convention:
#   diagonal        density of that trait
#   lower triangle  scatter + loess/lm smooth
#   upper triangle  Pearson r, text size scaled by |r|
#
# Above MAX_FACTORS traits the grid is unreadable and the render cost is
# quadratic, so a placeholder is written instead — PNG *and* SVG, both declared
# outputs, exit 0 (same contract as pregea_dbmem.R's write_skip() and
# ld_decay_analyze.R's .write_placeholder_plot()). The count itself is
# data-derived, so per CLAUDE.md rule 8 the *explanation* of the cap lives in
# the Shiny help note, not baked onto the image.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(svglite)
})

source("/pipeline/scripts/R/utils/theme_clinego.R")

args <- commandArgs(trailingOnly = TRUE)
META_FILE   <- args[1]
OUT_PNG     <- args[2]
MAX_FACTORS <- as.integer(args[3])

OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

save_both <- function(p, width, height) {
    suppressMessages(ggsave(OUT_PNG, p, width = width, height = height, dpi = 300, limitsize = FALSE))
    suppressMessages(ggsave(OUT_SVG, p, width = width, height = height,
                            device = svglite::svglite, bg = "transparent", limitsize = FALSE))
}

write_placeholder <- function(msg) {
    p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, color = "grey50") +
        theme_void()
    save_both(p, width = 10, height = 7)
    message("INFO: Wrote placeholder pairs plot: ", msg)
}

meta <- fread(META_FILE, sep = "\t", header = TRUE,
              colClasses = c("site" = "character", "sample" = "character"))

meta_cols <- c("sample", "site", "latitude", "longitude", "lat", "lon", "ID")
cand <- colnames(meta)[!tolower(colnames(meta)) %in% tolower(meta_cols)]
trait_cols <- cand[vapply(cand, function(x) is.numeric(meta[[x]]), logical(1))]

# Constant / all-NA traits carry no pairwise information and break cor()
keep <- vapply(trait_cols, function(x) {
    v <- suppressWarnings(as.numeric(meta[[x]]))
    !all(is.na(v)) && stats::sd(v, na.rm = TRUE) > 0
}, logical(1))
if (length(trait_cols) > 0 && any(!keep)) {
    message("INFO: dropped ", sum(!keep), " constant/all-NA trait(s): ",
            paste(trait_cols[!keep], collapse = ", "))
}
trait_cols <- trait_cols[keep]
n <- length(trait_cols)

if (n < 2) {
    write_placeholder(paste0("Not enough traits for a pairs plot (", n, " usable)"))
    quit(status = 0)
}
if (n > MAX_FACTORS) {
    write_placeholder(paste0(n, " traits exceeds the pairs-plot limit of ", MAX_FACTORS))
    quit(status = 0)
}

dt <- meta[, trait_cols, with = FALSE]

cell_theme <- theme_clinego(base_size = 9) +
    theme(plot.margin = margin(4, 4, 4, 4),
          axis.title = element_text(size = 8))

panel_density <- function(tr) {
    ggplot(dt, aes(x = .data[[tr]])) +
        geom_density(fill = CLINEGO_NEUTRAL, alpha = 0.5, color = CLINEGO_NEUTRAL) +
        labs(x = tr, y = NULL) +
        cell_theme
}

panel_scatter <- function(x, y) {
    ggplot(dt, aes(x = .data[[x]], y = .data[[y]])) +
        geom_point(color = CLINEGO_NEUTRAL, alpha = 0.7, size = 1.4) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                    color = CLINEGO_THRESHOLD, linewidth = 0.6) +
        labs(x = x, y = y) +
        cell_theme
}

panel_cor <- function(x, y) {
    r <- suppressWarnings(stats::cor(dt[[x]], dt[[y]], use = "pairwise.complete.obs"))
    lab <- if (is.na(r)) "NA" else sprintf("r = %.2f", r)
    size <- if (is.na(r)) 4 else 4 + 4 * abs(r)
    col  <- if (is.na(r)) "grey50" else if (r < 0) CLINEGO_REMOVED else CLINEGO_RETAINED
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = lab, size = size, color = col) +
        theme_void()
}

plots <- list()
for (i in seq_len(n)) {
    for (j in seq_len(n)) {
        # row i, column j — grid is filled row-major by ggarrange
        plots[[length(plots) + 1]] <- if (i == j) {
            panel_density(trait_cols[i])
        } else if (i > j) {
            panel_scatter(trait_cols[j], trait_cols[i])
        } else {
            panel_cor(trait_cols[i], trait_cols[j])
        }
    }
}

grid <- ggarrange(plotlist = plots, ncol = n, nrow = n)

cell_in <- 2.2
save_both(grid, width = n * cell_in, height = n * cell_in)
message("INFO: Wrote ", n, "x", n, " pairs grid for: ",
        paste(trait_cols, collapse = ", "))
