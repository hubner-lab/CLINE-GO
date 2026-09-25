# =============================================================================
# ms_style.R -- the shared look of every manuscript figure (simulation, Trifolium, ...).
# source() it; then use the MS_* colours, theme_ms() and ms_save(). Needs ggplot2 + svglite.
#
# PALETTE: ltc "minou" (Theodosiou), hex verbatim from ltc's R/ltc_functions.R, as in
# mvp_main_figure.R:45 (ltc is not in the image, so the six values are hardcoded).
# Scheme chosen by the user 2026-09-25: green/blue carry the focus, red/orange the contrast,
# grey the references.
#   focus     teal #00798c, sage #66a182, navy #2e4057
#   contrast  red #d1495b, amber #edae49, and MS_ORANGE halfway between them
#   reference grey #8d96a3 (zero lines, thresholds, "no effect")
#   sequential white -> sage -> teal -> navy (neutral magnitudes: counts, densities)
#   quality   red -> orange -> amber -> green (scores where higher is better; green = correct)
#   diverging red -> white -> teal (signed data only: bad -> neutral -> good)
#
# PRINT FORMAT: figures are drawn at final size in mm, 7 pt text, 6 pt tick labels, no titles
# (captions carry the explanation). ms_save() writes an EDITABLE SVG:
#   - fix_text_size = FALSE drops svglite's textLength/lengthAdjust, which otherwise re-stretch
#     every string to its original width when its font size is changed in an editor;
#   - svglite writes the font it resolved, and the image ships no Arial, only the metric-
#     compatible Nimbus Sans; font-family is rewritten to "Arial, Helvetica, sans-serif", so the
#     file opens in Arial with the layout measured on identical metrics.
# =============================================================================
MINOU <- c(teal = "#00798c", red = "#d1495b", amber = "#edae49",
           sage = "#66a182", navy = "#2e4057", grey = "#8d96a3")
MS_ORANGE    <- grDevices::colorRampPalette(c(MINOU[["red"]], MINOU[["amber"]]))(3)[2]
MS_REF       <- MINOU[["grey"]]
MS_INK       <- "#222222"
MS_DIVERGING  <- c(low = MINOU[["red"]], mid = "#FFFFFF", high = MINOU[["teal"]])
MS_SEQUENTIAL <- c("#FFFFFF", MINOU[["sage"]], MINOU[["teal"]], MINOU[["navy"]])
MS_QUALITY    <- c(MINOU[["red"]], MS_ORANGE, MINOU[["amber"]], MINOU[["sage"]])
# Categorical order for figures with no semantic grouping: focus hues first, then contrast.
MS_CATEGORICAL <- c(MINOU[["teal"]], MINOU[["red"]], MINOU[["amber"]], MINOU[["sage"]],
                    MINOU[["navy"]], MS_ORANGE, MINOU[["grey"]])

MS_BASE <- 7                       # pt: axis titles, strips
MS_TICK <- 6                       # pt: tick labels, legends, in-plot labels
MS_LAB  <- MS_TICK / ggplot2::.pt  # geom_text()/geom_label() size giving 6 pt

theme_ms <- function(base = MS_BASE, tick = MS_TICK) {
    ggplot2::theme_classic(base_size = base) + ggplot2::theme(
        text             = ggplot2::element_text(colour = MS_INK),
        axis.title       = ggplot2::element_text(size = base, colour = MS_INK),
        axis.text        = ggplot2::element_text(size = tick, colour = MS_INK),
        axis.line        = ggplot2::element_line(colour = MS_INK, linewidth = 0.3),
        axis.ticks       = ggplot2::element_line(colour = MS_INK, linewidth = 0.3),
        strip.background = ggplot2::element_blank(),
        strip.text       = ggplot2::element_text(size = base, colour = MS_INK),
        legend.position  = "bottom",
        legend.title     = ggplot2::element_text(size = tick, colour = MS_INK),
        legend.text      = ggplot2::element_text(size = tick, colour = MS_INK),
        legend.key.size  = ggplot2::unit(3, "mm"),
        legend.margin    = ggplot2::margin(0, 0, 0, 0),
        panel.spacing    = ggplot2::unit(3, "mm"),
        plot.margin      = ggplot2::margin(2, 3, 2, 2, "mm"))
}
#' White or MS_INK, whichever has the higher WCAG contrast on each background colour.
ms_text_on <- function(bg) {
    lum <- function(x) {
        v <- grDevices::col2rgb(x) / 255
        v <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
        colSums(v * c(0.2126, 0.7152, 0.0722))
    }
    l <- lum(bg)
    ifelse((1.05) / (l + 0.05) >= (l + 0.05) / (lum(MS_INK) + 0.05), "#FFFFFF", MS_INK)
}
scale_colour_ms <- function(...) ggplot2::scale_colour_manual(values = MS_CATEGORICAL, ...)
scale_fill_ms   <- function(...) ggplot2::scale_fill_manual(values = MS_CATEGORICAL, ...)

MS_SVG_FONT <- "font-family: Arial, Helvetica, sans-serif;"
#' Save `p` as <path_stem>.svg (editable, see header) and <path_stem>.png (600 dpi), w x h mm.
ms_save <- function(path_stem, p, w, h, dpi = 600) {
    svg <- paste0(path_stem, ".svg")
    ggplot2::ggsave(svg, p, width = w, height = h, units = "mm", device = svglite::svglite,
                    fix_text_size = FALSE, bg = "white")
    s <- readLines(svg, warn = FALSE)
    s <- gsub('font-family: "Nimbus Sans";', MS_SVG_FONT, s, fixed = TRUE)
    fams <- unique(unlist(regmatches(s, gregexpr("font-family:[^;]*;", s))))
    if (!identical(fams, MS_SVG_FONT) || any(grepl("textLength", s, fixed = TRUE)))
        stop(basename(svg), ": fonts [", paste(fams, collapse = " | "),
             "] or textLength left after the rewrite")
    writeLines(s, svg)
    ggplot2::ggsave(paste0(path_stem, ".png"), p, width = w, height = h, units = "mm",
                    dpi = dpi, bg = "white")
    invisible(svg)
}
