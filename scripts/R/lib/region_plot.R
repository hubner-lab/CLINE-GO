# region_plot.R — what plot_regionplot.R's regional Manhattan computes besides the topr call:
# each series' significance line, the window drawn around a region, and readable gene labels.
# Usage: source("/pipeline/scripts/R/lib/region_plot.R")
# Requires: compute_pval_threshold() (scripts/R/utils/pval_threshold.R), CLINEGO_SIGNIFICANCE
#   (scripts/R/utils/theme_clinego.R); ggplot2 and ggrepel installed (add_threshold_lines,
#   repel_gene_labels).

# The p-value threshold of one series (one method's p-values for one trait) under the rule given
# as "<rule>_<value>" (bonf_0.05, qval_0.05, top_100), from compute_pval_threshold(), the
# function behind the pipeline's own significance calls. plot_regionplot.R used a copy of its
# own that took every trait's threshold from the table's first trait column and counted NA rows
# in the Bonferroni denominator (c34204). NA when the rule yields no threshold.
series_threshold <- function(pvalues, adjust) {
    rule  <- sub("_.*$", "", adjust)
    value <- suppressWarnings(as.numeric(sub("^[^_]*_", "", adjust)))
    res <- compute_pval_threshold(pvalues, rule, value)
    if (identical(res$status, "ok")) res$threshold else NA_real_
}

# One line per distinct threshold (c34204: topr was given min(thresholds) only, so a method with a
# looser rule, RDA at q 0.05 beside EMMAX and LFMM at Bonferroni, had every hit drawn under the
# only line). Every line is drawn in one red, CLINEGO_SIGNIFICANCE, with a dash pattern of its own
# (F7: in the series' colours, a line shared by series of different colours took the plum of
# CLINEGO_THRESHOLD and no line was named). A line is named by the methods it carries: a method
# whose every series with a threshold is on this line by its name alone (the series legend
# already names each trait), otherwise by the series of it that are, in series order, wrapped
# at whole names. Nine EMMAX traits at one Bonferroni line read "EMMAX", not nine names in five
# rows, which pushed topr's Biotype legend over the last x tick label. Series without a threshold
# draw none.
# @param thresholds one p-value threshold per series, NA where its rule yields none
# @param labels the series' legend labels, in the same order
# @param methods the series' methods, in the same order; by default each label stands alone
# @return data.frame(thresh, label, linetype), one row per line, by increasing threshold (the
#   highest line first)
threshold_lines <- function(thresholds, labels, methods = labels) {
    ok <- !is.na(thresholds)
    th <- sort(unique(thresholds[ok]))
    lab <- vapply(th, function(t) {
        on <- ok & thresholds == t
        parts <- vapply(unique(methods[on]), function(m) {
            if (all(on[ok & methods == m])) m
            else paste(unique(labels[on & methods == m]), collapse = ", ")
        }, character(1), USE.NAMES = FALSE)
        wrap_series_names(parts)
    }, character(1))
    data.frame(thresh = th, label = lab, linetype = threshold_linetypes(length(th)),
               stringsAsFactors = FALSE)
}

# Series names joined by ", ", a new line wherever the next name would pass width characters:
# a Bonferroni line shared by nine EMMAX traits names all nine without one 130-character line.
# A name is never split.
wrap_series_names <- function(names, width = 30L) {
    out <- character(0)
    cur <- ""
    for (s in names) {
        cand <- if (nzchar(cur)) paste0(cur, ", ", s) else s
        if (nzchar(cur) && nchar(cand) > width) {
            out <- c(out, paste0(cur, ","))
            cur <- s
        } else {
            cur <- cand
        }
    }
    paste(c(out, cur), collapse = "\n")
}

# The dash patterns the threshold lines take in turn, highest line first: ggplot2's named dashes,
# then on/off hex patterns unlike them. Never "solid" (no dash) or "dotted" ("13": the region's
# edges, ad9f15). Past these, every two-digit hex pattern not already taken, so each line keeps a
# pattern of its own however many there are (above four or five they are hard to tell apart).
THRESHOLD_LINETYPES <- c("dashed", "dotdash", "longdash", "twodash", "F282", "131373", "22", "F3")

threshold_linetypes <- function(n) {
    hex  <- c(1:9, LETTERS[1:6])
    more <- setdiff(as.vector(outer(hex, hex, paste0)), c("44", "13", "73", THRESHOLD_LINETYPES))
    c(THRESHOLD_LINETYPES, more)[seq_len(n)]
}

# Draws threshold_lines()' lines on one of topr's panels, each in CLINEGO_SIGNIFICANCE with its
# own dash pattern. legend = TRUE (the main panel) keys every pattern to its line's label in a
# legend under the series legend; the overview panel draws the same lines without one.
# mirror = TRUE draws each line below zero as well, as topr does once it draws the series past
# its ntop below the axis.
# @param p topr's main or overview panel (regionplot(extract_plots = TRUE))
# @param lines threshold_lines()' data.frame; p is returned unchanged when it has no row
add_threshold_lines <- function(p, lines, linewidth, legend = TRUE, mirror = FALSE) {
    if (nrow(lines) == 0L) return(p)
    y <- -log10(lines$thresh)
    d <- data.frame(yintercept = if (mirror) c(y, -y) else y,
                    label = rep(lines$label, if (mirror) 2L else 1L), stringsAsFactors = FALSE)
    guide <- if (legend) {
        ggplot2::guide_legend(order = 2,
                              theme = ggplot2::theme(legend.key.width = grid::unit(2.5, "lines")))
    } else "none"
    # show.legend NA, not TRUE: TRUE also draws the line across every key of the series legend
    p <- p +
        ggplot2::geom_hline(data = d, ggplot2::aes(yintercept = yintercept, linetype = label),
                            colour = CLINEGO_SIGNIFICANCE, linewidth = linewidth,
                            show.legend = if (legend) NA else FALSE) +
        ggplot2::scale_linetype_manual(values = stats::setNames(lines$linetype, lines$label),
                                       breaks = lines$label, name = "Significance threshold",
                                       guide = guide)
    # The series legend first: a guide left at order 0 is placed after every numbered one.
    if (legend) p <- p + ggplot2::guides(colour = ggplot2::guide_legend(order = 1))
    p
}

# The window drawn for a region: flank bp on either side, starting at 1 at the latest (ad9f15).
region_window <- function(start, end, flank = 0) {
    c(max(1, start - flank), end + flank)
}

# The region's edges, wherever the window shows more than the region (ad9f15): drawn as vertical
# lines, so genes just outside the region can be judged against its bounds. NULL when the window
# is the region.
region_edge_lines <- function(region_start, region_end, window_start, window_end) {
    c(if (window_start < region_start) region_start, if (window_end > region_end) region_end)
}

# topr writes gene names with geom_text at each gene's start, on its own box, so close tandem
# copies print over each other (e9466a: two RCC1 labels in Trifolium 5_5557888). Each gene-name
# layer is redrawn with ggrepel: the label starts just right of its box, moves only up or down
# away from its neighbours, every label is kept (max.overlaps = Inf) and a thin white halo keeps
# it readable over a neighbouring box. Genes without a name do not push the labels.
# Port of the derived manuscript copy (manuscript/trifolium/.../fig3/plot_regionplot_topr_ms.R).
# @param gene_plot topr's gene panel (regionplot(extract_plots = TRUE)$gene_plot)
# @return the panel with every GeomText layer replaced; unchanged when it has none
repel_gene_labels <- function(gene_plot, colour, default_size = 3) {
    is_txt <- vapply(gene_plot$layers, function(l) inherits(l$geom, "GeomText"), logical(1))
    for (i in which(is_txt)) {
        l <- gene_plot$layers[[i]]
        if (!all(c("gene_symbol", "gene_start", "gene_end") %in% names(l$data)))
            stop("topr's gene-name layer has no gene_symbol / gene_start / gene_end columns")
        lab <- l$data[!is.na(l$data$gene_symbol) & nzchar(l$data$gene_symbol), , drop = FALSE]
        span <- diff(range(c(l$data$gene_start, l$data$gene_end), na.rm = TRUE))
        lab$gene_start <- lab$gene_end + 0.004 * span
        size <- l$aes_params$size
        gene_plot$layers[[i]] <- ggrepel::geom_text_repel(
            data = lab, mapping = l$mapping, inherit.aes = FALSE,
            size = if (is.null(size)) default_size else size, hjust = 0, vjust = 0,
            colour = colour, direction = "y", bg.colour = "white", bg.r = 0.12,
            box.padding = 0.05, point.padding = 0, min.segment.length = Inf,
            max.overlaps = Inf, seed = 1)
    }
    gene_plot
}
