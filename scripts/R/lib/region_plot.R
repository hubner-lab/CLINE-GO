# region_plot.R — what plot_regionplot.R's regional Manhattan computes besides the topr call:
# each series' significance line, the window drawn around a region, and readable gene labels.
# Usage: source("/pipeline/scripts/R/lib/region_plot.R")
# Requires: compute_pval_threshold() (scripts/R/utils/pval_threshold.R); ggplot2 and ggrepel
#   installed (repel_gene_labels).

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

# One line per distinct threshold, in the colour of the series it belongs to, or in shared_col
# where series of different colours share it (c34204: topr was given min(thresholds) only, so a
# method with a looser rule, RDA at q 0.05 beside EMMAX and LFMM at Bonferroni, had every hit
# drawn under the only line). Series without a threshold draw none.
# @return list(thresh, color), by increasing threshold, for topr's sign_thresh/sign_thresh_color
threshold_lines <- function(thresholds, series_colors, shared_col) {
    ok <- !is.na(thresholds)
    th <- sort(unique(thresholds[ok]))
    col <- vapply(th, function(t) {
        owners <- unique(series_colors[ok & thresholds == t])
        if (length(owners) == 1L) owners else shared_col
    }, character(1))
    list(thresh = th, color = unname(col))
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
