# factor_order.R — one display order for named factors in every figure that lists them.
# Usage: source("/pipeline/scripts/R/utils/factor_order.R")
#
# Names ending in _<N> (bio_1 ... bio_19, trait_2) come first, in numeric order; the others
# follow in the order given (climate_multivariate after the bioclim series, traits after the
# predictors). A plain sort puts bio_10 to bio_19 between bio_1 and bio_2 (9e5dd0: the region
# plots' series and legend; the climate correlogram's axes the same way).
#
# Callers: plot_regionplot.R (series), plot_correlation_heatmap.R (axes), plot_density.R (panels).

# @param x character vector of factor names
# @return integer permutation: x[numbered_order(x)] is the display order
numbered_order <- function(x) {
    x <- as.character(x)
    numbered <- grepl("_[0-9]+$", x)
    n <- rep(Inf, length(x))
    n[numbered] <- as.numeric(sub("^.*_([0-9]+)$", "\\1", x[numbered]))
    order(n, seq_along(x))
}
