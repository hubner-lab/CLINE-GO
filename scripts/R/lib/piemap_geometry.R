# piemap_geometry.R — where and how a geo piemap draws its pies (plot_piemap.R).
# Usage: source("/pipeline/scripts/R/lib/piemap_geometry.R")
# Base R only.
#
# site_centroids()  one point per sampling site, the mean of its samples' coordinates (cdb1a1).
#                   A site whose samples carry jittered GPS coordinates was drawn as one pie per
#                   distinct coordinate, each showing the site mean. pregea_dbmem.R and
#                   mantel_test.R collapse a site to the same centroid.
# pie_x_stretch()   coord_sf draws longitude/latitude with ONE aspect correction for the whole
#                   panel, 1 / cos(mid-latitude of the panel's y range) (ggplot2's CoordSf). A pie
#                   drawn with one radius in degrees therefore came out 1/cos(lat) times taller
#                   than wide: 1.19x at 33 N (f3fcce). Its longitude radius is stretched by that
#                   factor instead.
# pie_wedges()      the pies as polygons, slices clockwise from 12 o'clock as scatterpie drew
#                   them, the longitude radius stretched.

# @param samples data.frame with site, latitude, longitude (one row per sample)
# @return data.frame site, latitude, longitude: one row per site that has a coordinate
site_centroids <- function(samples) {
    site <- as.character(samples$site)
    lat  <- tapply(as.numeric(samples$latitude),  site, mean, na.rm = TRUE)
    lon  <- tapply(as.numeric(samples$longitude), site, mean, na.rm = TRUE)
    out <- data.frame(site = names(lat), latitude = as.numeric(lat),
                      longitude = as.numeric(lon[names(lat)]), stringsAsFactors = FALSE)
    out[is.finite(out$latitude) & is.finite(out$longitude), , drop = FALSE]
}

# @param y_range the latitude range the panel trains on: the raster's and the pies' extent.
#   coord_sf expands it symmetrically, so its middle is the panel's.
# @param lonlat whether the map's CRS is longitude/latitude (terra::is.lonlat); a projected
#   map gets no correction from coord_sf, so its pies need none either
# @return the factor to multiply a pie's longitude radius by
pie_x_stretch <- function(y_range, lonlat = TRUE) {
    if (!isTRUE(lonlat)) return(1)
    1 / cos(mean(range(y_range)) * pi / 180)
}

# @param pies data.frame, one row per pie: x, y (centre), r (radius in latitude units), group
#   (pie id) and one numeric column per slice. A pie's slices are its share of the row's sum;
#   NA and non-positive slices are left out, and a pie with nothing left is not drawn.
# @param slice_cols the slice columns, in drawing (and legend) order
# @param x_stretch pie_x_stretch()
# @param n_full vertices of a full circle; a slice gets its share, at least 2
# @return data.frame pie, slice (factor, levels = slice_cols), wedge (polygon id), x, y
pie_wedges <- function(pies, slice_cols, x_stretch = 1, n_full = 120L) {
    out <- vector("list", 0L)
    for (i in seq_len(nrow(pies))) {
        amount <- as.numeric(unlist(pies[i, slice_cols, drop = TRUE]))
        amount[!is.finite(amount) | amount < 0] <- 0
        if (sum(amount) <= 0 || !is.finite(pies$r[i]) ||
            !is.finite(pies$x[i]) || !is.finite(pies$y[i])) next
        share  <- amount / sum(amount)
        ends   <- cumsum(share) * 2 * pi
        starts <- c(0, ends[-length(ends)])
        for (j in which(share > 0)) {
            theta <- seq(starts[j], ends[j], length.out = max(2L, ceiling(n_full * share[j])) + 1L)
            x <- pies$x[i] + pies$r[i] * x_stretch * sin(theta)
            y <- pies$y[i] + pies$r[i] * cos(theta)
            if (share[j] < 1) {   # a slice is a wedge from the centre; a whole pie is the circle
                x <- c(pies$x[i], x)
                y <- c(pies$y[i], y)
            }
            out[[length(out) + 1L]] <- data.frame(
                pie = pies$group[i], slice = slice_cols[j],
                wedge = paste(pies$group[i], slice_cols[j], sep = "\r"),
                x = x, y = y, stringsAsFactors = FALSE)
        }
    }
    if (length(out) == 0L)
        return(data.frame(pie = character(0), slice = factor(character(0), levels = slice_cols),
                          wedge = character(0), x = numeric(0), y = numeric(0)))
    res <- do.call(rbind, out)
    res$slice <- factor(res$slice, levels = slice_cols)
    res
}
