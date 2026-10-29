# climate_cells.R — pair present and future raster cells by their cell ID.
#
# climate_present_all.tsv and climate_future_*_all.tsv hold one row per raster cell, keyed by
# `ID`, the terra cell index. A future raster with another NoData mask drops another set of
# cells, so the two tables need not list the same cells in the same rows. The three offset
# engines used to pair them by row position: two complete.cases() masks of different length
# were recycled, and an equal count with a different ID set paired one cell's present with
# another cell's future (2c9c67). Used by rda_offset.R, geometric_offset.R and
# gradient_forest_offset.R.

# Row indices of the cells present in BOTH tables with every predictor non-NA on both sides,
# in the present table's row order: list(pres = rows of `present`, fut = rows of `future`).
# Neither table is copied; on two tables listing the same cells in the same rows the result
# is exactly the old positional mask, which(complete(present) & complete(future)).
pair_climate_cells <- function(present, future, predictors) {
    for (tab in list(present, future)) {
        missing_cols <- setdiff(c("ID", predictors), names(tab))
        if (length(missing_cols) > 0) {
            stop("pair_climate_cells: a climate table lacks column(s) ",
                 paste(missing_cols, collapse = ", "))
        }
        if (anyDuplicated(tab$ID)) {
            stop("pair_climate_cells: a climate table lists raster cell ID ",
                 tab$ID[anyDuplicated(tab$ID)], " more than once")
        }
    }
    complete <- function(tab) do.call(stats::complete.cases, lapply(predictors, function(p) tab[[p]]))
    fut <- match(present$ID, future$ID)
    ok  <- !is.na(fut) & complete(present)
    ok[ok] <- complete(future)[fut[ok]]
    n_unmatched <- sum(is.na(fut))
    if (n_unmatched > 0 || nrow(future) != nrow(present)) {
        message("WARNING: present and future climate list different raster cells (present ",
                nrow(present), ", future ", nrow(future), ", ", n_unmatched,
                " present cell(s) without a future cell): cells are paired by ID and the ",
                "unpaired ones get no offset")
    }
    list(pres = which(ok), fut = fut[ok])
}
