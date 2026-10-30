# climate_cells.R — the raster cell tables: build one, pair present and future cells by ID.
#
# climate_present_all.tsv and climate_future_*_all.tsv hold one row per raster cell, keyed by
# `ID`, the terra cell index. climate_cell_table() builds them (download_climate_present.R,
# merge_climate_future.R). A future raster with another NoData mask drops another set of
# cells, so the two tables need not list the same cells in the same rows. The three offset
# engines used to pair them by row position: two complete.cases() masks of different length
# were recycled, and an equal count with a different ID set paired one cell's present with
# another cell's future (2c9c67). pair_climate_cells() is used by rda_offset.R,
# geometric_offset.R and gradient_forest_offset.R.

# One row per cell with a value in EVERY layer of `r` (na.omit()'s rule): `ID` (the cell index),
# then one column per layer, in cell order. The extractors used terra::extract(r, 1:ncell(r)),
# which built the whole grid, NoData cells included, and dropped them with na.omit() one line
# later (f11151); terra::as.data.frame(na.rm = TRUE) does the same internally (terra 1.9).
# Reading `cells_per_block` cells at a time keeps only land rows: on a 20M-cell crop at
# 30 arc-sec the step took 35 s and 6.7 GB at peak, against 117 s and 12.9 GB. ID stays an
# integer: fwrite() writes a double 100000 as 1e+05.
climate_cell_table <- function(r, cells_per_block = 1e6) {
    n_col  <- terra::ncol(r)
    n_row  <- terra::nrow(r)
    step   <- max(1L, as.integer(cells_per_block %/% n_col))
    terra::readStart(r)
    on.exit(terra::readStop(r))
    blocks <- lapply(seq(1L, n_row, by = step), function(row) {
        v <- terra::readValues(r, row = row, nrows = min(step, n_row - row + 1L), mat = TRUE)
        keep <- which(stats::complete.cases(v))
        data.table::data.table(ID = as.integer((row - 1L) * n_col + keep), v[keep, , drop = FALSE])
    })
    cells <- data.table::rbindlist(blocks)
    data.table::setnames(cells, c("ID", names(r)))
    cells
}

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
