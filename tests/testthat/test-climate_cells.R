# scripts/R/lib/climate_cells.R — pair_climate_cells() (finding 2c9c67).
#
# The offset engines paired present and future raster cells by ROW: a future table in another
# order, or with another NoData mask, gave each cell another cell's future climate.

cells_tab <- function(ID, bio_1 = ID * 1.5, bio_2 = ID - 0.25) {
    data.table::data.table(ID = as.integer(ID), bio_1 = bio_1, bio_2 = bio_2)
}

test_that("identical tables give the old positional mask exactly", {
    pres <- cells_tab(c(3, 7, 9, 12, 20))
    fut  <- cells_tab(c(3, 7, 9, 12, 20))
    fut$bio_2[4] <- NA
    old_ok <- stats::complete.cases(pres[, -1]) & stats::complete.cases(fut[, -1])
    got <- quiet(pair_climate_cells(pres, fut, c("bio_1", "bio_2")))
    expect_identical(got$pres, which(old_ok))
    expect_identical(got$fut, which(old_ok))
})

test_that("a future table in another row order is paired by ID", {
    pres <- cells_tab(c(3, 7, 9, 12))
    fut  <- cells_tab(c(12, 3, 9, 7), bio_1 = c(120, 30, 90, 70))
    got <- quiet(pair_climate_cells(pres, fut, c("bio_1", "bio_2")))
    expect_identical(got$pres, 1:4)
    expect_identical(fut$ID[got$fut], pres$ID)
    expect_identical(fut$bio_1[got$fut], c(30, 70, 90, 120))
})

test_that("an equal count with a different ID set pairs only the shared cells, in present order", {
    # The case merge_climate_future.R's row-count check cannot see.
    pres <- cells_tab(c(3, 7, 9, 12))
    fut  <- cells_tab(c(3, 9, 12, 15))
    expect_message(got <- pair_climate_cells(pres, fut, c("bio_1", "bio_2")),
                   "1 present cell\\(s\\) without a future cell")
    expect_identical(pres$ID[got$pres], c(3L, 9L, 12L))
    expect_identical(fut$ID[got$fut], c(3L, 9L, 12L))
})

test_that("a cell with an NA predictor on either side is left out", {
    pres <- cells_tab(c(1, 2, 3, 4))
    fut  <- cells_tab(c(4, 3, 2, 1))
    pres$bio_1[2] <- NA                 # cell 2 incomplete in the present
    fut$bio_2[fut$ID == 4] <- NA        # cell 4 incomplete in the future
    got <- quiet(pair_climate_cells(pres, fut, c("bio_1", "bio_2")))
    expect_identical(pres$ID[got$pres], c(1L, 3L))
    expect_identical(fut$ID[got$fut], c(1L, 3L))
})

test_that("only the named predictors decide completeness", {
    pres <- cells_tab(c(1, 2)); pres$bio_9 <- c(NA, 1)
    fut  <- cells_tab(c(1, 2))
    got <- quiet(pair_climate_cells(pres, fut, c("bio_1", "bio_2")))
    expect_identical(got$pres, 1:2)
})

test_that("a duplicated cell ID or a missing column stops", {
    expect_error(pair_climate_cells(cells_tab(c(1, 2, 2)), cells_tab(1:3), c("bio_1", "bio_2")),
                 "cell ID 2 more than once")
    expect_error(pair_climate_cells(cells_tab(1:3), cells_tab(1:3)[, !"bio_2"], c("bio_1", "bio_2")),
                 "lacks column\\(s\\) bio_2")
})
