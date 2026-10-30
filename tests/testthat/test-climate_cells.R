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

# ---- climate_cell_table() (f11151) ------------------------------------------------------------

# A 300 x 400 three-layer raster: 120,000 cells, so IDs reach 1e5; 5,000 cells NA in every layer
# (sea), 300 NA in one layer only, and cell 100000 kept on land.
fx_cell_raster <- function() {
    set.seed(11)
    r <- terra::rast(nrows = 300, ncols = 400, nlyrs = 3, xmin = 0, xmax = 40, ymin = 0, ymax = 30)
    v <- matrix(runif(terra::ncell(r) * 3) * 100, ncol = 3)
    v[sample(terra::ncell(r), 5000), ] <- NA
    v[sample(terra::ncell(r), 300), 2] <- NA
    v[100000, ] <- c(1.5, 2.5, 3.5)
    terra::values(r) <- v
    names(r) <- c("bio_1", "bio_2", "bio_12")
    r
}

test_that("the cell table equals the old whole-grid extract minus its NA rows", {
    r <- fx_cell_raster()
    old <- terra::extract(r, seq_len(terra::ncell(r)))
    old$ID <- seq_len(terra::ncell(r))
    old <- stats::na.omit(old)[, c("ID", names(r))]
    got <- climate_cell_table(r)
    expect_identical(names(got), c("ID", "bio_1", "bio_2", "bio_12"))
    expect_type(got$ID, "integer")
    expect_identical(nrow(got), nrow(old))
    for (col in names(old)) expect_identical(got[[col]], old[[col]], info = col)
})

test_that("the table does not depend on the block size, an all-NoData block included", {
    r <- fx_cell_raster()
    v <- terra::values(r)
    v[1:1200, ] <- NA                          # rows 1-3: one whole block of 2 rows is sea
    terra::values(r) <- v
    one_block <- climate_cell_table(r, cells_per_block = terra::ncell(r))
    for (per_block in c(800, 1000, 30001)) {   # 2 rows, 2 rows, 75 rows per block (last one short)
        expect_identical(climate_cell_table(r, cells_per_block = per_block), one_block, info = per_block)
    }
    expect_identical(one_block$ID, which(stats::complete.cases(v)))
})

test_that("a cell NA in one layer is dropped, and the written table matches the old bytes", {
    r <- fx_cell_raster()
    got <- climate_cell_table(r)
    v <- terra::values(r)
    expect_identical(got$ID, which(stats::complete.cases(v)))
    expect_true(100000L %in% got$ID)
    old <- terra::extract(r, seq_len(terra::ncell(r)))
    old$ID <- seq_len(terra::ncell(r))
    old <- stats::na.omit(old)[, c("ID", names(r))]
    f_old <- tempfile(fileext = ".tsv"); f_new <- tempfile(fileext = ".tsv")
    on.exit(unlink(c(f_old, f_new)))
    data.table::fwrite(old, f_old, sep = "\t")
    data.table::fwrite(got, f_new, sep = "\t")
    expect_identical(tools::md5sum(f_new)[[1]], tools::md5sum(f_old)[[1]])
    expect_true(any(readLines(f_new) == "100000\t1.5\t2.5\t3.5"))
})
