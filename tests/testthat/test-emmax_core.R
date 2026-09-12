# scripts/R/utils/emmax_core.R
#
# WHY THIS FILE SOURCES ITS OWN LIB INSTEAD OF USING helper-libs.R.
# run_emmax() dispatches through two rebindable string globals, EMMAX_BIN and
# EMMAX_RUN (:5-6), and the only hermetic way to exercise its stop-on-nonzero
# contract is to point them at /bin/false. test_dir() sources every test file
# into ONE shared parent environment, so rebinding them there would leak into
# every alphabetically-later file. sys.source() into a private env keeps it local.
#
# The headline target is load_pca_covariates(). Its k <= 0 guard exists to stop a
# silent-wrong-output bug documented in the code itself (:52-57): the LEA branch
# `1:k` -> `1:0` errors on index 0, but the PLINK branch `3:(2+k)` -> `3:2`
# does NOT — it quietly selects two wrong columns in reversed order. A guard that
# only prevents a crash on one of two branches is exactly the thing to pin.

emmax_env <- local({
    e <- new.env(parent = globalenv())
    sys.source(file.path(getOption("clinego.repo_root", "/pipeline"),
                         "scripts", "R", "utils", "emmax_core.R"),
               envir = e)
    e
})

# LEA projections: plain space-separated numeric matrix, no FID/IID prefix.
write_lea_pca <- function(d, n = 6, n_pc = 5, name = "projections.txt") {
    p <- file.path(d, name)
    m <- matrix(round(seq_len(n * n_pc) / 10, 3), nrow = n)
    write.table(m, p, sep = " ", row.names = FALSE, col.names = FALSE)
    p
}

# PLINK eigenvec: a constant FID column, a unique IID column, then the PCs.
write_plink_pca <- function(d, n = 6, n_pc = 5, name = "plink.eigenvec") {
    p <- file.path(d, name)
    m <- data.frame(FID = rep("FAM", n), IID = sprintf("ID%03d", seq_len(n)))
    for (i in seq_len(n_pc)) m[[paste0("V", i)]] <- round(seq_len(n) + i / 10, 3)
    write.table(m, p, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    p
}

test_that("the globals are plain strings, so nothing executes at source time", {
    expect_type(emmax_env$EMMAX_BIN, "character")
    expect_type(emmax_env$EMMAX_RUN, "character")
    # Sourcing this file must not require the binaries to exist or be runnable —
    # which is what made the "unsourceable" claim in helper-libs.R wrong.
    expect_length(emmax_env$EMMAX_BIN, 1L)
})

test_that("run_emmax stops when the wrapper exits non-zero", {
    e <- new.env(parent = emmax_env)
    e$EMMAX_RUN <- "/bin/false"
    e$EMMAX_BIN <- "ignored"
    f <- emmax_env$run_emmax
    environment(f) <- e
    expect_error(quiet(f("-v -d 10")), "EMMAX exited with status")
})

test_that("run_emmax returns quietly when the wrapper exits zero", {
    e <- new.env(parent = emmax_env)
    e$EMMAX_RUN <- "/bin/true"
    e$EMMAX_BIN <- "ignored"
    f <- emmax_env$run_emmax
    environment(f) <- e
    expect_identical(quiet(f("-v")), 0L)
})

test_that("load_pca_covariates returns NULL for k = 0 without reading the file", {
    # The kinship-only model in preGEA's #PC ladder. The path is deliberately
    # nonexistent: the guard must fire before any fread().
    got <- quiet(emmax_env$load_pca_covariates("/nonexistent/pca.txt", 0))
    expect_null(got)
})

test_that("load_pca_covariates returns NULL for a negative k too", {
    got <- quiet(emmax_env$load_pca_covariates("/nonexistent/pca.txt", -1))
    expect_null(got)
})

test_that("load_pca_covariates reads the LEA format from column 1", {
    d <- withr::local_tempdir()
    p <- write_lea_pca(d, n = 6, n_pc = 5)
    got <- quiet(emmax_env$load_pca_covariates(p, 3))

    expect_identical(names(got), c("PC1", "PC2", "PC3"))
    expect_identical(nrow(got), 6L)
    raw <- data.table::fread(p, sep = " ", header = FALSE)
    expect_equal(got$PC1, raw[[1]])
    expect_equal(got$PC3, raw[[3]])
})

test_that("load_pca_covariates autodetects PLINK eigenvec and skips FID/IID", {
    d <- withr::local_tempdir()
    p <- write_plink_pca(d, n = 6, n_pc = 5)
    got <- quiet(emmax_env$load_pca_covariates(p, 2))

    expect_identical(names(got), c("PC1", "PC2"))
    raw <- data.table::fread(p, sep = " ", header = FALSE)
    # Columns 3 and 4, in order — NOT columns 1-2 and NOT reversed.
    expect_equal(got$PC1, raw[[3]])
    expect_equal(got$PC2, raw[[4]])
})

test_that("load_pca_covariates subsets rows in DROP mode", {
    d <- withr::local_tempdir()
    p <- write_lea_pca(d, n = 6, n_pc = 5)
    all_samples <- sprintf("ID%03d", 1:6)
    order_file <- file.path(d, "samples_order.list")
    writeLines(all_samples, order_file)
    vcf_samples <- all_samples[c(2, 4, 5)]

    got <- quiet(emmax_env$load_pca_covariates(p, 2,
                                               samples_order_path = order_file,
                                               vcf_samples = vcf_samples))
    expect_identical(nrow(got), 3L)
    raw <- data.table::fread(p, sep = " ", header = FALSE)
    # Positional binding is the whole point: rows must follow the ORDER of the
    # samples_order file, not the order of vcf_samples.
    expect_equal(got$PC1, raw[[1]][c(2, 4, 5)])
})

test_that("load_pca_covariates errors when rows and VCF samples disagree", {
    d <- withr::local_tempdir()
    p <- write_lea_pca(d, n = 6, n_pc = 5)
    # 4 VCF samples, 6 PCA rows, and no samples_order_path to subset with.
    expect_error(
        quiet(emmax_env$load_pca_covariates(p, 2,
                                            vcf_samples = sprintf("ID%03d", 1:4))),
        "does not match VCF sample count")
})

test_that("load_pca_covariates accepts a matching row count with no subsetting", {
    d <- withr::local_tempdir()
    p <- write_lea_pca(d, n = 4, n_pc = 5)
    got <- quiet(emmax_env$load_pca_covariates(p, 2,
                                               vcf_samples = sprintf("ID%03d", 1:4)))
    expect_identical(nrow(got), 4L)
})
