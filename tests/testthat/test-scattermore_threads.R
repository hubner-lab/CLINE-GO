# scripts/R/utils/scattermore_threads.R (de7512). scattermore's point kernel defaulted to
# threads = 0, one thread per host core, past every cap a job sets: a Manhattan save ran 194
# threads in a 1-thread job on the 192-core lab. The util sets the default to OMP_NUM_THREADS,
# Snakemake's per-job thread count, and 1 when it is unset.
#
# Sourced here, not in helper-libs.R: each test sources it under its own OMP_NUM_THREADS. The
# change is to scattermore's namespace, so it holds for the rest of the run; the last source
# below leaves the default at 1.

.smt_util <- file.path(.clinego_R, "utils", "scattermore_threads.R")
.smt_default <- function() formals(get("apply_kernel_rgbwt", envir = asNamespace("scattermore")))$threads

test_that("the kernel still takes a threads argument (a scattermore upgrade must fail here)", {
    expect_true("threads" %in% names(formals(scattermore::apply_kernel_rgbwt)))
})

test_that("the kernel's default is OMP_NUM_THREADS, or 1 when it is unset or not a count", {
    withr::with_envvar(c(OMP_NUM_THREADS = "3"), source(.smt_util, local = new.env()))
    expect_identical(.smt_default(), 3L)
    withr::with_envvar(c(OMP_NUM_THREADS = "x"), source(.smt_util, local = new.env()))
    expect_identical(.smt_default(), 1L)
    withr::with_envvar(c(OMP_NUM_THREADS = NA), source(.smt_util, local = new.env()))
    expect_identical(.smt_default(), 1L)
})

test_that("the blurred raster is the same for any thread count", {
    set.seed(7)
    xy  <- cbind(runif(5000), runif(5000))
    pts <- scattermore::scatter_points_rgbwt(xy, xlim = c(0, 1), ylim = c(0, 1),
                                             out_size = c(300, 200), RGBA = c(10, 120, 90, 200))
    one <- scattermore::apply_kernel_rgbwt(pts, radius = 3, threads = 1)
    expect_identical(scattermore::apply_kernel_rgbwt(pts, radius = 3, threads = 0), one)
    expect_identical(scattermore::apply_kernel_rgbwt(pts, radius = 3, threads = 4), one)
})

test_that("every script that draws geom_scattermore() sources the util", {
    scripts <- list.files(file.path(.clinego_root, "scripts"), pattern = "\\.R$", full.names = TRUE)
    code <- lapply(scripts, function(f) grep("^\\s*#", readLines(f, warn = FALSE), value = TRUE, invert = TRUE))
    draws <- vapply(code, function(l) any(grepl("geom_scattermore\\(|add_scatter_layer\\(", l)), logical(1))
    caps  <- vapply(code, function(l) any(grepl("source\\(\"/pipeline/scripts/R/utils/scattermore_threads\\.R\"\\)", l)),
                    logical(1))
    expect_gt(sum(draws), 0)
    expect_identical(basename(scripts[draws & !caps]), character(0))
})
