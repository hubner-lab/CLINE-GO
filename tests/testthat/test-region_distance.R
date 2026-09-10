# scripts/R/lib/region_distance.R
#
# Turns Association.snp_clumping_distance into an actual bp distance, either fixed
# or inverted out of the LD-decay table written by mode=structure.

# --- fixtures --------------------------------------------------------------

# Builder, not a shared object: resolve_clumping_distance() reads the table with
# fread() every call, but the rows are also handed to resolve_row_distance() as
# lists and nothing here should be able to leak between tests.
ld_row <- function(group = "All", scope = "genome_wide", method = "hill_weir",
                   C_hat = 0.001, n_samples = 50, max_dist_bp = 500000,
                   half_decay_bp = NA_real_, r2_02_bp = NA_real_,
                   r2_intercept = NA_real_) {
    data.table::data.table(
        group = group, scope = scope, method = method,
        C_hat = C_hat, n_samples = n_samples, max_dist_bp = max_dist_bp,
        half_decay_bp = half_decay_bp, r2_02_bp = r2_02_bp,
        r2_intercept = r2_intercept
    )
}

write_ld <- function(dt) {
    path <- withr::local_tempfile(fileext = ".tsv", .local_envir = parent.frame())
    data.table::fwrite(dt, path, sep = "\t")
    path
}

# --- hill_weir -------------------------------------------------------------

test_that("hill_weir matches the closed form and decays with distance", {
    # Hand-evaluated at x = C*d = 0: ((10+0)/((2)(11))) * (1 + (3*12)/(n*2*11))
    n <- 50
    expect_equal(hill_weir(0, C = 0.001, n = n),
                 (10 / 22) * (1 + 36 / (n * 22)))

    d <- c(1, 1e3, 1e4, 1e5)
    r <- hill_weir(d, C = 1e-3, n = 50)
    expect_length(r, 4L)                       # vectorised over d
    expect_true(all(diff(r) < 0))              # monotone decreasing
    expect_true(all(r > 0 & r < 1))
})

# --- invert_hill_weir ------------------------------------------------------

test_that("invert_hill_weir round-trips through hill_weir", {
    d <- invert_hill_weir(C_hat = 1e-3, n_samples = 50, r2_target = 0.2,
                          max_dist_bp = 5e5)
    expect_false(is.na(d))
    expect_equal(hill_weir(d, C = 1e-3, n = 50), 0.2, tolerance = 1e-6)
})

test_that("invert_hill_weir returns NA on missing or zero-length input", {
    expect_true(is.na(invert_hill_weir(NA_real_, 50, 0.2, 5e5)))
    expect_true(is.na(invert_hill_weir(1e-3, NA_real_, 0.2, 5e5)))
    expect_true(is.na(invert_hill_weir(1e-3, 50, 0.2, NA_real_)))
    expect_true(is.na(invert_hill_weir(numeric(0), 50, 0.2, 5e5)))
})

test_that("invert_hill_weir warns and returns NA when the root is not bracketed", {
    # r2_target = 0.99 is above hill_weir() everywhere on [1, max_dist], so uniroot
    # cannot bracket a sign change.
    d <- expect_message(invert_hill_weir(1e-3, 50, 0.99, 5e5), "inversion failed")
    expect_true(is.na(d))
})

# --- loess_fallback --------------------------------------------------------

test_that("loess_fallback picks whichever precomputed column is nearer the target", {
    # r2_intercept 0.8 -> half_r2 0.4. Target 0.35 is nearer 0.4 than 0.2.
    row <- as.list(ld_row(half_decay_bp = 1000, r2_02_bp = 9000, r2_intercept = 0.8))
    expect_equal(suppressMessages(loess_fallback(row, 0.35)), 1000)

    # Target 0.2 is exactly r2_02's own r2, so that column wins.
    expect_equal(suppressMessages(loess_fallback(row, 0.2)), 9000)
})

test_that("loess_fallback defaults half_r2 to 0.3 when r2_intercept is NA", {
    row <- as.list(ld_row(half_decay_bp = 1000, r2_02_bp = 9000))
    expect_equal(suppressMessages(loess_fallback(row, 0.29)), 1000)  # nearer 0.3 than 0.2
})

test_that("loess_fallback passes a lone non-NA column through, else returns NA", {
    expect_equal(loess_fallback(as.list(ld_row(r2_02_bp = 4321)), 0.2), 4321)
    expect_equal(loess_fallback(as.list(ld_row(half_decay_bp = 1234)), 0.2), 1234)
    expect_true(is.na(loess_fallback(as.list(ld_row()), 0.2)))
})

# --- resolve_row_distance --------------------------------------------------

test_that("resolve_row_distance inverts when method is hill_weir and falls back otherwise", {
    hw <- as.list(ld_row(C_hat = 1e-3, n_samples = 50, max_dist_bp = 5e5))
    expect_equal(resolve_row_distance(hw, 0.2),
                 round(invert_hill_weir(1e-3, 50, 0.2, 5e5)))

    lo <- as.list(ld_row(method = "loess", r2_02_bp = 7777))
    expect_equal(resolve_row_distance(lo, 0.2), 7777)
})

test_that("resolve_row_distance falls back when the hill_weir fit is unusable", {
    row <- as.list(ld_row(C_hat = NA_real_, r2_02_bp = 8888))
    expect_equal(suppressMessages(resolve_row_distance(row, 0.2)), 8888)
})

# --- resolve_clumping_distance: fixed specs --------------------------------

test_that("a numeric spec is returned as an integer, unchanged", {
    expect_identical(suppressMessages(resolve_clumping_distance("50000")), 50000L)
    expect_identical(suppressMessages(resolve_clumping_distance(50000)), 50000L)
})

test_that("a non-numeric non-auto spec warns and falls back to 1e6", {
    d <- expect_message(resolve_clumping_distance("wide"), "not numeric")
    expect_identical(d, 1000000L)
})

# --- resolve_clumping_distance: auto specs ---------------------------------

test_that("auto specs are matched case-insensitively", {
    path <- write_ld(ld_row())
    a <- suppressMessages(resolve_clumping_distance("AUTO_GENOME_WIDE", path))
    b <- suppressMessages(resolve_clumping_distance("auto_genome_wide", path))
    expect_identical(a, b)
})

test_that("'auto' is the deprecated alias for auto_genome_wide", {
    path <- write_ld(ld_row())
    expect_message(resolve_clumping_distance("auto", path), "deprecated alias")
    expect_identical(suppressMessages(resolve_clumping_distance("auto", path)),
                     suppressMessages(resolve_clumping_distance("auto_genome_wide", path)))
})

test_that("an auto spec with no LD table errors instead of guessing", {
    expect_error(resolve_clumping_distance("auto_genome_wide", "NULL"),
                 "LD decay table not found")
    expect_error(resolve_clumping_distance("auto_genome_wide", "/nonexistent/ld.tsv"),
                 "LD decay table not found")
})

test_that("an empty LD table warns and falls back to 1e6", {
    path <- write_ld(ld_row()[0])
    d <- expect_message(resolve_clumping_distance("auto_genome_wide", path), "is empty")
    expect_identical(d, 1000000L)
})

test_that("auto_genome_wide errors when there is no genome_wide row", {
    path <- write_ld(ld_row(scope = "1"))
    expect_error(suppressMessages(resolve_clumping_distance("auto_genome_wide", path)),
                 "No genome-wide LD decay row")
})

test_that("auto_genome_wide returns the scalar distance implied by the row", {
    path <- write_ld(ld_row(C_hat = 1e-3, n_samples = 50, max_dist_bp = 5e5))
    d <- suppressMessages(resolve_clumping_distance("auto_genome_wide", path, 0.2))
    expect_identical(d, as.integer(round(invert_hill_weir(1e-3, 50, 0.2, 5e5))))
    expect_length(d, 1L)
    expect_null(names(d))
})

test_that("auto_per_chromosome returns a NAMED integer vector keyed by chromosome", {
    path <- write_ld(rbind(
        ld_row(scope = "genome_wide", C_hat = 1e-3),
        ld_row(scope = "1",  C_hat = 5e-4),
        ld_row(scope = "2",  C_hat = 2e-3)
    ))
    d <- suppressMessages(resolve_clumping_distance("auto_per_chromosome", path, 0.2))
    expect_identical(sort(names(d)), c("1", "2"))    # genome_wide is not a chromosome
    expect_true(is.integer(d))
    # Slower decay (smaller C) means LD extends further.
    expect_gt(d[["1"]], d[["2"]])
})

test_that("auto_per_chromosome falls back to the genome-wide distance per chromosome", {
    # chr 2's row cannot be inverted (C_hat NA, no precomputed columns), so it must
    # inherit the genome-wide value rather than dropping out of the map.
    path <- write_ld(rbind(
        ld_row(scope = "genome_wide", C_hat = 1e-3),
        ld_row(scope = "1", C_hat = 5e-4),
        ld_row(scope = "2", C_hat = NA_real_)
    ))
    d <- suppressMessages(resolve_clumping_distance("auto_per_chromosome", path, 0.2))
    gw <- as.integer(round(invert_hill_weir(1e-3, 50, 0.2, 5e5)))
    expect_identical(d[["2"]], gw)
})

test_that("auto_per_chromosome errors when a chromosome fails and there is no fallback", {
    path <- write_ld(ld_row(scope = "1", C_hat = NA_real_))
    expect_error(suppressMessages(resolve_clumping_distance("auto_per_chromosome", path)),
                 "no genome-wide fallback")
})
