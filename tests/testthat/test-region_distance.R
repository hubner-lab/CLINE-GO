# scripts/R/lib/region_distance.R
#
# Turns Association.snp_clumping_distance into an actual bp distance, either fixed
# or read from the LD-decay table written by mode=structure under one of three
# rules: critical (decay_critical_bp), half_decay (decay_half_bp) or r2_threshold
# (decay_r2thr_bp).

# --- fixtures --------------------------------------------------------------

# Builder, not a shared object: resolve_clumping_distance() reads the table with
# fread() every call, but the rows are also handed to resolve_row_distance() as
# lists and nothing here should be able to leak between tests.
ld_row <- function(group = "All", scope = "genome_wide", critical = 3100,
                   half = 1240, r2thr = 2600, r2_threshold = 0.2) {
    data.table::data.table(
        group = group, scope = scope, n_samples = 50L, n_pairs = 1e5,
        max_dist_bp = 5e5, short_r2 = 0.45, background_r2 = 0.06,
        background_mean_r2 = 0.02, decay_critical_bp = critical,
        decay_half_bp = half, decay_r2thr_bp = r2thr, r2_threshold = r2_threshold,
        decay_method = "loess_log"
    )
}

write_ld <- function(dt) {
    path <- withr::local_tempfile(fileext = ".tsv", .local_envir = parent.frame())
    data.table::fwrite(dt, path, sep = "\t")
    path
}

# --- rules -----------------------------------------------------------------

test_that("resolve_row_distance reads the column of each rule", {
    row <- as.list(ld_row())
    expect_identical(resolve_row_distance(row, "critical"), 3100)
    expect_identical(resolve_row_distance(row, "half_decay"), 1240)
    expect_identical(resolve_row_distance(row, "r2_threshold"), 2600)
})

test_that("each rule resolves auto_genome_wide to its own column", {
    path <- write_ld(ld_row())
    r <- function(rule) suppressMessages(resolve_clumping_distance("auto_genome_wide", path, rule))
    expect_identical(r("critical"), 3100L)
    expect_identical(r("half_decay"), 1240L)
    expect_identical(r("r2_threshold"), 2600L)
    expect_identical(r("CRITICAL"), 3100L)     # rule names are case-insensitive
})

test_that("the default rule is critical", {
    path <- write_ld(ld_row())
    expect_identical(suppressMessages(resolve_clumping_distance("auto_genome_wide", path)), 3100L)
})

test_that("a numeric rule is the r2_threshold rule, and only at the table's level", {
    path <- write_ld(ld_row(r2_threshold = 0.2))
    expect_identical(suppressMessages(resolve_clumping_distance("auto_genome_wide", path, 0.2)),
                     2600L)
    expect_identical(suppressMessages(resolve_clumping_distance("auto_genome_wide", path, "0.2")),
                     2600L)
    expect_error(suppressMessages(resolve_clumping_distance("auto_genome_wide", path, 0.1)),
                 "LDdecay.r2_threshold")
})

test_that("an unknown rule name is an error, not a fallback", {
    path <- write_ld(ld_row())
    expect_error(resolve_clumping_distance("auto_genome_wide", path, "extent"),
                 "unknown LD-decay rule")
})

test_that("a rule that is never reached stops and names the remedies", {
    path <- write_ld(ld_row(critical = NA_real_))
    expect_error(suppressMessages(resolve_clumping_distance("auto_genome_wide", path, "critical")),
                 "not reached.*LDdecay.max_distance")
    # the other rules of the same row still resolve
    expect_identical(suppressMessages(resolve_clumping_distance("auto_genome_wide", path, "half_decay")),
                     1240L)
})

test_that("a curve that starts below the rule's level is named, not blamed on max_distance", {
    # a small group whose unlinked background is as high as its short-range LD
    row <- ld_row(critical = NA_real_)
    row[, background_r2 := 0.6]
    path <- write_ld(row)
    expect_error(suppressMessages(resolve_clumping_distance("auto_genome_wide", path, "critical")),
                 "starts at or below that level")
})

test_that("a group missing from the table is named, with the groups present", {
    path <- write_ld(ld_row(group = "All"))
    expect_error(suppressMessages(resolve_clumping_distance("auto_genome_wide", path,
                                                            "critical", group = "Tel Aviv")),
                 "group 'Tel Aviv' is not in .*groups: All")
})

test_that("a table from before the PLINK estimator is refused, not misread", {
    old <- data.table::data.table(group = "All", scope = "genome_wide", method = "hill_weir",
                                  C_hat = 1e-3, n_samples = 50, max_dist_bp = 5e5,
                                  half_decay_bp = 1500, r2_02_bp = 3000, r2_intercept = 0.2)
    path <- write_ld(old)
    expect_error(resolve_clumping_distance("auto_genome_wide", path), "re-run mode=structure")
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

test_that("an empty LD table is an error, not a silent 1e6 bp", {
    path <- write_ld(ld_row()[0])
    expect_error(resolve_clumping_distance("auto_genome_wide", path), "is empty")
})

test_that("auto_genome_wide errors when there is no genome_wide row", {
    path <- write_ld(ld_row(scope = "1"))
    expect_error(suppressMessages(resolve_clumping_distance("auto_genome_wide", path)),
                 "No genome-wide LD decay row.*LDdecay.scope")
})

test_that("auto_genome_wide returns a scalar", {
    path <- write_ld(ld_row())
    d <- suppressMessages(resolve_clumping_distance("auto_genome_wide", path))
    expect_length(d, 1L)
    expect_null(names(d))
})

test_that("auto_per_chromosome returns a NAMED integer vector keyed by chromosome", {
    path <- write_ld(rbind(
        ld_row(scope = "genome_wide", critical = 3000),
        ld_row(scope = "1",  critical = 4000),
        ld_row(scope = "2",  critical = 2000)
    ))
    d <- suppressMessages(resolve_clumping_distance("auto_per_chromosome", path))
    expect_identical(sort(names(d)), c("1", "2"))    # genome_wide is not a chromosome
    expect_true(is.integer(d))
    expect_identical(d[["1"]], 4000L)
    expect_identical(d[["2"]], 2000L)
})

test_that("auto_per_chromosome falls back to the genome-wide distance per chromosome", {
    # chr 2 never reaches the critical r2, so it inherits the genome-wide value
    # rather than dropping out of the map.
    path <- write_ld(rbind(
        ld_row(scope = "genome_wide", critical = 3000),
        ld_row(scope = "1", critical = 4000),
        ld_row(scope = "2", critical = NA_real_)
    ))
    d <- suppressMessages(resolve_clumping_distance("auto_per_chromosome", path))
    expect_identical(d[["2"]], 3000L)
})

test_that("auto_per_chromosome errors when a chromosome fails and there is no fallback", {
    path <- write_ld(ld_row(scope = "1", critical = NA_real_))
    expect_error(suppressMessages(resolve_clumping_distance("auto_per_chromosome", path)),
                 "no genome-wide value to fall back to")
})

test_that("auto_per_chromosome without per-chromosome rows is an error, not NA distances", {
    # LDdecay.scope = genome_wide writes no chromosome rows; an empty map used to make
    # every region distance -Inf/NA downstream (regions.R took max(integer(0))).
    path <- write_ld(ld_row())
    expect_error(suppressMessages(resolve_clumping_distance("auto_per_chromosome", path)),
                 "no per-chromosome rows")
})

# --- resolve_clumping_distance_strict / format_clumping_distance -----------
# The association scripts (find_sig_snps.R, combine_selected_snps.R) resolve the
# shared distance with the strict variant: no silent 1e6 bp fallback.

test_that("the strict resolver passes numbers and auto specs through unchanged", {
    expect_identical(suppressMessages(resolve_clumping_distance_strict("50000")), 50000L)
    path <- write_ld(ld_row())
    expect_identical(suppressMessages(resolve_clumping_distance_strict("auto_genome_wide", path, "half_decay")),
                     suppressMessages(resolve_clumping_distance("auto_genome_wide", path, "half_decay")))
})

test_that("the strict resolver refuses a spec that is neither a number nor auto_*", {
    expect_error(resolve_clumping_distance_strict("wide"), "neither a number")
})

test_that("format_clumping_distance describes a scalar and a per-chromosome map", {
    expect_identical(format_clumping_distance(500L), "500 bp")
    expect_match(format_clumping_distance(c("1" = 100L, "2" = 300L)),
                 "per chromosome, 100-300 bp", fixed = TRUE)
})
