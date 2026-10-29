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
                   half = 1240, r2thr = 2600, r2_threshold = 0.2, status = "reliable") {
    data.table::data.table(
        group = group, scope = scope, n_samples = 50L, n_pairs = 1e5,
        max_dist_bp = 5e5, short_r2 = 0.45, background_r2 = 0.06,
        background_mean_r2 = 0.02, decay_critical_bp = critical,
        decay_half_bp = half, decay_r2thr_bp = r2thr, r2_threshold = r2_threshold,
        decay_method = "loess_log", ld_status = status,
        ld_status_reason = if (status == "unreliable") "short-range r2 0.45 is 1.50x the critical r2" else ""
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

test_that("auto_per_chromosome hands its genome-wide value on for chromosomes with no row", {
    # f503bd: a small contig with no anchor pair has no LD-decay row at all; regions.R's
    # distance_for_chr() gave it the map's maximum (4000 here), not the genome-wide 3000.
    path <- write_ld(rbind(
        ld_row(scope = "genome_wide", critical = 3000),
        ld_row(scope = "1", critical = 4000),
        ld_row(scope = "2", critical = 2000)
    ))
    d <- suppressMessages(resolve_clumping_distance("auto_per_chromosome", path))
    expect_identical(attr(d, "genome_wide"), 3000L)
    expect_identical(suppressMessages(distance_for_chr(d, c("1", "contig_9"))), c(4000L, 3000L))
    # Through the cross-method agreement, which coerces the spec first.
    expect_identical(attr(as_distance_spec(d), "genome_wide"), 3000L)

    no_gw <- write_ld(rbind(ld_row(scope = "1", critical = 4000),
                            ld_row(scope = "2", critical = 2000)))
    d2 <- suppressMessages(resolve_clumping_distance("auto_per_chromosome", no_gw))
    expect_null(attr(d2, "genome_wide"))
    expect_identical(suppressMessages(distance_for_chr(d2, "contig_9")), 4000L)
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

# --- reliability (ld_status) -----------------------------------------------
# An auto_* distance refuses a curve ld_decay_analyze.R marked unreliable, whatever
# rule is asked: the curve meets the unlinked background before it halves, so no
# distance read off it measures decay.

test_that("an unreliable genome-wide curve stops auto_genome_wide under every rule", {
    path <- write_ld(ld_row(status = "unreliable"))
    for (rule in c("critical", "half_decay", "r2_threshold")) {
        expect_error(suppressMessages(resolve_clumping_distance("auto_genome_wide", path, rule)),
                     "is unreliable: short-range r2 0.45 is 1.50x.*declare snp_clumping_distance")
    }
    # a declared distance is untouched by the table
    expect_identical(suppressMessages(resolve_clumping_distance("150000", path)), 150000L)
})

test_that("an unreliable chromosome falls back to a reliable genome-wide distance", {
    path <- write_ld(rbind(
        ld_row(scope = "genome_wide", critical = 3000),
        ld_row(scope = "1", critical = 4000),
        ld_row(scope = "2", critical = 900, status = "unreliable")
    ))
    msgs <- character(0)
    d <- withCallingHandlers(resolve_clumping_distance("auto_per_chromosome", path),
                             message = function(m) {
                                 msgs <<- c(msgs, conditionMessage(m))
                                 invokeRestart("muffleMessage")
                             })
    expect_identical(d[["1"]], 4000L)
    expect_identical(d[["2"]], 3000L)    # not its own 900 bp
    expect_true(any(grepl("chr 2: LD decay unreliable", msgs)))
})

test_that("an unreliable chromosome with an unreliable genome-wide curve stops", {
    path <- write_ld(rbind(
        ld_row(scope = "genome_wide", critical = 3000, status = "unreliable"),
        ld_row(scope = "1", critical = 4000),
        ld_row(scope = "2", critical = 900, status = "unreliable")
    ))
    expect_error(suppressMessages(resolve_clumping_distance("auto_per_chromosome", path)),
                 "chromosome 2, with no reliable genome-wide value")
})

# A table written before the check (no ld_status) used to read as reliable, whatever its
# curve (96e4fd). Its status is now derived from the row's own short_r2 and background_r2
# under ld_reliability(), the rule ld_decay_analyze.R writes it with.

test_that("a table without ld_status: a reliable curve resolves, saying its status was derived", {
    path <- write_ld(ld_row()[, !c("ld_status", "ld_status_reason")])   # 0.45 / 0.06 = 7.5x
    msgs <- character(0)
    d <- withCallingHandlers(resolve_clumping_distance("auto_genome_wide", path),
                             message = function(m) {
                                 msgs <<- c(msgs, conditionMessage(m))
                                 invokeRestart("muffleMessage")
                             })
    expect_identical(d, 3100L)
    expect_true(any(grepl("predates the reliability check.*derived", msgs)))
})

test_that("a table without ld_status: an unreliable curve stops as a written status does", {
    legacy <- ld_row()[, !c("ld_status", "ld_status_reason")]
    legacy[, background_r2 := 0.3]                                       # 0.45 / 0.3 = 1.5x
    path <- write_ld(legacy)
    expect_error(suppressMessages(resolve_clumping_distance("auto_genome_wide", path)),
                 "is unreliable: short-range r2 0.45 is 1.50x the critical r2")
    # per chromosome: an unreliable chromosome falls back to a reliable genome-wide value
    path <- write_ld(rbind(
        ld_row(scope = "genome_wide", critical = 3000),
        ld_row(scope = "2", critical = 900)[, background_r2 := 0.3]
    )[, !c("ld_status", "ld_status_reason")])
    d <- suppressMessages(resolve_clumping_distance("auto_per_chromosome", path))
    expect_identical(d[["2"]], 3000L)
})

test_that("ld_reliability: each branch", {
    expect_identical(ld_reliability(0.5, 0.05, "loess_log")$status, "reliable")
    expect_identical(ld_reliability(0.5, 0.05, "loess_log")$reason, "")
    r <- ld_reliability(0.464, 0.26936, "loess_log")          # Trifolium, ratio 1.72
    expect_identical(r$status, "unreliable")
    expect_match(r$reason, "1\\.72x")
    expect_identical(ld_reliability(0.4, 0.2, "loess_log")$status, "reliable")   # ratio = 2
    expect_identical(ld_reliability(0.4, 0.2, "loess_log", min_ratio = 3)$status, "unreliable")
    expect_match(ld_reliability(NA, 0.1, "insufficient_data")$reason, "too few distance bins")
    expect_match(ld_reliability(0.5, NA, "loess_log")$reason, "one chromosome")
    expect_identical(ld_reliability(0.8, 1, "loess_log")$status, "unreliable")   # starts below
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
