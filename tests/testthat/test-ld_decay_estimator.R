# Unit tests for scripts/R/lib/ld_decay.R: the LOESS-based LD-decay distances
# measured against a background level, and the PLINK pair aggregation that feeds them.
# Synthetic per-bp tables with a known decay, r2(d) = bg + (r0 - bg) * exp(-d / lambda),
# so every crossing has an analytic value.

synthetic_stat <- function(r0 = 0.3, bg = 0.02, lambda = 500, max_dist = 50000,
                           n_pairs = 100L) {
    d  <- seq_len(max_dist)
    r2 <- bg + (r0 - bg) * exp(-d / lambda)
    data.table::data.table(dist_bp = d, sum_r2 = r2 * n_pairs, n_pairs = n_pairs)
}
crossing_true <- function(target, r0 = 0.3, bg = 0.02, lambda = 500) {
    -lambda * log((target - bg) / (r0 - bg))
}

test_that("ld_log_breaks: 10-bp bins below 100 bp, then increasing to max_dist", {
    b <- ld_log_breaks(50000)
    expect_equal(b[1:11], seq(0, 100, 10))
    expect_true(all(diff(b) > 0))
    expect_gte(max(b), 50000)
    expect_equal(ld_log_breaks(60), c(seq(0, 60, 10)))
})

test_that("ld_bin_stat: pair-weighted mean distance and mean r2 per bin", {
    raw <- data.table::data.table(dist_bp = c(5, 6, 15), sum_r2 = c(2, 6, 3),
                                  n_pairs = c(10L, 30L, 10L))
    out <- ld_bin_stat(raw, c(0, 10, 20))
    expect_equal(nrow(out), 2L)
    expect_equal(out$dist_bp[1], (5 * 10 + 6 * 30) / 40)
    expect_equal(out$mean_r2[1], 8 / 40)
    expect_equal(out$n_pairs, c(40L, 10L))
    expect_null(ld_bin_stat(raw[0], c(0, 10)))
})

test_that("ld_critical_r2: quantile of sqrt(r2), squared, non-finite values dropped", {
    r2 <- c(seq(0, 0.5, length.out = 101), NA, Inf)
    expect_equal(ld_critical_r2(r2, 0.95),
                 stats::quantile(sqrt(seq(0, 0.5, length.out = 101)), 0.95, names = FALSE)^2)
    expect_true(is.na(ld_critical_r2(c(NA_real_, NaN))))
})

test_that("ld_decay_distances recovers the analytic crossings of a known decay", {
    raw <- synthetic_stat()
    binned <- ld_bin_stat(raw, ld_log_breaks(max(raw$dist_bp)))
    res <- ld_decay_distances(binned, background_r2 = 0.03, r2_threshold = 0.2)
    expect_equal(res$method, "loess_log")
    s <- res$short_r2
    expect_equal(s, 0.3, tolerance = 0.03)
    expect_equal(res$half_bp,     crossing_true(s / 2),            tolerance = 0.05)
    expect_equal(res$critical_bp, crossing_true(0.03),             tolerance = 0.05)
    expect_equal(res$r2thr_bp,    crossing_true(0.2),              tolerance = 0.05)
    # the background lies below half of the short-range LD, so it is reached last
    expect_lt(res$half_bp, res$critical_bp)
    expect_null(res$half_bg_bp)   # the halfway rule is not a published one and is gone
})

test_that("crossings are NA when the target is never reached or the curve starts below it", {
    raw <- synthetic_stat(r0 = 0.15, bg = 0.05)
    binned <- ld_bin_stat(raw, ld_log_breaks(max(raw$dist_bp)))
    res <- ld_decay_distances(binned, background_r2 = 0.02, r2_threshold = 0.2)
    expect_true(is.na(res$r2thr_bp))      # starts below 0.2
    expect_true(is.na(res$critical_bp))   # levels off at 0.05, never reaches 0.02
    expect_false(is.na(res$half_bp))      # half of 0.15 lies above the 0.05 plateau
    # half of the start at or below the plateau is never reached either
    flat <- ld_bin_stat(synthetic_stat(r0 = 0.1, bg = 0.05),
                        ld_log_breaks(50000))
    expect_true(is.na(ld_decay_distances(flat)$half_bp))
    res_hi <- ld_decay_distances(binned, background_r2 = 0.5)
    expect_true(is.na(res_hi$critical_bp)) # background above the short-range LD
})

test_that("ld_fit_curve keeps the smoothed curve within [0, 1]", {
    # Sparse, noisy bins make a degree-2 LOESS overshoot; r^2 cannot.
    set.seed(2)
    binned <- data.table::data.table(dist_bp = 10^seq(2, 4.7, length.out = 12),
                                     mean_r2 = c(0.9, 0, 1, 0.05, 1, 0, 0.95, 0, 1, 0.02, 0.9, 0.01),
                                     n_pairs = 3L)
    cv <- ld_fit_curve(binned)
    expect_true(all(cv$r2 >= 0 & cv$r2 <= 1))
    expect_true(all(diff(cv$r2_mono) <= 0))
})

# --- PLINK pairs -> per-distance sums ------------------------------------------

plink_pairs <- function() {
    data.table::data.table(
        CHR_A = c("1", "1", "1", "1", "2", "2", "1"),
        BP_A  = c(100L, 100L, 300L, 100L, 50L, 50L, 100L),
        SNP_A = c("v1", "v1", "v3", "v1", "v8", "v8", "v1"),
        BP_B  = c(300L, 100L, 100L, 900L, 60L, 90L, 600L),
        SNP_B = c("v3", "v2", "v1", "v9", "v9b", "v10", "v5"),
        R2    = c("0.5", "0.9", "0.5", "0.2", "nan", "0.4", "0.1"))
}

test_that("ld_pairs_to_stat sums r2 per distance, dropping distance 0, NaN and far pairs", {
    st <- ld_pairs_to_stat(plink_pairs(), max_dist = 600, dedupe = FALSE)
    # distance 0 (v1-v2), NaN (v8-v9b) and 800 bp (> max_dist) are dropped;
    # v1-v3 and v3-v1 both stay without dedupe (an exhaustive PLINK run never repeats one)
    expect_equal(st$gw$dist_bp, c(40L, 200L, 500L))
    expect_equal(st$gw$n_pairs, c(1L, 2L, 1L))
    expect_equal(st$gw$sum_r2, c(0.4, 1.0, 0.1))
    expect_equal(sort(unique(st$chr$chr)), c("1", "2"))
    expect_equal(st$chr[chr == "2", sum(n_pairs)], 1L)
})

test_that("ld_pairs_to_stat with dedupe counts an anchor-anchor pair once", {
    st <- ld_pairs_to_stat(plink_pairs(), max_dist = 600, dedupe = TRUE)
    expect_equal(st$gw[dist_bp == 200, n_pairs], 1L)   # v1-v3 reported from both anchors
    expect_equal(st$gw[dist_bp == 200, sum_r2], 0.5)
})

test_that("ld_anchor_plan is exhaustive under the budget and samples above it", {
    chr <- rep(c("1", "2"), each = 100)
    pos <- c(seq(100, 10000, by = 100), seq(100, 10000, by = 100))
    # within 1 kb every SNP has up to 10 partners on either side
    all_pairs <- ld_anchor_plan(chr, pos, max_dist = 1000, max_pairs = 1e6)
    expect_true(all_pairs$exhaustive)
    expect_equal(all_pairs$n_anchors, 200L)
    expect_equal(all_pairs$total_pairs, 2 * (100 * 10 - sum(1:10)))
    sampled <- ld_anchor_plan(chr, pos, max_dist = 1000, max_pairs = 500)
    expect_false(sampled$exhaustive)
    expect_equal(sampled$n_anchors, as.integer(ceiling(500 / sampled$mean_partners)))
    # SNPs sharing a position are not each other's partners
    dup <- ld_anchor_plan(c("1", "1", "1"), c(10, 10, 20), max_dist = 100, max_pairs = 10)
    expect_equal(dup$total_pairs, 2)
})

test_that("ld_write_stat writes PopLDdecay's six-column layout that the reader takes", {
    d <- withr::local_tempdir()
    st <- ld_pairs_to_stat(plink_pairs(), max_dist = 600)
    f <- file.path(d, "All.stat.gz")
    ld_write_stat(st$gw, f)
    back <- data.table::fread(cmd = paste("zcat", shQuote(f)), header = TRUE)
    expect_equal(names(back), c("#Dist", "Mean_r^2", "Mean_D'", "Sum_r^2", "Sum_D'", "NumberPairs"))
    expect_equal(back$`Sum_r^2`, st$gw$sum_r2)
    expect_equal(back$NumberPairs, st$gw$n_pairs)
    expect_equal(back$`Mean_r^2`, st$gw$sum_r2 / st$gw$n_pairs)
})

test_that("ld_background_markers keeps per_chr per chromosome and never exceeds the cap", {
    few <- ld_background_markers(rep(c("1", "2", "3"), each = 100), rep(1:100, 3), 50, 2000, 1)
    expect_equal(few$per_chr, 50)
    expect_identical(as.vector(table(few$markers$chr)), c(50L, 50L, 50L))
    again <- ld_background_markers(rep(c("1", "2", "3"), each = 100), rep(1:100, 3), 50, 2000, 1)
    expect_equal(few$markers, again$markers)   # seeded
    # a fragmented assembly: 3,000 contigs of 2 SNPs give one marker on each of 2,000
    frag <- ld_background_markers(rep(sprintf("ctg%04d", 1:3000), each = 2),
                                  rep(c(10L, 20L), 3000), 50, 2000, 1)
    expect_equal(frag$per_chr, 1)
    expect_identical(nrow(frag$markers), 2000L)
    expect_identical(data.table::uniqueN(frag$markers$chr), 2000L)
})

test_that("ld_scratch_dir makes a fresh folder beside a file or a folder", {
    d <- withr::local_tempdir()
    s <- ld_scratch_dir(file.path(d, "out.tsv"), "x_tmp")
    expect_identical(normalizePath(s), normalizePath(file.path(d, "x_tmp")))
    writeLines("left over", file.path(s, "old.txt"))
    s2 <- ld_scratch_dir(file.path(d, "out.tsv"), "x_tmp")   # a failed run's leftover goes
    expect_identical(list.files(s2), character(0))
    dir.create(file.path(d, "stat_gw"))
    s3 <- ld_scratch_dir(paste0(file.path(d, "stat_gw"), "/"), "y_tmp")   # a folder: its parent
    expect_identical(normalizePath(s3), normalizePath(file.path(d, "y_tmp")))
})

test_that("too few bins give insufficient_data, not an error", {
    raw <- synthetic_stat(max_dist = 40)
    binned <- ld_bin_stat(raw, ld_log_breaks(40))
    res <- ld_decay_distances(binned, background_r2 = 0.03)
    expect_equal(res$method, "insufficient_data")
    expect_true(is.na(res$half_bp))
})

# CONTENT check of the CLI integration (the wrapper row in test-cli-wrappers.R only
# asserts exit 0 and non-empty outputs): ld_decay_analyze.R with a background table
# writes the background-aware columns and the curve table.
source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

test_that("ld_decay_analyze.R writes background-aware distances and curves", {
    d <- withr::local_tempdir()
    gw <- file.path(d, "stat_gw"); dir.create(gw)
    fx_stat_gz(gw, "All")      # 0.6 * exp(-d / 20000) + 0.02, at 1-60 kb
    mf <- write_tsv(data.table::data.table(group = "All", n_samples = 50L, list_file = "All"),
                    file.path(d, "manifest.tsv"))
    bg <- write_tsv(data.table::data.table(group = "All", critical_r2 = 0.1, mean_r2 = 0.04),
                    file.path(d, "ld_decay_background.tsv"))
    out <- file.path(d, "half.tsv"); curve <- file.path(d, "curve.tsv")
    res <- run_wrapper("ld_decay_analyze.R",
                       c("genome_wide", gw, file.path(d, "none"), mf,
                         file.path(d, "chromosomes.txt"), out, "NULL", "NULL",
                         bg, "0.2", curve, "0.3"), wd = d)
    expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
    tab <- data.table::fread(out)
    expect_identical(names(tab), c("group", "scope", "n_samples", "n_pairs", "max_dist_bp",
                                   "short_r2", "background_r2", "background_mean_r2",
                                   "decay_critical_bp", "decay_half_bp", "decay_r2thr_bp",
                                   "r2_threshold", "decay_method"))
    expect_equal(tab$decay_method, "loess_log")
    expect_equal(tab$background_r2, 0.1)
    expect_equal(tab$background_mean_r2, 0.04)
    # analytic crossing of 0.1 on 0.6 * exp(-d / 20000) + 0.02: d = 20000 * log(0.6 / 0.08)
    expect_equal(tab$decay_critical_bp, 20000 * log(0.6 / 0.08), tolerance = 0.08)
    cv <- data.table::fread(curve)
    expect_gt(nrow(cv), 100L)
    expect_true(all(diff(cv$r2_mono) <= 0))
    expect_false(dir.exists(file.path(d, "analyze_tmp")))   # scratch removed, never in tempdir()
})

test_that("ld_excluded_contigs matches LDdecay.exclude_chromosomes to the VCF's contigs", {
    chrs <- c("1", "2", "X", "MT", "Pt", "Y")
    expect_identical(ld_excluded_contigs(chrs, "NULL"), character(0))
    expect_identical(ld_excluded_contigs(chrs, ""), character(0))
    expect_identical(ld_excluded_contigs(chrs, "X,MT"), c("X", "MT"))
    # case, a chr prefix and M for MT are ignored, as PLINK reads them; VCF order is kept
    expect_identical(ld_excluded_contigs(chrs, "chrM, pt,x"), c("X", "MT", "Pt"))
    expect_identical(ld_excluded_contigs(c("chr1", "chrM", "chr1"), "M"), "chrM")
    # a listed name with no contig is reported, never an error
    expect_message(out <- ld_excluded_contigs(chrs, "X,ChrC"),
                   "names no contig of the VCF: ChrC")
    expect_identical(out, "X")
})
