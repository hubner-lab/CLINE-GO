# HEAVY tier: CONTENT of the PLINK LD-decay producer (scripts/ld_decay_pairs.R) and the
# background (scripts/ld_decay_background.R). The wrapper rows only assert exit 0.
#
# The fixture (heavy_ld_fixture) has 26 numeric chromosomes, heterozygous calls and
# missing calls. In PLINK's default human mode chromosomes 23-26 would be X/Y/XY/MT and
# one of them would lose its pairs; the scripts run PLINK with --chr-set 95 no-x no-y
# no-xy no-mt, and every r^2 must equal R's pairwise-complete squared correlation.

read_stat <- function(path) {
    dt <- data.table::fread(cmd = paste("zcat", shQuote(path)), header = TRUE)
    data.table::setnames(dt, c("dist_bp", "mean_r2", "mean_dprime", "sum_r2", "sum_dprime", "n_pairs"))
    dt
}

test_that("ld_decay_pairs.R reproduces R's pairwise r2 per distance, on every chromosome", {
    require_tools(c("plink", "bcftools"))
    d  <- withr::local_tempdir()
    fx <- heavy_ld_fixture(d)
    # an earlier run's file for a chromosome that has no pairs now must not survive
    dir.create(file.path(d, "chr"))
    writeLines("stale", file.path(d, "chr", "All_99.stat.gz"))
    res <- run_wrapper("ld_decay_pairs.R",
                       c(fx$vcf, fx$sl, "20", "1e6", "1", "both", file.path(d, "gw"),
                         file.path(d, "chr"), file.path(d, "summary.tsv")), wd = d)
    expect_identical(res$status, 0L, info = paste(res$output, collapse = "\n"))
    expect_false(file.exists(file.path(d, "chr", "All_99.stat.gz")))

    got <- read_stat(file.path(d, "gw", "All.stat.gz"))
    exp <- heavy_ld_expected(fx, 20000)
    expect_identical(got$n_pairs, exp$n_pairs)        # every pair once, none lost
    expect_equal(got$dist_bp, exp$dist_bp)
    expect_equal(got$sum_r2, exp$sum_r2, tolerance = 1e-5)

    # the sub-group is computed on its own samples
    got_g <- read_stat(file.path(d, "gw", "G_1.stat.gz"))
    exp_g <- heavy_ld_expected(fx, 20000, rows = 1:15)
    expect_identical(sum(got_g$n_pairs), sum(exp_g$n_pairs))
    expect_equal(sum(got_g$sum_r2), sum(exp_g$sum_r2), tolerance = 1e-5)

    # per-chromosome files carry the VCF's chromosome names, all 26 of them
    chr_files <- list.files(file.path(d, "chr"), pattern = "^All_")
    expect_setequal(chr_files, paste0("All_", 1:26, ".stat.gz"))
    s <- data.table::fread(file.path(d, "summary.tsv"))
    expect_true(all(s$exhaustive))
    expect_false(dir.exists(file.path(d, "pairs_tmp")))   # PLINK's set and pair tables removed
})

test_that("a pair budget below the total samples anchors and counts each pair once", {
    require_tools(c("plink", "bcftools"))
    d  <- withr::local_tempdir()
    fx <- heavy_ld_fixture(d)
    res <- run_wrapper("ld_decay_pairs.R",
                       c(fx$vcf, fx$sl, "20", "500", "1", "genome_wide", file.path(d, "gw"),
                         file.path(d, "chr"), file.path(d, "summary.tsv")), wd = d)
    expect_identical(res$status, 0L, info = paste(res$output, collapse = "\n"))
    s <- data.table::fread(file.path(d, "summary.tsv"))
    expect_false(any(s$exhaustive))
    got <- read_stat(file.path(d, "gw", "All.stat.gz"))
    exp <- heavy_ld_expected(fx, 20000)
    m <- merge(got, exp, by = "dist_bp", suffixes = c("_got", "_all"))
    expect_identical(nrow(m), nrow(got))              # no distance without a real pair
    expect_true(all(m$n_pairs_got <= m$n_pairs_all))  # never more pairs than exist
    expect_lt(sum(got$n_pairs), sum(exp$n_pairs))
    expect_identical(length(list.files(file.path(d, "chr"))), 0L)  # genome_wide scope only
})

test_that("ld_decay_background.R gives the mean and critical r2 of unlinked pairs", {
    require_tools(c("plink", "bcftools"))
    d  <- withr::local_tempdir()
    fx <- heavy_ld_fixture(d)
    out <- file.path(d, "bg.tsv")
    res <- run_wrapper("ld_decay_background.R", c(fx$vcf, fx$sl, "15", "0.95", "1", out), wd = d)
    expect_identical(res$status, 0L, info = paste(res$output, collapse = "\n"))
    bg <- data.table::fread(out, colClasses = c(group = "character"))
    expect_identical(bg$group, c("All", "G 1"))
    # 15 markers per chromosome = every SNP: all between-chromosome pairs, in R
    X <- fx$X; chr <- fx$info$chr
    ij <- which(upper.tri(diag(ncol(X))) & outer(chr, chr, "!="), arr.ind = TRUE)
    r2 <- apply(ij, 1, function(p) suppressWarnings(stats::cor(X[, p[1]], X[, p[2]],
                                                               use = "pairwise.complete.obs")^2))
    r2 <- r2[is.finite(r2)]
    expect_identical(bg[group == "All", n_unlinked_pairs], length(r2))
    expect_equal(bg[group == "All", mean_r2], mean(r2), tolerance = 1e-4)
    expect_equal(bg[group == "All", critical_r2], ld_critical_r2(r2, 0.95), tolerance = 1e-4)
    expect_false(dir.exists(file.path(d, "background_tmp")))
})
