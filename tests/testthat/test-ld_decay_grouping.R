# LD-decay grouping and curve fitting — CONTENT assertions, quick tier.
#
# test-cli-wrappers.R asserts only exit status and non-empty outputs, which both
# of the defects below pass: ld_decay_prepare.R wrote an EMPTY All.txt and a
# manifest saying n_samples=0 at exit 0, and ld_decay_analyze.R then reported
# every curve as a LOESS fallback without saying why. Findings ef93c6 / 600acb.
#
# These run the two scripts as subprocesses (SCOPE="genome_wide", so no VCF is
# read) and read what they wrote.

source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

# A clusters_K{k}.tsv exactly as extract_clusters.R writes it: sample, site,
# then the K Q-matrix columns. The `site` column is the one that used to be
# dragged into max.col().
fx_clusters <- function(d, samples, sites, k = 3L) {
    q <- matrix(0.1, nrow = length(samples), ncol = k)
    # Give each sample an unambiguous winner, cycling through the K clusters.
    for (i in seq_along(samples)) q[i, ((i - 1L) %% k) + 1L] <- 0.8
    dt <- data.table::data.table(sample = samples, site = sites)
    for (j in seq_len(k)) dt[[paste0("C", j)]] <- q[, j]
    write_tsv(dt, file.path(d, paste0("clusters_K", k, ".tsv")))
}

prepare_cluster_run <- function(d, n = 12L, k = 3L, min_samples = 1L) {
    samples <- sprintf("ID%03d", seq_len(n))
    sites   <- rep(c("NEG", "TAV", "GAL"), length.out = n)
    meta <- write_tsv(
        data.table::data.table(site = sites, sample = samples,
                               latitude = 30 + seq_len(n) * 0.1,
                               longitude = 34 + seq_len(n) * 0.1),
        file.path(d, "metadata.tsv"))
    cl   <- fx_clusters(d, samples, sites, k = k)
    lists <- file.path(d, "sample_lists"); chrv <- file.path(d, "chr_vcfs")
    res <- run_wrapper("ld_decay_prepare.R",
                       c("NULL", meta, "cluster", as.character(min_samples),
                         "genome_wide", lists, chrv, cl))
    list(res = res, lists = lists, n = n, samples = samples)
}

test_that("group_by='cluster' assigns every sample a cluster, not NA", {
    d <- withr::local_tempdir()
    r <- prepare_cluster_run(d)

    expect_identical(r$res$status, 0L, info = r$res$output)
    # The old code coerced .SD to character via the `site` column.
    expect_false(grepl("NAs introduced by coercion", r$res$output, fixed = TRUE))
    expect_match(r$res$output, "Q-matrix columns: C1, C2, C3")

    manifest <- data.table::fread(file.path(r$lists, "manifest.tsv"))
    # 12 samples cycling over 3 clusters at min_samples=1 => All + C1..C3.
    expect_identical(sort(manifest$group), c("All", "C1", "C2", "C3"))
    expect_identical(manifest[group == "All", n_samples], r$n)
    expect_identical(sum(manifest[group != "All", n_samples]), r$n)
})

test_that("All.txt holds every sample, including any that lost a group", {
    d <- withr::local_tempdir()
    r <- prepare_cluster_run(d)

    all_txt <- readLines(file.path(r$lists, "All.txt"))
    expect_identical(sort(all_txt), sort(r$samples))
})

test_that("All.txt still covers the dataset when no group clears min_samples", {
    d <- withr::local_tempdir()
    # min_samples far above any cluster size: every per-group list is dropped,
    # but 'All' is the dataset and must survive.
    r <- prepare_cluster_run(d, min_samples = 99L)

    expect_identical(r$res$status, 0L, info = r$res$output)
    manifest <- data.table::fread(file.path(r$lists, "manifest.tsv"))
    expect_identical(manifest$group, "All")
    expect_identical(manifest[group == "All", n_samples], r$n)
    expect_length(readLines(file.path(r$lists, "All.txt")), r$n)
})

test_that("group_by must be 'site' or 'cluster'", {
    d <- withr::local_tempdir()
    meta <- fx_metadata(d)
    res <- run_wrapper("ld_decay_prepare.R",
                       c("NULL", meta, "population", "1", "genome_wide",
                         file.path(d, "sample_lists"), file.path(d, "chr_vcfs"),
                         "NULL"))
    # Previously `groups` was simply never created and the script died later with
    # "object 'groups' not found".
    expect_gt(res$status, 0L)
    expect_match(res$output, "Unsupported group_by")
})

test_that("a clusters table with no C* columns fails loudly", {
    d <- withr::local_tempdir()
    meta <- fx_metadata(d, n = 6)
    cl <- write_tsv(data.table::data.table(sample = sprintf("ID%03d", 1:6),
                                           site = "NEG",
                                           cluster1 = 0.9, cluster2 = 0.1),
                    file.path(d, "clusters_K2.tsv"))
    res <- run_wrapper("ld_decay_prepare.R",
                       c("NULL", meta, "cluster", "1", "genome_wide",
                         file.path(d, "sample_lists"), file.path(d, "chr_vcfs"), cl))
    expect_gt(res$status, 0L)
    expect_match(res$output, "No Q-matrix columns")
})

# --- 600acb: the Hill-Weir sample-size term ---------------------------------

analyze_with_n <- function(d, n_samples) {
    gw <- file.path(d, "stat_gw"); chrd <- file.path(d, "stat_chr")
    dir.create(gw, showWarnings = FALSE); dir.create(chrd, showWarnings = FALSE)
    fx_stat_gz(gw, "All")
    mf <- write_tsv(data.table::data.table(group = "All", n_samples = n_samples),
                    file.path(d, "manifest.tsv"))
    tbl <- file.path(d, "ld_decay_half_distances.tsv")
    res <- run_wrapper("ld_decay_analyze.R",
                       c("genome_wide", gw, chrd, mf,
                         file.path(d, "chromosomes.txt"), tbl,
                         file.path(d, "ld_decay_gw.png"), "NULL"))
    list(res = res, table = data.table::fread(tbl))
}

test_that("n_samples=0 is reported, not silently absorbed into a LOESS fallback", {
    d <- withr::local_tempdir()
    r <- analyze_with_n(d, 0L)

    expect_identical(r$res$status, 0L, info = r$res$output)
    expect_match(r$res$output, "n_samples=0 is not usable")
    # Still produces a usable curve, and says which model it came from.
    expect_identical(r$table$method, "loess")
})

test_that("a usable n_samples still fits Hill-Weir", {
    d <- withr::local_tempdir()
    r <- analyze_with_n(d, 50L)

    expect_identical(r$res$status, 0L, info = r$res$output)
    expect_false(grepl("is not usable", r$res$output, fixed = TRUE))
    expect_identical(r$table$method, "hill_weir")
})
