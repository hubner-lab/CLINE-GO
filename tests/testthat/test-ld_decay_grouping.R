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
                         "genome_wide", lists, chrv, cl, "NULL"))
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
                         "NULL", "NULL"))
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
                         file.path(d, "sample_lists"), file.path(d, "chr_vcfs"), cl, "NULL"))
    expect_gt(res$status, 0L)
    expect_match(res$output, "No Q-matrix columns")
})

# --- 600acb: n_samples no longer enters any fit ---------------------------------
# The Hill-Weir curve needed the group size and fell back to LOESS when it was 0.
# The LOESS estimator reads only the binned r^2, so n_samples is reported, not used.

analyze_with_n <- function(d, n_samples) {
    gw <- file.path(d, "stat_gw"); chrd <- file.path(d, "stat_chr")
    dir.create(gw, showWarnings = FALSE); dir.create(chrd, showWarnings = FALSE)
    fx_stat_gz(gw, "All")
    mf <- write_tsv(data.table::data.table(group = "All", n_samples = n_samples,
                                           list_file = "All"),
                    file.path(d, "manifest.tsv"))
    tbl <- file.path(d, "ld_decay_half_distances.tsv")
    res <- run_wrapper("ld_decay_analyze.R",
                       c("genome_wide", gw, chrd, mf,
                         file.path(d, "chromosomes.txt"), tbl,
                         file.path(d, "ld_decay_gw.png"), "NULL",
                         "NULL", "0.2", "NULL", "0.3"))
    list(res = res, table = data.table::fread(tbl))
}

test_that("the decay distances do not depend on the reported n_samples", {
    d0 <- withr::local_tempdir(); d50 <- withr::local_tempdir()
    r0 <- analyze_with_n(d0, 0L); r50 <- analyze_with_n(d50, 50L)

    expect_identical(r0$res$status, 0L, info = r0$res$output)
    expect_identical(r0$table$decay_method, "loess_log")
    expect_identical(r0$table$n_samples, 0L)
    expect_identical(r50$table$n_samples, 50L)
    expect_identical(r0$table$decay_half_bp, r50$table$decay_half_bp)
    expect_identical(r0$table$decay_r2thr_bp, r50$table$decay_r2thr_bp)
})

test_that("group names that are not file-safe get unique stems in the manifest", {
    d <- withr::local_tempdir()
    samples <- sprintf("ID%03d", 1:12)
    # "Tel Aviv" and "Tel_Aviv" sanitise to the same name, and a site named "All"
    # must not overwrite the whole-dataset list.
    sites <- rep(c("Tel Aviv", "Tel_Aviv", "All"), each = 4)
    sites[12] <- "all"   # one sample: below min_samples elsewhere, here a group of 1
    meta <- write_tsv(data.table::data.table(site = sites, sample = samples,
                                             latitude = 31, longitude = 34.8),
                      file.path(d, "metadata.tsv"))
    lists <- file.path(d, "sample_lists")
    res <- run_wrapper("ld_decay_prepare.R",
                       c("NULL", meta, "site", "1", "genome_wide", lists,
                         file.path(d, "chr_vcfs"), "NULL", "NULL"))
    expect_identical(res$status, 0L, info = res$output)
    mf <- data.table::fread(file.path(lists, "manifest.tsv"), colClasses = "character")
    # the site called "All" is renamed: "All" means the whole dataset everywhere downstream
    expect_identical(mf$group, c("All", "Tel Aviv", "Tel_Aviv", "All_site", "all_site"))
    expect_match(res$output, "reported as 'All_site'")
    expect_match(res$output, "reported as 'all_site'")
    # stems differ even ignoring case: a case-insensitive file system merges the rest
    expect_false(anyDuplicated(tolower(mf$list_file)) > 0)
    expect_identical(mf$list_file[1], "All")
    # every stem has its own sample list, and All.txt still holds every sample
    for (i in seq_len(nrow(mf))) {
        expect_true(file.exists(file.path(lists, paste0(mf$list_file[i], ".txt"))))
    }
    expect_identical(sort(readLines(file.path(lists, "All.txt"))), samples)
    expect_identical(readLines(file.path(lists, paste0(mf$list_file[2], ".txt"))), samples[1:4])
})

test_that("chromosomes.txt leaves out LDdecay.exclude_chromosomes", {
    d <- withr::local_tempdir()
    meta <- fx_metadata(d)
    vcf <- file.path(d, "tiny.vcf")
    row <- function(chr) paste(c(chr, "100", ".", "A", "G", ".", "PASS", "."), collapse = "\t")
    writeLines(c("##fileformat=VCFv4.2",
                 paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
                       collapse = "\t"),
                 row("1"), row("2"), row("MT"), row("Pt")), vcf)
    chrv <- file.path(d, "chr_vcfs")
    res <- run_wrapper("ld_decay_prepare.R",
                       c(vcf, meta, "site", "1", "per_chromosome",
                         file.path(d, "sample_lists"), chrv, "NULL", "chrM,Pt,ChrC"))
    expect_identical(res$status, 0L, info = res$output)
    expect_identical(readLines(file.path(chrv, "chromosomes.txt")), c("1", "2"))
    expect_match(res$output, "names no contig of the VCF: ChrC")
})
