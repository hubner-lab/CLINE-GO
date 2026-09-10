# scripts/R/utils/io_selected_snps.R, io_pvalues.R, and the parsing half of
# manhattan_utils.R.
#
# The recurring hazard here is fread()'s type inference: a chromosome named "1" or
# a sample named "88" becomes an integer and then silently fails every %in% and
# left_join against the character VCF header. Both readers pin colClasses for
# exactly that reason.

# --- io_selected_snps round-trip -------------------------------------------

test_that("selected_snps survives a write/read round-trip with chr and SNPID as character", {
    dir <- withr::local_tempdir()
    path <- file.path(dir, "selected.tsv")

    dt <- data.table::data.table(
        SNPID = c("1:100", "1:200"), chr = c("1", "1"), pos = c(100L, 200L),
        pvalue = c(1.234567e-09, 5e-3), pval_threshold = 1e-6,
        method = "EMMAX", trait = "bio_1"
    )
    write_selected_snps(dt, path)
    back <- read_selected_snps(path)

    expect_true(is.character(back$chr))     # would be integer without colClasses
    expect_true(is.character(back$SNPID))
    expect_identical(back$chr, dt$chr)
    expect_identical(back$SNPID, dt$SNPID)
    expect_identical(back$pos, dt$pos)
})

test_that("p-value precision survives the round-trip", {
    dir  <- withr::local_tempdir()
    path <- file.path(dir, "selected.tsv")
    dt <- data.table::data.table(SNPID = "1:100", chr = "1", pos = 100L,
                                 pvalue = 1.2345678901e-12)
    write_selected_snps(dt, path)
    expect_equal(read_selected_snps(path)$pvalue, 1.2345678901e-12, tolerance = 1e-20)
})

test_that("a purely numeric chromosome name is still read back as character", {
    dir  <- withr::local_tempdir()
    path <- file.path(dir, "selected.tsv")
    dt <- data.table::data.table(SNPID = c("88:1", "108:1"), chr = c("88", "108"),
                                 pos = c(1L, 1L))
    write_selected_snps(dt, path)
    back <- read_selected_snps(path)
    expect_identical(back$chr, c("88", "108"))
    expect_true(all(back$chr %in% c("88", "108")))   # the join that would break
})

test_that("write_selected_snps does not quote its output", {
    dir  <- withr::local_tempdir()
    path <- file.path(dir, "selected.tsv")
    write_selected_snps(data.table::data.table(SNPID = "1:100", chr = "1",
                                               trait = "bio_1,bio_2"), path)
    expect_false(any(grepl('"', readLines(path))))
})

test_that("an empty table round-trips as an empty table", {
    dir  <- withr::local_tempdir()
    path <- file.path(dir, "selected.tsv")
    dt <- data.table::data.table(SNPID = character(), chr = character(),
                                 pos = integer())
    write_selected_snps(dt, path)
    expect_identical(nrow(read_selected_snps(path)), 0L)
})

# --- read_pvalues_tsv ------------------------------------------------------

test_that("read_pvalues_tsv forces chr to character", {
    dir  <- withr::local_tempdir()
    path <- file.path(dir, "pvals.tsv")
    data.table::fwrite(data.table::data.table(SNPID = "1:100", chr = 1L, pos = 100L,
                                              bio_1 = 1e-9), path, sep = "\t")
    expect_true(is.character(read_pvalues_tsv(path)$chr))
})

test_that("write_pvalues and write_qvalues both emit unquoted TSV", {
    dir <- withr::local_tempdir()
    dt  <- data.table::data.table(SNPID = "1:100", chr = "1", pos = 100L, bio_1 = 1e-9)
    write_pvalues(dt, file.path(dir, "p.tsv"))
    write_qvalues(dt, file.path(dir, "q.tsv"))
    for (f in c("p.tsv", "q.tsv")) {
        lines <- readLines(file.path(dir, f))
        expect_false(any(grepl('"', lines)))
        expect_identical(strsplit(lines[1], "\t")[[1]],
                         c("SNPID", "chr", "pos", "bio_1"))
    }
})

# --- parse_assoc_files_str -------------------------------------------------

test_that("parse_assoc_files_str splits METHOD:ADJUST:PATH into a named list", {
    out <- parse_assoc_files_str("EMMAX:bonf_0.05:/a/e.tsv,LFMM:qval_0.1:/a/l.tsv")
    expect_named(out, c("EMMAX", "LFMM"))
    expect_identical(out$EMMAX$adjust, "bonf_0.05")
    expect_identical(out$LFMM$filepath, "/a/l.tsv")
    expect_identical(out$LFMM$method, "LFMM")
})

test_that("parse_assoc_files_str SILENTLY drops any item that is not 3 colon-parts", {
    # Documented contract, and the reason a path containing a colon cannot be
    # passed through this format at all.
    out <- parse_assoc_files_str("EMMAX:bonf_0.05:/a/e.tsv,BROKEN:only_two")
    expect_named(out, "EMMAX")

    expect_length(parse_assoc_files_str("A:b:c:d"), 0L)
    expect_length(parse_assoc_files_str(""), 0L)
})

test_that("a later entry for the same method replaces the earlier one", {
    out <- parse_assoc_files_str("EMMAX:bonf_0.05:/first.tsv,EMMAX:top_10:/second.tsv")
    expect_length(out, 1L)
    expect_identical(out$EMMAX$filepath, "/second.tsv")
})

# --- load_assoc_data -------------------------------------------------------

write_pval_file <- function(dir, name, dt) {
    path <- file.path(dir, name)
    data.table::fwrite(dt, path, sep = "\t")
    path
}

test_that("load_assoc_data returns tidy per-trait rows and -log10 thresholds", {
    dir <- withr::local_tempdir()
    dt  <- data.table::data.table(
        SNPID = paste0("1:", 1:5), chr = "1", pos = 1:5,
        bio_1 = c(1e-9, 1e-8, 0.5, 0.6, 0.7),
        bio_2 = c(0.5, 0.5, 0.5, 0.5, 1e-9)
    )
    path <- write_pval_file(dir, "emmax.tsv", dt)
    info <- parse_assoc_files_str(paste0("EMMAX:custom_1e-6:", path))

    res <- quiet(load_assoc_data(info, c("bio_1", "bio_2")))

    expect_identical(nrow(res$data), 10L)               # 5 SNPs x 2 traits
    expect_setequal(unique(res$data$trait), c("bio_1", "bio_2"))
    expect_identical(unique(res$data$method), "EMMAX")
    expect_true(all(c("log10p", "is_significant") %in% colnames(res$data)))

    # Thresholds are keyed "{method}_{trait}" and held as -log10, not raw p.
    expect_named(res$thresholds, c("EMMAX_bio_1", "EMMAX_bio_2"))
    expect_equal(res$thresholds$EMMAX_bio_1, -log10(1e-6))
})

test_that("significance uses the same INCLUSIVE boundary as sig_snps.R", {
    # log10p >= threshold_log10 is the plot-side spelling of p <= threshold, so a
    # SNP sitting exactly on the cutoff is significant in the table AND the plot.
    dir  <- withr::local_tempdir()
    dt   <- data.table::data.table(SNPID = c("1:1", "1:2"), chr = "1", pos = 1:2,
                                   bio_1 = c(1e-6, 0.5))
    path <- write_pval_file(dir, "emmax.tsv", dt)
    info <- parse_assoc_files_str(paste0("EMMAX:custom_1e-6:", path))

    res <- quiet(load_assoc_data(info, "bio_1"))
    expect_true(res$data$is_significant[res$data$SNPID == "1:1"])
    expect_false(res$data$is_significant[res$data$SNPID == "1:2"])
})

test_that("a trait missing from the file is warned about and skipped", {
    dir  <- withr::local_tempdir()
    dt   <- data.table::data.table(SNPID = "1:1", chr = "1", pos = 1L, bio_1 = 1e-9)
    path <- write_pval_file(dir, "emmax.tsv", dt)
    info <- parse_assoc_files_str(paste0("EMMAX:custom_1e-6:", path))

    expect_message(load_assoc_data(info, c("bio_1", "bio_99")), "bio_99 not found")
    res <- quiet(load_assoc_data(info, c("bio_1", "bio_99")))
    expect_named(res$thresholds, "EMMAX_bio_1")
})

test_that("a failed threshold keeps the data but flags nothing significant", {
    dir  <- withr::local_tempdir()
    # qval needs >= 10 tests; with 3 the status is too_few_tests.
    dt   <- data.table::data.table(SNPID = paste0("1:", 1:3), chr = "1", pos = 1:3,
                                   bio_1 = c(1e-9, 1e-9, 1e-9))
    path <- write_pval_file(dir, "emmax.tsv", dt)
    info <- parse_assoc_files_str(paste0("EMMAX:qval_0.1:", path))

    res <- quiet(load_assoc_data(info, "bio_1"))
    expect_identical(nrow(res$data), 3L)                # data retained for plotting
    expect_false(any(res$data$is_significant))
    expect_true(is.na(res$thresholds$EMMAX_bio_1))
})

test_that("panel_label adds a panel column, used by the Miami plot", {
    dir  <- withr::local_tempdir()
    dt   <- data.table::data.table(SNPID = "1:1", chr = "1", pos = 1L, bio_1 = 1e-9)
    path <- write_pval_file(dir, "emmax.tsv", dt)
    info <- parse_assoc_files_str(paste0("EMMAX:custom_1e-6:", path))

    plain <- quiet(load_assoc_data(info, "bio_1"))
    panel <- quiet(load_assoc_data(info, "bio_1", panel_label = "GWAS"))
    expect_false("panel" %in% colnames(plain$data))
    expect_identical(unique(panel$data$panel), "GWAS")
})
