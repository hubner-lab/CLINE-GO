# Phase 2c (ID binding) — CONTENT assertions on the scripts, quick tier.
#
# Each script used to pair two inputs by ROW POSITION. SIMDATA cannot show the fix: every
# table it writes is in the metadata's row order, so a join by ID gives the same pairing.
# These tests break that alignment on purpose — one input permuted, or an equal row count
# with a different ID set — and require either the unpermuted output, byte for byte, or a
# stop. Under the positional code every one of them fails.
#   lfmm2vcf.R               953fac  writes the output it is given
#   gradient_forest_model.R  ec0d05  climate rows matched to `samples` by id
#   gradient_forest_offset.R 2c9c67  present/future cells paired by ID
#   rda_offset.R             2c9c67 + 87eac5  cells by ID, site offsets keyed on sample
#   amova.R                  52289b  sites by sample id, LD VCF's 0_ prefix handled
# geometric_offset.R needs LEA, so its case is in tests/heavy/test-heavy-offsets.R.
# The amova.R fixture also carries Phase 4a's 563755 test: amova.tsv records the SNPs its
# distance was computed on.

source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

expect_ok <- function(res) expect_identical(res$status, 0L, info = res$output)

# ---- lfmm2vcf.R (953fac) -----------------------------------------------------

test_that("lfmm2vcf.R writes the output it is given, whatever the input path holds", {
    # '_lfmm' in a folder name matched the unanchored gsub('.lfmm', '.vcf', input), which sent
    # the output to <...>/run.vcf_x/... instead of the rule's {output}.
    d <- file.path(withr::local_tempdir(), "run_lfmm_x")
    dir.create(d)
    f <- fx_imputed_triple(d)
    out <- file.path(d, "restored.vcf")
    expect_ok(run_wrapper("lfmm2vcf.R", c(f$lfmm, f$vcf, "20000", out)))
    expect_identical(readLines(out)[1:2], readLines(f$vcf)[1:2])
    body <- data.table::fread(out, skip = "#CHROM", colClasses = "character")
    expect_identical(nrow(body), length(f$ac))
    expect_identical(body[["0_S01"]], c("0/0", "0/1", "1/1")[f$g[1, ] + 1L])
    expect_setequal(list.files(dirname(d)), "run_lfmm_x")

    old <- run_wrapper("lfmm2vcf.R", c(f$lfmm, f$vcf, "20000"))
    expect_false(identical(old$status, 0L))
    expect_match(old$output, "Usage: lfmm2vcf.R")
})

# ---- gradient_forest_model.R (ec0d05) ----------------------------------------

gf_fit <- function(w, clim, out) {
    run_wrapper("gradient_forest_model.R",
                c(w$lfmm, w$snps, w$vcfsnp, w$removed, w$meta, clim, "bio_1,bio_12",
                  "100", "0.5", "without", "adaptive", out, "NULL", "NULL"))
}

test_that("gradient_forest_model.R takes each climate row by sample id, not by position", {
    d <- withr::local_tempdir()
    w <- fx_offset_world(d)
    expect_ok(gf_fit(w, w$site_pres, file.path(d, "gf_a.qs")))
    rev_clim <- fx_reorder_rows(w$site_pres, w$n:1, file.path(d, "clim_rev.tsv"))
    expect_ok(gf_fit(w, rev_clim, file.path(d, "gf_b.qs")))
    a <- qs::qread(file.path(d, "gf_a.qs")); b <- qs::qread(file.path(d, "gf_b.qs"))
    expect_s3_class(a, "gradientForest")
    expect_identical(as.data.frame(a$X), as.data.frame(b$X))
    expect_identical(a$overall.imp, b$overall.imp)

    short <- fx_reorder_rows(w$site_pres, 2:w$n, file.path(d, "clim_short.tsv"))
    res <- gf_fit(w, short, file.path(d, "gf_c.qs"))
    expect_false(identical(res$status, 0L))
    expect_match(res$output, "1 sample\\(s\\) of .* have no row in .* \\(e.g. S01\\)")
})

# ---- gradient_forest_offset.R (2c9c67) ---------------------------------------

gf_offset <- function(d, w, gf, fut_all, tag) {
    o <- file.path(d, paste0(tag, c("_offset.tif", "_map.tsv", "_site.tsv")))
    res <- run_wrapper("gradient_forest_offset.R",
                       c(gf, "bio_1,bio_12", fut_all, w$all_pres, w$raster, w$meta,
                         o[1], o[2], o[3], "s1", "TRUE", "NULL"))
    list(res = res, map = o[2], site = o[3])
}

test_that("gradient_forest_offset.R pairs present and future cells by ID", {
    d <- withr::local_tempdir()
    w <- fx_offset_world(d)
    gf <- file.path(d, "gf.qs")
    expect_ok(gf_fit(w, w$site_pres, gf))
    base <- gf_offset(d, w, gf, w$all_fut, "base")
    expect_ok(base$res)

    shuffled <- fx_reorder_rows(w$all_fut, c(51:100, 1:50), file.path(d, "fut_shuffled.tsv"))
    perm <- gf_offset(d, w, gf, shuffled, "perm")
    expect_ok(perm$res)
    expect_identical(readLines(perm$map), readLines(base$map))
    expect_identical(readLines(perm$site), readLines(base$site))

    # One cell missing from the future: that cell has no offset, every other one is unchanged.
    holed <- fx_reorder_rows(w$all_fut, -37, file.path(d, "fut_holed.tsv"))
    hole <- gf_offset(d, w, gf, holed, "hole")
    expect_ok(hole$res)
    expect_match(hole$res$output, "1 present cell\\(s\\) without a future cell")
    m_base <- as.matrix(data.table::fread(base$map, header = FALSE))
    m_hole <- as.matrix(data.table::fread(hole$map, header = FALSE))
    cell <- which(t(m_hole) %in% NA)        # values() order is row-major, as the TSV
    expect_identical(cell, 37L)
    expect_identical(t(m_hole)[-37], t(m_base)[-37])
})

# ---- rda_offset.R (2c9c67, 87eac5) -------------------------------------------

rda_offset <- function(d, w, tag, meta = w$meta, site_fut = w$site_fut, all_fut = w$all_fut) {
    o <- file.path(d, paste0(tag, c("_site.tsv", "_map.tsv", "_offset.tif", "_imp.png",
                                    "_diag.tsv")))
    plots <- file.path(d, paste0(tag, "_plots"))
    res <- run_wrapper("rda_offset.R",
                       c(w$lfmm, w$vcfsnp, w$removed, w$snps, w$site_pres, site_fut,
                         w$all_pres, all_fut, w$raster, meta, "bio_1,bio_12", "2", "0.05",
                         "99", "0", "NULL", "NULL", "NULL", "42", "1",
                         o[1], o[2], o[3], o[4], o[5], plots, "s1"))
    list(res = res, site = o[1], map = o[2])
}

test_that("rda_offset.R keys site offsets on the sample and pairs raster cells by ID", {
    d <- withr::local_tempdir()
    w <- fx_offset_world(d)
    base <- rda_offset(d, w, "base")
    expect_ok(base$res)
    site <- data.table::fread(base$site, colClasses = c(site = "character", sample = "character"))
    meta <- data.table::fread(w$meta, colClasses = "character")
    expect_identical(site$site, meta$site[match(site$sample, meta$sample)])

    # The metadata in another row order: the same offset for the same sample (87eac5).
    meta_perm <- fx_reorder_rows(w$meta, c(12:7, 1:6), file.path(d, "meta_perm.tsv"))
    m <- rda_offset(d, w, "meta", meta = meta_perm)
    expect_ok(m$res)
    expect_identical(readLines(m$site), readLines(base$site))

    # The future cells in another row order: the same map (2c9c67).
    shuffled <- fx_reorder_rows(w$all_fut, c(100:1), file.path(d, "fut_shuffled.tsv"))
    p <- rda_offset(d, w, "cells", all_fut = shuffled)
    expect_ok(p$res)
    expect_identical(readLines(p$map), readLines(base$map))
    expect_identical(readLines(p$site), readLines(base$site))

    # The future site rows in another order would pair each sample with another's future.
    fut_perm <- fx_reorder_rows(w$site_fut, c(2, 1, 3:w$n), file.path(d, "site_fut_perm.tsv"))
    s <- rda_offset(d, w, "sitefut", site_fut = fut_perm)
    expect_false(identical(s$res$status, 0L))
    expect_match(s$res$output, "not in the same sample order")
})

# ---- amova.R (52289b) ------------------------------------------------------------

# A VCF named the way the LD-pruned VCFs are, 0_<sample>, with site-structured allele
# frequencies so AMOVA has a between-site component to find.
fx_amova_vcf <- function(d, meta, n_snp = 40L) {
    p <- c(A = 0.1, B = 0.45, C = 0.75, D = 0.95)
    g <- vapply(seq_len(n_snp), function(j)
        stats::rbinom(nrow(meta), 2, pmin(0.98, p[meta$site] * stats::runif(1, 0.7, 1))),
        numeric(nrow(meta)))
    gt <- matrix(c("0/0", "0/1", "1/1")[g + 1], nrow = nrow(meta))
    vcf <- file.path(d, "ld.vcf")
    writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000000>",
                 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                 paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                         "FORMAT", paste0("0_", meta$sample)), collapse = "\t"),
                 vapply(seq_len(n_snp), function(j) paste(
                     c("1", j * 1000L, paste0("snp", j), "A", "G", ".", "PASS", ".", "GT",
                       gt[, j]), collapse = "\t"), "")), vcf)
    vcf
}

amova <- function(d, vcf, meta, tag) {
    dirs <- vapply(c("plots", "tables", "inter"), function(x) fx_inter_dir(d, paste0(tag, "_", x)), "")
    res <- run_wrapper("amova.R", c(vcf, meta, "1", dirs[["plots"]], dirs[["tables"]], dirs[["inter"]]))
    list(res = res, table = file.path(dirs[["tables"]], "amova.tsv"))
}

test_that("amova.R gives each individual its own site, by id, in any metadata order", {
    d <- withr::local_tempdir()
    set.seed(5)
    meta_dt <- data.table::data.table(site = rep(c("A", "B", "C", "D"), each = 3),
                                      sample = sprintf("S%02d", 1:12),
                                      latitude = 31, longitude = 34)
    meta <- write_tsv(meta_dt, file.path(d, "metadata.tsv"))
    vcf <- fx_amova_vcf(d, meta_dt)
    base <- amova(d, vcf, meta, "base")
    expect_ok(base$res)
    tab <- data.table::fread(base$table)
    expect_identical(tab[Parameter == "df_between", as.numeric(Value)], 3)

    # Interleaved, not reversed: a reversal of 3-sample blocks keeps every site group intact
    # under positional binding, so only a mixing order can tell the two apart.
    perm <- fx_reorder_rows(meta, c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12),
                            file.path(d, "metadata_perm.tsv"))
    p <- amova(d, vcf, perm, "perm")
    expect_ok(p$res)
    expect_identical(readLines(p$table), readLines(base$table))

    short <- fx_reorder_rows(meta, 2:12, file.path(d, "metadata_short.tsv"))
    s <- amova(d, vcf, short, "short")
    expect_false(identical(s$res$status, 0L))
    expect_match(s$res$output, "1 VCF sample\\(s\\) have no row in .* \\(e.g. 0_S01\\)")
})

AMOVA_SUBSAMPLE_ROWS <- c("n_snps_vcf", "snp_cap", "n_snps_sampled", "n_loci_used", "seed")

test_that("amova.tsv records the SNPs its distance was computed on (563755)", {
    # amova.R draws RANDOM_N = 10000 SNPs (seed 42) before bitwise.dist(), and the table
    # recorded neither the cap nor the loci used, so its percentages could not be
    # reproduced from it. A VCF above the cap tells the sampled count from the VCF's.
    d <- withr::local_tempdir()
    set.seed(6)
    meta_dt <- data.table::data.table(site = rep(c("A", "B", "C", "D"), each = 3),
                                      sample = sprintf("S%02d", 1:12),
                                      latitude = 31, longitude = 34)
    meta <- write_tsv(meta_dt, file.path(d, "metadata.tsv"))
    r <- amova(d, fx_amova_vcf(d, meta_dt, n_snp = 10050L), meta, "cap")
    expect_ok(r$res)
    tab <- data.table::fread(r$table, colClasses = "character")
    expect_identical(tail(tab$Parameter, 5), AMOVA_SUBSAMPLE_ROWS)
    num <- function(p) as.numeric(tab[Parameter == p, Value])
    expect_identical(num("n_snps_vcf"), 10050)
    expect_identical(num("snp_cap"), 10000)
    expect_identical(num("n_snps_sampled"), 10000)
    expect_identical(num("n_loci_used"), 10000)
    expect_identical(num("seed"), 42)
    # The AMOVA rows themselves are still there, in their old order.
    expect_identical(head(tab$Parameter, 8),
                     c("Between_populations_var_percent", "Within_populations_var_percent",
                       "Phi-samples-total", "p-value", "Nperm",
                       "df_between", "df_within", "df_total"))
})

test_that("a single-site amova.tsv lists the same subsample rows, NA where no VCF was read", {
    d <- withr::local_tempdir()
    set.seed(7)
    meta_dt <- data.table::data.table(site = "A", sample = sprintf("S%02d", 1:4),
                                      latitude = 31, longitude = 34)
    meta <- write_tsv(meta_dt, file.path(d, "metadata.tsv"))
    r <- amova(d, fx_amova_vcf(d, meta_dt, n_snp = 5L), meta, "one")
    expect_ok(r$res)
    tab <- data.table::fread(r$table, colClasses = "character", na.strings = c("", "NA"))
    expect_identical(tab[Parameter == "status", Value], "skipped_single_site")
    expect_identical(tail(tab$Parameter, 5), AMOVA_SUBSAMPLE_ROWS)
    expect_identical(tab[Parameter == "snp_cap", Value], "10000")
    expect_identical(tab[Parameter == "seed", Value], "42")
    expect_true(all(is.na(tab[Parameter %in% c("n_snps_vcf", "n_snps_sampled", "n_loci_used"), Value])))
})
