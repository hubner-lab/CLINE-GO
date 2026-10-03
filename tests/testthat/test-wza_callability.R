# WZA floor reporting — CONTENT assertions, quick tier.
#
# The WZA rank transform maps the best per-SNP p to 1/(m+1), so the per-SNP z is
# capped at qnorm(1 - 1/(m+1)). A window holding ONE SNP has Z_W = z, so its WZA p
# cannot fall below that cap's p — and with most windows single-SNP, Bonferroni is
# unreachable however strong the association. A header-only sig-windows file then
# looks exactly like a genuine absence of signal (finding 6ab65e).
#
# compute_wza.R's existing wrapper row asserts exit 0 and a non-empty table, which
# the defect passes. These assert the LOG, which is where the finding's own evidence
# came from (wza_LFMM.log:28) and where a user looks at an empty result.

source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

# n SNPs spread over `spacing` bp so a window of `window` bp holds a predictable
# number of them. bio_1 carries a strong signal at the first SNP.
fx_wza_inputs <- function(d, n = 200L, spacing = 10000L) {
    set.seed(42)
    pos <- seq_len(n) * spacing
    pv <- write_tsv(data.table::data.table(
        SNPID = sprintf("snp%03d", seq_len(n)),
        chr   = "1",
        pos   = as.integer(pos),
        bio_1 = c(1e-30, 1e-25, 1e-20, runif(n - 3, 0.01, 1))),
        file.path(d, "pvalues.tsv"))
    maf <- write_tsv(data.table::data.table(
        CHR = "1", POS = as.integer(pos),
        MAF = seq(0.1, 0.5, length.out = n)),
        file.path(d, "maf_pos.tsv"))
    list(pv = pv, maf = maf)
}

run_wza <- function(d, window_bp) {
    f <- fx_wza_inputs(d)
    out <- file.path(d, paste0("wza_", window_bp, ".tsv"))
    res <- run_wrapper("compute_wza.R",
                       c(f$pv, f$maf, "NULL", as.character(window_bp),
                         as.character(window_bp), "All", out))
    list(res = res, out = out)
}

test_that("a single-SNP-dominated WZA reports its floor and names it uncallable", {
    d <- withr::local_tempdir()
    # 10 kb windows over 10 kb SNP spacing: every window holds exactly one SNP.
    r <- run_wza(d, 10000L)

    expect_identical(r$res$status, 0L, info = r$res$output)
    # The structural facts, all of them, in one line per trait.
    expect_match(r$res$output, "n_snps_ranked=200")
    expect_match(r$res$output, "z_cap=2\\.5")          # qnorm(1 - 1/201) = 2.576
    expect_match(r$res$output, "single_snp_windows=")
    expect_match(r$res$output, "min_achievable_p_single_snp=")
    expect_match(r$res$output, "UNCALLABLE BY CONSTRUCTION")

    # The floor is BINDING, which is the whole claim: the best observed window p
    # equals the best ACHIEVABLE p for a single-SNP window, so the 1e-30 p-value in
    # the input bought nothing at all.
    grab <- function(key) {
        m <- regmatches(r$res$output,
                        regexpr(paste0(key, "=[0-9.e+-]+"), r$res$output))
        expect_length(m, 1L)
        as.numeric(sub(paste0(key, "="), "", m, fixed = TRUE))
    }
    expect_equal(grab("observed_min_p"), grab("min_achievable_p_single_snp"),
                 tolerance = 1e-6)
    # ... and that floor is well above the Bonferroni threshold it is compared with.
    expect_gt(grab("observed_min_p"), grab("bonf_0.05"))
})

test_that("wide windows lift the floor and the uncallable warning goes away", {
    d <- withr::local_tempdir()
    # 1 Mb windows over 10 kb spacing: ~100 SNPs per window, so Z_W is a real sum
    # and is no longer bounded by one capped z.
    r <- run_wza(d, 1000000L)

    expect_identical(r$res$status, 0L, info = r$res$output)
    expect_false(grepl("UNCALLABLE BY CONSTRUCTION", r$res$output, fixed = TRUE))
})

test_that("the floor report does not alter the WZA p-values", {
    d <- withr::local_tempdir()
    r <- run_wza(d, 10000L)
    tbl <- data.table::fread(r$out)
    # Reporting is additive: the output schema is untouched and no new file appears.
    expect_identical(names(tbl), c("SNPID", "chr", "pos", "n_snps", "mean_maf", "bio_1"))
    expect_true(all(tbl$bio_1 >= 0 & tbl$bio_1 <= 1, na.rm = TRUE))
    expect_identical(list.files(d, pattern = "diagnostic"), character(0))
})

test_that("find_sig_snps --wza says 0 called and names the margin", {
    d <- withr::local_tempdir()
    r <- run_wza(d, 10000L)
    sig <- file.path(d, "sig.tsv")
    res <- run_wrapper("find_sig_snps.R",
                       c(r$out, "bonf_0.05", "5000", "LFMM", "1", sig, "--wza"))

    expect_identical(res$status, 0L, info = res$output)
    expect_match(res$output, "0 of [0-9]+ windows called")
    expect_match(res$output, "Best window p=")
    expect_match(res$output, "before reading this as absence of signal")
})

test_that("the per-SNP path is untouched by the WZA callability report", {
    d <- withr::local_tempdir()
    pv <- fx_pvalues(d)
    sig <- file.path(d, "sig.tsv")
    # No --wza: the same "nothing passed" situation must NOT get the WZA wording.
    res <- run_wrapper("find_sig_snps.R",
                       c(pv, "bonf_0.05", "5000", "LFMM", "1", sig))
    expect_identical(res$status, 0L, info = res$output)
    expect_false(grepl("windows called", res$output, fixed = TRUE))
    expect_false(grepl("absence of signal", res$output, fixed = TRUE))
})
