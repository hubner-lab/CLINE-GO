# scripts/R/utils/emmax_core.R
#
# WHY THIS FILE SOURCES ITS OWN LIB INSTEAD OF USING helper-libs.R.
# run_emmax() dispatches through two rebindable string globals, EMMAX_BIN and
# EMMAX_RUN (:5-6), and the only hermetic way to exercise its stop-on-nonzero
# contract is to point them at /bin/false. test_dir() sources every test file
# into ONE shared parent environment, so rebinding them there would leak into
# every alphabetically-later file. sys.source() into a private env keeps it local.
#
# The headline target is load_pca_covariates(). Its k <= 0 guard exists to stop a
# silent-wrong-output bug documented in the code itself (:52-57): the LEA branch
# `1:k` -> `1:0` errors on index 0, but the PLINK branch `3:(2+k)` -> `3:2`
# does NOT — it quietly selects two wrong columns in reversed order. A guard that
# only prevents a crash on one of two branches is exactly the thing to pin.

emmax_env <- local({
    e <- new.env(parent = globalenv())
    sys.source(file.path(getOption("clinego.repo_root", "/pipeline"),
                         "scripts", "R", "utils", "emmax_core.R"),
               envir = e)
    e
})

# LEA projections: plain space-separated numeric matrix, no FID/IID prefix.
write_lea_pca <- function(d, n = 6, n_pc = 5, name = "projections.txt") {
    p <- file.path(d, name)
    m <- matrix(round(seq_len(n * n_pc) / 10, 3), nrow = n)
    write.table(m, p, sep = " ", row.names = FALSE, col.names = FALSE)
    p
}

# PLINK eigenvec: a constant FID column, a unique IID column, then the PCs.
write_plink_pca <- function(d, n = 6, n_pc = 5, name = "plink.eigenvec") {
    p <- file.path(d, name)
    m <- data.frame(FID = rep("FAM", n), IID = sprintf("ID%03d", seq_len(n)))
    for (i in seq_len(n_pc)) m[[paste0("V", i)]] <- round(seq_len(n) + i / 10, 3)
    write.table(m, p, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    p
}

test_that("the globals are plain strings, so nothing executes at source time", {
    expect_type(emmax_env$EMMAX_BIN, "character")
    expect_type(emmax_env$EMMAX_RUN, "character")
    # Sourcing this file must not require the binaries to exist or be runnable —
    # which is what made the "unsourceable" claim in helper-libs.R wrong.
    expect_length(emmax_env$EMMAX_BIN, 1L)
})

test_that("run_emmax stops when the wrapper exits non-zero", {
    e <- new.env(parent = emmax_env)
    e$EMMAX_RUN <- "/bin/false"
    e$EMMAX_BIN <- "ignored"
    f <- emmax_env$run_emmax
    environment(f) <- e
    expect_error(quiet(f("-v -d 10")), "EMMAX exited with status")
})

test_that("run_emmax returns quietly when the wrapper exits zero", {
    e <- new.env(parent = emmax_env)
    e$EMMAX_RUN <- "/bin/true"
    e$EMMAX_BIN <- "ignored"
    f <- emmax_env$run_emmax
    environment(f) <- e
    expect_identical(quiet(f("-v")), 0L)
})

test_that("load_pca_covariates returns NULL for k = 0 without reading the file", {
    # The kinship-only model in preGEA's #PC ladder. The path is deliberately
    # nonexistent: the guard must fire before any fread().
    got <- quiet(emmax_env$load_pca_covariates("/nonexistent/pca.txt", 0))
    expect_null(got)
})

test_that("load_pca_covariates returns NULL for a negative k too", {
    got <- quiet(emmax_env$load_pca_covariates("/nonexistent/pca.txt", -1))
    expect_null(got)
})

test_that("load_pca_covariates reads the LEA format from column 1", {
    d <- withr::local_tempdir()
    p <- write_lea_pca(d, n = 6, n_pc = 5)
    got <- quiet(emmax_env$load_pca_covariates(p, 3))

    expect_identical(names(got), c("PC1", "PC2", "PC3"))
    expect_identical(nrow(got), 6L)
    raw <- data.table::fread(p, sep = " ", header = FALSE)
    expect_equal(got$PC1, raw[[1]])
    expect_equal(got$PC3, raw[[3]])
})

test_that("load_pca_covariates autodetects PLINK eigenvec and skips FID/IID", {
    d <- withr::local_tempdir()
    p <- write_plink_pca(d, n = 6, n_pc = 5)
    got <- quiet(emmax_env$load_pca_covariates(p, 2))

    expect_identical(names(got), c("PC1", "PC2"))
    raw <- data.table::fread(p, sep = " ", header = FALSE)
    # Columns 3 and 4, in order — NOT columns 1-2 and NOT reversed.
    expect_equal(got$PC1, raw[[3]])
    expect_equal(got$PC2, raw[[4]])
})

test_that("load_pca_covariates subsets rows in DROP mode", {
    d <- withr::local_tempdir()
    p <- write_lea_pca(d, n = 6, n_pc = 5)
    all_samples <- sprintf("ID%03d", 1:6)
    order_file <- file.path(d, "samples_order.list")
    writeLines(all_samples, order_file)
    vcf_samples <- all_samples[c(2, 4, 5)]

    got <- quiet(emmax_env$load_pca_covariates(p, 2,
                                               samples_order_path = order_file,
                                               vcf_samples = vcf_samples))
    expect_identical(nrow(got), 3L)
    raw <- data.table::fread(p, sep = " ", header = FALSE)
    # Positional binding is the whole point: rows must follow the ORDER of the
    # samples_order file, not the order of vcf_samples.
    expect_equal(got$PC1, raw[[1]][c(2, 4, 5)])
})

test_that("load_pca_covariates errors when rows and VCF samples disagree", {
    d <- withr::local_tempdir()
    p <- write_lea_pca(d, n = 6, n_pc = 5)
    # 4 VCF samples, 6 PCA rows, and no samples_order_path to subset with.
    expect_error(
        quiet(emmax_env$load_pca_covariates(p, 2,
                                            vcf_samples = sprintf("ID%03d", 1:4))),
        "does not match VCF sample count")
})

test_that("load_pca_covariates accepts a matching row count with no subsetting", {
    d <- withr::local_tempdir()
    p <- write_lea_pca(d, n = 4, n_pc = 5)
    got <- quiet(emmax_env$load_pca_covariates(p, 2,
                                               vcf_samples = sprintf("ID%03d", 1:4)))
    expect_identical(nrow(got), 4L)
})

# --- read_vcf_samples (finding ab75b3) ---------------------------------------

write_vcf_header_only <- function(d, samples, n_meta) {
    p <- file.path(d, "long_header.vcf")
    writeLines(c("##fileformat=VCFv4.2",
                 sprintf("##contig=<ID=scaffold_%d,length=1000>", seq_len(n_meta)),
                 paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                         "FORMAT", samples), collapse = "\t"),
                 paste(c("scaffold_1", "100", ".", "A", "G", ".", "PASS", ".", "GT",
                         rep("0/1", length(samples))), collapse = "\t")), p)
    p
}

test_that("read_vcf_samples finds #CHROM past header line 1000", {
    # head -1000 never reached a #CHROM pushed down by one ##contig line per scaffold,
    # so the function returned no samples and emmax.R bound phenotypes against nothing.
    d <- withr::local_tempdir()
    p <- write_vcf_header_only(d, c("S1", "S2", "S3"), n_meta = 1500)
    expect_identical(emmax_env$read_vcf_samples(p), c("S1", "S2", "S3"))
})

test_that("read_vcf_samples keeps numeric-looking ids as written", {
    d <- withr::local_tempdir()
    p <- write_vcf_header_only(d, c("001", "002", "010"), n_meta = 2)
    expect_identical(emmax_env$read_vcf_samples(p), c("001", "002", "010"))
})

test_that("read_vcf_samples ignores a header line that merely contains CHROM", {
    d <- withr::local_tempdir()
    p <- file.path(d, "chromosome_contig.vcf")
    writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=CHROMOSOME_I,length=1000>",
                 "##source=#CHROM-like text",
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB"), p)
    expect_identical(emmax_env$read_vcf_samples(p), c("A", "B"))
})

# --- tfam_phenotype (finding df2fe3) -----------------------------------------
# emmax.R used to cbind() the trait column onto the TFAM ids, so a trait table in another
# sample order gave EMMAX each sample's id with another sample's value, and no error.

tfam_fixture <- function() {
    data.table::data.table(FAMID = c("S1", "S2", "S3", "S4"),
                           INDID = c("S1", "S2", "S3", "S4"))
}

test_that("tfam_phenotype looks each trait value up by IID, whatever the trait table's order", {
    trait <- data.table::data.table(sample = c("S3", "S1", "S4", "S2"),
                                    bio_1 = c(30, 10, 40, 20), bio_2 = c(3, 1, 4, 2))
    got <- emmax_env$tfam_phenotype(tfam_fixture(), trait, "bio_1")
    expect_identical(names(got), c("FAMID", "INDID", "bio_1"))
    expect_identical(got$INDID, c("S1", "S2", "S3", "S4"))
    expect_identical(got$bio_1, c(10, 20, 30, 40))
})

test_that("tfam_phenotype equals the old positional cbind when the orders agree", {
    trait <- data.table::data.table(sample = c("S1", "S2", "S3", "S4"),
                                    bio_1 = c(0.25, NA, -1.5, 2))
    got <- emmax_env$tfam_phenotype(tfam_fixture(), trait, "bio_1")
    old <- cbind(tfam_fixture(), trait[, "bio_1"])
    expect_identical(capture.output(write.table(got, sep = "\t", col.names = FALSE,
                                                row.names = FALSE, quote = FALSE)),
                     capture.output(write.table(old, sep = "\t", col.names = FALSE,
                                                row.names = FALSE, quote = FALSE)))
})

test_that("tfam_phenotype stops on a TFAM individual the trait table does not list", {
    trait <- data.table::data.table(sample = c("S1", "S2", "S4", "S5"), bio_1 = 1:4)
    expect_error(emmax_env$tfam_phenotype(tfam_fixture(), trait, "bio_1"),
                 "1 TFAM individual\\(s\\) have no row in the trait table \\(e.g. S3\\)")
})

test_that("tfam_phenotype stops on a duplicated sample id or a missing key column", {
    dup <- data.table::data.table(sample = c("S1", "S2", "S3", "S4", "S2"), bio_1 = 1:5)
    expect_error(emmax_env$tfam_phenotype(tfam_fixture(), dup, "bio_1"),
                 "more than once \\(e.g. S2\\)")
    nokey <- data.table::data.table(id = c("S1", "S2", "S3", "S4"), bio_1 = 1:4)
    expect_error(emmax_env$tfam_phenotype(tfam_fixture(), nokey, "bio_1"),
                 "no `sample` column")
    expect_error(emmax_env$tfam_phenotype(tfam_fixture(), dup[1:4], "bio_9"),
                 "trait 'bio_9' is not a column")
})
