# Shared fixture builders and subprocess runner for the CLI wrapper smoke tests.
#
# Sourced by BOTH test roots:
#   tests/testthat/test-cli-wrappers.R   quick gate — no genomics tool needed
#   tests/heavy/test-heavy-wrappers.R    opt-in heavy gate — plink/LEA/vcftools/EMMAX
#
# It lives in tests/lib/ rather than tests/testthat/helper-*.R on purpose: the
# heavy root must be able to source it without reading a file named as the quick
# root's helper, and test_dir() must not source it twice.
#
# Side-effect free at source time apart from defining functions and the two path
# constants below.
#
# THE CWD SANDBOX, and why it is the real containment mechanism.
# run_wrapper() runs each script with its working directory set to an empty
# throwaway dir, and expect_wrapper_ok() then asserts that dir is STILL empty.
# Several wrappers resolve a path relative to the cwd — generate_simdata.R:8
# defaults OUTDIR to "data/", and base-graphics calls drop Rplots.pdf where they
# stand (geometric_offset.R:459, rda_offset.R:397). Without the sandbox those
# land in the mounted repo. A `git status` check cannot substitute: data/ and
# *_results/ are gitignored (.gitignore:4,23), so git is blind to exactly the
# writes that matter most.

REPO    <- getOption("clinego.repo_root", "/pipeline")
SCRIPTS <- file.path(REPO, "scripts")

# Wrappers that must NEVER appear in a spec table. Both take no argv at all and
# read-modify-write the tracked SIMDATA fixtures in place (add_related_samples.R:30-31,
# add_pregea_sites.R:38-39). Asserted by a test, not left to review discipline.
WRAPPER_DENYLIST <- c("add_related_samples.R", "add_pregea_sites.R")

# ------------------------------------------------------------------ fixtures

write_tsv <- function(dt, path) {
    data.table::fwrite(dt, path, sep = "\t")
    path
}

fx_metadata <- function(d, n = 12) {
    dt <- data.table::data.table(
        site      = rep(c("NEG", "TAV", "GAL"), length.out = n),
        sample    = sprintf("ID%03d", seq_len(n)),
        latitude  = 30 + seq_len(n) * 0.1,
        longitude = 34 + seq_len(n) * 0.1,
        height    = as.numeric(seq_len(n)),
        flowering_time = as.numeric(rev(seq_len(n)))
    )
    write_tsv(dt, file.path(d, "metadata.tsv"))
}

fx_climate_site <- function(d) {
    dt <- data.table::data.table(
        sample = sprintf("ID%03d", 1:8),
        bio_1  = c(1.0, 2.0, 3.0, 4.0, 2.5, 3.5, 1.5, 4.5),
        bio_2  = c(9.0, 7.0, 5.0, 3.0, 8.0, 4.0, 6.0, 2.0),
        bio_12 = rep(5.0, 8)                       # invariant on purpose
    )
    write_tsv(dt, file.path(d, "climate_present_site.tsv"))
}

fx_pvalues <- function(d, name = "pvalues.tsv") {
    n <- 20
    dt <- data.table::data.table(
        SNPID = sprintf("snp%02d", seq_len(n)),
        chr   = rep(c("1", "2"), each = n / 2),
        pos   = seq_len(n) * 10000L,
        bio_1 = c(1e-9, 1e-8, runif(n - 2, 0.05, 1)),
        bio_2 = c(runif(n - 2, 0.05, 1), 1e-9, 1e-8)
    )
    write_tsv(dt, file.path(d, name))
}

fx_sig_snps <- function(path, method, traits = c("bio_1", "bio_2")) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    dt <- data.table::data.table(
        SNPID           = c("snp01", "snp02", "snp11", "snp12"),
        chr             = c("1", "1", "2", "2"),
        pos             = c(10000L, 20000L, 110000L, 120000L),
        pvalue          = c(1e-9, 1e-8, 1e-9, 1e-8),
        pval_threshold  = 0.0025,
        method          = method,
        trait           = rep(traits, length.out = 4),
        overlap_traits  = "",
        overlap_snps    = "",
        overlap_distance = 10000L
    )
    write_tsv(dt, path)
}

fx_regions <- function(d) {
    dt <- data.table::data.table(
        region_id = c("1_5000-30000", "2_105000-130000"),
        trait     = c("bio_1", "bio_2"),
        chr       = c("1", "2"),
        start     = c(5000L, 105000L),
        end       = c(30000L, 130000L)
    )
    write_tsv(dt, file.path(d, "regions.tsv"))
}

fx_gff <- function(d) {
    p <- file.path(d, "genes.gff3")
    lines <- c(
        "##gff-version 3",
        paste("1", "src", "gene", "8000",  "15000", ".", "+", ".", "ID=g1;Name=ALPHA", sep = "\t"),
        paste("1", "src", "exon", "9000",  "10000", ".", "+", ".", "ID=e1;Parent=g1",  sep = "\t"),
        paste("2", "src", "gene", "110000", "125000", ".", "-", ".", "ID=g2;Name=BETA", sep = "\t")
    )
    writeLines(lines, p)
    p
}

fx_allsnps <- function(d) {
    p <- file.path(d, "all.vcfsnp")
    writeLines(c("1 9500 snp01 C A . . PR GT",
                 "1 12000 snp02 G T . . PR GT",
                 "2 115000 snp11 A C . . PR GT"), p)
    p
}

fx_sample_list <- function(d, name, samples) {
    p <- file.path(d, name)
    writeLines(paste(samples, samples), p)   # plink --keep: FID IID, no header
    p
}

# --- fixtures added 2026-09-12 with the second wrapper batch ---------------

# CHR/POS/MAF, the shape compute_wza.R expects (its own arg names, uppercase).
fx_maf <- function(d, name = "maf_pos.tsv") {
    n <- 20
    dt <- data.table::data.table(
        CHR = rep(c("1", "2"), each = n / 2),
        POS = seq_len(n) * 10000L,
        MAF = seq(0.05, 0.5, length.out = n)
    )
    write_tsv(dt, file.path(d, name))
}

# LFMM genotype matrix: space-separated, no header, no row names. 0/1/2 and 9 = NA.
fx_lfmm_matrix <- function(d, name = "geno.lfmm", n_ind = 6, n_snp = 8) {
    set.seed(42)
    m <- matrix(sample(0:2, n_ind * n_snp, replace = TRUE), nrow = n_ind)
    p <- file.path(d, name)
    write.table(m, p, sep = " ", row.names = FALSE, col.names = FALSE)
    p
}

# sNMF Q-matrix exactly as extract_clusters.R writes it: sample, site, C1..Ck.
# The `site` column is load-bearing — plot_structure.R:25 does select(-site) and
# errors without it.
fx_clusters <- function(d, k = 3, n = 8, name = NULL) {
    if (is.null(name)) name <- paste0("clusters_K", k, ".tsv")
    q <- matrix(runif(n * k), nrow = n)
    q <- q / rowSums(q)
    dt <- data.table::data.table(
        sample = sprintf("ID%03d", seq_len(n)),
        site   = rep(c("NEG", "TAV", "GAL"), length.out = n)
    )
    for (i in seq_len(k)) dt[[paste0("C", i)]] <- q[, i]
    write_tsv(dt, file.path(d, name))
}

# selected_snps.tsv: SNPID, chr, pos, one column PER METHOD, min_pvalue last.
# The per-method cells hold comma-separated TRAIT NAMES (or "" / literal `""`),
# not p-values — the schema misread that cost a session during Tier 4. Verified
# against SIMDATA_results/GEA/tables/selected_snps.tsv.
fx_selected_snps <- function(d, name = "selected_snps.tsv",
                             methods = c("EMMAX", "LFMM"),
                             traits = c("bio_1", "bio_2")) {
    dt <- data.table::data.table(
        SNPID = c("snp01", "snp02", "snp11", "snp12"),
        chr   = c("1", "1", "2", "2"),
        pos   = c(10000L, 20000L, 110000L, 120000L)
    )
    for (i in seq_along(methods)) {
        # Row 1 is hit by every method (a cross-method overlap), row 4 by none.
        dt[[methods[i]]] <- c(traits[1],
                              if (i == 1) paste(traits, collapse = ",") else "",
                              traits[min(i, length(traits))],
                              "")
    }
    dt$min_pvalue <- c(1e-9, 1e-8, 1e-7, 1e-6)
    write_tsv(dt, file.path(d, name))
}

# One number per line — what LEA writes and what plot_pregea_screeplot.R /
# plot_pca_structure.R read back with readLines() / fread(header = FALSE).
fx_eigenvalues <- function(d, name = "eigenvalues.txt", n = 10) {
    p <- file.path(d, name)
    writeLines(format(sort(runif(n, 1, 100), decreasing = TRUE), trim = TRUE), p)
    p
}

# Headerless space-separated PCA projections, n rows x n_pc columns.
fx_projections <- function(d, name = "projections.txt", n = 8, n_pc = 5) {
    p <- file.path(d, name)
    m <- matrix(round(rnorm(n * n_pc), 4), nrow = n)
    write.table(m, p, sep = " ", row.names = FALSE, col.names = FALSE)
    p
}

# A non-empty stand-in for a plot that a script only file.exists()-tests
# (write_summary.R:342 is the case). Deliberately NOT a real PNG: opening a
# graphics device at a small size raises "figure margins too large" in the
# container's default device, and no consumer of this fixture decodes the file.
fx_dummy_png <- function(d, name = "placeholder.png") {
    p <- file.path(d, name)
    writeLines("not a real PNG - existence-only fixture", p)
    p
}

# INTER_DIR arguments MUST carry a trailing slash. Several scripts build their
# intermediate path by string concatenation rather than file.path() — e.g.
# plot_density.R:19, plot_structure.R:18, plot_correlation_heatmap.R:37 — and
# work in production only because common.smk defines INTER with a trailing
# slash. Filed; until fixed, fixtures must honour the convention.
fx_inter_dir <- function(d, name = "inter") {
    p <- file.path(d, name)
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    paste0(p, "/")
}

# ------------------------------------------------------------------- runner

# Invoked exactly as the Snakefile does: `Rscript /pipeline/scripts/<name>.R ...`.
# The scripts are not executable (mode 644) and carry no reliable shebang, so
# executing the path directly gives "Permission denied" rather than running them.
run_wrapper <- function(script, args, wd = NULL) {
    # shQuote every argument. system2() builds a SHELL command line, so an
    # unquoted space-separated file list (combine_selected_snps.R's argv[1],
    # combine_pheno_pvalues.R's argv[1]) would be split into separate argv
    # entries and every later positional would shift by one — the script then
    # writes its output over one of its own inputs. The Snakefile quotes these
    # the same way ("{params.files_str}", gwas.smk:93, _assoc_downstream.smk:81).
    invoke <- function() {
        suppressWarnings(system2("Rscript",
                                 shQuote(c(file.path(SCRIPTS, script), args)),
                                 stdout = TRUE, stderr = TRUE))
    }
    out <- if (is.null(wd)) invoke() else withr::with_dir(wd, invoke())
    status <- attr(out, "status")
    list(status = if (is.null(status)) 0L else as.integer(status),
         output = paste(out, collapse = "\n"))
}

# A spec is a list with:
#   label    test-name suffix (two specs may share one script)
#   script   basename inside scripts/
#   build(d) writes fixtures into d, returns whatever args() needs
#   args(d, built)     character argv
#   outputs(d)         files asserted to exist and be non-empty
#   out_min(d)         OPTIONAL list(dir=, pattern=, n=) — assert at least n
#                      files match pattern in dir. For scripts that derive their
#                      own basenames (plot_manhattan.R:156-208): re-deriving them
#                      in the test would duplicate production logic.
#   tools    OPTIONAL character vector of executables/paths the heavy tier
#            requires; asserted present, never skipped (see require_tools()).
expect_wrapper_ok <- function(spec) {
    d <- withr::local_tempdir()
    sandbox <- file.path(d, "cwd")
    dir.create(sandbox, recursive = TRUE, showWarnings = FALSE)

    args <- spec$args(d, spec$build(d))
    res  <- run_wrapper(spec$script, args, wd = sandbox)

    testthat::expect_identical(
        res$status, 0L,
        info = paste0(spec$script, " exited ", res$status, "\nargs: ",
                      paste(args, collapse = " "), "\n", res$output))

    for (f in spec$outputs(d)) {
        testthat::expect_true(file.exists(f),
            info = paste0(spec$script, " declared ", f, " but did not write it\n",
                          res$output))
        testthat::expect_gt(file.size(f), 0)
    }

    if (!is.null(spec$out_min)) {
        om    <- spec$out_min(d)
        found <- list.files(om$dir, pattern = om$pattern)
        # expect_true, not expect_gte: testthat 3.2.3's expect_gte takes no
        # `info` argument, and the message is the whole point of this assertion.
        testthat::expect_true(length(found) >= om$n,
            info = paste0(spec$script, ": expected >= ", om$n, " files matching ",
                          om$pattern, " in ", om$dir, ", found ",
                          length(found), "\n", res$output))
    }

    # The cwd sandbox must be untouched: anything here is a path the script
    # resolved relative to its working directory, which in production is the
    # mounted repo.
    leaked <- list.files(sandbox, all.files = TRUE, no.. = TRUE, recursive = TRUE)
    testthat::expect_identical(
        leaked, character(0),
        info = paste0(spec$script, " wrote into its working directory: ",
                      paste(leaked, collapse = ", ")))

    invisible(d)
}

# Negative control: hand the script a first input that does not exist and require
# a non-zero exit with no output file written. This is the assertion that catches
# "silently wrote a stub and exited 0", and it is opt-in per spec because a few
# wrappers legitimately treat a declared input as optional — write_summary.R's
# read_opt() (:199,:261,:317) returns NULL for a missing file and skips the
# guarded block, which is itself filed as a defect.
expect_wrapper_fails_without_input <- function(spec) {
    d <- withr::local_tempdir()
    sandbox <- file.path(d, "cwd")
    dir.create(sandbox, recursive = TRUE, showWarnings = FALSE)

    args <- spec$args(d, spec$build(d))
    args[1] <- file.path(d, "definitely-not-here.tsv")
    outs <- spec$outputs(d)
    unlink(outs)

    res <- run_wrapper(spec$script, args, wd = sandbox)
    testthat::expect_gt(res$status, 0L)
    for (f in outs) {
        testthat::expect_false(file.exists(f),
            info = paste0(spec$script, " wrote ", f, " despite a missing input\n",
                          res$output))
    }
}

# Heavy tier only. An absent tool inside cline-go:latest means the IMAGE is
# wrong, which must be red — never a skip. A skip here would be the exact
# failure mode the quick/heavy split exists to avoid.
require_tools <- function(tools) {
    for (t in tools) {
        ok <- if (grepl("/", t)) file.exists(t) else nzchar(Sys.which(t))
        testthat::expect_true(ok, info = paste0(
            "required tool missing from the image: ", t,
            " — heavy tests fail rather than skip on this, see tests/heavy/"))
    }
}
