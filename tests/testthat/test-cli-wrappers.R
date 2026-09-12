# Smoke tests for the thin R CLI wrappers in scripts/.
#
# 40 of the 78 files in scripts/*.R define no functions at all: they read
# commandArgs(trailingOnly=TRUE) positionally, call a library in scripts/R/lib,
# and write a file. Tier 1 tested the libraries BEHIND them. Nothing tested the
# argument order IN FRONT of them, which is their entire failure mode — swap two
# args and the script still runs, just on the wrong file.
#
# Each case builds a fixture in a temp dir, invokes the script exactly as the
# Snakefile would, and asserts exit 0 plus every declared output present and
# non-empty. That is deliberately shallow: this file checks PLUMBING. What the
# numbers should be is asserted in the Tier 1 lib tests.
#
# Two constraints worth knowing before extending the table:
#   * Every wrapper source()s /pipeline/scripts/R/... with an ABSOLUTE path, so
#     these run only inside the container with the repo mounted at /pipeline.
#   * A wrapper's outputs are not always in argv. vcf2lfmm.R derives
#     <base>.lfmm + <base>.lfmm_nmissing from its INPUT path, and
#     plot_cross_entropy.R derives the SVG and QS siblings of the PNG it is
#     given (:19-20). `outputs` must list what is written, not what is passed.
#
# Out of scope on purpose: every wrapper needing a VCF, plink, LEA/sNMF, EMMAX,
# GAPIT or a WorldClim raster. Those are integration runs, not smoke tests.

REPO  <- getOption("clinego.repo_root", "/pipeline")
SCRIPTS <- file.path(REPO, "scripts")

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

# ------------------------------------------------------------------- runner

# Invoked exactly as the Snakefile does: `Rscript /pipeline/scripts/<name>.R ...`.
# The scripts are not executable (mode 644) and carry no reliable shebang, so
# executing the path directly gives "Permission denied" rather than running them.
run_wrapper <- function(script, args) {
    # shQuote every argument. system2() builds a SHELL command line, so an
    # unquoted space-separated file list (combine_selected_snps.R's argv[1],
    # combine_pheno_pvalues.R's argv[1]) would be split into separate argv
    # entries and every later positional would shift by one — the script then
    # writes its output over one of its own inputs. The Snakefile quotes these
    # the same way ("{params.files_str}", gwas.smk:93, _assoc_downstream.smk:81).
    out <- suppressWarnings(system2("Rscript",
                                    shQuote(c(file.path(SCRIPTS, script), args)),
                                    stdout = TRUE, stderr = TRUE))
    status <- attr(out, "status")
    list(status = if (is.null(status)) 0L else as.integer(status),
         output = paste(out, collapse = "\n"))
}

expect_wrapper_ok <- function(spec) {
    d <- withr::local_tempdir()
    args <- spec$args(d, spec$build(d))
    res  <- run_wrapper(spec$script, args)

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
    invisible(d)
}

# -------------------------------------------------------------------- table

WRAPPERS <- list(
    list(
        label   = "trait_summary.R",
        script  = "trait_summary.R",
        build   = function(d) fx_metadata(d),
        args    = function(d, meta) c(meta, file.path(d, "trait_summary.tsv")),
        outputs = function(d) file.path(d, "trait_summary.tsv")
    ),
    list(
        label   = "check_climate_variance.R (bio)",
        script  = "check_climate_variance.R",
        build   = function(d) fx_climate_site(d),
        args    = function(d, site) c(site, file.path(d, "invariant.tsv"), "bio"),
        outputs = function(d) file.path(d, "invariant.tsv")
    ),
    list(
        label   = "check_climate_variance.R (traits)",
        script  = "check_climate_variance.R",
        build   = function(d) fx_metadata(d),
        args    = function(d, meta) c(meta, file.path(d, "invariant_traits.tsv"), "traits"),
        outputs = function(d) file.path(d, "invariant_traits.tsv")
    ),
    list(
        # The only variadic wrapper: args[seq_len(n-1)] are inputs, args[n] is the
        # output. Swap the order and it writes over an input — the best
        # argument-order canary in the set.
        label   = "assemble_pvalues.R (variadic)",
        script  = "assemble_pvalues.R",
        build   = function(d) {
            a <- data.table::data.table(SNPID = c("s1", "s2"), chr = c("1", "1"),
                                        pos = c(100L, 200L), height = c(0.01, 0.2))
            b <- data.table::data.table(SNPID = c("s1", "s2"), chr = c("1", "1"),
                                        pos = c(100L, 200L), flowering_time = c(0.3, 0.04))
            c(write_tsv(a, file.path(d, "t_height.tsv")),
              write_tsv(b, file.path(d, "t_ft.tsv")))
        },
        args    = function(d, files) c(files, file.path(d, "wide.tsv")),
        outputs = function(d) file.path(d, "wide.tsv")
    ),
    list(
        label   = "filter_arrange_metadata.R",
        script  = "filter_arrange_metadata.R",
        build   = function(d) c(fx_metadata(d),
                                fx_sample_list(d, "vcf_samples.list",
                                               sprintf("ID%03d", c(3, 1, 2)))),
        args    = function(d, f) c(f[1], f[2], file.path(d, "metadata_ordered.tsv")),
        outputs = function(d) file.path(d, "metadata_ordered.tsv")
    ),
    list(
        label   = "filter_coord_samples.R",
        script  = "filter_coord_samples.R",
        build   = function(d) fx_metadata(d),
        args    = function(d, meta) c(meta,
                                      file.path(d, "coord_valid.list"),
                                      file.path(d, "metadata_climate.tsv"),
                                      file.path(d, "coord_missing_summary.tsv")),
        outputs = function(d) file.path(d, c("coord_valid.list", "metadata_climate.tsv",
                                             "coord_missing_summary.tsv"))
    ),
    list(
        label   = "filter_climate_valid_samples.R",
        script  = "filter_climate_valid_samples.R",
        build   = function(d) {
            meta <- fx_metadata(d)
            lst  <- fx_sample_list(d, "coord_valid.list", sprintf("ID%03d", 1:12))
            excl <- write_tsv(data.table::data.table(
                        sample = character(), site = character(),
                        latitude = numeric(), longitude = numeric(),
                        reason = character(), distance_km = numeric()),
                    file.path(d, "climate_na_excluded.tsv"))
            c(lst, meta, excl)
        },
        args    = function(d, f) c(f[1], f[2], f[3],
                                   file.path(d, "climate_valid.list"),
                                   file.path(d, "metadata_climate_out.tsv")),
        outputs = function(d) file.path(d, c("climate_valid.list",
                                             "metadata_climate_out.tsv"))
    ),
    list(
        label   = "find_sig_snps.R",
        script  = "find_sig_snps.R",
        build   = function(d) fx_pvalues(d),
        # cpu = 1 deliberately: sig_snps.R:44 superassigns inside mclapply, so the
        # diagnostics slot empties at cpu >= 2 (filed). Nothing here reads it, but
        # a smoke test should not depend on the core count either way.
        args    = function(d, pv) c(pv, "bonf_0.05", "10000", "EMMAX", "1",
                                    file.path(d, "sig_snps.tsv")),
        outputs = function(d) file.path(d, "sig_snps.tsv")
    ),
    list(
        label   = "create_regions.R",
        script  = "create_regions.R",
        build   = function(d) {
            fx_sig_snps(file.path(d, "EMMAX", "sig.tsv"), "EMMAX")
            file.path(d, "EMMAX", "sig.tsv")
        },
        args    = function(d, sig) c(sig, "10000",
                                     file.path(d, "regions_per_trait.tsv"),
                                     file.path(d, "regions_combined.tsv")),
        outputs = function(d) file.path(d, c("regions_per_trait.tsv",
                                             "regions_combined.tsv"))
    ),
    list(
        # combine_selected_snps.R:31 derives the method name from each file's
        # PARENT DIRECTORY, so the fixture must live at <d>/<METHOD>/<file>.
        label   = "combine_selected_snps.R",
        script  = "combine_selected_snps.R",
        build   = function(d) {
            a <- file.path(d, "EMMAX", "sig.tsv")
            b <- file.path(d, "LFMM",  "sig.tsv")
            fx_sig_snps(a, "EMMAX")
            fx_sig_snps(b, "LFMM")
            paste(a, b)
        },
        args    = function(d, files) c(files, "Union", "10000", "bio_1,bio_2",
                                       file.path(d, "selected_snps.tsv")),
        outputs = function(d) file.path(d, "selected_snps.tsv")
    ),
    list(
        label   = "combine_pheno_pvalues.R",
        script  = "combine_pheno_pvalues.R",
        build   = function(d) {
            n <- 30
            mk <- function(trait) {
                dt <- data.table::data.table(SNPID = sprintf("s%02d", seq_len(n)),
                                             chr = "1", pos = seq_len(n) * 100L)
                dt[[trait]] <- runif(n)
                write_tsv(dt, file.path(d, paste0("p_", trait, ".tsv")))
            }
            paste(mk("height"), mk("flowering_time"))
        },
        args    = function(d, files) c(files, file.path(d, "pvalues.tsv"),
                                       file.path(d, "qvalues.tsv")),
        outputs = function(d) file.path(d, c("pvalues.tsv", "qvalues.tsv"))
    ),
    list(
        label   = "find_genes_around_regions.R",
        script  = "find_genes_around_regions.R",
        build   = function(d) c(fx_gff(d), fx_regions(d), fx_allsnps(d)),
        args    = function(d, f) c(f[1], f[2], "gene", "1000", f[3], "1",
                                   file.path(d, "genes_per_region.tsv"),
                                   file.path(d, "genes_collapsed.tsv")),
        outputs = function(d) file.path(d, c("genes_per_region.tsv",
                                             "genes_collapsed.tsv"))
    )
)

# --------------------------------------------------------------------- tests

for (spec in WRAPPERS) {
    local({
        s <- spec
        test_that(paste0("wrapper runs and writes its outputs: ", s$label), {
            skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
            skip_if_not(nzchar(Sys.which("Rscript")), "Rscript not on PATH")
            expect_wrapper_ok(s)
        })
    })
}

test_that("find_genes_around_regions.R actually annotates the overlapping genes", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    gff <- fx_gff(d); reg <- fx_regions(d); snps <- fx_allsnps(d)
    out <- file.path(d, "genes_per_region.tsv")
    res <- run_wrapper("find_genes_around_regions.R",
                       c(gff, reg, "gene", "1000", snps, "1",
                         out, file.path(d, "genes_collapsed.tsv")))
    expect_identical(res$status, 0L, info = res$output)

    genes <- data.table::fread(out, colClasses = c(chr = "character"))
    # Both fixture genes sit inside their region, so a plumbing error that fed
    # the GFF and the regions table to each other's argument would show up here
    # as an empty table rather than as a non-zero exit.
    expect_gt(nrow(genes), 0)
    expect_true(all(c("region_id", "gene_id", "chr") %in% names(genes)))
})

test_that("a wrapper given a missing input fails loudly instead of writing a stub", {
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    out <- file.path(d, "trait_summary.tsv")
    res <- run_wrapper("trait_summary.R", c(file.path(d, "nope.tsv"), out))
    expect_gt(res$status, 0L)
    expect_false(file.exists(out))
})
