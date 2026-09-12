# Smoke tests for the thin R CLI wrappers in scripts/ — QUICK tier.
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
# The fixture builders and the runner live in tests/lib/wrapper_harness.R,
# shared with the heavy tier (tests/heavy/test-heavy-wrappers.R). Read that file
# before adding a row — in particular the cwd-sandbox rationale and the
# trailing-slash requirement on every INTER_DIR argument.
#
# Two constraints worth knowing before extending the table:
#   * Every wrapper source()s /pipeline/scripts/R/... with an ABSOLUTE path, so
#     these run only inside the container with the repo mounted at /pipeline.
#   * A wrapper's outputs are not always in argv. plot_density.R derives its SVG
#     and QS siblings from the PNG it is given (:18-19); plot_manhattan.R derives
#     every basename from PLOT_DIR (:156-208) — use `out_min` for those rather
#     than re-deriving production names here.
#
# THIS TIER TAKES NO GENOMICS TOOL. Anything needing a VCF, plink, LEA, EMMAX,
# GAPIT, vcftools, a raster or the network belongs in tests/heavy/ — which is
# expected GREEN, merely slow, and is run with `tests/run_all.sh --heavy`.

source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

# -------------------------------------------------------------------- table

WRAPPERS <- list(
    list(
        label   = "trait_summary.R",
        script  = "trait_summary.R",
        build   = function(d) fx_metadata(d),
        args    = function(d, meta) c(meta, file.path(d, "trait_summary.tsv")),
        outputs = function(d) file.path(d, "trait_summary.tsv"),
        fails_without_input = TRUE
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
        outputs = function(d) file.path(d, "metadata_ordered.tsv"),
        fails_without_input = TRUE
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
                                             "coord_missing_summary.tsv")),
        fails_without_input = TRUE
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
        outputs = function(d) file.path(d, "sig_snps.tsv"),
        fails_without_input = TRUE
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
                                             "regions_combined.tsv")),
        fails_without_input = TRUE
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
    ),

    # ===================================================================
    # Added 2026-09-12. All of these were previously listed as "out of
    # scope: needs a VCF / plink / LEA / raster" — measured false. None
    # of them touches a genomics tool: plot_pca_structure.R and
    # plot_pregea_screeplot.R read LEA *output text* with fread(header =
    # FALSE) / readLines() and never call library(LEA).
    # ===================================================================

    list(
        label   = "subset_lfmm_matrix.R",
        script  = "subset_lfmm_matrix.R",
        build   = function(d) {
            samples <- sprintf("ID%03d", 1:6)
            m     <- fx_lfmm_matrix(d, n_ind = 6, n_snp = 8)
            order <- file.path(d, "samples_order.list")
            writeLines(samples, order)               # one per line, NOT FID IID
            keep  <- fx_sample_list(d, "coord_valid.list", samples[1:4])
            c(m, order, keep)
        },
        args    = function(d, f) c(f[1], f[2], f[3], file.path(d, "geno_sub.lfmm")),
        outputs = function(d) file.path(d, "geno_sub.lfmm"),
        fails_without_input = TRUE
    ),
    list(
        label   = "prepare_phenotypes.R (MEAN)",
        script  = "prepare_phenotypes.R",
        build   = function(d) fx_metadata(d),
        args    = function(d, meta) c(meta, "MEAN", file.path(d, "pheno_mean"),
                                      file.path(d, "missing_mean.tsv")),
        outputs = function(d) c(file.path(d, "missing_mean.tsv"),
                                file.path(d, "pheno_mean", "all_phenotypes.tsv")),
        fails_without_input = TRUE
    ),
    list(
        label   = "prepare_phenotypes.R (MEDIAN)",
        script  = "prepare_phenotypes.R",
        build   = function(d) fx_metadata(d),
        args    = function(d, meta) c(meta, "MEDIAN", file.path(d, "pheno_median"),
                                      file.path(d, "missing_median.tsv")),
        outputs = function(d) c(file.path(d, "missing_median.tsv"),
                                file.path(d, "pheno_median", "all_phenotypes.tsv"))
    ),
    list(
        # DROP takes the per-trait branch (:110-131) instead of the single
        # all_phenotypes.tsv one, so it writes a different set of files. That
        # divergence is the reason all three strategies get a row.
        label   = "prepare_phenotypes.R (DROP)",
        script  = "prepare_phenotypes.R",
        build   = function(d) fx_metadata(d),
        args    = function(d, meta) c(meta, "DROP", file.path(d, "pheno_drop"),
                                      file.path(d, "missing_drop.tsv")),
        outputs = function(d) c(file.path(d, "missing_drop.tsv"),
                                file.path(d, "pheno_drop", "height_phenotype.tsv"),
                                file.path(d, "pheno_drop", "height_samples.list")),
        out_min = function(d) list(dir = file.path(d, "pheno_drop"),
                                   pattern = "_phenotype\\.tsv$", n = 2)
    ),
    list(
        label   = "promote_snp_set.R",
        script  = "promote_snp_set.R",
        build   = function(d) fx_selected_snps(d),
        args    = function(d, sel) c(sel, "testset", file.path(d, "snp_sets"), "GEA"),
        outputs = function(d) c(file.path(d, "snp_sets", "testset", "selected_snps.tsv"),
                                file.path(d, "snp_sets", "manifest.json")),
        fails_without_input = TRUE
    ),
    list(
        # write_summary.R is the one row that must NOT get the negative control:
        # read_opt() (:317) returns NULL for a missing declared input and the
        # guarded block is simply skipped, so it exits 0 with a thinner table.
        # That is filed as a defect, and quarantined in test-known-bugs.R.
        label   = "write_summary.R (traits)",
        script  = "write_summary.R",
        build   = function(d) {
            ts <- write_tsv(data.table::data.table(
                      trait = c("height", "flowering_time"), n = c(47L, 44L),
                      n_missing = c(0L, 3L), pct_missing = c(0, 6.38),
                      mean = c(36.1, 18.9), sd = c(5.9, 4.8),
                      min = c(26.3, 11.3), median = c(36.2, 18.0), max = c(47, 26)),
                  file.path(d, "trait_summary.tsv"))
            inv <- write_tsv(data.table::data.table(predictor = character(),
                                                    reason = character()),
                             file.path(d, "trait_invariant.tsv"))
            c(ts, inv, fx_dummy_png(d, "trait_pairs.png"))
        },
        args    = function(d, f) c("traits", file.path(d, "pipeline_summary.tsv"),
                                   f[1], f[2], f[3]),
        outputs = function(d) file.path(d, "pipeline_summary.tsv"),
        fails_without_input = FALSE
    ),
    # write_summary.R's OTHER two modes. Same OUTPUT contract (argv[2] is the
    # only file it ever writes, read-modify-write via update_summary at :31-37),
    # entirely different argv beyond it.
    list(
        # The gwas_only shape. summary.smk:152-158 always passes all NINE
        # positionals and uses the literal string "NULL" for the climate ones
        # when Climate.enabled is false — it never sends a short argv.
        #
        # That distinction matters and cost a red run to learn: only args 6-9 are
        # `length(args) >= n` guarded (:352-355). CLIMATE_SITE = args[4] and
        # PREDICTORS = args[5] are NOT, so a genuinely short argv makes args[4]
        # NA and `if (NA != "NULL")` dies with "missing value where TRUE/FALSE
        # needed". Filed; this row deliberately passes the shape the rule
        # actually emits rather than pinning the crash.
        label   = "write_summary.R (structure, gwas_only shape)",
        script  = "write_summary.R",
        build   = function(d) fx_ld_decay_table(d),
        args    = function(d, f) c("structure", file.path(d, "pipeline_summary.tsv"),
                                   "3", "NULL", "NULL", f, "site", "both", "NULL"),
        outputs = function(d) file.path(d, "pipeline_summary.tsv"),
        fails_without_input = FALSE
    ),
    list(
        label   = "write_summary.R (structure, all optional inputs)",
        script  = "write_summary.R",
        build   = function(d) {
            # NOT fx_climate_site: that fixture is 4 columns, and :362 computes
            # n_climate_vars as ncol - 4 assuming site/sample/lat/lon, which
            # would bake the known-wrong 0 into this row. The real shape is
            # those four identity columns PLUS one column per predictor.
            cs <- write_tsv(data.table::data.table(
                      site = c("NEG", "TAV"), sample = c("ID001", "ID002"),
                      latitude = c(30.85, 32.08), longitude = c(34.78, 34.78),
                      bio_1 = c(19.4, 20.1), bio_12 = c(90, 540)),
                  file.path(d, "climate_present_site.tsv"))
            na_ex <- write_tsv(data.table::data.table(sample = "ID003"),
                               file.path(d, "climate_na_excluded.tsv"))
            c(cs, fx_ld_decay_table(d), na_ex)
        },
        args    = function(d, f) c("structure", file.path(d, "pipeline_summary.tsv"),
                                   "3", f[1], "bio_1,bio_12", f[2], "site",
                                   "both", f[3]),
        outputs = function(d) file.path(d, "pipeline_summary.tsv"),
        fails_without_input = FALSE
    ),
    list(
        # The variance-partition table is args[5] of the CLIMATE mode. There is
        # no `varpart` mode — an unknown mode warns and quits 0 (:606-609),
        # which is exactly how a row aimed at one would pass while testing
        # nothing.
        label   = "write_summary.R (climate + varpart)",
        script  = "write_summary.R",
        build   = function(d) {
            inv <- write_tsv(data.table::data.table(
                       predictor = "bio_12", reason = "zero variance across sites"),
                   file.path(d, "climate_invariant_predictors.tsv"))
            conf <- write_tsv(data.table::data.table(
                        confounded = TRUE, r2_climate_on_structure = 0.41),
                    file.path(d, "climate_confounding.tsv"))
            px <- write_tsv(data.table::data.table(
                      variable = c("bio_1", "bio_12"), Px = c(0.38, 0.11)),
                  file.path(d, "px_per_variable.tsv"))
            c(inv, fx_dbmem_diag(d), fx_varpart(d), conf, px)
        },
        args    = function(d, f) c("climate", file.path(d, "pipeline_summary.tsv"),
                                   f[1], f[2], f[3], f[4], f[5]),
        outputs = function(d) file.path(d, "pipeline_summary.tsv"),
        # read_opt() (:261) returns NULL for a declared-but-missing input and the
        # guarded block is skipped: exit 0 with a thinner table. Filed and
        # quarantined, same as the traits row above.
        fails_without_input = FALSE
    ),
    list(
        label   = "pregea_ladder_stats.R",
        script  = "pregea_ladder_stats.R",
        build   = function(d) fx_pvalues(d),
        args    = function(d, pv) c(pv, "lfmm", "K", "3", "0.05", "0.05", "FALSE",
                                    file.path(d, "ladder_stats.tsv")),
        outputs = function(d) file.path(d, "ladder_stats.tsv"),
        fails_without_input = TRUE
    ),
    list(
        label   = "compute_pairwise_ondemand.R",
        script  = "compute_pairwise_ondemand.R",
        build   = function(d) {
            gea  <- file.path(d, "gea_sig.tsv")
            gwas <- file.path(d, "gwas_sig.tsv")
            fx_sig_snps(gea,  "EMMAX", traits = c("bio_1", "bio_2"))
            fx_sig_snps(gwas, "EMMAX", traits = c("height", "flowering_time"))
            c(gea, gwas)
        },
        args    = function(d, f) c(f[1], f[2], "10000", "10000", "1",
                                   file.path(d, "pw_collapsed.tsv"),
                                   file.path(d, "pw_table.tsv")),
        outputs = function(d) file.path(d, c("pw_collapsed.tsv", "pw_table.tsv"))
    ),
    list(
        label   = "compute_pairwise_overlaps.R (both sides)",
        script  = "compute_pairwise_overlaps.R",
        build   = function(d) c(fx_selected_snps(d, "gea_selected.tsv"),
                                fx_selected_snps(d, "gwas_selected.tsv",
                                                 traits = c("height", "flowering_time"))),
        args    = function(d, f) c(f[1], f[2], "10000", "1",
                                   file.path(d, "ov_collapsed.tsv"),
                                   file.path(d, "ov_pairwise.tsv")),
        outputs = function(d) file.path(d, c("ov_collapsed.tsv", "ov_pairwise.tsv"))
    ),
    list(
        # The literal string "NULL" is a supported value for either side
        # (:35) — a one-sided run must still write both tables.
        label   = "compute_pairwise_overlaps.R (GWAS = \"NULL\")",
        script  = "compute_pairwise_overlaps.R",
        build   = function(d) fx_selected_snps(d, "gea_only.tsv"),
        args    = function(d, sel) c(sel, "NULL", "10000", "1",
                                     file.path(d, "one_collapsed.tsv"),
                                     file.path(d, "one_pairwise.tsv")),
        outputs = function(d) file.path(d, c("one_collapsed.tsv", "one_pairwise.tsv"))
    ),
    list(
        label   = "compute_wza.R",
        script  = "compute_wza.R",
        build   = function(d) c(fx_pvalues(d), fx_maf(d)),
        args    = function(d, f) c(f[1], f[2], "NULL", "50000", "10000", "All",
                                   file.path(d, "wza.tsv")),
        outputs = function(d) file.path(d, "wza.tsv"),
        fails_without_input = TRUE
    ),
    list(
        # The poster child for the cwd sandbox: generate_simdata.R:8 falls back
        # to a RELATIVE "data/" when argv is empty, and data/ is gitignored, so
        # a no-arg run would overwrite the working SIMDATA fixtures invisibly.
        # Always pass an explicit dir — with a trailing slash (:301,:441 concat).
        label   = "generate_simdata.R",
        script  = "generate_simdata.R",
        build   = function(d) fx_inter_dir(d, "simdata"),
        args    = function(d, outdir) outdir,
        outputs = function(d) file.path(d, "simdata", c("SIMDATA.vcf", "SIMDATA.gff3"))
    ),

    # ---- plot smokes: exit 0 + declared PNG/SVG exists and is non-empty.
    # Plot CONTENT is never read (repo rule 1). What this catches is the actual
    # failure mode of a plotting wrapper: an argument-order change or an error
    # inside the ggplot chain.

    list(
        label   = "plot_density.R",
        script  = "plot_density.R",
        build   = function(d) c(fx_climate_site(d), fx_inter_dir(d)),
        args    = function(d, f) c(f[1], "bio_1,bio_2", file.path(d, "density.png"), f[2]),
        outputs = function(d) file.path(d, c("density.png", "density.svg"))
    ),
    list(
        label   = "plot_structure.R",
        script  = "plot_structure.R",
        build   = function(d) c(fx_clusters(d, k = 3), fx_inter_dir(d)),
        args    = function(d, f) c(f[1], "3", file.path(d, "structure_K3.png"), f[2]),
        outputs = function(d) file.path(d, c("structure_K3.png", "structure_K3.svg"))
    ),
    list(
        label   = "plot_trait_pairs.R",
        script  = "plot_trait_pairs.R",
        build   = function(d) fx_metadata(d),
        args    = function(d, meta) c(meta, file.path(d, "trait_pairs.png"), "8"),
        outputs = function(d) file.path(d, c("trait_pairs.png", "trait_pairs.svg")),
        fails_without_input = TRUE
    ),
    list(
        label   = "plot_correlation_heatmap.R (one block)",
        script  = "plot_correlation_heatmap.R",
        build   = function(d) c(fx_climate_site(d), fx_inter_dir(d)),
        args    = function(d, f) c(f[1], "NULL", file.path(d, "corr.png"), f[2]),
        outputs = function(d) file.path(d, c("corr.png", "corr.svg"))
    ),
    list(
        label   = "plot_correlation_heatmap.R (two blocks)",
        script  = "plot_correlation_heatmap.R",
        build   = function(d) c(fx_climate_site(d), fx_metadata(d, n = 8),
                                fx_inter_dir(d)),
        args    = function(d, f) c(f[1], f[2], file.path(d, "corr2.png"), f[3],
                                   "Traits x climate"),
        outputs = function(d) file.path(d, c("corr2.png", "corr2.svg"))
    ),
    list(
        # PLOT_DIR-style output: plot_manhattan.R:156-208 derives every basename
        # itself, so assert a count of matching files rather than re-deriving
        # production names in the test.
        label   = "plot_manhattan.R",
        script  = "plot_manhattan.R",
        build   = function(d) {
            dir.create(file.path(d, "plots"), showWarnings = FALSE)
            fx_pvalues(d)
        },
        args    = function(d, pv) c(pv, "bonf_0.05", "3", "EMMAX", "bio_1",
                                    file.path(d, "plots"), "bio_1,bio_2"),
        outputs = function(d) character(0),
        out_min = function(d) list(dir = file.path(d, "plots"),
                                   pattern = "\\.png$", n = 1)
    ),
    list(
        label   = "plot_manhattan_combined.R",
        script  = "plot_manhattan_combined.R",
        build   = function(d) {
            dir.create(file.path(d, "plots_comb"), showWarnings = FALSE)
            paste0("EMMAX:bonf_0.05:", fx_pvalues(d))
        },
        args    = function(d, files_str) c(files_str, "bio_1,bio_2", "3",
                                           file.path(d, "plots_comb")),
        outputs = function(d) character(0),
        out_min = function(d) list(dir = file.path(d, "plots_comb"),
                                   pattern = "\\.png$", n = 1)
    ),
    list(
        label   = "plot_miami.R",
        script  = "plot_miami.R",
        build   = function(d) {
            dir.create(file.path(d, "plots_miami"), showWarnings = FALSE)
            gea  <- fx_pvalues(d, "gea_pv.tsv")
            gwas <- data.table::data.table(
                SNPID = sprintf("snp%02d", 1:20), chr = rep(c("1", "2"), each = 10),
                pos = seq_len(20) * 10000L,
                height = c(1e-9, runif(19, 0.05, 1)),
                flowering_time = c(runif(19, 0.05, 1), 1e-9))
            c(paste0("EMMAX:bonf_0.05:", gea),
              paste0("EMMAX:bonf_0.05:", write_tsv(gwas, file.path(d, "gwas_pv.tsv"))))
        },
        args    = function(d, f) c(f[1], f[2], "bio_1,bio_2",
                                   "height,flowering_time", "3",
                                   file.path(d, "plots_miami")),
        outputs = function(d) character(0),
        out_min = function(d) list(dir = file.path(d, "plots_miami"),
                                   pattern = "\\.png$", n = 1)
    ),
    list(
        label   = "plot_pca_structure.R",
        script  = "plot_pca_structure.R",
        build   = function(d) c(fx_clusters(d, k = 3), fx_projections(d),
                                fx_eigenvalues(d), fx_inter_dir(d)),
        args    = function(d, f) c(f[1], f[2], f[3], "3",
                                   file.path(d, "pca_structure_K3.png"), f[4]),
        outputs = function(d) file.path(d, c("pca_structure_K3.png",
                                             "pca_structure_K3.svg"))
    ),
    list(
        label   = "plot_pregea_screeplot.R",
        script  = "plot_pregea_screeplot.R",
        build   = function(d) c(fx_eigenvalues(d), fx_inter_dir(d)),
        args    = function(d, f) c(f[1], "3", "2,3,4,5", "5",
                                   file.path(d, "scree.png"),
                                   file.path(d, "scree.tsv"), f[2]),
        outputs = function(d) file.path(d, c("scree.png", "scree.tsv")),
        fails_without_input = TRUE
    ),
    # check_invariants.R is the one wrapper that is a VALIDATOR rather than a
    # producer: it writes nothing (outputs = character(0)) and its exit status IS
    # the result. Its own header (:12-18) names the file -> checker shaping layer
    # as the untested part, and the two misreads found while writing it both
    # lived there. The checkers themselves are covered by test-invariants.R
    # against hand-built fixtures; these two rows cover the layer in between.
    list(
        label   = "check_invariants.R (clean tree)",
        script  = "check_invariants.R",
        build   = function(d) fx_results_tree(d, violations = FALSE),
        args    = function(d, res) c(res, "--modules", "GEA"),
        outputs = function(d) character(0),
        # argv[1] is normalizePath(mustWork = TRUE) at :41, so a missing results
        # dir is a hard error.
        fails_without_input = TRUE
    ),
    list(
        label   = "check_invariants.R (violating tree exits 1)",
        script  = "check_invariants.R",
        build   = function(d) fx_results_tree(d, violations = TRUE),
        args    = function(d, res) c(res, "--modules", "GEA"),
        outputs = function(d) character(0),
        # The load-bearing row. Exit 1 proves the shaping layer actually reached
        # the checkers: a schema misread would silently produce an empty or
        # NULL table, every checker would return no_violations(), and the run
        # would exit 0 while asserting nothing. Both mutations
        # (pvalue_out_of_range, chromosome_name_not_normalized) are severity
        # "error" — quit(status = 1) at :272-274 counts only those.
        expect_status = 1L
    )
)

# --------------------------------------------------------------------- tests

# One LOUD test for the mount, so a broken mount is a single red rather than 40
# silent skips. The per-row skip_if_not below then stays as the graceful path for
# anyone running the file outside the container.
test_that("the repo is mounted where the wrappers expect it", {
    expect_true(dir.exists(SCRIPTS),
                info = paste0("scripts/ not found at ", SCRIPTS,
                              " — run this suite with -v $PWD:/pipeline"))
})

for (spec in WRAPPERS) {
    local({
        s <- spec
        test_that(paste0("wrapper runs and writes its outputs: ", s$label), {
            skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
            skip_if_not(nzchar(Sys.which("Rscript")), "Rscript not on PATH")
            expect_wrapper_ok(s)
        })

        if (isTRUE(s$fails_without_input)) {
            test_that(paste0("wrapper fails loudly on a missing input: ", s$label), {
                skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
                expect_wrapper_fails_without_input(s)
            })
        }
    })
}

test_that("no spec invokes a wrapper that rewrites the tracked fixtures", {
    # add_related_samples.R and add_pregea_sites.R take NO argv and
    # read-modify-write data/SIMDATA.* in place. data/ is gitignored, so nothing
    # downstream would notice. Keep this a mechanism, not a review convention.
    scripts_used <- vapply(WRAPPERS, function(s) s$script, character(1))
    expect_length(intersect(scripts_used, WRAPPER_DENYLIST), 0L)
})

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
