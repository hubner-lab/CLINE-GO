# Phase 3 (table-changing fixes) — CONTENT assertions on the scripts, quick tier.
#
# Each block runs a script on a fixture that takes the branch its fix changed; every one of
# them fails on the script as it was before the fix.
#   create_regions.R, compute_pairwise_overlaps.R      dfb365  single-method combine_method
#   compute_pairwise_overlaps.R, _ondemand.R           ad7fce  a shared locus is counted once
#   pregea_recommend.R                                 3602c5  the RDA fallback says it is one
#   gradient_forest_model.R                            5ae1e3  random model drawn outside the set
#   write_summary.R                                    2e7ba2  n_climate_variables = predictors
# The library halves are in test-regions.R, test-region_distance.R, test-combine_sigsnps.R,
# test-gff_parsing.R and test-pregea_calibration.R.

source(file.path(getOption("clinego.repo_root", "/pipeline"),
                 "tests", "lib", "wrapper_harness.R"))

expect_ok <- function(res) expect_identical(res$status, 0L, info = res$output)
read_chr <- function(path) data.table::fread(path, colClasses = c(chr = "character"))

# ---- dfb365: a single-method combine_method selects by that method's calls alone ----------

# combine_method EMMAX: s1 and s2 are EMMAX's bio_1 hits, s3 its bio_2 hit. LFMM also called
# s1 for bio_2, a call that is no part of the selection.
fx_single_method_selection <- function(d) {
    sel <- write_tsv(data.table::data.table(
               SNPID = c("s1", "s2", "s3"), chr = c("1", "1", "2"),
               pos = c(10000L, 20000L, 110000L),
               EMMAX = c("bio_1", "bio_1", "bio_2"), LFMM = c("bio_2", "", "bio_2"),
               min_pvalue = c(1e-9, 1e-8, 1e-7)),
           file.path(d, "selected_snps.tsv"))
    per <- write_tsv(data.table::data.table(
               SNPID = c("s1", "s2", "s3"), chr = c("1", "1", "2"),
               pos = c(10000L, 20000L, 110000L), trait = c("bio_1", "bio_1", "bio_2"),
               min_pvalue = c(1e-9, 1e-8, 1e-7)),
           file.path(d, "selected_snps_per_trait.tsv"))
    list(sel = sel, per = per)
}

run_create_regions <- function(d, f, combine_method, tag) {
    out <- file.path(d, paste0(tag, c("_per_trait.tsv", "_combined.tsv")))
    res <- run_wrapper("create_regions.R", c(f$sel, "10000", out[1], out[2], "NULL", "critical",
                                             "All", f$per, combine_method))
    expect_ok(res)
    list(per_trait = read_chr(out[1]), combined = read_chr(out[2]), log = res$output)
}

test_that("create_regions.R takes a SNP's traits from the selecting method alone", {
    d <- withr::local_tempdir()
    f <- fx_single_method_selection(d)
    one <- run_create_regions(d, f, "EMMAX", "emmax")
    expect_identical(one$per_trait[, paste(chr, trait)], c("1 bio_1", "2 bio_2"))
    expect_identical(one$per_trait$methods, c("EMMAX", "EMMAX"))
    expect_true(all(is.na(one$per_trait$other_traits) | one$per_trait$other_traits == ""))
    expect_identical(one$combined$traits, c("bio_1", "bio_2"))
    expect_false("LFMM_snps" %in% names(one$combined))

    # A strategy is no method name: every column still gives traits, as before.
    all <- run_create_regions(d, f, "All", "all")
    expect_identical(all$per_trait[, paste(chr, trait)], c("1 bio_1", "1 bio_2", "2 bio_2"))
})

test_that("compute_pairwise_overlaps.R takes a SNP's traits from the selecting method alone", {
    d <- withr::local_tempdir()
    f <- fx_single_method_selection(d)
    run <- function(gea_combine, tag) {
        out <- file.path(d, paste0(tag, c("_collapsed.tsv", "_pairwise.tsv")))
        expect_ok(run_wrapper("compute_pairwise_overlaps.R",
                              c(f$sel, "NULL", "10000", "1", out[1], out[2], gea_combine, "All")))
        read_chr(out[1])
    }
    expect_identical(run("EMMAX", "emmax")[trait == "bio_2", SNPID], "s3")
    expect_setequal(run("All", "all")[trait == "bio_2", SNPID], c("s1", "s3"))
})

# ---- ad7fce: the overlap score counts a shared locus once ---------------------------------

# bio_1 and bio_2 share s1-s3 on chr 1 (one window, which holds them) and lie 1 kb apart, two
# SNPs each, on chr 2 (a window with no shared SNP). exact 3 + windows 2 - 1 holding exact = 4;
# the old exact + windows gave 5.
ad7fce_snps <- data.table::data.table(
    SNPID = c("s1", "s2", "s3", "t1", "t2", "u1", "u2"),
    chr   = c("1", "1", "1", "2", "2", "2", "2"),
    pos   = c(10000L, 12000L, 14000L, 500000L, 501000L, 502000L, 503000L),
    trait = c("bio_1,bio_2", "bio_1,bio_2", "bio_1,bio_2", "bio_1", "bio_1", "bio_2", "bio_2"))

expect_ad7fce_row <- function(pw) {
    expect_identical(nrow(pw), 1L)
    expect_identical(pw$exact_matches, 3L)
    expect_identical(pw$window_overlaps, 2L)
    expect_identical(pw$overlap_score, 4L)
}

test_that("compute_pairwise_overlaps.R counts each shared SNP once in the overlap score", {
    d <- withr::local_tempdir()
    sel <- write_tsv(ad7fce_snps[, .(SNPID, chr, pos, EMMAX = trait, min_pvalue = 1e-8)],
                     file.path(d, "gea_selected.tsv"))
    out <- file.path(d, c("collapsed.tsv", "pairwise.tsv"))
    expect_ok(run_wrapper("compute_pairwise_overlaps.R",
                          c(sel, "NULL", "10000", "2", out[1], out[2], "All", "All")))
    expect_ad7fce_row(data.table::fread(out[2]))
})

test_that("compute_pairwise_ondemand.R, the app's copy, scores alike", {
    d <- withr::local_tempdir()
    long <- ad7fce_snps[, .(trait = strsplit(trait, ",", fixed = TRUE)[[1]]),
                        by = .(SNPID, chr, pos)]
    sig <- write_tsv(long[, .(SNPID, chr, pos, pvalue = 1e-8, method = "EMMAX", trait)],
                     file.path(d, "gea_sig.tsv"))
    out <- file.path(d, c("collapsed.tsv", "table.tsv"))
    expect_ok(run_wrapper("compute_pairwise_ondemand.R",
                          c(sig, "NULL", "10000", "10000", "2", out[1], out[2])))
    expect_ad7fce_row(data.table::fread(out[2]))
})

# ---- 3602c5: the RDA fallback rung is labelled as one -------------------------------------

recommend <- function(d, anova_p) {
    w <- function(dt, name) write_tsv(dt, file.path(d, name))
    f <- c(w(data.table::data.table(trait = "__pooled__", rung_value_num = 2:3, hits_qval = 10,
                                    hist_flatness_ks = c(0.2, 0.1)), "lfmm.tsv"),
           w(data.table::data.table(trait = "__pooled__", rung_value_num = 0:1,
                                    lambda_gc = c(1.3, 1.05)), "emmax.tsv"),
           w(data.table::data.table(status = "ok", condition_pcs = 0:2,
                                    r2_adj = c(0.10, 0.14, 0.13), max_vif = c(1.5, 2, 3),
                                    n_aliased = 0L, anova_full_p = anova_p,
                                    n_axes_sig = 2L), "rda.tsv"),
           w(data.table::data.table(predictor = c("bio_1", "bio_2"), action = "kept"),
             "collin.tsv"),
           w(data.table::data.table(condition_pcs = rep(0:2, each = 2), axis = rep(1:2, 3),
                                    axis_eig = 0.3, axis_p = 0.01), "axis.tsv"))
    out <- file.path(d, "recs.tsv")
    expect_ok(run_wrapper("pregea_recommend.R",
                          c(f, "3", "0.15", "0.90", "10", "0.05", out)))
    data.table::fread(out, colClasses = "character")[method == "RDA" & param == "condition_pcs"]
}

test_that("pregea_recommend.R says when no RDA rung passed and the fallback was taken", {
    d <- withr::local_tempdir()
    # SIMDATA's case: no rung reaches anova_full_p < axis_alpha.
    fb <- recommend(d, anova_p = c(0.475, 0.42, 0.51))
    expect_identical(fb$recommended_value, "1")                 # the max-R2adj rung
    expect_match(fb$rule, "^max R2adj \\(no rung passed VIF/aliasing/significance")
    expect_match(fb$note, "no rung passed the rule; smallest full-model anova p 0.42")

    ok <- recommend(d, anova_p = 0.001)
    expect_identical(ok$rule, paste0("smallest rung passing VIF/aliasing/significance, ",
                                     "R2adj within 1 SE of ladder max"))
    expect_true(is.na(ok$note) || ok$note == "")
})

# ---- 5ae1e3: the random GF model is drawn outside the SNP set -----------------------------

gf_random <- function(w, out) {
    run_wrapper("gradient_forest_model.R",
                c(w$lfmm, w$snps, w$vcfsnp, w$removed, w$meta, w$site_pres, "bio_1,bio_12",
                  "100", "0.5", "without", "random", out, "NULL", "NULL"))
}

test_that("gradient_forest_model.R draws the random model from SNPs outside the set", {
    d <- withr::local_tempdir()
    w <- fx_offset_world(d)            # 30 SNPs (columns V1-V30); the set is the first 15
    out <- file.path(d, "setA_nospatial", "random_model.qs")
    dir.create(dirname(out))
    res <- gf_random(w, out)
    expect_ok(res)
    gf <- qs::qread(out)
    expect_s3_class(gf, "gradientForest")
    expect_setequal(colnames(gf$Y), paste0("V", 16:30))         # 15 drawn from the 15 outside
    expect_match(res$output, "Using 15 random SNPs for the 15-SNP adaptive set, drawn from the 15 SNPs outside it")

    # A set holding every SNP leaves nothing to draw from.
    vs <- data.table::fread(w$vcfsnp, header = FALSE)
    write_tsv(data.table::data.table(SNPID = paste0(vs$V1, ":", vs$V2)), w$snps)
    res <- gf_random(w, file.path(d, "setA_nospatial", "random_all.qs"))
    expect_false(identical(res$status, 0L))
    expect_match(res$output, "no SNP is left to draw a random model from")
})

# ---- 2e7ba2: n_climate_variables counts the predictor columns -----------------------------

test_that("write_summary.R counts the site table's predictor columns, not ncol - 4", {
    d <- withr::local_tempdir()
    # climate_present_site.tsv as download_climate_present.R writes it: sample + predictors.
    cs <- write_tsv(data.table::data.table(sample = c("ID001", "ID002"), bio_1 = c(19.4, 20.1),
                                           bio_2 = c(9, 11), bio_12 = c(90, 540)),
                    file.path(d, "climate_present_site.tsv"))
    out <- file.path(d, "pipeline_summary.tsv")
    expect_ok(run_wrapper("write_summary.R",
                          c("structure", out, "3", cs, "bio_1,bio_2,bio_12", "NULL", "NULL",
                            "NULL", "NULL", "NULL")))
    s <- data.table::fread(out, colClasses = "character")
    expect_identical(s[step == "structure" & metric == "n_climate_variables", value], "3")
    expect_identical(s[step == "structure" & metric == "climate_predictors", value],
                     "bio_1,bio_2,bio_12")
})
