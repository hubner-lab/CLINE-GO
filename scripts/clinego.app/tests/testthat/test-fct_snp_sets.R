# R/fct_snp_sets.R — the curated-SNP-set store the Shiny app writes and the
# Snakemake pipeline resolves by glob.
#
# It is a CONTRACT with the pipeline, not just app state: promote_snp_set.R:81
# asserts the four TSV columns written here, and common.smk:1423-1428 discovers
# sets with glob("{store}*/selected_snps.tsv"). Those cross-side assertions live
# in tests/testthat/test-equivalence-app-pipeline.R; this file pins the app half.
#
# The fixture shape matters and is easy to get wrong: save_snp_set() aggregates
# min(min_pvalue), so its input must be the EIGHT-column success shape of
# compute_interactive_sigsnps(), not the seven-column .empty_sigsnps_assoc().
# combine_sigsnps() stamps min_pvalue at fct_combine.R:135; the empty path never
# does.
#
# Writes land under project_base(), redirected to a session tempdir by
# helper-pipeline-path.R. Project names are unique per INVOCATION, not merely per
# block — CLINEGO_TEST_ROOT is stable for the whole session, so a fixed name
# makes a test that asserts absence unrepeatable.

# The 8-column long shape compute_interactive_sigsnps() returns on success.
fx_sigsnps <- function(snpid = c("1:100", "1:100", "2:300"),
                       chr   = c("1", "1", "2"),
                       pos   = c(100L, 100L, 300L),
                       pvalue = c(1e-8, 3e-9, 5e-7),
                       method = c("EMMAX", "LFMM", "EMMAX"),
                       trait  = c("bio_1", "bio_1", "bio_2")) {
    dt <- data.table::data.table(
        SNPID = snpid, chr = chr, pos = pos, pvalue = pvalue,
        method = method, trait = trait, region_id = NA_character_)
    dt[, min_pvalue := min(pvalue, na.rm = TRUE), by = "SNPID"]
    dt[]
}

new_project <- function(tag) basename(tempfile(paste0("SS_", tag, "_")))

# ---------------------------------------------------------------- manifest reads

test_that("read_snp_sets_manifest returns an empty list when nothing is stored", {
    expect_identical(read_snp_sets_manifest(new_project("EMPTY")), list())
})

test_that("read_snp_sets_manifest returns an empty list for a NULL project", {
    expect_identical(read_snp_sets_manifest(NULL), list())
})

test_that("read_snp_sets_manifest degrades to an empty list on corrupt JSON", {
    project <- new_project("CORRUPT")
    path <- snp_sets_manifest_path(project)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("{ not json", path)
    expect_identical(read_snp_sets_manifest(project), list())
})

test_that("list_snp_sets returns a typed zero-row frame when the manifest is empty", {
    got <- list_snp_sets(new_project("LIST_EMPTY"))
    expect_identical(nrow(got), 0L)
    expect_identical(names(got), c("name", "n_snps"))
    expect_type(got$name, "character")
    expect_type(got$n_snps, "integer")
})

test_that("set_exists is FALSE for an unknown project or name", {
    expect_false(set_exists(new_project("NOPE"), "anything"))
})

# ---------------------------------------------------------------- save

test_that("save_snp_set writes the four contract columns, unique by SNPID", {
    project <- new_project("SAVE")
    n <- save_snp_set(project, "setA", fx_sigsnps(), list(source_module = "GEA"))
    expect_identical(n, 2L)   # 1:100 appears twice in the long input

    got <- data.table::fread(snp_set_path(project, "setA"), sep = "\t")
    expect_identical(names(got), c("SNPID", "chr", "pos", "min_pvalue"))
    expect_identical(nrow(got), 2L)
    expect_setequal(got$SNPID, c("1:100", "2:300"))
})

test_that("save_snp_set keeps the MINIMUM p across a SNP's rows", {
    project <- new_project("MINP")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    got <- data.table::fread(snp_set_path(project, "setA"), sep = "\t")
    # 1:100 was seen at 1e-8 (EMMAX) and 3e-9 (LFMM).
    expect_equal(got[SNPID == "1:100"]$min_pvalue, 3e-9)
    expect_equal(got[SNPID == "2:300"]$min_pvalue, 5e-7)
})

test_that("save_snp_set orders rows by chr then pos", {
    project <- new_project("ORDER")
    save_snp_set(project, "setA", fx_sigsnps(
        snpid = c("2:50", "1:900", "1:100"), chr = c("2", "1", "1"),
        pos = c(50L, 900L, 100L), pvalue = c(1e-7, 1e-8, 1e-9),
        method = "EMMAX", trait = "bio_1"), list())
    got <- data.table::fread(snp_set_path(project, "setA"), sep = "\t",
                             colClasses = c(chr = "character"))
    expect_identical(got$SNPID, c("1:100", "1:900", "2:50"))
})

test_that("save_snp_set registers the set in the manifest with its params", {
    project <- new_project("MANIFEST")
    save_snp_set(project, "setA", fx_sigsnps(),
                 list(source_module = "GEA", strategy = "Union",
                      threshold_type = "bonf", threshold_value = 0.05))
    man <- read_snp_sets_manifest(project)
    expect_length(man, 1L)
    expect_identical(man[[1]]$name, "setA")
    expect_identical(man[[1]]$n_snps, 2L)
    expect_identical(man[[1]]$source_module, "GEA")
    expect_identical(man[[1]]$strategy, "Union")
    expect_true(nzchar(man[[1]]$created))
})

test_that("set_exists and list_snp_sets see a saved set", {
    project <- new_project("EXISTS")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    expect_true(set_exists(project, "setA"))
    expect_false(set_exists(project, "setB"))
    got <- list_snp_sets(project)
    expect_identical(got$name, "setA")
    expect_identical(got$n_snps, 2L)
})

test_that("save_snp_set leaves no tempfile beside the TSV or the manifest", {
    project <- new_project("TMP")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    expect_identical(list.files(dirname(snp_set_path(project, "setA")),
                                pattern = "^file.*\\.tsv$"), character(0))
    expect_identical(list.files(dirname(snp_sets_manifest_path(project)),
                                pattern = "^file.*\\.json$"), character(0))
})

# ---------------------------------------------------------------- upsert

test_that("re-saving a set replaces its manifest entry rather than duplicating it", {
    project <- new_project("UPSERT")
    save_snp_set(project, "setA", fx_sigsnps(), list(strategy = "Union"))
    save_snp_set(project, "setA",
                 fx_sigsnps(snpid = "9:1", chr = "9", pos = 1L, pvalue = 1e-6,
                            method = "EMMAX", trait = "bio_1"),
                 list(strategy = "Cross-method"))
    man <- read_snp_sets_manifest(project)
    expect_length(man, 1L)
    expect_identical(man[[1]]$strategy, "Cross-method")
    expect_identical(man[[1]]$n_snps, 1L)
})

test_that("re-saving MOVES a set to the end of manifest order", {
    # Remove-then-append, not replace-in-place (:83-85). Manifest order is what
    # list_snp_sets returns and therefore what the picker shows, so re-saving an
    # existing set reorders the UI.
    project <- new_project("REORDER")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    save_snp_set(project, "setB", fx_sigsnps(), list())
    expect_identical(list_snp_sets(project)$name, c("setA", "setB"))

    save_snp_set(project, "setA", fx_sigsnps(), list())
    expect_identical(list_snp_sets(project)$name, c("setB", "setA"))
})

test_that("saving a second set leaves the first one's rows untouched", {
    project <- new_project("TWO")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    save_snp_set(project, "setB",
                 fx_sigsnps(snpid = "9:1", chr = "9", pos = 1L, pvalue = 1e-6,
                            method = "EMMAX", trait = "bio_1"), list())
    a <- data.table::fread(snp_set_path(project, "setA"), sep = "\t")
    expect_identical(nrow(a), 2L)
    expect_identical(list_snp_sets(project)$name, c("setA", "setB"))
})

# ---------------------------------------------------------------- delete

test_that("delete_snp_set removes the set directory and its manifest entry", {
    project <- new_project("DELETE")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    save_snp_set(project, "setB", fx_sigsnps(), list())

    delete_snp_set(project, "setA", remove_gf_results = FALSE)
    expect_false(file.exists(snp_set_path(project, "setA")))
    expect_false(set_exists(project, "setA"))
    expect_true(set_exists(project, "setB"))
    expect_identical(list_snp_sets(project)$name, "setB")
})

test_that("delete_snp_set is a no-op for a name that was never saved", {
    project <- new_project("DEL_MISSING")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    expect_true(delete_snp_set(project, "ghost", remove_gf_results = FALSE))
    expect_true(set_exists(project, "setA"))
})

test_that("delete_snp_set does NOT over-delete a prefix sibling", {
    # The documented worry (:94-97): deleting "foo" must not take "foo_bar".
    # The suffix set is exact — {name}, {name}_spatial, {name}_nospatial — so
    # this holds.
    project <- new_project("PREFIX")
    save_snp_set(project, "foo", fx_sigsnps(), list())
    save_snp_set(project, "foo_bar", fx_sigsnps(), list())

    delete_snp_set(project, "foo", remove_gf_results = TRUE)
    expect_false(set_exists(project, "foo"))
    expect_true(set_exists(project, "foo_bar"))
    expect_true(file.exists(snp_set_path(project, "foo_bar")))
})

test_that("delete_snp_set removes the maladaptation results for its own suffixes", {
    project <- new_project("GF")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    targets <- c(
        file.path(mod_path(project, MOD_MALAD, "tables", "gradient_forest"), "setA_spatial"),
        file.path(mod_path(project, MOD_MALAD, "plots",  "gradient_forest"), "setA_nospatial"),
        file.path(mod_path(project, MOD_INTER, "geometric_offset"), "setA"))
    for (t in targets) dir.create(t, recursive = TRUE, showWarnings = FALSE)
    other <- file.path(mod_path(project, MOD_MALAD, "tables", "gradient_forest"), "setB_spatial")
    dir.create(other, recursive = TRUE, showWarnings = FALSE)

    delete_snp_set(project, "setA", remove_gf_results = TRUE)
    for (t in targets) expect_false(dir.exists(t))
    expect_true(dir.exists(other))
})

test_that("delete_snp_set leaves maladaptation results alone when asked not to", {
    project <- new_project("GF_KEEP")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    keep <- file.path(mod_path(project, MOD_MALAD, "tables", "gradient_forest"), "setA_spatial")
    dir.create(keep, recursive = TRUE, showWarnings = FALSE)

    delete_snp_set(project, "setA", remove_gf_results = FALSE)
    expect_true(dir.exists(keep))
})

test_that("delete_snp_set returns TRUE whether or not anything was deleted", {
    # No success signal: a caller cannot tell a real deletion from a no-op.
    project <- new_project("RET")
    expect_true(delete_snp_set(project, "never-existed", remove_gf_results = FALSE))
})
