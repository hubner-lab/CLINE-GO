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
# combine_sigsnps() stamps min_pvalue at lib/snp_set_select.R:118; the empty path never
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
    expect_identical(names(got), c("name", "n_snps", "origin"))
    expect_type(got$name, "character")
    expect_type(got$n_snps, "integer")
    expect_type(got$origin, "character")
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

# ---------------------------------------------------------------- order, names

test_that("save_snp_set breaks a chr/pos tie by SNPID, as the pipeline build does", {
    # collapse_snp_set() (lib/snp_set_select.R) is shared with build_snp_set.R:
    # without the SNPID tiebreak two IDs at one position kept their input order.
    project <- new_project("TIE")
    save_snp_set(project, "setA", fx_sigsnps(
        snpid = c("1:100_b", "1:100_a"), chr = "1", pos = 100L,
        pvalue = c(1e-7, 1e-8), method = "EMMAX", trait = "bio_1"), list())
    got <- data.table::fread(snp_set_path(project, "setA"), sep = "\t")
    expect_identical(got$SNPID, c("1:100_a", "1:100_b"))
})

test_that("save_snp_set refuses the gea. / gwas. prefixes of config sets", {
    project <- new_project("PREFIX_REFUSE")
    for (nm in c("gea.x", "gwas.x", "GEA.x")) {
        expect_error(save_snp_set(project, nm, fx_sigsnps(), list()),
                     "reserved for sets declared in config", info = nm)
    }
    expect_false(dir.exists(file.path(snp_sets_dir(project), "gea.x")))
    # Only the prefix is reserved: a dot elsewhere is still a valid app name.
    expect_silent(save_snp_set(project, "my.gea.x", fx_sigsnps(), list()))
})

test_that("snp_set_name_problem names each reason a name is refused", {
    expect_identical(snp_set_name_problem(""), "Enter a name.")
    expect_match(snp_set_name_problem("a b"), "Only letters")
    expect_match(snp_set_name_problem("gwas.t"), "reserved")
    expect_null(snp_set_name_problem("bonf05.union"))
})

# ---------------------------------------------------------------- record

fx_registry <- list(
    EMMAX = list(adjust = "bonf", threshold = "0.05", family = "univariate_pvalue"),
    LFMM  = list(adjust = "bonf", threshold = "0.05", family = "univariate_pvalue"),
    RDA   = list(adjust = "bonf", threshold = "0.01", family = "multivariate_pvalue"))

fx_record <- function(pairs = c("bio_2::LFMM", "bio_1::EMMAX", "climate_multivariate::RDA"),
                      overrides = list(), strategy = "All", regime = "snp", name = "setA") {
    snp_set_record_for_save(
        name = name, source_module = "GEA", pairs = pairs, overrides = overrides,
        master_type = "bonf", master_value = 0.05, registry_defaults = fx_registry,
        strategy = strategy, clumping_distance = 100000L, regime = regime, k_best = 3L)
}

test_that("the record holds the sorted pairs and each cell's effective rule", {
    rec <- fx_record(overrides = list("bio_2::LFMM" = list(type = "qval", value = 0.1)))
    expect_identical(rec$pairs, c("bio_1::EMMAX", "bio_2::LFMM", "climate_multivariate::RDA"))
    expect_identical(rec$cells[["bio_1::EMMAX"]], list(type = "bonf", value = 0.05))
    expect_identical(rec$cells[["bio_2::LFMM"]], list(type = "qval", value = 0.1))
    # RDA is pinned to its registry rule, exactly as the matrix shows it.
    expect_identical(rec$cells[["climate_multivariate::RDA"]], list(type = "bonf", value = 0.01))
    expect_identical(rec$strategy, "Union")   # the app's default "All", normalised
    expect_identical(rec$k_best, 3L)
    expect_type(rec$config_snippet, "character")
})

test_that("save_snp_set writes the record and adds pairs and k_best to the manifest", {
    project <- new_project("RECORD")
    rec <- fx_record(pairs = "bio_1::EMMAX")
    save_snp_set(project, "setA", fx_sigsnps(), list(source_module = "GEA"), record = rec)

    on_disk <- yaml::read_yaml(snp_set_record_path(project, "setA"))
    expect_identical(on_disk$origin, "shiny")
    expect_identical(on_disk$n_snps, 2L)
    expect_identical(on_disk$pairs, "bio_1::EMMAX")
    rec_txt <- paste(readLines(snp_set_record_path(project, "setA")), collapse = "\n")
    expect_match(rec_txt, "pairs:\n- bio_1::EMMAX", fixed = TRUE)   # a sequence, not a scalar
    expect_true(nzchar(on_disk$created))
    # The shared reader accepts it: build_snp_set.R takes a record as its DECL.
    expect_identical(read_snp_set_decl(snp_set_record_path(project, "setA"))$pairs,
                     "bio_1::EMMAX")

    raw <- paste(readLines(snp_sets_manifest_path(project)), collapse = "\n")
    expect_match(raw, '"pairs": \\[\\s*"bio_1::EMMAX"\\s*\\]')
    expect_identical(read_snp_sets_manifest(project)[[1]]$k_best, 3L)
})

test_that("re-saving without a record removes the old one rather than leave it stale", {
    project <- new_project("STALE")
    save_snp_set(project, "setA", fx_sigsnps(), list(), record = fx_record())
    expect_true(file.exists(snp_set_record_path(project, "setA")))
    save_snp_set(project, "setA", fx_sigsnps(), list())
    expect_false(file.exists(snp_set_record_path(project, "setA")))
})

# ---------------------------------------------------------------- list: both origins

fx_config_set <- function(project, name, n) {
    dir <- file.path(snp_sets_dir(project), name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    tsv <- data.table::data.table(SNPID = character(), chr = character(),
                                  pos = integer(), min_pvalue = numeric())
    if (n > 0) tsv <- data.table::data.table(SNPID = paste0("1:", seq_len(n)), chr = "1",
                                             pos = seq_len(n), min_pvalue = 1e-6)
    data.table::fwrite(tsv, file.path(dir, "selected_snps.tsv"), sep = "\t")
    write_snp_set_record(list(
        name = name, origin = "config", source_module = "GEA", regime = "snp",
        strategy = "Union", clumping_distance = 100000L, pairs = "bio_1::EMMAX",
        cells = list("bio_1::EMMAX" = list(type = "bonf", value = 0.01)),
        k_best = 3L, n_snps = n, created = "2026-10-14T00:00:00"),
        snp_set_record_path(project, name))
}

test_that("list_snp_sets merges app sets and config sets", {
    project <- new_project("MERGE")
    save_snp_set(project, "mysel", fx_sigsnps(), list())
    fx_config_set(project, "gea.strict", 5L)
    fx_config_set(project, "gea.empty", 0L)
    # A directory without a TSV is not a set (the pipeline's glob agrees).
    dir.create(file.path(snp_sets_dir(project), "half"), recursive = TRUE)

    got <- list_snp_sets(project)
    expect_identical(got$name, c("mysel", "gea.empty", "gea.strict"))   # manifest first
    expect_identical(got$origin, c("shiny", "config", "config"))
    expect_identical(got$n_snps, c(2L, 0L, 5L))
    expect_true(set_exists(project, "gea.strict"))
})

test_that("a manifest entry whose directory is gone is no longer listed", {
    project <- new_project("GONE")
    save_snp_set(project, "setA", fx_sigsnps(), list())
    unlink(file.path(snp_sets_dir(project), "setA"), recursive = TRUE)
    expect_identical(nrow(list_snp_sets(project)), 0L)
})

test_that("the picker offers only sets with SNPs, labelled by origin", {
    project <- new_project("PICK")
    save_snp_set(project, "mysel", fx_sigsnps(), list())
    fx_config_set(project, "gea.strict", 5L)
    fx_config_set(project, "gea.empty", 0L)
    sets <- pickable_snp_sets(project)
    expect_identical(sets$name, c("mysel", "gea.strict"))
    expect_identical(snp_set_label(sets),
                     c("mysel (2 SNPs) · app", "gea.strict (5 SNPs) · config"))
})

# ---------------------------------------------------------------- config snippet

snippet_entry <- function(sn, group = "GEA") yaml::yaml.load(sn$yaml)[[group]]$snp_sets[[1]]

test_that("a non-rectangular selection is written as explicit pairs", {
    sn <- snp_set_config_snippet(fx_record())
    e <- snippet_entry(sn)
    expect_null(sn$note)
    expect_identical(e$name, "setA")
    expect_identical(unlist(e$pairs),
                     c("bio_1::EMMAX", "bio_2::LFMM", "climate_multivariate::RDA"))
    expect_null(e$traits)
    expect_identical(e$threshold_type, "bonf")
    expect_identical(e$strategy, "Union")
    expect_identical(e$clumping_distance, 100000L)
})

test_that("a full rectangle is written as traits and methods, even with one of each", {
    e <- snippet_entry(snp_set_config_snippet(fx_record(
        pairs = c("bio_1::EMMAX", "bio_2::EMMAX", "bio_1::LFMM", "bio_2::LFMM"))))
    expect_identical(e$methods, c("EMMAX", "LFMM"))
    expect_identical(e$traits, c("bio_1", "bio_2"))
    expect_null(e$pairs)
    one <- snp_set_config_snippet(fx_record(pairs = "bio_1::EMMAX"))
    expect_match(one$yaml, "traits:\n    - bio_1\n", fixed = TRUE)   # a sequence
})

test_that("overrides cover exactly the cells config would otherwise get wrong", {
    rec <- fx_record(overrides = list(
        "bio_2::LFMM"               = list(type = "qval",   value = 0.1),
        "bio_1::EMMAX"              = list(type = "bonf",   value = 0.05),  # = master
        "climate_multivariate::RDA" = list(type = "bonf",   value = 0.05),  # = master, not pin
        "bio_3::EMMAX"              = list(type = "custom", value = 1e-3))) # not selected
    e <- snippet_entry(snp_set_config_snippet(rec, fx_registry))
    expect_setequal(names(e$overrides), c("bio_2::LFMM", "climate_multivariate::RDA"))
    expect_identical(e$overrides[["bio_2::LFMM"]], list(type = "qval", value = 0.1))
    expect_identical(e$overrides[["climate_multivariate::RDA"]],
                     list(type = "bonf", value = 0.05))
    # Without the registry the pin is unknown, so the RDA cell reads as an
    # override even when it holds the registry rule — never a silent omission.
    plain <- snippet_entry(snp_set_config_snippet(fx_record(), list()))
    expect_true("climate_multivariate::RDA" %in% names(plain$overrides))
})

test_that("a config record without a master rule takes its most common cell rule", {
    rec <- list(name = "gea.strict", origin = "config", source_module = "GEA",
                regime = "snp", strategy = "Cross-method", clumping_distance = 5e4,
                pairs = list("bio_1::EMMAX", "bio_2::EMMAX", "bio_3::EMMAX",
                             "climate_multivariate::RDA"),
                cells = list("bio_1::EMMAX" = list(type = "custom", value = 1e-6),
                             "bio_2::EMMAX" = list(type = "bonf", value = 0.2),
                             "bio_3::EMMAX" = list(type = "bonf", value = 0.2),
                             "climate_multivariate::RDA" = list(type = "bonf", value = 0.01)))
    e <- snippet_entry(snp_set_config_snippet(rec, fx_registry))
    expect_identical(e$name, "strict")                      # prefix stripped
    expect_identical(c(e$threshold_type, as.character(e$threshold_value)), c("bonf", "0.2"))
    expect_identical(names(e$overrides), "bio_1::EMMAX")
    expect_identical(e$clumping_distance, 50000L)
})

test_that("a name with dots becomes a valid config name", {
    e <- snippet_entry(snp_set_config_snippet(fx_record(name = "bonf.05.union")))
    expect_identical(e$name, "bonf_05_union")
})

test_that("a GWAS record yields a GWAS.snp_sets block", {
    rec <- fx_record(pairs = "height::EMMAX")
    rec$source_module <- "GWAS"
    expect_match(snp_set_config_snippet(rec)$yaml, "^GWAS:\n  snp_sets:")
})

test_that("a WZA set yields a note instead of a block", {
    sn <- snp_set_config_snippet(fx_record(regime = "wza"))
    expect_null(sn$yaml)
    expect_match(sn$note, "WZA")
    expect_null(fx_record(regime = "wza")$config_snippet)
})

test_that("a legacy manifest-only set gets a traits x methods block and a warning note", {
    project <- new_project("LEGACY")
    save_snp_set(project, "old.set", fx_sigsnps(), list(
        source_module = "GEA", threshold_type = "bonf", threshold_value = 0.05,
        threshold_overrides = list("bio_1::LFMM" = list(type = "custom", value = 1e-6),
                                   "bio_9::LFMM" = list(type = "custom", value = 1e-6)),
        regime = "snp", strategy = "All", clumping_distance = 100000L,
        traits = "bio_1", methods = c("EMMAX", "LFMM")))
    d <- snp_set_details(project, "old.set")
    expect_false(d$has_record)
    expect_match(d$snippet$note, "exact trait x method cells are\\s+unknown")
    e <- snippet_entry(d$snippet)
    expect_identical(e$name, "old_set")
    expect_identical(e$traits, "bio_1")
    expect_match(d$snippet$yaml, "traits:\n    - bio_1\n", fixed = TRUE)
    expect_identical(names(e$overrides), "bio_1::LFMM")    # bio_9 is outside the block
    expect_identical(e$strategy, "Union")
})

test_that("snp_set_details prefers the record and builds a snippet for a config set", {
    project <- new_project("DETAILS")
    fx_config_set(project, "gea.strict", 5L)
    d <- snp_set_details(project, "gea.strict", fx_registry)
    expect_true(d$has_record)
    expect_identical(d$origin, "config")
    expect_identical(snippet_entry(d$snippet)$name, "strict")

    save_snp_set(project, "mysel", fx_sigsnps(), list(), record = fx_record(name = "mysel"))
    d2 <- snp_set_details(project, "mysel")
    expect_identical(d2$origin, "shiny")
    expect_identical(d2$snippet$yaml, read_snp_set_record(project, "mysel")$config_snippet)

    expect_null(snp_set_details(project, "nosuch"))
})

test_that("the snippet block renders the YAML with a copy button, or only the note", {
    html <- as.character(snp_set_snippet_ui(list(yaml = "GEA:\n", note = NULL), "x1"))
    expect_match(html, '<pre id="x1"')
    expect_match(html, "isSecureContext")
    note_only <- as.character(snp_set_snippet_ui(list(yaml = NULL, note = "WZA"), "x2"))
    expect_false(grepl("<pre", note_only))
})
