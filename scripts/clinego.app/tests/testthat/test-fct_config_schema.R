# R/fct_config_schema.R — the dot-path config accessors and the input coercion
# that every config write goes through.
#
# Why this is worth pinning: config_set_by_path() is what the sidebar uses to edit
# a project's YAML, and input_to_config_value() decides the TYPE that lands in the
# file. Its numeric branch carries a comment explaining that a double where an
# integer belongs forks a second output directory (common.smk interpolates several
# of these into paths), so a regression here silently duplicates results trees
# rather than raising anything.
#
# Only the `numeric` branch had coverage before (test-fct_discovery.R:51); the five
# other branches are asserted here.

test_that("config_get_by_path walks a nested config", {
    cfg <- list(sNMF = list(k_start = 2L, k_end = 7L),
                Input = list(vcf = "x.vcf"))
    expect_identical(config_get_by_path(cfg, "sNMF.k_start"), 2L)
    expect_identical(config_get_by_path(cfg, "Input.vcf"), "x.vcf")
    # A single-segment path is the whole sub-list.
    expect_identical(config_get_by_path(cfg, "sNMF"), cfg$sNMF)
})

test_that("config_get_by_path returns NULL rather than erroring on a bad path", {
    cfg <- list(sNMF = list(k_start = 2L), Filter = list(maf = 0.05))
    expect_null(config_get_by_path(cfg, "sNMF.nope"))
    expect_null(config_get_by_path(cfg, "nope.k_start"))
    # Walking THROUGH a leaf: maf is a number, so maf.deeper cannot exist.
    expect_null(config_get_by_path(cfg, "Filter.maf.deeper"))
    expect_null(config_get_by_path(list(), "anything"))
})

test_that("config_set_by_path sets a leaf without disturbing its siblings", {
    cfg <- list(sNMF = list(k_start = 2L, k_end = 7L))
    got <- config_set_by_path(cfg, "sNMF.k_end", 9L)
    expect_identical(got$sNMF$k_end, 9L)
    expect_identical(got$sNMF$k_start, 2L)
})

test_that("config_set_by_path creates missing intermediate lists", {
    got <- config_set_by_path(list(), "Maladaptation.methods.gradient_forest.ntree", 500L)
    expect_identical(got$Maladaptation$methods$gradient_forest$ntree, 500L)
})

test_that("config_set_by_path is a pure function of its input", {
    cfg <- list(sNMF = list(k_start = 2L))
    invisible(config_set_by_path(cfg, "sNMF.k_start", 99L))
    # R copy-on-modify, but the sidebar relies on it: assert rather than assume.
    expect_identical(cfg$sNMF$k_start, 2L)
})

test_that("config_set_by_path REPLACES a non-list intermediate", {
    # Current behaviour, pinned deliberately (:553): descending through a leaf
    # clobbers it rather than erroring. The sidebar never does this, but a schema
    # edit that moves a key from leaf to group would hit it, and silently losing
    # the old value is the kind of thing worth having written down.
    cfg <- list(Climate = "enabled")
    got <- config_set_by_path(cfg, "Climate.predictors", list("bio_1"))
    expect_true(is.list(got$Climate))
    expect_identical(got$Climate$predictors, list("bio_1"))
    expect_false("enabled" %in% unlist(got$Climate))
})

test_that("config_set_by_path accepts a single-segment path", {
    got <- config_set_by_path(list(a = 1), "project_name", "SIMDATA")
    expect_identical(got$project_name, "SIMDATA")
})

test_that("input_to_config_value writes whole numbers as integers", {
    # The documented reason: yaml::write_yaml() renders an R double as "7.0", and
    # common.smk interpolates k_end / r2 / window straight into output PATHS, so a
    # double forks a second directory and orphans the first.
    got <- input_to_config_value(7, "numeric")
    expect_type(got, "integer")
    expect_identical(got, 7L)

    frac <- input_to_config_value(0.05, "numeric")
    expect_type(frac, "double")
    expect_identical(frac, 0.05)
})

test_that("input_to_config_value rejects a non-numeric numeric input", {
    expect_null(input_to_config_value("", "numeric"))
    expect_null(input_to_config_value("abc", "numeric"))
    expect_null(input_to_config_value(NA, "numeric"))
})

test_that("input_to_config_value coerces checkbox and checkbox_invert", {
    expect_true(input_to_config_value(TRUE, "checkbox"))
    expect_false(input_to_config_value(FALSE, "checkbox"))
    expect_false(input_to_config_value(NULL, "checkbox"))

    # checkbox_invert backs a UI switch whose sense is the opposite of the config
    # key (the config says "disable X", the switch says "enable X").
    expect_false(input_to_config_value(TRUE, "checkbox_invert"))
    expect_true(input_to_config_value(FALSE, "checkbox_invert"))
    expect_true(input_to_config_value(NULL, "checkbox_invert"))
})

test_that("input_to_config_value splits a textarea into a YAML list", {
    got <- input_to_config_value("bio_1, bio_12 ,bio_5", "textarea")
    expect_type(got, "list")
    expect_identical(got, list("bio_1", "bio_12", "bio_5"))
})

test_that("input_to_config_value coerces an all-numeric textarea to numbers", {
    got <- input_to_config_value("2, 3, 4", "textarea")
    expect_identical(got, list(2, 3, 4))
    # One non-numeric entry keeps the whole list character — mixing types in a
    # YAML sequence would be worse than keeping strings.
    expect_identical(input_to_config_value("2, x", "textarea"), list("2", "x"))
})

test_that("input_to_config_value drops empty textarea entries and blank input", {
    expect_identical(input_to_config_value("bio_1,,bio_2", "textarea"),
                     list("bio_1", "bio_2"))
    expect_null(input_to_config_value("", "textarea"))
    expect_null(input_to_config_value("   ", "textarea"))
    expect_null(input_to_config_value(",,", "textarea"))
    expect_null(input_to_config_value(NULL, "textarea"))
})

test_that("input_to_config_value passes method_table through untouched", {
    tbl <- list(list(method = "EMMAX", adjust = "bonf", threshold = "0.05"))
    expect_identical(input_to_config_value(tbl, "method_table"), tbl)
})

test_that("input_to_config_value treats text, select and unknown types alike", {
    for (ty in c("text", "select", "something_new")) {
        expect_identical(input_to_config_value("  SIMDATA  ", ty), "  SIMDATA  ")
        expect_null(input_to_config_value("   ", ty))
        expect_null(input_to_config_value(NULL, ty))
    }
})

test_that("schema_for_tab and schema_split partition the schema", {
    entries <- schema_for_tab("processing")
    expect_true(length(entries) > 0)
    expect_true(all(vapply(entries, function(e) e$tab == "processing", logical(1))))

    split <- schema_split(entries)
    expect_named(split, c("mandatory", "optional"))
    # Every entry lands in exactly one half.
    expect_identical(length(split$mandatory) + length(split$optional), length(entries))
})
