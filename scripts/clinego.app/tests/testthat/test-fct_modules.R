# R/fct_modules.R — the module bar's registry and the regime filter over it.
#
# Two things here fail silently rather than loudly.
#
# 1. The regime filter is `identical(regime, "gwas_only")` and nothing else. Any
#    other string — a typo, a different case, a NULL coerced upstream — returns
#    the FULL registry, so a gwas_only project would render Environmental, PreGEA,
#    GEA, GEAxGWAS and Maladaptation tabs for outputs it can never produce.
#
# 2. module_groups() relies on `levels = unique(reg$group)` to keep bar order. A
#    bare split() sorts alphabetically, which would silently reorder the module
#    bar into Association / Factors / Maladaptation / Setup / Structure — still
#    "working", just not the pipeline's order.

REG_IDS <- c("home", "processing", "prestructure", "structure",
             "climate", "traits", "pregea", "gea", "gwas",
             "gea_x_gwas", "maladaptation")

# The six that survive gwas_only. `traits` is in the list deliberately: it is the
# only factor-characterization module such a project has, since mode=climate
# raises without climate while mode=traits needs neither coordinates nor climate.
GWAS_ONLY_IDS <- c("home", "processing", "prestructure", "structure",
                   "traits", "gwas")

# ---------------------------------------------------------------- registry

test_that("module_registry lists the eleven modules in pipeline order", {
    reg <- module_registry()
    expect_identical(reg$id, REG_IDS)
    expect_identical(nrow(reg), 11L)
})

test_that("module_registry's parallel columns are all the same length", {
    reg <- module_registry()
    expect_identical(length(reg$label), 11L)
    expect_identical(length(reg$group), 11L)
    expect_identical(length(reg$gwas_only), 11L)
    # A short recycle here would mislabel modules without erroring.
    expect_false(any(is.na(reg$label)))
    expect_false(any(is.na(reg$group)))
})

test_that("module_registry keeps ids and labels as character, not factors", {
    reg <- module_registry()
    expect_type(reg$id, "character")
    expect_type(reg$label, "character")
    expect_type(reg$group, "character")
    expect_type(reg$gwas_only, "logical")
})

# ---------------------------------------------------------------- module_input_id

test_that("module_input_id prefixes the registry id", {
    expect_identical(module_input_id("gea"), "mb_gea")
    # Vectorized, because the bar builder maps it over the registry.
    expect_identical(module_input_id(c("home", "gwas")), c("mb_home", "mb_gwas"))
})

test_that("module_input_id covers every registry id without collision", {
    ids <- module_input_id(module_registry()$id)
    expect_identical(length(unique(ids)), 11L)
})

# ---------------------------------------------------------------- regime filter

test_that("module_ids_for_regime returns all eleven under standard", {
    expect_identical(module_ids_for_regime("standard"), REG_IDS)
    expect_identical(module_ids_for_regime(), REG_IDS)
})

test_that("module_ids_for_regime drops the six climate-dependent modules under gwas_only", {
    expect_identical(module_ids_for_regime("gwas_only"), GWAS_ONLY_IDS)
    expect_false("climate" %in% module_ids_for_regime("gwas_only"))
    expect_false("maladaptation" %in% module_ids_for_regime("gwas_only"))
})

test_that("module_ids_for_regime keeps traits under gwas_only", {
    # Deliberate, and the registry comment says so: the only factor
    # characterization a gwas_only project can run.
    expect_true("traits" %in% module_ids_for_regime("gwas_only"))
})

test_that("module_ids_for_regime silently falls back to ALL modules on an unrecognised regime", {
    # Characterisation of a real trap: the filter is identical(), so anything
    # that is not exactly "gwas_only" yields the full registry rather than an
    # error. A typo in a config therefore shows tabs the project cannot fill.
    expect_identical(module_ids_for_regime("GWAS_ONLY"), REG_IDS)
    expect_identical(module_ids_for_regime("gwas-only"), REG_IDS)
    expect_identical(module_ids_for_regime(""), REG_IDS)
    expect_identical(module_ids_for_regime(NULL), REG_IDS)
})

# ---------------------------------------------------------------- hidden sections

test_that("hidden_sidebar_sections hides the three map-dependent sections under gwas_only", {
    expect_identical(hidden_sidebar_sections("gwas_only"),
                     c("Map", "Piemap", "Pop Stats"))
})

test_that("hidden_sidebar_sections hides nothing under standard", {
    expect_identical(hidden_sidebar_sections("standard"), character(0))
    expect_identical(hidden_sidebar_sections(), character(0))
})

test_that("hidden_sidebar_sections leaves LD Decay visible under gwas_only", {
    # Deliberate: with one site it degrades to the genome-wide "All" curve,
    # which is still worth having.
    expect_false("LD Decay" %in% hidden_sidebar_sections("gwas_only"))
})

# ---------------------------------------------------------------- module_groups

test_that("module_groups preserves pipeline bar order, not alphabetical order", {
    got <- names(module_groups("standard"))
    expect_identical(got, c("Setup", "Structure", "Factors",
                            "Association", "Maladaptation"))
    # The bug this pins: a bare split() would sort these.
    expect_false(identical(got, sort(got)))
})

test_that("module_groups returns four groups under gwas_only, not five", {
    got <- module_groups("gwas_only")
    # Maladaptation loses its only member and disappears entirely; Factors
    # survives with just traits.
    expect_identical(names(got), c("Setup", "Structure", "Factors", "Association"))
    expect_identical(nrow(got$Factors), 1L)
    expect_identical(got$Factors$id, "traits")
    expect_identical(got$Association$id, "gwas")
})

test_that("module_groups partitions the registry without loss or duplication", {
    for (regime in c("standard", "gwas_only")) {
        got <- module_groups(regime)
        ids <- unlist(lapply(got, function(g) g$id), use.names = FALSE)
        expect_setequal(ids, module_ids_for_regime(regime))
        expect_identical(length(ids), length(module_ids_for_regime(regime)))
    }
})
