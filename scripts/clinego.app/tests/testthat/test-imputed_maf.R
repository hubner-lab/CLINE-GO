# sNMF.imputed_maf in the app: the sidebar field, the summary reader and the
# before -> after KPI box on the GEA / GWAS tabs.
#
# The pipeline side (common.smk) reads a blank / null value as "use Filter.maf" and
# 0 as "filter off". The sidebar must therefore write NOTHING for a blank field and
# a real 0 for zero — a schema or coercion change that wrote 0 for blank would turn
# the filter off without the user asking.

entry_for <- function(key) {
    hits <- Filter(function(e) identical(e$key, key), config_schema())
    expect_length(hits, 1L)
    hits[[1]]
}

test_that("sNMF.imputed_maf is an optional numeric on the Structure tab", {
    e <- entry_for("sNMF.imputed_maf")
    expect_identical(e$tab, "structure")
    expect_identical(e$section, "Structure")
    expect_identical(e$type, "numeric")
    expect_false(e$mandatory)
    expect_null(e$show_if)
    expect_equal(c(e$min, e$max, e$step), c(0, 0.49, 0.01))
    expect_match(e$help, "Filter.maf")
    expect_match(e$help, "0 disables")
    expect_true(any(vapply(schema_split(schema_for_tab("structure"))$optional,
                           function(x) identical(x$key, "sNMF.imputed_maf"), logical(1))))
})

test_that("the field stays visible in a gwas_only project", {
    # GWAS tests the same imputed matrix, so the filter applies in both regimes.
    expect_false("Structure" %in% hidden_sidebar_sections("gwas_only"))
})

test_that("a blank field writes nothing and zero writes zero", {
    e <- entry_for("sNMF.imputed_maf")
    expect_null(input_to_config_value(NA, e$type))
    expect_null(input_to_config_value("", e$type))
    expect_identical(input_to_config_value(0, e$type), 0L)
    expect_identical(input_to_config_value(0.1, e$type), 0.1)
})

test_that("the default config ships the key as null, i.e. Filter.maf", {
    path <- system.file("config_default.yaml", package = "clinego.app")
    skip_if(!nzchar(path), "config_default.yaml not installed")
    cfg <- yaml::read_yaml(path)
    expect_true("imputed_maf" %in% names(cfg$sNMF))
    expect_null(cfg$sNMF$imputed_maf)
})

summary_with_imputed_maf <- function() data.frame(
    step   = c("processing", "structure", "structure", "structure", "gea", "gea", "gea"),
    metric = c("snps_after_filtering", "imputed_maf_threshold", "snps_before_imputed_maf",
               "snps_after_imputed_maf", "imputed_maf_threshold", "snps_before_imputed_maf",
               "snps_after_imputed_maf"),
    value  = c("209466", "0.1", "6380", "6033", "0.1", "209466", "174568"),
    stringsAsFactors = FALSE
)

test_that("imputed_maf_counts reads one step's three rows as numbers", {
    s <- summary_with_imputed_maf()
    expect_identical(imputed_maf_counts(s, "gea"),
                     list(threshold = 0.1, before = 209466, after = 174568))
    expect_identical(imputed_maf_counts(s, "structure"),
                     list(threshold = 0.1, before = 6380, after = 6033))
})

test_that("imputed_maf_counts returns NULL for a run that predates the filter", {
    s <- summary_with_imputed_maf()
    expect_null(imputed_maf_counts(s, "gwas"))
    expect_null(imputed_maf_counts(s[s$metric != "snps_after_imputed_maf", ], "gea"))
    expect_null(imputed_maf_counts(NULL, "gea"))
    expect_null(imputed_maf_counts(s[0, ], "gea"))
})

test_that("imputed_maf_counts keeps a zero threshold (filter off) rather than dropping it", {
    s <- summary_with_imputed_maf()
    s$value[s$step == "gea" & s$metric == "imputed_maf_threshold"] <- "0"
    s$value[s$step == "gea" & s$metric == "snps_after_imputed_maf"] <- "209466"
    expect_identical(imputed_maf_counts(s, "gea")$threshold, 0)
})

box_html <- function(x) as.character(htmltools::renderTags(x)$html)

test_that("imputed_maf_value_box shows the threshold and before -> after counts", {
    html <- box_html(imputed_maf_value_box(list(threshold = 0.1, before = 209466,
                                                after = 174568)))
    expect_match(html, "Imputed MAF filter (≥ 0.1)", fixed = TRUE)
    expect_match(html, "209,466", fixed = TRUE)
    expect_match(html, "174,568", fixed = TRUE)
    expect_match(html, "→", fixed = TRUE)
})

test_that("imputed_maf_value_box labels a zero threshold as off and a missing run as a dash", {
    off <- box_html(imputed_maf_value_box(list(threshold = 0, before = 337, after = 337)))
    expect_match(off, "Imputed MAF filter (off)", fixed = TRUE)
    none <- box_html(imputed_maf_value_box(NULL))
    expect_match(none, "—", fixed = TRUE)
    expect_false(grepl("→", none, fixed = TRUE))
})
