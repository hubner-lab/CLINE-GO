# R/mod_home.R — the Home tab's KPI row (bffd38).
#
# The boxes asked pipeline_summary.tsv for metrics it never had (n_samples,
# n_snps_filtered, n_regions_combined), and the lookup ran inside a data.table's
# [ ], where its `step` and `metric` arguments name the table's own columns: every
# box showed the first row's value.

home_summary <- function(project, rows) {
    dir.create(project_base(project), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(rows, pipeline_summary_path(project), sep = "\t")
}

test_that("the Home KPIs read the summary's own metric names", {
    a <- paste0("home_", as.integer(Sys.time()))
    withr::defer(unlink(project_base(a), recursive = TRUE))
    home_summary(a, data.table::data.table(
        step   = c("traits", "processing", "processing", "gea", "gwas"),
        metric = c("n_traits", "samples_after_filtering", "snps_after_filtering",
                   "regions_combined", "regions_combined"),
        value  = c("3", "50", "337", "4", "2")))

    pd <- shiny::reactiveVal(list(name = a, config = list(), k_best = 3L))
    shiny::testServer(mod_home_server, args = list(project_data = pd,
                                                   has_project = shiny::reactive(TRUE)), {
        expect_identical(output$n_samples, "50")
        expect_identical(output$n_snps, "337")
        expect_identical(output$n_regions, "4")
        # a gwas_only project reports its GWAS regions
        pd(list(name = a, config = list(Regime = list(mode = "gwas_only")), k_best = 3L))
        session$flushReact()
        expect_identical(output$n_regions, "2")
    })
})

# Home's Pipeline Status (N3, 813d39): every module of the bar, Environmental included,
# and GEAxGWAS done when its own outputs exist (mode=gea_x_gwas writes no summary row).
test_that("the status list has Environmental and reads GEAxGWAS from its outputs", {
    a <- paste0("home_status_", as.integer(Sys.time()))
    withr::defer(unlink(project_base(a), recursive = TRUE))
    home_summary(a, data.table::data.table(step = c("processing", "climate"),
                                           metric = c("snps_raw", "n_predictors"),
                                           value = c("354", "3")))
    st <- check_module_status(a)
    expect_true(st[["climate"]])
    expect_false(st[["gea_x_gwas"]])
    plots <- mod_path(a, MOD_GEAXGWAS, "plots")
    dir.create(plots, recursive = TRUE)
    file.create(file.path(plots, "miami_combined_K3.png"))
    expect_true(check_module_status(a)[["gea_x_gwas"]])
    expect_identical(geaxgwas_outputs(a), "miami_combined_K3.png")

    pd <- shiny::reactiveVal(list(name = a, config = list(), k_best = 3L))
    shiny::testServer(mod_home_server, args = list(project_data = pd,
                                                   has_project = shiny::reactive(TRUE)), {
        html <- as.character(output$module_status$html)
        expect_match(html, "Environmental")
        expect_match(html, "GEAxGWAS")
    })
})

test_that("Home without a project says so instead of staying blank (I11)", {
    shiny::testServer(mod_home_server,
                      args = list(project_data = shiny::reactive(shiny::req(FALSE)),
                                  has_project = shiny::reactive(FALSE)), {
        expect_match(as.character(output$home_content$html), "No project yet")
    })
})

# R/mod_gwas.R — the GWAS missing-data alert's colour (154411): filling a missing
# outcome with its MEAN / MEDIAN is the case to warn about, not dropping samples.

test_that("the missing-data alert warns on a fill, not on DROP", {
    tbl <- function(strategy, n_available, n_total = 50L)
        data.table::data.table(trait = "height", n_total = n_total,
                               n_available = n_available, strategy = strategy)
    expect_false(pheno_missing_alert_warns(tbl("DROP", 44L)))
    expect_true(pheno_missing_alert_warns(tbl("MEAN", 44L)))
    expect_true(pheno_missing_alert_warns(tbl("median", 44L)))
    expect_false(pheno_missing_alert_warns(tbl("MEAN", 50L)))     # nothing was filled
    expect_true(pheno_missing_alert_warns(tbl("DROP", 20L)))      # under half its samples
})
