# R/app_config.R — get_pipeline_path()'s four-tier priority chain.
#
# This function decides where EVERY {PROJECT}_results/ read and write lands, so
# the tier ORDER is the contract, not the value any one tier happens to return.
# The tests assert ordering and the guard on each tier; they deliberately do not
# assert "/pipeline", which resolves differently inside and outside the container.
#
# The asymmetry worth knowing: tier 1 is a bare is.null() check with no nzchar()
# guard, unlike tiers 2 and 3. options(clinego.pipeline_path = "") therefore wins
# and returns "", where the same empty string in the golem config or the env var
# would fall through. That is also the tier helper-pipeline-path.R uses to
# redirect the whole suite, so the suite's hermeticity rests on it.
#
# Tier 2 is NOT dormant in the shipped package: inst/golem-config.yml sets
# pipeline_path in both `default:` and `production:`, so it always returns a
# non-empty value and tiers 3-4 are unreachable unless that file is missing or
# config::get errors. Every test below that needs a lower tier therefore has to
# neutralise tier 2 explicitly rather than assume it is quiet.

# Tier 2 reads through get_golem_config(); forcing it to error is the only way to
# reach tiers 3 and 4 with the config file present, and it exercises the
# tryCatch at :37 at the same time.
with_golem_config_failing <- function(code) {
    testthat::local_mocked_bindings(
        get_golem_config = function(...) stop("no config"),
        .package = "clinego.app"
    )
    force(code)
}

test_that("tier 1: a runtime option beats everything downstream", {
    withr::local_options(clinego.pipeline_path = "/tier1")
    withr::local_envvar(PIPELINE_PATH = "/tier3")
    expect_identical(get_pipeline_path(), "/tier1")
})

test_that("tier 1 has no nzchar guard, so an empty option wins and returns \"\"", {
    # Characterisation of a real asymmetry: tiers 2 and 3 both require nzchar,
    # tier 1 only requires non-NULL. An empty string here silently redirects
    # every project path to the working directory.
    withr::local_options(clinego.pipeline_path = "")
    withr::local_envvar(PIPELINE_PATH = "/tier3")
    expect_identical(get_pipeline_path(), "")
})

test_that("tier 2: with the option unset, the golem config supplies the path", {
    withr::local_options(clinego.pipeline_path = NULL)
    withr::local_envvar(PIPELINE_PATH = "/tier3")
    # The env var must NOT win here — the shipped config sets pipeline_path, so
    # tier 2 answers first.
    expect_false(identical(get_pipeline_path(), "/tier3"))
    expect_true(nzchar(get_pipeline_path()))
})

test_that("tier 3: the env var is used only once the option and config are gone", {
    withr::local_options(clinego.pipeline_path = NULL)
    withr::local_envvar(PIPELINE_PATH = "/tier3")
    with_golem_config_failing({
        expect_identical(get_pipeline_path(), "/tier3")
    })
})

test_that("tier 3 is skipped when the env var is empty", {
    withr::local_options(clinego.pipeline_path = NULL)
    withr::local_envvar(PIPELINE_PATH = "")
    with_golem_config_failing({
        expect_identical(get_pipeline_path(), "/pipeline")
    })
})

test_that("tier 4: the literal default is the last resort", {
    withr::local_options(clinego.pipeline_path = NULL)
    withr::local_envvar(PIPELINE_PATH = NULL)
    with_golem_config_failing({
        expect_identical(get_pipeline_path(), "/pipeline")
    })
})

test_that("a config error is caught rather than propagated", {
    withr::local_options(clinego.pipeline_path = NULL)
    withr::local_envvar(PIPELINE_PATH = "/tier3")
    with_golem_config_failing({
        # The tryCatch at :37 is what keeps a malformed golem-config.yml from
        # taking down every path lookup in the app.
        expect_no_error(get_pipeline_path())
    })
})

test_that("the test helper's redirect is in force for the rest of the suite", {
    # Not a tautology: this is the assertion that the suite is hermetic. If the
    # helper's plain options() call were ever changed to withr::local_options(),
    # it would unwind when the helper finished sourcing and every project write
    # would land in the real /pipeline tree instead.
    expect_identical(getOption("clinego.pipeline_path"), CLINEGO_TEST_ROOT)
    expect_identical(get_pipeline_path(), CLINEGO_TEST_ROOT)
    expect_true(dir.exists(CLINEGO_TEST_ROOT))
})

# --------------------------------------------------- www/ asset resource prefix

# app_ui() references two scripts as src = "www/...", which resolve only through
# a registered "www" resource prefix. Until 2026-09-15 the only addResourcePath("www", ...)
# in the whole package lived in dev.R, so clinego.app::run_app() — the documented
# production path — served 404s for both, losing the config-dirty handler and the
# module bar that drives tab switching. These tests read the sources, because the
# defect is a relationship between two files rather than a value.

.app_r <- function(f) file.path(testthat::test_path("..", ".."), "R", f)

test_that("every www/ asset referenced by the package exists in inst/app/www", {
    code <- unlist(lapply(c("app_ui.R", "utils_ui.R"), function(f) {
        src <- readLines(.app_r(f), warn = FALSE)
        src[!grepl("^\\s*#", src)]            # comments mention the pattern too
    }))
    refs <- regmatches(code, gregexpr('src\\s*=\\s*"www/[^"]+"', code))
    refs <- sub('^.*"www/([^"]+)"$', "\\1", unlist(refs))

    expect_gt(length(refs), 0)
    for (r in refs) expect_true(file.exists(app_sys("app", "www", r)), info = r)
})

test_that("app_ui registers the www resource prefix itself, not only dev.R", {
    ui <- paste(readLines(.app_r("app_ui.R"), warn = FALSE), collapse = "\n")
    expect_match(ui, 'addResourcePath\\(\\s*"www"')
})

test_that("the registration is skipped when a www prefix is already set", {
    # dev.R registers the BIND MOUNT before app_ui() ever runs; re-registering
    # the installed copy there would shadow the file being edited.
    ui <- paste(readLines(.app_r("app_ui.R"), warn = FALSE), collapse = "\n")
    expect_match(ui, '"www" %in% names\\(shiny::resourcePaths\\(\\)\\)')
})
