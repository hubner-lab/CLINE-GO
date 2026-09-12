# Quick-tier CANARY for the heavy tier. Runs in the DEFAULT gate.
#
# The heavy tests live in their own root (tests/heavy/) and only run with
# `tests/run_all.sh --heavy`, which means nobody would notice if they rotted: a
# renamed script, a deleted file or an emptied test file would leave the default
# gate green and the heavy gate quietly covering less than it claims.
#
# So the default gate asserts the heavy tier's SHAPE without running it — cheap,
# hermetic, no tool touched. This is the same instinct as
# test-equivalence-app-pipeline.R:570's structural qvalue assertion: when the
# behaviour cannot be checked here, check the cause that would break it.

heavy_dir <- file.path(getOption("clinego.repo_root", "/pipeline"), "tests", "heavy")

test_that("the heavy tier still exists and is wired to a runner", {
    skip_if_not(dir.exists(heavy_dir), "tests/heavy not mounted")

    expect_true(file.exists(file.path(heavy_dir, "heavy_wrappers.R")))
    expect_true(file.exists(file.path(heavy_dir, "helper-heavy.R")))
    expect_true(file.exists(file.path(
        getOption("clinego.repo_root", "/pipeline"), "tests", "run_heavy.R")))

    test_files <- list.files(heavy_dir, pattern = "^test-.*\\.R$")
    expect_gt(length(test_files), 0L)
})

test_that("every heavy test file contains at least one test_that block", {
    # Catches the failure mode a separate root makes invisible: a heavy file
    # emptied or commented out while the gate keeps reporting PASS.
    skip_if_not(dir.exists(heavy_dir), "tests/heavy not mounted")

    for (f in list.files(heavy_dir, pattern = "^test-.*\\.R$", full.names = TRUE)) {
        src <- readLines(f, warn = FALSE)
        expect_true(any(grepl("test_that\\(", src)),
                    info = paste0(basename(f), " contains no test_that block"))
    }
})

test_that("the heavy spec table is well formed and names scripts that exist", {
    skip_if_not(dir.exists(heavy_dir), "tests/heavy not mounted")
    # heavy_wrappers.R is side-effect free by design, so the quick tier can read
    # it. Source it into a private env: it defines REPO-derived helpers that must
    # not leak into the shared test_dir() environment.
    e <- new.env(parent = environment())
    e$REPO <- getOption("clinego.repo_root", "/pipeline")
    e$SCRIPTS <- file.path(e$REPO, "scripts")
    sys.source(file.path(heavy_dir, "heavy_wrappers.R"), envir = e)

    expect_true(length(e$HEAVY_WRAPPERS) > 0L)
    for (s in e$HEAVY_WRAPPERS) {
        expect_true(all(c("label", "script", "build", "args", "outputs") %in% names(s)),
                    info = paste0("incomplete heavy spec: ", s$label))
        expect_true(file.exists(file.path(e$SCRIPTS, s$script)),
                    info = paste0("heavy spec names a missing script: ", s$script))
    }
})

test_that("no heavy spec targets a wrapper on the denylist", {
    skip_if_not(dir.exists(heavy_dir), "tests/heavy not mounted")
    e <- new.env(parent = environment())
    e$REPO <- getOption("clinego.repo_root", "/pipeline")
    e$SCRIPTS <- file.path(e$REPO, "scripts")
    sys.source(file.path(heavy_dir, "heavy_wrappers.R"), envir = e)

    used <- vapply(e$HEAVY_WRAPPERS, function(s) s$script, character(1))
    # WRAPPER_DENYLIST comes from the shared harness, sourced by
    # test-cli-wrappers.R into the same parent environment.
    expect_length(intersect(used, WRAPPER_DENYLIST), 0L)
})
