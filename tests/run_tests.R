#!/usr/bin/env Rscript
# Unit tests for scripts/R/lib and scripts/R/utils.
#
# This is a SECOND testthat root, separate from scripts/clinego.app/tests. That one
# is package-scoped (test_check("clinego.app")) and can only test code that lives
# inside the package; the pipeline libraries are plain .R files with no DESCRIPTION
# and no NAMESPACE, so they need test_dir() and a helper that source()s them.
#
# Run (the -v mount is required — the image ships no pipeline code):
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -v $PWD:/pipeline \
#     cline-go:latest Rscript /pipeline/tests/run_tests.R
#
# Exits non-zero on any failure or error, so it can gate a merge.
#
# Do NOT route this through scripts/clinego.app/dev.R and do not prepend
# .R_libs_dev to .libPaths(): that directory carries testthat 3.3.2 while the
# image pins 3.2.3, and the two disagree about test_dir()'s return shape.

suppressPackageStartupMessages(library(testthat))

.here <- local({
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) == 0) getwd()
    else dirname(normalizePath(sub("^--file=", "", file_arg)))
})

# Published for helper-libs.R, which must locate scripts/R/ without depending on
# test_path() semantics (they differ between a package run and a bare test_dir()).
options(clinego.repo_root = normalizePath(file.path(.here, "..")))

results <- test_dir(file.path(.here, "testthat"),
                    reporter = "summary", stop_on_failure = FALSE)

df <- as.data.frame(results)
if (any(df$failed > 0 | df$error > 0)) quit(status = 1)
