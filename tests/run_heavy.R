#!/usr/bin/env Rscript
# HEAVY tier runner — the wrapper tests that need a genomics toolchain.
#
# A THIRD testthat root, deliberately separate from tests/testthat:
#   * the quick gate must stay fast, and test_dir() sources every file's top level
#     even for tests it will skip;
#   * a heavy test parked behind a skip in the shared root would inflate the
#     recorded "N skipped" count, whose entire meaning in this repo is
#     "correct-behaviour assertions for known, quarantined defects"
#     (tests/testthat/test-known-bugs.R). Two different meanings behind one number
#     is how a real regression gets read as a quarantine.
#
# CONTRACT: expected GREEN. Separated by RUNTIME and tool dependency, not by
# correctness. Contrast --invariants, which is expected RED by design (see
# tests/run_all.sh's header).
#
# Run (the -v mount is required — the image ships no pipeline code):
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -v $PWD:/pipeline \
#     cline-go:latest Rscript /pipeline/tests/run_heavy.R
# or, with the rest of the gate:
#   tests/run_all.sh --heavy
#
# Exits non-zero on any failure or error.

suppressPackageStartupMessages(library(testthat))

.here <- local({
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) == 0) getwd()
    else dirname(normalizePath(sub("^--file=", "", file_arg)))
})

# MUST be set before test_dir(): helper-heavy.R and helper-libs.R both resolve
# scripts/R via this option, and their test_path() fallback would otherwise
# resolve relative to tests/heavy/ instead of the repo root.
options(clinego.repo_root = normalizePath(file.path(.here, "..")))

results <- test_dir(file.path(.here, "heavy"),
                    reporter = "summary", stop_on_failure = FALSE)

df <- as.data.frame(results)
if (any(df$failed > 0 | df$error > 0)) quit(status = 1)
