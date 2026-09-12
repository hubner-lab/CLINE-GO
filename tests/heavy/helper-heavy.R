# Helper for the HEAVY test root. Sourced automatically by test_dir() because of
# the helper- prefix.
#
# Reuses the quick root's helper verbatim (package attaches + the shared lib
# source order + quiet()) and then the shared wrapper harness. Nothing is
# duplicated between the two roots: that is exactly why tests/lib/wrapper_harness.R
# exists as a third file rather than living in tests/testthat/.

.heavy_root <- getOption(
    "clinego.repo_root",
    normalizePath(file.path(testthat::test_path(), ".."), mustWork = TRUE)
)

source(file.path(.heavy_root, "tests", "testthat", "helper-libs.R"))
source(file.path(.heavy_root, "tests", "lib", "wrapper_harness.R"))
source(file.path(.heavy_root, "tests", "heavy", "heavy_wrappers.R"))
