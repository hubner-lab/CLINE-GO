# THE GATE. Runs unconditionally, and is the reason the seven
# skip_if_not(exists("compute_pval_threshold", ...)) guards elsewhere in this suite
# are safe to leave alone.
#
# The trap it closes: test-fct_threshold_rules.R:170,202,230,244,275,291,310 each
# skip when the shared pipeline libs are not in the package namespace. Those seven
# blocks are the app suite's ONLY real science assertions. If .onLoad() ever stops
# delivering the libs — a broken mount, a moved path, a silent sys.source() no-op —
# every one of them turns into a skip and the suite still reports green. A suite
# that cannot tell "passed" from "did not run" asserts nothing.
#
# So the CAUSE is asserted here instead, actively and without a skip, in two
# separately-diagnosing blocks: one says "the files are not there", the other says
# "the files are there but the namespace did not get them". This is the same
# instinct as the structural qvalue assertion at
# tests/testthat/test-equivalence-app-pipeline.R:570 — when the behaviour cannot be
# checked directly, check the thing whose breakage would hide it.
#
# Consequence, accepted deliberately: running this suite with no /pipeline mount
# now produces ONE clear failure instead of a green run with silent skips. That is
# the point.

test_that("CLINEGO_SHARED_LIBS is non-empty and every path exists", {
    expect_gt(length(CLINEGO_SHARED_LIBS), 0L)
    missing <- CLINEGO_SHARED_LIBS[!file.exists(CLINEGO_SHARED_LIBS)]
    expect_identical(
        missing, character(0),
        info = paste0(
            "shared pipeline libraries not found:\n  ",
            paste(missing, collapse = "\n  "),
            "\nMount the pipeline root at /pipeline (-v $PWD:/pipeline). ",
            "Without them the threshold/region tests in this suite SKIP rather ",
            "than fail, so the suite would go green having asserted nothing."))
})

test_that(".onLoad actually put the shared functions in the package namespace", {
    # file.exists() passing is not enough: sys.source() into asNamespace() is the
    # step that can silently do nothing (zzz.R:59-61), and inherits = FALSE is
    # essential — exists() would otherwise find these via the search path in a
    # session that happened to attach them.
    ns <- asNamespace("clinego.app")
    for (fn in c("compute_pval_threshold", "cluster_snps_to_regions")) {
        expect_true(exists(fn, envir = ns, inherits = FALSE),
                    info = paste0(fn, " is not in the clinego.app namespace — ",
                                  "zzz.R's load_shared_libs() did not deliver it"))
    }
})

test_that("load_shared_libs reports the paths it could not find", {
    e <- new.env(parent = baseenv())
    missing <- load_shared_libs(e, paths = c("/nonexistent/a.R", "/nonexistent/b.R"))
    expect_identical(missing, c("/nonexistent/a.R", "/nonexistent/b.R"))
    # Nothing was defined from files that do not exist.
    expect_length(ls(e), 0L)
})

test_that("load_shared_libs sources what it finds and returns only the misses", {
    d <- withr::local_tempdir()
    good <- file.path(d, "good.R")
    writeLines("helper_from_good <- function() 42", good)

    # baseenv(), not emptyenv(): sourced code needs `<-` to be resolvable.
    e <- new.env(parent = baseenv())
    missing <- load_shared_libs(e, paths = c(good, "/nonexistent/b.R"))

    expect_identical(missing, "/nonexistent/b.R")
    expect_true(exists("helper_from_good", envir = e, inherits = FALSE))
    expect_identical(get("helper_from_good", envir = e)(), 42)
})

test_that("load_shared_libs returns an empty character vector when all paths load", {
    d <- withr::local_tempdir()
    f <- file.path(d, "a.R")
    writeLines("x_from_a <- 1", f)
    expect_identical(load_shared_libs(new.env(parent = baseenv()), paths = f),
                     character(0))
})
