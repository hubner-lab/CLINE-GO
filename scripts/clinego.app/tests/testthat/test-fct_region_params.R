# R/fct_region_params.R — the exploratory-parameter store, kept separate from the
# pipeline config YAML. It backs the UI rule in CLAUDE.md: never silently discard
# a user's chosen parameter value, and show "modified" against the SAVED params
# rather than the config default.
#
# Two behaviours here are easy to assume wrongly and both are pinned below.
#
# 1. save_region_params() returns invisible(params_list) whether or not the write
#    succeeded — the failure path only message()s. A caller cannot use the return
#    value to decide anything.
# 2. The JSON round trip is NOT identity for vectors. write_json(auto_unbox =
#    TRUE) unboxes length-1 vectors to scalars and read_json(simplifyVector =
#    FALSE) brings everything back as lists, so a stored c(1, 2) returns
#    list(1, 2). Tests assert the shape that actually comes back.
#
# Writes go under project_base(), which resolves through get_pipeline_path() —
# helper-pipeline-path.R has already redirected that to a session tempdir. Each
# block uses a unique project string so nothing collides.

# ---------------------------------------------------------------- path + empty

test_that("region_params_path lands in the project's _intermediate dir", {
    p <- region_params_path("RP_PATH")
    expect_true(grepl("RP_PATH_results", p, fixed = TRUE))
    expect_identical(basename(p), "region_params.json")
    expect_identical(basename(dirname(p)), "_intermediate")
})

test_that("read_region_params returns the empty skeleton when nothing is stored", {
    got <- read_region_params("RP_MISSING")
    expect_identical(got, list(global = list(), regions = list()))
})

test_that("read_region_params returns the skeleton for corrupt JSON rather than erroring", {
    path <- region_params_path("RP_CORRUPT")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("{not valid json", path)
    expect_identical(read_region_params("RP_CORRUPT"),
                     list(global = list(), regions = list()))
})

test_that("read_region_params treats an empty file as absent", {
    # file_ok() requires size > 0, so a zero-byte file never reaches read_json.
    path <- region_params_path("RP_EMPTY")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    file.create(path)
    expect_identical(read_region_params("RP_EMPTY"),
                     list(global = list(), regions = list()))
})

# ---------------------------------------------------------------- set/get, in memory

test_that("set_region_param builds every missing level of the nesting", {
    p <- set_region_param(.empty_region_params(), "gea", "5_100-200", "hap_scan",
                          list(epsilon = 0.6))
    expect_identical(p$regions$gea$`5_100-200`$hap_scan, list(epsilon = 0.6))
})

test_that("set_region_param does not disturb a sibling region or module", {
    p <- .empty_region_params()
    p <- set_region_param(p, "gea",  "r1", "hap_scan", list(epsilon = 0.6))
    p <- set_region_param(p, "gea",  "r2", "hap_scan", list(epsilon = 0.9))
    p <- set_region_param(p, "gwas", "r1", "hap_scan", list(epsilon = 0.1))
    expect_identical(get_region_param(p, "gea",  "r1", "hap_scan"), list(epsilon = 0.6))
    expect_identical(get_region_param(p, "gea",  "r2", "hap_scan"), list(epsilon = 0.9))
    expect_identical(get_region_param(p, "gwas", "r1", "hap_scan"), list(epsilon = 0.1))
})

test_that("set_region_param overwrites the same slot rather than appending", {
    p <- .empty_region_params()
    p <- set_region_param(p, "gea", "r1", "hap_scan", list(epsilon = 0.6))
    p <- set_region_param(p, "gea", "r1", "hap_scan", list(epsilon = 0.8))
    expect_identical(get_region_param(p, "gea", "r1", "hap_scan"), list(epsilon = 0.8))
    expect_length(p$regions$gea$r1, 1L)
})

test_that("set_region_param keeps two computation types on one region apart", {
    p <- .empty_region_params()
    p <- set_region_param(p, "gea", "r1", "hap_scan", list(epsilon = 0.6))
    p <- set_region_param(p, "gea", "r1", "hap_viz",  list(epsilon = 0.2))
    expect_identical(get_region_param(p, "gea", "r1", "hap_scan"), list(epsilon = 0.6))
    expect_identical(get_region_param(p, "gea", "r1", "hap_viz"),  list(epsilon = 0.2))
})

test_that("get_region_param returns NULL for any missing level", {
    p <- set_region_param(.empty_region_params(), "gea", "r1", "hap_scan", list(a = 1))
    expect_null(get_region_param(p, "gwas", "r1", "hap_scan"))   # module
    expect_null(get_region_param(p, "gea",  "r9", "hap_scan"))   # region
    expect_null(get_region_param(p, "gea",  "r1", "hap_viz"))    # type
})

test_that("get_region_param returns NULL when any argument is NULL", {
    p <- set_region_param(.empty_region_params(), "gea", "r1", "hap_scan", list(a = 1))
    expect_null(get_region_param(NULL, "gea", "r1", "hap_scan"))
    expect_null(get_region_param(p, NULL, "r1", "hap_scan"))
    expect_null(get_region_param(p, "gea", NULL, "hap_scan"))
    expect_null(get_region_param(p, "gea", "r1", NULL))
})

test_that("global params are keyed by module and independent of regions", {
    p <- set_global_param(.empty_region_params(), "gea", "region_distance", 50000)
    expect_identical(get_global_param(p, "gea", "region_distance"), 50000)
    expect_null(get_global_param(p, "gwas", "region_distance"))
    expect_null(get_global_param(p, "gea", "other_key"))
    expect_identical(p$regions, list())
})

test_that("get_global_param returns NULL for NULL params", {
    expect_null(get_global_param(NULL, "gea", "region_distance"))
})

# ---------------------------------------------------------------- disk round trip

test_that("save_region_params then read_region_params preserves a scalar param", {
    p <- set_region_param(.empty_region_params(), "gea", "r1", "hap_scan",
                          list(epsilon = 0.6))
    save_region_params("RP_RT", p)
    got <- read_region_params("RP_RT")
    expect_identical(get_region_param(got, "gea", "r1", "hap_scan"), list(epsilon = 0.6))
})

test_that("a length-1 vector survives the round trip as a scalar, a longer one as a list", {
    # auto_unbox = TRUE on write + simplifyVector = FALSE on read. This is the
    # reason the round trip is not identical() for vector values, and a caller
    # that stored c(1, 2) gets list(1, 2) back.
    p <- .empty_region_params()
    p <- set_region_param(p, "gea", "r1", "hap_scan",
                          list(one = 5, many = c(1, 2, 3)))
    save_region_params("RP_VEC", p)
    got <- get_region_param(read_region_params("RP_VEC"), "gea", "r1", "hap_scan")
    expect_length(got$many, 3L)
    expect_type(got$many, "list")
    expect_equal(unlist(got$many), c(1, 2, 3))
    expect_equal(got$one, 5)
})

test_that("a whole-number double comes back as an INTEGER, a fractional one stays double", {
    # JSON has one number type. write_json emits 5 for both 5L and 5, and
    # read_json picks the narrowest R type that fits — so a double with no
    # fractional part changes type across the round trip while 0.6 does not.
    #
    # This is the sharp edge for the "modified" badge in CLAUDE.md's UI rule:
    # that badge compares the current input against the SAVED params, and an
    # identical() comparison of a whole-number parameter against its stored copy
    # is FALSE purely because of storage type. Use == or all.equal there.
    p <- .empty_region_params()
    p <- set_region_param(p, "gea", "r1", "hap_scan",
                          list(whole = 5, fractional = 0.6, already_int = 7L))
    save_region_params("RP_NUMTYPE", p)
    got <- get_region_param(read_region_params("RP_NUMTYPE"), "gea", "r1", "hap_scan")

    expect_type(got$whole, "integer")
    expect_false(identical(got$whole, 5))       # the trap
    expect_equal(got$whole, 5)                  # but numerically intact

    expect_type(got$fractional, "double")
    expect_identical(got$fractional, 0.6)       # unchanged, so the rule is type-driven

    expect_type(got$already_int, "integer")
    expect_identical(got$already_int, 7L)
})

test_that("save_region_params creates the _intermediate directory if absent", {
    # The project name must be unique per INVOCATION, not merely per test block:
    # CLINEGO_TEST_ROOT is stable for the whole session, so a fixed name would
    # make this the one test in the file that cannot be run twice in one session
    # — the first run creates the directory the second asserts is absent.
    # Harmless under test_check(), which runs each file once, but it reports as a
    # spurious failure under any re-run harness (it showed up as a bogus
    # "mutation leaked" verdict during this file's falsification run).
    project <- basename(tempfile("RP_MKDIR_"))
    path <- region_params_path(project)
    expect_false(dir.exists(dirname(path)))
    save_region_params(project, .empty_region_params())
    expect_true(file.exists(path))
})

test_that("save_region_params leaves no tempfile beside the target", {
    # It writes to a tempfile in the same dir and renames. A leaked .json temp
    # would be picked up by nothing, but it would accumulate in a real results
    # tree on every save.
    save_region_params("RP_TMP", .empty_region_params())
    leftovers <- list.files(dirname(region_params_path("RP_TMP")),
                            pattern = "^file.*\\.json$")
    expect_identical(leftovers, character(0))
})

test_that("save_region_params returns its input whether or not the write worked", {
    p <- set_global_param(.empty_region_params(), "gea", "region_distance", 1000)
    expect_identical(save_region_params("RP_RET", p), p)
    # The return value carries no success signal: the error path only message()s,
    # so a caller cannot branch on it.
})

test_that("a later save replaces the stored params rather than merging them", {
    a <- set_global_param(.empty_region_params(), "gea", "region_distance", 111)
    save_region_params("RP_REPLACE", a)
    b <- set_global_param(.empty_region_params(), "gea", "other", 222)
    save_region_params("RP_REPLACE", b)
    got <- read_region_params("RP_REPLACE")
    expect_null(get_global_param(got, "gea", "region_distance"))
    # expect_equal, not expect_identical: 222 is stored as a double and read back
    # as an integer — see the number-type round-trip test above.
    expect_equal(get_global_param(got, "gea", "other"), 222)
})
