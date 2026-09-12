# R/fct_pipeline.R — the Snakemake run-log readers.
#
# pipeline_read_progress() is the largest untested parser in the app: a ~110-line
# single-pass state machine over Snakemake's log, driving the Run tab's progress
# bar and rule list. It is NOT extracted from its file reader, so these tests
# write a log file rather than passing lines — the extraction is filed as a
# request rather than done here.
#
# The invariant that does NOT hold, and is pinned below so nobody "fixes" it:
# `done` and `failed` are LOG-EVENT COUNTERS while `rules` is a JOBID REGISTRY.
# jobid 0 (the `all` DAG root) is counted in `done` but never registered, and a
# block that ends before its jobid line is dropped entirely — so
# done != sum(status == "done") by design.
#
# Two parsing rules are asymmetric and easy to reverse:
#   "N of M steps" is UNANCHORED and last-wins, taking nums[2] positionally over
#   the whole line; "^total N" is first-wins, guarded by is.na(total).

write_log <- function(lines) {
    project <- basename(tempfile("PL_"))
    path <- pipeline_run_log_path(project)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path)
    project
}

# A complete, well-formed rule block.
blk <- function(rule, jobid, log = NULL) {
    c(paste0("rule ", rule, ":"),
      "    input: a.vcf",
      if (!is.null(log)) paste0("    log: ", log),
      paste0("    jobid: ", jobid),
      "")
}

# ---------------------------------------------------------------- empty states

test_that("pipeline_read_progress returns the empty shape when the log is missing", {
    got <- pipeline_read_progress(basename(tempfile("PL_NONE_")))
    expect_identical(got$total, NA_integer_)
    expect_identical(got$done, 0L)
    expect_identical(got$failed, 0L)
    expect_identical(got$current_rule, "")
    expect_identical(got$rules, list())
})

test_that("pipeline_read_progress returns the empty shape for a zero-line log", {
    got <- pipeline_read_progress(write_log(character(0)))
    expect_identical(got$total, NA_integer_)
    expect_identical(got$rules, list())
})

# ---------------------------------------------------------------- block parsing

test_that("both 'rule' and 'localrule' open a block", {
    got <- pipeline_read_progress(write_log(c(blk("filter_vcf", 1), blk("copy_meta", 2))))
    expect_identical(length(got$rules), 2L)
    expect_identical(vapply(got$rules, function(r) r$rule, character(1)),
                     c("filter_vcf", "copy_meta"))

    got2 <- pipeline_read_progress(write_log(c("localrule copy_meta:", "    jobid: 7", "")))
    expect_identical(length(got2$rules), 1L)
    expect_identical(got2$rules[[1]]$rule, "copy_meta")
})

test_that("a registered rule starts as running and sets current_rule", {
    got <- pipeline_read_progress(write_log(blk("filter_vcf", 1)))
    expect_identical(got$rules[[1]]$status, "running")
    expect_identical(got$rules[[1]]$jobid, 1L)
    expect_identical(got$current_rule, "filter_vcf")
})

test_that("jobid 0 is NOT registered — it is the all-rule DAG root", {
    got <- pipeline_read_progress(write_log(c(blk("all", 0), blk("filter_vcf", 1))))
    expect_identical(length(got$rules), 1L)
    expect_identical(got$rules[[1]]$rule, "filter_vcf")
})

test_that("a block ending before its jobid line is silently DROPPED", {
    # Terminated by a blank line...
    got <- pipeline_read_progress(write_log(c(
        "rule orphan:", "    input: a.vcf", "", blk("filter_vcf", 1))))
    expect_identical(length(got$rules), 1L)
    expect_identical(got$rules[[1]]$rule, "filter_vcf")

    # ...and by a timestamp line, the other terminator.
    got2 <- pipeline_read_progress(write_log(c(
        "rule orphan:", "    input: a.vcf",
        "[Fri Sep 12 10:00:00 2026]", blk("filter_vcf", 1))))
    expect_identical(length(got2$rules), 1L)
})

test_that("the log: line keeps the first path and drops the rest after a comma", {
    got <- pipeline_read_progress(write_log(
        blk("filter_vcf", 1, log = "logs/a.log, logs/b.log")))
    expect_identical(got$rules[[1]]$log_file, "logs/a.log")
})

test_that("a rule with no log: line carries NA rather than an empty string", {
    got <- pipeline_read_progress(write_log(blk("filter_vcf", 1)))
    expect_true(is.na(got$rules[[1]]$log_file))
})

test_that("a repeated jobid updates in place instead of appending", {
    got <- pipeline_read_progress(write_log(c(blk("filter_vcf", 1), blk("filter_vcf", 1))))
    expect_identical(length(got$rules), 1L)
})

test_that("rules come back in first-seen order, not numeric jobid order", {
    got <- pipeline_read_progress(write_log(c(blk("c_rule", 9), blk("a_rule", 2))))
    expect_identical(vapply(got$rules, function(r) r$jobid, integer(1)), c(9L, 2L))
})

# ---------------------------------------------------------------- finished / error

test_that("Finished jobid marks the matching rule done", {
    got <- pipeline_read_progress(write_log(c(
        blk("filter_vcf", 1), "Finished jobid: 1 (Rule: filter_vcf)")))
    expect_identical(got$done, 1L)
    expect_identical(got$rules[[1]]$status, "done")
})

test_that("Finished jobid counts even for an id that was never registered", {
    # done is a log-event counter: jobid 0 and any unknown id still increment it.
    got <- pipeline_read_progress(write_log(c(
        blk("filter_vcf", 1),
        "Finished jobid: 0 (Rule: all)",
        "Finished jobid: 99 (Rule: ghost)")))
    expect_identical(got$done, 2L)
    expect_identical(got$rules[[1]]$status, "running")
})

test_that("done is deliberately NOT the count of rules marked done", {
    # The invariant that does not hold. Pinned so a future "fix" is a decision
    # rather than an accident.
    got <- pipeline_read_progress(write_log(c(
        blk("filter_vcf", 1),
        "Finished jobid: 0 (Rule: all)",
        "Finished jobid: 1 (Rule: filter_vcf)")))
    n_done <- sum(vapply(got$rules, function(r) r$status == "done", logical(1)))
    expect_identical(got$done, 2L)
    expect_identical(n_done, 1L)
    expect_false(identical(got$done, n_done))
})

test_that("Error in rule marks the LAST running instance of that rule", {
    got <- pipeline_read_progress(write_log(c(
        blk("emmax", 1), blk("emmax", 2), "Error in rule emmax:")))
    expect_identical(got$failed, 1L)
    expect_identical(got$rules[[1]]$status, "running")   # jobid 1 untouched
    expect_identical(got$rules[[2]]$status, "failed")    # jobid 2, the latest
})

test_that("Error in rule skips instances that are already done", {
    got <- pipeline_read_progress(write_log(c(
        blk("emmax", 1), blk("emmax", 2),
        "Finished jobid: 2 (Rule: emmax)",
        "Error in rule emmax:")))
    expect_identical(got$rules[[2]]$status, "done")
    expect_identical(got$rules[[1]]$status, "failed")
})

test_that("Error in rule counts even when it names an unregistered rule", {
    got <- pipeline_read_progress(write_log(c(blk("filter_vcf", 1), "Error in rule ghost:")))
    expect_identical(got$failed, 1L)
    expect_identical(got$rules[[1]]$status, "running")
})

# ---------------------------------------------------------------- total

test_that("'N of M steps' takes M and is LAST-wins", {
    got <- pipeline_read_progress(write_log(c(
        "1 of 12 steps (8%) done", "4 of 12 steps (33%) done")))
    expect_identical(got$total, 12L)

    # A later line with a different M overwrites the earlier one.
    got2 <- pipeline_read_progress(write_log(c(
        "1 of 12 steps (8%) done", "1 of 20 steps (5%) done")))
    expect_identical(got2$total, 20L)
})

test_that("'total N' is FIRST-wins and only applies when no steps line has landed", {
    got <- pipeline_read_progress(write_log(c("total       7", "total       9")))
    expect_identical(got$total, 7L)

    # A steps line already set total, so the later 'total' line is ignored.
    got2 <- pipeline_read_progress(write_log(c("2 of 12 steps (16%) done", "total       7")))
    expect_identical(got2$total, 12L)
})

test_that("'total N' seen FIRST is overwritten by a later steps line", {
    # Asymmetric on purpose: the steps branch has no is.na() guard.
    got <- pipeline_read_progress(write_log(c("total       7", "2 of 12 steps (16%) done")))
    expect_identical(got$total, 12L)
})

test_that("a timestamp prefix on the steps line corrupts total", {
    # The regex is UNANCHORED and nums[2] is positional over the WHOLE line, so
    # leading digits shift the window. Characterisation of a real fragility:
    # here total becomes 12 (from the time) rather than 20.
    got <- pipeline_read_progress(write_log("[Fri Sep 12 10:00:00 2026] 3 of 20 steps (15%) done"))
    expect_false(identical(got$total, 20L))
})

# ---------------------------------------------------------------- pipeline_cmd

test_that("pipeline_cmd builds the documented argument vector", {
    args <- pipeline_cmd("gea", "PROJ", 4, pipeline_path = "/pipeline")
    expect_identical(args[1:2], c("-c", "4"))
    expect_identical(args[3:4], c("-s", "/pipeline/Snakefile"))
    expect_identical(args[5:6], c("--config", "mode=gea"))
    expect_identical(args[9:10], c("--scheduler", "greedy"))
    expect_identical(args[11], "--rerun-incomplete")
})

test_that("pipeline_cmd prefers a per-project config when one exists on disk", {
    # Not a pure function: config_file_path() hits the filesystem, so both
    # branches need a real directory.
    root <- file.path(tempdir(), basename(tempfile("cfg_")))
    dir.create(root, recursive = TRUE, showWarnings = FALSE)

    args <- pipeline_cmd("gea", "PROJ", 2, pipeline_path = root)
    expect_identical(args[8], file.path(root, "config.yaml"))

    file.create(file.path(root, "config_PROJ.yaml"))
    args2 <- pipeline_cmd("gea", "PROJ", 2, pipeline_path = root)
    expect_identical(args2[8], file.path(root, "config_PROJ.yaml"))
})

# ---------------------------------------------------------------- log tail

test_that("pipeline_read_log_tail reports a missing log as a length-1 string", {
    # Indistinguishable in TYPE from a successful read — a caller cannot branch
    # on the return value alone.
    expect_identical(pipeline_read_log_tail(NULL), "Log file not found.")
    expect_identical(pipeline_read_log_tail(NA_character_), "Log file not found.")
    expect_identical(pipeline_read_log_tail(""), "Log file not found.")
    expect_identical(pipeline_read_log_tail(file.path(tempdir(), "nope.log")),
                     "Log file not found.")
})

test_that("pipeline_read_log_tail returns every line when under the limit", {
    p <- tempfile(fileext = ".log"); writeLines(c("a", "b", "c"), p)
    expect_identical(pipeline_read_log_tail(p, n = 10L), c("a", "b", "c"))
})

test_that("pipeline_read_log_tail prepends a notice and returns n + 1 lines when truncating", {
    p <- tempfile(fileext = ".log"); writeLines(as.character(1:10), p)
    got <- pipeline_read_log_tail(p, n = 3L)
    expect_length(got, 4L)
    expect_match(got[1], "^\\.\\.\\. \\(showing last 3 of 10 lines\\)")
    expect_identical(got[2:4], c("8", "9", "10"))
})

# ---------------------------------------------------------------- warnings

test_that("pipeline_read_warnings returns a typed zero-row frame when there is no log", {
    got <- pipeline_read_warnings(basename(tempfile("PL_NOWARN_")))
    expect_identical(nrow(got), 0L)
    expect_identical(names(got), c("message", "line"))
})

test_that("pipeline_read_warnings strips the prefix and reports whole-file line numbers", {
    project <- write_log(c("starting", "WARNING: low sample count", "rule x:",
                           "WARN: second one", "done"))
    got <- pipeline_read_warnings(project)
    expect_identical(got$message, c("low sample count", "second one"))
    expect_identical(got$line, c(2L, 4L))
})

test_that("pipeline_read_warnings matches case-insensitively but only at line start", {
    project <- write_log(c("warning: lowercase counts",
                           "a trailing WARNING: does not"))
    got <- pipeline_read_warnings(project)
    expect_identical(nrow(got), 1L)
    expect_identical(got$message, "lowercase counts")
})

test_that("pipeline_read_warnings returns character columns, not factors", {
    project <- write_log(c("WARNING: one"))
    expect_type(pipeline_read_warnings(project)$message, "character")
})
