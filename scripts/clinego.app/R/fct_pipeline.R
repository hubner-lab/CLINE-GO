#' Pipeline execution utilities
#'
#' Functions for launching Snakemake as a background process from within the
#' Shiny app, tracking progress by parsing .run.log, and handling cleanup.
#'
#' Both the Shiny app and Snakemake run inside the same Docker container at
#' /pipeline. Snakemake is installed via pipx at /root/.local/bin/snakemake
#' (or on PATH as snakemake).
#'
#' Progress is tracked by parsing the Snakemake 9 log format written to
#' {project}_logs/.run.log. Key patterns:
#'   "rule name:"           — rule started (followed by log: and jobid: fields)
#'   "Finished jobid: N"    — rule completed
#'   "Error in rule name:"  — rule failed
#'   "N of M steps"         — overall progress

#' Build the snakemake command argument vector
#'
#' @param mode pipeline mode string (e.g. "association")
#' @param project project name
#' @param cpu number of CPU cores
#' @param pipeline_path pipeline root path
#' @return character vector suitable for processx::process$new(command, args)
#' @noRd
pipeline_cmd <- function(mode, project, cpu,
                          pipeline_path = get_pipeline_path()) {
    configfile <- config_file_path(project, pipeline_path)

    c(
        "-c", as.character(cpu),
        "-s", file.path(pipeline_path, "Snakefile"),
        "--config", paste0("mode=", mode),
        "--configfile", configfile,
        "--scheduler", "greedy",
        "--rerun-incomplete"
    )
}

#' Launch Snakemake as a background process
#'
#' @param mode pipeline mode string
#' @param project project name
#' @param cpu number of CPU cores
#' @param pipeline_path pipeline root path
#' @return processx::process object
#' @noRd
pipeline_run <- function(mode, project, cpu,
                          pipeline_path = get_pipeline_path()) {
    log_dir <- file.path(pipeline_path, paste0(project, "_logs"))
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(log_dir, ".run.log")

    args <- pipeline_cmd(mode, project, cpu, pipeline_path)

    # Auto-clear any stale lock left by a prior crash / docker stop / kill signal.
    # Safe: the app's global pipeline_running lock ensures no other in-app snakemake
    # is running by the time we reach here, so any on-disk lock is necessarily stale.
    if (pipeline_is_locked(project, pipeline_path)) {
        pipeline_unlock(project, mode, pipeline_path)
    }

    # Find snakemake on PATH (pipx installs to ~/.local/bin, also /root/.local/bin)
    smk <- Sys.which("snakemake")
    if (!nzchar(smk)) {
        for (candidate in c("/root/.local/bin/snakemake",
                            "/home/shiny/.local/bin/snakemake",
                            "/usr/local/bin/snakemake")) {
            if (file.exists(candidate)) { smk <- candidate; break }
        }
    }
    if (!nzchar(smk)) stop("snakemake not found on PATH")

    processx::process$new(
        command = smk,
        args    = args,
        stdout  = log_file,
        stderr  = "2>&1",
        wd      = pipeline_path,
        cleanup = FALSE   # keep process alive after Shiny session ends
    )
}

#' Check whether a pipeline process is still running
#'
#' @param proc processx::process object (or NULL)
#' @return logical
#' @noRd
pipeline_is_running <- function(proc) {
    !is.null(proc) && inherits(proc, "process") && proc$is_alive()
}

#' Run snakemake --unlock for a project's working directory
#'
#' Safe to call even when no lock exists. Used both by pipeline_kill() (Stop button)
#' and pipeline_run() (auto-clear stale lock from a prior crash/docker stop).
#' mode is required by Snakemake to parse the Snakefile config block.
#'
#' @param project project name
#' @param mode pipeline mode string (e.g. "gea") — required for Snakefile config parsing
#' @param pipeline_path pipeline root path
#' @noRd
pipeline_unlock <- function(project, mode, pipeline_path = get_pipeline_path()) {
    tryCatch({
        smk <- Sys.which("snakemake")
        if (!nzchar(smk)) smk <- "/root/.local/bin/snakemake"
        if (nzchar(smk) && file.exists(smk)) {
            system2(smk,
                    args    = c("--unlock", "-s",
                                file.path(pipeline_path, "Snakefile"),
                                "--config", paste0("mode=", mode),
                                "--configfile",
                                config_file_path(project, pipeline_path)),
                    stdout  = FALSE,
                    stderr  = FALSE,
                    timeout = 10)
        }
    }, error = function(e) NULL)
    invisible(NULL)
}

#' Kill a running pipeline and unlock the Snakemake workdir
#'
#' @param proc processx::process object
#' @param project project name
#' @param mode pipeline mode string (passed through to pipeline_unlock)
#' @param pipeline_path pipeline root path
#' @noRd
pipeline_kill <- function(proc, project, mode, pipeline_path = get_pipeline_path()) {
    if (pipeline_is_running(proc)) {
        proc$kill()
        proc$wait(timeout = 3000)
    }
    pipeline_unlock(project, mode, pipeline_path)
    invisible(NULL)
}

#' Check whether a Snakemake lock exists for the project workdir
#'
#' @param project project name
#' @param pipeline_path pipeline root path
#' @return logical
#' @noRd
pipeline_is_locked <- function(project, pipeline_path = get_pipeline_path()) {
    results_dir <- file.path(pipeline_path, paste0(project, "_results"))
    lock_dir    <- file.path(results_dir, ".snakemake", "locks")
    if (!dir.exists(lock_dir)) return(FALSE)
    length(list.files(lock_dir, all.files = TRUE, no.. = TRUE)) > 0
}

#' Read and parse progress from the Snakemake run log (.run.log)
#'
#' Parses the Snakemake 9 log format written to {project}_logs/.run.log.
#'
#' Returns a list with:
#'   $total        — total job count (from "N of M steps" or "total N" job stats)
#'   $done         — jobs finished
#'   $failed       — jobs failed
#'   $current_rule — most recently started rule name
#'   $rules        — ordered list of rule summaries, each:
#'                   list(rule, jobid, status, log_file)
#'                   status: "running" | "done" | "failed"
#'
#' @param project project name
#' @param pipeline_path pipeline root path
#' @noRd
pipeline_read_progress <- function(project, pipeline_path = get_pipeline_path()) {
    log_path <- pipeline_run_log_path(project, pipeline_path)

    empty <- list(total = NA_integer_, done = 0L, failed = 0L,
                  current_rule = "", rules = list())

    if (!file.exists(log_path)) return(empty)

    lines <- tryCatch(readLines(log_path, warn = FALSE), error = function(e) character(0))
    if (length(lines) == 0) return(empty)

    total        <- NA_integer_
    done         <- 0L
    failed       <- 0L
    current_rule <- ""

    # Rule tracking: key = jobid (character), value = list(rule, jobid, status, log_file)
    rule_keys <- character(0)
    rule_map  <- list()

    # Parser state for multi-line rule blocks
    in_block      <- FALSE
    cur_rule      <- ""
    cur_jobid_str <- ""
    cur_log_file  <- NA_character_

    for (line in lines) {

        # ── Rule block start: "rule name:" or "localrule name:" ─────────────────
        if (grepl("^(rule|localrule) \\w+:", line)) {
            cur_rule      <- sub("^(?:rule|localrule) (\\w+):.*", "\\1", line)
            cur_jobid_str <- ""
            cur_log_file  <- NA_character_
            in_block      <- TRUE
            next
        }

        if (in_block) {
            # "    log: /path/to/log" — capture first path before any comma
            if (grepl("^\\s+log:", line)) {
                log_raw  <- trimws(sub("^\\s+log:\\s*", "", line))
                cur_log_file <- trimws(strsplit(log_raw, ",")[[1]][1])
            }
            # "    jobid: N" — register the rule once we know its id
            if (grepl("^\\s+jobid:", line)) {
                cur_jobid_str <- trimws(sub("^\\s+jobid:\\s*", "", line))
                # Skip "all" rule (jobid 0) — it's the DAG root, not a real task
                if (cur_jobid_str != "0" && nzchar(cur_rule)) {
                    key <- cur_jobid_str
                    if (!key %in% rule_keys) rule_keys <- c(rule_keys, key)
                    rule_map[[key]] <- list(
                        rule     = cur_rule,
                        jobid    = as.integer(cur_jobid_str),
                        status   = "running",
                        log_file = cur_log_file
                    )
                    current_rule <- cur_rule
                }
                in_block <- FALSE
            }
            # Blank line or new timestamp ends block without finding jobid
            if (grepl("^\\s*$", line) || grepl("^\\[", line)) {
                in_block <- FALSE
            }
            next
        }

        # ── "Finished jobid: N (Rule: name)" ────────────────────────────────────
        if (grepl("^Finished jobid:", line)) {
            key <- sub("^Finished jobid: (\\d+).*", "\\1", line)
            done <- done + 1L          # always count, including jobid 0 (all rule)
            if (key %in% rule_keys) {
                rule_map[[key]]$status <- "done"
            }
        }

        # ── "Error in rule name:" ────────────────────────────────────────────────
        if (grepl("^Error in rule ", line)) {
            err_rule <- sub("^Error in rule (\\w+):.*", "\\1", line)
            failed <- failed + 1L      # always count
            for (key in rev(rule_keys)) {
                if (!is.null(rule_map[[key]]) &&
                    rule_map[[key]]$rule == err_rule &&
                    rule_map[[key]]$status == "running") {
                    rule_map[[key]]$status <- "failed"
                    break
                }
            }
        }

        # ── "N of M steps (P%) done" → update total ─────────────────────────────
        if (grepl("\\d+ of \\d+ steps", line)) {
            nums <- regmatches(line, gregexpr("\\d+", line))[[1]]
            if (length(nums) >= 2) total <- as.integer(nums[2])
        }

        # ── "total  N" from Job stats table ─────────────────────────────────────
        if (grepl("^total\\s+\\d+", line) && is.na(total)) {
            total <- as.integer(trimws(sub("^total\\s+", "", line)))
        }
    }

    rules <- lapply(rule_keys, function(k) rule_map[[k]])
    list(total = total, done = done, failed = failed,
         current_rule = current_rule, rules = rules)
}

#' Read the last N lines of a pipeline log file
#'
#' @param log_path absolute path to the log file
#' @param n number of lines to show from the end (default 200)
#' @return character vector of lines (may include a leading truncation notice)
#' @noRd
pipeline_read_log_tail <- function(log_path, n = 200L) {
    if (is.null(log_path) || is.na(log_path) || !nzchar(log_path) ||
        !file.exists(log_path)) {
        return("Log file not found.")
    }
    lines <- tryCatch(
        readLines(log_path, warn = FALSE),
        error = function(e) paste("Error reading log:", conditionMessage(e))
    )
    if (length(lines) > n) {
        c(paste0("... (showing last ", n, " of ", length(lines), " lines) ..."),
          tail(lines, n))
    } else {
        lines
    }
}

#' Path to the main .run.log for a project
#'
#' @param project project name
#' @param pipeline_path pipeline root path
#' @noRd
pipeline_run_log_path <- function(project, pipeline_path = get_pipeline_path()) {
    file.path(pipeline_path, paste0(project, "_logs"), ".run.log")
}

#' Extract WARNING lines from the main run log
#'
#' Scans the .run.log for lines starting with WARNING: or WARN:.
#' Returns a data.frame with $message (stripped prefix) and $line (line number).
#' Used by mod_pipeline_runner to show toast notifications and the warning summary.
#'
#' @param project project name
#' @param pipeline_path pipeline root path
#' @return data.frame with columns: message, line
#' @noRd
pipeline_read_warnings <- function(project, pipeline_path = get_pipeline_path()) {
    log_path <- pipeline_run_log_path(project, pipeline_path)
    if (!file.exists(log_path)) return(data.frame(message = character(0), line = integer(0)))

    lines <- tryCatch(readLines(log_path, warn = FALSE), error = function(e) character(0))
    if (length(lines) == 0) return(data.frame(message = character(0), line = integer(0)))

    warn_idx <- grep("^WARNING:|^WARN:", lines, ignore.case = TRUE)
    if (length(warn_idx) == 0) return(data.frame(message = character(0), line = integer(0)))

    data.frame(
        message = sub("^WARNING:\\s*|^WARN:\\s*", "", lines[warn_idx], ignore.case = TRUE),
        line    = warn_idx,
        stringsAsFactors = FALSE
    )
}
