#' Pipeline execution utilities
#'
#' Functions for launching Snakemake as a background process from within the
#' Shiny app, tracking progress by parsing .run.log, and handling cleanup.
#'
#' Both the Shiny app and Snakemake run inside the same Docker container. Two roots:
#' the code root (get_pipeline_path(): Snakefile, scripts/, workflow/) and the project
#' root (get_project_root(): config_<p>.yaml, <p>_results/, <p>_logs/, inputs). They
#' are one folder in the developer form (the repo mounted at /pipeline) and two in the
#' image, where the code is baked in read-only and the project is mounted at /data.
#' Every Snakemake call passes --config project_root=<project root> and runs with the
#' project root as its working directory. Snakemake is installed via pipx at
#' /root/.local/bin/snakemake (or on PATH as snakemake).
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
#' @param code_root code root (the Snakefile's folder)
#' @param project_root project root (the config, results and logs folder)
#' @return character vector suitable for processx::process$new(command, args)
#' @noRd
pipeline_cmd <- function(mode, project, cpu,
                          code_root = get_pipeline_path(),
                          project_root = get_project_root()) {
    configfile <- config_file_path(project, project_root)

    c(
        "-c", as.character(cpu),
        "-s", file.path(code_root, "Snakefile"),
        "--config", paste0("mode=", mode), pipeline_root_arg(project_root),
        "--configfile", configfile,
        "--scheduler", "greedy",
        "--rerun-incomplete"
    )
}

#' The `project_root=<dir>` word of a Snakemake --config, as an absolute path
#'
#' common.smk refuses a relative project_root, so a relative option or env value is
#' made absolute here (against the app's working directory) rather than failing there.
#' @noRd
pipeline_root_arg <- function(project_root) {
    paste0("project_root=", normalizePath(project_root, winslash = "/", mustWork = FALSE))
}

#' Launch Snakemake as a background process
#'
#' @param mode pipeline mode string
#' @param project project name
#' @param cpu number of CPU cores
#' @param code_root code root (the Snakefile's folder)
#' @param project_root project root (config, results, logs; the working directory)
#' @return processx::process object
#' @noRd
pipeline_run <- function(mode, project, cpu,
                          code_root = get_pipeline_path(),
                          project_root = get_project_root()) {
    log_dir <- file.path(project_root, paste0(project, "_logs"))
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(log_dir, ".run.log")

    args <- pipeline_cmd(mode, project, cpu, code_root, project_root)

    # Never clear a lock here: it may belong to a run this session cannot see -- a
    # command-line run, another browser tab, or this app's own run from before a page
    # reload (ac6ab8). The runner asks the user before pipeline_unlock() removes one.
    if (pipeline_is_locked(project, project_root))
        stop("the results folder is locked by another Snakemake run")

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
        wd      = project_root,
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
#' Safe to call even when no lock exists. Used by pipeline_kill() (Stop button) and
#' by the runner's "Remove lock" confirmation (a lock a prior crash or docker stop
#' left behind). mode is required by Snakemake to parse the Snakefile config block.
#'
#' @param project project name
#' @param mode pipeline mode string (e.g. "gea") — required for Snakefile config parsing
#' @param code_root code root (the Snakefile's folder)
#' @param project_root project root (config, results, logs; the working directory)
#' @noRd
pipeline_unlock <- function(project, mode, code_root = get_pipeline_path(),
                            project_root = get_project_root()) {
    tryCatch({
        smk <- Sys.which("snakemake")
        if (!nzchar(smk)) smk <- "/root/.local/bin/snakemake"
        if (nzchar(smk) && file.exists(smk)) {
            # In the project root like every other Snakemake call, never in the
            # code root (read-only in the image).
            processx::run(smk,
                          args    = pipeline_unlock_args(project, mode, code_root,
                                                         project_root),
                          wd      = project_root,
                          error_on_status = FALSE,
                          timeout = 10)
        }
    }, error = function(e) NULL)
    invisible(NULL)
}

#' The argument vector of pipeline_unlock()'s snakemake --unlock call
#' @noRd
pipeline_unlock_args <- function(project, mode, code_root = get_pipeline_path(),
                                 project_root = get_project_root()) {
    c("--unlock", "-s", file.path(code_root, "Snakefile"),
      "--config", paste0("mode=", mode), pipeline_root_arg(project_root),
      "--configfile", config_file_path(project, project_root))
}

#' Kill a running pipeline and unlock the Snakemake workdir
#'
#' @param proc processx::process object
#' @param project project name
#' @param mode pipeline mode string (passed through to pipeline_unlock)
#' @param code_root code root (the Snakefile's folder)
#' @param project_root project root (config, results, logs)
#' @noRd
pipeline_kill <- function(proc, project, mode, code_root = get_pipeline_path(),
                          project_root = get_project_root()) {
    if (pipeline_is_running(proc)) {
        proc$kill()
        proc$wait(timeout = 3000)
    }
    pipeline_unlock(project, mode, code_root, project_root)
    invisible(NULL)
}

#' Check whether a Snakemake lock exists for the project workdir
#'
#' @param project project name
#' @param project_root project root
#' @return logical
#' @noRd
pipeline_is_locked <- function(project, project_root = get_project_root()) {
    results_dir <- file.path(project_root, paste0(project, "_results"))
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
#' @param project_root project root
#' @noRd
pipeline_read_progress <- function(project, project_root = get_project_root()) {
    log_path <- pipeline_run_log_path(project, project_root)
    lines <- if (file.exists(log_path))
        tryCatch(readLines(log_path, warn = FALSE), error = function(e) character(0))
    else character(0)
    pipeline_parse_progress(lines)
}

#' Parse the lines of a Snakemake 9 run log into progress (see pipeline_read_progress()).
#'
#' Pure: no file access, so a test passes the lines themselves (f8a481).
#'
#' @param lines character vector, one log line each
#' @return list(total, done, failed, current_rule, rules)
#' @noRd
pipeline_parse_progress <- function(lines) {
    empty <- list(total = NA_integer_, done = 0L, failed = 0L,
                  current_rule = "", rules = list())
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
        # M is captured from the phrase itself: taking the line's second number
        # read a timestamp prefix's digits instead (a8501b).
        steps <- regmatches(line, regexec("(\\d+) of (\\d+) steps", line))[[1]]
        if (length(steps) == 3) total <- as.integer(steps[3])

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
#' @param project_root project root
#' @noRd
pipeline_run_log_path <- function(project, project_root = get_project_root()) {
    file.path(project_root, paste0(project, "_logs"), ".run.log")
}

#' Extract WARNING lines from the main run log
#'
#' Scans the .run.log for lines starting with WARNING: or WARN:.
#' Returns a data.frame with $message (stripped prefix) and $line (line number).
#' Used by mod_pipeline_runner to show toast notifications and the warning summary.
#'
#' @param project project name
#' @param project_root project root
#' @return data.frame with columns: message, line
#' @noRd
pipeline_read_warnings <- function(project, project_root = get_project_root()) {
    log_path <- pipeline_run_log_path(project, project_root)
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
