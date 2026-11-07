#!/usr/bin/env Rscript
# collect_resources.R: wall time and peak memory of every pipeline mode run of the three Arabidopsis
# projects, as manuscript/arabidopsis/run_mode.sh recorded them.
#
# Reads ANALYSIS/Arabidopsis1135/resources/modes.tsv, one row per finished run_mode.sh call, and for
# each row the files run_mode.sh wrote under the same <config stem>_<mode>_<stamp> name: the memory
# trace (resources/*.tsv, one sample every 15 s), the processes above 1 GiB (resources/*_procs.tsv)
# and the Snakemake log (runlogs/*.log). Sample and SNP counts come from each project's
# Processing/tables/filtering_summary.tsv. Writes ANALYSIS/Arabidopsis1135/time_memory/:
#   time_memory.tsv  one row per real run, completed or failed: steps, wall time, peak memory, the
#                    script that held most of it at the peak, and the largest single process
#   skipped.tsv      the other rows of modes.tsv and why: dry runs (an -n token in the arguments, the
#                    failed `-- -n` call included) and the pilot's pregea, which the user stopped
#   coverage.tsv     each project's planned modes, measured or not, and why not
#   totals.tsv       wall time and largest peak per project, over its chain (each mode once) and over all
#                    its runs
#   time_memory.png  draft figure: wall time and peak memory per mode, memory over each project's runs
#
# Units and limits of the record:
#   - Memory is GiB: run_mode.sh divides bytes by 2^30, though modes.tsv names its columns _gb. anon is
#     the container's anonymous memory, what its jobs hold; cgroup is memory.current, the whole charge,
#     page cache included.
#   - A peak is the largest of samples taken every 15 s, so it is a lower bound. A run shorter than the
#     interval may have no sample: its peak is NA, not 0.
#   - Wall time runs from the launcher's start to the first poll after the run exits, so it can be up
#     to 15 s long.
#   - host_used_peak_gib is the whole host's used memory (free -g), other users' jobs included: context,
#     not this project's cost.
#
# Run from the repo root, in the image:
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp -v $PWD:/data -w /data \
#     cline-go:6cdedcc Rscript manuscript/arabidopsis/collect_resources.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(here)
})
source("/pipeline/scripts/R/utils/theme_clinego.R")

analysis_dir <- here("ANALYSIS", "Arabidopsis1135")
resources_dir <- file.path(analysis_dir, "resources")
runlogs_dir <- file.path(analysis_dir, "runlogs")
out_dir <- file.path(analysis_dir, "time_memory")
dir.create(out_dir, showWarnings = FALSE)

GIB <- 2^30
KIB_PER_GIB <- 2^20
SAMPLING_SECONDS <- 15
MAX_STAMP_LAG <- 60  # seconds from a file's stamp to the start modes.tsv records for its run

MODES <- c("processing", "prestructure", "structure", "climate", "traits", "pregea", "gea", "gwas",
           "gea_x_gwas", "maladaptation")

# The three projects, in order of filtered matrix size. The stress test plans five modes (Phase 16).
projects <- data.table(
  config = c("config_Arabidopsis1135chr5", "config_Arabidopsis1135", "config_Arabidopsis1135stress"),
  project = c("Arabidopsis1135chr5", "Arabidopsis1135", "Arabidopsis1135stress"),
  role = c("pilot, chromosome 5", "main", "stress test"))
planned <- rbind(
  data.table(project = "Arabidopsis1135chr5", mode = MODES),
  data.table(project = "Arabidopsis1135", mode = MODES),
  data.table(project = "Arabidopsis1135stress",
             mode = c("processing", "prestructure", "structure", "gea", "gwas")))

# Facts from the plan's run log that no file states. Every key must match exactly one row.
# A run the user stopped: left out of the table, with the reason.
stopped <- data.table(
  config = "config_Arabidopsis1135chr5", mode = "pregea", start = "2026-11-05T21:59:53+02:00",
  reason = "stopped by the user in pregea_rda_setup on 2026-11-06 to move to the whole dataset; not a failure")
# What a run covered, where a mode ran twice for real or a setting shaped its time. A run outside the
# chain is a second pass, not part of running each mode once.
run_notes <- data.table(
  config = c("config_Arabidopsis1135chr5", "config_Arabidopsis1135chr5", "config_Arabidopsis1135",
             "config_Arabidopsis1135", "config_Arabidopsis1135", "config_Arabidopsis1135stress"),
  mode = c("processing", "processing", "pregea", "maladaptation", "maladaptation", "gea"),
  start = c("2026-11-05T19:59:57+02:00", "2026-11-05T20:19:53+02:00", "2026-11-06T14:19:26+02:00",
            "2026-11-07T13:11:19+02:00", "2026-11-07T15:36:04+02:00", "2026-11-07T21:50:49+02:00"),
  note = c("whole mode at Filter.relatedness 0.99",
           "rerun from relatedness_filter at Filter.relatedness 0.998 (attempt 1)",
           "RDA ladder at 99 permutations on 32 forks (attempt 3)",
           "3 offset methods x 6 SNP sets x 30 gardens",
           "20 MAF-matched random draws of the union set; the 18 real runs did not rerun",
           "LEA geno2lfmm() segfaulted (exit 139) on the full matrix"),
  label = c("processing", "processing, rerun", "pregea", "maladaptation", "maladaptation, random draws",
            "gea"),
  in_chain = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE))
# Planned modes with no real run.
not_measured <- data.table(
  project = c(rep("Arabidopsis1135chr5", 5), "Arabidopsis1135stress"),
  mode = c("pregea", "gea", "gwas", "gea_x_gwas", "maladaptation", "gwas"),
  why = c("stopped by the user in the RDA ladder",
          rep("not run: the user stopped the pilot at pregea and moved to the whole dataset", 4),
          "not run: it needs vcf_to_lfmm_full, which failed in gea"))

# ------------------------------------------------------------------ modes.tsv and its files
modes <- fread(file.path(resources_dir, "modes.tsv"), sep = "\t", colClasses = "character")
modes[, `:=`(wall_seconds = as.numeric(wall_seconds), peak_anon_gib = as.numeric(peak_anon_gb),
             peak_cgroup_gib = as.numeric(peak_current_gb), exit = as.integer(exit))]
modes[, c("peak_anon_gb", "peak_current_gb") := NULL]
if (!all(modes$config %in% projects$config)) stop("collect_resources.R: modes.tsv names an unknown config")
modes <- merge(modes, projects, by = "config", sort = FALSE)
modes[, dry_run := vapply(strsplit(trimws(snakemake_args), "\\s+"), function(tok) "-n" %in% tok, logical(1))]

offsets <- unique(sub("^.*([+-][0-9]{2}:[0-9]{2})$", "\\1", modes$start))
if (length(offsets) != 1) stop("collect_resources.R: starts carry more than one UTC offset")
# Stamps are the host's local time, as starts are, so both compare as clock readings.
clock_seconds <- function(x, format) as.numeric(as.POSIXct(x, format = format, tz = "UTC"))
modes[, start_clock := clock_seconds(substr(start, 1, 19), "%Y-%m-%dT%H:%M:%S")]
modes[, start_epoch := as.numeric(as.POSIXct(sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", start),
                                             format = "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"))]

trace_files <- list.files(resources_dir, pattern = "^config_.*_[0-9]{8}_[0-9]{6}\\.tsv$")
traces <- data.table(base = sub("\\.tsv$", "", trace_files))
traces[, `:=`(prefix = sub("_[0-9]{8}_[0-9]{6}$", "", base),
              stamp_clock = clock_seconds(sub("^.*_([0-9]{8}_[0-9]{6})$", "\\1", base), "%Y%m%d_%H%M%S"))]
# A row takes the latest stamp at or before its start; one file per row and one row per file.
modes[, base := vapply(seq_len(.N), function(i) {
  hit <- traces[prefix == paste0(modes$config[i], "_", modes$mode[i]) &
                  modes$start_clock[i] - stamp_clock >= 0 &
                  modes$start_clock[i] - stamp_clock <= MAX_STAMP_LAG]
  if (nrow(hit) == 0) {
    stop(sprintf("collect_resources.R: no trace file for %s %s at %s", modes$config[i], modes$mode[i],
                 modes$start[i]))
  }
  hit[which.max(stamp_clock), base]
}, character(1))]
if (anyDuplicated(modes$base)) stop("collect_resources.R: one trace file matches two rows of modes.tsv")
if (!setequal(modes$base, traces$base)) {
  stop("collect_resources.R: trace files with no row in modes.tsv: ",
       paste(setdiff(traces$base, modes$base), collapse = ", "))
}

# Steps from the Snakemake log: the first job table's total, the last "k of n steps" line, failed rules.
log_steps <- function(base) {
  lines <- readLines(file.path(runlogs_dir, paste0(base, ".log")), warn = FALSE)
  total <- grep("^total\\s+[0-9]+$", lines, value = TRUE)
  done <- grep("^[0-9]+ of [0-9]+ steps", lines, value = TRUE)
  failed <- unique(sub("^Error in rule (\\S+):$", "\\1", grep("^Error in rule \\S+:$", lines, value = TRUE)))
  list(steps_planned = if (length(total)) as.integer(sub("^total\\s+", "", total[1])) else 0L,
       steps_done = if (length(done)) as.integer(sub(" of .*$", "", done[length(done)])) else 0L,
       failed_rule = paste(failed, collapse = ","))
}
modes[, c("steps_planned", "steps_done", "failed_rule") := rbindlist(lapply(base, log_steps))]

# ------------------------------------------------------------------ which rows are real runs
key_rows <- function(keys, what) {
  idx <- vapply(seq_len(nrow(keys)), function(i) {
    hit <- which(modes$config == keys$config[i] & modes$mode == keys$mode[i] & modes$start == keys$start[i])
    if (length(hit) != 1) stop(sprintf("collect_resources.R: %s key %s %s %s matches %d rows", what,
                                       keys$config[i], keys$mode[i], keys$start[i], length(hit)))
    hit
  }, integer(1))
  idx
}
stopped_rows <- key_rows(stopped, "stopped")
modes[, skip_reason := fifelse(dry_run, fifelse(exit == 0L, "dry run", sprintf("dry run, exit %d", exit)), "")]
modes[stopped_rows, skip_reason := sprintf("%s (exit %d)", stopped$reason, exit)]
if (any(modes$dry_run[stopped_rows])) stop("collect_resources.R: a stopped run is a dry run")

skipped <- modes[skip_reason != "", .(project, mode, start, wall_seconds, exit, steps_planned, steps_done,
                                     peak_anon_gib, snakemake_args, reason = skip_reason)]
kept <- modes[skip_reason == ""]
kept[, status := fifelse(exit == 0L, "completed", "failed")]
if (!all(paste(kept$project, kept$mode) %in% paste(planned$project, planned$mode))) {
  stop("collect_resources.R: a real run of a mode outside its project's plan")
}

note_rows <- match(key_rows(run_notes, "run note"), which(modes$skip_reason == ""))
if (anyNA(note_rows)) stop("collect_resources.R: a run note names a skipped row")
kept[, `:=`(note = "", label = mode, in_chain = TRUE)]
kept[note_rows, `:=`(note = run_notes$note, label = run_notes$label, in_chain = run_notes$in_chain)]
if (kept[, anyDuplicated(paste(project, label))]) stop("collect_resources.R: two runs share one label")
if (kept[status == "failed" & note == "", .N] > 0) stop("collect_resources.R: a failed run has no note")

# ------------------------------------------------------------------ memory: trace and processes
read_trace <- function(base) {
  fread(file.path(resources_dir, paste0(base, ".tsv")), sep = "\t",
        colClasses = c(time = "numeric", anon_bytes = "numeric", current_bytes = "numeric",
                       host_used_gb = "numeric"))
}
process_name <- function(command) {
  vapply(strsplit(trimws(command), "\\s+"), function(tok) {
    file_arg <- grep("^--file=", tok, value = TRUE)
    if (length(file_arg)) return(basename(sub("^--file=", "", file_arg[1])))
    first <- basename(tok[1])
    if (grepl("^(python[0-9.]*|bash|sh)$", first) && length(tok) > 1) return(basename(tok[2]))
    first
  }, character(1))
}
memory_of <- function(base, anon_recorded, cgroup_recorded) {
  trace <- read_trace(base)
  if (nrow(trace) == 0) {
    if (anon_recorded != 0 || cgroup_recorded != 0) stop("collect_resources.R: peaks with no sample, ", base)
    anon <- NA_real_
    cgroup <- NA_real_
    host <- NA_real_
  } else {
    anon <- max(trace$anon_bytes) / GIB
    cgroup <- max(trace$current_bytes) / GIB
    host <- max(trace$host_used_gb)
    if (abs(round(anon, 1) - anon_recorded) > 0.051 || abs(round(cgroup, 1) - cgroup_recorded) > 0.051) {
      stop("collect_resources.R: trace peaks differ from modes.tsv, ", base)
    }
  }
  procs <- fread(file.path(resources_dir, paste0(base, "_procs.tsv")), sep = "\t", quote = "",
                 colClasses = c(time = "numeric", pid = "integer", rss_kb = "numeric", command = "character"))
  procs[, name := process_name(command)]
  largest <- if (nrow(procs)) procs[which.max(rss_kb)] else NULL
  # What held the memory at the peak: the script with the largest summed RSS in the poll that took the
  # peak sample (run_mode.sh stamps a poll's trace row and process rows with one time). Forks share
  # pages, so a sum of their RSS can pass the container's anon.
  holder <- NULL
  if (nrow(trace) > 0 && nrow(procs) > 0) {
    peak_poll <- trace$time[which.max(trace$anon_bytes)]
    at_peak <- procs[time == peak_poll]
    if (nrow(at_peak)) holder <- at_peak[, .(n = .N, gib = sum(rss_kb) / KIB_PER_GIB), by = name][which.max(gib)]
  }
  list(peak_anon_gib = anon, peak_cgroup_gib = cgroup, memory_samples = nrow(trace),
       host_used_peak_gib = host,
       peak_holder = if (is.null(holder)) NA_character_ else holder$name,
       peak_holder_processes = if (is.null(holder)) NA_integer_ else holder$n,
       peak_holder_gib = if (is.null(holder)) NA_real_ else holder$gib,
       largest_process = if (is.null(largest)) NA_character_ else largest$name,
       largest_process_gib = if (is.null(largest)) NA_real_ else largest$rss_kb / KIB_PER_GIB)
}
kept[, c("peak_anon_gib", "peak_cgroup_gib", "memory_samples", "host_used_peak_gib", "peak_holder",
         "peak_holder_processes", "peak_holder_gib", "largest_process", "largest_process_gib") :=
       rbindlist(Map(memory_of, base, peak_anon_gib, peak_cgroup_gib))]

# ------------------------------------------------------------------ project sizes
size_of <- function(project) {
  filtering <- fread(here(paste0(project, "_results"), "Processing", "tables", "filtering_summary.tsv"))
  row <- filtering[stage == "After MAF + SNP missingness filter"]
  if (nrow(row) != 1) stop("collect_resources.R: no filtered row in ", project, "'s filtering_summary.tsv")
  list(n_samples = row$n_samples, n_snps = row$n_snps)
}
projects[, c("n_samples", "n_snps") := rbindlist(lapply(project, size_of))]
projects[, role_label := sprintf("%s: %s samples x %s SNPs", role, formatC(n_samples, format = "d", big.mark = ","),
                                 formatC(n_snps, format = "d", big.mark = ","))]
kept <- merge(kept, projects[, .(project, n_samples, n_snps, role_label)], by = "project", sort = FALSE)

# ------------------------------------------------------------------ tables
hms <- function(s) sprintf("%d:%02d:%02d", as.integer(s %/% 3600), as.integer(s %% 3600 %/% 60),
                           as.integer(round(s %% 60)))
kept[, project_order := match(project, projects$project)]
kept[, label_order := match(mode, MODES) + fifelse(in_chain, 0, 0.5)]
setorder(kept, project_order, label_order, start)
kept[, wall := hms(wall_seconds)]

time_memory <- kept[, .(project, role, n_samples, n_snps, mode, label, note, in_chain, start, end,
                        wall_seconds, wall, steps_planned, steps_done, exit, status, failed_rule,
                        peak_anon_gib = round(peak_anon_gib, 1), peak_cgroup_gib = round(peak_cgroup_gib, 1),
                        memory_samples, host_used_peak_gib, peak_holder, peak_holder_processes,
                        peak_holder_gib = round(peak_holder_gib, 1), largest_process,
                        largest_process_gib = round(largest_process_gib, 1), snakemake_args)]
fwrite(time_memory, file.path(out_dir, "time_memory.tsv"), sep = "\t")
fwrite(skipped, file.path(out_dir, "skipped.tsv"), sep = "\t")

coverage <- merge(planned, kept[, .(completed = sum(status == "completed"), failed = sum(status == "failed")),
                                by = .(project, mode)], by = c("project", "mode"), all.x = TRUE)
coverage[is.na(completed), `:=`(completed = 0L, failed = 0L)]
coverage <- merge(coverage, not_measured, by = c("project", "mode"), all.x = TRUE)
if (coverage[completed + failed == 0 & is.na(why), .N] > 0) stop("collect_resources.R: an unmeasured mode has no reason")
if (coverage[completed + failed > 0 & !is.na(why), .N] > 0) stop("collect_resources.R: a reason names a measured mode")
coverage[, status := fifelse(completed > 0, "completed", fifelse(failed > 0, "failed", "not measured"))]
coverage[is.na(why), why := ""]
coverage[, `:=`(project_order = match(project, projects$project), mode_order = match(mode, MODES))]
setorder(coverage, project_order, mode_order)
fwrite(coverage[, .(project, mode, status, completed_runs = completed, failed_runs = failed, why)],
       file.path(out_dir, "coverage.tsv"), sep = "\t")

totals <- kept[, .(runs = .N, chain_runs = sum(in_chain), chain_wall_seconds = sum(wall_seconds[in_chain]),
                   all_wall_seconds = sum(wall_seconds),
                   chain_max_peak_anon_gib = round(max(peak_anon_gib[in_chain], na.rm = TRUE), 1),
                   chain_max_peak_run = label[in_chain][which.max(peak_anon_gib[in_chain])],
                   max_peak_anon_gib = round(max(peak_anon_gib, na.rm = TRUE), 1),
                   max_peak_run = label[which.max(peak_anon_gib)]),
               by = .(project, role, project_order)]
totals[, `:=`(chain_wall = hms(chain_wall_seconds), all_wall = hms(all_wall_seconds))]
setorder(totals, project_order)
totals[, project_order := NULL]
fwrite(totals, file.path(out_dir, "totals.tsv"), sep = "\t")

# ------------------------------------------------------------------ figure (draft)
label_levels <- unique(kept[order(label_order), label])
kept[, y := match(label, rev(label_levels)) + (2 - project_order) * 0.27]
kept[, role_label := factor(role_label, levels = projects$role_label)]
project_colors <- setNames(CLINEGO_CATEGORICAL[seq_len(nrow(projects))], projects$role_label)
bar_half <- 0.12
y_axis <- scale_y_continuous(breaks = seq_along(label_levels), labels = rev(label_levels),
                             limits = c(0.5, length(label_levels) + 0.5), expand = c(0, 0))

time_floor <- 5
time_breaks <- c(10, 60, 600, 3600, 36000)
plot_time <- ggplot(kept) +
  geom_rect(aes(xmin = time_floor, xmax = wall_seconds, ymin = y - bar_half, ymax = y + bar_half,
                fill = role_label)) +
  geom_text(data = kept[status == "failed"],
            aes(x = wall_seconds, y = y, label = sprintf("failed, %d of %d steps", steps_done, steps_planned)),
            hjust = -0.08, size = 2.8, color = CLINEGO_REMOVED) +
  scale_x_log10(breaks = time_breaks, labels = c("10 s", "1 min", "10 min", "1 h", "10 h"),
                limits = c(time_floor, 2e5), expand = c(0, 0)) +
  y_axis +
  scale_fill_manual(values = project_colors, name = NULL) +
  labs(x = "Wall time", y = NULL) +
  theme_clinego()

memory_max <- max(kept$peak_cgroup_gib, na.rm = TRUE)
sampled <- kept[!is.na(peak_anon_gib)]
plot_memory <- ggplot(sampled) +
  geom_rect(aes(xmin = 0, xmax = peak_anon_gib, ymin = y - bar_half, ymax = y + bar_half, fill = role_label)) +
  geom_point(aes(x = peak_cgroup_gib, y = y, color = role_label), shape = 21, fill = "white", size = 1.6,
             stroke = 0.7) +
  geom_text(data = sampled[!is.na(peak_holder)],
            aes(x = pmax(peak_cgroup_gib, peak_anon_gib), y = y,
                label = sprintf("%s%s %.0f", peak_holder,
                                fifelse(peak_holder_processes > 1, sprintf(" x%d", peak_holder_processes), ""),
                                peak_holder_gib)),
            hjust = -0.12, size = 2.4, color = CLINEGO_COL$secondary) +
  geom_text(data = kept[is.na(peak_anon_gib)], aes(x = 0, y = y, label = "not sampled"),
            hjust = -0.1, size = 2.4, color = CLINEGO_COL$muted) +
  scale_x_continuous(limits = c(0, memory_max * 1.35), expand = c(0, 0)) +
  y_axis +
  scale_fill_manual(values = project_colors, name = NULL) +
  scale_color_manual(values = project_colors, guide = "none") +
  labs(x = "Peak memory, GiB", y = NULL) +
  theme_clinego()

traces_long <- rbindlist(lapply(seq_len(nrow(kept)), function(i) {
  trace <- read_trace(kept$base[i])
  if (nrow(trace) == 0) return(NULL)
  trace[, .(project = kept$project[i], label = kept$label[i], start = kept$start[i],
            elapsed = time - kept$start_epoch[i], anon_gib = anon_bytes / GIB,
            cgroup_gib = current_bytes / GIB)]
}))
bands <- kept[order(project_order, start), .(label, start, wall_seconds, role_label,
                                              band_start = cumsum(shift(wall_seconds, fill = 0))),
              by = project]
bands[, `:=`(band_end = band_start + wall_seconds, shade = seq_len(.N) %% 2 == 0), by = project]
bands[, show_label := wall_seconds / sum(wall_seconds) >= 0.03, by = project]
bands[show_label == TRUE, label_vjust := fifelse(seq_len(.N) %% 2 == 1, 1.4, 3.2), by = project]
traces_long <- merge(traces_long, bands[, .(project, start, band_start, role_label)], by = c("project", "start"))
traces_long[, hours := (band_start + elapsed) / 3600]
plot_trace <- ggplot() +
  geom_rect(data = bands[shade == TRUE],
            aes(xmin = band_start / 3600, xmax = band_end / 3600, ymin = -Inf, ymax = Inf),
            fill = "grey93") +
  geom_line(data = traces_long, aes(hours, cgroup_gib, group = start, color = role_label),
            linewidth = 0.3, alpha = 0.45) +
  geom_line(data = traces_long, aes(hours, anon_gib, group = start, color = role_label), linewidth = 0.5) +
  geom_text(data = bands[show_label == TRUE],
            aes(x = (band_start + band_end) / 7200, y = Inf, label = label, vjust = label_vjust),
            size = 2.5, color = CLINEGO_COL$secondary) +
  facet_wrap(~role_label, ncol = 1, scales = "free_x") +
  scale_color_manual(values = project_colors, guide = "none") +
  scale_y_continuous(limits = c(0, memory_max * 1.12), expand = c(0, 0)) +
  labs(x = "Hours of run time, the project's runs end to end (idle time between runs left out)",
       y = "Memory, GiB (thick: anonymous; thin: cgroup)") +
  theme_clinego_grid()

caption <- paste(
  "Modes top to bottom in pipeline order; within a mode, projects in order of filtered matrix size.",
  "b: bar = peak anonymous memory, circle = peak cgroup charge with page cache, text = the script holding most",
  "memory in the poll of the peak (process count, summed RSS in GiB).",
  "Memory is sampled every 15 s, so a peak is a lower bound and a run under 15 s may have no sample.",
  "c: bands alternate by run; runs under 3 % of a project's time carry no label.",
  sep = "\n")
figure <- ((plot_time | plot_memory) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) /
  plot_trace + plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", caption = caption,
                  theme = theme(plot.caption = element_text(hjust = 0, size = 9, color = CLINEGO_COL$secondary)))
ggsave(file.path(out_dir, "time_memory.png"), figure, width = 14, height = 13, dpi = 200)

message(sprintf("collect_resources.R: %d real runs, %d skipped rows, %d of %d planned modes measured; tables and figure in %s",
                nrow(kept), nrow(skipped), coverage[status != "not measured", .N], nrow(coverage), out_dir))
