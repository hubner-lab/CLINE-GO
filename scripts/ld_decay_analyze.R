#!/usr/bin/env Rscript
# LD-decay distances per group and scope from the per-distance r^2 sums of
# ld_decay_pairs.R (PLINK --r2, PopLDdecay's .stat.gz layout), against the background
# of ld_decay_background.R (the same estimator between chromosomes), and the plots.
#
# For each group (manifest.tsv) and scope (genome-wide, or one chromosome) the mean
# r^2 is binned (10 bp below 100 bp, log-spaced beyond), smoothed by LOESS, and three
# distances are read off the curve (scripts/R/lib/ld_decay.R):
#   decay_critical_bp -- where it meets the critical r^2 of unlinked markers
#                        (Breseghello & Sorrells 2006; Laido et al. 2014);
#   decay_half_bp     -- where it falls to half its short-range value
#                        (Zhao et al. 2011; Vos et al. 2017);
#   decay_r2thr_bp    -- where it falls to the fixed r^2 R2_THRESHOLD.
# NA means the curve never reaches that level within the measured distance.
#
# .stat.gz layout (tab-separated, with header):
#   #Dist  Mean_r^2  Mean_D'  Sum_r^2  Sum_D'  NumberPairs

library(data.table)
library(ggplot2)
library(svglite)

args <- commandArgs(trailingOnly = TRUE)
################
SCOPE            <- args[1]   # "genome_wide", "per_chromosome", or "both"
STAT_GW_DIR      <- args[2]   # Directory with genome-wide .stat.gz files
STAT_CHR_DIR     <- args[3]   # Directory with per-chromosome .stat.gz files
MANIFEST_PATH    <- args[4]   # manifest.tsv (group, n_samples, list_file)
CHR_LIST_PATH    <- args[5]   # chromosomes.txt
OUT_TABLE        <- args[6]   # ld_decay_half_distances.tsv
PLOT_GW_PATH     <- args[7]   # genome-wide plot path (or "NULL")
PLOT_CHR_PATH    <- args[8]   # per-chromosome plot path (or "NULL")
BACKGROUND_TABLE <- args[9]   # ld_decay_background.tsv (or "NULL")
R2_THRESHOLD     <- as.numeric(args[10])
OUT_CURVE        <- args[11]  # ld_decay_curve.tsv (or "NULL")
LOESS_SPAN       <- as.numeric(args[12])
################
stopifnot(length(args) == 12)
source("/pipeline/scripts/R/lib/ld_decay.R")
source("/pipeline/scripts/R/utils/theme_clinego.R")
tmp <- ld_scratch_dir(STAT_GW_DIR, "analyze_tmp")   # fread(cmd = "zcat ...")'s copies

message("INFO: Analyzing LD decay results")
message(paste0("INFO: Scope: ", SCOPE))

#=============================================================================
# FUNCTIONS
#=============================================================================

read_stat_gz <- function(filepath) {
    if (!file.exists(filepath)) {
        message(paste0("WARNING: File not found: ", filepath))
        return(NULL)
    }
    dt <- tryCatch(
        fread(cmd = paste("zcat", shQuote(filepath)), header = TRUE, tmpdir = tmp),
        error = function(e) {
            message(paste0("WARNING: Failed to read ", filepath, ": ", e$message))
            return(NULL)
        }
    )
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    setnames(dt, colnames(dt), c("dist_bp", "mean_r2", "mean_dprime", "sum_r2", "sum_dprime", "n_pairs"))
    dt[dist_bp > 0]
}

# Placeholder PNG+SVG when no group has data
.write_placeholder_plot <- function(path, msg) {
    if (is.null(path) || path == "NULL" || path == "") return(invisible(NULL))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, color = "grey50") +
        theme_void()
    suppressMessages(ggsave(path, p, width = 10, height = 7, dpi = 150))
    suppressMessages(ggsave(sub("\\.png$", ".svg", path), p, width = 10, height = 7,
                            device = svglite))
}

.kb <- function(bp) {
    ifelse(is.na(bp), "not reached",
           ifelse(bp >= 1e6, paste0(signif(bp / 1e6, 3), " Mb"),
                  ifelse(bp >= 1e3, paste0(signif(bp / 1e3, 3), " kb"),
                         paste0(round(bp), " bp"))))
}
.dist_labels <- function(x) .kb(x)

#=============================================================================
# READ MANIFEST AND BACKGROUND
#=============================================================================
manifest <- fread(MANIFEST_PATH, colClasses = "character")
manifest[, n_samples := as.integer(n_samples)]
if (!"list_file" %in% names(manifest)) manifest[, list_file := group]
message(paste0("INFO: Groups to analyze: ", paste(manifest$group, collapse = ", ")))

bg <- if (BACKGROUND_TABLE != "NULL" && file.exists(BACKGROUND_TABLE)) {
    fread(BACKGROUND_TABLE, colClasses = c(group = "character"))
} else data.table(group = character(), critical_r2 = numeric(), mean_r2 = numeric())
if (!"mean_r2" %in% names(bg)) bg[, mean_r2 := NA_real_]
bg_of <- function(grp, col) {
    v <- bg[group == grp][[col]]
    if (length(v) == 0) NA_real_ else v[1]
}

chr_list <- character(0)
if (SCOPE %in% c("per_chromosome", "both") && file.exists(CHR_LIST_PATH)) {
    chr_list <- readLines(CHR_LIST_PATH)
    chr_list <- chr_list[chr_list != ""]
    message(paste0("INFO: Chromosomes: ", paste(chr_list, collapse = ", ")))
}

#=============================================================================
# DISTANCES PER GROUP AND SCOPE
#=============================================================================
results <- list(); curve_rows <- list(); point_rows <- list()

analyze_one <- function(raw, grp, n_samp, scope) {
    if (is.null(raw)) {
        message(paste0("WARNING: No data for group '", grp, "' (", scope, ")"))
        return(NULL)
    }
    binned <- ld_bin_stat(raw, ld_log_breaks(max(raw$dist_bp)))
    dd <- ld_decay_distances(binned, background_r2 = bg_of(grp, "critical_r2"),
                             r2_threshold = R2_THRESHOLD, span = LOESS_SPAN)
    if (!is.null(binned)) {
        point_rows[[length(point_rows) + 1]] <<- binned[, .(group = grp, scope = scope,
                                                            dist_bp, mean_r2)]
    }
    if (!is.null(dd$curve)) {
        curve_rows[[length(curve_rows) + 1]] <<- dd$curve[, .(group = grp, scope = scope,
            dist_bp = round(dist_bp, 1), r2 = round(r2, 6), r2_mono = round(r2_mono, 6))]
    }
    data.table(group = grp, scope = scope, n_samples = n_samp, n_pairs = sum(raw$n_pairs),
               max_dist_bp = max(raw$dist_bp), short_r2 = round(dd$short_r2, 4),
               background_r2 = dd$background_r2,
               background_mean_r2 = bg_of(grp, "mean_r2"),
               decay_critical_bp = round(dd$critical_bp), decay_half_bp = round(dd$half_bp),
               decay_r2thr_bp = round(dd$r2thr_bp), r2_threshold = R2_THRESHOLD,
               decay_method = dd$method)
}

if (SCOPE %in% c("genome_wide", "both")) {
    message("INFO: Processing genome-wide results")
    for (i in seq_len(nrow(manifest))) {
        raw <- read_stat_gz(file.path(STAT_GW_DIR, paste0(manifest$list_file[i], ".stat.gz")))
        row <- analyze_one(raw, manifest$group[i], manifest$n_samples[i], "genome_wide")
        if (!is.null(row)) {
            results[[length(results) + 1]] <- row
            message(paste0("INFO: ", row$group, " genome-wide: critical ", .kb(row$decay_critical_bp),
                           ", half-decay ", .kb(row$decay_half_bp), ", r2 = ", R2_THRESHOLD,
                           " at ", .kb(row$decay_r2thr_bp)))
        }
    }
}

if (SCOPE %in% c("per_chromosome", "both") && length(chr_list) > 0) {
    message("INFO: Processing per-chromosome results")
    for (i in seq_len(nrow(manifest))) {
        for (chr in chr_list) {
            raw <- read_stat_gz(file.path(STAT_CHR_DIR,
                                          paste0(manifest$list_file[i], "_", chr, ".stat.gz")))
            row <- analyze_one(raw, manifest$group[i], manifest$n_samples[i], chr)
            if (!is.null(row)) results[[length(results) + 1]] <- row
        }
    }
}

#=============================================================================
# WRITE OUTPUT TABLES
#=============================================================================
result_table <- if (length(results) > 0) rbindlist(results) else
    data.table(group = character(), scope = character(), n_samples = integer(),
               n_pairs = numeric(), max_dist_bp = numeric(), short_r2 = numeric(),
               background_r2 = numeric(), background_mean_r2 = numeric(),
               decay_critical_bp = numeric(), decay_half_bp = numeric(),
               decay_r2thr_bp = numeric(), r2_threshold = numeric(),
               decay_method = character())
fwrite(result_table, OUT_TABLE, sep = "\t")
message(paste0("INFO: Written results table: ", OUT_TABLE, " (", nrow(result_table), " rows)"))

if (OUT_CURVE != "NULL") {
    curve_table <- if (length(curve_rows) > 0) rbindlist(curve_rows) else
        data.table(group = character(), scope = character(), dist_bp = numeric(),
                   r2 = numeric(), r2_mono = numeric())
    fwrite(curve_table, OUT_CURVE, sep = "\t")
    message(paste0("INFO: Written LD-decay curves: ", OUT_CURVE, " (", nrow(curve_table), " rows)"))
}

#=============================================================================
# PLOTS
#=============================================================================
# Binned means (points: bins, not SNP pairs), the LOESS curve per group, and for the
# 'All' group the critical r^2 (horizontal) and the critical and half-decay distances
# (vertical), the lines the clumping distance is read from.
points_all <- if (length(point_rows)) rbindlist(point_rows) else NULL
curves_all <- if (length(curve_rows)) rbindlist(curve_rows) else NULL

group_palette <- function(groups) {
    others <- setdiff(groups, "All")
    cols <- c(All = CLINEGO_COL$fg, setNames(clinego_cluster_palette(length(others)), others))
    cols[groups]
}

decay_plot <- function(pts, crv, res, facet = FALSE) {
    groups <- c("All", setdiff(unique(res$group), "All"))
    groups <- groups[groups %in% res$group]
    lab <- res[, .(n = n_samples[1]), by = group]
    lab_map <- setNames(paste0(lab$group, " (n = ", lab$n, ")"), lab$group)
    cols <- group_palette(groups); names(cols) <- lab_map[groups]
    widths <- setNames(ifelse(groups == "All", 1.1, 0.5), lab_map[groups])
    pts <- pts[group %in% groups][, group_label := lab_map[group]]
    crv <- crv[group %in% groups][, group_label := lab_map[group]]
    all_res <- res[group == "All"]
    vlines <- rbind(all_res[!is.na(decay_critical_bp), .(scope, dist_bp = decay_critical_bp,
                                                          rule = "critical r²")],
                    all_res[!is.na(decay_half_bp), .(scope, dist_bp = decay_half_bp,
                                                      rule = "half-decay")])
    hline <- all_res[!is.na(background_r2), .(scope, r2 = background_r2)]
    p <- ggplot() +
        geom_point(data = pts, aes(dist_bp, mean_r2, color = group_label), size = 0.6, alpha = 0.5) +
        geom_line(data = crv, aes(dist_bp, r2, color = group_label, linewidth = group_label)) +
        scale_color_manual(values = cols, name = "Group") +
        scale_linewidth_manual(values = widths, guide = "none") +
        scale_x_log10(labels = .dist_labels) +
        labs(x = "Distance between SNPs", y = expression(Mean~r^2)) +
        theme_clinego() +
        theme(legend.direction = "vertical")
    if (nrow(hline)) {
        p <- p + geom_hline(data = hline, aes(yintercept = r2), color = CLINEGO_THRESHOLD,
                            linetype = "dashed", linewidth = 0.4)
    }
    if (nrow(vlines)) {
        p <- p + geom_vline(data = vlines, aes(xintercept = dist_bp, linetype = rule),
                            color = CLINEGO_THRESHOLD, linewidth = 0.5) +
            scale_linetype_manual(values = c("critical r²" = "solid", "half-decay" = "dotted"),
                                  name = "'All' decay distance")
    }
    if (facet) p <- p + facet_wrap(~ scope, scales = "free_x")
    p
}

save_plot <- function(p, path, width, height) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    ggsave(path, p, width = width, height = height, dpi = 300)
    ggsave(sub("\\.png$", ".svg", path), p, width = width, height = height, device = svglite)
    message(paste0("INFO: Saved plot: ", path))
}

gw_res <- result_table[scope == "genome_wide"]
if (PLOT_GW_PATH != "NULL") {
    if (nrow(gw_res) == 0 || is.null(curves_all) || !nrow(curves_all[scope == "genome_wide"])) {
        message("WARNING: No genome-wide LD decay data — writing placeholder plot")
        .write_placeholder_plot(PLOT_GW_PATH, "No LD decay data (insufficient group sizes)")
    } else {
        p <- decay_plot(points_all[scope == "genome_wide"], curves_all[scope == "genome_wide"], gw_res) +
            labs(title = "LD decay (genome-wide)")
        save_plot(p, PLOT_GW_PATH, 10, 7)
    }
}

chr_res <- result_table[scope != "genome_wide"]
if (PLOT_CHR_PATH != "NULL") {
    if (nrow(chr_res) == 0 || is.null(curves_all) || !nrow(curves_all[scope != "genome_wide"])) {
        message("WARNING: No per-chromosome LD decay data — writing placeholder plot")
        .write_placeholder_plot(PLOT_CHR_PATH, "No LD decay data (insufficient group sizes)")
    } else {
        p <- decay_plot(points_all[scope != "genome_wide"], curves_all[scope != "genome_wide"],
                        chr_res, facet = TRUE) +
            labs(title = "LD decay (per chromosome)")
        n_chr <- uniqueN(chr_res$scope)
        save_plot(p, PLOT_CHR_PATH, 12, 3 + ceiling(n_chr / min(3, n_chr)) * 3)
    }
}

unlink(tmp, recursive = TRUE)
message("INFO: LD decay analysis complete")
