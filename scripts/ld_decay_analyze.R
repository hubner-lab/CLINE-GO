#!/usr/bin/env Rscript
# Analyze PopLDdecay output: read .stat.gz files, fit Hill-Weir curves,
# extract half-decay distances, generate ggplot2 plots
#
# PopLDdecay .stat.gz format (tab-separated, no header):
#   Dist(bp)  Mean_r^2  NumberPairs  Sum_r^2
# One row per bp distance from 0 to MaxDist*1000

library(data.table)
library(ggplot2)
library(svglite)

args <- commandArgs(trailingOnly = TRUE)
################
SCOPE          <- args[1]  # "genome_wide", "per_chromosome", or "both"
STAT_GW_DIR    <- args[2]  # Directory with genome-wide .stat.gz files
STAT_CHR_DIR   <- args[3]  # Directory with per-chromosome .stat.gz files
MANIFEST_PATH  <- args[4]  # manifest.tsv (group, n_samples)
CHR_LIST_PATH  <- args[5]  # chromosomes.txt
OUT_TABLE      <- args[6]  # ld_decay_half_distances.tsv
PLOT_GW_PATH   <- args[7]  # genome-wide plot path (or "NULL")
PLOT_CHR_PATH  <- args[8]  # per-chromosome plot path (or "NULL")
################

message("INFO: Analyzing LD decay results")
message(paste0("INFO: Scope: ", SCOPE))

#=============================================================================
# FUNCTIONS
#=============================================================================

read_stat_gz <- function(filepath) {
    # PopLDdecay .stat.gz format (with header):
    # #Dist  Mean_r^2  Mean_D'  Sum_r^2  Sum_D'  NumberPairs
    if (!file.exists(filepath)) {
        message(paste0("WARNING: File not found: ", filepath))
        return(NULL)
    }
    dt <- tryCatch(
        fread(cmd = paste0("zcat ", filepath), header = TRUE),
        error = function(e) {
            message(paste0("WARNING: Failed to read ", filepath, ": ", e$message))
            return(NULL)
        }
    )
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    # Standardize column names (PopLDdecay uses #Dist, Mean_r^2, etc.)
    cols <- colnames(dt)
    setnames(dt, cols, c("dist_bp", "mean_r2", "mean_dprime", "sum_r2", "sum_dprime", "n_pairs"))
    # Remove distance 0 (self-comparisons)
    dt <- dt[dist_bp > 0]
    dt
}

bin_stat <- function(dt, n_bins = 100) {
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    max_dist <- max(dt$dist_bp)
    bin_size <- max(1000, ceiling(max_dist / n_bins))  # minimum 1kb bins
    dt[, bin := floor(dist_bp / bin_size) * bin_size + bin_size / 2]
    binned <- dt[, .(mean_r2 = sum(sum_r2) / sum(n_pairs),
                     n_pairs = sum(n_pairs),
                     dist_kb = bin[1] / 1000),
                 by = bin]
    binned[order(bin)]
}

hill_weir <- function(d, C, n) {
    # Hill & Weir (1988) expectation of r^2
    x <- C * d
    ((10 + x) / ((2 + x) * (11 + x))) *
    (1 + ((3 + x) * (12 + 12 * x + x^2)) / (n * (2 + x) * (11 + x)))
}

fit_ld_curve <- function(binned, n_samples) {
    # Fit Hill-Weir curve, fallback to LOESS
    if (is.null(binned) || nrow(binned) < 5) {
        return(list(method = "insufficient_data", fit_fn = NULL,
                    r2_intercept = NA, half_decay_bp = NA, r2_02_bp = NA,
                    C_hat = NA, max_dist_bp = NA))
    }

    r2_intercept <- binned[which.min(bin), mean_r2]
    max_dist_bp <- max(binned$bin)

    # Hill-Weir carries a sample-size correction term 1/(n*(2+x)*(11+x)), so an
    # absent or degenerate n makes the model meaningless: n = 0 sends the term to
    # Inf and nls always fails, which used to surface only as an INFO line and a
    # silent LOESS fallback (finding 600acb). Say so, and skip straight to LOESS.
    n_usable <- !is.na(n_samples) && n_samples >= 2
    if (!n_usable) {
        message(paste0("WARNING: n_samples=", n_samples, " is not usable in the ",
                       "Hill-Weir sample-size term (need >= 2) — fitting LOESS ",
                       "instead. Check the LD-decay manifest."))
    }

    # Try Hill-Weir nls fit
    fit_result <- if (!n_usable) NULL else tryCatch({
        fit <- nls(mean_r2 ~ hill_weir(bin, C, n_samples),
                   data = binned,
                   start = list(C = 0.001),
                   control = nls.control(maxiter = 500, warnOnly = TRUE))
        C_hat <- coef(fit)[["C"]]

        # Create prediction function
        fit_fn <- function(d) hill_weir(d, C_hat, n_samples)

        # Half-decay distance: where curve crosses r2_intercept / 2
        half_target <- r2_intercept / 2
        half_decay_bp <- tryCatch(
            uniroot(function(d) fit_fn(d) - half_target, c(1, max_dist_bp))$root,
            error = function(e) NA
        )

        # r2 = 0.2 threshold distance
        r2_02_bp <- tryCatch(
            uniroot(function(d) fit_fn(d) - 0.2, c(1, max_dist_bp))$root,
            error = function(e) NA
        )

        list(method = "hill_weir", fit_fn = fit_fn,
             r2_intercept = r2_intercept, half_decay_bp = half_decay_bp,
             r2_02_bp = r2_02_bp, C_hat = C_hat, max_dist_bp = max_dist_bp)
    }, error = function(e) {
        message(paste0("INFO: Hill-Weir fit failed, using LOESS: ", e$message))
        NULL
    })

    if (!is.null(fit_result)) return(fit_result)

    # Fallback: LOESS
    tryCatch({
        lo <- loess(mean_r2 ~ bin, data = binned, span = 0.3)
        fit_fn <- function(d) predict(lo, newdata = data.frame(bin = d))

        half_target <- r2_intercept / 2
        pred_grid <- data.frame(bin = seq(1, max_dist_bp, length.out = 1000))
        pred_grid$r2 <- predict(lo, newdata = pred_grid)
        pred_grid <- pred_grid[!is.na(pred_grid$r2), ]

        half_decay_bp <- if (any(pred_grid$r2 <= half_target)) {
            pred_grid$bin[which(pred_grid$r2 <= half_target)[1]]
        } else NA

        r2_02_bp <- if (any(pred_grid$r2 <= 0.2)) {
            pred_grid$bin[which(pred_grid$r2 <= 0.2)[1]]
        } else NA

        list(method = "loess", fit_fn = fit_fn,
             r2_intercept = r2_intercept, half_decay_bp = half_decay_bp,
             r2_02_bp = r2_02_bp, C_hat = NA, max_dist_bp = max_dist_bp)
    }, error = function(e) {
        message(paste0("WARNING: LOESS fallback also failed: ", e$message))
        list(method = "failed", fit_fn = NULL,
             r2_intercept = r2_intercept, half_decay_bp = NA, r2_02_bp = NA,
             C_hat = NA, max_dist_bp = max_dist_bp)
    })
}

#=============================================================================
# READ MANIFEST
#=============================================================================
manifest <- fread(MANIFEST_PATH)
message(paste0("INFO: Groups to analyze: ", paste(manifest$group, collapse = ", ")))

# Read chromosome list if needed
chr_list <- character(0)
if (SCOPE %in% c("per_chromosome", "both") && file.exists(CHR_LIST_PATH)) {
    chr_list <- readLines(CHR_LIST_PATH)
    chr_list <- chr_list[chr_list != ""]
    message(paste0("INFO: Chromosomes: ", paste(chr_list, collapse = ", ")))
}

#=============================================================================
# PROCESS GENOME-WIDE RESULTS
#=============================================================================
results <- list()
gw_binned_all <- list()  # For plotting

if (SCOPE %in% c("genome_wide", "both")) {
    message("INFO: Processing genome-wide results")
    for (i in seq_len(nrow(manifest))) {
        grp <- manifest$group[i]
        n_samp <- manifest$n_samples[i]
        stat_file <- file.path(STAT_GW_DIR, paste0(grp, ".stat.gz"))
        raw <- read_stat_gz(stat_file)
        binned <- bin_stat(raw)

        if (is.null(binned)) {
            message(paste0("WARNING: No data for group '", grp, "' genome-wide"))
            next
        }

        fit <- fit_ld_curve(binned, n_samp)
        binned[, group := grp]
        gw_binned_all[[grp]] <- binned

        results[[length(results) + 1]] <- data.table(
            group = grp, scope = "genome_wide",
            half_decay_bp = round(fit$half_decay_bp),
            r2_02_bp = round(fit$r2_02_bp),
            n_samples = n_samp,
            n_pairs = sum(binned$n_pairs),
            r2_intercept = round(fit$r2_intercept, 4),
            method = fit$method,
            C_hat = fit$C_hat,
            max_dist_bp = fit$max_dist_bp
        )
        message(paste0("INFO: ", grp, " genome-wide: half-decay=",
                        round(fit$half_decay_bp / 1000, 1), "kb, method=", fit$method))
    }
}

#=============================================================================
# PROCESS PER-CHROMOSOME RESULTS
#=============================================================================
chr_binned_all <- list()  # For plotting

if (SCOPE %in% c("per_chromosome", "both") && length(chr_list) > 0) {
    message("INFO: Processing per-chromosome results")
    for (i in seq_len(nrow(manifest))) {
        grp <- manifest$group[i]
        n_samp <- manifest$n_samples[i]
        for (chr in chr_list) {
            stat_file <- file.path(STAT_CHR_DIR, paste0(grp, "_", chr, ".stat.gz"))
            raw <- read_stat_gz(stat_file)
            binned <- bin_stat(raw)

            if (is.null(binned)) {
                message(paste0("WARNING: No data for group '", grp, "' chr ", chr))
                next
            }

            fit <- fit_ld_curve(binned, n_samp)
            binned[, `:=`(group = grp, chr = chr)]
            chr_binned_all[[paste0(grp, "_", chr)]] <- binned

            results[[length(results) + 1]] <- data.table(
                group = grp, scope = chr,
                half_decay_bp = round(fit$half_decay_bp),
                r2_02_bp = round(fit$r2_02_bp),
                n_samples = n_samp,
                n_pairs = sum(binned$n_pairs),
                r2_intercept = round(fit$r2_intercept, 4),
                method = fit$method,
                C_hat = fit$C_hat,
                max_dist_bp = fit$max_dist_bp
            )
        }
    }
}

#=============================================================================
# WRITE OUTPUT TABLE
#=============================================================================
if (length(results) > 0) {
    result_table <- rbindlist(results)
} else {
    result_table <- data.table(group = character(), scope = character(),
                               half_decay_bp = numeric(), r2_02_bp = numeric(),
                               n_samples = integer(), n_pairs = integer(),
                               r2_intercept = numeric(), method = character())
}
fwrite(result_table, OUT_TABLE, sep = "\t")
message(paste0("INFO: Written results table: ", OUT_TABLE, " (", nrow(result_table), " rows)"))

# Helper: create a blank placeholder PNG+SVG when no data is available
.write_placeholder_plot <- function(path, msg) {
    if (is.null(path) || path == "NULL" || path == "") return(invisible(NULL))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, color = "grey50") +
        theme_void()
    suppressMessages(ggsave(path, p, width = 10, height = 7, dpi = 150))
    svg_path <- sub("\\.png$", ".svg", path)
    suppressMessages(ggsave(svg_path, p, width = 10, height = 7, device = svglite))
}

#=============================================================================
# GENOME-WIDE PLOT
#=============================================================================
if (PLOT_GW_PATH != "NULL" && length(gw_binned_all) == 0) {
    message("WARNING: No genome-wide LD decay data — writing placeholder plot")
    .write_placeholder_plot(PLOT_GW_PATH, "No LD decay data (insufficient group sizes)")
}
if (PLOT_GW_PATH != "NULL" && length(gw_binned_all) > 0) {
    dir.create(dirname(PLOT_GW_PATH), recursive = TRUE, showWarnings = FALSE)
    message("INFO: Generating genome-wide LD decay plot")
    gw_data <- rbindlist(gw_binned_all)

    gw_results <- result_table[scope == "genome_wide"]

    # Build legend labels (r2<0.2 only — half-decay kept in table only)
    gw_results[, label := paste0(group, " (n=", n_samples,
                                  ", r2<0.2=", round(r2_02_bp / 1000, 0), "kb)")]
    label_map <- setNames(gw_results$label, gw_results$group)
    gw_data[, group_label := label_map[group]]

    # Color scheme: "All" in black, others in distinct colors
    groups_ordered <- c("All", setdiff(unique(gw_data$group), "All"))
    n_groups <- length(groups_ordered)
    if (n_groups > 1) {
        colors <- c("black", scales::hue_pal()(n_groups - 1))
    } else {
        colors <- "black"
    }
    names(colors) <- sapply(groups_ordered, function(g) label_map[g])

    # Size: "All" thicker
    sizes <- c(1.2, rep(0.8, n_groups - 1))
    names(sizes) <- sapply(groups_ordered, function(g) label_map[g])

    # Curves only — no raw points (WGS data produces too many)
    p <- ggplot(gw_data, aes(x = dist_kb, y = mean_r2, color = group_label)) +
        geom_smooth(method = "loess", span = 0.3, se = FALSE,
                    aes(linewidth = group_label)) +
        geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey40") +
        scale_color_manual(values = colors) +
        scale_linewidth_manual(values = sizes) +
        labs(x = "Distance (kb)", y = expression(Mean~r^2),
             title = "LD Decay (Genome-wide)",
             color = "Group", linewidth = "Group") +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.direction = "vertical",
              legend.text = element_text(size = 8))

    # Add r2<0.2 vertical lines (prominent — the practical LD extent)
    for (j in seq_len(nrow(gw_results))) {
        r2_02 <- gw_results$r2_02_bp[j]
        if (!is.na(r2_02)) {
            col_idx <- match(gw_results$label[j], names(colors))
            p <- p + geom_vline(xintercept = r2_02 / 1000, linetype = "dashed",
                                color = colors[col_idx], alpha = 0.8)
        }
    }

    ggsave(PLOT_GW_PATH, p, width = 10, height = 7, dpi = 300)
    svg_path <- sub("\\.png$", ".svg", PLOT_GW_PATH)
    ggsave(svg_path, p, width = 10, height = 7, device = svglite)
    message(paste0("INFO: Saved genome-wide plot: ", PLOT_GW_PATH))
}

#=============================================================================
# PER-CHROMOSOME PLOT
#=============================================================================
if (PLOT_CHR_PATH != "NULL" && length(chr_binned_all) == 0) {
    message("WARNING: No per-chromosome LD decay data — writing placeholder plot")
    .write_placeholder_plot(PLOT_CHR_PATH, "No LD decay data (insufficient group sizes)")
}
if (PLOT_CHR_PATH != "NULL" && length(chr_binned_all) > 0) {
    dir.create(dirname(PLOT_CHR_PATH), recursive = TRUE, showWarnings = FALSE)
    message("INFO: Generating per-chromosome LD decay plot")
    chr_data <- rbindlist(chr_binned_all)

    # Color by group with per-chromosome r2<0.2 metrics in legend
    groups_ordered <- c("All", setdiff(unique(chr_data$group), "All"))
    n_groups <- length(groups_ordered)
    if (n_groups > 1) {
        colors <- c("black", scales::hue_pal()(n_groups - 1))
    } else {
        colors <- "black"
    }

    # Build legend labels with per-chromosome r2<0.2 distances
    chr_results <- result_table[scope != "genome_wide"]
    chr_order <- sort(unique(chr_results$scope))
    label_map <- sapply(groups_ordered, function(grp) {
        grp_rows <- chr_results[group == grp]
        if (nrow(grp_rows) == 0) return(grp)
        # Build per-chr r2<0.2 string: "chr1=Xkb, chr2=Ykb, ..."
        chr_metrics <- sapply(chr_order, function(ch) {
            val <- grp_rows[scope == ch, r2_02_bp]
            if (length(val) == 0 || is.na(val[1])) return(paste0(ch, "=NA"))
            paste0(ch, "=", round(val[1] / 1000, 0), "kb")
        })
        paste0(grp, " (r2<0.2: ", paste(chr_metrics, collapse = ", "), ")")
    })
    names(colors) <- label_map
    chr_data[, group_label := label_map[group]]

    # Curves only — no raw points (WGS data produces too many)
    p <- ggplot(chr_data, aes(x = dist_kb, y = mean_r2, color = group_label)) +
        geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 0.8) +
        geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey40") +
        facet_wrap(~ chr, scales = "free_x") +
        scale_color_manual(values = colors) +
        labs(x = "Distance (kb)", y = expression(Mean~r^2),
             title = "LD Decay (Per Chromosome)",
             color = "Group") +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.direction = "vertical",
              legend.text = element_text(size = 7),
              strip.text = element_text(size = 9))

    # Determine plot height based on number of chromosomes
    n_chr <- length(unique(chr_data$chr))
    n_cols <- min(3, n_chr)
    n_rows <- ceiling(n_chr / n_cols)
    plot_height <- 3 + n_rows * 3

    ggsave(PLOT_CHR_PATH, p, width = 12, height = plot_height, dpi = 300)
    svg_path <- sub("\\.png$", ".svg", PLOT_CHR_PATH)
    ggsave(svg_path, p, width = 12, height = plot_height, device = svglite)
    message(paste0("INFO: Saved per-chromosome plot: ", PLOT_CHR_PATH))
}

message("INFO: LD decay analysis complete")
