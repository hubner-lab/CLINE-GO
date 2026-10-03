# region_distance.R — LD-decay based snp_clumping_distance resolver
# Usage: source("/pipeline/scripts/R/lib/region_distance.R")
# Requires: data.table

# Hill-Weir (1988) expected r^2 at distance d given parameter C and sample size n
hill_weir <- function(d, C, n) {
    x <- C * d
    ((10 + x) / ((2 + x) * (11 + x))) *
    (1 + ((3 + x) * (12 + 12 * x + x^2)) / (n * (2 + x) * (11 + x)))
}

# Invert Hill-Weir: find distance d where r^2 = r2_target
invert_hill_weir <- function(C_hat, n_samples, r2_target, max_dist_bp) {
    if (length(C_hat) == 0 || length(n_samples) == 0 || length(max_dist_bp) == 0 ||
        is.na(C_hat) || is.na(n_samples) || is.na(max_dist_bp)) return(NA_real_)
    tryCatch(
        uniroot(function(d) hill_weir(d, C_hat, n_samples) - r2_target,
                c(1, max_dist_bp))$root,
        error = function(e) {
            message(paste0("WARNING: Hill-Weir inversion failed (r2=", r2_target,
                           "): ", e$message))
            NA_real_
        }
    )
}

# Fall back to nearest precomputed column when Hill-Weir fit is unavailable
loess_fallback <- function(row, r2_target) {
    half_val   <- as.numeric(row[["half_decay_bp"]])
    r2_02_val  <- as.numeric(row[["r2_02_bp"]])
    r2_intercept <- as.numeric(row[["r2_intercept"]])
    half_r2    <- if (!is.na(r2_intercept)) r2_intercept / 2 else 0.3

    if (!is.na(half_val) && !is.na(r2_02_val)) {
        if (abs(r2_target - half_r2) < abs(r2_target - 0.2)) {
            message(paste0("WARNING: Using half_decay_bp (r2=", round(half_r2, 3),
                           ") as fallback for r2=", r2_target))
            return(half_val)
        } else {
            message(paste0("WARNING: Using r2_02_bp as fallback for r2=", r2_target))
            return(r2_02_val)
        }
    } else if (!is.na(r2_02_val)) {
        return(r2_02_val)
    } else if (!is.na(half_val)) {
        return(half_val)
    }
    return(NA_real_)
}

# Resolve distance for one LD-decay TSV row (as named list / data.table row)
resolve_row_distance <- function(row, r2_target) {
    method <- row[["method"]]
    if (!is.null(method) && method == "hill_weir") {
        C_hat     <- as.numeric(row[["C_hat"]])
        n_samples <- as.numeric(row[["n_samples"]])
        max_dist  <- as.numeric(row[["max_dist_bp"]])

        if (length(C_hat) == 0 || length(max_dist) == 0 || is.na(C_hat) || is.na(max_dist)) {
            message("WARNING: C_hat/max_dist_bp missing — using precomputed fallback. Re-run mode=structure.")
            return(loess_fallback(row, r2_target))
        }
        dist <- invert_hill_weir(C_hat, n_samples, r2_target, max_dist)
        if (!is.na(dist)) return(round(dist))
        return(loess_fallback(row, r2_target))
    }
    return(loess_fallback(row, r2_target))
}

# Resolve snp_clumping_distance from config spec + optional LD decay table.
#
# @param spec         character: "auto_per_chromosome" | "auto_genome_wide" | "auto" | integer string
# @param ld_decay_path  path to LD decay TSV (or "NULL" if not available)
# @param r2_threshold numeric: r2 cutoff for inversion (default 0.2)
# @param group        character: group name in LD decay table (default "All")
# @return either:
#   - scalar integer          (genome-wide or fixed)
#   - named integer vector    (per-chromosome, names = chr names)
resolve_clumping_distance <- function(spec, ld_decay_path = "NULL",
                                      r2_threshold = 0.2, group = "All") {
    # `group` is also a column of the LD-decay table, so a bare `group == group`
    # inside data.table's i is a self-comparison that matches EVERY row rather
    # than the requested one (finding 240d18). Filter against a differently
    # named local copy instead.
    want_group <- as.character(group)
    spec_lower <- tolower(as.character(spec))
    auto_mode  <- spec_lower %in% c("auto_per_chromosome", "auto_genome_wide", "auto")

    if (!auto_mode) {
        dist <- suppressWarnings(as.numeric(spec))
        if (is.na(dist)) {
            message(paste0("WARNING: snp_clumping_distance='", spec,
                           "' is not numeric — using 1e6 bp fallback"))
            return(1000000L)
        }
        message(paste0("INFO: Fixed snp_clumping_distance: ", dist, " bp"))
        return(as.integer(dist))
    }

    if (ld_decay_path == "NULL" || !file.exists(ld_decay_path)) {
        stop("snp_clumping_distance='", spec, "' but LD decay table not found: ", ld_decay_path)
    }
    ld_table <- data.table::fread(ld_decay_path)

    if (nrow(ld_table) == 0) {
        message(paste0("WARNING: LD decay table is empty (", ld_decay_path,
                       ") — using 1e6 bp fallback for '", spec, "'"))
        return(1000000L)
    }

    if (spec_lower %in% c("auto_genome_wide", "auto")) {
        if (spec_lower == "auto") {
            message("INFO: 'auto' treated as 'auto_genome_wide' (deprecated alias)")
        }
        row <- ld_table[group == want_group & scope == "genome_wide"]
        if (nrow(row) == 0) stop("No genome-wide LD decay row for group='", group, "'")
        dist_val <- resolve_row_distance(as.list(row[1]), r2_threshold)
        if (is.na(dist_val)) stop("Could not resolve genome-wide distance from LD decay table")
        message(paste0("INFO: AUTO genome-wide distance (r2<", r2_threshold, "): ",
                       round(dist_val), " bp"))
        return(as.integer(round(dist_val)))
    }

    # auto_per_chromosome
    chr_rows <- ld_table[group == want_group & scope != "genome_wide"]
    gw_row   <- ld_table[group == want_group & scope == "genome_wide"]

    if (nrow(gw_row) == 0 && nrow(chr_rows) == 0)
        stop("No LD decay rows for group='", group, "'")

    gw_dist <- if (nrow(gw_row) > 0) {
        d <- resolve_row_distance(as.list(gw_row[1]), r2_threshold)
        if (!is.na(d)) as.integer(round(d)) else NULL
    } else NULL

    dist_map <- setNames(integer(0), character(0))
    for (i in seq_len(nrow(chr_rows))) {
        row      <- as.list(chr_rows[i])
        chr_name <- row[["scope"]]
        d        <- resolve_row_distance(row, r2_threshold)
        if (is.na(d)) {
            if (!is.null(gw_dist)) {
                message(paste0("WARNING: chr ", chr_name, " inversion failed; using gw fallback: ", gw_dist, " bp"))
                d <- gw_dist
            } else {
                stop("r2 inversion failed for chr ", chr_name, " and no genome-wide fallback available")
            }
        }
        dist_map[chr_name] <- as.integer(round(d))
    }

    message(paste0("INFO: AUTO per-chromosome distances (r2<", r2_threshold, "):"))
    for (chr in names(dist_map)) message(paste0("  chr ", chr, ": ", dist_map[chr], " bp"))
    if (!is.null(gw_dist)) message(paste0("  genome-wide fallback: ", gw_dist, " bp"))

    dist_map
}
