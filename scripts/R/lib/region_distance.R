# region_distance.R — LD-decay based snp_clumping_distance resolver
# Usage: source("/pipeline/scripts/R/lib/region_distance.R")
# Requires: data.table
#
# An auto_* distance is read from the LD-decay table (ld_decay_analyze.R) under one of
# three rules, each a column of that table:
#   critical     -- decay_critical_bp: the LD-decay curve meets the critical r^2 of
#                   unlinked markers (Breseghello & Sorrells 2006; Laido et al. 2014);
#   half_decay   -- decay_half_bp: r^2 falls to half its short-range value
#                   (Zhao et al. 2011; Vos et al. 2017);
#   r2_threshold -- decay_r2thr_bp: r^2 falls to LDdecay.r2_threshold.

LD_RULE_COLUMNS <- c(critical = "decay_critical_bp", half_decay = "decay_half_bp",
                     r2_threshold = "decay_r2thr_bp")

# The rule name for a rule argument. A number (or numeric string) is the r^2 level of
# the r2_threshold rule from callers that predate the rules; it must match the level
# the table was computed at (`table_threshold`, NA when unknown).
.ld_rule <- function(rule, table_threshold = NA_real_) {
    r <- tolower(trimws(as.character(rule)))
    if (r %in% names(LD_RULE_COLUMNS)) return(r)
    num <- suppressWarnings(as.numeric(r))
    if (is.na(num)) {
        stop("unknown LD-decay rule '", rule, "': expected one of ",
             paste(names(LD_RULE_COLUMNS), collapse = ", "))
    }
    if (!is.na(table_threshold) && abs(num - table_threshold) > 1e-9) {
        stop("r2 threshold ", num, " requested, but the LD-decay table was computed at r2 = ",
             table_threshold, "; set LDdecay.r2_threshold instead")
    }
    "r2_threshold"
}

# The distance of one LD-decay table row under `rule` (NA when the curve never reaches
# the rule's level within the measured distance).
resolve_row_distance <- function(row, rule) {
    as.numeric(row[[LD_RULE_COLUMNS[[rule]]]])
}

# Stop for a rule that gives NA on `row` (a table row as a list), naming why: too few
# bins to fit a curve, no level to decay to, a curve that starts at or below the level
# (a small group whose unlinked background is as high as its short-range LD), or a
# curve that never falls to it within LDdecay.max_distance.
.not_reached <- function(rule, group, where, row) {
    s <- suppressWarnings(as.numeric(row[["short_r2"]]))
    target <- switch(rule, critical = suppressWarnings(as.numeric(row[["background_r2"]])),
                     half_decay = s / 2,
                     r2_threshold = suppressWarnings(as.numeric(row[["r2_threshold"]])))
    why <- if (identical(row[["decay_method"]], "insufficient_data")) {
        "too few distance bins to fit a curve"
    } else if (length(target) == 0 || is.na(target)) {
        "there is no level to decay to (no critical r2 without unlinked pairs)"
    } else if (!is.na(s) && s <= target) {
        paste0("the curve starts at or below that level (short-range r2 ", signif(s, 3),
               ", level ", signif(target, 3), "), so there is no decay to read")
    } else {
        "LD never decays to that level within LDdecay.max_distance; raising it can help"
    }
    stop("LD-decay rule '", rule, "' is not reached for group '", group, "' (", where,
         "): ", why, ". Choose another clumping_ld_rule or set a fixed distance in bp")
}

# Resolve snp_clumping_distance from config spec + optional LD decay table.
#
# @param spec           character: "auto_per_chromosome" | "auto_genome_wide" | "auto" | integer string
# @param ld_decay_path  path to LD decay TSV (or "NULL" if not available)
# @param rule           "critical" | "half_decay" | "r2_threshold" (a number = r2_threshold)
# @param group          character: group name in LD decay table (default "All")
# @return either:
#   - scalar integer          (genome-wide or fixed)
#   - named integer vector    (per-chromosome, names = chr names)
resolve_clumping_distance <- function(spec, ld_decay_path = "NULL",
                                      rule = "critical", group = "All") {
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
    ld_table <- data.table::fread(ld_decay_path, colClasses = c(group = "character",
                                                               scope = "character"))
    if (nrow(ld_table) == 0) {
        stop("snp_clumping_distance='", spec, "' but the LD decay table is empty: ",
             ld_decay_path, " (no group had LD-decay data)")
    }
    missing_cols <- setdiff(LD_RULE_COLUMNS, names(ld_table))
    if (length(missing_cols)) {
        stop("LD decay table ", ld_decay_path, " lacks ", paste(missing_cols, collapse = ", "),
             ": it predates the PLINK LD-decay estimator; re-run mode=structure")
    }
    table_thr <- if ("r2_threshold" %in% names(ld_table)) ld_table$r2_threshold[1] else NA_real_
    rule <- .ld_rule(rule, table_thr)
    if (!want_group %in% ld_table$group) {
        stop("LD-decay group '", group, "' is not in ", ld_decay_path, " (groups: ",
             paste(unique(ld_table$group), collapse = ", "), "): check ld_decay_group, ",
             "LDdecay.group_by and LDdecay.min_samples")
    }

    if (spec_lower %in% c("auto_genome_wide", "auto")) {
        if (spec_lower == "auto") {
            message("INFO: 'auto' treated as 'auto_genome_wide' (deprecated alias)")
        }
        row <- ld_table[group == want_group & scope == "genome_wide"]
        if (nrow(row) == 0) {
            stop("No genome-wide LD decay row for group='", group, "': auto_genome_wide ",
                 "needs LDdecay.scope genome_wide or both")
        }
        dist_val <- resolve_row_distance(as.list(row[1]), rule)
        if (is.na(dist_val)) .not_reached(rule, group, "genome-wide", as.list(row[1]))
        message(paste0("INFO: AUTO genome-wide distance (", rule, "): ", round(dist_val), " bp"))
        return(as.integer(round(dist_val)))
    }

    # auto_per_chromosome
    chr_rows <- ld_table[group == want_group & scope != "genome_wide"]
    gw_row   <- ld_table[group == want_group & scope == "genome_wide"]

    if (nrow(chr_rows) == 0) {
        stop("snp_clumping_distance='auto_per_chromosome' but the LD decay table has no ",
             "per-chromosome rows for group '", group, "': set LDdecay.scope to ",
             "per_chromosome or both, or use auto_genome_wide")
    }

    gw_dist <- if (nrow(gw_row) > 0) {
        d <- resolve_row_distance(as.list(gw_row[1]), rule)
        if (!is.na(d)) as.integer(round(d)) else NULL
    } else NULL

    dist_map <- setNames(integer(0), character(0))
    for (i in seq_len(nrow(chr_rows))) {
        row      <- as.list(chr_rows[i])
        chr_name <- row[["scope"]]
        d        <- resolve_row_distance(row, rule)
        if (is.na(d)) {
            if (!is.null(gw_dist)) {
                message(paste0("WARNING: chr ", chr_name, ": rule '", rule,
                               "' not reached; using the genome-wide distance ", gw_dist, " bp"))
                d <- gw_dist
            } else {
                .not_reached(rule, group, paste0("chromosome ", chr_name,
                                                 ", with no genome-wide value to fall back to"),
                             row)
            }
        }
        dist_map[chr_name] <- as.integer(round(d))
    }

    message(paste0("INFO: AUTO per-chromosome distances (", rule, "):"))
    for (chr in names(dist_map)) message(paste0("  chr ", chr, ": ", dist_map[chr], " bp"))
    if (!is.null(gw_dist)) message(paste0("  genome-wide fallback: ", gw_dist, " bp"))

    dist_map
}

# Strict variant for the association scripts (find_sig_snps.R, combine_selected_snps.R):
# a spec that is neither a number of bp nor an auto_* mode is a configuration error,
# never the silent 1e6 bp fallback resolve_clumping_distance() keeps for old callers.
resolve_clumping_distance_strict <- function(spec, ld_decay_path = "NULL",
                                             rule = "critical", group = "All") {
    spec_lower <- tolower(trimws(as.character(spec)))
    is_auto    <- spec_lower %in% c("auto_per_chromosome", "auto_genome_wide", "auto")
    if (!is_auto && is.na(suppressWarnings(as.numeric(spec)))) {
        stop("clumping distance '", spec, "' is neither a number of bp nor an auto_* mode")
    }
    resolve_clumping_distance(spec, ld_decay_path, rule, group)
}

# One-line description of a resolved distance for the logs.
format_clumping_distance <- function(dist) {
    if (length(dist) == 1 && is.null(names(dist))) return(paste0(dist, " bp"))
    paste0("per chromosome, ", min(dist), "-", max(dist), " bp (",
           paste0(names(dist), "=", dist, collapse = ", "), ")")
}
