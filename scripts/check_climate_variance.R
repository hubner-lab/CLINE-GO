#!/usr/bin/env Rscript
# Detect invariant factors at sample sites: all values NA (all_na), one non-NA
# value (single_value) or zero variance (zero_variance).
# Writes a TSV with columns: predictor, reason
# Header-only file when all factors are valid.
#
# Serves two modules (Phase 3):
#   COLUMNS=traits         phenotypic traits — every numeric column that is not
#                          sample/site/coordinate metadata (mode=traits), the same
#                          exclusion rule as plot_density.R:27-31
#   COLUMNS=bio_1,bio_2,…  climate predictors, named: the predictor universe
#                          (ALL_BIO: the 19 WorldClim layers, or the custom columns),
#                          not only Climate.predictors, because the app strikes every
#                          invariant chip from this table (mode=climate). Each must be
#                          a column of SITE_FILE.

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
SITE_FILE  <- args[1]
OUT_FILE   <- args[2]
if (length(args) < 3 || !nzchar(args[3]))
    stop("usage: check_climate_variance.R SITE_FILE OUT_FILE traits|<predictor,predictor,...>")
COLUMNS    <- args[3]

site <- fread(SITE_FILE, sep = '\t', header = TRUE)

meta_cols <- c("sample", "site", "latitude", "longitude", "lat", "lon", "ID")

if (identical(COLUMNS, "traits")) {
    cand <- colnames(site)[!tolower(colnames(site)) %in% tolower(meta_cols)]
    target_cols <- cand[vapply(cand, function(x) is.numeric(site[[x]]), logical(1))]
    label <- "traits"
} else {
    wanted  <- trimws(strsplit(COLUMNS, ",", fixed = TRUE)[[1]])
    wanted  <- wanted[nzchar(wanted)]
    missing <- setdiff(wanted, colnames(site))
    if (length(wanted) == 0 || length(missing) > 0)
        stop("check_climate_variance.R: predictor(s) not in ", SITE_FILE, ": ",
             paste(if (length(wanted) == 0) "<none named>" else missing, collapse = ", "))
    # A name list, not grep("^bio_"): for custom climate the pattern matched nothing,
    # and an empty table read as "every predictor varies". Table order is kept.
    target_cols <- intersect(colnames(site), wanted)
    label <- "climate predictors"
}

invariant <- data.table(predictor = character(), reason = character())

for (col in target_cols) {
    vals <- site[[col]]
    # Coerce to numeric — empty strings become NA
    vals <- suppressWarnings(as.numeric(as.character(vals)))

    n_obs <- sum(!is.na(vals))
    if (n_obs == 0) {
        invariant <- rbind(invariant, data.table(predictor = col, reason = "all_na"))
        message("WARNING: ", col, " — all values are NA/empty")
    } else if (n_obs == 1) {
        # sd() of one value is NA, and `else if (NA)` used to abort the mode.
        invariant <- rbind(invariant, data.table(predictor = col, reason = "single_value"))
        message("WARNING: ", col, " — one non-NA value only")
    } else if (sd(vals, na.rm = TRUE) == 0) {
        invariant <- rbind(invariant, data.table(predictor = col, reason = "zero_variance"))
        message("WARNING: ", col, " — zero variance across all samples")
    }
}

fwrite(invariant, OUT_FILE, sep = '\t', quote = FALSE)

if (nrow(invariant) == 0) {
    message("INFO: All ", length(target_cols), " ", label, " have valid variance")
} else {
    message("INFO: Found ", nrow(invariant), " invariant ", label, ": ",
            paste(invariant$predictor, collapse = ", "))
}
