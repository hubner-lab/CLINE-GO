# factor_correlation.R — the correlation matrix of a wide numeric factor table, optionally
# joined with a second block of factors: what every correlogram draws and every correlation
# table lists (bffcfc).
# Usage: source("/pipeline/scripts/R/lib/factor_correlation.R")
# Requires: data.table
#
# Callers:
#   plot_correlation_heatmap.R  draws factor_correlation()$r
#   correlation_table.R         writes correlation_pairs(factor_correlation(...))
# on the same tables, so a correlation table holds exactly the figure's numbers:
#   Environmental (mode=climate)  A = climate_present_site.tsv,     B = NULL
#   Phenotypic    (mode=traits)   A = metadata.tsv,                 B = NULL
#   Phenotypic + climate          A = metadata_climate_valid.tsv,   B = climate_present_site.tsv
#
# The two blocks are joined on `sample` whenever both carry it; a row-count mismatch without a
# shared key is a hard error instead of a silent misalignment.

FACTOR_META_COLS <- c("site", "sample", "latitude", "longitude", "lat", "lon", "ID")

# site and sample read as character (numeric IDs must not become integers), named only when
# the table has them: fread warns about every colClasses name it cannot find, and the climate
# site table has no site column.
read_factor_block <- function(path) {
    keys <- intersect(c("site", "sample"), names(data.table::fread(path, nrows = 0L)))
    if (length(keys) == 0L) return(data.table::fread(path))
    data.table::fread(path, colClasses = stats::setNames(rep("character", length(keys)), keys))
}

# Numeric factor columns only (metadata columns dropped) — same exclusion rule as
# plot_density.R:27-31, so all three factor-characterization scripts agree on what counts as
# a factor.
factor_cols <- function(dt) {
    cols <- setdiff(colnames(dt), colnames(dt)[tolower(colnames(dt)) %in% tolower(FACTOR_META_COLS)])
    cols[vapply(cols, function(x) is.numeric(dt[[x]]), logical(1))]
}

# @param table_a path of the required wide table (numeric factor columns + metadata columns)
# @param table_b path of an optional second block, "NULL" to skip
# @return list:
#   r        correlation matrix of the non-constant factors (the figure's matrix)
#   n_obs    matrix: rows with both values present, over the rows each pair was correlated on
#   basis    matrix: "sites" for pairs correlated over one row per site, else "samples"
#   dropped  constant or all-NA factors left out
factor_correlation <- function(table_a, table_b = "NULL") {
    a <- read_factor_block(table_a)
    factors <- a[, factor_cols(a), with = FALSE]
    # Site key, carried in lockstep with `factors` — see the BLOCK-WISE WEIGHTING note
    # below. NULL when table A has no 'site' column (e.g. climate_present_site.tsv,
    # which is keyed by sample only); a fallback is derived there.
    site_vec <- if ("site" %in% colnames(a)) as.character(a$site) else NULL

    if (!identical(table_b, "NULL") && nzchar(table_b)) {
        b <- read_factor_block(table_b)
        b_factors <- b[, factor_cols(b), with = FALSE]

        if ("sample" %in% colnames(a) && "sample" %in% colnames(b)) {
            # Keyed join — the only alignment guarantee that survives either table
            # being filtered independently (e.g. climate-NA-excluded samples).
            joined <- merge(cbind(sample = a$sample,
                                  site = if (is.null(site_vec)) NA_character_ else site_vec,
                                  factors),
                            cbind(sample = b$sample, b_factors),
                            by = "sample", sort = FALSE)
            if (nrow(joined) == 0)
                stop("correlation heatmap: no samples shared between the two factor tables")
            if (nrow(joined) < nrow(a))
                message("INFO: ", nrow(a) - nrow(joined), " sample(s) present in table A dropped by the join")
            site_vec <- if (!all(is.na(joined$site))) as.character(joined$site) else NULL
            factors <- joined[, setdiff(colnames(joined), c("sample", "site")), with = FALSE]
        } else {
            # No shared key: fall back to positional binding, but only when the row
            # counts prove the two blocks describe the same samples in the same order.
            if (nrow(a) != nrow(b))
                stop("correlation heatmap: tables have no 'sample' column to join on and ",
                     "different row counts (", nrow(a), " vs ", nrow(b), ") — refusing to ",
                     "bind them positionally")
            message("INFO: no shared 'sample' key — binding ", nrow(a), " rows positionally")
            factors <- cbind(factors, b_factors)
        }
    }

    if (ncol(factors) < 2)
        stop("correlation heatmap: need at least 2 numeric factor columns, found ", ncol(factors))

    # Remove constant columns (filter NA from sd on non-numeric columns)
    std_devs <- apply(factors, 2, sd, na.rm = TRUE)
    non_constant <- names(std_devs)[!is.na(std_devs) & std_devs != 0]
    dropped <- setdiff(colnames(factors), non_constant)
    if (length(dropped) > 0)
        message("INFO: dropped ", length(dropped), " constant/all-NA factor(s): ",
                paste(dropped, collapse = ", "))
    factors <- factors[, non_constant, with = FALSE]

    message("INFO: correlation matrix ", nrow(factors), " rows x ", ncol(factors),
            " factors: ", paste(colnames(factors), collapse = ", "))

    # BLOCK-WISE WEIGHTING. Some factor columns are properties of the SITE, not of the
    # sample: every sample sequenced at a site carries that site's climate row. Over
    # per-sample rows a site with 30 individuals then counts 30x against a site with 1
    # in every |r| a user reads off this figure when curating predictors — and those |r|
    # would disagree with the site-weighted collinearity screen in pregea_rda_setup.R.
    # Measured on test_data (9 sites, 46 samples, only mildly unbalanced): max
    # |r_sample - r_site| = 0.335, enough to flip which predictor the screen keeps.
    #
    # A column is treated as site-level when it is CONSTANT WITHIN EVERY SITE, which is
    # checked from the data rather than assumed, so this stays correct for the joint
    # traits x climate table (traits vary within a site and stay per-sample) and for a
    # traits-only table (nothing is site-level, nothing changes). Where table A has no
    # 'site' column, identical rows ARE the sites — climate values are constant within
    # one — so row identity is the fallback key.
    site_key <- if (is.null(site_vec)) do.call(paste, c(as.list(factors), sep = "\r")) else site_vec
    n_sites  <- length(unique(site_key))
    site_level <- vapply(colnames(factors), function(cn) {
        all(vapply(split(factors[[cn]], site_key),
                   function(v) length(unique(v[!is.na(v)])) <= 1, logical(1)))
    }, logical(1))

    cor_matrix <- cor(factors, use = 'pairwise.complete.obs')
    n_obs <- crossprod(!is.na(as.matrix(factors)))
    basis <- matrix("samples", ncol(factors), ncol(factors),
                    dimnames = dimnames(cor_matrix))
    sl_cols <- names(site_level)[site_level]
    if (length(sl_cols) > 1 && n_sites >= 3 && n_sites < nrow(factors)) {
        first_of_site <- !duplicated(site_key)
        site_rows <- as.data.frame(factors)[first_of_site, sl_cols, drop = FALSE]
        cor_matrix[sl_cols, sl_cols] <- cor(site_rows, use = 'pairwise.complete.obs')
        n_obs[sl_cols, sl_cols] <- crossprod(!is.na(as.matrix(site_rows)))
        basis[sl_cols, sl_cols] <- "sites"
        message("INFO: ", length(sl_cols), " site-level factor(s) correlated over ",
                n_sites, " distinct sites (", nrow(factors), " samples)",
                if (any(!site_level)) paste0("; ", sum(!site_level),
                    " per-sample factor(s) over samples") else "")
    } else if (length(sl_cols) > 1 && n_sites < 3) {
        message("WARNING: only ", n_sites, " distinct site(s) — correlations left ",
                "sample-weighted (a site-level correlation needs >= 3 sites)")
    }

    list(r = cor_matrix, n_obs = n_obs, basis = basis, dropped = dropped)
}

# Long table of a factor_correlation() result: one row per pair of factors (i < j, in the
# matrix's column order), with r, the rows it was computed over and their basis.
correlation_pairs <- function(fc) {
    f <- colnames(fc$r)
    idx <- which(upper.tri(fc$r), arr.ind = TRUE)
    idx <- idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
    data.table::data.table(factor_1 = f[idx[, "row"]],
                           factor_2 = f[idx[, "col"]],
                           r        = fc$r[idx],
                           basis    = fc$basis[idx],
                           n_obs    = as.integer(fc$n_obs[idx]))
}
