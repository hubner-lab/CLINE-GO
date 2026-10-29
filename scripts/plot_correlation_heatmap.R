#!/usr/bin/env Rscript
# Correlogram of a wide numeric factor table, optionally joined with a second
# block of factors.
#
# Generalized (Phase 3) from the old climate-only signature, which took a
# climate table as arg 1 and positionally cbind()ed the trait columns of a
# metadata table onto it. Two callers now need different combinations:
#
#   Environmental (mode=climate)  A = climate_present_site.tsv, B = NULL
#   Phenotypic    (mode=traits)   A = metadata.tsv,             B = NULL
#   Phenotypic + climate          A = metadata_climate_valid.tsv, B = climate_present_site.tsv
#
# The positional cbind() it replaces was only correct because the climate rule
# was careful to pass the climate-valid metadata subset (same row count AND
# same row order). That invariant lived in a Snakefile comment, not in code.
# Here the two blocks are joined on `sample` whenever both carry it, and a
# row-count mismatch without a shared key is a hard error instead of a silent
# misalignment.

library(ggplot2)
library(data.table)
library(dplyr)
library(qs)
library(ggcorrplot)

args = commandArgs(trailingOnly=TRUE)
#################################
TABLE_A   = args[1]    # required wide table (numeric factor columns + metadata columns)
TABLE_B   = args[2]    # optional second block, "NULL" to skip
OUT_PNG   = args[3]
INTER_DIR = args[4]
TITLE     = if (length(args) >= 5 && nzchar(args[5])) args[5] else "Correlogram of factors"
#################################

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- file.path(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Plot theme
plot_theme <-
  theme_classic(base_size = 8, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = 'bold', size = 10, angle = -90),
        axis.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 12))

#################################

META_COLS <- c("site", "sample", "latitude", "longitude", "lat", "lon", "ID")

`%||%` <- function(x, y) if (is.null(x)) y else x

read_block <- function(path) {
  fread(path, colClasses = c("site" = "character", "sample" = "character"))
}

# Numeric factor columns only (metadata columns dropped) — same exclusion rule
# as plot_density.R:27-31, so all three factor-characterization scripts agree
# on what counts as a factor.
factor_cols <- function(dt) {
  cols <- setdiff(colnames(dt), colnames(dt)[tolower(colnames(dt)) %in% tolower(META_COLS)])
  cols[vapply(cols, function(x) is.numeric(dt[[x]]), logical(1))]
}

a <- read_block(TABLE_A)
factors <- a[, factor_cols(a), with = FALSE]
# Site key, carried in lockstep with `factors` — see the BLOCK-WISE WEIGHTING note
# below. NULL when table A has no 'site' column (e.g. climate_present_site.tsv,
# which is keyed by sample only); a fallback is derived there.
site_vec <- if ("site" %in% colnames(a)) as.character(a$site) else NULL

if (!identical(TABLE_B, "NULL") && nzchar(TABLE_B)) {
  b <- read_block(TABLE_B)
  b_factors <- b[, factor_cols(b), with = FALSE]

  if ("sample" %in% colnames(a) && "sample" %in% colnames(b)) {
    # Keyed join — the only alignment guarantee that survives either table
    # being filtered independently (e.g. climate-NA-excluded samples).
    joined <- merge(cbind(sample = a$sample, site = site_vec %||% NA_character_, factors),
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
factors <- factors %>% dplyr::select(all_of(non_constant))

message("INFO: correlation matrix ", nrow(factors), " rows x ", ncol(factors),
        " factors: ", paste(colnames(factors), collapse = ", "))

# Compute correlation matrix.
#
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
site_key <- site_vec %||% do.call(paste, c(as.list(factors), sep = "\r"))
n_sites  <- length(unique(site_key))
site_level <- vapply(colnames(factors), function(cn) {
    all(vapply(split(factors[[cn]], site_key),
               function(v) length(unique(v[!is.na(v)])) <= 1, logical(1)))
}, logical(1))

cor_matrix <- cor(factors, use = 'pairwise.complete.obs')
sl_cols <- names(site_level)[site_level]
if (length(sl_cols) > 1 && n_sites >= 3 && n_sites < nrow(factors)) {
    first_of_site <- !duplicated(site_key)
    cor_matrix[sl_cols, sl_cols] <-
        cor(as.data.frame(factors)[first_of_site, sl_cols, drop = FALSE],
            use = 'pairwise.complete.obs')
    message("INFO: ", length(sl_cols), " site-level factor(s) correlated over ",
            n_sites, " distinct sites (", nrow(factors), " samples)",
            if (any(!site_level)) paste0("; ", sum(!site_level),
                " per-sample factor(s) over samples") else "")
} else if (length(sl_cols) > 1 && n_sites < 3) {
    message("WARNING: only ", n_sites, " distinct site(s) — correlations left ",
            "sample-weighted (a site-level correlation needs >= 3 sites)")
}

# Build the heatmap using ggcorrplot
gHM <- ggcorrplot(cor_matrix,
                  type = "lower",
                  lab = TRUE,
                  digits = 1,
                  lab_size = 4) +
  ggtitle(TITLE) +
  plot_theme

# Determine size based on number of factors
scale_factor <- 0.7
width <- max(4, ncol(factors) * scale_factor)
height <- max(4, ncol(factors) * scale_factor)

# Save
ggsave(OUT_PNG, gHM, width = width, height = height, units = "in")
ggsave(OUT_SVG, gHM, width = width, height = height, units = "in",
       device = svglite::svglite, bg = 'transparent')
qsave(gHM, OUT_QS)
