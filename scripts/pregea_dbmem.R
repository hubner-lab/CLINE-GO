#!/usr/bin/env Rscript
# pregea_dbmem.R — SHARED spatial-vector artifact (dbMEM eigenvectors),
# reused by preGEA's varpart block AND (follow-up, see docs/rda_research.md
# A20/C4) Gradient Forest, replacing gradient_forest_model.R's ad hoc
# raw-Euclidean-dist + "keep half the positive eigenvalues" heuristic.
#
# Ports rda_varpart_lasky.Rmd's PCNM block (geodesic distance via terra,
# MST-truncation threshold, keep-all-positive-eigenvalues), upgraded from
# vegan::pcnm() to adespatial::dbmem() per the spec.
#
# Inputs are COORDINATES ONLY (deliberate) — cheap enough to be pulled into
# maladaptation mode on demand once the Gradient Forest swap lands, without
# dragging genotype work into that DAG.
#
# ALWAYS site level, unconditionally — climate/spatial predictors are
# constant within a site, so sample-level dbMEM would replicate each value
# ~N times per site and artificially inflate the shared variance fraction
# (the pseudoreplication warned about in rda_varpart_lasky.Rmd). Sites are
# collapsed to their coordinate CENTROID (mean lat/lon) so a site is one
# point regardless of how many samples it has, down to n=1. Site-level MEMs
# are BROADCAST back to samples by site match, so the output artifact is
# ALWAYS sample-indexed — consumers (varpart, future GF) must match() on
# `sample`, never assume row order (a wider coord-valid sample set may be
# joining against this climate-valid one). Below 3 sites (or if dbmem()
# degenerates to rank 0) the script writes a skip record instead of crashing
# — see write_skip() — surfaced to the user as a warning badge in the Shiny
# Climate tab, not as a config knob.

suppressPackageStartupMessages({
    library(data.table)
    library(vegan)
    library(adespatial)
    library(terra)
    library(ggplot2)
})

source("/pipeline/scripts/R/utils/theme_clinego.R")

args <- commandArgs(trailingOnly = TRUE)
################################################################################
METADATA_CLIMATE_VALID <- args[1]
CLIMATE_VALID_SAMPLES  <- args[2]
OUT_VECTORS_TSV        <- args[3]
OUT_DIAG_TSV           <- args[4]
OUT_PNG                <- args[5]
INTER_DIR              <- args[6]
MIN_SITES              <- 3L                 # hard floor — spantree()/dbmem() need >=3 points
################################################################################

for (d in c(dirname(OUT_VECTORS_TSV), dirname(OUT_DIAG_TSV), dirname(OUT_PNG), INTER_DIR))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)

diag <- list()
add_diag <- function(k, v) diag[[k]] <<- v

meta <- fread(METADATA_CLIMATE_VALID, colClasses = c(site = "character", sample = "character"))
sample_order <- fread(CLIMATE_VALID_SAMPLES, header = FALSE, colClasses = "character")$V2
meta <- meta[match(sample_order, sample)]   # enforce CLIMATE_VALID_SAMPLES order
n_rows <- nrow(meta)

# Site centroid (mean lat/lon per site) — a site is ONE point regardless of
# how many samples it has (down to n=1), and regardless of per-sample GPS
# jitter within a site. n_unique_coords is recorded alongside for diagnostics
# only; the guard below is on n_sites, since the two can diverge in either
# direction (two sites sharing one coordinate; one site with jittered coords).
site_coords <- meta[, .(latitude = mean(latitude), longitude = mean(longitude)), by = site]
n_sites <- nrow(site_coords)
n_unique_coords <- length(unique(paste(meta$latitude, meta$longitude)))
add_diag("n_rows", n_rows)
add_diag("n_sites", n_sites)
add_diag("n_unique_coords", n_unique_coords)
add_diag("spatial_level", "site")
message("INFO: preGEA dbMEM — level=site n_rows=", n_rows, " n_sites=", n_sites,
       " n_unique_coords=", n_unique_coords)

write_skip <- function(status) {
    add_diag("status", status)
    add_diag("mst_threshold_m", NA_real_); add_diag("mst_threshold_km", NA_real_)
    add_diag("n_mem_total", 0L); add_diag("n_mem_positive", 0L)
    out <- meta[, .(sample, site)]
    fwrite(out, OUT_VECTORS_TSV, sep = "\t", quote = FALSE)
    diag_dt <- data.table(names(diag), sapply(diag, as.character)); setnames(diag_dt, c("key", "value"))
    fwrite(diag_dt, OUT_DIAG_TSV, sep = "\t", quote = FALSE)
    # Empty-state placeholder — this IS the entire plot content (CLAUDE.md Rule 8 exception).
    g <- ggplot() + annotate("text", x = 0, y = 0, label = paste0(
            "dbMEM skipped: ", status, "\n(", n_sites, " site(s), need >= ", MIN_SITES, ")")) +
        theme_void()
    ggsave(OUT_PNG, g, width = 6, height = 4, dpi = 150)
    ggsave(OUT_SVG, g, width = 6, height = 4, device = svglite::svglite, bg = "transparent")
    message("INFO: dbMEM skipped (", status, ") — wrote 0-MEM vectors file")
}

if (n_sites < MIN_SITES) {
    write_skip("too_few_sites")
    quit(status = 0)
}

################################################################################
# Point set = site centroids (always)
################################################################################
pt_ids <- site_coords$site
lon <- site_coords$longitude; lat <- site_coords$latitude
n_pts <- length(pt_ids)

################################################################################
# Geodesic distance (terra::distance on EPSG:4326 — NOT raw Euclidean, unlike
# gradient_forest_model.R's current heuristic) + MST truncation
################################################################################
pts_sv <- terra::vect(cbind(lon, lat), crs = "EPSG:4326")
d_mat  <- as.matrix(terra::distance(pts_sv))
diag(d_mat) <- 0
rownames(d_mat) <- colnames(d_mat) <- pt_ids

sp_tree   <- vegan::spantree(as.dist(d_mat))
threshold <- max(sp_tree$dist)
add_diag("mst_threshold_m", threshold)
add_diag("mst_threshold_km", threshold / 1000)
message("INFO: MST truncation threshold = ", round(threshold / 1000, 2), " km")

################################################################################
# adespatial::dbmem() — keep ALL positive-eigenvalue MEMs
################################################################################
mem_obj <- tryCatch(
    adespatial::dbmem(as.dist(d_mat), thresh = threshold, MEM.autocor = "positive", silent = TRUE),
    error = function(e) { message("WARNING: dbmem() failed: ", e$message); NULL }
)

if (is.null(mem_obj) || ncol(as.data.frame(mem_obj)) == 0) {
    write_skip("degenerate_rank")
    quit(status = 0)
}

mem_df <- as.data.frame(mem_obj)
n_mem <- ncol(mem_df)
mem_eigs <- attr(mem_obj, "values")
n_mem_positive <- if (!is.null(mem_eigs)) sum(mem_eigs > 0) else n_mem
colnames(mem_df) <- paste0("MEM", seq_len(n_mem))
mem_df$.pt_id <- pt_ids

add_diag("n_mem_total", n_mem)
add_diag("n_mem_positive", n_mem_positive)
add_diag("status", "ok")
message("INFO: dbMEM produced ", n_mem, " positive-autocorrelation MEM(s)")

################################################################################
# Broadcast site-level MEMs back to samples (always sample-indexed output)
################################################################################
out <- merge(meta[, .(sample, site)], mem_df, by.x = "site", by.y = ".pt_id", all.x = TRUE)
setcolorder(out, c("sample", "site", paste0("MEM", seq_len(n_mem))))
out <- out[match(sample_order, sample)]   # restore CLIMATE_VALID_SAMPLES order
fwrite(out, OUT_VECTORS_TSV, sep = "\t", quote = FALSE)

diag_dt <- data.table(names(diag), sapply(diag, as.character)); setnames(diag_dt, c("key", "value"))
fwrite(diag_dt, OUT_DIAG_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote ", nrow(out), " sample rows x ", n_mem, " MEM(s) to ", OUT_VECTORS_TSV)

################################################################################
# Screeplot
################################################################################
eig_dt <- data.table(mem = seq_len(n_mem), eigenvalue = if (!is.null(mem_eigs)) mem_eigs[seq_len(n_mem)] else rep(NA_real_, n_mem))
g <- ggplot(eig_dt, aes(x = factor(mem), y = eigenvalue)) +
    geom_col(fill = CLINEGO_RETAINED) +
    labs(x = "MEM", y = "Eigenvalue (Moran's I)",
        title = "dbMEM spatial eigenvector screeplot",
        subtitle = sprintf("site-level | %d positive-autocorrelation MEM(s) | MST threshold=%.1f km",
                           n_mem_positive, threshold / 1000)) +
    theme_clinego()
ggsave(OUT_PNG, g, width = 7, height = 5, dpi = 300)
ggsave(OUT_SVG, g, width = 7, height = 5, device = svglite::svglite, bg = "transparent")

message("INFO: preGEA dbMEM complete.")
