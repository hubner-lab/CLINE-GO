# Shared fixture builders and subprocess runner for the CLI wrapper smoke tests.
#
# Sourced by BOTH test roots:
#   tests/testthat/test-cli-wrappers.R   quick gate — no genomics tool needed
#   tests/heavy/test-heavy-wrappers.R    opt-in heavy gate — plink/LEA/vcftools/EMMAX
#
# It lives in tests/lib/ rather than tests/testthat/helper-*.R on purpose: the
# heavy root must be able to source it without reading a file named as the quick
# root's helper, and test_dir() must not source it twice.
#
# Side-effect free at source time apart from defining functions and the two path
# constants below.
#
# THE CWD SANDBOX, and why it is the real containment mechanism.
# run_wrapper() runs each script with its working directory set to an empty
# throwaway dir, and expect_wrapper_ok() then asserts that dir is STILL empty.
# Several wrappers resolve a path relative to the cwd — generate_simdata.R:8
# defaults OUTDIR to "data/", and base-graphics calls drop Rplots.pdf where they
# stand (geometric_offset.R:459, rda_offset.R:397). Without the sandbox those
# land in the mounted repo. A `git status` check cannot substitute: data/ and
# *_results/ are gitignored (.gitignore:4,23), so git is blind to exactly the
# writes that matter most.

REPO    <- getOption("clinego.repo_root", "/pipeline")
SCRIPTS <- file.path(REPO, "scripts")

# Wrappers that must NEVER appear in a spec table. Both take no argv at all and
# read-modify-write the tracked SIMDATA fixtures in place (add_related_samples.R:30-31,
# add_pregea_sites.R:38-39). Asserted by a test, not left to review discipline.
WRAPPER_DENYLIST <- c("add_related_samples.R", "add_pregea_sites.R")

# ------------------------------------------------------------------ fixtures

write_tsv <- function(dt, path) {
    data.table::fwrite(dt, path, sep = "\t")
    path
}

fx_metadata <- function(d, n = 12) {
    dt <- data.table::data.table(
        site      = rep(c("NEG", "TAV", "GAL"), length.out = n),
        sample    = sprintf("ID%03d", seq_len(n)),
        latitude  = 30 + seq_len(n) * 0.1,
        longitude = 34 + seq_len(n) * 0.1,
        height    = as.numeric(seq_len(n)),
        flowering_time = as.numeric(rev(seq_len(n)))
    )
    write_tsv(dt, file.path(d, "metadata.tsv"))
}

fx_climate_site <- function(d) {
    dt <- data.table::data.table(
        sample = sprintf("ID%03d", 1:8),
        bio_1  = c(1.0, 2.0, 3.0, 4.0, 2.5, 3.5, 1.5, 4.5),
        bio_2  = c(9.0, 7.0, 5.0, 3.0, 8.0, 4.0, 6.0, 2.0),
        bio_12 = rep(5.0, 8)                       # invariant on purpose
    )
    write_tsv(dt, file.path(d, "climate_present_site.tsv"))
}

fx_pvalues <- function(d, name = "pvalues.tsv") {
    n <- 20
    dt <- data.table::data.table(
        SNPID = sprintf("snp%02d", seq_len(n)),
        chr   = rep(c("1", "2"), each = n / 2),
        pos   = seq_len(n) * 10000L,
        bio_1 = c(1e-9, 1e-8, runif(n - 2, 0.05, 1)),
        bio_2 = c(runif(n - 2, 0.05, 1), 1e-9, 1e-8)
    )
    write_tsv(dt, file.path(d, name))
}

fx_sig_snps <- function(path, method, traits = c("bio_1", "bio_2")) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    dt <- data.table::data.table(
        SNPID           = c("snp01", "snp02", "snp11", "snp12"),
        chr             = c("1", "1", "2", "2"),
        pos             = c(10000L, 20000L, 110000L, 120000L),
        pvalue          = c(1e-9, 1e-8, 1e-9, 1e-8),
        pval_threshold  = 0.0025,
        method          = method,
        trait           = rep(traits, length.out = 4),
        overlap_traits  = "",
        overlap_snps    = "",
        overlap_distance = 10000L
    )
    write_tsv(dt, path)
}

fx_regions <- function(d) {
    dt <- data.table::data.table(
        region_id = c("1_5000-30000", "2_105000-130000"),
        trait     = c("bio_1", "bio_2"),
        chr       = c("1", "2"),
        start     = c(5000L, 105000L),
        end       = c(30000L, 130000L)
    )
    write_tsv(dt, file.path(d, "regions.tsv"))
}

fx_gff <- function(d) {
    p <- file.path(d, "genes.gff3")
    lines <- c(
        "##gff-version 3",
        paste("1", "src", "gene", "8000",  "15000", ".", "+", ".", "ID=g1;Name=ALPHA", sep = "\t"),
        paste("1", "src", "exon", "9000",  "10000", ".", "+", ".", "ID=e1;Parent=g1",  sep = "\t"),
        paste("2", "src", "gene", "110000", "125000", ".", "-", ".", "ID=g2;Name=BETA", sep = "\t")
    )
    writeLines(lines, p)
    p
}

fx_allsnps <- function(d) {
    p <- file.path(d, "all.vcfsnp")
    writeLines(c("1 9500 snp01 C A . . PR GT",
                 "1 12000 snp02 G T . . PR GT",
                 "2 115000 snp11 A C . . PR GT"), p)
    p
}

fx_sample_list <- function(d, name, samples) {
    p <- file.path(d, name)
    writeLines(paste(samples, samples), p)   # plink --keep: FID IID, no header
    p
}

# Everything write_summary.R's processing branch reads (args 3-11), with the
# sample counts chosen so the accounting identity total == after + removed +
# het + related closes (10 = 8 + 2 + 0 + 0). `truncate_keep` shortens ONLY the
# keep-list to the shape plink --keep leaves behind when it cannot find an ID
# (audit 2026-09-13 B41): samples_total then reads 7 while 8 samples passed,
# and the script must refuse to write that.
fx_processing_summary_inputs <- function(d, truncate_keep = FALSE) {
    vcf_lines <- c("##fileformat=VCFv4.2",
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tID001",
                   paste0("1\t", 1:5 * 100, "\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"))
    vcf_filt <- file.path(d, "filt.vcf"); writeLines(vcf_lines, vcf_filt)
    vcf_ld   <- file.path(d, "ld.vcf");   writeLines(vcf_lines[1:4], vcf_ld)
    ids  <- sprintf("ID%03d", 1:10)
    keep <- if (truncate_keep) ids[1:7] else ids
    samples_list     <- fx_sample_list(d, "samples.list", keep)
    samples_filtered <- fx_sample_list(d, "samples_filtered.list", ids[1:8])
    samples_removed  <- file.path(d, "samples_removed.list")
    writeLines(paste(ids[9:10], ids[9:10], c(0.31, 0.44)), samples_removed)
    stats <- file.path(d, "raw_stats.txt")
    writeLines(c("SN\t0\tnumber of SNPs:\t5", "TSTV\t0\t3\t2\t1.50"), stats)
    filt <- write_tsv(data.table::data.table(
                stage = c("Raw VCF", "After sample missingness filter"),
                n_samples = c(10L, 8L), n_snps = c(5L, 5L)),
            file.path(d, "filtering_summary.tsv"))
    c(vcf_filt, vcf_ld, samples_list, samples_filtered, samples_removed, stats, filt)
}

# Everything compare_offsets.R reads, on a 2 x 3 raster (6 cells). Returns a
# named list of paths. Hand-placed so the expected ExDet classes are known
# WITHOUT re-deriving them with the script's own formulas: six reference sites
# (bio_1 10..20, bio_2 100..180), then future cells
#   1 (15,140) analog   2 (5,140) type 1 (bio_1 below min: NT1 = -0.5)
#   3 (20,100) type 2   4 (10,180) type 2 (inside both ranges, Mahalanobis
#   ratio 2.3 / 2.7 > 1 against the reference maximum)   5 (12,120) analog
#   6 (18,170) analog  -> 50 % novel, 16.7 % type 1, 33.3 % type 2.
# The two offset models rank the six sites identically, so Kendall's W is 1
# by definition (the audit's B2 fixture: the transposed call gave 0).
fx_compare_offsets <- function(d) {
    sites <- paste0("S", 1:6)
    meta <- data.table::data.table(
        site = rep(sites, each = 2), sample = sprintf("ID%03d", 1:12),
        latitude = rep(c(30.1, 30.9, 31.6, 32.2, 32.8, 33.3), each = 2),
        longitude = rep(c(34.5, 34.8, 35.0, 35.2, 35.3, 35.5), each = 2))
    meta_p <- write_tsv(meta, file.path(d, "metadata_climate_valid.tsv"))
    ref_b1 <- c(10, 12, 14, 16, 18, 20); ref_b2 <- c(100, 150, 120, 180, 110, 160)
    clim_site <- data.table::data.table(sample = meta$sample,
        bio_1 = rep(ref_b1, each = 2), bio_2 = rep(ref_b2, each = 2))
    clim_site_p <- write_tsv(clim_site, file.path(d, "climate_present_site.tsv"))
    pres_all <- data.table::data.table(ID = 1:6, bio_1 = ref_b1, bio_2 = ref_b2)
    fut_all  <- data.table::data.table(ID = 1:6,
        bio_1 = c(15, 5, 20, 10, 12, 18), bio_2 = c(140, 140, 100, 180, 120, 170))
    pres_all_p <- write_tsv(pres_all, file.path(d, "climate_present_all.tsv"))
    fut_all_p  <- write_tsv(fut_all,  file.path(d, "climate_future_all.tsv"))
    r <- terra::rast(nrows = 2, ncols = 3, xmin = 34, xmax = 37, ymin = 30, ymax = 32,
                     crs = "EPSG:4326", vals = 1:6)
    ras_p <- file.path(d, "present.tif")
    terra::writeRaster(r, ras_p, overwrite = TRUE)
    site_tbl <- function(off) data.table::data.table(
        site = meta$site, sample = meta$sample, genetic_offset = rep(off, each = 2))
    a_site <- write_tsv(site_tbl(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)), file.path(d, "a_site.tsv"))
    b_site <- write_tsv(site_tbl(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)), file.path(d, "b_site.tsv"))
    a_map <- file.path(d, "a_map.tsv"); b_map <- file.path(d, "b_map.tsv")
    writeLines(c("1\t2\t3", "4\t5\t6"), a_map)
    writeLines(c("2\t4\t6", "8\t10\t12"), b_map)
    list(cache = file.path(d, "cache"), novelty = file.path(d, "novelty"),
         a_site = a_site, a_map = a_map, b_site = b_site, b_map = b_map,
         clim_site = clim_site_p, pres_all = pres_all_p, fut_all = fut_all_p,
         raster = ras_p, meta = meta_p)
}

# Every input of pregea_varpart.R (25 positionals, in argv order — see the
# script header). n_sites sites x n_per_site samples, one coordinate and one
# climate row per site (site properties, as the pipeline writes them), dbMEM
# vectors broadcast per site the way pregea_dbmem.R emits them, a Q-matrix
# and a LEA-shaped projections file in the same sample order. Returns the
# argv vector minus the five output paths + plot/inter dirs, which the spec
# row supplies. With 8 sites and 2 predictors the site-level residual df is 5
# (status "ok"); with 3 sites it is 0 (status "insufficient_sites").
fx_varpart_inputs <- function(d, n_sites = 8L, n_per_site = 2L, k = 3L, seed = 7L) {
    set.seed(seed)
    sites   <- sprintf("S%02d", seq_len(n_sites))
    samples <- sprintf("ID%03d", seq_len(n_sites * n_per_site))
    site_of <- rep(sites, each = n_per_site)
    lat <- 30 + seq_len(n_sites) * 0.3; lon <- 34 + seq_len(n_sites) * 0.2
    meta <- write_tsv(data.table::data.table(
        site = site_of, sample = samples,
        latitude = rep(lat, each = n_per_site), longitude = rep(lon, each = n_per_site)),
        file.path(d, "metadata_climate_valid.tsv"))
    clim_site <- data.table::data.table(site = sites,
        bio_1 = round(rnorm(n_sites), 4), bio_12 = round(rnorm(n_sites), 4))
    clim <- write_tsv(data.table::data.table(sample = samples,
        bio_1  = clim_site$bio_1[match(site_of, sites)],
        bio_12 = clim_site$bio_12[match(site_of, sites)]),
        file.path(d, "climate_present_site_scaled.tsv"))
    mem_site <- matrix(round(rnorm(n_sites * 2), 4), ncol = 2)
    dbmem <- write_tsv(data.table::data.table(sample = samples, site = site_of,
        MEM1 = mem_site[match(site_of, sites), 1], MEM2 = mem_site[match(site_of, sites), 2]),
        file.path(d, "dbmem_vectors.tsv"))
    q <- matrix(runif(length(samples) * k), ncol = k); q <- q / rowSums(q)
    cl <- data.table::data.table(sample = samples, site = site_of)
    for (i in seq_len(k)) cl[[paste0("C", i)]] <- round(q[, i], 4)
    clusters <- write_tsv(cl, file.path(d, "clusters_K3.tsv"))
    valid <- file.path(d, "climate_valid_samples.list")
    writeLines(paste(0, samples), valid)
    order_p <- file.path(d, "samples_order.list")
    writeLines(samples, order_p)
    proj <- fx_projections(d, n = length(samples), n_pc = 5)
    eig  <- fx_eigenvalues(d, n = 5)
    c(proj, eig, dbmem, clusters, clim, valid, order_p, "bio_1,bio_12",
      "pcs", "0.8", "5", "2", "NULL", "0.05", "99", "99", "42", meta)
}

# --- fixtures added 2026-09-12 with the second wrapper batch ---------------

# CHR/POS/MAF, the shape compute_wza.R expects (its own arg names, uppercase).
fx_maf <- function(d, name = "maf_pos.tsv") {
    n <- 20
    dt <- data.table::data.table(
        CHR = rep(c("1", "2"), each = n / 2),
        POS = seq_len(n) * 10000L,
        MAF = seq(0.05, 0.5, length.out = n)
    )
    write_tsv(dt, file.path(d, name))
}

# LFMM genotype matrix: space-separated, no header, no row names. 0/1/2 and 9 = NA.
fx_lfmm_matrix <- function(d, name = "geno.lfmm", n_ind = 6, n_snp = 8) {
    set.seed(42)
    m <- matrix(sample(0:2, n_ind * n_snp, replace = TRUE), nrow = n_ind)
    p <- file.path(d, name)
    write.table(m, p, sep = " ", row.names = FALSE, col.names = FALSE)
    p
}

# sNMF Q-matrix exactly as extract_clusters.R writes it: sample, site, C1..Ck.
# The `site` column is load-bearing — plot_structure.R:25 does select(-site) and
# errors without it.
fx_clusters <- function(d, k = 3, n = 8, name = NULL) {
    if (is.null(name)) name <- paste0("clusters_K", k, ".tsv")
    q <- matrix(runif(n * k), nrow = n)
    q <- q / rowSums(q)
    dt <- data.table::data.table(
        sample = sprintf("ID%03d", seq_len(n)),
        site   = rep(c("NEG", "TAV", "GAL"), length.out = n)
    )
    for (i in seq_len(k)) dt[[paste0("C", i)]] <- q[, i]
    write_tsv(dt, file.path(d, name))
}

# selected_snps.tsv: SNPID, chr, pos, one column PER METHOD, min_pvalue last.
# The per-method cells hold comma-separated TRAIT NAMES (or "" / literal `""`),
# not p-values — the schema misread that cost a session during Tier 4. Verified
# against SIMDATA_results/GEA/tables/selected_snps.tsv.
fx_selected_snps <- function(d, name = "selected_snps.tsv",
                             methods = c("EMMAX", "LFMM"),
                             traits = c("bio_1", "bio_2")) {
    dt <- data.table::data.table(
        SNPID = c("snp01", "snp02", "snp11", "snp12"),
        chr   = c("1", "1", "2", "2"),
        pos   = c(10000L, 20000L, 110000L, 120000L)
    )
    for (i in seq_along(methods)) {
        # Row 1 is hit by every method (a cross-method overlap), row 4 by none.
        dt[[methods[i]]] <- c(traits[1],
                              if (i == 1) paste(traits, collapse = ",") else "",
                              traits[min(i, length(traits))],
                              "")
    }
    dt$min_pvalue <- c(1e-9, 1e-8, 1e-7, 1e-6)
    write_tsv(dt, file.path(d, name))
}

# One number per line — what LEA writes and what plot_pregea_screeplot.R /
# plot_pca_structure.R read back with readLines() / fread(header = FALSE).
fx_eigenvalues <- function(d, name = "eigenvalues.txt", n = 10) {
    p <- file.path(d, name)
    writeLines(format(sort(runif(n, 1, 100), decreasing = TRUE), trim = TRUE), p)
    p
}

# Headerless space-separated PCA projections, n rows x n_pc columns.
fx_projections <- function(d, name = "projections.txt", n = 8, n_pc = 5) {
    p <- file.path(d, name)
    m <- matrix(round(rnorm(n * n_pc), 4), nrow = n)
    write.table(m, p, sep = " ", row.names = FALSE, col.names = FALSE)
    p
}

# A non-empty stand-in for a plot that a script only file.exists()-tests
# (write_summary.R:342 is the case). Deliberately NOT a real PNG: opening a
# graphics device at a small size raises "figure margins too large" in the
# container's default device, and no consumer of this fixture decodes the file.
fx_dummy_png <- function(d, name = "placeholder.png") {
    p <- file.path(d, name)
    writeLines("not a real PNG - existence-only fixture", p)
    p
}

# INTER_DIR arguments MUST carry a trailing slash. Several scripts build their
# intermediate path by string concatenation rather than file.path() — e.g.
# plot_density.R:19, plot_structure.R:18, plot_correlation_heatmap.R:37 — and
# work in production only because common.smk defines INTER with a trailing
# slash. Filed; until fixed, fixtures must honour the convention.
fx_inter_dir <- function(d, name = "inter") {
    p <- file.path(d, name)
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    paste0(p, "/")
}

# Climate/varpart summary inputs. write_summary.R's `climate` branch reads five
# tables, each with a different shape.
#
# NOTE ON NAMING: there is NO `varpart` MODE. write_summary.R dispatches on
# processing|prestructure|pregea|climate|traits|structure|gea|gwas|
# maladaptation|gea_x_gwas, and an unknown mode warns and quits 0 (:606-609).
# The variance-partition table is args[5] INSIDE the `climate` branch (:257).
fx_varpart <- function(d, name = "variance_partition.tsv") {
    # The REAL schema pregea_varpart.R writes (verified against
    # SIMDATA_results/climate/tables/varpart/variance_partition.tsv): percents, a
    # p_value per unique fraction, a `group`, and "Unexplained" (not "Residual")
    # — consumers key on that exact component name (mod_climate.R,
    # mod_dashboard_outputs.R, write_summary.R). unit / n_units / df_env were
    # added 2026-09-13 with the site-level fit. Before that date this fixture
    # carried fractions, "Residual" and no p_value, so the write_summary.R row
    # that consumes it never exercised the real keys.
    write_tsv(data.table::data.table(
        group        = c("Unique to one factor", "Unique to one factor",
                         "Shared between factors", "Unexplained"),
        component    = c("Climate", "Structure", "Climate \u2229 Structure", "Unexplained"),
        variance_pct = c(12, 31, 7, 50),
        p_value      = c(0.012, 0.001, NA, NA),
        model        = "2-way (climate + structure)", status = "ok",
        unit = "site", n_units = 8L, df_env = 7L), file.path(d, name))
}

fx_dbmem_diag <- function(d, name = "dbmem_diagnostics.tsv") {
    # Long key/value, not one row per metric: write_summary.R:282-286 does
    # setNames(dbmem$value, dbmem$key) and picks six keys out by name.
    #
    # `key` cannot be passed to data.table() as a column name — it is that
    # function's own KEY argument, so data.table(key = c(...)) tries to set a
    # key from those strings and dies with "some columns are not in the
    # data.table". Build it under another name and rename.
    dt <- data.table::data.table(
        k = c("spatial_level", "n_sites", "n_unique_coords",
              "mst_threshold_km", "n_mem_positive", "status"),
        value = c("site", "9", "9", "48.2", "4", "ok"))
    data.table::setnames(dt, "k", "key")
    write_tsv(dt, file.path(d, name))
}

# The eight columns of ld_decay_analyze.R's EMPTY skeleton (:247-250), which are
# also what write_summary.R's structure branch reads back (:383-397) — one
# fixture for both sides, so a column rename cannot pass one and fail the other
# silently. The group == "All" & scope == "genome_wide" row is what :385 picks.
#
# NOTE the populated table has TEN columns: the per-result data.table (:183-193)
# adds C_hat and max_dist_bp, which the empty skeleton omits. Filed. Every
# current consumer reads only the four columns above, so it is latent, but a
# reader that globs columns sees a different schema depending on whether the run
# found any data.
fx_ld_decay_table <- function(d, name = "ld_decay_half_distances.tsv") {
    write_tsv(data.table::data.table(
        group = c("All", "All", "NEG"),
        scope = c("genome_wide", "per_chromosome", "genome_wide"),
        half_decay_bp = c(1240, 1180, 990),
        r2_02_bp = c(3100, 2980, 2450),
        n_samples = c(50L, 50L, 17L),
        n_pairs = c(12000L, 4100L, 3900L),
        r2_intercept = c(0.62, 0.60, 0.66),
        method = c("hill_weir", "hill_weir", "loess")),
        file.path(d, name))
}

# An ALIGNED metadata / clusters / climate triple for mantel_test.R. Returns
# c(samples, clusters, env).
#
# No combination of the existing fixtures can drive that script, for three
# independent reasons — each of which is a hard stop(), not a warning:
#   1. fx_metadata and fx_clusters both cycle NEG/TAV/GAL = 3 sites, and
#      mantel_test.R:285-289 stop()s below 4 sites after site aggregation
#      (3 sites give 3 pairs and no meaningful permutation test). Below 8 it
#      only warns (:290-295), so 8 is the smallest quiet default.
#   2. ENV is read POSITIONALLY, not joined by sample (:219-220), and
#      :270-271 stopifnot()s that geo, env and clust all have the metadata's row
#      count. fx_metadata(n = 12) against fx_climate_site's 8 rows fails that.
#   3. Clusters are aligned by match() on sample with anyNA -> stop() (:231-236),
#      so the cluster table must name every metadata sample.
# Predictors must also vary ACROSS SITES: scale() turns a zero-variance column
# into NaN, those get dropped, and ncol(env) == 0 is another stop() (:299-311).
# fx_climate_site's bio_12 is deliberately invariant, so it cannot be reused here.
fx_geo_triple <- function(d, n_sites = 8L, n_per_site = 2L, k = 3L) {
    n <- n_sites * n_per_site
    site <- rep(sprintf("S%02d", seq_len(n_sites)), each = n_per_site)
    samp <- sprintf("ID%03d", seq_len(n))
    # One coordinate pair per site, shared by that site's samples — the real
    # shape, and what makes site aggregation meaningful.
    lat <- rep(30.5 + seq_len(n_sites) * 0.35, each = n_per_site)
    lon <- rep(34.5 + seq_len(n_sites) * 0.20, each = n_per_site)

    samples <- write_tsv(data.table::data.table(
        site = site, sample = samp, latitude = lat, longitude = lon),
        file.path(d, "metadata_climate_valid.tsv"))

    q <- matrix(0, nrow = n, ncol = k)
    for (i in seq_len(n)) {
        w <- (seq_len(k) + i) %% k + 1L
        q[i, ] <- w / sum(w)
    }
    cl <- data.table::data.table(sample = samp, site = site)
    for (j in seq_len(k)) data.table::set(cl, j = paste0("C", j), value = q[, j])
    clusters <- write_tsv(cl, file.path(d, paste0("clusters_K", k, ".tsv")))

    # Both predictors vary between sites and stay constant within one.
    env <- write_tsv(data.table::data.table(
        sample = samp,
        bio_1  = rep(18 + seq_len(n_sites) * 0.8, each = n_per_site),
        bio_12 = rep(600 - seq_len(n_sites) * 45, each = n_per_site)),
        file.path(d, "climate_present_site.tsv"))

    c(samples, clusters, env)
}

# A PopLDdecay .stat.gz as read_stat_gz() actually parses it (:39-49): gzipped,
# WITH a header, and SIX columns renamed positionally to
# dist_bp, mean_r2, mean_dprime, sum_r2, sum_dprime, n_pairs.
#
# The script's own file header (:5-6) documents FOUR columns in a different
# order ("Dist(bp) Mean_r^2 NumberPairs Sum_r^2"). That is stale and a fixture
# built from it dies inside setnames(). Filed. The in-function comment at
# :33-34 is the correct one.
#
# n_rows spans enough distance to clear two thresholds at once: bin_stat() uses
# bins of max(1000, ceiling(max_dist/100)) and fit_ld_curve() returns
# "insufficient_data" below FIVE bins (:75). dist_bp == 0 rows are dropped (:50).
fx_stat_gz <- function(d, group, chr = NULL, n_rows = 60L) {
    nm <- if (is.null(chr)) paste0(group, ".stat.gz") else paste0(group, "_", chr, ".stat.gz")
    p  <- file.path(d, nm)
    dist <- seq_len(n_rows) * 1000L
    # A decaying r^2 so the Hill-Weir nls has something to converge on.
    mean_r2 <- 0.6 * exp(-dist / 20000) + 0.02
    n_pairs <- rep(500L, n_rows)
    dt <- data.table::data.table(
        `#Dist`   = dist,
        `Mean_r^2` = round(mean_r2, 6),
        `Mean_D'`  = round(mean_r2 + 0.1, 6),
        `Sum_r^2`  = round(mean_r2 * n_pairs, 4),
        `Sum_D'`   = round((mean_r2 + 0.1) * n_pairs, 4),
        NumberPairs = n_pairs)
    con <- gzfile(p, "w")
    on.exit(close(con), add = TRUE)
    utils::write.table(dt, con, sep = "\t", quote = FALSE, row.names = FALSE)
    p
}

# A minimal but INTERNALLY CONSISTENT {PROJECT}_results/ tree for
# check_invariants.R. Returns the results directory.
#
# Only the GEA module is populated: check_invariants.R's read_tsv() returns NULL
# for a missing file and every caller skips, so an absent module costs nothing
# and a second one would only duplicate the same code path.
#
# Consistency is the whole point. The default tree must produce ZERO
# error-severity violations, which pins every cross-table relationship the
# shaping layer builds:
#   * summary counts equal the row counts of the tables they summarise
#   * a region's snp_ids name SNPs that exist in selected_snps.tsv and fall
#     inside its own bounds, and snp_count agrees with that list
#   * length == end - start (NOT end - start + 1; regions.R:90)
#   * min_pvalue equals the smallest p for that SNP across the *_sig_snps_* tables
#   * every chr in every table is in the canonical set, which is built from the
#     *_pvalues_K*.tsv tables alone
#   * genes reference regions that exist in regions_per_trait.tsv
#
# Two schema traps this fixture exists to honour, both recorded at
# check_invariants.R:12-18 as the misreads found while writing it:
#   * selected_snps.tsv's per-method columns hold TRAIT NAMES, not p-values
#   * the sig tables' overlap_traits must name OTHER traits, never the row's own
#
# violations = TRUE applies exactly two independent error-severity mutations —
# one p-value of 1.5 (check_pvalues_table) and one gene row on "chr2"
# (check_chromosome_names, both the ^chr-prefix and the not-canonical halves).
# Neither cascades into another checker, so the exit status is traceable.
fx_results_tree <- function(d, project = "FX", violations = FALSE) {
    res <- file.path(d, paste0(project, "_results"))
    gea <- file.path(res, "GEA", "tables")
    mth <- file.path(gea, "methods", "EMMAX")
    prc <- file.path(res, "Processing", "tables")
    for (p in c(mth, prc)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

    write_tsv(data.table::data.table(
        step = c(rep("processing", 7L), rep("gea", 4L)),
        metric = c("samples_total", "samples_removed", "samples_het_outliers_removed",
                   "samples_removed_relatedness", "samples_after_filtering",
                   "samples_with_coordinates", "samples_dropped_missing_coordinates",
                   "selected_snps_total", "regions_per_trait", "regions_combined",
                   "genes_found"),
        # 10 - 2 - 0 - 0 = 8, and 8 + 0 = 8: both accounting identities close.
        value = c("10", "2", "0", "0", "8", "8", "0", "3", "2", "2", "2")),
        file.path(res, "pipeline_summary.tsv"))

    write_tsv(data.table::data.table(
        stage = c("raw", "maf", "missingness"),
        n_samples = c(10L, 10L, 8L),
        n_snps    = c(100L, 80L, 70L)),
        file.path(prc, "filtering_summary.tsv"))

    write_tsv(data.table::data.table(
        site   = rep(c("NEG", "TAV"), each = 4L),
        sample = sprintf("ID%03d", seq_len(8L))),
        file.path(prc, "metadata.tsv"))

    snpid <- c("1:100", "1:200", "2:100")
    minp  <- c(1e-8, 2e-8, 3e-8)

    # EMMAX holds trait names, not numbers. min_pvalue is the only numeric column.
    write_tsv(data.table::data.table(
        SNPID = snpid, chr = c("1", "1", "2"), pos = c(100L, 200L, 100L),
        EMMAX = c("bio_1", "bio_1", "bio_2"), min_pvalue = minp),
        file.path(gea, "selected_snps.tsv"))

    write_tsv(data.table::data.table(
        region_id = c("1_50-250_bio_1", "2_50-150_bio_2"),
        trait = c("bio_1", "bio_2"), chr = c("1", "2"),
        start = c(50L, 50L), end = c(250L, 150L), length = c(200L, 100L),
        snp_count = c(2L, 1L), snp_ids = c("1:100,1:200", "2:100")),
        file.path(gea, "regions_per_trait.tsv"))

    write_tsv(data.table::data.table(
        region_id = c("1_50-250", "2_50-150"), chr = c("1", "2"),
        start = c(50L, 50L), end = c(250L, 150L), length = c(200L, 100L),
        snp_count = c(2L, 1L), snp_ids = c("1:100,1:200", "2:100")),
        file.path(gea, "regions_combined.tsv"))

    genes <- data.table::data.table(
        region_id = c("1_50-250_bio_1", "2_50-150_bio_2"),
        gene_id = c("g1", "g2"), chr = c("1", "2"),
        gene_start = c(60L, 60L), gene_end = c(120L, 120L))
    if (violations) {
        # A gene table chr is checked by check_chromosome_names but takes part in
        # no membership check, so the prefix mutation stays contained.
        genes <- rbind(genes, data.table::data.table(
            region_id = "2_50-150_bio_2", gene_id = "g3", chr = "chr2",
            gene_start = 70L, gene_end = 130L))
    }
    write_tsv(genes, file.path(gea, "genes_per_region.tsv"))
    write_tsv(data.table::data.table(
        region_id = c("1_50-250_bio_1", "2_50-150_bio_2"),
        gene_ids = c("g1", "g2")),
        file.path(gea, "genes_per_region_collapsed.tsv"))
    write_tsv(data.table::data.table(
        gene_id = c("g1", "g2"), chr = c("1", "2")),
        file.path(gea, "genes_combined.tsv"))

    # The canonical chromosome set is derived from THIS table alone, so it must
    # carry every chr the region/gene tables use.
    pv <- data.table::data.table(
        SNPID = snpid, chr = c("1", "1", "2"), pos = c(100L, 200L, 100L),
        bio_1 = c(1e-8, 2e-8, 0.4), bio_2 = c(0.5, 0.6, 3e-8))
    if (violations) pv[1L, bio_1 := 1.5]
    write_tsv(pv, file.path(mth, "EMMAX_pvalues_K3.tsv"))

    # One threshold variant only: a second would fire
    # check_single_threshold_variant, which is a real staleness defect and not
    # what either row is measuring.
    write_tsv(data.table::data.table(
        SNPID = snpid, chr = c("1", "1", "2"), pos = c(100L, 200L, 100L),
        pvalue = minp, trait = c("bio_1", "bio_1", "bio_2"),
        overlap_traits = c("bio_2", "", "bio_1"),
        overlap_snps   = c("2:100", "", "1:100")),
        file.path(mth, "EMMAX_pvalues_K3_sig_snps_bonf_0.05.tsv"))

    res
}

# ------------------------------------------------------------------- runner

# Invoked exactly as the Snakefile does: `Rscript /pipeline/scripts/<name>.R ...`.
# The scripts are not executable (mode 644) and carry no reliable shebang, so
# executing the path directly gives "Permission denied" rather than running them.
run_wrapper <- function(script, args, wd = NULL) {
    # shQuote every argument. system2() builds a SHELL command line, so an
    # unquoted space-separated file list (combine_selected_snps.R's argv[1],
    # combine_pheno_pvalues.R's argv[1]) would be split into separate argv
    # entries and every later positional would shift by one — the script then
    # writes its output over one of its own inputs. The Snakefile quotes these
    # the same way ("{params.files_str}", gwas.smk:93, _assoc_downstream.smk:81).
    invoke <- function() {
        suppressWarnings(system2("Rscript",
                                 shQuote(c(file.path(SCRIPTS, script), args)),
                                 stdout = TRUE, stderr = TRUE))
    }
    out <- if (is.null(wd)) invoke() else withr::with_dir(wd, invoke())
    status <- attr(out, "status")
    list(status = if (is.null(status)) 0L else as.integer(status),
         output = paste(out, collapse = "\n"))
}

# A spec is a list with:
#   label    test-name suffix (two specs may share one script)
#   script   basename inside scripts/
#   build(d) writes fixtures into d, returns whatever args() needs
#   args(d, built)     character argv
#   outputs(d)         files asserted to exist and be non-empty
#   out_min(d)         OPTIONAL list(dir=, pattern=, n=) — assert at least n
#                      files match pattern in dir. For scripts that derive their
#                      own basenames (plot_manhattan.R:156-208): re-deriving them
#                      in the test would duplicate production logic.
#   tools    OPTIONAL character vector of executables/paths the heavy tier
#            requires; asserted present, never skipped (see require_tools()).
#   expect_status OPTIONAL integer exit status, default 0L. Only one wrapper is
#            not a "runs clean" smoke: check_invariants.R is a VALIDATOR, and
#            quit(status = 1) on an error-severity violation (:272-274) is its
#            contract, not a failure. Asserting exit 1 on a tree built to
#            violate is what proves the file -> checker shaping layer — the
#            layer its own header (:12-18) names as untested — actually reaches
#            the checkers.
expect_wrapper_ok <- function(spec) {
    d <- withr::local_tempdir()
    sandbox <- file.path(d, "cwd")
    dir.create(sandbox, recursive = TRUE, showWarnings = FALSE)

    args <- spec$args(d, spec$build(d))
    res  <- run_wrapper(spec$script, args, wd = sandbox)

    want_status <- if (is.null(spec$expect_status)) 0L else as.integer(spec$expect_status)
    testthat::expect_identical(
        res$status, want_status,
        info = paste0(spec$script, " exited ", res$status, ", expected ", want_status,
                      "\nargs: ", paste(args, collapse = " "), "\n", res$output))

    for (f in spec$outputs(d)) {
        testthat::expect_true(file.exists(f),
            info = paste0(spec$script, " declared ", f, " but did not write it\n",
                          res$output))
        testthat::expect_gt(file.size(f), 0)
    }

    if (!is.null(spec$out_min)) {
        om    <- spec$out_min(d)
        found <- list.files(om$dir, pattern = om$pattern)
        # expect_true, not expect_gte: testthat 3.2.3's expect_gte takes no
        # `info` argument, and the message is the whole point of this assertion.
        testthat::expect_true(length(found) >= om$n,
            info = paste0(spec$script, ": expected >= ", om$n, " files matching ",
                          om$pattern, " in ", om$dir, ", found ",
                          length(found), "\n", res$output))
    }

    # The cwd sandbox must be untouched: anything here is a path the script
    # resolved relative to its working directory, which in production is the
    # mounted repo.
    leaked <- list.files(sandbox, all.files = TRUE, no.. = TRUE, recursive = TRUE)
    testthat::expect_identical(
        leaked, character(0),
        info = paste0(spec$script, " wrote into its working directory: ",
                      paste(leaked, collapse = ", ")))

    invisible(d)
}

# Negative control: hand the script a first input that does not exist and require
# a non-zero exit with no output file written. This is the assertion that catches
# "silently wrote a stub and exited 0", and it is opt-in per spec because a few
# wrappers legitimately treat a declared input as optional — write_summary.R's
# read_opt() (:199,:261,:317) returns NULL for a missing file and skips the
# guarded block, which is itself filed as a defect.
expect_wrapper_fails_without_input <- function(spec) {
    d <- withr::local_tempdir()
    sandbox <- file.path(d, "cwd")
    dir.create(sandbox, recursive = TRUE, showWarnings = FALSE)

    args <- spec$args(d, spec$build(d))
    args[1] <- file.path(d, "definitely-not-here.tsv")
    outs <- spec$outputs(d)
    unlink(outs)

    res <- run_wrapper(spec$script, args, wd = sandbox)
    testthat::expect_gt(res$status, 0L)
    for (f in outs) {
        testthat::expect_false(file.exists(f),
            info = paste0(spec$script, " wrote ", f, " despite a missing input\n",
                          res$output))
    }
}

# Heavy tier only. An absent tool inside cline-go:latest means the IMAGE is
# wrong, which must be red — never a skip. A skip here would be the exact
# failure mode the quick/heavy split exists to avoid.
require_tools <- function(tools) {
    for (t in tools) {
        ok <- if (grepl("/", t)) file.exists(t) else nzchar(Sys.which(t))
        testthat::expect_true(ok, info = paste0(
            "required tool missing from the image: ", t,
            " — heavy tests fail rather than skip on this, see tests/heavy/"))
    }
}
