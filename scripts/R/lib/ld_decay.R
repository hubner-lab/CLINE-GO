# =============================================================================
# ld_decay.R -- LD-decay distances from per-distance r^2 sums.
#
# r^2 is PLINK 1.9's --r2, the squared correlation of allele counts with missing
# calls excluded pairwise, for both the decay curve (pairs on one chromosome) and its
# background (pairs on different chromosomes), so the two come from one estimator.
# ld_pairs_to_stat() sums PLINK's pair table per bp distance, in PopLDdecay's .stat.gz
# layout, which ld_decay_analyze.R reads.
#
# The curve is the pair-weighted mean r^2 by distance, smoothed by LOESS (degree 2)
# on log10(distance) over bins that are fine near zero (10 bp below 100 bp) and
# log-spaced beyond. Three decay distances are read off it, each the first distance
# at which the running minimum of the smoothed curve falls to a target:
#   critical -- the background LD of unlinked markers: the 95th percentile of
#               sqrt(r^2) between markers on different chromosomes, squared
#               (Breseghello & Sorrells 2006, as applied by Laido et al. 2014);
#   half     -- half of the short-range LD (Zhao et al. 2011; Vos et al. 2017);
#   r2thr    -- a fixed r^2 level.
#
# Pure functions (data.tables and numbers in, numbers out), plus two file helpers,
# ld_write_stat() and ld_scratch_dir(). Sourced by scripts/ld_decay_prepare.R,
# scripts/ld_decay_pairs.R, scripts/ld_decay_background.R, scripts/ld_decay_analyze.R
# and tests/testthat/helper-libs.R.
# =============================================================================

# Bin breaks: fine bins up to `fine_until`, then `per_decade` log-spaced bins up to
# `max_dist`. Returned ascending and unique; the last break is >= max_dist.
ld_log_breaks <- function(max_dist, fine_bp = 10, fine_until = 100, per_decade = 20) {
    stopifnot(is.numeric(max_dist), max_dist > 0)
    fine <- seq(0, min(fine_until, max_dist), by = fine_bp)
    if (max_dist <= fine_until) return(unique(c(fine, max_dist)))
    n_log <- max(2L, ceiling(per_decade * log10(max_dist / fine_until)) + 1L)
    logb  <- 10^seq(log10(fine_until), log10(max_dist), length.out = n_log)
    sort(unique(round(c(fine, logb))))
}

# raw: data.table with dist_bp, sum_r2, n_pairs (distance 0 already removed). Returns
# one row per non-empty bin: the pair-weighted mean distance, the mean r^2 and the
# pair count, ordered by distance.
ld_bin_stat <- function(raw, breaks) {
    if (is.null(raw) || nrow(raw) == 0) return(NULL)
    dt <- data.table::as.data.table(raw)[dist_bp > 0 & n_pairs > 0]
    if (nrow(dt) == 0) return(NULL)
    dt[, bin := findInterval(dist_bp, breaks, left.open = TRUE)]
    out <- dt[, .(dist_bp = sum(dist_bp * n_pairs) / sum(n_pairs),
                  mean_r2 = sum(sum_r2) / sum(n_pairs),
                  n_pairs = sum(n_pairs)), by = bin]
    out[order(dist_bp)][, bin := NULL][]
}

# LOESS of mean r^2 on log10(distance), weighted by pair counts, predicted on a grid
# of `n_grid` log-spaced distances. r2_mono is the running minimum, so a first
# crossing is not undone by a small upward wiggle of the smooth.
ld_fit_curve <- function(binned, span = 0.3, n_grid = 1000L, min_bins = 8L) {
    if (is.null(binned) || nrow(binned) < min_bins) return(NULL)
    d <- data.frame(lx = log10(binned$dist_bp), y = binned$mean_r2, w = binned$n_pairs)
    # With few bins loess cannot estimate its residual scale and warns ("NaNs
    # produced"); the fitted values, which are all that is used here, are unaffected.
    fit <- tryCatch(
        suppressWarnings(stats::loess(y ~ lx, data = d, weights = w, span = span, degree = 2)),
        error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    lx <- seq(min(d$lx), max(d$lx), length.out = n_grid)
    y  <- stats::predict(fit, newdata = data.frame(lx = lx))
    if (anyNA(y)) return(NULL)
    y <- pmin(pmax(y, 0), 1)   # a sparse curve can overshoot the range of r^2
    data.table::data.table(dist_bp = 10^lx, r2 = y, r2_mono = cummin(y))
}

# First distance at which the running-minimum curve reaches `target`, interpolated
# linearly in log10(distance). NA when the target is NA, when the curve starts at or
# below it (nothing to decay to) or when it is not reached within the curve.
ld_first_crossing <- function(curve, target) {
    if (is.null(curve) || is.na(target)) return(NA_real_)
    y <- curve$r2_mono
    if (y[1] <= target) return(NA_real_)
    i <- which(y <= target)[1]
    if (is.na(i)) return(NA_real_)
    lx <- log10(curve$dist_bp[c(i - 1L, i)])
    f  <- (y[i - 1L] - target) / (y[i - 1L] - y[i])
    10^(lx[1] + f * (lx[2] - lx[1]))
}

# Critical r^2 of unlinked marker pairs: the `quantile` of sqrt(r^2), squared.
ld_critical_r2 <- function(r2, quantile = 0.95) {
    r2 <- r2[is.finite(r2)]
    if (length(r2) == 0) return(NA_real_)
    stats::quantile(sqrt(r2), probs = quantile, names = FALSE)^2
}

# All decay distances for one binned curve. short_r2 is the smoothed r^2 at the
# shortest binned distance; background_r2 comes from ld_critical_r2() (NA when the
# dataset has a single chromosome or no background table was given).
ld_decay_distances <- function(binned, background_r2 = NA_real_, r2_threshold = 0.2,
                               span = 0.3) {
    curve <- ld_fit_curve(binned, span = span)
    if (is.null(curve)) {
        return(list(method = "insufficient_data", short_r2 = NA_real_,
                    background_r2 = background_r2, critical_bp = NA_real_,
                    half_bp = NA_real_, r2thr_bp = NA_real_, curve = NULL))
    }
    s <- curve$r2[1]
    list(method = "loess_log", short_r2 = s, background_r2 = background_r2,
         critical_bp = ld_first_crossing(curve, background_r2),
         half_bp     = ld_first_crossing(curve, s / 2),
         r2thr_bp    = ld_first_crossing(curve, r2_threshold),
         curve = curve)
}

# Whether an LD-decay curve can give a clumping distance at all: "reliable" or
# "unreliable", with the reason. The test is the curve's room to decay, short_r2 over
# the critical r^2 of unlinked markers. Below `min_ratio` (default 2) the unlinked
# background lies above half the short-range LD, so the curve reaches noise before it
# has halved and every distance read off it is set by the background, not by decay.
# A small sample (unlinked r^2 rises roughly as 1/n) or a few strongly differentiated,
# inbred lineages compress the ratio. Also unreliable: no curve (too few bins) and no
# background (a single chromosome), where the ratio cannot be measured.
ld_reliability <- function(short_r2, background_r2, method, min_ratio = 2) {
    out <- function(ok, why = "") list(status = if (ok) "reliable" else "unreliable",
                                       reason = why)
    if (identical(method, "insufficient_data") || is.na(short_r2)) {
        return(out(FALSE, "too few distance bins to fit a decay curve"))
    }
    if (is.na(background_r2)) {
        return(out(FALSE, paste0("no critical r2 of unlinked markers (one chromosome), so ",
                                 "the curve's room to decay cannot be measured")))
    }
    ratio <- short_r2 / background_r2
    if (ratio < min_ratio) {
        return(out(FALSE, sprintf(paste0(
            "short-range r2 %.3g is %.2fx the critical r2 of unlinked markers (%.3g), below ",
            "%g: the unlinked background lies above half the short-range LD, so the curve ",
            "reaches noise before it halves (typical of a small sample or a few inbred, ",
            "strongly differentiated lineages)"), short_r2, ratio, background_r2, min_ratio)))
    }
    out(TRUE)
}

# How many anchor SNPs PLINK should pair with every SNP in the window, given a budget
# of pairs. chr/pos describe every SNP; max_dist is in bp. When all pairs within
# max_dist fit the budget the run is exhaustive (no anchors); otherwise anchors are
# drawn so that anchors x mean partners is about max_pairs.
ld_anchor_plan <- function(chr, pos, max_dist, max_pairs) {
    stopifnot(length(chr) == length(pos), max_dist > 0, max_pairs > 0)
    partners <- unlist(lapply(split(pos, factor(chr, unique(chr))), function(p) {
        p <- sort(p)
        hi <- findInterval(p + max_dist, p)
        lo <- findInterval(p - max_dist, p, left.open = TRUE) + 1L
        same <- findInterval(p, p) - findInterval(p, p, left.open = TRUE)
        (hi - lo + 1L) - same     # partners at distance 1..max_dist, either side
    }), use.names = FALSE)
    total <- sum(as.numeric(partners)) / 2
    mean_partners <- if (length(partners)) mean(partners) else 0
    if (total <= max_pairs || mean_partners == 0) {
        return(list(exhaustive = TRUE, n_anchors = length(pos), total_pairs = total,
                    mean_partners = mean_partners))
    }
    list(exhaustive = FALSE,
         n_anchors = as.integer(min(length(pos), ceiling(max_pairs / mean_partners))),
         total_pairs = total, mean_partners = mean_partners)
}

# Per-distance sums from a PLINK --r2 table (CHR_A, BP_A, SNP_A, BP_B, SNP_B, R2).
# Pairs at distance 0 or beyond max_dist and non-finite r^2 are dropped; with
# dedupe = TRUE an unordered pair reported twice (both SNPs anchors under
# --ld-snp-list) is counted once. Returns list(gw, chr): data.tables of
# dist_bp, sum_r2, n_pairs, the second also keyed by chr.
ld_pairs_to_stat <- function(pairs, max_dist, dedupe = FALSE) {
    p <- data.table::as.data.table(pairs)
    p[, R2 := suppressWarnings(as.numeric(R2))]   # PLINK writes "nan" for a monomorphic SNP
    p <- p[is.finite(R2)]
    p[, dist_bp := abs(BP_B - BP_A)]
    p <- p[dist_bp > 0 & dist_bp <= max_dist]
    if (dedupe && nrow(p)) {
        p[, key := ifelse(SNP_A < SNP_B, paste(SNP_A, SNP_B), paste(SNP_B, SNP_A))]
        p <- unique(p, by = "key")
    }
    gw  <- p[, .(sum_r2 = sum(R2), n_pairs = .N), keyby = dist_bp]
    chr <- p[, .(sum_r2 = sum(R2), n_pairs = .N), keyby = .(chr = as.character(CHR_A), dist_bp)]
    list(gw = gw, chr = chr)
}

# Seeded markers for the unlinked background: up to `per_chr` on each chromosome and
# never more than `max_markers` in all, since unlinked pairs grow with the square of
# the marker count. A fragmented assembly gets fewer markers per contig, down to one,
# and beyond `max_markers` contigs a seeded sample of those. Returns the markers
# (chr, pos) and the number per chromosome that was allowed.
ld_background_markers <- function(chr, pos, per_chr, max_markers, seed) {
    dt <- data.table::data.table(chr = as.character(chr), pos = pos)
    k <- min(per_chr, max(1L, max_markers %/% data.table::uniqueN(dt$chr)))
    set.seed(seed)
    picked <- dt[, .SD[sort(sample.int(.N, min(.N, k)))], by = chr]
    if (nrow(picked) > max_markers) picked <- picked[sort(sample.int(.N, max_markers))]
    list(markers = picked, per_chr = k)
}

# The contigs among `chrs` (VCF names) that `exclude` names: LDdecay.exclude_chromosomes
# as one comma-separated word ("X,Y,MT,Pt"; "NULL" or "" for none), for contigs that do
# not recombine like autosomes. Case and a leading "chr" are ignored and M is MT, as
# PLINK reads them, so "chrM" names the filtered VCF's MT. A listed name with no contig
# in `chrs` is reported, so a typo shows in the log.
ld_excluded_contigs <- function(chrs, exclude) {
    ex <- trimws(strsplit(as.character(exclude), ",", fixed = TRUE)[[1]])
    ex <- ex[nzchar(ex) & ex != "NULL"]
    if (!length(ex)) return(character(0))
    canon <- function(x) {
        x <- sub("^chr", "", tolower(x))
        ifelse(x == "m", "mt", x)
    }
    chrs <- unique(as.character(chrs))
    absent <- ex[!canon(ex) %in% canon(chrs)]
    if (length(absent)) {
        message("INFO: LDdecay.exclude_chromosomes names no contig of the VCF: ",
                paste(absent, collapse = ", "))
    }
    chrs[canon(chrs) %in% canon(ex)]
}

# Write per-distance sums in PopLDdecay's .stat.gz layout (header row, six columns:
# #Dist Mean_r^2 Mean_D' Sum_r^2 Sum_D' NumberPairs; D' is not computed, so NA).
ld_write_stat <- function(stat, path) {
    out <- data.table::data.table(stat$dist_bp, stat$sum_r2 / stat$n_pairs, NA_real_,
                                  stat$sum_r2, NA_real_, stat$n_pairs)
    data.table::setnames(out, c("#Dist", "Mean_r^2", "Mean_D'", "Sum_r^2", "Sum_D'",
                                "NumberPairs"))
    data.table::fwrite(out, path, sep = "\t", na = "NA", compress = "gzip")
    invisible(path)
}

# Scratch directory `name` beside `path` (a file or a folder), in the LD-decay work
# folder on the results volume. In a container tempdir() sits on the host's root disk,
# which 48 parallel runs filled (2026-10-26), so the LD-decay scripts keep every
# temporary file here, fread(cmd = )'s copies too (tmpdir =). A leftover of a failed
# run is removed first. The caller removes the directory when done: a top-level
# on.exit() never runs under Rscript.
ld_scratch_dir <- function(path, name) {
    d <- file.path(dirname(normalizePath(path, mustWork = FALSE)), name)
    unlink(d, recursive = TRUE)
    dir.create(d, recursive = TRUE)
    d
}
