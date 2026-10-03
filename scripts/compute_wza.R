#!/usr/bin/env Rscript
# Weighted-Z Analysis (WZA) for GEA/GWAS p-values
#
# Implements Booker et al. (2024) Mol Ecol Resour 24:e13768
# DOI: 10.1111/1755-0998.13768
#
# Aggregates per-SNP p-values into fixed-width genomic windows using
# heterozygosity-weighted Z-scores, then corrects for window SNP-count
# variation via polynomial regression.
#
# Steps (Booker et al. §2):
#   1. Rank-transform per-SNP p-values → empirical p (genome-wide)
#   2. Convert empirical p → one-sided z: z_i = qnorm(1 - empirical_p_i)
#   3. Per window k: Z_W,k = Σ(p̄q̄·z_i) / sqrt(Σ((p̄q̄)²))
#      where p̄q̄ = 2·MAF·(1−MAF) is expected heterozygosity
#   4. Correct for SNP-count variation: fit linear model mean(Z_W) ~ n_snps
#      and sd(Z_W) ~ n_snps using a sliding bin of up to 50 windows.
#      Per-window WZA p = 1 - pnorm(Z_W, predicted_mean, predicted_sd)
#
# Output: TSV in same wide format as pvalues input (SNPID, chr, pos, n_snps,
#         mean_maf, <trait_1>, ..., <trait_n>)
# SNPID = "{chr}:{start}-{end}", pos = window center

library(data.table)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
options(scipen = 99999)
################
PVALUES_TSV    = args[1]  # per-SNP pvalues (wide: SNPID, chr, pos, trait_cols)
MAF_TSV        = args[2]  # maf_pos.tsv: CHR, POS, MAF
LD_DECAY_PATH  = if (length(args) >= 3) args[3] else "NULL"
WINDOW_SIZE    = if (length(args) >= 4) args[4] else "auto_genome_wide"
FALLBACK_BP    = if (length(args) >= 5) as.integer(args[5]) else 10000L
LD_DECAY_GROUP = if (length(args) >= 6) args[6] else "All"
OUTPUT         = args[7]
################

# No WZA-specific MAF cut. The SNP set is whatever `Filter.maf` already admitted — WZA does
# not get to apply a second, invisible filter on top of the user's.
#
# This used to be `MIN_MAF_WZA <- 0.05` (Booker et al. §3), applied before windowing. The
# recommendation is sound in its original setting, but hardcoding it here made the WZA arm
# silently disagree with every other output of the same run: on the MVP benchmark at
# `Filter.maf: 0.01` it dropped 2 502 of 10 325 SNPs and 28 of 70 causal loci, so WZA recall
# carried a 0.60 ceiling that looked like a method property and was actually a hidden filter.
# A user who wants the Booker cut sets `Filter.maf: 0.05` and gets it everywhere, consistently.
#
# SNPs with no MAF record are still dropped: the weight h = 2·MAF·(1−MAF) is undefined
# without one, so they cannot enter a window sum at all.
SLIDING_BIN   <- 50L    # windows per sliding bin for n_snps correction

message("INFO: WZA — reading inputs")

# ─── Resolve window size ──────────────────────────────────────────────────────

hill_weir <- function(d, C, n) {
    x <- C * d
    ((10 + x) / ((2 + x) * (11 + x))) *
    (1 + ((3 + x) * (12 + 12 * x + x^2)) / (n * (2 + x) * (11 + x)))
}

invert_hill_weir <- function(C_hat, n_samples, r2_target, max_dist_bp) {
    if (is.na(C_hat) || is.na(n_samples) || is.na(max_dist_bp)) return(NA_real_)
    tryCatch(
        uniroot(function(d) hill_weir(d, C_hat, n_samples) - r2_target,
                c(1, max_dist_bp))$root,
        error = function(e) NA_real_
    )
}

resolve_window_size <- function(spec, ld_path, group, fallback_bp) {
    spec_lower <- tolower(trimws(spec))
    if (spec_lower %in% c("auto_genome_wide", "auto")) {
        if (ld_path == "NULL" || !file.exists(ld_path)) {
            message(paste0("WARNING: window_size='", spec, "' but LD decay table not found. ",
                           "Using fallback: ", fallback_bp, " bp"))
            return(fallback_bp)
        }
        ld <- fread(ld_path)
        row <- ld[group == LD_DECAY_GROUP & scope == "genome_wide"]
        if (nrow(row) == 0) {
            message(paste0("WARNING: No genome-wide LD decay row for group='", LD_DECAY_GROUP,
                           "'. Using fallback: ", fallback_bp, " bp"))
            return(fallback_bp)
        }
        row <- as.list(row[1])
        if (row$method == "hill_weir") {
            d <- invert_hill_weir(as.numeric(row$C_hat), as.numeric(row$n_samples),
                                  0.2, as.numeric(row$max_dist_bp))
        } else {
            d <- as.numeric(row$r2_02_bp)
        }
        if (is.na(d)) {
            message(paste0("WARNING: LD-decay inversion failed. Using fallback: ", fallback_bp, " bp"))
            return(fallback_bp)
        }
        d <- max(round(d), fallback_bp)  # never smaller than fallback
        message(paste0("INFO: WZA window size from LD decay (r²=0.2, group=", LD_DECAY_GROUP,
                       "): ", d, " bp"))
        return(as.integer(d))
    } else if (spec_lower == "auto_per_chromosome") {
        # Use genome-wide as approximation for WZA (per-chr adds complexity for little gain)
        message("INFO: WZA window_size=auto_per_chromosome treated as auto_genome_wide")
        return(resolve_window_size("auto_genome_wide", ld_path, group, fallback_bp))
    } else {
        d <- as.integer(spec)
        message(paste0("INFO: WZA fixed window size: ", d, " bp"))
        return(d)
    }
}

window_bp <- resolve_window_size(WINDOW_SIZE, LD_DECAY_PATH, LD_DECAY_GROUP, FALLBACK_BP)

# ─── Load inputs ─────────────────────────────────────────────────────────────

pval_dt <- fread(PVALUES_TSV)
pval_dt$chr <- as.character(pval_dt$chr)
pval_dt$pos <- as.integer(pval_dt$pos)
trait_cols <- setdiff(colnames(pval_dt), c("SNPID", "chr", "pos"))
message(paste0("INFO: ", nrow(pval_dt), " SNPs, ", length(trait_cols), " traits"))

if (length(trait_cols) == 0) {
    stop("No trait columns found in pvalues TSV")
}

maf_dt <- fread(MAF_TSV, colClasses = c("CHR" = "character"))
colnames(maf_dt) <- tolower(colnames(maf_dt))
# normalise chr column name
if ("chr" %in% colnames(maf_dt) && !"pos" %in% colnames(maf_dt)) {
    # maf_pos.tsv is CHR  POS  MAF
    setnames(maf_dt, old = colnames(maf_dt)[1:3], new = c("chr", "pos", "maf"))
}
maf_dt$chr <- as.character(maf_dt$chr)
maf_dt$pos <- as.integer(maf_dt$pos)

# Join MAF into pval table
pval_dt <- merge(pval_dt, maf_dt[, .(chr, pos, maf)], by = c("chr", "pos"), all.x = TRUE)

# Drop only SNPs with no usable MAF — see the note at MIN_MAF_WZA's removal above. The MAF
# threshold itself belongs to Filter.maf and has already been applied upstream.
n_before <- nrow(pval_dt)
pval_dt  <- pval_dt[!is.na(maf) & maf > 0]
message(paste0("INFO: SNPs with a usable MAF: ", nrow(pval_dt), "/", n_before,
               " (no WZA-specific MAF cut; Filter.maf governs)"))

if (nrow(pval_dt) == 0) {
    message("WARNING: No SNPs carry a usable MAF. Writing empty WZA output.")
    empty <- data.table(SNPID = character(), chr = character(), pos = integer(),
                        n_snps = integer(), mean_maf = numeric())
    for (tr in trait_cols) empty[[tr]] <- numeric()
    fwrite(empty, OUTPUT, sep = "\t")
    quit(save = "no", status = 0)
}

# ─── Assign windows ──────────────────────────────────────────────────────────

pval_dt[, window_start := (floor((pos - 1L) / window_bp)) * window_bp + 1L]
pval_dt[, window_end   := window_start + window_bp - 1L]
pval_dt[, window_id    := paste0(chr, ":", window_start, "-", window_end)]

# ─── WZA per trait ───────────────────────────────────────────────────────────

# Rank-transform p-values genome-wide → empirical p ∈ (0,1)
# Use mid-rank to handle ties; exclude NAs
empirical_p <- function(p) {
    is_valid <- !is.na(p)
    n <- sum(is_valid)
    ep <- rep(NA_real_, length(p))
    ep[is_valid] <- rank(p[is_valid], ties.method = "average") / (n + 1)
    ep
}

# Weighted-Z per window: Z_W = Σ(h·z) / sqrt(Σ(h²))  where h = 2·maf·(1−maf)
wza_window <- function(h, z) {
    valid <- !is.na(h) & !is.na(z) & h > 0
    if (sum(valid) == 0) return(NA_real_)
    hv <- h[valid]; zv <- z[valid]
    sum(hv * zv) / sqrt(sum(hv^2))
}

# N-snp polynomial correction: predict mean+sd of Z_W from n_snps in window,
# using a sliding bin of SLIDING_BIN windows ordered by n_snps.
poly_correction <- function(n_snps_vec, zw_vec) {
    ord <- order(n_snps_vec)
    ns <- n_snps_vec[ord]
    zw <- zw_vec[ord]
    N  <- length(ns)

    bin_mean <- numeric(N)
    bin_sd   <- numeric(N)
    half <- floor(SLIDING_BIN / 2L)
    for (i in seq_len(N)) {
        idx <- max(1L, i - half):min(N, i + half)
        bin_mean[i] <- mean(zw[idx], na.rm = TRUE)
        bin_sd[i]   <- sd(zw[idx],   na.rm = TRUE)
    }
    bin_sd[is.na(bin_sd) | bin_sd == 0] <- sd(zw, na.rm = TRUE)  # global fallback

    if (N >= 3) {
        fit_mean <- suppressWarnings(lm(bin_mean ~ poly(ns, min(2L, N - 1L), raw = TRUE)))
        fit_sd   <- suppressWarnings(lm(bin_sd   ~ poly(ns, min(2L, N - 1L), raw = TRUE)))
        pred_mean <- predict(fit_mean, newdata = data.frame(ns = ns))
        pred_sd   <- pmax(predict(fit_sd, newdata = data.frame(ns = ns)), 1e-6)
    } else {
        pred_mean <- rep(mean(zw, na.rm = TRUE), N)
        pred_sd   <- rep(max(sd(zw, na.rm = TRUE), 1e-6), N)
    }

    result        <- data.table(ord_idx = ord, pred_mean = pred_mean, pred_sd = pred_sd)
    result[order(ord_idx)]
}

# How small a WZA p CAN this trait produce, and is that small enough to ever be
# called? The rank transform maps the best per-SNP p to 1/(m+1), so the per-SNP z is
# capped at qnorm(1 - 1/(m+1)) — 2.768 at m = 354. A window holding a single SNP has
# Z_W = h*z/sqrt(h^2) = z, so its WZA p cannot go below the capped z's p under the
# n_snps correction, however strong the underlying association. With most windows
# single-SNP that makes Bonferroni unreachable BY CONSTRUCTION, and a header-only
# sig-windows file used to look identical to a genuine absence of signal (finding
# 6ab65e). This reports the floor so the two are distinguishable.
#
# Nothing here changes a computed p-value, and no new output file is written: the
# finding's own evidence came from the rule log, and a diagnostics TSV nobody reads
# is its own defect (cf. rda_offset_diagnostics.tsv, finding 0214e1).
report_wza_floor <- function(trait, p_vec, win_zw, pred_mean, pred_sd) {
    m <- sum(!is.na(p_vec))
    if (m == 0L) {
        message(paste0("WARNING: WZA trait=", trait, " has no non-NA p-values"))
        return(invisible(NULL))
    }
    z_cap <- qnorm(1 - 1 / (m + 1))

    single <- which(win_zw$n_snps == 1L & !is.na(pred_mean) & !is.na(pred_sd))
    n_win  <- nrow(win_zw)
    floor_single <- if (length(single) > 0) {
        min(pnorm(z_cap, pred_mean[single], pred_sd[single], lower.tail = FALSE))
    } else NA_real_

    bonf <- 0.05 / sum(!is.na(win_zw$wza_p))
    message(paste0("INFO: WZA trait=", trait,
                   " n_snps_ranked=", m,
                   " z_cap=", signif(z_cap, 4),
                   " single_snp_windows=", length(single), "/", n_win,
                   " min_achievable_p_single_snp=", signif(floor_single, 4),
                   " observed_min_p=", signif(min(win_zw$wza_p, na.rm = TRUE), 4),
                   " bonf_0.05=", signif(bonf, 4)))

    if (!is.na(floor_single) && floor_single > bonf) {
        message(paste0("WARNING: WZA trait=", trait,
                       " — no single-SNP window can reach Bonferroni 0.05/",
                       sum(!is.na(win_zw$wza_p)), " = ", signif(bonf, 4),
                       ": the rank transform caps the per-SNP z at ", signif(z_cap, 4),
                       " (m=", m, " SNPs), so the smallest p such a window can take is ",
                       signif(floor_single, 4), ". ", length(single), " of ", n_win,
                       " windows hold one SNP, so an empty sig-windows table here means ",
                       "UNCALLABLE BY CONSTRUCTION, not absence of signal. Widen the ",
                       "window (GEA.WZA.window_size) or raise the threshold."))
    }
    invisible(NULL)
}

# Compute per-window WZA for all traits
windows_base <- pval_dt[, .(
    chr       = chr[1],
    start     = window_start[1],
    end       = window_end[1],
    pos       = as.integer((window_start[1] + window_end[1]) / 2L),
    n_snps    = .N,
    mean_maf  = mean(maf, na.rm = TRUE)
), by = window_id]
setorder(windows_base, chr, start)

wza_results <- copy(windows_base)

for (tr in trait_cols) {
    message(paste0("INFO:   WZA for trait: ", tr))

    # empirical p genome-wide
    pval_dt[, emp_p := empirical_p(get(tr))]
    pval_dt[, z_i   := qnorm(1 - emp_p)]

    # heterozygosity weight h = 2·MAF·(1−MAF)
    pval_dt[, maf_h := 2 * maf * (1 - maf)]

    # per-window Z_W
    win_zw <- pval_dt[, .(Z_W = wza_window(maf_h, z_i)), by = window_id]

    # merge window metadata for correction
    win_zw <- merge(win_zw, windows_base[, .(window_id, n_snps)], by = "window_id")
    valid  <- !is.na(win_zw$Z_W)

    pred_mean <- rep(NA_real_, nrow(win_zw))
    pred_sd   <- rep(NA_real_, nrow(win_zw))

    if (sum(valid) >= 2) {
        corr <- poly_correction(win_zw$n_snps[valid], win_zw$Z_W[valid])
        pred_mean[valid] <- corr$pred_mean
        pred_sd[valid]   <- corr$pred_sd
        wza_p <- rep(NA_real_, nrow(win_zw))
        wza_p[valid] <- pnorm(win_zw$Z_W[valid], corr$pred_mean, corr$pred_sd, lower.tail = FALSE)
    } else {
        message(paste0("WARNING: Too few windows with Z_W for trait '", tr,
                       "' (n=", sum(valid), "). Reporting uncorrected p-values."))
        global_sd <- max(sd(win_zw$Z_W, na.rm = TRUE), 1e-6)
        wza_p     <- pnorm(win_zw$Z_W, mean(win_zw$Z_W, na.rm = TRUE), global_sd, lower.tail = FALSE)
    }

    # Upper-tail pnorm() avoids the catastrophic cancellation of `1 - pnorm()`, which
    # loses precision for small p and rounds to exactly 0 for very large Z_W. Floor at
    # the smallest positive double so a window p is never exactly 0 — otherwise
    # -log10(p) becomes +Inf downstream and crashes scattermore in the Manhattan
    # background. NA windows (no valid SNPs for the trait) are preserved as NA.
    wza_p <- pmax(wza_p, .Machine$double.xmin)

    win_zw[, wza_p := wza_p]

    report_wza_floor(tr, pval_dt[[tr]], win_zw, pred_mean, pred_sd)

    # merge into results keyed by window_id
    wza_results <- merge(
        wza_results,
        win_zw[, .(window_id, wza_p = wza_p)],
        by = "window_id", all.x = TRUE
    )
    setnames(wza_results, "wza_p", tr)
}

# ─── Format output ───────────────────────────────────────────────────────────

# Output mirrors pvalues TSV format: SNPID, chr, pos, n_snps, mean_maf, <traits>
setnames(wza_results, "window_id", "SNPID")
col_order <- c("SNPID", "chr", "pos", "n_snps", "mean_maf", trait_cols)
wza_results <- wza_results[, ..col_order]
setorder(wza_results, chr, pos)

fwrite(wza_results, OUTPUT, sep = "\t")

# Diagnostic summary
n_windows_total <- nrow(wza_results)
for (tr in trait_cols) {
    n_non_na <- sum(!is.na(wza_results[[tr]]))
    message(paste0("INFO: WZA trait=", tr,
                   " windows_with_p=", n_non_na,
                   " windows_na=", n_windows_total - n_non_na))
}
message(paste0("INFO: WZA complete — n_windows=", n_windows_total,
               " n_dropped_no_maf=", n_before - nrow(pval_dt),
               " output=", OUTPUT))
