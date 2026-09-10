#!/usr/bin/env Rscript
# =============================================================================
# mvp_wza_gate.R -- is window-based calling (WZA) viable on this benchmark, and if so,
# how does it score?
#
# WHY A GATE AND NOT JUST A SCORE
# -------------------------------
# MVP_README.md:102 declared WZA unusable on this dataset: "~10k SNPs over 20 LGs gives
# order 10^1-10^2 windows, an order of magnitude short of what window statistics need."
# That was measured under the pipeline default `window_size: auto_genome_wide`. It is a
# statement about the WINDOW SIZE, not about the dataset: the simulated genome is 1 Mb,
# and the same rescaling already applied to LD.window (100 -> 10 kb) applies here. At a
# fixed 1 kb window the genome yields ~1000 windows at ~8-10 SNPs each.
#
# So the viability question has to be re-asked with numbers rather than inherited. The
# binding constraint is NOT the window count -- it is how many windows actually contain a
# causal locus that WZA can see, because WZA can only show signal in those. This script
# computes that, applies the gate, and refuses to score if the gate fails.
#
# THE MAF SUBTLETY -- a pipeline defect this script surfaced, now fixed
# ---------------------------------------------------------------------
# scripts/compute_wza.R used to hardcode MIN_MAF_WZA <- 0.05 (Booker et al. sec.3) and drop
# every SNP below it BEFORE windowing, on top of whatever Filter.maf had already done. On
# MVP1232548 at Filter.maf 0.01 that silently removed 2 502 of 10 325 SNPs and 28 of the 70
# causal loci, so WZA carried a recall ceiling of 0.60 that read as a method property and was
# actually a hidden second filter. Removed 2026-08-03: WZA now windows exactly the SNP set
# Filter.maf admitted, like every other output of the same run.
#
# --wza-maf therefore defaults to 0 (no extra cut) and exists only to reproduce the old
# behaviour for comparison. Whenever it is non-zero, the causal denominator is restricted to
# match, and both counts are emitted so the gap stays visible.
#
# WINDOW-LEVEL SCORING CONTRACT -- the SNP-level conventions, lifted one unit up:
#   TP window  a called window containing >= 1 testable causal locus
#   FP window  a called window on a BACKGROUND linkage group (11-20) -- the same rule as
#              precision_strict, where only background_neutral counts as a false positive
#   ignored    a called window on LG 1-10 with no causal locus inside: that is the window
#              analogue of linked_neutral and is not charged as an error
#   recall     testable causal loci covered by called windows / testable causal loci
#              (LOCUS-level, so it stays comparable to the SNP-level recall_testable)
#
# Usage:
#   Rscript mvp_wza_gate.R --wza="EMMAX=/p/a.tsv,LFMM=/p/b.tsv" --truth=FILE \
#       --loci-freq=FILE --tag=STR [--outdir=DIR] [--min-causal-windows=30]
#       [--wza-maf=0] [--bg-lg-min=11]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(qvalue)
    library(ggplot2)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/pval_threshold.R"))
source(file.path(PIPELINE_ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args <- parse_kv_args(commandArgs(trailingOnly = TRUE))
req  <- function(k) {
    if (is.null(args[[k]]) || !nzchar(args[[k]])) stop("Missing required --", gsub("_", "-", k))
    args[[k]]
}
opt  <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]

TAG        <- req("tag")
TRUTH_F    <- req("truth")
FREQ_F     <- req("loci_freq")
OUTDIR     <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/wza"))
MIN_CW     <- as.integer(opt("min_causal_windows", "30"))
WZA_MAF    <- as.numeric(opt("wza_maf", "0"))
BG_LG_MIN  <- as.integer(opt("bg_lg_min", "11"))
TOP_MAX    <- as.numeric(opt("top_max", "100"))

spec <- strsplit(req("wza"), ",", fixed = TRUE)[[1]]
tabs <- list()
for (s in spec) {
    kv <- strsplit(s, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Bad --wza entry (expected NAME=PATH): ", s)
    if (!file.exists(kv[2])) stop("WZA table not found: ", kv[2])
    tabs[[kv[1]]] <- kv[2]
}
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

truth  <- load_truth(TRUTH_F)
freq   <- fread(FREQ_F, colClasses = c(chr = "character"))
freq[, key := paste(chr, pos, sep = ":")]
causal <- merge(truth$causal[, .(key, chr, pos)], freq[, .(key, maf)], by = "key", all.x = TRUE)
if (anyNA(causal$maf)) stop("FATAL: ", sum(is.na(causal$maf)), " causal loci missing from ", FREQ_F)
causal_wza <- if (WZA_MAF > 0) causal[maf > WZA_MAF] else causal
message("INFO: causal loci: ", nrow(causal), " total, ", nrow(causal_wza),
        if (WZA_MAF > 0) paste0(" above the extra WZA MAF cut (", WZA_MAF, ")")
        else " (no extra WZA MAF cut; Filter.maf governs)")

# Window bounds come from the SNPID the pipeline itself wrote ("{chr}:{start}-{end}"),
# not from a window size re-derived here. If compute_wza.R ever changes its binning,
# this follows it instead of silently disagreeing.
parse_windows <- function(path) {
    w <- fread(path, colClasses = c(chr = "character"))
    stopifnot(all(c("SNPID", "chr", "pos", "n_snps") %in% names(w)))
    w[, c("wstart", "wend") := {
        se <- tstrsplit(sub("^[^:]*:", "", SNPID), "-", fixed = TRUE)
        list(as.numeric(se[[1]]), as.numeric(se[[2]]))
    }]
    w
}

# ------------------------------------------------------------------- the gate
ref  <- parse_windows(tabs[[1]])

# STALE-TABLE GUARD. Every method's WZA table is built by the same rule from the same
# GEA.wza.window_size, so the window sets must be identical. They are not when a table is
# a leftover from an earlier run at a different window size -- observed exactly once, with
# 999-window tables (1 kb, this run) sitting next to 100-window tables (auto_genome_wide,
# journal 05) in the same directory. That silently produced "recall 1.000 at top-100",
# which was not a detection result at all: top-100 of a 100-window table calls the entire
# genome. Fatal, because the failure mode looks like a spectacular success.
for (m in names(tabs)[-1]) {
    o <- parse_windows(tabs[[m]])
    if (!identical(sort(o$SNPID), sort(ref$SNPID)))
        stop("FATAL: ", m, " has a different window set than ", names(tabs)[1], " (",
             nrow(o), " vs ", nrow(ref), " windows). One of these tables is stale -- ",
             "re-run mode=gea so every method is windowed at the same GEA.wza.window_size.")
}
cov  <- function(win, loci) {
    # windows x loci containment, per chromosome. Small tables -- a join is not worth it.
    vapply(seq_len(nrow(win)), function(i)
        sum(loci$chr == win$chr[i] & loci$pos >= win$wstart[i] & loci$pos <= win$wend[i]),
        integer(1))
}
ref[, n_causal      := cov(ref, causal)]
ref[, n_causal_wza  := cov(ref, causal_wza)]
ref[, is_background := as.integer(chr) >= BG_LG_MIN]

gate <- data.table(
    tag = TAG,
    n_windows              = nrow(ref),
    median_snps_per_window = median(ref$n_snps),
    min_snps_per_window    = min(ref$n_snps),
    max_snps_per_window    = max(ref$n_snps),
    n_snps_windowed        = sum(ref$n_snps),
    n_causal_total         = nrow(causal),
    n_causal_above_maf     = nrow(causal_wza),
    n_windows_with_causal  = sum(ref$n_causal_wza > 0),
    min_causal_windows     = MIN_CW)
gate[, passed := n_windows_with_causal >= MIN_CW]
fwrite(gate, file.path(OUTDIR, paste0(TAG, "_wza_gate.tsv")), sep = "\t")
message("\n=== WZA viability gate ===")
print(t(gate))

if (!gate$passed) {
    message("\nGATE FAILED: only ", gate$n_windows_with_causal,
            " windows contain a testable causal locus (need >= ", MIN_CW,
            "). WZA is not viable on this replicate at this window size; no scores written.")
    quit(save = "no", status = 0)
}

# ------------------------------------------------------------------ scoring
# Denominator is the LOCUS count WZA can reach, not the window count: a window with three
# causal loci inside is one call but three recoverable loci, and recall has to stay
# comparable to the SNP-level recall_testable it will be printed beside.
n_testable <- nrow(causal_wza)
GRID <- rbind(
    data.table(adjust = "bonf", value = c("0.01", "0.05", "0.1", "0.5", "1")),
    data.table(adjust = "qval", value = c("0.05", "0.1", "0.2")),
    data.table(adjust = "top",  value = as.character(c(10, 20, 30, 50, 100)[c(10, 20, 30, 50, 100) <= TOP_MAX])),
    data.table(adjust = "custom", value = c("1e-6", "1e-5", "1e-4", "1e-3")))

rows <- list()
for (m in names(tabs)) {
    w <- parse_windows(tabs[[m]])
    w[, n_causal_wza  := cov(w, causal_wza)]
    w[, is_background := as.integer(chr) >= BG_LG_MIN]
    tc <- setdiff(names(w), c(STRUCTURAL_COLS, "wstart", "wend", "n_causal", "n_causal_wza",
                              "is_background"))
    for (g in seq_len(nrow(GRID))) {
        adj <- GRID$adjust[g]; val <- GRID$value[g]
        sel <- rep(FALSE, nrow(w))
        for (t in tc) {
            th <- suppressMessages(compute_pval_threshold(w[[t]], adj, as.numeric(val)))
            if (identical(th$status, "ok") && !is.na(th$threshold))
                sel <- sel | (!is.na(w[[t]]) & w[[t]] < th$threshold)
        }
        cw <- w[sel]
        tp_windows <- sum(cw$n_causal_wza > 0)
        fp_windows <- sum(cw$is_background == 1L)
        loci_found <- sum(cw$n_causal_wza)
        prec <- if (tp_windows + fp_windows > 0) tp_windows / (tp_windows + fp_windows) else NA_real_
        rec  <- loci_found / n_testable
        rows[[length(rows) + 1L]] <- data.table(
            tag = TAG, method = m, adjust = adj, value = val,
            n_called_windows = nrow(cw), tp_windows = tp_windows, fp_windows = fp_windows,
            ignored_linked_windows = nrow(cw) - tp_windows - fp_windows,
            causal_loci_found = loci_found, n_testable = n_testable,
            # A rule that calls every window is not an operating point; recall is 1.0 by
            # construction. Flagged rather than dropped so the grid stays complete.
            calls_whole_genome = nrow(cw) == nrow(w),
            precision = prec, recall = rec,
            f1 = if (!is.na(prec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA_real_,
            windows_per_tp = if (tp_windows > 0) nrow(cw) / tp_windows else NA_real_)
    }
}
out <- rbindlist(rows)
fwrite(out, file.path(OUTDIR, paste0(TAG, "_wza_scores.tsv")), sep = "\t")
message("\n=== best F1 per method (window-level) ===")
print(out[!is.na(f1)][order(-f1, n_called_windows), .SD[1], by = method][
    , .(method, adjust, value, n_called_windows, tp_windows, fp_windows,
        loci = causal_loci_found, precision = round(precision, 3),
        recall = round(recall, 3), f1 = round(f1, 3))])

p <- ggplot(out[!is.na(precision) & n_called_windows > 0],
            aes(recall, precision, colour = method)) +
    geom_point(size = 2.4, alpha = 0.85) +
    scale_colour_n(uniqueN(out$method)) +
    scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
    labs(title = paste0("WZA window-level detection (", TAG, ")"),
         x = "Recall (causal loci testable by WZA)",
         y = "Precision (TP windows / (TP + background windows))", colour = NULL) +
    theme_clinego()
ggsave(file.path(OUTDIR, paste0(TAG, "_wza_pr.png")), p, width = 9, height = 6, dpi = 150)
ggsave(file.path(OUTDIR, paste0(TAG, "_wza_pr.svg")), p, width = 9, height = 6)
message("INFO: WZA outputs in ", OUTDIR)
