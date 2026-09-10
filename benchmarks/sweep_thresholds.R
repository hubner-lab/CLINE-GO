#!/usr/bin/env Rscript
# =============================================================================
# sweep_thresholds.R — the FREE half of the GEA hyperparameter search.
#
# Model fitting (LFMM K, EMMAX n_pcs, RDA condition_pcs) is expensive and each
# knob gets its own ladder. Significance CALLING is not: adjust x threshold, and
# every cross-method combine rule, are pure post-processing on p-value tables
# that already exist. This script sweeps all of that in one pass, which is what
# turns a combinatorial grid into a sum of independent ladders plus free
# post-hoc scoring.
#
# Sweeps, per method (journal-06 grid; see TOP_MAX below for why `top` is capped):
#   bonf   0.01 0.05 0.1 0.5 1
#   qval   0.05 0.1 0.2
#   top    10 20 30 50 100
#   custom 1e-6 1e-5 1e-4 1e-3
#
# The grid still runs past conventional significance. The premise is that a user reads
# the operating point off the Manhattan plot: climate-association signal is noisy, so the
# useful question is the shape of the whole precision/recall tradeoff, not the single point
# FDR 0.05 lands on. bonf 1 is the "one expected false positive genome-wide" line
# (cutoff = 1/n_tested).
#
# Combine rules across methods, at each shared (adjust, value) and each combine window:
#   union, intersection, at_least_2, and a threshold-free rank_sum top-N.
# Combining is now WINDOWED (--combine-window-kb): two methods calling SNPs within W bp of
# each other count as agreeing, which is what the pipeline's own combine_sigsnps.R does over
# GEA.snp_clumping_distance. W=0 is the old exact-key behaviour and is asserted identical.
#
# HETEROGENEOUS per-method operating points are NOT done here -- see
# benchmarks/sweep_portfolio.R, which reads this script's output to choose each method's
# candidate points and then combines across method subsets.
#
# Outputs: long-format rows appended to --out, a wide summary TSV, and a
# precision-recall scatter with the Pareto frontier drawn through it.
#
# Usage:
#   Rscript sweep_thresholds.R --methods="EMMAX=/p/a.tsv,LFMM=/p/b.tsv,RDA=/p/c.tsv" \
#       --truth=FILE --tag=STR [--tp-window-kb=0] [--combine-window-kb=0,1,2.5,5] \
#       [--top-max=100] [--out=FILE] [--outdir=DIR]
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

TAG       <- req("tag")
TRUTH_F   <- req("truth")
OUTDIR    <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/laruson_eval"))
OUT_TSV   <- opt("out", file.path(OUTDIR, "benchmark_eval.tsv"))
WINDOW    <- as.numeric(opt("tp_window_kb", "0"))

spec <- strsplit(req("methods"), ",", fixed = TRUE)[[1]]
methods <- list()
for (s in spec) {
    kv <- strsplit(s, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Bad --methods entry (expected NAME=PATH): ", s)
    if (!file.exists(kv[2])) stop("Method p-value table not found: ", kv[2])
    methods[[kv[1]]] <- kv[2]
}
if (!length(methods)) stop("--methods produced no entries")

# TOP_MAX caps every rank-based rule (`top` and the rank_sum top-N vector).
#
# Journal 05 ran `top` to 5000 because the highly-polygenic replicates (1231798, 395 causal;
# 1232923, 363) were still rising monotonically at 1000, and a censored frontier reported as
# a frontier understates every rule. Journal 06 restricts the benchmark to ONE replicate
# (1232548, 70 causal) and caps at 100 by instruction: calling 5000 of ~10 300 SNPs is half
# the genome and is not an operating point any user would accept, whatever its F1.
#
# On this replicate the cap costs nothing measurable -- 70 causal loci, so top-100 can still
# reach recall 1.0. It WOULD censor the polygenic seeds, which is one of the reasons the
# single-replicate choice and the cap go together. Raise it only alongside a seed that needs it.
TOP_MAX  <- as.numeric(opt("top_max", "100"))
TOP_GRID <- c(10, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000)
TOP_GRID <- TOP_GRID[TOP_GRID <= TOP_MAX]

GRID <- rbind(
    data.table(adjust = "bonf",   value = c("0.01", "0.05", "0.1", "0.5", "1")),
    # qval above 0.2 is not a threshold anyone reports; it survived in journal 05 only to
    # extend the frontier into the high-recall tail that TOP_MAX now bounds anyway.
    data.table(adjust = "qval",   value = c("0.05", "0.1", "0.2")),
    data.table(adjust = "top",    value = as.character(TOP_GRID)),
    data.table(adjust = "custom", value = c("1e-6", "1e-5", "1e-4", "1e-3"))
)

# Cross-method agreement windows, in kb. 0 reproduces journal 05's exact-key set operations
# (asserted below). The non-zero rungs are calibrated on this replicate's own truth table --
# linked_neutral -> nearest causal is median 2 351 bp, p75 4 496, p90 6 819 -- and 5 kb is
# GEA.snp_clumping_distance from the config, i.e. what the pipeline itself would use.
COMBINE_WINDOWS <- as.numeric(strsplit(opt("combine_window_kb", "0"), ",", fixed = TRUE)[[1]])
if (any(is.na(COMBINE_WINDOWS))) stop("--combine-window-kb must be a comma-separated numeric list")

truth <- load_truth(TRUTH_F)
message("INFO: truth = ", truth$n_simulated, " simulated causal loci")

# ------------------------------------------------------------------ load once
loaded <- list()
for (m in names(methods)) {
    lp <- load_pvalues(methods[[m]], "all")
    tc <- testable_causal(lp$pv, lp$trait_cols, truth)
    loaded[[m]] <- list(pv = lp$pv, trait_cols = lp$trait_cols,
                        testable = tc, n_testable = length(tc),
                        rank = rank_keys(lp$pv, lp$trait_cols))
    message(sprintf("INFO: %-8s %d SNPs, traits=[%s], testable causal=%d",
                    m, nrow(lp$pv), paste(lp$trait_cols, collapse = ","), length(tc)))
}

# The combine denominator is the union of what any method could have found —
# using a single method's testable set would misreport combined recall.
testable_union <- unique(unlist(lapply(loaded, `[[`, "testable")))
n_testable_union <- length(testable_union)
message("INFO: union testable causal across methods = ", n_testable_union)

E   <- new_emitter()
sum_rows <- list()
add_summary <- function(config, kind, adjust, value, s, testable_keys, combine_window_kb = 0) {
    sum_rows[[length(sum_rows) + 1L]] <<- data.table(
        tag = TAG, config = config, kind = kind, adjust = adjust, value = value,
        combine_window_kb = combine_window_kb,
        n_called = s$n_called, tp = s$tp, tp_any_causal = s$tp_any_causal,
        expected_linked = s$expected_linked,
        fp_background = s$fp_background, fn_testable = s$fn_testable,
        precision = s$precision_strict, precision_all = s$precision_all,
        recall_testable = s$recall_testable, recall_simulated = s$recall_simulated,
        f1 = s$f1, calls_per_tp = s$calls_per_tp, n_testable = length(testable_keys))
}

# ------------------------------------------------------- per-method threshold
called_cache <- list()   # [[adjust_value]][[method]] = called key vector
for (g in seq_len(nrow(GRID))) {
    adj <- GRID$adjust[g]; val <- GRID$value[g]
    key <- paste0(adj, "_", val)
    called_cache[[key]] <- list()
    for (m in names(loaded)) {
        L   <- loaded[[m]]
        res <- call_by_threshold(L$pv, L$trait_cols, adj, val, quiet = TRUE)
        called_cache[[key]][[m]] <- res$called
        s   <- score_calls(res$called, WINDOW, truth, L$testable)
        cfg <- paste0(TAG, "|", m, "|", key)
        E$emit_scored(cfg, s, WINDOW)
        add_summary(cfg, m, adj, val, s, L$testable)
    }
}

# ------------------------------------------------------------- combine rules
# Homogeneous rules: every method at the SAME (adjust, value). Heterogeneous per-method
# operating points are the job of benchmarks/sweep_portfolio.R, which reads this script's
# output to pick each method's candidate points first.
for (key in names(called_cache)) {
    parts <- called_cache[[key]]
    if (length(parts) < 2) next
    adj <- sub("_[^_]*$", "", key); val <- sub("^.*_", "", key)

    for (cw in COMBINE_WINDOWS) {
        support <- combine_support(parts, cw * 1000)
        combos  <- list(
            union        = support$snp,
            intersection = support[n_methods >= length(parts), snp],
            at_least_2   = support[n_methods >= 2, snp]
        )
        for (cname in names(combos)) {
            s   <- score_calls(combos[[cname]], WINDOW, truth, testable_union)
            cfg <- paste0(TAG, "|COMBINE_", cname, "|", key,
                          if (cw > 0) paste0("|cw", cw, "kb") else "")
            E$emit_scored(cfg, s, WINDOW)
            add_summary(cfg, paste0("combine_", cname), adj, val, s, testable_union, cw)
        }
    }
}

# ---- unit test: at window 0, combine_support() must reproduce the exact-key set logic
# journal 05 used. Asserted rather than assumed -- this is the one place where a silent
# change would rewrite every combine number in the series without any error surfacing.
local({
    key <- names(called_cache)[1]
    parts <- called_cache[[key]]
    if (length(parts) >= 2) {
        tab <- table(unlist(parts, use.names = FALSE))
        sup <- combine_support(parts, 0)
        chk <- function(a, b, what) if (!setequal(a, b))
            stop("combine_support(window=0) disagrees with exact-key ", what,
                 ": ", length(a), " vs ", length(b))
        chk(sup$snp,                              names(tab),                       "union")
        chk(sup[n_methods >= 2, snp],             names(tab)[tab >= 2],             "at_least_2")
        chk(sup[n_methods >= length(parts), snp], names(tab)[tab == length(parts)], "intersection")
        message("INFO: window-0 combine identity check passed on ", key)
    }
})

# ------------------------------------------------- threshold-free rank_sum
# Each method ranks every SNP it tested; ranks are summed across methods and the
# top N taken. Methods that did not test a SNP contribute a worst-case rank, so
# a SNP missing from one method is penalised rather than silently favoured.
rank_tables <- lapply(names(loaded), function(m) {
    r <- loaded[[m]]$rank$keys
    # Column is `snp`, not `key`: data.table(key = ...) would be taken as the
    # constructor's key= argument and create no column.
    data.table(snp = r, rank = seq_along(r), method = m)
})
rt    <- rbindlist(rank_tables)
worst <- rt[, .(worst = max(rank)), by = method]
wide  <- dcast(rt, snp ~ method, value.var = "rank")
for (m in names(loaded)) {
    w <- worst[method == m, worst]
    set(wide, which(is.na(wide[[m]])), m, w + 1L)
}
wide[, rank_sum := rowSums(.SD), .SDcols = names(loaded)]
setorder(wide, rank_sum)
for (n in TOP_GRID) {
    if (n > nrow(wide)) next
    s   <- score_calls(wide$snp[seq_len(n)], WINDOW, truth, testable_union)
    cfg <- paste0(TAG, "|COMBINE_rank_sum|top_", n)
    E$emit_scored(cfg, s, WINDOW)
    add_summary(cfg, "combine_rank_sum", "top", as.character(n), s, testable_union)
}

# ------------------------------------------------------------------- outputs
summary_dt <- rbindlist(sum_rows)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
sum_path <- file.path(OUTDIR, paste0(TAG, "_threshold_sweep.tsv"))
fwrite(summary_dt, sum_path, sep = "\t")
message("INFO: wrote ", sum_path)
append_eval(E$collect(), OUT_TSV)

# Pareto frontier over (recall_testable, precision): a config is on the frontier
# when no other config beats it on both axes at once.
pf <- summary_dt[!is.na(precision) & !is.na(recall_testable) & n_called > 0]
if (nrow(pf)) {
    setorder(pf, -recall_testable, -precision)
    front <- pf[0]; best_p <- -Inf
    for (i in seq_len(nrow(pf))) {
        if (pf$precision[i] > best_p) { front <- rbind(front, pf[i]); best_p <- pf$precision[i] }
    }
    fwrite(front, file.path(OUTDIR, paste0(TAG, "_pareto.tsv")), sep = "\t")

    p <- ggplot(pf, aes(x = recall_testable, y = precision, colour = kind)) +
        geom_point(size = 2.4, alpha = 0.85) +
        geom_step(data = front[order(recall_testable)],
                  aes(x = recall_testable, y = precision),
                  inherit.aes = FALSE, direction = "vh",
                  colour = CLINEGO_THRESHOLD, linetype = "dashed", linewidth = 0.5) +
        scale_colour_n(uniqueN(pf$kind)) +
        scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
        labs(title = paste0("Causal-locus detection: precision vs recall (", TAG, ")"),
             x = "Recall (causal loci in the tested SNP set)",
             y = "Precision (TP / (TP + background-neutral FP))",
             colour = "Method / rule") +
        theme_clinego()
    ggsave(file.path(OUTDIR, paste0(TAG, "_precision_recall.png")), p,
           width = 9, height = 6, dpi = 150)
    ggsave(file.path(OUTDIR, paste0(TAG, "_precision_recall.svg")), p, width = 9, height = 6)
    message("INFO: wrote precision-recall plot + Pareto frontier (", nrow(front), " configs)")
}

# Console report: best by F1 and best recall at >=0.9 precision, per kind.
message("\n=== best F1 per method/rule ===")
best_f1 <- summary_dt[!is.na(f1)][order(-f1, n_called), .SD[1], by = kind]
print(best_f1[, .(kind, adjust, value, cw = combine_window_kb, n_called, tp, fp_background,
                  expected_linked, precision = round(precision, 3),
                  precision_all = round(precision_all, 3),
                  recall = round(recall_testable, 3), f1 = round(f1, 3),
                  per_tp = round(calls_per_tp, 1))])
message("\n=== best recall at precision >= 0.9 ===")
hp <- summary_dt[!is.na(precision) & precision >= 0.9][order(-recall_testable), .SD[1], by = kind]
print(hp[, .(kind, adjust, value, n_called, tp, fp_background,
             precision = round(precision, 3), recall = round(recall_testable, 3))])
