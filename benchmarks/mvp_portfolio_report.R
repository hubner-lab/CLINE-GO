#!/usr/bin/env Rscript
# =============================================================================
# mvp_portfolio_report.R -- aggregate journal 06's per-cell outputs into the tables
# and figures the lab journal actually prints.
#
# Reads (all optional -- whatever is present is aggregated, missing pieces are reported
# as absent rather than silently skipped):
#   sweep/{PROJ}/{TAG}_threshold_sweep.tsv   homogeneous grid, per cell x axis
#   sweep/{PROJ}/rank.tsv                    AUC-PR per method x cell x axis
#   portfolio/{TAG}_portfolio.tsv            heterogeneous subsets x combine windows
#   redundancy/{TAG}_method_clusters.tsv     method families
#   profile/{TAG}_maf_cut.tsv                recall below vs above the MAF cut
#   wza/{TAG}_wza_gate.tsv                   window-calling viability
#
# Writes benchmarks/mvp_eval/report06/.
#
# Every table here answers one of the questions journal 06 was opened to settle:
#   ladder_aucpr        which structure-correction rung each method wants
#   solo_vs_combine     does combining beat the best single method, and by how much
#   combine_window      does agreement-within-a-window beat exact-key agreement
#   gapit_marginal      does adding GAPIT models past the family representative help
#   headline            the single best operating point per configuration group
#
# Usage:
#   Rscript mvp_portfolio_report.R --proj=MVP1232548 [--arm-tag=_m05] [--call-cap=100]
#       [--eval-dir=DIR] [--outdir=DIR]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt  <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
PROJ     <- opt("proj", "MVP1232548")
EVAL_DIR <- opt("eval_dir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval"))
# --arm-tag suffixes every input subdirectory AND the default output directory together, so a
# parallel arm cannot half-read the other arm's tables: "" is the base arm, "_m05" the
# Filter.maf 0.05 arm. Splitting these into separate flags would make a mismatch expressible.
ARM_TAG  <- opt("arm_tag", "")
OUTDIR   <- opt("outdir", file.path(EVAL_DIR, paste0("report06", ARM_TAG)))
subdir   <- function(n) file.path(EVAL_DIR, paste0(n, ARM_TAG))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# GAPIT models, so "did the extra models help" can be asked without hardcoding the answer.
GAPIT_MODELS <- c("GLM", "MLM", "CMLM", "ECMLM", "SUPER", "MLMM", "FarmCPU", "BLINK")

say <- function(...) message("INFO: ", ...)
absent <- function(what, path) message("ABSENT: ", what, " -- ", path,
                                       " not found; that section is omitted, not empty")

# ------------------------------------------------------------------ 1. ladders
rank_f <- file.path(subdir("sweep06"), PROJ, "rank.tsv")
if (file.exists(rank_f)) {
    R <- fread(rank_f)
    # eval_detection.R writes long rows keyed by a "PROJ|METHOD|CELL|AXIS" config label.
    R <- R[metric %in% c("auc_pr", "best_f1", "n_causal_testable")]
    parts <- tstrsplit(R$config, "|", fixed = TRUE)
    R[, `:=`(proj = parts[[1]], method = parts[[2]], cell = parts[[3]], axis = parts[[4]])]
    L <- dcast(R, proj + method + cell + axis ~ metric, value.var = "value")
    setorder(L, axis, method, cell)
    fwrite(L, file.path(OUTDIR, "ladder_aucpr.tsv"), sep = "\t")

    best <- L[!is.na(auc_pr), .(best_cell = cell[which.max(auc_pr)],
                                best_auc  = max(auc_pr),
                                default_auc = auc_pr[cell == "c3"][1]), by = .(method, axis)]
    best[, gain_over_default := best_auc / default_auc]
    fwrite(best, file.path(OUTDIR, "ladder_best_rung.tsv"), sep = "\t")
    say("ladder: ", nrow(L), " method x cell x axis rows")

    p <- ggplot(L[!is.na(auc_pr)], aes(cell, auc_pr, group = method, colour = method)) +
        geom_line(linewidth = 0.7) + geom_point(size = 2) +
        facet_wrap(~axis) + scale_colour_n(uniqueN(L$method)) +
        labs(title = paste0("Structure-correction ladder, AUC-PR (", PROJ, ")"),
             x = "cell (n_pcs / LFMM K rung)", y = "AUC-PR", colour = NULL) +
        theme_clinego()
    ggsave(file.path(OUTDIR, "ladder_aucpr.png"), p, width = 11, height = 5.5, dpi = 150)
    ggsave(file.path(OUTDIR, "ladder_aucpr.svg"), p, width = 11, height = 5.5)
} else absent("ladder", rank_f)

# ------------------------------------------------- 2-4. portfolio-derived tables
port_fs <- list.files(subdir("portfolio"), pattern = "_portfolio\\.tsv$",
                      full.names = TRUE)
if (length(port_fs)) {
    P <- rbindlist(lapply(port_fs, fread), fill = TRUE)
    say("portfolio: ", nrow(P), " operating points from ", length(port_fs), " file(s)")

    # ---- solo vs combine. Compared INSIDE one table so both sides use the same truth
    # denominator; comparing across two files is how denominators drift apart unnoticed.
    solo_best <- P[rule == "solo" & !is.na(f1)][order(-f1, n_called), .SD[1], by = tag][
        , .(tag, solo_subset = subset, solo_f1 = f1, solo_prec = precision,
            solo_recall = recall_testable, solo_called = n_called)]
    comb_best <- P[rule != "solo" & !is.na(f1)][order(-f1, n_called), .SD[1], by = tag][
        , .(tag, combine_rule = rule, combine_subset = subset, combine_cw = combine_window_kb,
            combine_assignment = assignment, combine_f1 = f1, combine_prec = precision,
            combine_recall = recall_testable, combine_called = n_called)]
    SC <- merge(solo_best, comb_best, by = "tag", all = TRUE)
    SC[, f1_gain := combine_f1 - solo_f1]
    fwrite(SC, file.path(OUTDIR, "solo_vs_combine.tsv"), sep = "\t")

    # ---- combine window. The paired comparison: hold (tag, subset, rule, assignment)
    # fixed and vary ONLY the window, so the delta cannot be a subset effect in disguise.
    CW <- P[!is.na(combine_window_kb) & rule %in% c("union", "intersection") |
            grepl("^at_least_", rule)]
    CW <- CW[!is.na(f1) & !is.na(combine_window_kb)]
    if (nrow(CW)) {
        base <- CW[combine_window_kb == 0,
                   .(tag, subset, rule, assignment, f1_w0 = f1, prec_w0 = precision,
                     rec_w0 = recall_testable, called_w0 = n_called)]
        CWD <- merge(CW[combine_window_kb > 0], base,
                     by = c("tag", "subset", "rule", "assignment"))
        CWD[, `:=`(f1_delta = f1 - f1_w0, prec_delta = precision - prec_w0,
                   recall_delta = recall_testable - rec_w0)]
        fwrite(CWD, file.path(OUTDIR, "combine_window_effect.tsv"), sep = "\t")
        summ <- CWD[, .(comparisons = .N, median_f1_delta = median(f1_delta),
                        wins = sum(f1_delta > 0), ties = sum(f1_delta == 0),
                        losses = sum(f1_delta < 0)), by = .(rule, combine_window_kb)]
        setorder(summ, rule, combine_window_kb)
        fwrite(summ, file.path(OUTDIR, "combine_window_summary.tsv"), sep = "\t")
        message("\n=== does a combine window beat exact-key agreement? ===")
        print(summ)

        p <- ggplot(CWD, aes(factor(combine_window_kb), f1_delta, fill = rule)) +
            geom_hline(yintercept = 0, colour = CLINEGO_THRESHOLD, linetype = "dashed") +
            geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
            scale_fill_n(uniqueN(CWD$rule)) +
            labs(title = paste0("F1 change from windowed cross-method agreement (", PROJ, ")"),
                 x = "combine window (kb)", y = "F1 - F1 at exact-key agreement", fill = NULL) +
            theme_clinego()
        ggsave(file.path(OUTDIR, "combine_window_effect.png"), p, width = 9, height = 5.5, dpi = 150)
        ggsave(file.path(OUTDIR, "combine_window_effect.svg"), p, width = 9, height = 5.5)
    } else absent("combine-window comparison", "portfolio rows with a window")

    # ---- marginal value of extra GAPIT models. A subset is "GAPIT-heavy" when it carries
    # more than one GAPIT model: if the family really is redundant, those should not beat
    # subsets carrying exactly one.
    P[, n_gapit := vapply(strsplit(subset, "+", fixed = TRUE),
                          function(v) sum(v %in% GAPIT_MODELS), integer(1))]
    GM <- P[rule != "solo" & !is.na(f1), .(best_f1 = max(f1), n = .N),
            by = .(tag, n_methods, n_gapit)][order(tag, n_methods, n_gapit)]
    fwrite(GM, file.path(OUTDIR, "gapit_marginal.tsv"), sep = "\t")

    # ---- headline. TWO rankings, because on this benchmark they disagree and the
    # disagreement is the finding, not a nuisance.
    #
    # `f1`         max F1 on (precision_strict, recall_testable) -- the metric of record.
    # `capped`     max F1 subject to n_called <= CALL_CAP.
    #
    # Why both: precision_strict charges ONLY background_neutral calls as false positives
    # (Lotterhos's convention), and roughly half this genome is linked_neutral. So a rule can
    # call 1 458 SNPs, hit 31 causal loci, put every other call on a linked_neutral site, and
    # score precision_strict = 1.000 -- while the user actually reads 47 SNPs per real hit.
    # That is a correct number under the published convention and a useless operating point.
    #
    # The cap is also where the `top <= 100` instruction actually binds. It constrains the
    # `top` grid only; a badly calibrated method at bonf_0.05 (GLM calls 1 184 SNPs here,
    # 11% of the genome) blows past it through a rule the cap never touched, and union /
    # at_least_k across methods can exceed any single method's cap by construction.
    CALL_CAP <- as.numeric(opt("call_cap", "100"))
    head_cols <- function(d) d[, .(tag, subset, rule, combine_window_kb, assignment, n_called,
                                   tp, tp_any_causal, precision, precision_all,
                                   recall_testable, f1, calls_per_tp)]
    HEAD <- head_cols(P[!is.na(f1)][order(-f1, n_called), .SD[1], by = tag])[, ranking := "max_f1"]
    CAP  <- head_cols(P[!is.na(f1) & n_called <= CALL_CAP][order(-f1, n_called), .SD[1], by = tag])
    if (nrow(CAP)) CAP[, ranking := paste0("max_f1_at_n_called_le_", CALL_CAP)]
    HEAD <- rbind(HEAD, CAP, fill = TRUE)
    fwrite(HEAD, file.path(OUTDIR, "headline.tsv"), sep = "\t")
    message("\n=== headline operating point per tag, both rankings ===")
    print(HEAD[order(tag, ranking),
               .(tag, ranking, subset, rule, cw = combine_window_kb, n_called, tp,
                 prec_strict = round(precision, 3), prec_all = round(precision_all, 3),
                 rec = round(recall_testable, 3), f1 = round(f1, 3),
                 per_tp = round(calls_per_tp, 1))])
} else absent("portfolio", subdir("portfolio"))

# ------------------------------------------------------- 5. side tables, verbatim
for (nm in list(c("redundancy", "_method_clusters.tsv", "method_families.tsv"),
                c("profile",    "_maf_cut.tsv",         "maf_cut.tsv"),
                c("profile",    "_maf_recall.tsv",      "maf_recall.tsv"),
                c("wza",        "_wza_gate.tsv",        "wza_gate.tsv"),
                c("wza",        "_wza_scores.tsv",      "wza_scores.tsv"))) {
    fs <- list.files(subdir(nm[1]), pattern = paste0(nm[2], "$"), full.names = TRUE)
    if (!length(fs)) { absent(nm[3], file.path(subdir(nm[1]), paste0("*", nm[2]))); next }
    fwrite(rbindlist(lapply(fs, fread), fill = TRUE), file.path(OUTDIR, nm[3]), sep = "\t")
    say("copied ", length(fs), " file(s) -> ", nm[3])
}

message("\nINFO: journal-06 report tables in ", OUTDIR)
