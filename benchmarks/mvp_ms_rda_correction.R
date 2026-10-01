#!/usr/bin/env Rscript
# =============================================================================
# mvp_ms_rda_correction.R -- supplementary figure + numbers brief for the manuscript's
# "corrected vs uncorrected RDA" comparison on the SS-Clines corpus (600 replicates), and the
# detection numbers at the operating point chosen at the Phase 3b gate (2026-10-01).
#
# WHY. Journal 18 (Steps 9-11) measured RDA two ways on the same 600 replicates: corrected
# (partial RDA, Condition(PC1..PC_k_best), the one-fit arm params_onefit/) and uncorrected
# (condition_pcs = 0, params_rdaunc/). The user decided (2026-10-01): the manuscript shows that
# comparison in one supplementary figure, and from then on uses UNCORRECTED RDA as the third
# method beside LFMM and EMMAX, with the rules 1/3, 2/3, 3/3, every method at top 0.25 % of its
# SNPs and a 5 kb agreement window. Plan: ~/.claude/plans/twinkly-purring-sparkle.md.
#
# PLOTTING + NUMBERS ONLY. Every input is a table the re-measure driver, the uncorrected-arm
# check or the harvest summaries already wrote (compute once, reuse). Nothing is re-scored.
#
# DATA SOURCES (fixed; the same rung name exists in two directories):
#   RDA comparison   remeasure600/{onefit,rdaunc}/rank_metrics.tsv, onefit/covariates.tsv,
#                    onefit_diag/rda_diag.tsv, rdaunc_check/{check,binding}.tsv
#   detection        remeasure600/rdaunc/ -- default 20-rung grid, whose top_0.0025 rung IS
#                    0.25 % and which carries the 5 kb window. NOT ledger_rdaunc/ (windows 1/2.5/4).
#                    Rules from calls_combine.tsv, single methods from calls_per_method.tsv.
#   window table     remeasure600/{onefit,rdaunc}/calls_combine.tsv, top_0.0025, windows 0-10 kb
#   LD / geometry    onefit_diag/ld_decay.tsv, params_rdaunc p-table keys, config_MVP{seed}_c1.yaml
#
# SCORING CONVENTIONS -- every number in the brief carries one tag:
#   [rank]   AUC-PR / R-precision (lib_detection.R auc_pr_from_rank): causal = TP, background-
#            neutral = the only FP, linked-neutral excluded; SNP rank = min p over traits.
#   [causal] precision = causal / (causal + background), recall = testable causal hit / testable
#            causal; F1 = 0 for a panel with no causal hit; precision medianed over panels with at
#            least one causal or background hit.
#   [A1]     the manuscript detection plane (mvp_j16_gallery.R:183): TP = causal + linked,
#            FP = background. Rebuilt here from tp_any_causal + expected_linked vs fp_background;
#            Phase 4a's rebuilt panels must reproduce these counts.
#
# OUTPUTS (FIG_OUT, default benchmarks/mvp_eval/figures_ssclines_rda/):
#   S_rda_correction.{svg,png}   AUC-PR by genetic architecture, RDA corrected / uncorrected,
#                                LFMM reference; free x-axis per architecture. (A second panel,
#                                AUC-PR gain vs adj R2 kept, was dropped 2026-10-01 by the user:
#                                the mechanism stays in the brief as numbers.)
#   numbers.tsv                  key, value, shown, unit, convention, source -- every brief number
#   window_table.tsv             agreement-window table (also read by journal 18 Step 12)
#   detection_table.tsv          detection at the operating point, rule x architecture
#   RDA_CORRECTION_BRIEF.md      the brief, generated from the tables above
#
#   FIG_OUT  BOOT_SEED (20261001)  N_BOOT (2000)
#   MVP_ARM / MVP_ADDED / MVP_N_EXPECT via benchmarks/mvp_arm.R (mandatory, 600)
# =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(ggdist)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
RM   <- file.path(EVAL, "remeasure600")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_ssclines_rda"))
BOOT_SEED   <- as.integer(Sys.getenv("BOOT_SEED", "20261001"))
N_BOOT      <- as.integer(Sys.getenv("N_BOOT", "2000"))
JITTER_SEED <- 16L                       # as mvp_ms_sim_figures.R
RUNG   <- "top_0.0025"                   # top 0.25 % of each method's SNPs
WIN_KB <- 5                              # agreement window
source(file.path(ROOT, "benchmarks/ms_style.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

man  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(man)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))
SEEDS <- sort(PRIM$seed)
N <- length(SEEDS)

rd <- function(f, ...) fread(f, colClasses = c(seed = "character"), ...)
has_all <- function(D, what) {
    miss <- setdiff(SEEDS, D$seed)
    if (length(miss)) stop(what, ": ", length(miss), " replicate(s) missing, e.g. ", miss[1])
    D[seed %in% SEEDS]
}

# ------------------------------------------------------------------ vocabulary
ARCH <- c("oliogenic" = "oligogenic", "mod-polygenic" = "moderately polygenic",
          "highly-polygenic" = "highly polygenic")
ARCH_WRAP <- c("oligogenic" = "oligogenic", "moderately polygenic" = "moderately\npolygenic",
               "highly polygenic" = "highly\npolygenic")
SERIES <- c("RDA uncorrected", "RDA corrected", "LFMM")
SER_COL <- c("RDA uncorrected" = MINOU[["amber"]], "RDA corrected" = MINOU[["navy"]],
             "LFMM" = MS_REF)
RULES <- c("1/3 methods", "2/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")

# ------------------------------------------------------------------ number store
NUM <- list()
put <- function(key, value, shown, unit = "", conv = "", src = "") {
    stopifnot(length(value) == 1, !key %in% names(NUM))
    row <- data.table(k_ = key, value = format(value, digits = 15), shown = shown,
                      unit = unit, convention = conv, source = src)
    NUM[[key]] <<- setnames(row, "k_", "key")      # data.table(key = ) would set a sort key
    invisible(shown)
}
s <- function(key) { if (!key %in% names(NUM)) stop("brief asks for unknown number: ", key); NUM[[key]]$shown }
f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
pc <- function(x) paste0(formatC(100 * x, format = "f", digits = 0), " %")
pc1 <- function(x) paste0(formatC(100 * x, format = "f", digits = 1), " %")
fp <- function(p) if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(p, format = "f", digits = 3)
md_table <- function(D) {
    D <- as.data.table(D)[, lapply(.SD, function(x) trimws(as.character(x)))]
    h <- paste0("| ", paste(names(D), collapse = " | "), " |")
    r <- paste0("|", paste(rep("---", ncol(D)), collapse = "|"), "|")
    b <- apply(D, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
    paste(c(h, r, b), collapse = "\n")
}

# =============================================================================
# 1. RDA corrected vs uncorrected
# =============================================================================
cv <- has_all(rd(file.path(RM, "onefit/covariates.tsv")), "covariates")
cv[, arch := factor(ARCH[arch_level], levels = ARCH)]
stopifnot(!anyNA(cv$arch), all(table(cv$arch) == N / 3))
rk_c <- has_all(rd(file.path(RM, "onefit/rank_metrics.tsv")), "rank_metrics onefit")
rk_u <- has_all(rd(file.path(RM, "rdaunc/rank_metrics.tsv")), "rank_metrics rdaunc")
# LFMM / EMMAX are symlinked into the uncorrected arm: their rows must be the same numbers
chk <- merge(rk_c[method != "RDA", .(seed, method, a = aucpr)], rk_u[method != "RDA", .(seed, method, b = aucpr)],
             by = c("seed", "method"))
stopifnot(nrow(chk) == 2 * N, isTRUE(all.equal(chk$a, chk$b, tolerance = 0)))
AU <- rbind(rk_c[method == "RDA", .(seed, series = "RDA corrected", aucpr, r_precision)],
            rk_u[method == "RDA", .(seed, series = "RDA uncorrected", aucpr, r_precision)],
            rk_c[method == "LFMM", .(seed, series = "LFMM", aucpr, r_precision)],
            rk_c[method == "EMMAX", .(seed, series = "EMMAX", aucpr, r_precision)])
stopifnot(all(table(AU$series) == N), nrow(AU) == 4 * N, !anyNA(AU$aucpr), all(AU$aucpr > 0))
AU <- merge(AU, cv[, .(seed, arch)], by = "seed")
SRC_RK <- "remeasure600/{onefit,rdaunc}/rank_metrics.tsv"

for (a in c("pooled", ARCH)) {
    D <- if (a == "pooled") AU else AU[arch == a]
    for (sr in c("RDA corrected", "RDA uncorrected", "LFMM", "EMMAX")) {
        x <- D[series == sr]
        k <- paste0("aucpr_", gsub(" ", "_", sr), "_", gsub(" ", "_", a))
        put(paste0(k, "_median"), median(x$aucpr), f3(median(x$aucpr)), "", "[rank]", SRC_RK)
        put(paste0(k, "_q25"), quantile(x$aucpr, 0.25), f3(quantile(x$aucpr, 0.25)), "", "[rank]", SRC_RK)
        put(paste0(k, "_q75"), quantile(x$aucpr, 0.75), f3(quantile(x$aucpr, 0.75)), "", "[rank]", SRC_RK)
        put(sub("^aucpr_", "rprec_", paste0(k, "_median")), median(x$r_precision), f3(median(x$r_precision)),
            "", "[rank]", SRC_RK)
        put(sub("^aucpr_", "rprec_", paste0(k, "_q25")), quantile(x$r_precision, 0.25),
            f3(quantile(x$r_precision, 0.25)), "", "[rank]", SRC_RK)
        put(sub("^aucpr_", "rprec_", paste0(k, "_q75")), quantile(x$r_precision, 0.75),
            f3(quantile(x$r_precision, 0.75)), "", "[rank]", SRC_RK)
    }
}

# recall denominators: testable causal loci per replicate (same set for every method)
NT_ <- merge(rk_c[method == "RDA", .(seed, n_testable)], cv[, .(seed, arch)], by = "seed")
stopifnot(rk_c[, uniqueN(n_testable), by = seed][, all(V1 == 1)])
for (a in ARCH) {
    x <- NT_[arch == a]$n_testable; k <- gsub(" ", "_", a)
    put(paste0("ntest_min_", k), min(x), f0(min(x)), "causal loci", "", SRC_RK)
    put(paste0("ntest_median_", k), median(x), f0(median(x)), "causal loci", "", SRC_RK)
    put(paste0("ntest_max_", k), max(x), f0(max(x)), "causal loci", "", SRC_RK)
}

# paired: uncorrected vs corrected, per replicate
W <- dcast(AU, seed + arch ~ series, value.var = "aucpr")
setnames(W, c("RDA corrected", "RDA uncorrected"), c("corr", "unc"))
W[, d := unc - corr][, lr := log2(unc / corr)]
stopifnot(all(is.finite(W$lr)))
set.seed(BOOT_SEED)
put("boot_seed", BOOT_SEED, as.character(BOOT_SEED), "", "", "this script")
put("n_boot", N_BOOT, f0(N_BOOT), "resamples", "", "this script")
for (a in c("pooled", ARCH)) {
    D <- if (a == "pooled") W else W[arch == a]
    k <- gsub(" ", "_", a)
    wt <- wilcox.test(D$unc, D$corr, paired = TRUE, exact = FALSE)
    bt <- replicate(N_BOOT, median(sample(D$d, replace = TRUE)))
    put(paste0("n_", k), nrow(D), f0(nrow(D)), "replicates", "", "covariates.tsv")
    put(paste0("wins_unc_", k), sum(D$unc > D$corr), f0(sum(D$unc > D$corr)), "replicates", "[rank]", SRC_RK)
    put(paste0("wins_lfmm_over_unc_", k), sum(D$LFMM > D$unc), f0(sum(D$LFMM > D$unc)), "replicates", "[rank]", SRC_RK)
    put(paste0("wilcox_V_", k), unname(wt$statistic), f0(unname(wt$statistic)), "", "[rank]", SRC_RK)
    put(paste0("wilcox_p_", k), wt$p.value, fp(wt$p.value), "", "[rank]", SRC_RK)
    put(paste0("dmed_", k), median(D$d), f3(median(D$d)), "AUC-PR", "[rank]", SRC_RK)
    put(paste0("dmed_lo_", k), quantile(bt, 0.025), f3(quantile(bt, 0.025)), "AUC-PR", "[rank]", "bootstrap")
    put(paste0("dmed_hi_", k), quantile(bt, 0.975), f3(quantile(bt, 0.975)), "AUC-PR", "[rank]", "bootstrap")
    rm_ <- median(D$unc) / median(D$corr)
    put(paste0("ratio_med_", k), rm_, formatC(rm_, format = "f", digits = 1), "x", "[rank]", SRC_RK)
}

# ---- mechanism
dg <- rd(file.path(RM, "onefit_diag/rda_diag.tsv"))
dg <- has_all(dcast(dg[key %in% c("adj_r_squared", "max_vif", "vif_flagged_predictors", "gif_lambda",
                                  "rda_axes", "condition_pcs")], seed ~ key, value.var = "value"), "rda_diag")
ck <- has_all(rd(file.path(RM, "rdaunc_check/check.tsv")), "rdaunc_check")
stopifnot(all(ck$u1_ok), all(dg$rda_axes == "2"))
MC <- Reduce(function(x, y) merge(x, y, by = "seed"),
             list(W, dg, ck[, .(seed, adj_r2_unc, gif_lambda_unc, n_snps, u2_max_rel_dev, u2_n_above_tol,
                                u2_binding_disagree)],
                  cv[, .(seed, k_best, r2_pc1_temp, K_authors)]))
stopifnot(nrow(MC) == N, all(as.integer(MC$condition_pcs) == MC$k_best))
MC[, `:=`(r2c = as.numeric(adj_r_squared), vif = as.numeric(max_vif), gifc = as.numeric(gif_lambda))]
MC[, kept := r2c / adj_r2_unc]
SRC_MC <- "onefit_diag/rda_diag.tsv + rdaunc_check/check.tsv"
put("u1_pass", sum(ck$u1_ok), f0(sum(ck$u1_ok)), "replicates", "", "rdaunc_check/check.tsv")
put("u2_max_rel_dev", max(MC$u2_max_rel_dev), formatC(max(MC$u2_max_rel_dev), format = "e", digits = 1), "", "", "rdaunc_check/check.tsv")
put("u2_above_tol", sum(MC$u2_n_above_tol), f0(sum(MC$u2_n_above_tol)), "SNPs", "", "rdaunc_check/check.tsv")
put("u2_binding_disagree", sum(MC$u2_binding_disagree), f0(sum(MC$u2_binding_disagree)), "SNPs", "", "rdaunc_check/check.tsv")
put("n_snps_total", sum(MC$n_snps), f0(sum(MC$n_snps)), "SNPs", "", "rdaunc_check/check.tsv")
put("adjr2_corr_median", median(MC$r2c), f3(median(MC$r2c)), "", "", SRC_MC)
put("adjr2_unc_median", median(MC$adj_r2_unc), f3(median(MC$adj_r2_unc)), "", "", SRC_MC)
put("kept_median", median(MC$kept), pc(median(MC$kept)), "", "", SRC_MC)
put("kept_q25", quantile(MC$kept, 0.25), pc(quantile(MC$kept, 0.25)), "", "", SRC_MC)
put("kept_q75", quantile(MC$kept, 0.75), pc(quantile(MC$kept, 0.75)), "", "", SRC_MC)
put("vif_median", median(MC$vif), formatC(median(MC$vif), format = "f", digits = 1), "", "", SRC_MC)
put("vif_ge10", sum(MC$vif >= 10), f0(sum(MC$vif >= 10)), "replicates", "", SRC_MC)
put("bio_flagged", sum(grepl("bio", MC$vif_flagged_predictors)), f0(sum(grepl("bio", MC$vif_flagged_predictors))),
    "replicates", "", SRC_MC)
put("gif_corr_median", median(MC$gifc), f3(median(MC$gifc)), "", "", SRC_MC)
put("gif_unc_median", median(MC$gif_lambda_unc), f3(median(MC$gif_lambda_unc)), "", "", SRC_MC)
put("r2pc1temp_median", median(MC$r2_pc1_temp), f2(median(MC$r2_pc1_temp)), "", "", "covariates.tsv")
put("r2pc1temp_min", min(MC$r2_pc1_temp), f2(min(MC$r2_pc1_temp)), "", "", "covariates.tsv")
put("r2pc1temp_max", max(MC$r2_pc1_temp), f2(max(MC$r2_pc1_temp)), "", "", "covariates.tsv")
for (kb in sort(unique(MC$k_best))) put(paste0("kbest_n_", kb), sum(MC$k_best == kb), f0(sum(MC$k_best == kb)),
                                        "replicates", "", "covariates.tsv")
RHO <- rbindlist(lapply(ARCH, function(a) {
    D <- MC[arch == a]
    rbindlist(lapply(c("kept", "vif", "k_best", "r2_pc1_temp"), function(v) {
        ct <- suppressWarnings(cor.test(D$lr, D[[v]], method = "spearman", exact = FALSE))
        data.table(arch = a, covariate = v, rho = unname(ct$estimate), p = ct$p.value, n = nrow(D))
    }))
}))
for (i in seq_len(nrow(RHO))) {
    k <- paste0("rho_", RHO$covariate[i], "_", gsub(" ", "_", RHO$arch[i]))
    put(k, RHO$rho[i], f2(RHO$rho[i]), "", "[rank]", SRC_MC)
    put(sub("^rho_", "rhop_", k), RHO$p[i], fp(RHO$p[i]), "", "[rank]", SRC_MC)
}
bd <- has_all(rd(file.path(RM, "rdaunc_check/binding.tsv")), "binding")
BD <- bd[, .(n = sum(n), share = sum(unc_larger) / sum(n), med_log10 = median(median_log10_unc_over_corr)), by = cls]
stopifnot(setequal(BD$cls, c("causal", "background", "linked")))
for (c0 in BD$cls) {
    put(paste0("bind_share_", c0), BD[cls == c0, share], pc1(BD[cls == c0, share]), "", "", "rdaunc_check/binding.tsv")
    put(paste0("bind_log10_", c0), BD[cls == c0, med_log10], f3(BD[cls == c0, med_log10]), "log10", "", "rdaunc_check/binding.tsv")
    put(paste0("bind_n_", c0), BD[cls == c0, n], f0(BD[cls == c0, n]), "SNPs", "", "rdaunc_check/binding.tsv")
}

# =============================================================================
# 2. detection at the operating point: top 0.25 % per method, 5 kb window
# =============================================================================
score_rows <- function(D, what) {
    D <- copy(D)
    bad <- D[n_called != tp_any_causal + expected_linked + fp_background]
    if (nrow(bad)) stop(what, ": ", nrow(bad), " row(s) with untracked SNPs (n_called != causal + linked + background)")
    D[, .(seed, n = n_called, causal = tp_any_causal, linked = expected_linked, background = fp_background,
          P = precision_strict, R = recall_testable, F1 = fifelse(is.na(f1), 0, f1))]
}
detection_for <- function(arm) {
    cc <- rd(file.path(RM, arm, "calls_combine.tsv"))[rung == RUNG]
    pm <- rd(file.path(RM, arm, "calls_per_method.tsv"))[rung == RUNG]
    un <- cc[combine == "union"]
    stopifnot(un[, uniqueN(n_called), by = seed][, all(V1 == 1)])     # union is window-invariant
    rbind(score_rows(un[window_kb == 0], "union")[, rule := "1/3 methods"],
          score_rows(cc[combine == "ge2" & window_kb == WIN_KB], "ge2")[, rule := "2/3 methods"],
          score_rows(cc[combine == "all3" & window_kb == WIN_KB], "all3")[, rule := "3/3 methods"],
          score_rows(pm[method == "LFMM"], "LFMM")[, rule := "LFMM"],
          score_rows(pm[method == "RDA"], "RDA")[, rule := "RDA"],
          score_rows(pm[method == "EMMAX"], "EMMAX")[, rule := "EMMAX"])
}
DET <- detection_for("rdaunc")
stopifnot(all(table(DET$rule) == N), nrow(DET) == 6 * N, setequal(unique(DET$seed), SEEDS))
DET <- merge(DET, cv[, .(seed, arch)], by = "seed")
det_sum <- function(D) D[, .(n_median = median(n), causal_median = median(causal), linked_median = median(linked),
                             background_median = median(background), TP_A1_median = median(causal + linked),
                             P_median = median(P, na.rm = TRUE), R_median = median(R), F1_median = median(F1),
                             empty = sum(n == 0), lt3 = sum(n < 3), n_rep = .N), by = rule]
DT <- rbind(det_sum(DET)[, arch := "pooled"], DET[, det_sum(.SD), by = arch][, arch := as.character(arch)], fill = TRUE)
DT[, rule := factor(rule, levels = RULES)][, arch := factor(arch, levels = c("pooled", ARCH))]
setorder(DT, arch, rule)
fwrite(DT, file.path(OUT, "detection_table.tsv"), sep = "\t")
SRC_DET <- "remeasure600/rdaunc/{calls_combine,calls_per_method}.tsv, top_0.0025, window 5 kb"
for (i in seq_len(nrow(DT))) {
    k <- paste0("det_", gsub("[ /]", "_", DT$rule[i]), "_", gsub(" ", "_", DT$arch[i]))
    put(paste0(k, "_n"), DT$n_median[i], f0(DT$n_median[i]), "SNPs", "[A1]", SRC_DET)
    put(paste0(k, "_P"), DT$P_median[i], f2(DT$P_median[i]), "", "[causal]", SRC_DET)
    put(paste0(k, "_R"), DT$R_median[i], f3(DT$R_median[i]), "", "[causal]", SRC_DET)
    put(paste0(k, "_F1"), DT$F1_median[i], f3(DT$F1_median[i]), "", "[causal]", SRC_DET)
    put(paste0(k, "_empty"), DT$empty[i], f0(DT$empty[i]), "replicates", "", SRC_DET)
}
# the same rules with CORRECTED RDA as the third method -- what the switch changed (pooled)
DC <- det_sum(detection_for("onefit"))[, rule := factor(rule, levels = RULES)][order(rule)]
for (i in seq_len(nrow(DC))) {
    k <- paste0("detcorr_", gsub("[ /]", "_", DC$rule[i]))
    put(paste0(k, "_n"), DC$n_median[i], f0(DC$n_median[i]), "SNPs", "[A1]", "remeasure600/onefit")
    put(paste0(k, "_P"), DC$P_median[i], f2(DC$P_median[i]), "", "[causal]", "remeasure600/onefit")
    put(paste0(k, "_R"), DC$R_median[i], f3(DC$R_median[i]), "", "[causal]", "remeasure600/onefit")
    put(paste0(k, "_F1"), DC$F1_median[i], f3(DC$F1_median[i]), "", "[causal]", "remeasure600/onefit")
    put(paste0(k, "_empty"), DC$empty[i], f0(DC$empty[i]), "replicates", "", "remeasure600/onefit")
}

# =============================================================================
# 3. agreement-window table + LD / genome geometry + config facts
# =============================================================================
WT <- rbindlist(lapply(c("rdaunc", "onefit"), function(arm) {
    cc <- rd(file.path(RM, arm, "calls_combine.tsv"))[rung == RUNG & combine %in% c("ge2", "all3")]
    cc[, .(empty = sum(n_called == 0), lt3 = sum(n_called < 3), size_median = median(n_called),
           P_median = median(precision_strict, na.rm = TRUE), R_median = median(recall_testable),
           F1_median = median(fifelse(is.na(f1), 0, f1)),
           linked_share_median = median(expected_linked[n_called > 0] / n_called[n_called > 0]),
           n_rep = .N), by = .(combine, window_kb)][, third_method := c(rdaunc = "RDA uncorrected",
                                                                         onefit = "RDA corrected")[[arm]]]
}))
stopifnot(all(WT$n_rep == N), setequal(WT$window_kb, c(0, 1, 2.5, 5, 10)))
WT[, rule := c(ge2 = "2/3 methods", all3 = "3/3 methods")[combine]]
setcolorder(WT, c("third_method", "rule", "window_kb"))[, combine := NULL]
setorder(WT, third_method, rule, window_kb)
fwrite(WT, file.path(OUT, "window_table.tsv"), sep = "\t")
for (i in seq_len(nrow(WT))) {
    k <- paste0("win_", c("RDA uncorrected" = "unc", "RDA corrected" = "corr")[[WT$third_method[i]]], "_",
                c("2/3 methods" = "ge2", "3/3 methods" = "all3")[[WT$rule[i]]], "_", WT$window_kb[i])
    put(paste0(k, "_empty"), WT$empty[i], f0(WT$empty[i]), "replicates", "", "window_table.tsv")
}

ld <- has_all(rd(file.path(RM, "onefit_diag/ld_decay.tsv"))[group == "All" & scope == "genome_wide"], "ld_decay")
stopifnot(nrow(ld) == N, all(ld$method == "hill_weir"))
put("ld_half_min", min(ld$half_decay_bp), f2(min(ld$half_decay_bp) / 1000), "kb", "", "onefit_diag/ld_decay.tsv")
put("ld_half_median", median(ld$half_decay_bp), f2(median(ld$half_decay_bp) / 1000), "kb", "", "onefit_diag/ld_decay.tsv")
put("ld_half_max", max(ld$half_decay_bp), f2(max(ld$half_decay_bp) / 1000), "kb", "", "onefit_diag/ld_decay.tsv")
put("ld_r2_02_median", median(ld$r2_02_bp), f0(median(ld$r2_02_bp)), "bp", "", "onefit_diag/ld_decay.tsv")
put("ld_r2_02_max", max(ld$r2_02_bp), f0(max(ld$r2_02_bp)), "bp", "", "onefit_diag/ld_decay.tsv")
put("ld_r2_intercept_median", median(ld$r2_intercept), f3(median(ld$r2_intercept)), "", "", "onefit_diag/ld_decay.tsv")

GEO <- rbindlist(lapply(SEEDS, function(sd) {
    k <- fread(file.path(EVAL, "params_rdaunc", paste0("MVP", sd), "c1", "RDA_pvalues.tsv"), select = 1L,
               colClasses = "character")[[1]]
    ch <- sub(":.*", "", k); ps <- as.numeric(sub(".*:", "", k))
    data.table(seed = sd, n_lg = uniqueN(ch), max_pos = max(ps), n_snps = length(k))
}))
put("lg_n_min", min(GEO$n_lg), f0(min(GEO$n_lg)), "linkage groups", "", "params_rdaunc p-table keys")
put("lg_n_max", max(GEO$n_lg), f0(max(GEO$n_lg)), "linkage groups", "", "params_rdaunc p-table keys")
put("lg_maxpos", max(GEO$max_pos), f0(max(GEO$max_pos)), "bp", "", "params_rdaunc p-table keys")
put("snps_min", min(GEO$n_snps), f0(min(GEO$n_snps)), "SNPs", "", "params_rdaunc p-table keys")
put("snps_median", median(GEO$n_snps), f0(median(GEO$n_snps)), "SNPs", "", "params_rdaunc p-table keys")
put("snps_max", max(GEO$n_snps), f0(max(GEO$n_snps)), "SNPs", "", "params_rdaunc p-table keys")
put("top025_n_min", ceiling(0.0025 * min(GEO$n_snps)), f0(ceiling(0.0025 * min(GEO$n_snps))), "SNPs per trait", "", "derived")
put("top025_n_max", ceiling(0.0025 * max(GEO$n_snps)), f0(ceiling(0.0025 * max(GEO$n_snps))), "SNPs per trait", "", "derived")

clump <- vapply(SEEDS, function(sd) {
    l <- readLines(file.path(ROOT, paste0("config_MVP", sd, "_c1.yaml")))
    v <- as.numeric(sub(".*snp_clumping_distance:\\s*([0-9]+).*", "\\1",
                        grep("^\\s*snp_clumping_distance:", l, value = TRUE)))
    if (!length(v) || anyNA(v) || any(v != v[1])) stop("config_MVP", sd, "_c1.yaml: snp_clumping_distance unreadable or mixed")
    v[1]
}, numeric(1))
stopifnot(all(clump == 5000))
put("clump_5000_configs", sum(clump == 5000), f0(sum(clump == 5000)), "configs", "", "config_MVP{seed}_c1.yaml")
dflt <- readLines(file.path(ROOT, "scripts/clinego.app/inst/config_default.yaml"))
dv <- unique(as.numeric(sub(".*snp_clumping_distance:\\s*([0-9]+).*", "\\1",
                            grep("^\\s*snp_clumping_distance:", dflt, value = TRUE))))
stopifnot(length(dv) == 1)
put("clump_pipeline_default", dv, f0(dv / 1000), "kb", "", "scripts/clinego.app/inst/config_default.yaml")

# =============================================================================
# 4. figure
# =============================================================================
AF <- AU[series %in% SERIES]
AF[, series := factor(series, levels = SERIES)][, ypos := as.numeric(factor(series, levels = rev(SERIES)))]
pa <- ggplot(AF, aes(aucpr, ypos)) +
    stat_halfeye(aes(fill = series, group = series), orientation = "horizontal", adjust = 0.8,
                 height = 0.55, justification = -0.2, .width = 0, point_colour = NA, slab_alpha = 0.75) +
    geom_boxplot(aes(group = series), orientation = "y", width = 0.14, outlier.shape = NA,
                 fill = "white", colour = MS_INK, linewidth = 0.25) +
    geom_point(aes(y = ypos - 0.25, colour = series),
               position = position_jitter(height = 0.08, width = 0, seed = JITTER_SEED),
               size = 0.6, alpha = 0.6, stroke = 0) +
    # free x per architecture: AUC-PR spans ~0-1 (oligogenic) but ~0-0.15 (highly polygenic),
    # so a shared axis would flatten the highly polygenic panel
    facet_grid(cols = vars(arch), scales = "free_x", labeller = labeller(arch = ARCH_WRAP)) +
    scale_y_continuous(breaks = seq_along(SERIES), labels = rev(SERIES)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
    scale_fill_manual(values = SER_COL, guide = "none") + scale_colour_manual(values = SER_COL, guide = "none") +
    labs(x = "AUC-PR (causal loci vs background)", y = NULL) + theme_ms()
ms_save(file.path(OUT, "S_rda_correction"), pa, 180, 60)
message("  S_rda_correction            180 x  60 mm")

# =============================================================================
# 5. numbers.tsv + brief
# =============================================================================
NT <- rbindlist(NUM)
fwrite(NT, file.path(OUT, "numbers.tsv"), sep = "\t")

a3 <- function(stem) paste(vapply(c("oligogenic", "moderately_polygenic", "highly_polygenic"),
                                  function(a) s(paste0(stem, a)), ""), collapse = " / ")
tab_auc <- rbindlist(lapply(c("RDA corrected", "RDA uncorrected", "LFMM", "EMMAX"), function(sr) {
    k <- gsub(" ", "_", sr)
    data.table(series = sr,
               oligogenic = sprintf("%s (%s-%s)", s(paste0("aucpr_", k, "_oligogenic_median")),
                                    s(paste0("aucpr_", k, "_oligogenic_q25")), s(paste0("aucpr_", k, "_oligogenic_q75"))),
               `moderately polygenic` = sprintf("%s (%s-%s)", s(paste0("aucpr_", k, "_moderately_polygenic_median")),
                                    s(paste0("aucpr_", k, "_moderately_polygenic_q25")), s(paste0("aucpr_", k, "_moderately_polygenic_q75"))),
               `highly polygenic` = sprintf("%s (%s-%s)", s(paste0("aucpr_", k, "_highly_polygenic_median")),
                                    s(paste0("aucpr_", k, "_highly_polygenic_q25")), s(paste0("aucpr_", k, "_highly_polygenic_q75"))),
               pooled = sprintf("%s (%s-%s)", s(paste0("aucpr_", k, "_pooled_median")),
                                s(paste0("aucpr_", k, "_pooled_q25")), s(paste0("aucpr_", k, "_pooled_q75"))),
               `R-precision, pooled` = sprintf("%s (%s-%s)", s(paste0("rprec_", k, "_pooled_median")),
                                               s(paste0("rprec_", k, "_pooled_q25")), s(paste0("rprec_", k, "_pooled_q75"))))
}))
tab_pair <- rbindlist(lapply(c("oligogenic", "moderately_polygenic", "highly_polygenic", "pooled"), function(k)
    data.table(architecture = gsub("_", " ", k), n = s(paste0("n_", k)),
               `uncorrected > corrected` = s(paste0("wins_unc_", k)),
               `median paired diff (95 % CI)` = sprintf("%s (%s to %s)", s(paste0("dmed_", k)), s(paste0("dmed_lo_", k)),
                                                       s(paste0("dmed_hi_", k))),
               `ratio of medians` = paste0(s(paste0("ratio_med_", k)), "x"),
               `Wilcoxon V` = s(paste0("wilcox_V_", k)), `Wilcoxon p` = s(paste0("wilcox_p_", k)),
               `LFMM > uncorrected` = s(paste0("wins_lfmm_over_unc_", k)))))
tab_rho <- RHO[, .(architecture = arch,
                   covariate = c(kept = "adj R2 kept after correction", vif = "max VIF, corrected fit",
                                 k_best = "k_best (PCs conditioned on)", r2_pc1_temp = "R2(PC1 ~ temperature)")[covariate],
                   `Spearman rho` = f2(rho), p = vapply(p, fp, ""), n = n)]
tab_bind <- BD[, .(class = c(causal = "causal", background = "background neutral", linked = "linked neutral")[cls],
                   SNPs = f0(n), `uncorrected p larger` = pc1(share), `median log10(p_unc / p_corr)` = f3(med_log10))]
tab_det <- DT[, .(architecture = as.character(arch), rule = as.character(rule), `median n` = f0(n_median),
                  `causal` = f0(causal_median), `linked` = f0(linked_median), `background` = f0(background_median),
                  `TP [A1]` = f0(TP_A1_median), `P [causal]` = f2(P_median), `R [causal]` = f3(R_median),
                  `F1 [causal]` = f3(F1_median), empty = empty, `n < 3` = lt3)]
tab_detc <- DC[, .(rule = as.character(rule), `median n` = f0(n_median), `P [causal]` = f2(P_median),
                   `R [causal]` = f3(R_median), `F1 [causal]` = f3(F1_median), empty = empty, `n < 3` = lt3)]
tab_win <- WT[, .(`third method` = third_method, rule, `window, kb` = window_kb, empty, `n < 3` = lt3,
                  `median n` = f0(size_median), `P [causal]` = f2(P_median), `R [causal]` = f3(R_median),
                  `F1 [causal]` = f3(F1_median), `linked share` = pc(linked_share_median))]

brief <- c(
"# Corrected vs uncorrected RDA, and detection at the operating point — numbers brief",
"",
"Base for writing the simulation section (SS-Clines, 600 replicates). Generated by",
"`benchmarks/mvp_ms_rda_correction.R`; every number below is in `numbers.tsv` with its source table",
"and scoring convention. No number in this file is typed by hand. Lab record: journal 18",
"(`work/journal/18_ssclines_onefit_rule.Rmd`), Steps 9-12.",
"",
"## 1. Decisions (user, Phase 3b gate)",
"",
"- 2026-09-30: RDA uncorrected (`condition_pcs = 0`) is run as a GEA method beside corrected RDA, on all",
"  600 replicates (plan `~/.claude/plans/glittery-fluttering-seahorse.md`).",
"- 2026-10-01: the manuscript first compares corrected and uncorrected RDA (one supplementary figure,",
"  `S_rda_correction`); from then on **uncorrected RDA is the third method** beside LFMM and EMMAX, under",
"  the rules 1/3, 2/3, 3/3 and the three single methods. Supersedes the 2026-09-30 decision to keep",
"  corrected (one-fit) RDA in the panels.",
"- 2026-10-01: operating point **top 0.25 %** of each method's SNPs; agreement window **5 kb**",
sprintf("  (= `snp_clumping_distance` in all %s replicate configs; section 6).", s("clump_5000_configs")),
"- pmax / intersection-union combination of the two RDA fits is **not reported** (not a GEA-literature",
"  method). It exists only in the lab record.",
"",
"## 2. Methods facts",
"",
"- RDA: `scripts/rda.R`, one fit per replicate, `vegan::rda`, predictors `bio_1`, `bio_2` (temperature,",
"  salinity optima, `benchmarks/convert_mvp.R:20`). **Corrected** = partial RDA, `Condition(PC1..PC_k)`, k = `k_best`",
sprintf("  (k = 2: %s, 3: %s, 4: %s, 5: %s, 7: %s, 8: %s replicates); **uncorrected** = the same call with",
        s("kbest_n_2"), s("kbest_n_3"), s("kbest_n_4"), s("kbest_n_5"), s("kbest_n_7"), s("kbest_n_8")),
"  `condition_pcs = 0`. SNP p-values from `rdadapt` on the retained axes (2 in every replicate, both fits),",
sprintf("  divided by the genomic inflation factor (median GIF corrected %s, uncorrected %s).",
        s("gif_corr_median"), s("gif_unc_median")),
"- The uncorrected arm re-issued each replicate's exact Snakemake `rda.R` command with `condition_pcs` set",
"  to 0 and outputs redirected (`benchmarks/mvp_rda_uncorrected_run.sh`; same image `cline-go:latest`",
sprintf("  93bd6025, same inputs, seed 42, 99 permutations). Checks: U1 passed on %s / 600; U2 —", s("u1_pass")),
sprintf("  max(p_corrected, p_uncorrected) reproduces the earlier two-fit arm on %s SNPs, max relative",
        s("n_snps_total")),
sprintf("  deviation %s, %s SNPs above 1e-9, %s disagreements about which fit gave the larger p.",
        s("u2_max_rel_dev"), s("u2_above_tol"), s("u2_binding_disagree")),
"- LFMM (K = `k_best`) and EMMAX (`k_best` PCs) p-values are unchanged from the one-fit arm (byte-identical).",
sprintf("- Calling rule: top 0.25 %% of SNPs per method **per predictor** (LFMM, EMMAX test `bio_1` and `bio_2`"),
sprintf("  separately), unioned over the two; RDA has one multivariate p. N = ceiling(0.0025 x tested SNPs) ="),
sprintf("  %s-%s SNPs per predictor (%s-%s tested SNPs,",
        s("top025_n_min"), s("top025_n_max"), s("snps_min"), s("snps_max")),
sprintf("  median %s).", s("snps_median")),
"- Rules: **1/3** = called by at least one method (union; window-independent by construction); **2/3** =",
"  a called SNP with a call from at least one other method within 5 kb (the pipeline's `Overlap` /",
"  `Cross-method` strategy, `workflow/rules/common.smk:299`, verified equal on 12 replicates x 100 checks);",
"  **3/3** = calls from all three methods within 5 kb (not a pipeline strategy; strictest reference).",
"- Scoring conventions (tags in the tables):",
"  - `[rank]` AUC-PR = mean precision at each causal locus in the ranked SNP list, causal loci not",
"    retrieved count 0; background-neutral SNPs are the only false positives, linked-neutral SNPs are",
"    excluded; SNP rank = min p over the two predictors. R-precision = share of causal loci among the top R SNPs,",
"    R = number of testable causal loci.",
"  - `[causal]` precision = causal / (causal + background); recall = testable causal loci hit / testable",
"    causal loci; F1 = 0 when a panel has no causal hit; precision medianed over panels with >= 1 causal",
"    or background hit.",
"  - `[A1]` the manuscript's detection-plane convention: TP = causal + linked-neutral hits, FP =",
"    background-neutral hits (linkage groups 11-20).",
"",
"## 3. Result 1 — corrected vs uncorrected RDA: ranking of causal loci",
"",
"AUC-PR `[rank]`, median (IQR), n = 200 replicates per architecture:",
"",
md_table(tab_auc),
"",
"Paired, uncorrected minus corrected RDA, per replicate (bootstrap: seed `BOOT_SEED` =",
sprintf("%s, %s resamples):", s("boot_seed"), s("n_boot")),
"",
md_table(tab_pair),
"",
"## 4. Result 2 — why: structure correction removes the climate signal on this landscape",
"",
sprintf("- adj R2 of the RDA model: corrected median %s, uncorrected %s; the corrected fit keeps a median %s",
        s("adjr2_corr_median"), s("adjr2_unc_median"), s("kept_median")),
sprintf("  (IQR %s-%s) of the uncorrected model's adj R2.", s("kept_q25"), s("kept_q75")),
sprintf("- Corrected fit: max VIF median %s; max VIF >= 10 in %s / 600 replicates; a climate predictor (`bio_*`)",
        s("vif_median"), s("vif_ge10")),
sprintf("  among the VIF-flagged terms in %s / 600 — the conditioning PCs are collinear with the predictors.",
        s("bio_flagged")),
"- Spearman correlation of log2(AUC-PR uncorrected / corrected) with each covariate, within architecture.",
"  Note: `adj R2 kept` and the AUC-PR ratio both measure the loss caused by the correction, so their",
"  correlation describes rather than proves; `k_best` is the independent covariate (more PCs conditioned",
"  on, larger gain without correction). PC1 alone does not explain it.",
"",
md_table(tab_rho),
"",
sprintf("- R2(PC1 ~ temperature): median %s, range %s-%s.", s("r2pc1temp_median"), s("r2pc1temp_min"), s("r2pc1temp_max")),
"- Per SNP (all 600 replicates): where the uncorrected fit gives the larger p (relative tolerance 1e-9):",
"",
md_table(tab_bind),
"",
"  On this corpus the correction did not preferentially suppress background SNPs: dropping it lowers",
"  causal p-values more than background ones. This is a per-SNP p-value comparison, not a false-positive",
"  rate at a threshold, and it does not carry over to landscapes where structure does not follow climate",
"  (section 8).",
"",
"## 5. Result 3 — detection at the operating point (uncorrected RDA as third method)",
"",
"top 0.25 % per method, 5 kb window; medians over replicates (n = 600 pooled, 200 per architecture);",
"`empty` / `n < 3` = number of replicates. Source: `detection_table.tsv`.",
"",
"**Recall is relative to each replicate's testable causal loci, which differ by orders of magnitude",
sprintf("between architectures** (median, range): oligogenic %s (%s-%s), moderately polygenic %s (%s-%s),",
        s("ntest_median_oligogenic"), s("ntest_min_oligogenic"), s("ntest_max_oligogenic"),
        s("ntest_median_moderately_polygenic"), s("ntest_min_moderately_polygenic"), s("ntest_max_moderately_polygenic")),
sprintf("highly polygenic %s (%s-%s). Compare recall (and F1) between rules **within** an architecture; the",
        s("ntest_median_highly_polygenic"), s("ntest_min_highly_polygenic"), s("ntest_max_highly_polygenic")),
"pooled rows mix the three scales and describe no single architecture.",
"",
md_table(tab_det),
"",
"Same rules with **corrected** RDA as the third method (pooled, for the record of what the switch changed):",
"",
md_table(tab_detc),
"",
"## 6. Agreement window — why 5 kb",
"",
sprintf("- Pipeline default `snp_clumping_distance` = %s kb (`config_default.yaml`) — larger than a whole",
        s("clump_pipeline_default")),
sprintf("  linkage group here (%s linkage groups per replicate, positions up to %s bp), so every hit on a",
        if (s("lg_n_min") == s("lg_n_max")) s("lg_n_min") else paste0(s("lg_n_min"), "-", s("lg_n_max")), s("lg_maxpos")),
sprintf("  linkage group would \"agree\". Every replicate config sets it to 5000 bp (%s / 600 checked); the",
        s("clump_5000_configs")),
"  pipeline uses this same parameter as the window of its Overlap (>= 2 methods) strategy.",
sprintf("- LD (Hill-Weir fit, all 600): half-decay %s-%s kb (median %s kb); the pipeline's `auto` mode (distance",
        s("ld_half_min"), s("ld_half_max"), s("ld_half_median")),
sprintf("  at r2 = 0.2) gives a median %s bp (max %s bp), because r2 at distance 0 is only ~%s — not usable.",
        s("ld_r2_02_median"), s("ld_r2_02_max"), s("ld_r2_intercept_median")),
"- 5 kb is ~3x the LD half-decay. Effect of the window at top 0.25 % (`window_table.tsv`):",
"",
md_table(tab_win),
"",
sprintf("- With uncorrected RDA: 2/3 is never empty at any window; 3/3 is empty in %s / %s / %s / %s / %s replicates",
        s("win_unc_all3_0_empty"), s("win_unc_all3_1_empty"), s("win_unc_all3_2.5_empty"),
        s("win_unc_all3_5_empty"), s("win_unc_all3_10_empty")),
"  at 0 / 1 / 2.5 / 5 / 10 kb — 3/3 is not \"never empty\"; at the chosen 5 kb it is empty in 2. A wider",
"  window raises 3/3's recall and lowers 2/3's precision.",
"",
"## 7. Figure — `S_rda_correction.svg` / `.png` (180 x 60 mm, one panel) — ready for the supplement",
"",
"- AUC-PR `[rank]` per replicate by genetic architecture (n = 200 each): RDA uncorrected (amber),",
"  RDA corrected (navy), LFMM (grey, reference). Half-violin = density; box = median and IQR; dots =",
"  replicates. **Each architecture has its own x-axis range** (AUC-PR spans ~0-1 in oligogenic but",
"  ~0-0.15 in highly polygenic replicates) — say so in the legend; compare within a panel, not across.",
"- The mechanism (section 4: adj R2 kept, VIF, Spearman rho) is text-only; no figure panel.",
"",
"## 8. Do not claim / open",
"",
"- Do not generalise the uncorrected advantage: it follows from structure being collinear with the",
"  climate gradients in SS-Clines (section 4). Where structure does not follow climate, uncorrected RDA is",
"  expected to inflate false positives; the pipeline keeps both options.",
"- \"Lotterhos 2023 ran RDA uncorrected\" appears only as a code comment",
"  (`benchmarks/mvp_write_sweep_configs.R:40-42`). **Unverified** — check the paper (`shelf grep`) before it",
"  enters the text.",
"- pmax (max of the two fits' p) is not a reported method (section 1).",
"- Offset accuracy with the new panels is **pending** (Phase 4: panels rebuilt with uncorrected RDA,",
"  top 0.25 %, 5 kb; offsets re-run). Detection figures (A1/D5/D6) are redrawn in Phase 5 from those",
"  panels, which must reproduce the `[A1]` counts in section 5."
)
stopifnot(!any(grepl("NA", brief[grepl("^\\|", brief)], fixed = TRUE)))
writeLines(brief, file.path(OUT, "RDA_CORRECTION_BRIEF.md"))
message(sprintf("wrote numbers.tsv (%d numbers), window_table.tsv (%d rows), detection_table.tsv (%d rows), brief (%d lines)",
                nrow(NT), nrow(WT), nrow(DT), length(brief)))
message("OK: ", N, " replicates -> ", OUT)
