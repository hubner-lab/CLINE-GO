# detection_table_land.R -- detection_table.tsv for one landscape (plan nifty-prancing-grove Phase 8).
# The block "2. detection at the operating point" of work/ld_window_lab/scratch_figs/rda_correction_ld.R,
# copied verbatim (score_rows, detection_for, the provenance + top-N budget checks, det_sum): the
# remeasure harness's own medians, which fig2_s2_land.R ties the gallery's detection counts to (two code
# paths, one operating point). The rest of that script (Figure S1 and its SS-Clines-only diagnostics)
# is not run: the user chose Figure 2 + S2 only (2026-11-06).
# Covariates (architecture) from the uncorrected arm's own covariates.tsv instead of remeasure600/onefit/.
# Env: REM_ROOT (holds rdaunc_rda2x_ldb/, onefit_rda2x_ldb/), FIG_OUT, MVP_SEEDS_FILE + MVP_N_EXPECT.
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
RM   <- Sys.getenv("REM_ROOT"); OUT <- Sys.getenv("FIG_OUT")
stopifnot(nzchar(RM), nzchar(OUT), dir.exists(RM))
RUNG   <- "top_0.0025"                   # top 0.25 % per predictor (RDA: x RDA_MULT)
RDA_MULT <- 2                            # RDA share = 0.25 % x number of predictors (2)
ARM_UNC  <- "rdaunc_rda2x_ldb"
ARM_CORR <- "onefit_rda2x_ldb"
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
ARCH <- c("oliogenic" = "oligogenic", "mod-polygenic" = "moderately polygenic",
          "highly-polygenic" = "highly polygenic")
RULES <- c("1/3 methods", "2/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")
cv <- has_all(rd(file.path(RM, ARM_UNC, "covariates.tsv")), "covariates")
cv[, arch := factor(ARCH[arch_level], levels = ARCH)]
stopifnot(!anyNA(cv$arch), all(table(cv$arch) == N / 3))

# ---- verbatim from rda_correction_ld.R, "2. detection at the operating point"
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
          score_rows(cc[combine == "ge2" & window_kb > 0], "ge2")[, rule := "2/3 methods"],
          score_rows(cc[combine == "all3" & window_kb > 0], "all3")[, rule := "3/3 methods"],
          score_rows(pm[method == "LFMM"], "LFMM")[, rule := "LFMM"],
          score_rows(pm[method == "RDA"], "RDA")[, rule := "RDA"],
          score_rows(pm[method == "EMMAX"], "EMMAX")[, rule := "EMMAX"])
}
for (arm in c(ARM_UNC, ARM_CORR)) {
    pv <- readLines(file.path(RM, arm, "provenance.tsv"))
    if (!any(pv == sprintf("share_mult\tLFMM=1,EMMAX=1,RDA=%g", RDA_MULT)))
        stop(arm, ": provenance does not record share_mult RDA=", RDA_MULT)
    pm <- rd(file.path(RM, arm, "calls_per_method.tsv"))[rung == RUNG]
    nt <- rd(file.path(RM, arm, "rank_metrics.tsv"))[, .(seed, method, n_snps)]
    chk <- merge(pm, nt, by = c("seed", "method"))
    chk[, want := ceiling(0.0025 * fifelse(method == "RDA", RDA_MULT, 1) * n_snps)]
    if (nrow(chk) != 3L * N || any(chk$value_applied != chk$want)) stop(arm, ": top N is not the RDA x ", RDA_MULT, " budget")
}
DET <- detection_for(ARM_UNC)
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
message(sprintf("wrote detection_table.tsv (%d rows, %d replicates) -> %s", nrow(DT), N, OUT))
