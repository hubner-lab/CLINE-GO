#!/usr/bin/env Rscript
# =============================================================================
# mvp_detection_600.R -- threshold-free detection (AUC-PR) for the SS-Clines corpus.
#
# WHY. The detection arm (AUC-PR per GEA method) was only ever scored on the 14
# legacy SS-Mtn seeds (journal 07, report07/aucpr_by_arm.tsv). The 600 SS-Clines
# replicates carry rule-level panel composition only (panel_pr_recomputed.tsv, fixed
# operating points). The raw p-value tables survived the 2026-09-23 cleanup, so the
# single-method ranking quality can be scored on the same corpus the offset arm
# reports, without re-running any GEA.
#
# CONVENTION -- verbatim the one journal 07 used (mvp_method_redundancy.R:75-86):
#   * ranking score = min p over the method's trait columns (rank_keys); LFMM/EMMAX
#     have bio_1 + bio_2, RDA has ONE multivariate column -- a different ranking basis,
#     stated wherever the numbers are shown;
#   * truth = truth_any.tsv (causal on either axis); positives = causal loci,
#     false positives = background_neutral only; linked hits are neither;
#   * denominator = causal loci present in the tested SNP set (testable_causal);
#   * no genomic control (the harvested tables are what production mode=gea wrote).
# SS-Clines ran one structure-correction rung (--cells=1 -> c1 = k_best everywhere),
# which is the legacy default rung c3.
#
# GATE. Before scoring the corpus, one legacy seed is re-scored from params07/c3 and
# must reproduce redundancy07/<seed>_c3_method_clusters.tsv to 1e-12 for all three
# methods. If it does not, the script stops: an AUC-PR on an unconfirmed convention is
# not shipped.
#
#   OUT_DIR   default benchmarks/mvp_eval/detection600
#   CELL      default c1
#   NCORES    default 16
#   MVP_ARM / MVP_ADDED / MVP_N_EXPECT   via benchmarks/mvp_arm.R (mandatory for SS-Clines)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OUT  <- Sys.getenv("OUT_DIR", file.path(EVAL, "detection600"))
CELL <- Sys.getenv("CELL", "c1")
NCOR <- as.integer(Sys.getenv("NCORES", "16"))
METHODS <- c("LFMM", "RDA", "EMMAX")

source(file.path(ROOT, "scripts/R/utils/pval_threshold.R"))   # lib_detection.R expects it
source(file.path(ROOT, "benchmarks/lib_detection.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

score_seed <- function(seed, pdir, truth_f) {
    truth <- load_truth(truth_f)
    rbindlist(lapply(METHODS, function(m) {
        f <- file.path(pdir, paste0(m, "_pvalues.tsv"))
        if (!file.exists(f)) stop("missing p-value table: ", f)
        lp <- load_pvalues(f, "all")
        rk <- rank_keys(lp$pv, lp$trait_cols)
        tc <- testable_causal(lp$pv, lp$trait_cols, truth)
        data.table(seed = seed, method = m,
                   aucpr = auc_pr_from_rank(rk$keys, truth, length(tc)),
                   n_testable = length(tc), n_snps = nrow(lp$pv),
                   n_traits = length(lp$trait_cols))
    }))
}

# ---------------------------------------------------------------- 1. the gate
GATE_SEED <- "1231288"
ref <- fread(file.path(EVAL, "redundancy07", sprintf("MVP%s_c3_method_clusters.tsv", GATE_SEED)))
got <- score_seed(GATE_SEED, file.path(EVAL, "params07", paste0("MVP", GATE_SEED), "c3"),
                  file.path(ROOT, "data/mvp", paste0("MVP", GATE_SEED), "truth_any.tsv"))
cmp <- merge(got[, .(method, aucpr, n_testable)],
             ref[method %in% METHODS, .(method, ref_aucpr = auc_pr, ref_testable = n_testable)],
             by = "method")
cmp[, abs_diff := abs(aucpr - ref_aucpr)]
fwrite(cmp, file.path(OUT, "gate_legacy_reproduction.tsv"), sep = "\t")
print(cmp)
if (nrow(cmp) != length(METHODS) || any(cmp$abs_diff > 1e-12) ||
    any(cmp$n_testable != cmp$ref_testable))
    stop("GATE FAILED: legacy seed ", GATE_SEED, " does not reproduce journal-07 AUC-PR")
message("GATE PASSED: journal-07 AUC-PR reproduced on seed ", GATE_SEED, " (max |diff| = ",
        format(max(cmp$abs_diff)), ")")

# ------------------------------------------------------------ 2. the corpus
man  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(man)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))

res <- mclapply(PRIM$seed, function(s)
    tryCatch(score_seed(s, file.path(EVAL, "params", paste0("MVP", s), CELL),
                        file.path(ROOT, "data/mvp", paste0("MVP", s), "truth_any.tsv")),
             error = function(e) data.table(seed = s, method = NA_character_,
                                            aucpr = NA_real_, error = conditionMessage(e))),
    mc.cores = NCOR, mc.preschedule = FALSE)
A <- rbindlist(res, fill = TRUE)
if ("error" %in% names(A) && any(!is.na(A$error))) {
    print(A[!is.na(error)])
    stop("scoring failed on ", A[!is.na(error), uniqueN(seed)], " seed(s)")
}
if ("error" %in% names(A)) A[, error := NULL]

# Covariates, parsed exactly as mvp_pleiotropy_report.R:92-97 does.
cov <- PRIM[, .(seed, arch_level, architecture, ispleiotropy, group, added, final_LA,
                n_causal_maf01, k_best)]
cov[, regime := fifelse(grepl("unequal-S", architecture), "unequal-S", "equal-S")]
cov[, pleiotropy := fifelse(ispleiotropy == 1L, "pleiotropy", "no pleiotropy")]
cov[, block := sub("^ssclines_", "", group)]
A <- merge(A, cov[, .(seed, arch_level, regime, pleiotropy, block, final_LA,
                      n_causal_maf01, k_best)], by = "seed")

# ------------------------------------------------------------- 3. the gates
stopifnot(nrow(A) == length(METHODS) * nrow(PRIM),
          uniqueN(A$seed) == nrow(PRIM),
          !anyNA(A$aucpr), all(A$aucpr >= 0 & A$aucpr <= 1),
          all(A[method == "RDA", n_traits] == 1L),
          all(A[method != "RDA", n_traits] == 2L))
setorder(A, seed, method)
fwrite(A, file.path(OUT, "aucpr_per_seed.tsv"), sep = "\t")

S <- A[, .(n = .N, median = median(aucpr), q25 = quantile(aucpr, .25), q75 = quantile(aucpr, .75)),
       by = .(method, arch_level)]
setorder(S, arch_level, -median)
fwrite(S, file.path(OUT, "aucpr_summary.tsv"), sep = "\t")
print(S)
message(sprintf("OK: %d rows, %d seeds, arm [%s] -> %s", nrow(A), uniqueN(A$seed),
                mvp_arm_label(), OUT))
