#!/usr/bin/env Rscript
# =============================================================================
# mvp_guard_cost_600.R -- what did RDA's intersection-union guard cost, SNP by SNP?
#
# WHY. Until Phase 2 (commit 3d4ccab) scripts/rda.R reported pmax(p_partial, p_unconstrained):
# the intersection-union test (Berger 1982), valid and conservative. The one-fit change reports
# p_partial alone -- a deliberate tradeoff (docs/gea-simulation-reanalysis.md), whose price this
# script measures on all 600 SS-Clines replicates. Both arms are on disk: the frozen pmax arm
# (benchmarks/mvp_eval/params/, chmod a-w) and the one-fit re-run (params_onefit/). The partial
# fit is unchanged between them (Phase 2a anchor), so per SNP p_new = p_partial and
# pmax_old = max(p_partial, p_unconstrained).
#
# WHAT IT MEASURES, per replicate:
#   guard_binding.tsv  which fit was binding, per SNP, split by truth class:
#                        lower  pmax_old > p_new beyond TOL -> the unconstrained fit was binding
#                        equal  |pmax_old - p_new| <= TOL   -> the partial fit was binding
#                      This is the share of the GENOME the guard touched. It is NOT the share of
#                      calls it removed (dossier #correction 2026-09-29) -- that is the next table.
#                      Also AUC-PR / R-precision of both definitions.
#   guard_rungs.tsv    at every rung of the re-measure grid: the RDA call set of each definition,
#                        removed = called by one-fit, not by pmax   (what the guard blocked)
#                        added   = called by pmax, not by one-fit
#                      each split into causal (TP) / background_neutral (FP) / linked_neutral /
#                      untracked. `added` is empty BY CONSTRUCTION only where the cutoff is fixed
#                      (bonf, custom: p_new <= pmax_old moves no SNP out). For top N the cutoff is
#                      the N-th smallest p of each definition and for qval pi0 moves with the
#                      distribution, so both differences are real there and both are reported.
#   assert_guard.tsv   what the assertions below checked
#
# SCORING CONVENTION as mvp_remeasure_600.R / lib_detection.R: exact position; causal = TP,
# background_neutral = the only FP, linked_neutral neither; calls via call_by_threshold()
# (the pipeline's per-trait `p <= cutoff`). RDA has one multivariate column.
#
# ASSERTIONS -- each stops the run, none warns:
#   G1 same SNP keys in the same order in both arms, same NA pattern.
#   G2 p_new <= pmax_old * (1 + TOL) on every SNP. TOL (default 1e-9) is relative: block
#      ssclines_ncline_ctredge was frozen with snakemake -c4 while the re-run used -c2, and its
#      RDA p differs at ~1e-13 relative (Phase 2b, dossier 2026-09-30). A bitwise test would
#      count those ~315k SNPs as "unconstrained fit binding"; n_violation_exact / n_lower_exact
#      keep the bitwise counts visible.
#   G3 bonf and custom: added is empty (fixed cutoff, nested by construction).
#   G4 both arms' n_called, tp and fp_background at every rung equal the RDA rows of
#      remeasure600/{pmax,onefit}/calls_per_method.tsv, and both AUC-PRs equal their
#      rank_metrics.tsv to 1e-12 -- this script and the re-measure driver score the same thing.
#   G5 600 replicates, no replicate missing, no scoring error.
#
#   OLD_DIR    default benchmarks/mvp_eval/params          (frozen pmax arm, read-only)
#   NEW_DIR    default benchmarks/mvp_eval/params_onefit
#   OLD_REMEASURE / NEW_REMEASURE  default benchmarks/mvp_eval/remeasure600/{pmax,onefit}
#                                  (grid.tsv is read from NEW_REMEASURE, so the rungs line up)
#   OUT_DIR    default benchmarks/mvp_eval/remeasure600/guard_cost
#   CELL       default c1;  NCORES default 16;  TOL default 1e-9
#   MVP_ARM / MVP_ADDED / MVP_N_EXPECT   via benchmarks/mvp_arm.R (mandatory for SS-Clines)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OLD  <- Sys.getenv("OLD_DIR", file.path(EVAL, "params"))
NEW  <- Sys.getenv("NEW_DIR", file.path(EVAL, "params_onefit"))
RM_OLD <- Sys.getenv("OLD_REMEASURE", file.path(EVAL, "remeasure600", "pmax"))
RM_NEW <- Sys.getenv("NEW_REMEASURE", file.path(EVAL, "remeasure600", "onefit"))
OUT  <- Sys.getenv("OUT_DIR", file.path(EVAL, "remeasure600", "guard_cost"))
CELL <- Sys.getenv("CELL", "c1")
NCOR <- as.integer(Sys.getenv("NCORES", "16"))
TOL  <- as.numeric(Sys.getenv("TOL", "1e-9"))

source(file.path(ROOT, "scripts/R/utils/pval_threshold.R"))
source(file.path(ROOT, "benchmarks/lib_detection.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GRID <- fread(file.path(RM_NEW, "grid.tsv"))
GRID_OLD <- fread(file.path(RM_OLD, "grid.tsv"))
if (!isTRUE(all.equal(GRID, GRID_OLD)))
    stop("grid.tsv differs between ", RM_OLD, " and ", RM_NEW, " -- the arms were not re-measured on one grid")
NESTED <- c("bonf", "custom")

classify <- function(keys, truth) {
    c(causal    = sum(keys %in% truth$causal$key),
      background = sum(keys %in% truth$bg_keys),
      linked    = sum(keys %in% truth$linked_keys),
      untracked = sum(!keys %in% c(truth$causal$key, truth$bg_keys, truth$linked_keys)))
}

# ------------------------------------------------------------------ one replicate
guard_seed <- function(seed) {
    truth <- load_truth(file.path(ROOT, "data/mvp", paste0("MVP", seed), "truth_any.tsv"))
    rd <- function(dir) {
        f <- file.path(dir, paste0("MVP", seed), CELL, "RDA_pvalues.tsv")
        if (!file.exists(f)) stop("missing RDA p-value table: ", f)
        lp <- load_pvalues(f, "all")
        if (length(lp$trait_cols) != 1L) stop("RDA table with ", length(lp$trait_cols), " p columns: ", f)
        lp
    }
    o <- rd(OLD); n <- rd(NEW)
    if (!identical(o$pv$key, n$pv$key) || !identical(o$trait_cols, n$trait_cols))       # G1
        stop(sprintf("G1 FAILED seed %s: SNP keys or column differ between arms", seed))
    po <- o$pv[[o$trait_cols]]; pn <- n$pv[[n$trait_cols]]
    if (!identical(is.na(po), is.na(pn)))
        stop(sprintf("G1 FAILED seed %s: NA pattern differs between arms", seed))
    ok <- !is.na(po)
    viol  <- ok & pn > po * (1 + TOL)
    if (any(viol))                                                                      # G2
        stop(sprintf("G2 FAILED seed %s: %d SNP(s) with p_new > pmax_old beyond TOL (max rel %.3g)",
                     seed, sum(viol), max((pn[viol] - po[viol]) / po[viol])))
    lower <- ok & po > pn * (1 + TOL)
    equal <- ok & !lower
    keys  <- o$pv$key
    cls   <- fifelse(keys %in% truth$causal$key, "causal",
             fifelse(keys %in% truth$bg_keys, "background",
             fifelse(keys %in% truth$linked_keys, "linked", "untracked")))
    testable <- testable_causal(n$pv, n$trait_cols, truth)
    ro <- rank_keys(o$pv, o$trait_cols); rn <- rank_keys(n$pv, n$trait_cols)
    rpo <- r_precision_from_rank(ro$keys, truth, length(testable))
    rpn <- r_precision_from_rank(rn$keys, truth, length(testable))

    bind <- data.table(seed = seed, n_snps = sum(ok),
        n_violation_exact = sum(ok & pn > po), n_lower = sum(lower), n_equal = sum(equal),
        n_lower_exact = sum(ok & pn < po),
        causal_n = sum(ok & cls == "causal"),         causal_lower = sum(lower & cls == "causal"),
        background_n = sum(ok & cls == "background"), background_lower = sum(lower & cls == "background"),
        linked_n = sum(ok & cls == "linked"),         linked_lower = sum(lower & cls == "linked"),
        median_p_old = median(po[ok]), median_p_new = median(pn[ok]),
        min_p_old = min(po[ok]), min_p_new = min(pn[ok]),
        # causal loci: median -log10 p gain where the unconstrained fit was binding
        causal_median_log10_gain = if (any(lower & cls == "causal"))
            median(log10(po[lower & cls == "causal"]) - log10(pn[lower & cls == "causal"])) else NA_real_,
        n_testable = length(testable),
        aucpr_old = auc_pr_from_rank(ro$keys, truth, length(testable)),
        aucpr_new = auc_pr_from_rank(rn$keys, truth, length(testable)),
        r_precision_old = rpo$r_precision, r_precision_new = rpn$r_precision)

    rows <- vector("list", nrow(GRID))
    for (g in seq_len(nrow(GRID))) {
        adj <- GRID$adjust[g]
        val <- if (adj == "top") ceiling(GRID$value[g] * nrow(n$pv)) else GRID$value[g]
        co  <- call_by_threshold(o$pv, o$trait_cols, adj, val, quiet = TRUE)
        cn  <- call_by_threshold(n$pv, n$trait_cols, adj, val, quiet = TRUE)
        rem <- setdiff(cn$called, co$called); add <- setdiff(co$called, cn$called)
        if (adj %in% NESTED && length(add))                                              # G3
            stop(sprintf("G3 FAILED seed %s rung %s: %d SNP(s) called by pmax but not by one-fit under a fixed cutoff",
                         seed, GRID$rung[g], length(add)))
        so <- score_calls(co$called, 0, truth, testable)
        sn <- score_calls(cn$called, 0, truth, testable)
        io <- co$info[[1]]; inn <- cn$info[[1]]
        # as.data.table(c(list, list)): data.table(..., as.list(x)) would NOT spread x into
        # columns, it makes one list column
        rows[[g]] <- as.data.table(c(list(seed = seed, adjust = adj, rung = GRID$rung[g],
            value_applied = as.numeric(val),
            status_old = io$status, status_new = inn$status,
            cutoff_old = io$threshold, cutoff_new = inn$threshold,
            n_old = length(co$called), n_new = length(cn$called),
            n_both = length(intersect(co$called, cn$called)),
            tp_old = so$tp, tp_new = sn$tp,
            fp_old = so$fp_background, fp_new = sn$fp_background,
            removed_n = length(rem), added_n = length(add)),
            as.list(setNames(classify(rem, truth), paste0("removed_", c("causal", "background", "linked", "untracked")))),
            as.list(setNames(classify(add, truth), paste0("added_", c("causal", "background", "linked", "untracked"))))))
    }
    list(bind = bind, rungs = rbindlist(rows))
}

# ------------------------------------------------------------------- the corpus
man  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(man)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))
SEEDS <- sort(PRIM$seed)
message(sprintf("INFO: %d replicate(s), %d rungs, TOL %g, %d cores, corpus [%s]",
                length(SEEDS), nrow(GRID), TOL, NCOR, mvp_arm_label()))

t0  <- Sys.time()
res <- mclapply(SEEDS, function(s)
    tryCatch(guard_seed(s),
             error = function(e) list(error = data.table(seed = s, error = conditionMessage(e)))),
    mc.cores = NCOR, mc.preschedule = FALSE)
message(sprintf("INFO: scored in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

bad <- vapply(res, function(r) !is.list(r) || !is.null(r$error), logical(1))
if (any(bad)) {
    errs <- rbindlist(lapply(res[bad], function(r)
        if (is.list(r)) r$error else data.table(seed = NA, error = as.character(r))))
    print(errs)
    fwrite(errs, file.path(OUT, "errors.tsv"), sep = "\t")
    stop("G5 FAILED: ", nrow(errs), " replicate(s) errored (see errors.tsv)")
}
B  <- rbindlist(lapply(res, `[[`, "bind"))
RG <- rbindlist(lapply(res, `[[`, "rungs"))
stopifnot(nrow(B) == length(SEEDS), uniqueN(B$seed) == length(SEEDS),
          nrow(RG) == nrow(GRID) * length(SEEDS))                                        # G5

# ------------------------------------------------------------------ G4
g4 <- function(dir, arm) {
    M <- fread(file.path(dir, "calls_per_method.tsv"), colClasses = c(seed = "character"))[method == "RDA"]
    R <- fread(file.path(dir, "rank_metrics.tsv"), colClasses = c(seed = "character"))[method == "RDA"]
    x <- merge(RG[, .(seed, rung, n = get(paste0("n_", arm)), tp = get(paste0("tp_", arm)),
                      fp = get(paste0("fp_", arm)))],
               M[, .(seed, rung, n_ref = n_called, tp_ref = tp, fp_ref = fp_background)],
               by = c("seed", "rung"), all.x = TRUE)
    if (anyNA(x$n_ref) || any(x$n != x$n_ref | x$tp != x$tp_ref | x$fp != x$fp_ref))
        stop(sprintf("G4 FAILED (%s): %d rung row(s) disagree with %s/calls_per_method.tsv",
                     arm, sum(is.na(x$n_ref) | x$n != x$n_ref | x$tp != x$tp_ref | x$fp != x$fp_ref), dir))
    a <- merge(B[, .(seed, aucpr = get(paste0("aucpr_", arm)))], R[, .(seed, ref = aucpr)], by = "seed", all.x = TRUE)
    if (anyNA(a$ref) || any(abs(a$aucpr - a$ref) > 1e-12))
        stop(sprintf("G4 FAILED (%s): AUC-PR disagrees with %s/rank_metrics.tsv", arm, dir))
    data.table(arm = arm, remeasure_dir = dir, rung_rows = nrow(x), aucpr_rows = nrow(a),
               max_abs_aucpr_diff = max(abs(a$aucpr - a$ref)))
}
AS <- rbind(g4(RM_OLD, "old"), g4(RM_NEW, "new"))
message("G4 PASSED: both arms reproduce the re-measure driver's RDA rows and AUC-PR")

# ------------------------------------------------------------------ write
w <- function(d, f) { fwrite(d, file.path(OUT, f), sep = "\t"); message("wrote ", f, " (", nrow(d), " rows)") }
w(B, "guard_binding.tsv"); w(RG, "guard_rungs.tsv"); w(AS, "assert_guard.tsv")
writeLines(c(sprintf("old_dir\t%s", OLD), sprintf("new_dir\t%s", NEW),
             sprintf("old_remeasure\t%s", RM_OLD), sprintf("new_remeasure\t%s", RM_NEW),
             sprintf("cell\t%s", CELL), sprintf("tol\t%g", TOL),
             sprintf("corpus\t%s", mvp_arm_label()), sprintf("n_replicates\t%d", length(SEEDS)),
             sprintf("run_at\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))),
           file.path(OUT, "provenance.tsv"))
message(sprintf("OK: %d replicates -> %s", length(SEEDS), OUT))
