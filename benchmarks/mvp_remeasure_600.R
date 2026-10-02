#!/usr/bin/env Rscript
# =============================================================================
# mvp_remeasure_600.R -- re-measure every GEA analysis parameter on the SS-Clines corpus.
#
# WHY. Every operating point, rule type and agreement window the simulation arm uses was
# inherited from one SS-Mtn anchor replicate and never measured on SS-Clines
# (docs/gea-simulation-reanalysis.md, Phase 1). The harvested p-value tables of all 600
# replicates survived the 2026-09-23 cleanup, so the whole ladder can be measured post hoc,
# with no pipeline compute. This script writes the measurements; journal 17 reads them.
#
# WHAT IT MEASURES, per replicate x method (LFMM, EMMAX, RDA):
#   calls_per_trait.tsv   every grid rung x trait: realised N (top), cutoff, rule STATUS
#                         (ok / no_hits / too_few_tests / engine_error -- never conflated),
#                         calls on that trait
#   calls_per_method.tsv  every grid rung: SNPs called (union over traits), P/R/F1 vs truth
#   calls_combine.tsv     every grid rung x agreement window: union / >=2 of 3 / 3 of 3,
#                         sizes and P/R/F1
#   rank_metrics.tsv      AUC-PR, R-precision (threshold-free)
#   pr_curves.tsv         the AUC-PR step points (one row per retrieved causal locus)
#   lambda_pi0.tsv        genomic-control lambda and qvalue pi0 per method x trait
#   rda_{ARM}_dist.tsv    the RDA p distribution of the arm read. ARM=pmax (the frozen arm): the
#                         combined p, pmax of the partial and unconstrained fit -- the "before"
#                         baseline. ARM=onefit (Phase 2): the one fit's rdadapt p -- the "after".
#   covariates.tsv        design cell, demography, K_authors, k_best, manifest structure stats
#   assert_*.tsv          what the assertions below checked
#
# GRID (per method, applied PER TRAIT and unioned over traits, as the pipeline does):
#   bonf   0.01 0.05 0.1 0.5 1          qval  0.01 0.02 0.05 0.1 0.15 0.2
#   top    0.1 0.25 0.5 1 2 % of tested SNPs, N = ceiling(share x n_snps) -- the corpus spans
#          5,141-32,130 SNPs, so a fixed N would not be a fixed share. N is recorded.
#   custom 1e-6 1e-5 1e-4 1e-3
# COMBINE: union / >=2 of 3 / 3 of 3, agreement windows 0 / 1 / 2.5 / 5 / 10 kb; every method
# at the SAME rung (homogeneous). >=2 of 3 is the pipeline's Cross-method strategy; 3 of 3 is
# not a pipeline strategy and is measured as the strictest reference.
#
# SCORING CONVENTION (lib_detection.R, as journals 07 and 16): truth = truth_any.tsv; exact
# position match (window 0); causal = TP, background_neutral = the only FP, linked_neutral
# neither; recall against causal loci present in the tested set. No genomic control -- the
# tables are what production mode=gea wrote. Threshold calls use `p <= cutoff` (the 2026-09-28
# harness fix), so threshold counts are NOT one-to-one comparable with journal 16's.
#
# ASSERTIONS -- each stops the run, none warns:
#   A1 window 0: combine_support() equals the exact-key tally at EVERY rung, all three rules
#      (sweep_thresholds.R:186-203 checked only the first rung).
#   A2 >=2 of 3 at every window and rung equals the pipeline's own
#      scripts/R/lib/combine_sigsnps.R .overlap_cross_method() SNP set, on 12 replicates (one per
#      genic-level x architecture design cell).
#   A3 AUC-PR equals benchmarks/mvp_eval/detection600/aucpr_per_seed.tsv to 1e-12 on every row
#      of the A3_METHODS: same tables, same function, and rank-based metrics are untouched by the
#      `<=` fix. detection600 was scored on the FROZEN pmax arm, so A3 can only assert a method
#      whose p-table is unchanged in the arm being read: all three on ARM=pmax, LFMM and EMMAX on
#      ARM=onefit (byte-identical to the frozen arm, Phase 2b gate). One-fit RDA's AUC-PR is
#      expected to differ -- that difference is a measured result, written to the same table
#      with asserted = FALSE, never a relaxed assertion.
#   A4 600 replicates x 3 methods, no replicate missing, no scoring error.
# journal 16's legacy-seed reproduction gate (mvp_detection_600.R) does NOT apply here: the
# `<=` fix deliberately changes threshold counts, so reproducing journal 07's would be a failure.
#
#   PARAMS_DIR  default benchmarks/mvp_eval/params         (the frozen pmax arm, read-only)
#   ARM         default pmax -- label of the arm read; names OUT_DIR's default and the RDA file
#   OUT_DIR     default benchmarks/mvp_eval/remeasure600/{ARM}
#   A3_METHODS  default LFMM,EMMAX,RDA -- methods A3 asserts (ARM=onefit: LFMM,EMMAX)
#   TOP_SHARES  default unset = the 20-rung GRID below. A comma list of shares (e.g.
#               0.0015,0.002,0.0025) replaces the GRID with those `top` rungs only -- the Phase 3b
#               detection ledger (plan ~/.claude/plans/glittery-fluttering-seahorse.md).
#   WINDOWS_KB  default unset = 0,1,2.5,5,10. A comma list replaces it; 0 is always added, because
#               A1 (window 0 == exact-key tally) is checked there and nowhere else.
#   SHARE_MULT  default unset = 1 for every method. `METHOD=x[,METHOD=y]` multiplies that method's
#               `top` share only (rung names keep the base share). Phase 4a (2026-10-01): RDA=2 --
#               RDA tests both predictors in ONE multivariate p, so it gets the same total budget
#               as LFMM/EMMAX's 0.25 % per predictor x 2 predictors (plan
#               ~/.claude/plans/we-are-now-starting-radiant-pillow.md). Non-`top` rungs ignore it.
#   CELL        default c1
#   NCORES      default 16
#   N_SEEDS     default all (a smaller number runs the first N replicates -- timing only;
#               A3 and A4 are then checked on that subset)
#   MVP_ARM / MVP_ADDED / MVP_N_EXPECT   via benchmarks/mvp_arm.R (mandatory for SS-Clines)
# =============================================================================
# dplyr first so data.table keeps between()/first()/last(); dplyr is only for combine_sigsnps.R
suppressPackageStartupMessages({ library(dplyr); library(data.table); library(parallel) })

ROOT   <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL   <- file.path(ROOT, "benchmarks/mvp_eval")
PDIR   <- Sys.getenv("PARAMS_DIR", file.path(EVAL, "params"))
ARM    <- Sys.getenv("ARM", "pmax")
if (!grepl("^[a-z0-9_]+$", ARM)) stop("ARM must match ^[a-z0-9_]+$, got: ", ARM)
OUT    <- Sys.getenv("OUT_DIR", file.path(EVAL, "remeasure600", ARM))
CELL   <- Sys.getenv("CELL", "c1")
NCOR   <- as.integer(Sys.getenv("NCORES", "16"))
NSEEDS <- Sys.getenv("N_SEEDS", "all")
METHODS <- c("LFMM", "EMMAX", "RDA")
A3_METHODS <- trimws(strsplit(Sys.getenv("A3_METHODS", paste(METHODS, collapse = ",")), ",")[[1]])
if (!length(A3_METHODS) || !all(A3_METHODS %in% METHODS))
    stop("A3_METHODS must be a non-empty subset of ", paste(METHODS, collapse = ","))
WINDOWS_KB <- c(0, 1, 2.5, 5, 10)
SHARE_MULT <- setNames(rep(1, length(METHODS)), METHODS)
if (nzchar(Sys.getenv("SHARE_MULT"))) {
    kv <- strsplit(strsplit(Sys.getenv("SHARE_MULT"), ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
    if (any(lengths(kv) != 2L)) stop("SHARE_MULT must be METHOD=x[,METHOD=y]")
    mv <- setNames(as.numeric(vapply(kv, `[`, "", 2L)), vapply(kv, `[`, "", 1L))
    if (anyNA(mv) || any(mv <= 0) || !all(names(mv) %in% METHODS))
        stop("SHARE_MULT: unknown method or non-positive multiplier: ", Sys.getenv("SHARE_MULT"))
    SHARE_MULT[names(mv)] <- mv
}
if (nzchar(Sys.getenv("WINDOWS_KB"))) {
    WINDOWS_KB <- sort(unique(c(0, as.numeric(strsplit(Sys.getenv("WINDOWS_KB"), ",")[[1]]))))
    if (anyNA(WINDOWS_KB) || any(WINDOWS_KB < 0)) stop("WINDOWS_KB must be non-negative numbers")
}

source(file.path(ROOT, "scripts/R/utils/pval_threshold.R"))
source(file.path(ROOT, "benchmarks/lib_detection.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
source(file.path(ROOT, "scripts/R/lib/combine_sigsnps.R"))   # .overlap_cross_method(), for A2
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GRID <- rbind(
    data.table(adjust = "bonf",   value = c(0.01, 0.05, 0.1, 0.5, 1)),
    data.table(adjust = "qval",   value = c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2)),
    data.table(adjust = "top",    value = c(0.001, 0.0025, 0.005, 0.01, 0.02)),   # shares
    data.table(adjust = "custom", value = c(1e-6, 1e-5, 1e-4, 1e-3)))
if (nzchar(Sys.getenv("TOP_SHARES"))) {
    shares <- as.numeric(strsplit(Sys.getenv("TOP_SHARES"), ",")[[1]])
    if (anyNA(shares) || any(shares <= 0 | shares >= 1)) stop("TOP_SHARES must be shares in (0, 1)")
    GRID <- data.table(adjust = "top", value = sort(unique(shares)))
}
GRID[, rung := paste0(adjust, "_", vapply(value, format, "", scientific = FALSE, drop0trailing = TRUE))]
stopifnot(!anyDuplicated(GRID$rung))

# ------------------------------------------------------------------ one replicate
remeasure_seed <- function(seed, check_pipeline) {
    truth <- load_truth(file.path(ROOT, "data/mvp", paste0("MVP", seed), "truth_any.tsv"))
    L <- setNames(lapply(METHODS, function(m) {
        f <- file.path(PDIR, paste0("MVP", seed), CELL, paste0(m, "_pvalues.tsv"))
        if (!file.exists(f)) stop("missing p-value table: ", f)
        lp <- load_pvalues(f, "all")
        lp$testable <- testable_causal(lp$pv, lp$trait_cols, truth)
        lp$rank     <- rank_keys(lp$pv, lp$trait_cols)
        lp
    }), METHODS)
    testable_union <- unique(unlist(lapply(L, `[[`, "testable")))

    # ---- threshold-free
    rank_rows <- list(); pr_rows <- list(); lam_rows <- list()
    for (m in METHODS) {
        lp <- L[[m]]; nt <- length(lp$testable)
        rp <- r_precision_from_rank(lp$rank$keys, truth, nt)
        rank_rows[[m]] <- data.table(seed = seed, method = m, n_snps = nrow(lp$pv),
            n_traits = length(lp$trait_cols), n_testable = nt,
            aucpr = auc_pr_from_rank(lp$rank$keys, truth, nt),
            r_precision = rp$r_precision, r_precision_strict = rp$r_precision_strict)
        pr_rows[[m]] <- pr_curve_from_rank(lp$rank$keys, truth, nt)[, `:=`(seed = seed, method = m)]
        for (tc in lp$trait_cols) {
            p  <- lp$pv[[tc]]; pn <- p[!is.na(p)]
            q  <- tryCatch(list(pi0 = qvalue::qvalue(pn)$pi0, status = "ok"),
                           error = function(e) list(pi0 = NA_real_, status = conditionMessage(e)))
            lam_rows[[paste(m, tc)]] <- data.table(seed = seed, method = m, trait = tc,
                n_tested = length(pn), n_na = sum(is.na(p)), lambda_gc = gc_lambda(p),
                qvalue_pi0 = q$pi0, qvalue_status = q$status)
        }
    }

    # ---- per-method grid. call_by_threshold() is the scoring core's own mirror of the
    # pipeline's per-trait selection; its $info carries each trait's cutoff and STATUS.
    called <- list(); trait_rows <- list(); meth_rows <- list()
    for (g in seq_len(nrow(GRID))) {
        adj <- GRID$adjust[g]; rung <- GRID$rung[g]
        called[[rung]] <- list()
        for (m in METHODS) {
            lp  <- L[[m]]
            val <- if (adj == "top") ceiling(GRID$value[g] * SHARE_MULT[[m]] * nrow(lp$pv)) else GRID$value[g]
            res <- call_by_threshold(lp$pv, lp$trait_cols, adj, val, quiet = TRUE)
            called[[rung]][[m]] <- res$called
            for (tc in lp$trait_cols) {
                th <- res$info[[tc]]
                ntc <- if (identical(th$status, "ok") && !is.na(th$threshold))
                           sum(lp$pv[[tc]] <= th$threshold, na.rm = TRUE) else 0L
                trait_rows[[length(trait_rows) + 1L]] <- data.table(seed = seed, method = m,
                    trait = tc, adjust = adj, rung = rung, value_applied = as.numeric(val),
                    status = th$status, cutoff = th$threshold, n_called_trait = ntc)
            }
            s <- score_calls(res$called, 0, truth, lp$testable)
            meth_rows[[length(meth_rows) + 1L]] <- data.table(seed = seed, method = m,
                adjust = adj, rung = rung, value_applied = as.numeric(val),
                status = paste(vapply(lp$trait_cols, function(tc) res$info[[tc]]$status, ""),
                               collapse = ";"),
                n_called = s$n_called, tp = s$tp, tp_any_causal = s$tp_any_causal,
                expected_linked = s$expected_linked, fp_background = s$fp_background,
                precision_strict = s$precision_strict, recall_testable = s$recall_testable,
                f1 = s$f1, n_testable = length(lp$testable))
        }
    }

    # ---- combine, with A1 at every rung and A2 on the check replicates
    comb_rows <- list(); a1 <- 0L; a2 <- 0L
    for (rk in names(called)) {
        gi    <- match(rk, GRID$rung)
        parts <- called[[rk]]
        adj   <- GRID$adjust[gi]
        tab   <- table(unlist(parts, use.names = FALSE))
        for (w in WINDOWS_KB) {
            sup  <- combine_support(parts, w * 1000)
            sets <- list(union = sup$snp, ge2 = sup[n_methods >= 2, snp],
                         all3 = sup[n_methods >= length(METHODS), snp])
            if (w == 0) {                                                   # A1
                chk <- function(a, b, what) if (!setequal(a, b))
                    stop(sprintf("A1 FAILED seed %s rung %s %s: combine_support %d vs exact-key %d",
                                 seed, rk, what, length(a), length(b)))
                chk(sets$union, names(tab), "union")
                chk(sets$ge2,   names(tab)[tab >= 2], "ge2")
                chk(sets$all3,  names(tab)[tab == length(METHODS)], "all3")
                a1 <- a1 + 1L
            }
            if (check_pipeline) {                                           # A2
                sig <- setNames(lapply(METHODS, function(m) {
                    lp <- L[[m]]
                    rbindlist(lapply(lp$trait_cols, function(tc) {
                        th <- call_by_threshold(lp$pv, tc, adj,
                                  if (adj == "top") ceiling(GRID$value[gi] * SHARE_MULT[[m]] * nrow(lp$pv))
                                  else GRID$value[gi])$info[[tc]]
                        if (!identical(th$status, "ok") || is.na(th$threshold)) return(NULL)
                        hit <- !is.na(lp$pv[[tc]]) & lp$pv[[tc]] <= th$threshold
                        data.table(SNPID = lp$pv$key[hit], chr = lp$pv$chr[hit],
                                   pos = as.integer(lp$pv$pos[hit]), trait = tc, method = m,
                                   pvalue = lp$pv[[tc]][hit])
                    }))
                }), METHODS)
                sig <- lapply(sig, function(d) if (is.null(d) || !nrow(d))
                    data.table(SNPID = character(), chr = character(), pos = integer(),
                               trait = character(), method = character(), pvalue = numeric())
                    else d)
                pipe <- .overlap_cross_method(sig, METHODS, as.integer(w * 1000), per_trait = FALSE)
                pipe_set <- if (is.null(pipe)) character(0) else unique(pipe$SNPID)
                if (!setequal(pipe_set, sets$ge2))
                    stop(sprintf("A2 FAILED seed %s rung %s window %s kb: harness >=2 has %d SNPs, pipeline Cross-method %d",
                                 seed, rk, w, length(sets$ge2), length(pipe_set)))
                a2 <- a2 + 1L
            }
            for (cn in names(sets)) {
                s <- score_calls(sets[[cn]], 0, truth, testable_union)
                comb_rows[[length(comb_rows) + 1L]] <- data.table(seed = seed, adjust = adj,
                    rung = rk, window_kb = w, combine = cn, n_called = s$n_called,
                    tp = s$tp, tp_any_causal = s$tp_any_causal,
                    expected_linked = s$expected_linked, fp_background = s$fp_background,
                    precision_strict = s$precision_strict, recall_testable = s$recall_testable,
                    f1 = s$f1, n_testable = length(testable_union))
            }
        }
    }

    # ---- the RDA p of this arm, as harvested (pmax: pmax(p_partial, p_unconstrained); onefit: one fit)
    prda <- L$RDA$pv[[L$RDA$trait_cols[1]]]; prda <- prda[!is.na(prda)]
    qq <- quantile(prda, c(0.001, 0.01, 0.05, 0.25, 0.5, 0.75))
    rda_row <- data.table(seed = seed, n = length(prda), min_p = min(prda),
        q001 = qq[[1]], q01 = qq[[2]], q05 = qq[[3]], q25 = qq[[4]], median = qq[[5]],
        q75 = qq[[6]], share_lt_005 = mean(prda < 0.05), share_ge_099 = mean(prda >= 0.99),
        lambda_gc = gc_lambda(prda))

    list(rank = rbindlist(rank_rows), pr = rbindlist(pr_rows), lam = rbindlist(lam_rows),
         trait = rbindlist(trait_rows), meth = rbindlist(meth_rows), comb = rbindlist(comb_rows),
         rda = rda_row, assert = data.table(seed = seed, a1_rungs = a1, a2_checks = a2))
}

# ------------------------------------------------------------------- the corpus
man  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(man)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))

cov <- PRIM[, .(seed, block = sub("^ssclines_", "", added), arch_level, architecture,
                ispleiotropy, K_authors, k_best, meanFst, r2_pc1_temp, r2_pc1_sal,
                n_causal_maf01, final_LA)]
cov[, regime := fifelse(grepl("unequal-S", architecture), "unequal-S", "equal-S")]
cov[, pleiotropy := fifelse(ispleiotropy == 1L, "pleiotropy", "no pleiotropy")]
cov[, design_cell := paste(arch_level, pleiotropy, regime, sep = " | ")]
cov[, k_floored := K_authors < 2]
stopifnot(uniqueN(cov$design_cell) == 12L, all(cov$k_best == pmax(2L, cov$K_authors)))
setorder(cov, seed)

# One replicate per design cell for A2: the lowest seed id in the cell, i.e. chosen on
# design, never on outcome.
CHECK <- cov[, .(seed = min(seed)), by = design_cell]$seed
stopifnot(length(CHECK) == 12L)

SEEDS <- cov$seed
if (NSEEDS != "all") SEEDS <- unique(c(head(SEEDS, as.integer(NSEEDS))))
message(sprintf("INFO: %d replicate(s), %d with the pipeline cross-check, %d cores, arm [%s]",
                length(SEEDS), sum(SEEDS %in% CHECK), NCOR, mvp_arm_label()))

t0  <- Sys.time()
res <- mclapply(SEEDS, function(s)
    tryCatch(remeasure_seed(s, s %in% CHECK),
             error = function(e) list(error = data.table(seed = s, error = conditionMessage(e)))),
    mc.cores = NCOR, mc.preschedule = FALSE)
message(sprintf("INFO: scored in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

bad <- vapply(res, function(r) !is.list(r) || !is.null(r$error), logical(1))
if (any(bad)) {
    errs <- rbindlist(lapply(res[bad], function(r)
        if (is.list(r)) r$error else data.table(seed = NA, error = as.character(r))))
    print(errs)
    fwrite(errs, file.path(OUT, "errors.tsv"), sep = "\t")
    stop("A4 FAILED: ", nrow(errs), " replicate(s) errored (see errors.tsv)")
}
grab <- function(k) rbindlist(lapply(res, `[[`, k))
R  <- grab("rank"); PR <- grab("pr"); LAM <- grab("lam"); TR <- grab("trait")
M  <- grab("meth"); C  <- grab("comb"); RDA <- grab("rda"); AS <- grab("assert")

# ------------------------------------------------------------------ A3, A4
stopifnot(nrow(R) == length(METHODS) * length(SEEDS), uniqueN(R$seed) == length(SEEDS),
          !anyNA(R$aucpr),
          nrow(M) == nrow(GRID) * length(METHODS) * length(SEEDS),
          nrow(C) == nrow(GRID) * length(WINDOWS_KB) * 3L * length(SEEDS),
          all(AS$a1_rungs == nrow(GRID)),
          all(AS[seed %in% CHECK, a2_checks] == nrow(GRID) * length(WINDOWS_KB)))
ref <- fread(file.path(EVAL, "detection600", "aucpr_per_seed.tsv"), colClasses = c(seed = "character"))
a3  <- merge(R[, .(seed, method, aucpr, n_testable)],
             ref[, .(seed, method, ref_aucpr = aucpr, ref_testable = n_testable)],
             by = c("seed", "method"), all.x = TRUE)
a3[, abs_diff := abs(aucpr - ref_aucpr)]
a3[, asserted := method %in% A3_METHODS]
a3c <- a3[asserted == TRUE]
if (anyNA(a3c$ref_aucpr) || any(a3c$abs_diff > 1e-12) || any(a3c$n_testable != a3c$ref_testable))
    stop("A3 FAILED: AUC-PR does not reproduce detection600 (max |diff| = ",
         format(max(a3c$abs_diff, na.rm = TRUE)), ")")
message(sprintf("A3 PASSED: AUC-PR reproduces detection600 on %d rows (%s); %d rows measured, not asserted",
                nrow(a3c), paste(A3_METHODS, collapse = ","), sum(!a3$asserted)))

# ------------------------------------------------------------------ write
w <- function(d, f) { fwrite(d, file.path(OUT, f), sep = "\t"); message("wrote ", f, " (", nrow(d), " rows)") }
w(cov[seed %in% SEEDS], "covariates.tsv")
w(R, "rank_metrics.tsv");  w(PR, "pr_curves.tsv");      w(LAM, "lambda_pi0.tsv")
w(TR, "calls_per_trait.tsv"); w(M, "calls_per_method.tsv"); w(C, "calls_combine.tsv")
w(RDA, sprintf("rda_%s_dist.tsv", ARM))
w(AS[, .(seed, a1_rungs, a2_checks, check_replicate = seed %in% CHECK)], "assert_combine.tsv")
w(a3, "assert_aucpr_detection600.tsv")
w(GRID, "grid.tsv")
writeLines(c(sprintf("arm\t%s", ARM), sprintf("params_dir\t%s", PDIR), sprintf("cell\t%s", CELL),
             sprintf("corpus\t%s", mvp_arm_label()), sprintf("a3_methods\t%s", paste(A3_METHODS, collapse = ",")),
             sprintf("n_replicates\t%d", length(SEEDS)),
             sprintf("grid\t%s", paste(GRID$rung, collapse = ",")),
             sprintf("windows_kb\t%s", paste(WINDOWS_KB, collapse = ",")),
             sprintf("share_mult\t%s", paste(names(SHARE_MULT), SHARE_MULT, sep = "=", collapse = ",")),
             sprintf("check_replicates\t%s", paste(CHECK, collapse = ",")),
             sprintf("run_at\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))),
           file.path(OUT, "provenance.tsv"))
message(sprintf("OK: %d replicates -> %s", length(SEEDS), OUT))
