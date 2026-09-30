#!/usr/bin/env Rscript
# mvp_onefit_pilot_check.R -- did the one-fit RDA re-run change ONLY what it was meant to?
#
#   Rscript benchmarks/mvp_onefit_pilot_check.R --seeds=1231418,1231310,1231131 \
#       [--ref=benchmarks/mvp_eval/params] [--new=benchmarks/mvp_eval/params_onefit] \
#       [--anchor=1231418:benchmarks/mvp_eval/onefit_pilot/pmax_diag_1231418] \
#       [--runlogs=benchmarks/mvp_eval/onefit_pilot/runlogs] \
#       [--out=benchmarks/mvp_eval/onefit_pilot]
#
# SS-Clines re-analysis, Phase 2a step 7 (plan: ~/.claude/plans/we-are-now-starting-radiant-
# pillow.md). LFMM / EMMAX byte-identity is NOT checked here -- run mvp_repro_compare.R with
# --new=benchmarks/mvp_eval/params_onefit --methods=LFMM,EMMAX for that. This script asserts
# the RDA side and the harvest:
#
#   A1 harvest complete: 3 p-value tables + RDA {candidates,diagnostics,anova} + LD decay +
#      pca.{eigenvalues,tracywidom}, and no BLINK table.
#   A2 RDA p_new <= pmax_old on EVERY SNP, same keys in the same order. The frozen arm wrote
#      pmax(p_partial, p_unconstrained); the partial fit is unchanged, so p_new = p_partial.
#      Compared with a relative tolerance (--tol, default 1e-9): block ssclines_ncline_ctredge
#      was frozen with snakemake -c4 (runbook default) while every other block and this re-run
#      used -c2, and its RDA p differs from the -c2 fit at ~1e-13 relative (max 6.3e-13 over
#      its 120 seeds) with LFMM/EMMAX still byte-identical. n_violations_exact keeps the bitwise
#      count so that difference stays visible.
#   A3 the removed keys / columns are gone: no *_partial / *_unconstrained diagnostics keys,
#      no split candidate columns, no anova `model` column -- i.e. the NEW rda.R ran.
#   A4 anchor seed (exact): its pmax-era side tables, stashed before the re-run, carry the
#      partial fit's own scalars. gif_lambda must equal gif_lambda_partial, rda_axes must equal
#      k_partial, and adj_r_squared / n_markers_fitted / max_vif / qvalue_method must match;
#      the anova table must equal the pmax table's model == "partial" rows. A2 alone would
#      pass even if the partial fit had moved; this is the witness that it did not.
#
# Writes pilot_check.tsv (one row per seed: call counts before/after on the pipeline's own
# thresholds, binding-fit split) and pilot_walltime.tsv (per seed x mode, from the runlogs).
# Exits 1 on any failed assertion, AFTER writing both tables.

suppressPackageStartupMessages(library(data.table))

A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
    a <- sub("^--", "", a)
    A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}
ROOT    <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
opt     <- function(k, d) if (is.null(A[[k]])) d else A[[k]]
SEEDS   <- strsplit(opt("seeds", stop("--seeds is required")), ",", fixed = TRUE)[[1]]
REF     <- file.path(ROOT, opt("ref", "benchmarks/mvp_eval/params"))
NEW     <- file.path(ROOT, opt("new", "benchmarks/mvp_eval/params_onefit"))
ANCHOR  <- strsplit(opt("anchor", "1231418:benchmarks/mvp_eval/onefit_pilot/pmax_diag_1231418"),
                    ":", fixed = TRUE)[[1]]
RUNLOGS <- file.path(ROOT, opt("runlogs", "benchmarks/mvp_eval/onefit_pilot/runlogs"))
OUT     <- file.path(ROOT, opt("out", "benchmarks/mvp_eval/onefit_pilot"))
TOL     <- as.numeric(opt("tol", "1e-9"))

source(file.path(ROOT, "scripts/R/utils/pval_threshold.R"))  # compute_pval_threshold()
suppressPackageStartupMessages(library(qvalue))                # its qval branch needs it

HARVEST <- c("EMMAX_pvalues.tsv", "LFMM_pvalues.tsv", "RDA_pvalues.tsv",
             "RDA_candidates.tsv", "RDA_diagnostics.tsv", "RDA_anova.tsv",
             "ld_decay_half_distances.tsv", "pca.eigenvalues", "pca.tracywidom")
REMOVED_KEYS <- c("gif_lambda_partial", "gif_lambda_unconstrained", "k_partial",
                  "k_unconstrained", "k_pin_capped", "unconstrained_fit_status")
REMOVED_COLS <- c("mahalanobis_partial", "mahalanobis_unconstrained",
                  "p_value_partial", "p_value_unconstrained")

fails <- character(0)
check <- function(ok, msg) if (!isTRUE(ok)) fails <<- c(fails, msg)

read_diag <- function(f) {
    d <- fread(f, sep = "\t", colClasses = "character")
    setNames(d$value, d$key)
}
n_called <- function(p, adjust, value) {
    r <- compute_pval_threshold(p, adjust, value)
    if (!identical(r$status, "ok")) return(0L)
    sum(p <= r$threshold, na.rm = TRUE)            # inclusive, as sig_snps.R
}

check_seed <- function(seed) {
    dn <- file.path(NEW, paste0("MVP", seed), "c1")
    dr <- file.path(REF, paste0("MVP", seed), "c1")
    tag <- paste0("MVP", seed)

    # ---- A1
    miss <- HARVEST[!file.exists(file.path(dn, HARVEST)) | file.size(file.path(dn, HARVEST)) == 0]
    check(length(miss) == 0, paste0(tag, " A1 missing/empty: ", paste(miss, collapse = ",")))
    check(!file.exists(file.path(dn, "BLINK_pvalues.tsv")), paste0(tag, " A1 BLINK harvested"))
    if (length(miss) > 0) return(data.table(seed = seed, harvest_ok = FALSE))

    # ---- A2
    po <- fread(file.path(dr, "RDA_pvalues.tsv"), colClasses = c(chr = "character"))
    pn <- fread(file.path(dn, "RDA_pvalues.tsv"), colClasses = c(chr = "character"))
    same_keys <- identical(po$SNPID, pn$SNPID)
    check(same_keys, paste0(tag, " A2 RDA SNP keys differ from the frozen table"))
    old <- po$climate_multivariate; new <- pn$climate_multivariate
    check(identical(is.na(old), is.na(new)), paste0(tag, " A2 NA pattern differs"))
    n_viol_exact <- sum(new > old, na.rm = TRUE)
    n_viol <- sum(new > old * (1 + TOL), na.rm = TRUE)
    max_rel_excess <- if (n_viol_exact > 0) max(((new - old) / old)[which(new > old)]) else 0
    check(n_viol == 0, paste0(tag, " A2 ", n_viol, " SNPs with p_new > pmax_old * (1 + ", TOL, ")"))

    # ---- A3
    dg <- read_diag(file.path(dn, "RDA_diagnostics.tsv"))
    left <- intersect(REMOVED_KEYS, names(dg))
    check(length(left) == 0, paste0(tag, " A3 removed diagnostics keys present: ", paste(left, collapse = ",")))
    cand <- names(fread(file.path(dn, "RDA_candidates.tsv"), nrows = 0))
    check(length(intersect(REMOVED_COLS, cand)) == 0, paste0(tag, " A3 split candidate columns present"))
    an <- fread(file.path(dn, "RDA_anova.tsv"))
    check(!"model" %in% names(an), paste0(tag, " A3 anova still has a model column"))
    check(as.numeric(dg[["condition_pcs"]]) > 0, paste0(tag, " A3 condition_pcs is 0 -- the pilot tests nothing"))

    # ---- A4
    anchored <- identical(seed, ANCHOR[1])
    if (anchored) {
        ad <- file.path(ROOT, ANCHOR[2])
        dp <- read_diag(file.path(ad, "RDA_diagnostics.tsv"))
        eq <- function(new_key, old_key = new_key)
            check(identical(dg[[new_key]], dp[[old_key]]),
                  paste0(tag, " A4 ", new_key, " = ", dg[[new_key]], " but pmax-era ", old_key,
                         " = ", dp[[old_key]]))
        eq("gif_lambda", "gif_lambda_partial"); eq("rda_axes", "k_partial")
        for (k in c("adj_r_squared", "n_markers_fitted", "max_vif", "qvalue_method",
                    "rda_axes_max", "condition_pcs")) eq(k)
        ap <- fread(file.path(ad, "RDA_anova.tsv"))
        ap <- ap[model == "partial"][, model := NULL]
        check(isTRUE(all.equal(as.data.frame(ap), as.data.frame(an), check.attributes = FALSE)),
              paste0(tag, " A4 anova differs from the pmax-era partial rows"))
    }

    m <- sum(!is.na(new))
    data.table(
        seed = seed, harvest_ok = TRUE, anchored = anchored, n_snps = m,
        k_best = dg[["k_best"]], condition_pcs = dg[["condition_pcs"]], rda_axes = dg[["rda_axes"]],
        gif_lambda = as.numeric(dg[["gif_lambda"]]), n_violations = n_viol,
        n_violations_exact = n_viol_exact, max_rel_excess = max_rel_excess,
        n_equal = sum(new == old, na.rm = TRUE),            # partial fit was the binding one
        n_lower = sum(new < old, na.rm = TRUE),             # unconstrained fit was binding
        median_p_old = median(old, na.rm = TRUE), median_p_new = median(new, na.rm = TRUE),
        min_p_old = min(old, na.rm = TRUE), min_p_new = min(new, na.rm = TRUE),
        bonf005_old = n_called(old, "bonf", 0.05), bonf005_new = n_called(new, "bonf", 0.05),
        qval01_old  = n_called(old, "qval", 0.1),  qval01_new  = n_called(new, "qval", 0.1),
        top05pct_new = n_called(new, "top", ceiling(0.005 * m)),
        applied_rule = dg[["candidate_rule"]], n_candidates_total = dg[["n_candidates_total"]])
}

res <- rbindlist(lapply(SEEDS, check_seed), fill = TRUE)

# ---- wall time per mode, from the driver's own timestamps
wall <- rbindlist(lapply(SEEDS, function(seed) {
    f <- file.path(RUNLOGS, paste0("MVP", seed, ".log"))
    if (!file.exists(f)) return(NULL)
    L <- readLines(f)
    L <- L[grepl(paste0("^\\[MVP", seed, "\\] "), L)]
    ts <- regmatches(L, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9:]{8}[+-][0-9:]{5}", L))
    has <- grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}T", L)
    ev  <- sub("^\\[MVP[0-9]+\\] ", "", L[has])
    ev  <- trimws(sub("[0-9]{4}-[0-9]{2}-[0-9]{2}T.*$", "", ev))
    # the interval each timestamp OPENS: "k_best=3  cells=1" -> seed setup, "c1" -> mode=gea,
    # "c1 gea finished" -> harvest
    ev  <- sub("^k_best=.*$", "setup", ev)
    ev  <- sub("^mode=", "", ev)
    ev  <- sub("^c1$", "gea", ev)
    ev  <- sub("^c1 gea finished$", "harvest", ev)
    t   <- as.POSIXct(sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", ts), format = "%Y-%m-%dT%H:%M:%S%z")
    if (length(t) < 2) return(NULL)
    data.table(seed = seed, step = ev[-length(ev)],
               minutes = round(as.numeric(diff(t), units = "mins"), 1))
}))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fwrite(res,  file.path(OUT, "pilot_check.tsv"), sep = "\t")
fwrite(wall, file.path(OUT, "pilot_walltime.tsv"), sep = "\t")
print(res); print(wall)
if (length(fails) > 0) {
    message("PILOT CHECK FAILED (", length(fails), "):\n  ", paste(fails, collapse = "\n  "))
    quit(status = 1)
}
message("PILOT CHECK PASSED: ", length(SEEDS), " seeds, A1-A4")
