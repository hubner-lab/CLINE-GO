#!/usr/bin/env Rscript
# =============================================================================
# mvp_rda_unc_check.R -- is the RDA-uncorrected arm (benchmarks/mvp_rda_uncorrected_run.sh) what it
# claims to be, and how does it relate to the two RDA definitions already on disk?
#
# Three arms of the same SNPs, same image, same partial-fit inputs:
#   corrected    params_onefit/MVP{seed}/c1/RDA_pvalues.tsv   p_corr = partial RDA, Condition(PC1..PCk)
#   uncorrected  params_rdaunc/MVP{seed}/c1/RDA_pvalues.tsv   p_unc  = RDA without Condition()
#   pmax         params/MVP{seed}/c1/RDA_pvalues.tsv          pmax_old = max(p_partial, p_unconstrained),
#                                                             the frozen arm (read-only)
#
# U1 (hard, stops the run): DONE + the four tables present; diagnostics say condition_pcs = 0,
#    structure_proxy "none ...", rda_axes = 2; SNP keys and NA pattern identical to the corrected table.
# U2 (reconstruction witness): relative deviation of pmax(p_corr, p_unc) from pmax_old, and how many
#    SNPs the two disagree on about WHICH fit gave the larger p. The deleted rda.R section 14 ran the
#    second rdadapt() with K pinned to the partial fit's K; this run selects K itself (2 = rank on
#    this corpus either way) and computes its own GIF, which the frozen arm never recorded. So U2 is
#    MEASURED by default and asserted only with U2_ASSERT=1 (Phase 3b plan: promote it if the pilot
#    seed is exact). Tolerance TOL (1e-9, relative) -- block ssclines_ncline_ctredge's pmax arm
#    carries ~1e-13 thread noise (-c4 vs -c2).
#
# Writes OUT_DIR/check.tsv (one row per seed) and OUT_DIR/binding.tsv (seed x truth class: where the
# uncorrected fit gave the larger p). Exits 1 on a U1 failure, or on U2 when asserted, AFTER writing.
#
#   SEEDS  comma list (e.g. the pilot 1231418), or unset = the corpus via benchmarks/mvp_arm.R
#          (MVP_ARM=primary_ssclines MVP_ADDED=<5 blocks> MVP_N_EXPECT=600)
#   CORR_DIR params_onefit | UNC_DIR params_rdaunc | OLD_DIR params | OUT_DIR remeasure600/rdaunc_check
#   TOL 1e-9 | U2_ASSERT 0 | NCORES 16
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
CORR <- Sys.getenv("CORR_DIR", file.path(EVAL, "params_onefit"))
UNC  <- Sys.getenv("UNC_DIR",  file.path(EVAL, "params_rdaunc"))
OLD  <- Sys.getenv("OLD_DIR",  file.path(EVAL, "params"))
OUT  <- Sys.getenv("OUT_DIR",  file.path(EVAL, "remeasure600", "rdaunc_check"))
TOL  <- as.numeric(Sys.getenv("TOL", "1e-9"))
U2_ASSERT <- identical(Sys.getenv("U2_ASSERT", "0"), "1")
NCOR <- as.integer(Sys.getenv("NCORES", "16"))
source(file.path(ROOT, "benchmarks/lib_detection.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

check_seed <- function(seed) {
    d <- file.path(UNC, paste0("MVP", seed), "c1")
    need <- file.path(d, c("DONE", "RDA_pvalues.tsv", "RDA_candidates.tsv", "RDA_diagnostics.tsv", "RDA_anova.tsv"))
    u1 <- character(0)
    miss <- basename(need)[!file.exists(need)]
    if (length(miss)) return(list(row = data.table(seed = seed, u1_ok = FALSE,
                                                   u1_problem = paste("missing:", paste(miss, collapse = ","))), bind = NULL))
    dg <- fread(need[4], colClasses = "character", sep = "\t")
    kv <- setNames(dg$value, dg$key)
    if (!identical(unname(kv["condition_pcs"]), "0")) u1 <- c(u1, paste0("condition_pcs=", kv["condition_pcs"]))
    if (!startsWith(unname(kv["structure_proxy"]) %||% "", "none")) u1 <- c(u1, "structure_proxy not 'none'")
    if (!identical(unname(kv["rda_axes"]), "2")) u1 <- c(u1, paste0("rda_axes=", kv["rda_axes"]))

    rd <- function(f) { lp <- load_pvalues(f, "all")
        if (length(lp$trait_cols) != 1L) stop("RDA table with ", length(lp$trait_cols), " p columns: ", f); lp }
    lu <- rd(need[2])
    lc <- rd(file.path(CORR, paste0("MVP", seed), "c1", "RDA_pvalues.tsv"))
    lo <- rd(file.path(OLD,  paste0("MVP", seed), "c1", "RDA_pvalues.tsv"))
    if (!identical(lu$pv$key, lc$pv$key)) u1 <- c(u1, "SNP keys/order differ from corrected table")
    pu <- lu$pv[[lu$trait_cols]]; pc <- lc$pv[[lc$trait_cols]]; po <- lo$pv[[lo$trait_cols]]
    if (!identical(is.na(pu), is.na(pc))) u1 <- c(u1, "NA pattern differs from corrected table")
    if (!identical(lo$pv$key, lc$pv$key)) u1 <- c(u1, "frozen pmax table keys differ")
    if (length(u1)) return(list(row = data.table(seed = seed, u1_ok = FALSE, u1_problem = paste(u1, collapse = "; ")),
                                bind = NULL))

    ok  <- !is.na(pc)
    rec <- pmax(pc, pu)
    rel <- abs(rec[ok] - po[ok]) / po[ok]
    unc_larger_new <- ok & pu > pc * (1 + TOL)       # this run: uncorrected fit gives the larger p
    unc_larger_old <- ok & po > pc * (1 + TOL)       # frozen arm: unconstrained fit was binding
    truth <- load_truth(file.path(ROOT, "data/mvp", paste0("MVP", seed), "truth_any.tsv"))
    keys <- lc$pv$key
    cls <- fifelse(keys %in% truth$causal$key, "causal",
           fifelse(keys %in% truth$bg_keys, "background",
           fifelse(keys %in% truth$linked_keys, "linked", "untracked")))
    bind <- data.table(cls = cls[ok], unc = unc_larger_new[ok], lr = log10(pu[ok]) - log10(pc[ok]))[
        , .(n = .N, unc_larger = sum(unc), median_log10_unc_over_corr = median(lr)), by = cls][, seed := seed]
    list(row = data.table(seed = seed, u1_ok = TRUE, u1_problem = "",
            n_snps = sum(ok), gif_lambda_unc = as.numeric(kv["gif_lambda"]),
            adj_r2_unc = as.numeric(kv["adj_r_squared"]),
            median_p_unc = median(pu[ok]), median_p_corr = median(pc[ok]),
            u2_max_rel_dev = max(rel), u2_n_above_tol = sum(rel > TOL),
            u2_binding_disagree = sum(unc_larger_new != unc_larger_old),
            share_unc_larger = mean(unc_larger_new[ok]), share_unc_binding_old = mean(unc_larger_old[ok])),
         bind = bind)
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

if (nzchar(Sys.getenv("SEEDS"))) {
    SEEDS <- trimws(strsplit(Sys.getenv("SEEDS"), ",")[[1]])
} else {
    source(file.path(ROOT, "benchmarks/mvp_arm.R"))
    man <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
    PRIM <- mvp_prim(man)
    stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))
    SEEDS <- sort(PRIM$seed)
}
message(sprintf("INFO: %d seed(s), TOL %g, U2 %s", length(SEEDS), TOL, if (U2_ASSERT) "ASSERTED" else "measured"))

res <- mclapply(SEEDS, function(s) tryCatch(check_seed(s), error = function(e)
    list(row = data.table(seed = s, u1_ok = FALSE, u1_problem = paste("error:", conditionMessage(e))), bind = NULL)),
    mc.cores = NCOR)
CK <- rbindlist(lapply(res, `[[`, "row"), fill = TRUE)
BD <- rbindlist(lapply(res, `[[`, "bind"))
fwrite(CK, file.path(OUT, "check.tsv"), sep = "\t")
if (nrow(BD)) fwrite(setcolorder(BD, "seed"), file.path(OUT, "binding.tsv"), sep = "\t")
message(sprintf("wrote check.tsv (%d rows), binding.tsv (%d rows)", nrow(CK), nrow(BD)))

bad <- CK[u1_ok == FALSE]
if (nrow(bad)) { print(bad[, .(seed, u1_problem)]); message("U1 FAILED on ", nrow(bad), " seed(s)") }
good <- CK[u1_ok == TRUE]
if (nrow(good)) message(sprintf("U2: max rel dev %.3g over %d seed(s); SNPs above TOL %d; binding disagreements %d",
                                max(good$u2_max_rel_dev), nrow(good), sum(good$u2_n_above_tol), sum(good$u2_binding_disagree)))
u2_fail <- U2_ASSERT && nrow(good) && any(good$u2_n_above_tol > 0)
if (u2_fail) message("U2 FAILED (asserted): pmax(p_corr, p_unc) does not reproduce the frozen pmax within TOL")
if (nrow(bad) || u2_fail) quit(status = 1)
message("OK: U1 passed", if (U2_ASSERT) " and U2 passed (asserted)" else "; U2 measured")
