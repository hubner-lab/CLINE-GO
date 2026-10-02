#!/usr/bin/env Rscript
# =============================================================================
# mvp_build_snp_sets.R -- build the SNP panels that feed mode=maladaptation, per replicate.
#
# The offset arm asks ONE question: which SNPs should be handed to the offset model? This
# script writes the candidate answers, plus the references that put them on an absolute
# scale, into {PROJECT}_results/_intermediate/snp_sets/{name}/selected_snps.tsv -- the store
# the pipeline's resolve_active_snp_sets() (workflow/rules/common.smk) reads at parse time.
#
# [rewritten 2026-10-01, SS-Clines re-analysis Phase 4a -- plan
#  ~/.claude/plans/we-are-now-starting-radiant-pillow.md]. The previous version (git 97f37ce)
# built scheme S1 (LFMM top_100, RDA top_100, EMMAX custom 1e-04, pmax RDA), whose operating
# points came from one SS-Mtn anchor replicate. Every value below is instead chosen from
# measurements on the 600 SS-Clines replicates themselves (journals 17/18, the Phase 3 and
# Phase 3b gates, benchmarks/mvp_eval/figures_ssclines_rda/RDA_CORRECTION_BRIEF.md).
#
# PANELS
#
#   solo_lfmm    LFMM alone                               <- candidate
#   solo_rda     RDA (uncorrected) alone                  <- candidate
#   solo_emmax   EMMAX alone                              <- candidate
#   union        1/3: called by >= 1 method               <- candidate
#   best         2/3: a call with another method's call within 5 kb   <- candidate
#   intersect3   3/3: calls of all three methods within 5 kb          <- candidate
#   truth        the replicate's causal loci              <- CEILING (not a candidate)
#   all          every tested SNP                         <- SNP list only (PR denominator)
#   neutral_all  background-neutral loci (QTN-free LGs)   <- SNP list only
#
# Names are kept from the previous builder (`best` = 2/3, `intersect3` = 3/3) because every
# downstream scorer, table and gallery script keys on them.
#
# THE THREE METHODS AND THEIR P-VALUES. LFMM (K = k_best) and EMMAX (k_best PCs) are the
# one-fit re-run (params_onefit/). RDA is the UNCORRECTED fit (condition_pcs = 0,
# params_rdaunc/): the Phase 3b gate (2026-10-01) made uncorrected RDA the third method after
# it out-ranked the structure-corrected fit on 557/600 replicates (AUC-PR, brief section 3).
# params_rdaunc/MVP{seed}/c1/ holds exactly that trio -- its LFMM/EMMAX tables are symlinks
# into params_onefit/ -- and it is the directory remeasure600/rdaunc/ was scored from, which is
# what makes the regression tie below an identity rather than a coincidence. This script
# asserts the trio before writing anything; there is no fallback root (the old PARAMS_ANCHOR
# fallthrough could silently read the frozen pmax arm).
#
# OPERATING POINT: top 0.25 % per method, per predictor, unioned over predictors -- exactly as
# the pipeline's find_significant_snps_per_trait() applies a `top N` rule, with
# N = ceiling(share x tested SNPs of that method's table). The share -> N conversion is ours
# (the corpus spans 5,141-32,130 tested SNPs, so a fixed N is not a fixed share); the rule
# applied with that N is the pipeline's own `top`. LFMM and EMMAX test bio_1 and bio_2
# separately, so each calls up to 2 x N SNPs (median 127).
# [changed 2026-10-01, user, Phase 4a] RDA has ONE multivariate p for both predictors, so at
# 0.25 % it called half as many SNPs (median 64). It now gets the same total budget: share x
# number of predictors = 0.5 % (median 127). remeasure600/rdaunc_rda2x/ is scored at exactly
# this mixed operating point (mvp_remeasure_600.R SHARE_MULT=RDA=2).
#
# AGREEMENT WINDOW: 5 kb for every combine rule (= snp_clumping_distance in all 600 replicate
# configs, the distance the pipeline's own Cross-method overlap uses; ~3x the measured LD
# half-decay). `union` is window-independent by construction, so one combine_support() call at
# 5 kb yields all three rules. 2/3 is the pipeline's Cross-method strategy; 3/3 is not a
# pipeline strategy and is the strictest reference.
#
# TRUTH ENTERS ONLY THE `truth` CEILING AND THE EVALUATION COUNTS. No candidate panel's size
# or membership depends on the truth table (hard constraint 3 of
# docs/gea-simulation-reanalysis.md): the rule is a fixed share applied identically to every
# replicate, so it can be built on real data.
#
# NO RANDOM FLOORS. Dropped by the user at Phase 4a (2026-10-01). The rand_* and the legacy
# `solo` (byte-identical to solo_lfmm) panels of the previous build are removed from disk --
# a stale pmax-era panel left in place would be swept and scored as if it were current. They
# are deleted only after --archive names an existing archive of the old snp_sets trees.
#
# `truth` is a pure function of truth_any.tsv (byte-identical across the Phase 2b reconversion),
# so the rebuilt panel must equal the one on disk: write_set() never rewrites a file whose
# content is unchanged (its mtime is kept), and the run FAILS unless every replicate's truth
# panel came out `unchanged`. [changed 2026-10-01, user] truth is nevertheless RE-SWEPT in
# Phase 4b rather than reused from offset12 -- offset12's offsets came from the pre-reconversion
# upstream matrices, and re-running it is cheap once geometric_offset fits LFMM2 once per job.
#
# EMPTY AND UNDERFILLED PANELS ARE RESULTS. A candidate with 0 SNPs gets no file (and any stale
# file is removed); it is recorded with n = 0. A panel with 1-2 SNPs IS written (its GEA
# precision/recall is real) but marked usable_for_offset = FALSE: scripts/rda_offset.R FATALs
# below 3 SNPs and that failure aborts the whole replicate's Snakemake run, so such a panel must
# be left out of that replicate's sweep config. Both are listed per cohort in
# panel_underfilled_{cohort}.tsv (MVP{seed} <tab> panel <tab> n, the runbook step-6 format).
#
# SIZE CURVE (Phase 4b item 2, first variant -- user 2026-10-01): on ONE block (--size_cohort,
# default ssclines_nvar_mvar, 120 replicates) the 1/3 rule is also built at the other shares of
# the remeasure grid, RDA again x number of predictors: union_top0.1pct, union_top0.5pct,
# union_top1pct, union_top2pct (the 0.25 % point is `union` itself). The rule is held fixed and
# only the share moves, so it shows whether offset accuracy keeps improving as a GEA panel grows;
# it does NOT show whether any panel of that size would do as well (that needed the random
# floors, which were dropped).
#
# REGRESSION TIE (fails the run): per replicate, every candidate's n / causal / linked /
# background counts equal remeasure600/rdaunc_rda2x (the same p-tables, RDA x 2) -- at rung
# top_0.0025 calls_combine.tsv window 5 (union / ge2 / all3) and calls_per_method.tsv
# (LFMM / RDA / EMMAX); size-curve panels against calls_combine `union` at their own rung.
#
# Usage (inside the image; seeds default to the whole arm via benchmarks/mvp_arm.R):
#   MVP_ARM=primary_ssclines MVP_ADDED=<5 block tags> MVP_N_EXPECT=600 \
#   Rscript mvp_build_snp_sets.R --archive=PATH [--seeds=CSV] [--outdir=DIR] [--ncores=16]
#     [--params=DIR] [--check=DIR] [--sets=CSV]
#   --seeds  a subset of the arm (never `all`: that would mean the 692-row manifest)
#   --sets   rebuild only these panels; others on disk are left alone (no stale pruning)
#   --size_cohort  block that gets the size-curve panels (`none` = no size curve)
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
    library(parallel)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(PIPELINE_ROOT, "benchmarks/mvp_eval")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/pval_threshold.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/mvp_arm.R"))

args    <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt     <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
SEEDS_A <- opt("seeds", "")
OUTDIR  <- opt("outdir", file.path(EVAL, "offset13"))
PARAMS  <- opt("params", file.path(EVAL, "params_rdaunc"))
CHECK   <- opt("check",  file.path(EVAL, "remeasure600/rdaunc_rda2x"))
ARCHIVE <- opt("archive", "")
NCORES  <- as.integer(opt("ncores", "16"))
SETS_A  <- opt("sets", "")
SIZE_COHORT <- opt("size_cohort", "ssclines_nvar_mvar")
WANT    <- if (nzchar(SETS_A)) strsplit(SETS_A, ",", fixed = TRUE)[[1]] else NULL

CELL      <- "c1"
METHODS   <- c("LFMM", "RDA", "EMMAX")
TOP_SHARE <- 0.0025
RUNG      <- "top_0.0025"          # the same rung's name in remeasure600/
WINDOW_KB <- 5
RDA_OFFSET_MIN <- 3L               # scripts/rda_offset.R FATAL floor
CANDIDATES <- c("solo_lfmm", "solo_rda", "solo_emmax", "union", "best", "intersect3")
REFERENCES <- c("truth", "all", "neutral_all")
# size-curve shares (the 0.25 % point is `union`); names carry the share in percent
SIZE_SHARES <- c(0.001, 0.005, 0.01, 0.02)
SIZE_SETS   <- setNames(paste0("union_top", c("0.1", "0.5", "1", "2"), "pct"), SIZE_SHARES)
MANAGED    <- c(CANDIDATES, REFERENCES, SIZE_SETS)
# Panels of the previous build that this one no longer writes. Removed (after --archive) so
# none of them can be swept or scored as current. Anything on disk that is in NEITHER list
# stops the run: an unknown panel is never deleted on a guess.
RETIRED    <- c("solo", "rand_best1", "rand_best2", "rand_best3", "rand_solo1", "rand_union1")

if (!is.null(WANT) && length(setdiff(WANT, MANAGED)))
    stop("--sets names unknown panels: ", paste(setdiff(WANT, MANAGED), collapse = ", "))
PRUNE <- is.null(WANT)
if (PRUNE && (!nzchar(ARCHIVE) || !file.exists(ARCHIVE) || file.size(ARCHIVE) == 0))
    stop("--archive must name an existing, non-empty archive of the current snp_sets trees ",
         "before retired panels are deleted (got: '", ARCHIVE, "')")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------- seeds
MAN  <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(MAN)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))
SEEDS <- if (nzchar(SEEDS_A)) strsplit(SEEDS_A, ",", fixed = TRUE)[[1]] else PRIM$seed
if (identical(SEEDS, "all")) stop("--seeds=all is not accepted: it would mean the whole manifest")
if (length(setdiff(SEEDS, PRIM$seed)))
    stop("seeds outside the arm [", mvp_arm_label(), "]: ",
         paste(head(setdiff(SEEDS, PRIM$seed), 5), collapse = ", "))
COHORT <- setNames(PRIM$added, PRIM$seed)
if (SIZE_COHORT != "none" && !(SIZE_COHORT %in% PRIM$added))
    stop("--size_cohort ", SIZE_COHORT, " is not a block of this arm")
message(sprintf("INFO: %d replicate(s) of [%s]; params %s; outdir %s",
                length(SEEDS), mvp_arm_label(), PARAMS, OUTDIR))

# ------------------------------------------------- preflight: the method trio
# Every (seed, method) must resolve INSIDE PARAMS, LFMM/EMMAX must be the one-fit tables and
# RDA must be a real file (the uncorrected fit), before a single panel is written.
pv_path <- function(s, m) file.path(PARAMS, paste0("MVP", s), CELL, paste0(m, "_pvalues.tsv"))
pre <- rbindlist(lapply(SEEDS, function(s) rbindlist(lapply(METHODS, function(m) {
    f <- pv_path(s, m)
    data.table(seed = s, method = m, exists = file.exists(f) && file.size(f) > 0,
               link = if (nzchar(Sys.readlink(f))) Sys.readlink(f) else "")
}))))
bad <- pre[!exists |
           (method %in% c("LFMM", "EMMAX") & !grepl("/params_onefit/", link)) |
           (method == "RDA" & nzchar(link))]
if (nrow(bad)) { print(head(bad, 20)); stop("PREFLIGHT FAILED: ", nrow(bad), " (seed, method) ",
                                            "pairs are not the one-fit LFMM/EMMAX + uncorrected RDA trio") }
message("PREFLIGHT PASSED: ", nrow(pre), " p-value tables, one-fit LFMM/EMMAX + uncorrected RDA")

# ----------------------------------------------------------------- writing
# Writes only when the content differs from what is on disk, so an unchanged panel keeps its
# mtime (Snakemake would otherwise rebuild its cached models) and `unchanged` is evidence that
# the panel is byte-identical to the one an earlier sweep consumed.
write_set <- function(proj, name, keys, pmin_by_key) {
    dir <- file.path(PIPELINE_ROOT, paste0(proj, "_results"), "_intermediate", "snp_sets", name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    parts <- tstrsplit(keys, ":", fixed = TRUE)
    dt <- data.table(SNPID = keys, chr = parts[[1]], pos = as.integer(parts[[2]]),
                     min_pvalue = if (is.null(pmin_by_key)) NA_real_ else unname(pmin_by_key[keys]))
    setorder(dt, chr, pos)
    target <- file.path(dir, "selected_snps.tsv")
    tmp <- tempfile("selected_snps_", tmpdir = dir, fileext = ".tsv")
    fwrite(dt, tmp, sep = "\t")
    if (file.exists(target) && unname(tools::md5sum(tmp)) == unname(tools::md5sum(target))) {
        unlink(tmp); return("unchanged")
    }
    status <- if (file.exists(target)) "replaced" else "new"
    if (!file.rename(tmp, target)) stop("cannot move ", tmp, " -> ", target)
    status
}

remove_set <- function(proj, name) {
    dir <- file.path(PIPELINE_ROOT, paste0(proj, "_results"), "_intermediate", "snp_sets", name)
    if (!dir.exists(dir)) return(FALSE)
    if (unlink(dir, recursive = TRUE) != 0 || dir.exists(dir)) stop("cannot remove ", dir)
    TRUE
}

# --------------------------------------------------------------- one seed
build_seed <- function(s) {
    proj <- paste0("MVP", s)
    truth_f <- file.path(PIPELINE_ROOT, "data/mvp", proj, "truth_any.tsv")
    if (!file.exists(truth_f)) stop(proj, ": missing ", truth_f)
    TR <- fread(truth_f, colClasses = c(chr = "character"))
    TR[, key := paste(chr, pos, sep = ":")]
    causal <- TR[category == "causal", key]
    linked <- TR[category == "linked_neutral", key]
    bgk    <- TR[category == "background_neutral", key]

    # ---- per-method calls: top N per predictor, N from that method's own table
    L <- setNames(lapply(METHODS, function(m) load_pvalues(pv_path(s, m), "all")), METHODS)
    keysets <- lapply(L, function(l) sort(l$pv$key))
    if (!all(vapply(keysets, identical, logical(1), keysets[[1]])))
        stop(proj, ": the three methods tested different SNP sets")
    universe <- keysets[[1]]

    # RDA's one multivariate p gets share x (number of predictors the per-predictor methods test)
    n_pred <- unique(c(length(L$LFMM$trait_cols), length(L$EMMAX$trait_cols)))
    if (length(n_pred) != 1L) stop(proj, ": LFMM and EMMAX test different numbers of predictors")
    if (length(L$RDA$trait_cols) != 1L) stop(proj, ": RDA table is not one multivariate p column")
    MULT <- c(LFMM = 1, EMMAX = 1, RDA = n_pred)
    call_method <- function(m, share) {
        l  <- L[[m]]
        n  <- ceiling(share * MULT[[m]] * nrow(l$pv))
        cl <- call_by_threshold(l$pv, l$trait_cols, "top", n)
        st <- vapply(l$trait_cols, function(tc) cl$info[[tc]]$status, character(1))
        if (any(st != "ok")) stop(proj, " ", m, ": top ", n, " rule status ", paste(st, collapse = ";"))
        list(called = cl$called, n = n)
    }
    parts <- list(); call_rows <- list(); pmin <- list()
    for (m in METHODS) {
        l  <- L[[m]]
        cm <- call_method(m, TOP_SHARE)
        parts[[m]] <- cm$called
        sc <- do.call(pmin.int, c(lapply(l$trait_cols, function(x) l$pv[[x]]), list(na.rm = TRUE)))
        pmin[[m]] <- setNames(sc, l$pv$key)
        call_rows[[m]] <- data.table(seed = s, method = m, n_tests = nrow(l$pv),
                                     n_traits = length(l$trait_cols), share = TOP_SHARE * MULT[[m]],
                                     top_n = cm$n, n_called = length(cm$called))
    }
    pmin_all <- Reduce(function(a, b) { k <- union(names(a), names(b)); pmin(a[k], b[k], na.rm = TRUE) }, pmin)

    sup <- combine_support(parts, WINDOW_KB * 1000)
    sets <- list(
        solo_lfmm   = parts$LFMM,
        solo_rda    = parts$RDA,
        solo_emmax  = parts$EMMAX,
        union       = sup$snp,
        best        = sup[n_methods >= 2L, snp],
        intersect3  = sup[n_methods >= 3L, snp],
        truth       = intersect(causal, universe),
        all         = universe,
        neutral_all = intersect(bgk, universe))
    # union at 5 kb must equal the plain union of the three call sets (window-independence)
    if (!setequal(sets$union, unique(unlist(parts)))) stop(proj, ": union is not window-independent")
    # size curve: the 1/3 rule at the other shares, this block only. union needs no window.
    if (identical(COHORT[[s]], SIZE_COHORT)) {
        for (sh in names(SIZE_SETS))
            sets[[SIZE_SETS[[sh]]]] <- unique(unlist(lapply(METHODS, function(m) call_method(m, as.numeric(sh))$called)))
    }

    rows <- list(); man <- list()
    for (nm in names(sets)) {
        if (!is.null(WANT) && !(nm %in% WANT)) next
        keys <- unique(sets[[nm]])
        n <- length(keys)
        status <- if (n == 0L) {
            if (remove_set(proj, nm)) "empty_removed_stale" else "empty"
        } else write_set(proj, nm, keys, if (nm %in% c(CANDIDATES, SIZE_SETS)) pmin_all else NULL)
        rows[[nm]] <- data.table(seed = s, cohort = COHORT[[s]], set = nm, n_snps = n,
            n_causal = sum(keys %in% causal), n_linked = sum(keys %in% linked),
            n_background = sum(keys %in% bgk), status = status,
            usable_for_offset = n >= RDA_OFFSET_MIN)
        if (n == 0L) next
        man[[length(man) + 1L]] <- list(
            name = nm, n_snps = n, created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
            source_module = "benchmark_offset13",
            threshold_type = if (nm %in% c(CANDIDATES, SIZE_SETS))
                sprintf("top_share_%s_per_predictor_rda_x%d",
                        if (nm %in% SIZE_SETS) names(SIZE_SETS)[SIZE_SETS == nm] else TOP_SHARE, n_pred)
                else "none",
            threshold_value = if (nm %in% SIZE_SETS) as.numeric(names(SIZE_SETS)[SIZE_SETS == nm])
                              else if (nm %in% CANDIDATES) TOP_SHARE else 0,
            regime = "snp",
            strategy = if (nm %in% SIZE_SETS) "union" else
                       switch(nm, union = "union", best = "at_least_2_5kb", intersect3 = "all_3_5kb",
                              truth = "causal_loci", all = "all_tested", neutral_all = "background_neutral",
                              "single_method"),
            traits = L$LFMM$trait_cols,
            methods = if (nm %in% SIZE_SETS) METHODS else
                      switch(nm, union = , best = , intersect3 = METHODS,
                             solo_lfmm = "LFMM", solo_rda = "RDA", solo_emmax = "EMMAX", character(0)))
    }

    # ---- retired panels: deleted, never left to be swept as current
    removed <- character(0)
    if (PRUNE) {
        store <- file.path(PIPELINE_ROOT, paste0(proj, "_results"), "_intermediate", "snp_sets")
        on_disk <- basename(list.dirs(store, recursive = FALSE))
        unknown <- setdiff(on_disk, c(MANAGED, RETIRED))
        if (length(unknown)) stop(proj, ": unknown panel dir(s) on disk, not deleting on a guess: ",
                                  paste(unknown, collapse = ", "))
        for (nm in intersect(on_disk, RETIRED)) if (remove_set(proj, nm)) removed <- c(removed, nm)
    }

    man_f <- file.path(PIPELINE_ROOT, paste0(proj, "_results"), "_intermediate", "snp_sets", "manifest.json")
    if (!is.null(WANT) && file.exists(man_f)) {
        old  <- tryCatch(fromJSON(man_f, simplifyDataFrame = FALSE), error = function(e) list())
        keep <- Filter(function(e) !(e$name %in% WANT), old)
        man  <- c(keep, man)
    }
    write_json(man, man_f, auto_unbox = TRUE, pretty = TRUE)

    list(rows = rbindlist(rows), calls = rbindlist(call_rows),
         removed = data.table(seed = rep(s, length(removed)), panel = removed))
}

t0  <- Sys.time()
res <- mclapply(SEEDS, function(s) tryCatch(build_seed(s),
                error = function(e) list(error = data.table(seed = s, error = conditionMessage(e)))),
                mc.cores = NCORES, mc.preschedule = FALSE)
message(sprintf("INFO: built in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
failed <- vapply(res, function(r) !is.list(r) || !is.null(r$error), logical(1))
if (any(failed)) {
    errs <- rbindlist(lapply(res[failed], function(r)
        if (is.list(r)) r$error else data.table(seed = NA_character_, error = as.character(r))))
    print(errs); fwrite(errs, file.path(OUTDIR, "errors.tsv"), sep = "\t")
    stop(nrow(errs), " replicate(s) failed (see errors.tsv)")
}
SUM   <- rbindlist(lapply(res, `[[`, "rows"))
CALLS <- rbindlist(lapply(res, `[[`, "calls"))
REM   <- rbindlist(lapply(res, `[[`, "removed"))
setorder(SUM, seed, set); setorder(CALLS, seed, method)

# --------------------------------------------------------------- contracts
SUM[, untracked := n_snps - n_causal - n_linked - n_background]
w <- dcast(SUM, seed ~ set, value.var = "n_snps")
chk <- function(ok, msg) if (!isTRUE(all(ok))) stop("CONTRACT FAILED: ", msg)
chk(SUM[set %in% CANDIDATES, untracked == 0L], "a candidate holds SNPs outside causal/linked/background")
if (all(c("union", "best", "intersect3") %in% names(w))) {
    chk(w$union >= w$best & w$best >= w$intersect3, "union >= 2/3 >= 3/3 violated")
    for (m in c("solo_lfmm", "solo_rda", "solo_emmax")) if (m %in% names(w))
        chk(w$union >= w[[m]], paste("union smaller than", m))
}
if (all(c("all", "union") %in% names(w))) chk(w$all >= w$union, "`all` smaller than `union`")
if (all(SIZE_SETS %in% names(w))) {
    ws <- w[!is.na(get(SIZE_SETS[[1]]))]
    chk(ws[[SIZE_SETS[[1]]]] <= ws$union & ws$union <= ws[[SIZE_SETS[[2]]]] &
        ws[[SIZE_SETS[[2]]]] <= ws[[SIZE_SETS[[3]]]] & ws[[SIZE_SETS[[3]]]] <= ws[[SIZE_SETS[[4]]]],
        "size-curve panels do not grow with the share")
}
chk(SUM[set %in% SIZE_SETS, untracked == 0L], "a size-curve panel holds SNPs outside causal/linked/background")
if ("truth" %in% SUM$set) {
    tr <- SUM[set == "truth"]
    chk(tr$n_snps == tr$n_causal, "truth holds non-causal SNPs")
    notsame <- tr[status != "unchanged"]
    if (nrow(notsame)) {
        print(notsame[, .(seed, n_snps, status)])
        stop("TRUTH PANEL CHANGED on ", nrow(notsame), " replicate(s): it is a pure function of ",
             "truth_any.tsv, which Phase 2b verified byte-identical, so a change means an input moved")
    }
}

# ---- regression tie against remeasure600/rdaunc
CMB <- fread(file.path(CHECK, "calls_combine.tsv"), colClasses = c(seed = "character"))[
    rung == RUNG & window_kb == WINDOW_KB]
PMT <- fread(file.path(CHECK, "calls_per_method.tsv"), colClasses = c(seed = "character"))[rung == RUNG]
ref <- rbind(
    CMB[, .(seed, set = c(union = "union", ge2 = "best", all3 = "intersect3")[combine],
            ref_n = n_called, ref_causal = tp_any_causal, ref_linked = expected_linked,
            ref_background = fp_background)],
    PMT[, .(seed, set = paste0("solo_", tolower(method)),
            ref_n = n_called, ref_causal = tp_any_causal, ref_linked = expected_linked,
            ref_background = fp_background)])
CMB_ALL <- fread(file.path(CHECK, "calls_combine.tsv"), colClasses = c(seed = "character"))
size_ref <- rbindlist(lapply(names(SIZE_SETS), function(sh) {
    rk <- paste0("top_", format(as.numeric(sh), scientific = FALSE, drop0trailing = TRUE))
    CMB_ALL[rung == rk & window_kb == WINDOW_KB & combine == "union",
            .(seed, set = SIZE_SETS[[sh]], ref_n = n_called, ref_causal = tp_any_causal,
              ref_linked = expected_linked, ref_background = fp_background)]
}))
ref <- rbind(ref, size_ref)
tie <- merge(SUM[set %in% c(CANDIDATES, SIZE_SETS), .(seed, set, n_snps, n_causal, n_linked, n_background)],
             ref, by = c("seed", "set"), all.x = TRUE)
tie[, match := !is.na(ref_n) & n_snps == ref_n & n_causal == ref_causal &
               n_linked == ref_linked & n_background == ref_background]
fwrite(tie, file.path(OUTDIR, "regression_tie.tsv"), sep = "\t")
if (!all(tie$match)) {
    print(head(tie[match == FALSE], 20))
    stop("REGRESSION TIE FAILED: ", sum(!tie$match), " of ", nrow(tie),
         " (seed, panel) rows differ from ", CHECK, " at ", RUNG, " / ", WINDOW_KB, " kb")
}
message(sprintf("REGRESSION TIE PASSED: %d (seed, panel) rows equal %s at %s / %g kb (+ size curve at its own rungs)",
                nrow(tie), basename(CHECK), RUNG, WINDOW_KB))

# ------------------------------------------------------------------ write
SUM_F <- file.path(OUTDIR, "snp_sets_summary.tsv")
if (!is.null(WANT) && file.exists(SUM_F)) {
    OLD <- fread(SUM_F, colClasses = c(seed = "character"))
    SUM <- rbind(OLD[!paste(seed, set) %in% SUM[, paste(seed, set)]], SUM, fill = TRUE)
    setorder(SUM, seed, set)
}
fwrite(SUM, SUM_F, sep = "\t")
fwrite(CALLS, file.path(OUTDIR, "method_calls.tsv"), sep = "\t")
fwrite(REM, file.path(OUTDIR, "retired_panels_removed.tsv"), sep = "\t")

# Pre-screen, runbook step-6 format, one file per cohort plus the pooled one. Candidates and
# truth only: all/neutral_all are never swept.
UF <- SUM[set %in% c(CANDIDATES, SIZE_SETS, "truth") & n_snps < RDA_OFFSET_MIN,
          .(project = paste0("MVP", seed), set, n_snps, cohort)]
setorder(UF, cohort, project, set)
fwrite(UF[, .(project, set, n_snps)], file.path(OUTDIR, "panel_underfilled.tsv"),
       sep = "\t", col.names = FALSE)
for (co in unique(SUM$cohort))
    fwrite(UF[cohort == co, .(project, set, n_snps)],
           file.path(OUTDIR, sprintf("panel_underfilled_%s.tsv", co)), sep = "\t", col.names = FALSE)

writeLines(c(sprintf("params_dir\t%s", PARAMS), sprintf("cell\t%s", CELL),
             sprintf("top_share\t%s", TOP_SHARE), sprintf("window_kb\t%s", WINDOW_KB),
             sprintf("rda_share_mult\tnumber of predictors (2)"),
             sprintf("size_cohort\t%s", SIZE_COHORT),
             sprintf("size_shares\t%s", paste(names(SIZE_SETS), collapse = ",")),
             sprintf("check_dir\t%s", CHECK), sprintf("archive\t%s", ARCHIVE),
             sprintf("corpus\t%s", mvp_arm_label()), sprintf("n_replicates\t%d", length(SEEDS)),
             sprintf("sets\t%s", if (is.null(WANT)) paste(MANAGED, collapse = ",") else SETS_A),
             sprintf("run_at\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))),
           file.path(OUTDIR, "provenance.tsv"))

message("\nPanel sizes (median [min-max]) and status counts:")
print(SUM[, .(median_n = median(n_snps), min_n = min(n_snps), max_n = max(n_snps),
              empty = sum(n_snps == 0L), lt3 = sum(n_snps < RDA_OFFSET_MIN)), by = set])
print(dcast(SUM, set ~ status, fun.aggregate = length))
message(sprintf("retired panel dirs removed: %d; underfilled (n<3) candidate/truth panels: %d on %d replicate(s)",
                nrow(REM), nrow(UF), uniqueN(UF$project)))
message("OK -> ", OUTDIR)
