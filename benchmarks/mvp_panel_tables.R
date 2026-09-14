#!/usr/bin/env Rscript
# =============================================================================
# mvp_panel_tables.R -- the two tier-A input tables every figure script reads.
#
# WHY THIS EXISTS. `phase1_seed_medians_solo.tsv` and `panel_pr_recomputed.tsv` are read by
# mvp_manuscript_figures.R, mvp_absolute_figures.R, mvp_panel_distributions.R and
# mvp_panel_explainer.R -- but nothing in the repo WROTE them. They were produced ad hoc in
# the session that built the 32-replicate figure set, which means the figures could not be
# regenerated from a clean checkout and the definitions lived only in a transcript.
#
# Both definitions below were recovered by re-deriving them against the existing
# offset09 files until they reproduced exactly (see --check), not by guessing:
#
#   phase1_seed_medians_solo.tsv   median tau per (seed, marker_set, method_label), taken
#                                  over the 100 LANDSCAPE gardens only -- the 12 climate-
#                                  novelty gardens are excluded -- and over PRIMARY seeds
#                                  only, dropping the 2 degenerate Est-Clines controls.
#                                  Verified: all-112-garden medians do NOT reproduce the
#                                  published file (-0.6736 vs -0.6804 on the first cell);
#                                  landscape-only does, to 15 significant digits.
#
#   panel_pr_recomputed.tsv        per (seed, panel) detection precision/recall of the
#                                  marker panel itself against that replicate's truth table,
#                                  keyed on chr:pos. Recomputed from the panel files rather
#                                  than read from snp_sets_summary.tsv, because that file is
#                                  overwritten by whichever --sets run wrote it last and so
#                                  covers a different replicate count per panel.
#                                  Controls are KEPT here (it is a property of the panel, not
#                                  of the offset), matching the published file.
#
# `_solo` in the first name is historical: it marks the version that includes the three
# single-method panels (solo_lfmm / solo_rda / solo_emmax). It is the file the figures read.
#
# Usage:
#   Rscript mvp_panel_tables.R [--outdir=DIR] [--seeds=all|CSV] [--check=DIR]
#
#   --check  regenerate into --outdir and diff against the tables already in DIR, cell by
#            cell. This is the regression gate: run it against offset09 before trusting a
#            new cohort's numbers.
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))   # parse_kv_args()

args   <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt    <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
OUTDIR <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset10"))
SEEDS_A <- opt("seeds", "all")
CHECK  <- opt("check", "")

MAN <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
             colClasses = c("seed" = "character"))
SEEDS <- if (SEEDS_A == "all") MAN$seed else strsplit(SEEDS_A, ",", fixed = TRUE)[[1]]

# ---------------------------------------------------------------- medians ------
GP_F <- file.path(OUTDIR, "garden_performance.tsv")
if (!file.exists(GP_F)) stop("Missing scorer output: ", GP_F,
                             " -- run eval_offset_lind.R first")
GP <- fread(GP_F, colClasses = c("seed" = "character"))

med <- GP[garden_type == "landscape" & control == FALSE & seed %in% SEEDS & !is.na(tau),
          .(tau = median(tau)), by = .(seed, marker_set, method_label)]
setorder(med, seed, marker_set, method_label)
fwrite(med, file.path(OUTDIR, "phase1_seed_medians_solo.tsv"), sep = "\t")
message(sprintf("wrote phase1_seed_medians_solo.tsv -- %d rows, %d seeds, %d panels, %d methods",
                nrow(med), uniqueN(med$seed), uniqueN(med$marker_set),
                uniqueN(med$method_label)))

# ---------------------------------------------------------------- panel PR -----
# Which panel directories to score. Explicit rather than a glob, and two exclusions matter:
#
#   rand_best2 / rand_best3  extra draws at the `best` size, there to show the floor's own
#                            sampling noise. Including them would put three rows of the same
#                            panel into every per-panel median.
#   solo                     BYTE-IDENTICAL to solo_lfmm -- both are LFMM alone at the same
#                            operating point (mvp_build_snp_sets.R's SOLO <- "LFMM"). It is
#                            kept on disk because the offset arm harvests it under its own
#                            name, but scoring both would double-count one panel.
#
# Verified against the published 32-replicate table: 12 panels x 32 seeds = 384 rows.
PANELS <- c("all", "neutral_all", "truth", "union", "best", "intersect3",
            "solo_lfmm", "solo_rda", "solo_emmax", "rand_best1",
            "rand_solo1", "rand_union1")

# AN EMPTY PANEL IS A RESULT, NOT MISSING DATA. `mvp_build_snp_sets.R` does not write a
# set file when the combination rule yields nothing ("intersect3 EMPTY -- set not written"),
# so a panel directory can be absent because the three methods agreed on ZERO loci on that
# replicate -- which is exactly what the GEA comparison should show. Emitting n = 0 keeps
# that visible; skipping the row would silently turn "the method found nothing" into "we
# did not measure this", and would bias every per-panel median upward by dropping precisely
# the replicates where the panel failed.
#
# RDA_OFFSET_MIN is the pipeline's own floor for the offset arm: scripts/rda_offset.R:174
# stops with FATAL when fewer than 3 candidate SNPs survive zero-variance filtering, and
# because that rule failing kills the whole seed's Snakemake run, such panels have to be
# dropped from the maladaptation config. `usable_for_offset` records that split so the two
# modules can be read side by side: the GEA table keeps every panel, the offset table
# carries only the usable ones, and the exclusions are enumerated rather than inferred.
RDA_OFFSET_MIN <- 3L

pr_rows <- list()
for (s in SEEDS) {
    proj    <- paste0("MVP", s)
    truth_f <- file.path(PIPELINE_ROOT, "data/mvp", proj, "truth_any.tsv")
    if (!file.exists(truth_f)) { message("  ", proj, ": no truth table -- skipped"); next }
    TR <- fread(truth_f, colClasses = c("chr" = "character"))
    TR[, key := paste(chr, pos, sep = ":")]
    causal <- TR[category == "causal"]$key
    linked <- TR[category == "linked_neutral"]$key

    for (p in PANELS) {
        f <- file.path(PIPELINE_ROOT, paste0(proj, "_results"), "_intermediate",
                       "snp_sets", p, "selected_snps.tsv")
        if (file.exists(f)) {
            S <- fread(f, colClasses = c("chr" = "character"))
            k <- paste(S$chr, S$pos, sep = ":")
        } else {
            k <- character(0)     # combination rule produced an empty set
        }
        n_causal <- sum(k %in% causal)
        n_linked <- sum(k %in% linked)
        pr_rows[[length(pr_rows) + 1L]] <- data.table(
            seed = s, set = p, n = length(k),
            n_causal = n_causal, n_linked = n_linked,
            causal_total = length(causal),
            # Precision counts ONLY causal loci as hits. Linked-neutral markers are reported
            # separately (n_linked) rather than folded in either direction -- they are
            # physically linked to true QTNs, so they are neither a clean hit nor a clean
            # false positive, and the composition figures need them as their own column.
            # An empty panel has no precision (0/0), but its recall IS 0 -- it found none
            # of the causal loci, which is a measurement, not an undefined quantity.
            precision = if (length(k)) n_causal / length(k) else NA_real_,
            recall    = if (length(causal)) n_causal / length(causal) else NA_real_,
            usable_for_offset = length(k) >= RDA_OFFSET_MIN)
    }
}
PR <- rbindlist(pr_rows)
setorder(PR, seed, set)
fwrite(PR, file.path(OUTDIR, "panel_pr_recomputed.tsv"), sep = "\t")
message(sprintf("wrote panel_pr_recomputed.tsv -- %d rows, %d seeds, %d panels",
                nrow(PR), uniqueN(PR$seed), uniqueN(PR$set)))

# ---------------------------------------------------- offset exclusions --------
# One row per (seed, panel) that the GEA side measured but the offset side could not use.
# This is the manuscript's own note: those replicates are absent from the offset medians
# because the panel was empty or near-empty, NOT because the offset failed to converge.
EXC <- PR[usable_for_offset == FALSE,
          .(seed, set, n_snps = n, n_causal,
            reason = fifelse(n == 0L, "panel empty -- combination rule selected no loci",
                     sprintf("panel too small for rda_offset (needs >= %d SNPs, rda_offset.R:174)",
                             RDA_OFFSET_MIN)))]
setorder(EXC, set, seed)
fwrite(EXC, file.path(OUTDIR, "panel_offset_exclusions.tsv"), sep = "\t")
if (nrow(EXC)) {
    message(sprintf("wrote panel_offset_exclusions.tsv -- %d (seed, panel) excluded from the offset arm:",
                    nrow(EXC)))
    print(EXC[, .N, by = set])
} else {
    message("wrote panel_offset_exclusions.tsv -- none; every panel is usable for the offset arm")
}

# Per-panel denominators, so no figure has to infer n from row counts.
COV <- PR[, .(n_seeds_gea = .N,
              n_empty = sum(n == 0L),
              n_seeds_offset = sum(usable_for_offset),
              median_n_snps = as.integer(median(n))), by = set]
setorder(COV, set)
fwrite(COV, file.path(OUTDIR, "panel_coverage.tsv"), sep = "\t")
message("wrote panel_coverage.tsv:")
print(COV)

# ---------------------------------------------------------------- check --------
if (nzchar(CHECK)) {
    message("\n== regression check against ", CHECK, " ==")
    # Gate outcome, accumulated across cmp_one() calls. THE TRAP THIS CLOSES: the
    # gate used to `return(invisible())` on every non-comparison and never set an
    # exit status, so BOTH "0 rows were compared" and "the numbers moved" printed
    # a line and exited 0. On ssclines_b1 it printed
    #   phase1_seed_medians_solo.tsv: ROW COUNT 3580 (ref) vs 0 (new)
    # and carried on -- zero rows compared, reported in a shape that reads like a
    # pass. Exit is now 0 = compared and identical, 1 = FAIL (a real regression),
    # 2 = NOT APPLICABLE (nothing was comparable, e.g. disjoint seed sets, which
    # is the expected state for every fresh SS-Clines block).
    N_FAIL <- 0L; N_NA <- 0L; N_OK <- 0L
    cmp_one <- function(fn, keys) {
        old_f <- file.path(CHECK, fn)
        if (!file.exists(old_f)) {
            message("  ", fn, ": NOT APPLICABLE -- absent in reference")
            N_NA <<- N_NA + 1L; return(invisible())
        }
        old <- fread(old_f, colClasses = c("seed" = "character"))
        new <- fread(file.path(OUTDIR, fn), colClasses = c("seed" = "character"))
        n_seed_shared <- length(intersect(unique(old$seed), unique(new$seed)))
        new <- new[seed %in% old$seed]                       # reference may cover fewer seeds
        setkeyv(old, keys); setkeyv(new, keys)
        if (nrow(new) == 0L) {
            message(sprintf(paste0("  %s: NOT APPLICABLE -- 0 of %d reference seeds present ",
                                   "in the new output (disjoint seed sets). NOTHING WAS COMPARED."),
                            fn, length(unique(old$seed))))
            N_NA <<- N_NA + 1L; return(invisible())
        }
        if (nrow(old) != nrow(new)) {
            message(sprintf(paste0("  %s: FAIL -- ROW COUNT %d (ref) vs %d (new, restricted to ",
                                   "ref seeds); %d seeds in common"),
                            fn, nrow(old), nrow(new), n_seed_shared))
            N_FAIL <<- N_FAIL + 1L; return(invisible())
        }
        j <- merge(old, new, by = keys, suffixes = c(".ref", ".new"))
        vals <- setdiff(names(old), keys)
        worst <- 0
        for (v in vals) {
            a <- j[[paste0(v, ".ref")]]; b <- j[[paste0(v, ".new")]]
            d <- if (is.numeric(a)) max(abs(a - b), na.rm = TRUE) else as.numeric(!all(a == b))
            worst <- max(worst, d)
            if (d > 1e-9) message(sprintf("    %s: max |diff| = %.3g", v, d))
        }
        message(sprintf("  %s: %d rows compared, max |diff| = %.3g %s",
                        fn, nrow(j), worst, if (worst <= 1e-9) "-- IDENTICAL" else "-- FAIL, DIFFERS"))
        if (worst <= 1e-9) N_OK <<- N_OK + 1L else N_FAIL <<- N_FAIL + 1L
    }
    cmp_one("phase1_seed_medians_solo.tsv", c("seed", "marker_set", "method_label"))
    cmp_one("panel_pr_recomputed.tsv",      c("seed", "set"))

    # One unambiguous verdict line, then an exit status that distinguishes the
    # three outcomes. A caller under `set -e` now stops on a real regression and
    # can special-case 2 ("this block shares no seeds with the reference").
    message(sprintf("\n== gate: %d identical, %d FAIL, %d not applicable ==",
                    N_OK, N_FAIL, N_NA))
    if (N_FAIL > 0L) {
        message("GATE FAILED -- the legacy numbers moved. Stop and diagnose.")
        quit(status = 1L)
    }
    # ANY uncompared table exits 2, not just all of them. The two tables have
    # different seed coverage -- panel_pr_recomputed is rebuilt for every seed it
    # can find while phase1_seed_medians_solo follows --outdir's
    # garden_performance -- so a block run against the legacy reference compares
    # one and skips the other. Reporting that as PASSED is the same trap in a
    # smaller box: a table nobody checked must never be inside a green verdict.
    if (N_NA > 0L) {
        message(sprintf(paste0("GATE NOT APPLICABLE -- %d of %d table(s) were never compared. ",
                               "This is NOT a pass: that table got no regression check at all."),
                        N_NA, N_NA + N_OK + N_FAIL))
        quit(status = 2L)
    }
    message("GATE PASSED.")
}
