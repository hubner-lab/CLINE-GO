#!/usr/bin/env Rscript
# =============================================================================
# mvp_build_snp_sets.R -- build the SNP sets that feed mode=maladaptation, per seed.
#
# The offset arm asks ONE question: which SNPs should be handed to the offset model?
# This script writes the candidate answers, plus the two references that put them on an
# absolute scale, into {PROJECT}_results/_intermediate/snp_sets/{name}/selected_snps.tsv --
# the store resolve_active_snp_sets() (workflow/rules/common.smk:1207) globs at parse time.
#
# THREE CANDIDATES, ONE CEILING, ONE FLOOR (size-matched)
#
#   best        LFMM+RDA+EMMAX, at_least_2, 5 kb agreement window   <- candidate
#   union       LFMM+RDA+EMMAX, union (at_least_1), 0 kb            <- candidate
#   solo        LFMM alone                                          <- candidate
#   truth       the seed's causal loci                              <- CEILING (not a candidate)
#   rand_*      draws from background-neutral loci only             <- FLOOR (not a candidate)
#
# Why these three candidates. `best` is scheme S1 of benchmarks/mvp_anchor_schemes.tsv,
# which journal 07's cross-seed vote ranks first on the base arm for highly-polygenic seeds
# and statistically indistinguishable from the union schemes for mod-polygenic ones
# (med_f1 0.3157 vs 0.328/0.3304, benchmarks/mvp_eval/report07/recommended_schemes.tsv).
# It is applied UNIFORMLY to all 14 seeds rather than per architecture: a per-seed scheme
# choice would make the cross-seed comparison a comparison of schemes as well as of set
# sizes. `union` is the same three methods with the agreement requirement dropped -- the
# maximise-recall arm. `solo` is scheme B1, the single best method, which separates "does
# combining help" from "does selecting help".
#
# CALIBRATION IS TRANSFERRED, NOT RE-DERIVED FROM TRUTH. All three candidates use S1's
# per-method operating points VERBATIM -- LFMM top_100, RDA top_100, EMMAX custom 1e-04
# (benchmarks/mvp_anchor_schemes.tsv, fixed on the anchor seed 1232548 in journal 06). This
# is journal 07's `verbatim` transfer mode: applying these points to a new replicate needs
# no truth, so the candidate sets stay reproducible on real data. Journal 07's `recal` mode
# re-picks the points using each replicate's own truth and is deliberately NOT used --
# it would leak the answer into the candidates.
#
# Why not the pipeline's own bonf 0.05, which would have been the more natural default:
# measured on all 14 seeds at their default rungs, it calls 0-34 SNPs per method and RDA
# calls 0 on 8 of 14, so `best` (agreement of >= 2 methods) comes out EMPTY on most seeds --
# there is nothing to hand a Gradient Forest. qval 0.1 is worse: it calls exactly 0 on every
# seed and method (the qvalue pi0 fit degenerates on these p-value distributions). custom
# 1e-04 calls 0-62. Only the top_N points produce sets of a size an offset model can use.
# Consequence to keep in mind when reading the results: set size is roughly constant across
# seeds by construction, while the number of causal loci is not (6 to 395), so `best` is a
# far larger share of the truth on the oligogenic seeds than on the polygenic ones.
#
# The three candidates therefore differ ONLY in the combination rule, which is exactly the
# question asked. Truth enters this script in two places only: the `truth` ceiling set, and
# the background-neutral pool the floor is drawn from.
#
# OPERATING POINT PER METHOD is read from each method's own is_default rung in
# benchmarks/mvp_sweep_cells*.tsv (LFMM K = k_best, EMMAX/RDA n_pcs = k_best). The ladders
# desync on seeds whose k_best is small enough for the ladder to hit its floor -- LFMM's
# default lands in c2 there while EMMAX/RDA stay in c3 -- so this is looked up per (seed,
# method) and never assumed to be one cell.
#
# THE FLOOR IS DRAWN FROM background_neutral ONLY. On seed 1232548 the truth table is
# 70 causal / 5207 linked_neutral / 5048 background_neutral, so sample(all_snps) would be
# ~50% linked loci -- physically linked to true QTNs, and therefore not a null at all.
# Three independent draws at the `best` size (a single draw's sampling noise would read as
# a method effect), plus one draw each at the `union` and `solo` sizes so every candidate
# has a size-matched floor and set size cannot be confounded with set quality.
#
# Usage:
#   Rscript mvp_build_snp_sets.R [--seeds=all|CSV] [--outdir=DIR] [--adjust=bonf] [--value=0.05]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/pval_threshold.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args    <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt     <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
SEEDS_A <- opt("seeds", "all")
OUTDIR  <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset08"))
# Which sets to (re)write. Default: every set. Pass e.g. --sets=all to add the `all` panel
# without touching the others -- rewriting a set file changes its mtime, which would make
# Snakemake rebuild that set's cached Gradient Forest model for no reason.
SETS_A  <- opt("sets", "")
WANT    <- if (nzchar(SETS_A)) strsplit(SETS_A, ",", fixed = TRUE)[[1]] else NULL
PARAMS  <- opt("params", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/params07"))
PARAMS_ANCHOR <- opt("params_anchor", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/params"))

# S1's per-method operating points, verbatim from benchmarks/mvp_anchor_schemes.tsv.
POINTS <- list(LFMM  = c("top",    "100"),
               RDA   = c("top",    "100"),
               EMMAX = c("custom", "1e-04"))
SUBSET    <- c("LFMM", "RDA", "EMMAX")   # scheme S1 / union arm
SOLO      <- "LFMM"                      # single-method arm, same operating point
WINDOW_KB <- 5                           # S1's agreement window
N_RAND_BEST <- 3

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MAN   <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
               colClasses = c("seed" = "character"))
SEEDS <- if (SEEDS_A == "all") MAN$seed else strsplit(SEEDS_A, ",", fixed = TRUE)[[1]]

# ---- default rung per (seed, method), unioned over the two cell manifests
# Group 1's 18 new seeds are written to their own manifest so a seed-scoped
# regeneration cannot truncate the one describing the runs that produced
# benchmarks/mvp_eval/sweep/. It must be listed here too, or default_cell()
# silently resolves nothing for every new seed and the panels come out empty --
# after the upstream hours have already been spent.
#
# The list is EXPLICIT, not a glob over mvp_sweep_cells_*.tsv: the _m05 / _p11ibs / _repair
# manifests describe different arms and carry their own is_default rows for the same
# (seed, method), which would make default_cell()'s uniqueness check fail. Each new cohort
# adds one line here (or is passed via MVP_CELL_FILES for a one-off).
cell_files <- c(file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_j07.tsv"),
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_p11.tsv"),
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_group1.tsv"),
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_group2.tsv"),
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_group3.tsv"),
                # SS-Clines complete-block arm: one manifest per demography block. Listed
                # up front rather than added block by block -- the vector is filtered by
                # file.exists() below, so naming a manifest before its block has run is a
                # no-op, while forgetting to add one after the sweep is a silent empty-panel
                # failure. Still EXPLICIT, not a glob, for the reason given above.
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_ssclines_nvar_mvar.tsv"),
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_ssclines_ncline_ns.tsv"),
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_ssclines_nequal_mconst.tsv"),
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_ssclines_ncline_ctredge.tsv"),
                file.path(PIPELINE_ROOT, "benchmarks/mvp_sweep_cells_ssclines_nequal_mbreaks.tsv"))
extra_cells <- Sys.getenv("MVP_CELL_FILES", "")
if (nzchar(extra_cells)) {
    cell_files <- c(cell_files, strsplit(extra_cells, ",", fixed = TRUE)[[1]])
}
message("cell manifests: ",
        paste(basename(cell_files[file.exists(cell_files)]), collapse = ", "))
CELLS <- rbindlist(lapply(cell_files[file.exists(cell_files)],
                          fread, colClasses = c("seed" = "character")), fill = TRUE)
DEFCELL <- unique(CELLS[is_default == TRUE, .(seed, method, cell, param, value)])

default_cell <- function(s, m) {
    r <- DEFCELL[seed == s & method == m]
    if (nrow(r) != 1L) stop("No unique default rung for seed ", s, " method ", m)
    r$cell[1]
}

# The anchor seed 1232548 was harvested by journal 06 into mvp_eval/params/, the other 13 by
# journal 07 into mvp_eval/params07/. Resolve per seed rather than assuming one root.
pv_path <- function(s, m) {
    cell <- default_cell(s, m)
    for (root in c(PARAMS, PARAMS_ANCHOR)) {
        p <- file.path(root, paste0("MVP", s), cell, paste0(m, "_pvalues.tsv"))
        if (file.exists(p)) return(p)
    }
    NA_character_
}

write_set <- function(proj, name, keys, pmin_by_key) {
    dir <- file.path(PIPELINE_ROOT, paste0(proj, "_results"), "_intermediate", "snp_sets", name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    parts <- tstrsplit(keys, ":", fixed = TRUE)
    dt <- data.table(SNPID = keys, chr = parts[[1]], pos = as.integer(parts[[2]]),
                     min_pvalue = if (is.null(pmin_by_key)) NA_real_ else pmin_by_key[keys])
    setorder(dt, chr, pos)
    fwrite(dt, file.path(dir, "selected_snps.tsv"), sep = "\t")
    nrow(dt)
}

summary_rows <- list()
manifest_rows <- list()

for (s in SEEDS) {
    proj <- paste0("MVP", s)
    message("\n=== ", proj, " ===")

    truth_f <- file.path(PIPELINE_ROOT, "data/mvp", proj, "truth_any.tsv")
    if (!file.exists(truth_f)) stop(proj, ": missing ", truth_f)
    TR <- fread(truth_f, colClasses = c("chr" = "character"))
    TR[, key := paste(chr, pos, sep = ":")]

    # ---- per-method calls at the pipeline's own threshold
    parts <- list(); pmin <- list(); universe <- NULL
    for (m in unique(c(SUBSET, SOLO))) {
        p <- pv_path(s, m)
        if (is.na(p)) { message("  ", m, ": NO p-value table -- method dropped for this seed"); next }
        L  <- load_pvalues(p, "all")
        pt <- POINTS[[m]]
        cl <- call_by_threshold(L$pv, L$trait_cols, pt[1], pt[2])
        parts[[m]] <- cl$called
        sc <- do.call(pmin.int, c(lapply(L$trait_cols, function(x) L$pv[[x]]), list(na.rm = TRUE)))
        pmin[[m]] <- setNames(sc, L$pv$key)
        if (is.null(universe)) universe <- L$pv$key
        message(sprintf("  %-6s cell=%s  %s_%s  called=%d",
                        m, default_cell(s, m), pt[1], pt[2], length(cl$called)))
    }
    if (is.null(universe)) stop(proj, ": no p-value tables at all")

    # min p across methods, for the min_pvalue column of the candidate sets
    pmin_all <- Reduce(function(a, b) {
        k <- union(names(a), names(b))
        pmin(a[k], b[k], na.rm = TRUE)
    }, pmin)

    sub_parts <- parts[intersect(SUBSET, names(parts))]

    sets <- list()
    # best: agreement of >= 2 of the three, within a 5 kb window
    cs <- combine_support(sub_parts, WINDOW_KB * 1000)
    sets[["best"]]  <- if (nrow(cs)) cs[n_methods >= 2, snp] else character(0)
    # union: any of the three, exact position
    cs0 <- combine_support(sub_parts, 0)
    sets[["union"]] <- if (nrow(cs0)) cs0[n_methods >= 1, snp] else character(0)
    # solo: LFMM alone (kept under its historical name so journals 08/09 still resolve)
    sets[["solo"]]  <- if (!is.null(parts[[SOLO]])) parts[[SOLO]] else character(0)
    # One panel per GEA method at its OWN verbatim operating point. These are what
    # make "combining methods helps" a measured claim rather than an assumption:
    # each is the same scan, same calibration, just without the agreement step, so
    # the only thing separating them from `best`/`intersect3` is the combination rule.
    for (m in SUBSET) {
        sets[[paste0("solo_", tolower(m))]] <-
            if (!is.null(parts[[m]])) parts[[m]] else character(0)
    }
    # ceiling: causal loci that survived filtering (i.e. are in the tested universe)
    sets[["truth"]] <- intersect(TR[category == "causal", key], universe)
    # `all` -- every SNP the GEA actually tested. This is the "all" marker panel of
    # Lind & Lotterhos 2025 and the practitioner default; it is a candidate, not a reference.
    sets[["all"]]   <- universe
    # `neutral_all` -- THEIR neutral panel, verbatim: "SNPs on linkage groups without any QTNs"
    # (mean N = 16,520 in their sims). Our background_neutral category is exactly that set --
    # convert_mvp.R:464 assigns it to LG 11-20, the QTN-free linkage groups, quoting the source
    # paper's own Methods. The size-matched `rand_best1` draw stays as OUR floor; this panel is
    # what makes the adaptive-vs-neutral comparison commensurable with their published numbers.
    sets[["neutral_all"]] <- intersect(TR[category == "background_neutral", key], universe)
    # `intersect3` -- maximum-precision GEA panel: agreement of ALL three methods. Complements
    # `union` (maximum recall) and `best` (>= 2 of 3). Written only if non-empty.
    if (length(sub_parts) == 3L) {
        cs3 <- combine_support(sub_parts, WINDOW_KB * 1000)
        sets[["intersect3"]] <- if (nrow(cs3)) cs3[n_methods >= 3, snp] else character(0)
    }

    # floor: background-neutral only, size-matched to each candidate
    bg_pool <- intersect(TR[category == "background_neutral", key], universe)
    draw <- function(n, tag, i) {
        if (n == 0 || length(bg_pool) == 0) return(character(0))
        n <- min(n, length(bg_pool))
        set.seed(as.integer(substr(s, 4, 7)) * 100L + i)   # reproducible, distinct per draw
        sample(bg_pool, n)
    }
    for (i in seq_len(N_RAND_BEST))
        sets[[paste0("rand_best", i)]] <- draw(length(sets[["best"]]),  "best",  i)
    sets[["rand_union1"]] <- draw(length(sets[["union"]]), "union", 11L)
    sets[["rand_solo1"]]  <- draw(length(sets[["solo"]]),  "solo",  21L)

    man <- list()
    for (nm in names(sets)) {
        if (!is.null(WANT) && !(nm %in% WANT)) next
        keys <- unique(sets[[nm]])
        if (!length(keys)) {
            message(sprintf("  %-12s EMPTY -- set not written", nm))
            summary_rows[[length(summary_rows) + 1L]] <- data.table(
                seed = s, set = nm, n_snps = 0L, n_causal = 0L, n_linked = 0L,
                n_background = 0L, written = FALSE)
            next
        }
        n <- write_set(proj, nm, keys,
                       if (nm %in% c("best", "union", "solo") || startsWith(nm, "solo_"))
                           pmin_all else NULL)
        nc <- sum(keys %in% TR[category == "causal", key])
        nl <- sum(keys %in% TR[category == "linked_neutral", key])
        nb <- sum(keys %in% TR[category == "background_neutral", key])
        message(sprintf("  %-12s n=%5d  causal=%4d  linked=%5d  background=%5d", nm, n, nc, nl, nb))
        summary_rows[[length(summary_rows) + 1L]] <- data.table(
            seed = s, set = nm, n_snps = n, n_causal = nc, n_linked = nl,
            n_background = nb, written = TRUE)
        man[[length(man) + 1L]] <- list(
            name = nm, n_snps = n, created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
            source_module = "benchmark_offset08",
            threshold_type = if (nm %in% c("best", "union", "solo") || startsWith(nm, "solo_"))
                paste(vapply(POINTS[intersect(SUBSET, names(parts))],
                             function(x) paste(x, collapse = "_"), character(1)),
                      collapse = ";") else "none",
            threshold_value = 0,
            regime = "snp",
            strategy = if (startsWith(nm, "solo_")) "single_method" else
                       switch(nm, best = "at_least_2_5kb", union = "union", solo = "solo",
                              truth = "causal_loci", "random_background"),
            traits = c("bio_1", "bio_2"),
            methods = if (nm %in% c("best", "union")) names(sub_parts)
                      else if (nm == "solo") SOLO
                      else if (startsWith(nm, "solo_")) toupper(sub("^solo_", "", nm))
                      else character(0))
    }
    man_f <- file.path(PIPELINE_ROOT, paste0(proj, "_results"),
                       "_intermediate", "snp_sets", "manifest.json")
    # When only some sets were rebuilt, keep the entries of the ones left untouched --
    # overwriting the manifest with a partial list would hide sets that still exist on disk.
    if (!is.null(WANT) && file.exists(man_f)) {
        old <- tryCatch(fromJSON(man_f, simplifyDataFrame = FALSE), error = function(e) list())
        keep <- Filter(function(e) !(e$name %in% vapply(man, `[[`, character(1), "name")), old)
        man <- c(keep, man)
    }
    write_json(man, man_f, auto_unbox = TRUE, pretty = TRUE)
    manifest_rows[[length(manifest_rows) + 1L]] <- data.table(seed = s, manifest = man_f,
                                                              n_sets = length(man))
}

SUM <- rbindlist(summary_rows)
SUM_F <- file.path(OUTDIR, "snp_sets_summary.tsv")
# Partial rebuild: merge into the existing summary instead of truncating it to the rebuilt sets.
if (!is.null(WANT) && file.exists(SUM_F)) {
    OLD <- fread(SUM_F, colClasses = c("seed" = "character"))
    SUM <- rbind(OLD[!paste(seed, set) %in% SUM[, paste(seed, set)]], SUM, fill = TRUE)
    setorder(SUM, seed, set)
}
fwrite(SUM, SUM_F, sep = "\t")
message("\nWrote ", file.path(OUTDIR, "snp_sets_summary.tsv"))
print(dcast(SUM, seed ~ set, value.var = "n_snps"))

# Structural checks -- these are contracts, not diagnostics.
bad <- SUM[grepl("^rand", set) & (n_causal > 0 | n_linked > 0)]
if (nrow(bad)) { print(bad); stop("Random sets contain causal/linked loci -- the floor is not a null") }
w <- dcast(SUM[written == TRUE], seed ~ set, value.var = "n_snps")
if ("best" %in% names(w) && "union" %in% names(w) && any(w$union < w$best, na.rm = TRUE))
    stop("union smaller than best on some seed -- the agreement rule is inverted")
if ("all" %in% names(w) && any(w$all < w$union, na.rm = TRUE))
    stop("`all` smaller than `union` on some seed -- the universe is not the full tested set")
message("Structural checks passed.")
