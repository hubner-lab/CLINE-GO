#!/usr/bin/env Rscript
# =============================================================================
# mvp_panel_rank.R -- how often each GEA panel is among the best / the worst, per engine.
#
# WHY THIS EXISTS. The claim the simulation has to carry is that the >=2-of-3 panel is the
# STABLE choice: a user who does not know which GEA method suits their data loses least by
# taking it. Distributions of delta (mvp_oracle_stats.R) show the average; they do not show
# how often a panel is the one you would have regretted. A within-replicate RANK does:
#
#   rank the six GEA panels inside each (replicate, offset engine) by accuracy (-tau)
#   -> share of replicates where a panel is in the top k / bottom k, k = 1, 2, 3
#
# Under random ranking every panel sits at k/6 (16.7 / 33.3 / 50 %); that reference is
# drawn on every figure, so separation is read against chance, not as a raw percentage.
# With six panels top-3 and bottom-3 are complements (exactly so under the tie rule below).
#
# Engines are kept APART (GFoffset, LFMM2offset, RDA-uncorrected): they may disagree, and
# which engine favours which panel is itself a result. "median of 3" -- the median accuracy
# over the three engines inside a replicate, as in mvp_oracle_stats.R -- is one more engine
# column, not the base. RDA-corrected is out of the manuscript (Phase 4b gate, 2026-10-03).
#
# Ranking within a replicate does not need the causal-loci oracle (it is constant inside a
# (seed, engine)), so ranks on accuracy == ranks on delta. delta is still written for reference.
#
# EMPTY PANEL = WORST. A panel with no SNPs yields no offset; it is ranked last (decided at
# the gate). Only panels listed as not usable in snp_sets_summary.tsv may be absent -- any
# other gap stops the script.
#
# NOISE. Adjacent panels can differ by less than any meaningful amount, and then who is
# "worst" is a coin flip. rank_gaps.tsv records, for every top-1 / bottom-1 call, the gap to
# the next panel, and rank_summary.tsv the share of those calls decided by a gap < 0.001 and
# < 0.005 (-tau units; panel-vs-oracle deltas are ~0.03).
#
# TIES [changed 2026-10-14, user]. Tied panels share the average rank (frank); every top-k /
# bottom-k share and V4 then count them by the one tie rule of mvp_rank_ties.R: a panel tied
# for ranks k and k + 1 counts 0.5 at each. Before, top-k / bottom-k counted rank <= k /
# rank >= 7 - k, so a tie at 4.5 was never "bottom 2" and one at 5.5 never "bottom 1", while
# V4 drew round(rank). Units without a tie give exactly the counts they gave before (checked).
#
# Outputs (FIG_OUT): rank_long.tsv, rank_summary.tsv (+ bootstrap 95 % CI over replicates),
# rank_gaps.tsv, winrate.tsv, and figure variants V1-V5 (ms_save, 180 mm wide).
#
# OFFSET_DIR and FIG_OUT have no defaults (older arms exist on disk).
#
# Usage:
#   OFFSET_DIR=offset13 FIG_OUT=/pipeline/benchmarks/mvp_eval/figures_ssclines_offset13/rank \
#   MVP_ARM=primary_ssclines MVP_ADDED=<5 block tags> MVP_N_EXPECT=600 \
#     Rscript /pipeline/benchmarks/mvp_panel_rank.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "")
OUT  <- Sys.getenv("FIG_OUT", "")
if (!nzchar(OFF) || !nzchar(OUT)) stop("set OFFSET_DIR and FIG_OUT explicitly (no defaults)")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
source(file.path(ROOT, "benchmarks/ms_style.R"))
source(file.path(ROOT, "benchmarks/mvp_rank_ties.R"))
B_BOOT <- as.integer(Sys.getenv("N_BOOT", "2000"))

rd <- function(...) {
    f <- file.path(...)
    if (!file.exists(f)) stop("MISSING: ", f)
    fread(f, colClasses = c("seed" = "character"))
}
emit <- function(dt, stem) {
    fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
    message(sprintf("  wrote %s.tsv (%d rows)", stem, nrow(dt)))
    invisible(dt)
}

# ------------------------------------------------------------------ vocabulary
# Labels, order and colours as mvp_ms_sim_figures.R:47-51.
PANELS <- c(gea_best = "2/3 methods", gea_union = "1/3 methods", gea_strict = "3/3 methods",
            gea_lfmm_only = "LFMM", gea_rda_only = "RDA", gea_emmax_only = "EMMAX")
SET_OF <- c(gea_best = "best", gea_union = "union", gea_strict = "intersect3",
            gea_lfmm_only = "solo_lfmm", gea_rda_only = "solo_rda", gea_emmax_only = "solo_emmax")
RULES <- unname(PANELS)
RULE_COL <- c("2/3 methods" = MINOU[["teal"]], "1/3 methods" = MINOU[["sage"]],
              "3/3 methods" = MINOU[["navy"]], "LFMM" = MINOU[["red"]],
              "RDA" = MINOU[["amber"]], "EMMAX" = MS_ORANGE)
ENGINES <- c(GFoffset = "gradientForest", LFMM2offset = "LFMM2", `RDA-uncorrected` = "RDA",
             median3 = "median of 3")
ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")   # corpus ships the typo
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
N_PANEL <- length(PANELS)

# ------------------------------------------------------------------ inputs
seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
PRIM  <- mvp_prim(seeds)
stopifnot(nrow(PRIM) == mvp_n_expect())
PRIM[, arch_lab := factor(arch_level, levels = ARCH_LEVELS, labels = ARCH_LABELS)]
stopifnot(!any(is.na(PRIM$arch_lab)))
message(sprintf("arm %s, OFFSET_DIR=%s: %d replicates", mvp_arm_label(), OFF, nrow(PRIM)))

acc <- rd(EVAL, OFF, "phase1_seed_medians_solo.tsv")
acc <- acc[seed %in% PRIM$seed & method_label %in% names(ENGINES)[1:3]]
acc[, accuracy := -tau]                      # NEVER abs(): backwards must stay negative
oracle <- acc[marker_set == "adaptive", .(seed, method_label, oracle = accuracy)]
stopifnot(nrow(oracle) == 3L * nrow(PRIM))
acc <- acc[marker_set %in% names(PANELS)]

# complete grid; an absent panel must be one the builder marked unusable
grid <- CJ(seed = PRIM$seed, method_label = names(ENGINES)[1:3], marker_set = names(PANELS))
acc  <- merge(grid, acc[, .(seed, method_label, marker_set, accuracy)],
              by = c("seed", "method_label", "marker_set"), all.x = TRUE)
miss <- unique(acc[is.na(accuracy), .(seed, set = SET_OF[marker_set])])
if (nrow(miss)) {
    SS <- rd(EVAL, OFF, "snp_sets_summary.tsv")
    chk <- merge(miss, SS[, .(seed, set, usable_for_offset)], by = c("seed", "set"), all.x = TRUE)
    if (any(is.na(chk$usable_for_offset) | chk$usable_for_offset))
        stop("panel absent from the scores but usable per snp_sets_summary:\n",
             paste(capture.output(print(chk)), collapse = "\n"))
    message(sprintf("empty panels ranked worst: %s",
                    paste(sprintf("%s/%s", chk$seed, chk$set), collapse = ", ")))
}
acc[, empty := is.na(accuracy)]

# median-of-3 engine: computed over present panels only; an empty panel stays empty
m3 <- acc[, .(accuracy = if (all(empty)) NA_real_ else median(accuracy), empty = all(empty)),
          by = .(seed, marker_set)][, method_label := "median3"]
acc <- rbind(acc, m3, use.names = TRUE)
oracle <- rbind(oracle, oracle[, .(oracle = median(oracle)), by = seed][, method_label := "median3"],
                use.names = TRUE)
acc <- merge(acc, oracle, by = c("seed", "method_label"))
acc[, delta := accuracy - oracle]
acc <- merge(acc, PRIM[, .(seed, arch_lab, block = added)], by = "seed")
acc[, panel := factor(PANELS[marker_set], levels = RULES)]
acc[, engine := factor(ENGINES[method_label], levels = ENGINES)]

# ------------------------------------------------------------------ ranks
# rank 1 = most accurate; an empty panel gets -Inf, i.e. last. Ties share the average rank here;
# every count below splits them evenly (mvp_rank_ties.R).
acc[, key := fifelse(empty, -Inf, accuracy)]
acc[, rank := frank(-key, ties.method = "average"), by = .(seed, engine)]
n_tied <- acc[, .(tied = any(duplicated(key))), by = .(seed, engine)][, sum(tied)]
message(sprintf("(seed, engine) units with a tie: %d of %d", n_tied, uniqueN(acc[, .(seed, engine)])))
# Tie groups from the raw accuracies (key) must be the groups the tie rule finds from the ranks
# (sets sharing an average rank): the same partition of every unit. Their sizes are reported,
# because a tie of three or more sets averages to a whole number.
acc[, `:=`(g_key = .GRP), by = .(seed, engine, key)][, `:=`(g_rank = .GRP), by = .(seed, engine, rank)]
stopifnot(uniqueN(acc[, .(g_key, g_rank)]) == uniqueN(acc$g_key),
          uniqueN(acc$g_key) == uniqueN(acc$g_rank))
TG <- acc[, .N, by = .(seed, engine, key)][N > 1L]
message(sprintf("tie groups from the raw accuracies: %d (two sets: %d, three or more: %d); non-integer ranks: %d",
                nrow(TG), sum(TG$N == 2L), sum(TG$N > 2L), acc[, sum(rank != round(rank))]))
acc[, c("g_key", "g_rank") := NULL]
RX <- rank_split(acc, c("seed", "engine"))      # (rank, weight) rows, one per tied rank

# gap from each top-1 / bottom-1 panel to the next one (-tau units); empty panels excluded
gaps <- acc[!(empty), {
    o <- order(-key); k <- key[o]
    if (length(k) < 2L) NULL else
        .(position = c("best", "worst"), panel = panel[o][c(1L, length(k))],
          gap = c(k[1] - k[2], k[length(k) - 1L] - k[length(k)]))
}, by = .(seed, engine, arch_lab, block)]
emit(gaps, "rank_gaps")
emit(acc[, .(seed, block, arch_lab, engine, panel, empty, accuracy, oracle, delta, rank)], "rank_long")

# ------------------------------------------------------------------ summary + bootstrap
# [changed 2026-10-14] Indicators under the tie rule: one fractional row per (seed, engine,
# panel), computed once and joined onto every stratum row of LONG in LONG's own row order --
# so each summary row keeps its n and its bootstrap resampling stream, and only tied units move.
METRICS <- c("top1", "top2", "top3", "bottom1", "bottom2", "bottom3", "mean_rank")
EXPECT  <- c(top1 = 1, top2 = 2, top3 = 3, bottom1 = 1, bottom2 = 2, bottom3 = 3) / N_PANEL
IND  <- rank_indicators(RX, c("seed", "engine"), n_set = N_PANEL)
LONG <- rbind(acc[, stratum := "pooled"][],
              copy(acc)[, stratum := as.character(arch_lab)][],
              copy(acc)[, stratum := block][], use.names = TRUE)
I <- IND[LONG[, .(seed, engine, panel, stratum)], on = .(seed, engine, panel)]
stopifnot(nrow(I) == nrow(LONG), !anyNA(I[, ..METRICS]))
setcolorder(I, c("seed", "engine", "panel", "stratum"))

# Check against the previous counts (rank <= k / rank >= 7 - k): units without a tie give exactly
# the same indicator; in a tied unit a panel moves by 0.5 where the cutoff falls inside its tie
# (the old count put the average rank on the wrong side), never otherwise; the mean rank is the
# average rank either way.
OLD <- acc[, .(seed, engine, panel,
               top1 = as.numeric(rank <= 1), top2 = as.numeric(rank <= 2), top3 = as.numeric(rank <= 3),
               bottom1 = as.numeric(rank >= N_PANEL), bottom2 = as.numeric(rank >= N_PANEL - 1),
               bottom3 = as.numeric(rank >= N_PANEL - 2), mean_rank = rank)]
CMP <- merge(OLD, IND, by = c("seed", "engine", "panel"), suffixes = c(".old", ".new"))
stopifnot(nrow(CMP) == nrow(acc))
CMP <- merge(CMP, RX[, .(tied = any(tie_n > 1L)), by = .(seed, engine)], by = c("seed", "engine"))
MOVED <- rbindlist(lapply(METRICS, function(v) {
    d <- CMP[[paste0(v, ".new")]] - CMP[[paste0(v, ".old")]]
    if (any(d[!CMP$tied] != 0)) stop(v, ": a unit without a tie changed its count")
    if (!all(abs(d[CMP$tied]) %in% c(0, 0.5))) stop(v, ": a tied panel moved by other than 0 or 0.5")
    if (v == "mean_rank" && any(abs(d) > 1e-12)) stop("mean rank changed")
    CMP[d != 0, .(metric = v, change = d[d != 0], seed, engine, panel)]
}))
message(sprintf("tie rule vs previous counts: untied units identical; %d panel x metric values moved by 0.5 (all in tied units)",
                nrow(MOVED)))
if (nrow(MOVED)) print(dcast(MOVED[, .N, by = .(engine, metric)], engine ~ metric, value.var = "N", fill = 0L))

set.seed(20261003)
boot_ci <- function(m) {
    # m: replicates x 1 numeric; resample replicates
    n <- length(m)
    b <- vapply(seq_len(B_BOOT), function(i) mean(m[sample.int(n, n, replace = TRUE)]), 0)
    quantile(b, c(.025, .975), names = FALSE)
}
SUMM <- I[, {
    out <- lapply(METRICS, function(v) {
        x <- as.numeric(.SD[[v]]); ci <- boot_ci(x)
        data.table(metric = v, value = mean(x), lo = ci[1], hi = ci[2])
    })
    rbindlist(out)
}, by = .(engine, panel, stratum), .SDcols = METRICS]
SUMM[, n := I[, .N, by = .(engine, panel, stratum)][SUMM, on = .(engine, panel, stratum), N]]
SUMM[metric != "mean_rank", `:=`(value = 100 * value, lo = 100 * lo, hi = 100 * hi,
                                 expected = 100 * EXPECT[metric])]
SUMM[metric == "mean_rank", expected := (N_PANEL + 1) / 2]
# Invariants only the split rule satisfies exactly (each panel spreads one unit of weight over the
# six ranks, and each rank holds one unit of weight): top 3 + bottom 3 = 100 % on every row, and
# per engine x stratum the six panels' top-k (bottom-k) shares sum to k x 100 %.
TB <- dcast(SUMM[metric %in% c("top3", "bottom3")], engine + panel + stratum ~ metric, value.var = "value")
KS <- SUMM[metric %like% "^(top|bottom)[123]$", .(s = sum(value)), by = .(engine, stratum, metric)]
stopifnot(max(abs(TB$top3 + TB$bottom3 - 100)) < 1e-9,
          max(abs(KS$s - 100 * as.integer(substring(KS$metric, nchar(KS$metric))))) < 1e-9)
message("invariants: top3 + bottom3 = 100 % on every row; top-k and bottom-k sum to k x 100 % per engine x stratum")

# share of top-1 / bottom-1 calls decided by a hair
hair <- gaps[, .(calls = .N, gap_median = median(gap),
                 lt_0.001_pct = 100 * mean(gap < 0.001), lt_0.005_pct = 100 * mean(gap < 0.005)),
             by = .(engine, panel, position)]
emit(hair, "rank_gap_summary")
setorder(SUMM, stratum, engine, metric, panel)
emit(SUMM, "rank_summary")

message("\n=== % replicates in top-k / bottom-k, pooled (expected by chance in brackets) ===")
for (v in c("top1", "bottom1", "top2", "bottom2", "top3")) {
    message(sprintf("-- %s [%.1f %%]", v, 100 * EXPECT[[v]]))
    print(dcast(SUMM[stratum == "pooled" & metric == v], panel ~ engine, value.var = "value"),
          digits = 3)
}
message("\n=== top-1 / bottom-1 calls decided by a gap < 0.005 (% of that panel's calls) ===")
print(dcast(hair[engine == "median of 3"], panel ~ position, value.var = "lt_0.005_pct"), digits = 3)

# ------------------------------------------------------------------ pairwise win rate
WIN <- rbindlist(lapply(RULES, function(a) rbindlist(lapply(setdiff(RULES, a), function(b) {
    x <- dcast(acc[panel %in% c(a, b)], seed + engine ~ panel, value.var = "key")
    x[, .(panel_a = a, panel_b = b, n = .N,
          win_pct = 100 * mean((get(a) > get(b)) + 0.5 * (get(a) == get(b)))), by = engine]
}))))
emit(WIN, "winrate")

# ------------------------------------------------------------------ figures
P <- SUMM[stratum == "pooled"]
P[, panel := factor(panel, levels = rev(RULES))]
P[, engine := factor(engine, levels = ENGINES)]
PCT <- function(x) paste0(abs(x), "%")

# V1-V3: diverging bars, bottom-k to the left, top-k to the right, one figure per k
diverging <- function(k) {
    d <- rbind(P[metric == paste0("top", k)][, side := "top"],
               P[metric == paste0("bottom", k)][, side := "bottom"])
    d[side == "bottom", `:=`(value = -value, lo = -hi, hi = -lo)]
    e <- 100 * k / N_PANEL
    lim <- max(abs(c(d$lo, d$hi)), e) * 1.12
    ggplot(d, aes(value, panel, fill = panel, alpha = side)) +
        geom_vline(xintercept = c(-e, e), linetype = "dashed", colour = MS_REF, linewidth = 0.3) +
        geom_vline(xintercept = 0, colour = MS_INK, linewidth = 0.3) +
        geom_col(width = 0.72) +
        geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.25, colour = MS_INK,
                      alpha = 1) +
        geom_text(aes(x = ifelse(value < 0, lo, hi), label = sprintf("%.0f", abs(value)),
                      hjust = ifelse(value < 0, 1.25, -0.25)),
                  size = MS_LAB, colour = MS_INK, alpha = 1) +
        facet_wrap(~ engine, nrow = 1) +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        scale_alpha_manual(values = c(top = 1, bottom = 0.45), guide = "none") +
        scale_x_continuous(labels = PCT, limits = c(-lim, lim)) +
        labs(x = sprintf("%% replicates in the bottom %d  |  top %d of 6 panels (dashed: chance, %.0f %%)",
                         k, k, e), y = NULL) +
        theme_ms() + theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25))
}
for (k in 1:3) ms_save(file.path(OUT, sprintf("V%d_rank_top%d", k, k)), diverging(k), 180, 62)

# V4: full rank distribution, ranks 1-6 stacked; a tied panel counts half at each of its ranks
# [changed 2026-10-14: was round(rank), half to even]
R4 <- RX[, .(n = sum(w)), by = .(engine, panel, r)]
R4[, pct := 100 * n / sum(n), by = .(engine, panel)]
R4[, r := factor(r, levels = N_PANEL:1)]
R4[, panel := factor(panel, levels = rev(RULES))]
RANK_COL <- setNames(colorRampPalette(c(MINOU[["red"]], "#F2F2F2", MINOU[["teal"]]))(N_PANEL),
                     N_PANEL:1)
v4 <- ggplot(R4, aes(pct, panel, fill = r)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    facet_wrap(~ engine, nrow = 1) +
    scale_fill_manual(values = RANK_COL, name = "rank (1 = most accurate)",
                      breaks = as.character(1:N_PANEL)) +
    scale_x_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    labs(x = "% replicates at each rank", y = NULL) + theme_ms()
ms_save(file.path(OUT, "V4_rank_distribution"), v4, 180, 66)

# V5: pairwise win rate, row panel beats column panel
W5 <- copy(WIN)
W5[, panel_a := factor(panel_a, levels = rev(RULES))][, panel_b := factor(panel_b, levels = RULES)]
W5[, engine := factor(engine, levels = ENGINES)]
v5 <- ggplot(W5, aes(panel_b, panel_a, fill = win_pct)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f", win_pct), colour = abs(win_pct - 50) > 12),
              size = MS_LAB) +
    facet_wrap(~ engine, nrow = 1) +
    scale_fill_gradient2(low = MS_DIVERGING[["low"]], mid = MS_DIVERGING[["mid"]],
                         high = MS_DIVERGING[["high"]], midpoint = 50, limits = c(20, 80),
                         oob = scales::squish, name = "% replicates row beats column") +
    scale_colour_manual(values = c(`TRUE` = "#FFFFFF", `FALSE` = MS_INK), guide = "none") +
    labs(x = NULL, y = NULL) + theme_ms() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(),
          axis.ticks = element_blank())
ms_save(file.path(OUT, "V5_winrate"), v5, 180, 62)

message("\nAll tables and figures written to ", OUT)
