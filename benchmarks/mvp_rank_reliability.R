#!/usr/bin/env Rscript
# =============================================================================
# mvp_rank_reliability.R -- is the within-replicate rank of the GEA panels signal or noise?
#
# WHY THIS EXISTS. The six panels' offset-accuracy distributions over 600 replicates overlap
# almost completely (medians within 0.012), so a reader can suspect that ranking panels inside a
# replicate (mvp_panel_rank.R, the main offset figure) orders noise. Three tests, reported in
# the manuscript TEXT only (Methods + Results; user decision 2026-10-03, no figure):
#
#   1. SPLIT-HALF RELIABILITY. The 100 landscape gardens of a replicate are split at random into
#      two halves of 50; each panel's accuracy (median -tau) is recomputed on each half and
#      centred on the replicate's mean over its six panels. Pearson r of the centred values
#      between halves (pooled over replicates x panels), Spearman-Brown-corrected to 100 gardens,
#      and the share of replicates whose LEAST (and most) accurate panel is the same in both
#      halves (chance 1/6). N_SPLITS random splits; mean and 2.5-97.5 % over splits.
#   2. ENGINE AGREEMENT. Kendall's W of the six panels' order across the three offset engines
#      (GF, LFMM2, RDA-uncorrected) per replicate, full gardens; permutation null = each engine's
#      panel values shuffled independently within the replicate (N_PERM). Plus pairwise engine
#      correlation of the centred values and agreement on the least accurate panel.
#   3. VARIANCE PARTITION. Accuracy (median of 3 engines) ~ replicate + panel, sums of squares:
#      the panel main effect vs the replicate effect vs the replicate x panel residual.
#
# Test 1 rules out garden-sampling noise, test 2 engine-specific noise (LFMM2 and RDA are
# deterministic given the panel; only GF is stochastic). Neither tests whether 2/3 ranks mid-field
# because it is built from the other panels' SNPs -- that is a framing question (a hedge), not one
# these data can settle.
#
# Replicates: those with all six panels scored (an empty panel has no garden rows); landscape
# gardens only, as every accuracy in the manuscript. Reads offset13 tables only, fits nothing.
#
# Outputs (FIG_OUT): reliability_split_half.tsv, reliability_engines.tsv, reliability_variance.tsv
# OFFSET_DIR and FIG_OUT have no defaults (older arms exist on disk).
#
# Usage (work/p5_figures.sh reliability):
#   OFFSET_DIR=offset13 FIG_OUT=/pipeline/benchmarks/mvp_eval/figures_ssclines_offset13/rank \
#   MVP_ARM=primary_ssclines MVP_ADDED=<5 block tags> MVP_N_EXPECT=600 \
#     Rscript /pipeline/benchmarks/mvp_rank_reliability.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "")
OUT  <- Sys.getenv("FIG_OUT", "")
if (!nzchar(OFF) || !nzchar(OUT)) stop("set OFFSET_DIR and FIG_OUT explicitly (no defaults)")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
N_SPLITS <- as.integer(Sys.getenv("N_SPLITS", "100"))
N_PERM   <- as.integer(Sys.getenv("N_PERM", "1000"))
SEED     <- 20261003L

PANELS  <- c(gea_best = "2/3 methods", gea_union = "1/3 methods", gea_strict = "3/3 methods",
             gea_lfmm_only = "LFMM", gea_rda_only = "RDA", gea_emmax_only = "EMMAX")
ENGINES <- c(GFoffset = "gradientForest", LFMM2offset = "LFMM2", `RDA-uncorrected` = "RDA")
N_PANEL <- length(PANELS)
emit <- function(dt, stem) {
    fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
    message(sprintf("  wrote %s.tsv (%d rows)", stem, nrow(dt)))
}

# ------------------------------------------------------------------ inputs
PRIM <- mvp_prim(fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character")))
stopifnot(nrow(PRIM) == mvp_n_expect())
G <- fread(file.path(EVAL, OFF, "garden_performance.tsv"),
           select = c("seed", "garden_id", "garden_type", "method_label", "marker_set", "tau"),
           colClasses = c(seed = "character"))
G <- G[seed %in% PRIM$seed & garden_type == "landscape" & method_label %in% names(ENGINES) &
       marker_set %in% names(PANELS)]
G[, acc := -tau]                                       # never abs(): backwards stays negative
full <- G[, uniqueN(marker_set), by = seed][V1 == N_PANEL, seed]
G <- G[seed %in% full]
NG <- G[, uniqueN(garden_id), by = seed]
stopifnot(length(full) >= nrow(PRIM) - 10L, all(NG$V1 == 100L),
          G[, .N, by = .(seed, method_label, marker_set)][, all(N == 100L)])
message(sprintf("%d of %d replicates have all %d panels scored; 100 landscape gardens each",
                length(full), nrow(PRIM), N_PANEL))

# Tie: medians over all 100 gardens equal the accuracy table every figure uses.
F <- G[, .(acc = median(acc)), by = .(seed, method_label, marker_set)]
P1 <- fread(file.path(EVAL, OFF, "phase1_seed_medians_solo.tsv"), colClasses = c(seed = "character"))
tie <- merge(F, P1[, .(seed, method_label, marker_set, a1 = -tau)],
             by = c("seed", "method_label", "marker_set"), all.x = TRUE)
stopifnot(!anyNA(tie$a1), max(abs(tie$acc - tie$a1)) < 1e-12)
message(sprintf("tie vs phase1_seed_medians_solo.tsv: %d rows, max |diff| %.1e", nrow(tie),
                max(abs(tie$acc - tie$a1))))

# ------------------------------------------------------------------ 1. split-half reliability
add_m3 <- function(A, by) rbind(A, A[, .(acc = median(acc)), by = by][, method_label := "median3"],
                                use.names = TRUE)
set.seed(SEED)
GD <- unique(G[, .(seed, garden_id)])
SH <- rbindlist(lapply(seq_len(N_SPLITS), function(s) {
    GD[, half := sample(rep(c("A", "B"), length.out = .N)), by = seed]
    H <- G[GD, on = .(seed, garden_id)]
    A <- add_m3(H[, .(acc = median(acc)), by = .(seed, method_label, marker_set, half)],
                c("seed", "marker_set", "half"))
    A[, c := acc - mean(acc), by = .(seed, method_label, half)]
    W <- dcast(A, seed + method_label + marker_set ~ half, value.var = "c")
    E <- W[, .(worst = which.min(A) == which.min(B), best = which.max(A) == which.max(B)),
           by = .(seed, method_label)]
    merge(W[, .(r = cor(A, B)), by = method_label],
          E[, .(worst_same = mean(worst), best_same = mean(best)), by = method_label],
          by = "method_label")[, split := s]
}))
SH[, r_full := 2 * r / (1 + r)]                        # Spearman-Brown, 50 -> 100 gardens
q <- function(x, p) unname(quantile(x, p))
SHS <- SH[, .(n_rep = length(full), n_splits = .N,
              r_half = mean(r), r_half_lo = q(r, .025), r_half_hi = q(r, .975),
              r_full = mean(r_full),
              worst_same_pct = 100 * mean(worst_same), worst_same_lo = 100 * q(worst_same, .025),
              worst_same_hi = 100 * q(worst_same, .975),
              best_same_pct = 100 * mean(best_same), chance_pct = 100 / N_PANEL), by = method_label]
SHS[, engine := c(ENGINES, median3 = "median of 3")[method_label]]
emit(SHS, "reliability_split_half")
print(SHS[, .(engine, r_half, r_full, worst_same_pct, best_same_pct, chance_pct)], digits = 3)

# ------------------------------------------------------------------ 2. engine agreement
# array replicate x panel x engine of full-garden accuracy
F[, c := acc - mean(acc), by = .(seed, method_label)]
ARR <- array(NA_real_, c(length(full), N_PANEL, length(ENGINES)),
             dimnames = list(full, names(PANELS), names(ENGINES)))
ARR[cbind(match(F$seed, full), match(F$marker_set, names(PANELS)), match(F$method_label, names(ENGINES)))] <- F$acc
stopifnot(!anyNA(ARR))
kendall_w <- function(M) {                              # rows = panels, columns = engines
    R <- apply(M, 2, rank); m <- ncol(R); n <- nrow(R)
    12 * sum((rowSums(R) - m * (n + 1) / 2)^2) / (m^2 * (n^3 - n))
}
W_obs <- apply(ARR, 1, kendall_w)
set.seed(SEED)
W_null <- vapply(seq_len(N_PERM), function(i)
    mean(apply(ARR, 1, function(M) kendall_w(apply(M, 2, sample)))), 0)
p_perm <- (1 + sum(W_null >= mean(W_obs))) / (1 + N_PERM)
pairs <- combn(names(ENGINES), 2, simplify = FALSE)
PW <- rbindlist(lapply(pairs, function(p) {
    a <- ARR[, , p[1]]; b <- ARR[, , p[2]]
    ca <- a - rowMeans(a); cb <- b - rowMeans(b)
    data.table(engine_a = ENGINES[[p[1]]], engine_b = ENGINES[[p[2]]],
               r_centred = cor(as.vector(ca), as.vector(cb)),
               worst_same_pct = 100 * mean(apply(a, 1, which.min) == apply(b, 1, which.min)),
               best_same_pct = 100 * mean(apply(a, 1, which.max) == apply(b, 1, which.max)))
}))
ENG <- rbind(
    data.table(test = "kendall_w", engine_a = "all three", engine_b = NA_character_,
               value = mean(W_obs), median = median(W_obs), null_mean = mean(W_null),
               null_lo = q(W_null, .025), null_hi = q(W_null, .975), expected_no_agreement = 1 / length(ENGINES),
               p_perm = p_perm, n_rep = length(full), n_perm = N_PERM),
    PW[, .(test = "pair", engine_a, engine_b, value = r_centred, worst_same_pct, best_same_pct,
           chance_pct = 100 / N_PANEL, n_rep = length(full))], use.names = TRUE, fill = TRUE)
emit(ENG, "reliability_engines")
message(sprintf("Kendall W: mean %.3f (median %.3f); permutation null mean %.3f [%.3f, %.3f]; p = %.2g",
                mean(W_obs), median(W_obs), mean(W_null), q(W_null, .025), q(W_null, .975), p_perm))
print(PW, digits = 3)

# ------------------------------------------------------------------ 3. variance partition
M3 <- F[, .(acc = median(acc)), by = .(seed, marker_set)]
fit <- aov(acc ~ factor(seed) + factor(marker_set), data = M3)
ss <- summary(fit)[[1]][["Sum Sq"]]
VP <- data.table(source = c("replicate", "panel", "replicate x panel (residual)"),
                 ss = ss, share_pct = 100 * ss / sum(ss), n_rep = length(full),
                 engine = "median of 3")
pr <- M3[, .(med = median(acc)), by = marker_set]
VP[, `:=`(panel_median_range = diff(range(pr$med)),
          within_rep_range_median = M3[, diff(range(acc)), by = seed][, median(V1)])]
emit(VP, "reliability_variance")
print(VP, digits = 3)
message("\nAll tables written to ", OUT)
