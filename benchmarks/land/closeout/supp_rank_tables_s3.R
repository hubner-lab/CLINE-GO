#!/usr/bin/env Rscript
# =============================================================================
# mvp_ms_supp_rank_tables.R -- supplementary rank tables of the simulation section (SS-Clines, 600
# replicates), built from Data S2 only (supplement build 2026-10-25, item A3):
#   t_ranks            % of replicates at each within-replicate rank 1..6 (1 = most accurate offset),
#                      ranks 1-2 and ranks 5-6, per offset engine and candidate set, 95 % percentile
#                      bootstrap CI
#   t_ranks_strata     the same under the median of the three engines, per genetic architecture (3),
#                      demographic scenario (5) and design cell (12); bootstrap within the stratum
#   t_offset_accuracy  offset accuracy of the set minus that of the causal loci of the same
#                      replicate (median of the three engines): median, min, max, % >= 0
#                      [audit3b] plus the set's own accuracy (median, min, max) and a last row with the
#                      causal loci's accuracy (one `oracle` per replicate), difference cells empty
#   t_rank_tests       Kruskal-Wallis tests of each set's rank across design cells and across
#                      demographic scenarios, Holm over the 12 tests
# Every value is recomputed from the Data S2 members the way mvp_ms_sim_figures.R (SIM) computes it,
# then asserted against SIM's own outputs (Table_S1.tsv, numbers.tsv); no number is copied from them.
#   * ties: benchmarks/mvp_rank_ties.R, rank_split(RL, c("seed", "engine")) as SIM:565 -- a set tied
#     for ranks k and k + 1 counts its weight 1/2 at each
#   * bootstrap as SIM:955-973: replicates in sort(unique(seed), method = "radix") order,
#     set.seed(20261005), 10,000 iterations, each drawing i <- sample.int(600, 600, replace = TRUE)
#     and no other random number; every statistic of an iteration from the same i
#   * strata as SIM:854-869 (Figure R3_least_accurate_strata): set.seed(20261005) before each
#     stratum, strata in the order architectures, scenarios, design cells, each in SIM's level order
#   * rows in the Figure 2b order SIM computes (PANEL_ORDER, recomputed and asserted)
#
# Usage: Rscript benchmarks/mvp_ms_supp_rank_tables.R <Data_S2.zip> <Table_S1.tsv> <numbers.tsv> <out_dir>
# Writes <out_dir>/{t_ranks,t_ranks_strata,t_offset_accuracy,t_rank_tests}.{tsv,md}.
# Driver (2026-10-25): cline-go:latest, OPENBLAS_NUM_THREADS=2 (as work/p5_figures.sh).
# =============================================================================
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: mvp_ms_supp_rank_tables.R <Data_S2.zip> <Table_S1.tsv> <numbers.tsv> <out_dir>")
zip <- args[1]; ts1_path <- args[2]; num_path <- args[3]; out_dir <- args[4]
for (f in c(zip, ts1_path, num_path)) if (!file.exists(f)) stop("not found: ", f)
if (!dir.exists(out_dir)) stop("output directory not found: ", out_dir)
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(ROOT, "benchmarks/mvp_rank_ties.R"))

# ------------------------------------------------------------------ vocabulary (SIM:127-146)
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
BLOCK_LAB <- c("N variable, m variable", "N cline north-south", "N equal, m constant",
               "N cline centre-edge", "N equal, m breaks")
RULES <- c("2/3 methods", "1/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")
ENG_LAB <- c(gradientForest = "Gradient Forest", LFMM2 = "Geometric offset", RDA = "RDA-based offset",
             `median of 3` = "Median of the three")
SET_LAB <- c("2/3 methods" = "2/3 set", "1/3 methods" = "1/3 set", "3/3 methods" = "3/3 set",
             "LFMM" = "LFMM2", "RDA" = "RDA", "EMMAX" = "EMMAX")   # [closeout] display label LFMM2; the data key stays LFMM
set_lab <- function(x) unname(SET_LAB[as.character(x)])
N_PANEL <- length(RULES)
NREP <- 600L
S1_SEED <- 20261005L
S1_B <- 10000L
cell_of <- function(a, p, r) paste(a, p, r, sep = ", ")                      # SIM:862

# ------------------------------------------------------------------ inputs (Data S2 members)
read_member <- function(member, ...) fread(cmd = paste("unzip -p", shQuote(zip), shQuote(member)), sep = "\t", ...)
RL  <- read_member("figures_ld_offset14/rank/rank_long.tsv", colClasses = c(seed = "character"))
DET <- read_member("figures_ld_gallery14/detection_per_seed.tsv", colClasses = c(seed = "character"))[label %in% RULES]
OFF <- read_member("figures_ld_gallery14/offset_delta_per_seed.tsv", colClasses = c(seed = "character"))[label %in% RULES]
fac <- function(D) {                                                          # SIM:213-221
    D[, arch := factor(arch, levels = ARCH_LABELS)]
    if ("regime" %in% names(D)) D[, regime := factor(regime, levels = c("equal-S", "unequal-S"))]
    if ("pleio"  %in% names(D)) D[, pleio  := factor(pleio,  levels = c("no pleiotropy", "pleiotropy"))]
    if ("block"  %in% names(D)) D[, block  := factor(block,  levels = BLOCK_LAB)]
    D[, label := factor(label, levels = RULES)]
    stopifnot(!anyNA(D$arch), !anyNA(D$label))
    D[]
}
DET <- fac(DET); OFF <- fac(OFF)
stopifnot(setequal(unique(RL$engine), names(ENG_LAB)), setequal(unique(RL$panel), RULES),
          uniqueN(RL$seed) == NREP, nrow(RL) == NREP * length(ENG_LAB) * N_PANEL,
          !anyDuplicated(RL[, .(seed, engine, panel)]), setequal(RL$seed, DET$seed))
# tie groups found from the ranks are the ties of the raw accuracies (SIM:552, 562-564)
RL[, key := fifelse(empty, -Inf, accuracy)]
RL[, g_key := .GRP, by = .(seed, engine, key)][, g_rank := .GRP, by = .(seed, engine, rank)]
stopifnot(uniqueN(RL[, .(g_key, g_rank)]) == uniqueN(RL$g_key), uniqueN(RL$g_key) == uniqueN(RL$g_rank))
RL[, c("g_key", "g_rank") := NULL]

# ------------------------------------------------------------------ ranks as (rank, weight) rows
RX <- rank_split(RL, c("seed", "engine"))
stopifnot(RX[tie_n == 1L, all(r == rank & w == 1)])
TI <- rank_indicators(RX, c("seed", "engine"), n_set = N_PANEL)
ORD <- TI[engine == "median of 3", .(bottom1 = 100 * mean(bottom1), bottom2 = 100 * mean(bottom2)), by = panel]
setorder(ORD, bottom1, bottom2)
PANEL_ORDER <- ORD$panel
message("[LD-lab] PANEL_ORDER: ", paste(PANEL_ORDER, collapse = " < "))
message("set order (Figure 2b): ", paste(set_lab(PANEL_ORDER), collapse = " < "))

STATS <- c(r1 = "rank 1", r2 = "rank 2", r3 = "rank 3", r4 = "rank 4", r5 = "rank 5", r6 = "rank 6",
           r12 = "ranks 1–2", r56 = "ranks 5–6")
ST <- RX[, .(r1 = sum(w[r == 1L]), r2 = sum(w[r == 2L]), r3 = sum(w[r == 3L]), r4 = sum(w[r == 4L]),
             r5 = sum(w[r == 5L]), r6 = sum(w[r == 6L]), r12 = sum(w[r <= 2L]), r56 = sum(w[r >= 5L])),
         by = .(seed, engine, panel)]
stopifnot(nrow(ST) == NREP * length(ENG_LAB) * N_PANEL)

# One matrix: rows = replicates in radix seed order, columns = every (engine, statistic, set)
build_big <- function(S, engines, seeds) {
    idx <- CJ(engine = engines, stat = names(STATS), panel = RULES, sorted = FALSE)
    M <- do.call(cbind, lapply(engines, function(e) do.call(cbind, lapply(names(STATS), function(m) {
        w <- dcast(S[engine == e], seed ~ panel, value.var = m)
        stopifnot(identical(w$seed, seeds), !anyNA(w))
        as.matrix(w[, RULES, with = FALSE])
    }))))
    stopifnot(ncol(M) == nrow(idx))
    list(M = M, idx = idx)
}
boot_ci <- function(M, n) {
    set.seed(S1_SEED)
    B <- matrix(NA_real_, S1_B, ncol(M))
    for (b in seq_len(S1_B)) {
        i <- sample.int(n, n, replace = TRUE)
        B[b, ] <- colMeans(M[i, , drop = FALSE])
    }
    data.table(value = 100 * colMeans(M),
               lo = apply(100 * B, 2, quantile, 0.025, names = FALSE),
               hi = apply(100 * B, 2, quantile, 0.975, names = FALSE))
}

# ==================================================================== t_ranks
ENGS <- names(ENG_LAB)
SEEDS <- sort(unique(ST$seed), method = "radix")
stopifnot(length(SEEDS) == NREP)
BG <- build_big(ST, ENGS, SEEDS)
TR <- cbind(BG$idx, boot_ci(BG$M, NREP))
TR[, n := NREP]
# shares: each set's six ranks sum to 100, and so does each rank over the six sets
s_set  <- TR[stat %in% sprintf("r%d", 1:6), sum(value), by = .(engine, panel)]
s_rank <- TR[stat %in% sprintf("r%d", 1:6), sum(value), by = .(engine, stat)]
chk <- c(sum_per_set = max(abs(s_set$V1 - 100)) < 1e-9, sum_per_rank = max(abs(s_rank$V1 - 100)) < 1e-9)
# against Table_S1.tsv section a (value, lo and hi), |difference| < 1e-9
TS1 <- fread(ts1_path, sep = "\t")
stopifnot(TS1[set == "LFMM", .N] == 22L)   # [closeout] the original Table_S1.tsv: LFMM in sections a (16), b (4), e (2)
TS1[set == "LFMM", set := "LFMM2"]   # [closeout] its set labels follow SET_LAB, so every tie below still pairs all rows
A_MET <- c(r6 = "% replicates least accurate", r56 = "% replicates among the two least accurate",
           r1 = "% replicates most accurate", r12 = "% replicates among the two most accurate")
xa <- merge(TR[stat %in% names(A_MET), .(set = set_lab(panel), engine = unname(ENG_LAB[engine]),
                                          metric = unname(A_MET[stat]), value, lo, hi)],
            TS1[section == "a", .(set, engine, metric, v1 = value, l1 = lo, h1 = hi)], by = c("set", "engine", "metric"))
chk["table_s1_a_rows"] <- nrow(xa) == 4L * 4L * N_PANEL && TS1[section == "a", .N] == nrow(xa)
dA <- xa[, max(abs(c(value - v1, lo - l1, hi - h1)))]
chk["table_s1_a_values"] <- dA < 1e-9
message(sprintf("t_ranks vs Table_S1.tsv (a): %d rows, max |difference| over value/lo/hi %.2e", nrow(xa), dA))
ex <- function(s, st) TR[engine == "median of 3" & panel == s & stat == st]
near3 <- function(x, v, l, h) abs(x$value - v) < 0.005 && abs(x$lo - l) < 0.005 && abs(x$hi - h) < 0.005
chk["ex_23_least"]  <- near3(ex("2/3 methods", "r6"), 2.33, 1.17, 3.67)   # [closeout] Phase 5 values (t_ranks.tsv)
chk["ex_emmax_least"] <- near3(ex("EMMAX", "r6"), 29.17, 25.50, 33.00)   # [closeout] Phase 5 values
chk["ex_23_last2"] <- near3(ex("2/3 methods", "r56"), 13.17, 10.50, 15.83)   # [closeout] Phase 5 values
print(chk); if (!all(chk)) warning("[LD-lab] t_ranks assertion failed: ", paste(names(chk)[!chk], collapse = ", "))
message("t_ranks asserts passed: ", paste(names(chk), collapse = ", "))

# ==================================================================== t_ranks_strata
COV <- unique(DET[, .(seed, arch, pleio, regime, block)])
stopifnot(nrow(COV) == NREP, !anyNA(COV), !anyDuplicated(COV$seed))
S3 <- merge(ST[engine == "median of 3"], COV, by = "seed")
stopifnot(nrow(S3) == NREP * N_PANEL)
S3[, cell := cell_of(arch, pleio, regime)]
CELLS <- CJ(arch = ARCH_LABELS, pleio = c("no pleiotropy", "pleiotropy"), regime = c("equal-S", "unequal-S"),
            sorted = FALSE)[, cell_of(arch, pleio, regime)]
STRATA <- rbind(data.table(kind = "genetic architecture", stratum = ARCH_LABELS, col = "arch"),
                data.table(kind = "demographic scenario", stratum = BLOCK_LAB, col = "block"),
                data.table(kind = "design cell", stratum = CELLS, col = "cell"))
TRS <- rbindlist(lapply(seq_len(nrow(STRATA)), function(k) {
    s <- STRATA[k]
    D <- S3[as.character(get(s$col)) == s$stratum]
    sd <- sort(unique(D$seed), method = "radix")
    G <- build_big(D, "median of 3", sd)
    cbind(kind = s$kind, stratum = s$stratum, G$idx[, .(stat, panel)], boot_ci(G$M, length(sd)), n = length(sd))
}))
stopifnot(nrow(TRS) == nrow(STRATA) * length(STATS) * N_PANEL)
message("strata bootstrap order: ", paste(sprintf("%s: %s", STRATA$kind, STRATA$stratum), collapse = " | "))
NUM <- fread(num_path, sep = "\t")
N8 <- NUM[startsWith(figure, "Figure S8 / Table S1 (c)") & engine == "median of 3",
          .(stratum_key = stratum, panel = rule, stat = c(least_accurate_pct = "r6", last_two_pct = "r56")[metric], v8 = value, n8 = n)]
x8 <- merge(TRS[stat %in% c("r6", "r56"), .(stratum_key = paste0(kind, ": ", stratum), panel, stat, value, n)],
            N8, by = c("stratum_key", "panel", "stat"))
chk2 <- c(numbers_rows = nrow(x8) == nrow(STRATA) * N_PANEL * 2L && nrow(N8) == nrow(x8),
          numbers_values = max(abs(x8$value - x8$v8)) < 1e-9, numbers_n = all(x8$n == x8$n8),
          n_per_stratum = TRS[, all(n == c(`genetic architecture` = 200, `demographic scenario` = 120, `design cell` = 50)[kind])])
message(sprintf("t_ranks_strata vs numbers.tsv (Figure S8 / Table S1 (c)): %d rows, max |difference| %.2e",
                nrow(x8), max(abs(x8$value - x8$v8))))
half_up <- function(v, d) floor(v * 10^d + 0.5) / 10^d                       # SIM:1018
RNG <- TRS[kind == "design cell" & stat == "r56" & panel %in% c("1/3 methods", "2/3 methods"),
           .(lo = min(value), hi = max(value)), by = .(regime = sub(".*, ", "", stratum), panel)]
WANT <- data.table(regime = c("unequal-S", "unequal-S", "equal-S", "equal-S"),
                   panel = c("1/3 methods", "2/3 methods", "1/3 methods", "2/3 methods"),
                   want_lo = c(27, 8, 6, 10), want_hi = c(52, 14, 12, 26))   # [closeout] Phase 5 values (t_ranks_strata.tsv)
RNG <- merge(RNG, WANT, by = c("regime", "panel"))
RNG[, ok := half_up(lo, 0) == want_lo & half_up(hi, 0) == want_hi]
for (i in seq_len(nrow(RNG)))
    message(sprintf("two least accurate by design cell, %s, %s: %.4f-%.4f %% (stated %g-%g), n = 50 per cell",
                    RNG$regime[i], set_lab(RNG$panel[i]), RNG$lo[i], RNG$hi[i], RNG$want_lo[i], RNG$want_hi[i]))
chk2["design_cell_ranges"] <- nrow(RNG) == 4L && all(RNG$ok)
print(chk2); if (!all(chk2)) warning("[LD-lab] t_ranks_strata assertion failed: ", paste(names(chk2)[!chk2], collapse = ", "))
message("t_ranks_strata asserts passed: ", paste(names(chk2), collapse = ", "))

# ==================================================================== t_offset_accuracy
DB <- OFF[, .(n = .N, median = median(delta), min = min(delta), max = max(delta), ge0_pct = 100 * mean(delta >= 0)),
          by = .(panel = as.character(label))]
B_MET <- c(median = "median delta accuracy", min = "min delta accuracy", max = "max delta accuracy",
           ge0_pct = "% replicates delta >= 0")
xb <- merge(melt(DB, id.vars = c("panel", "n"), variable.name = "stat", variable.factor = FALSE)[
                , .(set = set_lab(panel), metric = unname(B_MET[stat]), n, value)],
            TS1[section == "b" & set %in% SET_LAB, .(set, metric, n1 = n, v1 = value)], by = c("set", "metric"))
want_b <- data.table(set = c("2/3 set", "1/3 set", "3/3 set", "LFMM2", "RDA", "EMMAX"),   # [closeout] label, Phase 5 values (2/3, 3/3)
                     med = c(-0.030, -0.035, -0.045, -0.036, -0.026, -0.039),
                     ge0 = c(30.8, 26.3, 30.9, 25.8, 37.7, 25.5))
xw <- merge(DB[, .(set = set_lab(panel), median, ge0_pct)], want_b, by = "set")
chk3 <- c(table_s1_b_rows = nrow(xb) == 4L * N_PANEL, table_s1_b_values = max(abs(xb$value - xb$v1)) < 1e-9,
          table_s1_b_n = all(xb$n == xb$n1), stated = nrow(xw) == N_PANEL &&
              all(abs(xw$median - xw$med) < 0.0005 & abs(xw$ge0_pct - xw$ge0) < 0.05))
message(sprintf("t_offset_accuracy vs Table_S1.tsv (b): %d values, max |difference| %.2e", nrow(xb), max(abs(xb$value - xb$v1))))
print(chk3); if (!all(chk3)) warning("[LD-lab] t_offset_accuracy assertion failed: ", paste(names(chk3)[!chk3], collapse = ", "))
message("t_offset_accuracy asserts passed: ", paste(names(chk3), collapse = ", "))
# [audit3b] Table S3: each set's own offset accuracy and the causal loci's (`oracle`), from the rows already in OFF
d_id <- OFF[, max(abs(delta - (accuracy - oracle)))]
stopifnot(!anyNA(OFF$accuracy), !anyNA(OFF$oracle), d_id < 1e-12)          # delta = accuracy - oracle, every row
message(sprintf("[audit3b] OFF: %d rows; max |delta - (accuracy - oracle)| = %.2e", nrow(OFF), d_id))
ACC <- OFF[, .(median_accuracy = median(accuracy), min_accuracy = min(accuracy), max_accuracy = max(accuracy),
               n_acc = .N), by = .(panel = as.character(label))]
stopifnot(nrow(ACC) == N_PANEL, setequal(ACC$panel, DB$panel), all(ACC[match(DB$panel, panel), n_acc] == DB$n))
ORC <- OFF[, .(oracle = oracle[1L], n_values = uniqueN(oracle)), by = seed]
stopifnot(all(ORC$n_values == 1L), nrow(ORC) == NREP, setequal(ORC$seed, DET$seed))   # one causal-loci value per replicate
CAUSAL <- ORC[, .(n = .N, median_accuracy = median(oracle), min_accuracy = min(oracle), max_accuracy = max(oracle))]
for (i in seq_len(nrow(ACC)))
    message(sprintf("[audit3b] %s: median accuracy %.4f (lowest %.4f, highest %.4f), n = %d", set_lab(ACC$panel[i]),
                    ACC$median_accuracy[i], ACC$min_accuracy[i], ACC$max_accuracy[i], ACC$n_acc[i]))
message(sprintf("[audit3b] causal loci: median accuracy %.4f (lowest %.4f, highest %.4f), n = %d replicates",
                CAUSAL$median_accuracy, CAUSAL$min_accuracy, CAUSAL$max_accuracy, CAUSAL$n))

# ==================================================================== t_rank_tests (SIM:879-888)
M3 <- merge(RL[engine == "median of 3", .(seed, panel, rank)], COV, by = "seed")
stopifnot(nrow(M3) == NREP * N_PANEL)
KWT <- rbindlist(lapply(c("design cell", "demographic scenario"), function(f) rbindlist(lapply(PANEL_ORDER, function(s) {
    D <- M3[panel == s]
    g <- if (f == "design cell") factor(cell_of(D$arch, D$pleio, D$regime)) else factor(D$block)
    kt <- kruskal.test(D$rank, g)
    data.table(factor = f, panel = s, groups = nlevels(g), n = nrow(D), H = unname(kt$statistic),
               df = unname(kt$parameter), p = kt$p.value)
}))))
KWT[, p_holm := p.adjust(p, method = "holm")]
stopifnot(nrow(KWT) == 12L)
xe <- merge(KWT[, .(set = set_lab(panel), stratum = factor, H, df, p, p_holm, n)],
            TS1[section == "e", .(set, stratum, H1 = value, df1 = df, p1 = p, ph1 = p_holm, n1 = n)], by = c("set", "stratum"))
rel <- function(a, b) abs(a - b) / pmax(abs(b), .Machine$double.xmin)
NG <- NUM[startsWith(figure, "Figure S8 / Table S1 (e)") & metric == "groups", .(panel = rule, factor = stratum, g1 = value)]
xg <- merge(KWT[, .(panel, factor, groups)], NG, by = c("panel", "factor"))
e23 <- function(f) KWT[panel == "2/3 methods" & factor == f]
chk4 <- c(table_s1_e_rows = nrow(xe) == 12L,
          table_s1_e_values = all(rel(xe$H, xe$H1) < 1e-9 & xe$df == xe$df1 & rel(xe$p, xe$p1) < 1e-9 &
                                  rel(xe$p_holm, xe$ph1) < 1e-9 & xe$n == xe$n1),
          groups = nrow(xg) == 12L && all(xg$groups == xg$g1),
          ex_23_cells = with(e23("design cell"), abs(H - 6.75) < 0.005 && df == 11 && abs(p - 0.819) < 0.0005),   # [closeout] Phase 5 values (t_rank_tests.tsv)
          ex_23_scen  = with(e23("demographic scenario"), abs(H - 29.10) < 0.005 && df == 4 &&   # [closeout] Phase 5 values
                             abs(p - 7.44e-06) < 0.05e-06 && abs(p_holm - 2.23e-05) < 0.05e-05))
message(sprintf("t_rank_tests vs Table_S1.tsv (e): 12 tests, max relative |difference| H %.2e, p %.2e, Holm p %.2e",
                max(rel(xe$H, xe$H1)), max(rel(xe$p, xe$p1)), max(rel(xe$p_holm, xe$ph1))))
print(chk4); if (!all(chk4)) warning("[LD-lab] t_rank_tests assertion failed: ", paste(names(chk4)[!chk4], collapse = ", "))
message("t_rank_tests asserts passed: ", paste(names(chk4), collapse = ", "))

# ==================================================================== [closeout] Table S4's two tests
# (1) Per set: the 96 two-sided permutation p-values of Table_S1.tsv section a, as mvp_ms_sim_figures.R wrote them (sets
#     x rank 6 / ranks 5-6 / rank 1 / ranks 1-2 x rankings; 10,000 permutations of the set labels within each replicate,
#     seed 20261005, P = (k + 1) / 10,001), Holm-adjusted within each ranking's 24 tests. Not re-derived.
A2S <- setNames(names(A_MET), unname(A_MET))
PS <- TS1[section == "a", .(engine_lab = engine, set, stat = unname(A2S[metric]), p)]
stopifnot(nrow(PS) == 96L, !anyNA(PS$stat), !anyNA(PS$p), all(PS$p > 0 & PS$p <= 1), !anyDuplicated(PS[, .(engine_lab, set, stat)]),
          setequal(PS$engine_lab, unname(ENG_LAB)), PS[, .N, by = engine_lab][, all(N == 24L)])
PS[, p_holm := p.adjust(p, method = "holm"), by = engine_lab]
TR[, `:=`(engine_lab = unname(ENG_LAB[engine]), set_l = set_lab(panel))]
TR <- merge(TR, PS[, .(engine_lab, set_l = set, stat, p_perm = p, p_holm)], by = c("engine_lab", "set_l", "stat"),
            all.x = TRUE, sort = FALSE)
stopifnot(nrow(TR) == length(ENGS) * length(STATS) * N_PANEL, TR[!is.na(p_holm), .N] == 96L,
          TR[stat %in% names(A_MET), !anyNA(p_holm)], TR[!stat %in% names(A_MET), all(is.na(p_holm))])
message(sprintf("[closeout] per-set Holm: %d tests, %d with Holm p < 0.05 (%d raw)", nrow(PS), PS[p_holm < 0.05, .N], PS[p < 0.05, .N]))
# (2) Whole table, one test per ranking: the chi-square of its 6 x 6 table of set by rank (ties as the table counts them,
#     weight 1/2 at each tied rank, rank_split(); expected 100 per cell) against the same statistic after shuffling the
#     set labels within every replicate (each set's rank row moves whole, so ties are kept). 10,000 permutations; the
#     generator is set exactly as for the per-set tests (set.seed(20261005), then per permutation one runif() over the
#     replicate x set grid, replicates in radix seed order, sets in RULES order), so permutation b is the per-set tests'
#     permutation b, applied to all four rankings. p = (1 + k) / 10,001; Holm over the four rankings.
OMNI_B <- S1_B
E_CELL <- NREP / N_PANEL
VM <- lapply(setNames(ENGS, ENGS), function(e) {
    X <- RX[engine == e]
    V <- array(0, c(NREP, N_PANEL, N_PANEL))                                  # replicate x set x rank
    si <- match(X$seed, SEEDS); pj <- match(X$panel, RULES)
    stopifnot(!anyNA(si), !anyNA(pj), all(X$r %in% seq_len(N_PANEL)))
    for (j in seq_len(nrow(X))) V[si[j], pj[j], X$r[j]] <- V[si[j], pj[j], X$r[j]] + X$w[j]
    stopifnot(max(abs(apply(V, c(1, 2), sum) - 1)) < 1e-12, max(abs(apply(V, c(1, 3), sum) - 1)) < 1e-12)
    matrix(V, nrow = NREP * N_PANEL)                                          # row (s - 1) * NREP + i = V[i, s, ]
})
GRP <- rep(seq_len(N_PANEL), each = NREP)
chi_of <- function(O) sum((O - E_CELL)^2 / E_CELL)
OBS_T <- lapply(VM, function(M) rowsum(M, GRP))                               # set x rank counts
for (e in ENGS) {   # tie: the counts are t_ranks.tsv's rank 1..6 percentages x 6
    tr6 <- TR[engine == e & stat %in% sprintf("r%d", 1:6)]
    stopifnot(nrow(tr6) == 36L, max(abs(OBS_T[[e]][cbind(match(tr6$panel, RULES), as.integer(sub("r", "", tr6$stat)))] -
                                    tr6$value * NREP / 100)) < 1e-9)
}
CHI_OBS <- vapply(OBS_T, chi_of, 0)
# the per-set statistics on the same permuted labels (rank 6, ranks 5-6, rank 1, ranks 1-2), only to show that these are
# the per-set tests' permutations: their recomputed P must equal Table_S1.tsv's in all 96 cells (the P used is TS1's)
IND <- lapply(VM, function(M) cbind(r6 = M[, 6], r56 = M[, 5] + M[, 6], r1 = M[, 1], r12 = M[, 1] + M[, 2]))
OBS_I <- lapply(IND, function(M) 100 * rowsum(M, GRP) / NREP)                 # set x statistic, %
CH_I <- 100 * c(r6 = 1, r56 = 2, r1 = 1, r12 = 2) / N_PANEL
KK_I <- lapply(OBS_I, function(M) M * 0)
BLK <- rep(seq_len(NREP), each = N_PANEL)
PERM_CHI <- matrix(NA_real_, OMNI_B, length(ENGS), dimnames = list(NULL, ENGS))
set.seed(S1_SEED)
for (b in seq_len(OMNI_B)) {
    o <- order(BLK, runif(NREP * N_PANEL))
    P <- matrix((o - 1L) %% N_PANEL + 1L, nrow = NREP, byrow = TRUE)           # row i: a random order of the six sets
    idx <- as.vector((P - 1L) * NREP + seq_len(NREP))                          # set s of replicate i takes set P[i, s]'s row
    for (e in ENGS) {
        PERM_CHI[b, e] <- chi_of(rowsum(VM[[e]][idx, , drop = FALSE], GRP))
        pm <- 100 * rowsum(IND[[e]][idx, , drop = FALSE], GRP) / NREP
        KK_I[[e]] <- KK_I[[e]] + (abs(sweep(pm, 2, CH_I)) >= abs(sweep(OBS_I[[e]], 2, CH_I)) - 1e-9)
    }
}
OMNI <- data.table(engine = ENGS, chi2 = unname(CHI_OBS[ENGS]), permutations = OMNI_B,
                   k = vapply(ENGS, function(e) sum(PERM_CHI[, e] >= CHI_OBS[[e]] - 1e-9), 0L),
                   max_permuted_chi2 = apply(PERM_CHI, 2, max))
OMNI[, p := (1 + k) / (OMNI_B + 1)][, p_holm := p.adjust(p, method = "holm")][, seed := S1_SEED]
print(OMNI)
RECOMP <- rbindlist(lapply(ENGS, function(e) data.table(engine_lab = unname(ENG_LAB[e]), set = set_lab(RULES),
    stat = rep(colnames(KK_I[[e]]), each = N_PANEL), p_re = as.vector((KK_I[[e]] + 1) / (OMNI_B + 1)))))
RECOMP <- merge(RECOMP, PS[, .(engine_lab, set, stat, p)], by = c("engine_lab", "set", "stat"))
stopifnot(nrow(RECOMP) == 96L)
message(sprintf("[closeout] the omnibus permutations are the per-set tests': recomputed per-set P vs Table_S1.tsv, %d of 96 equal, max |difference| %.2e",
                RECOMP[abs(p_re - p) < 1e-12, .N], RECOMP[, max(abs(p_re - p))]))
fwrite(OMNI[, .(offset_method = unname(ENG_LAB[engine]), chi2, permutations, k, p, p_holm, seed)],
       file.path(out_dir, "t_ranks_omnibus.tsv"), sep = "\t")
fwrite(RECOMP, file.path(out_dir, "check_perset_p_recomputed.tsv"), sep = "\t")
fwrite(OMNI, file.path(out_dir, "check_omnibus_full.tsv"), sep = "\t")

# ==================================================================== TSVs (full precision, display labels)
ord_rows <- function(D) D[, `:=`(set_o = match(panel, PANEL_ORDER), stat_o = match(stat, names(STATS)))]
TR <- ord_rows(TR)[, eng_o := match(engine, ENGS)][order(eng_o, set_o, stat_o)]
fwrite(TR[, .(offset_method = unname(ENG_LAB[engine]), set = set_lab(panel), statistic = unname(STATS[stat]), value, lo, hi, n, p_holm)],   # [closeout] column engine -> offset_method; p_holm on the 96 tested rows
       file.path(out_dir, "t_ranks.tsv"), sep = "\t")
TRS[, str_o := match(paste(kind, stratum), paste(STRATA$kind, STRATA$stratum))]
TRS <- ord_rows(TRS)[order(str_o, set_o, stat_o)]
fwrite(TRS[, .(stratum_kind = kind, stratum = sub("centre-edge", "center-edge", stratum, fixed = TRUE), set = set_lab(panel), statistic = unname(STATS[stat]), value, lo, hi, n)],   # [closeout] center-edge
       file.path(out_dir, "t_ranks_strata.tsv"), sep = "\t")
DB <- DB[match(PANEL_ORDER, panel)]
# [audit3b] the set's own accuracy after n, then the differences; last row the causal loci, difference cells empty
DB <- merge(DB, ACC[, .(panel, median_accuracy, min_accuracy, max_accuracy)], by = "panel", sort = FALSE)[
    match(PANEL_ORDER, panel)]
stopifnot(identical(DB$panel, PANEL_ORDER), !anyNA(DB$median_accuracy))
TS3 <- rbind(DB[, .(set = set_lab(panel), n, median_accuracy, min_accuracy, max_accuracy, median_difference = median,
                    min_difference = min, max_difference = max, pct_replicates_difference_ge_0 = ge0_pct)],
             CAUSAL[, .(set = "Causal loci", n, median_accuracy, min_accuracy, max_accuracy,
                        median_difference = NA_real_, min_difference = NA_real_, max_difference = NA_real_,
                        pct_replicates_difference_ge_0 = NA_real_)])
stopifnot(nrow(TS3) == N_PANEL + 1L,
          identical(names(TS3), c("set", "n", "median_accuracy", "min_accuracy", "max_accuracy", "median_difference",
                                  "min_difference", "max_difference", "pct_replicates_difference_ge_0")))
fwrite(TS3, file.path(out_dir, "t_offset_accuracy.tsv"), sep = "\t")
KW_LAB <- c(`design cell` = "across the 12 design cells", `demographic scenario` = "across the 5 demographic scenarios")
fwrite(KWT[, .(test = unname(KW_LAB[factor]), set = set_lab(panel), H, df, p, p_holm, n, groups)],
       file.path(out_dir, "t_rank_tests.tsv"), sep = "\t")

# ==================================================================== markdown
minus <- function(s) sub("^-", "−", s)
f_pct <- function(x) minus(formatC(x, format = "f", digits = 1))
f_dif <- function(x) minus(formatC(x, format = "f", digits = 3))
f_ci  <- function(v, l, h) sprintf("%s [%s–%s]", f_pct(v), f_pct(l), f_pct(h))
SUP <- c(`-` = "⁻", `0` = "⁰", `1` = "¹", `2` = "²", `3` = "³", `4` = "⁴",
         `5` = "⁵", `6` = "⁶", `7` = "⁷", `8` = "⁸", `9` = "⁹")
f_p <- function(p) vapply(p, function(x) {
    if (x >= 0.001) return(formatC(signif(x, 3), format = "fg", digits = 3, flag = "#"))
    e <- floor(log10(x)); m <- x / 10^e
    if (round(m, 2) >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("%s × 10%s", formatC(m, format = "f", digits = 2), paste(SUP[strsplit(as.character(e), "")[[1]]], collapse = ""))
}, character(1))
md_row <- function(x) paste0("| ", paste(x, collapse = " | "), " |")
RANK_COLS <- c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Ranks 1–2", "Ranks 5–6")
rank_rows <- function(D, lead) {
    D <- copy(D)[, cell := f_ci(value, lo, hi)]
    W <- dcast(D, set_o + panel ~ stat_o, value.var = "cell")
    vapply(seq_len(nrow(W)), function(i) md_row(c(if (!is.null(lead)) lead(i), set_lab(W$panel[i]), unlist(W[i, -(1:2)]))), character(1))
}
# [closeout] Table S4's tests in t_ranks.md: Holm p after the CI in the four tested columns; one line per ranking for the
# whole-table test
f_ph <- function(p) trimws(formatC(signif(p, 2), format = "fg", digits = 2))   # formatC pads "1" to the digits width
rank_rows_p <- function(D) {
    D <- copy(D)[, cell := paste0(f_ci(value, lo, hi), fifelse(is.na(p_holm), "", paste0(", p = ", f_ph(p_holm))))]
    W <- dcast(D, set_o + panel ~ stat_o, value.var = "cell")
    vapply(seq_len(nrow(W)), function(i) md_row(c(set_lab(W$panel[i]), unlist(W[i, -(1:2)]))), character(1))
}
omni_line <- function(e) with(OMNI[engine == e], sprintf(
    "Whole table (6 × 6 counts of candidate set by rank): χ² = %s; permutation p = %s (k = %d of %s permutations), Holm-adjusted over the four rankings p = %s.",
    formatC(chi2, format = "f", digits = 1, big.mark = ","), f_ph(p), k, format(permutations, big.mark = ","), f_ph(p_holm)))
writeLines(c("# t_ranks", "",
             paste0("Cells: % of replicates [95% bootstrap CI]. In Rank 1, Rank 6, Ranks 1–2 and Ranks 5–6, p: the two-sided permutation test of the percentage against chance ",
                    "(16.7% for one rank, 33.3% for two; candidate-set labels shuffled within each replicate, 10,000 permutations, P = (k + 1) / 10,001), ",
                    "Holm-adjusted over the 24 tests of each ranking. Last line of each ranking: the χ² of its table of candidate set by rank (tied sets counted one half at each tied rank; ",
                    "expected 100 per cell) against the same permutations of the set labels; values in `t_ranks_omnibus.tsv`. ",
                    "Values: `t_ranks.tsv`; source: `Data_S2.zip` (`rank_long.tsv`); built by `benchmarks/mvp_ms_supp_rank_tables.R`."), "",
             unlist(lapply(ENGS, function(e) c(sprintf("## %s", ENG_LAB[[e]]), "",
                                               md_row(c("Candidate set", RANK_COLS)), md_row(rep("---", 9)),
                                               rank_rows_p(TR[engine == e]), "", omni_line(e), "")))),
           file.path(out_dir, "t_ranks.md"))
KIND_TITLE <- c(`genetic architecture` = "Genetic architecture", `demographic scenario` = "Demographic scenario",
                `design cell` = "Design cell")
writeLines(c("# t_ranks_strata", "",
             "Median of the three offset methods. Cells: % of the stratum's replicates [95% bootstrap CI]; n = replicates in the stratum. Values: `t_ranks_strata.tsv`; source: `Data_S2.zip` (`rank_long.tsv`, `detection_per_seed.tsv`); built by `benchmarks/mvp_ms_supp_rank_tables.R`.", "",
             unlist(lapply(names(KIND_TITLE), function(k) {
                 D <- TRS[kind == k]
                 c(sprintf("## %s", KIND_TITLE[[k]]), "", md_row(c(KIND_TITLE[[k]], "Candidate set", RANK_COLS)), md_row(rep("---", 10)),
                   unlist(lapply(STRATA[kind == k, stratum], function(s) {
                       Ds <- D[stratum == s]
                       rank_rows(Ds, function(i) if (i == 1L) sprintf("%s (n = %d)", sub("centre-edge", "center-edge", s, fixed = TRUE), Ds$n[1]) else "")   # [closeout] center-edge
                   })), "")
             }))),
           file.path(out_dir, "t_ranks_strata.md"))
writeLines(c("# t_offset_accuracy", "",
             "Accuracy = offset accuracy of the set, or of the causal loci in the last row, each the median of the three offset methods; difference = the set's accuracy − the causal loci's accuracy in the same replicate. Values: `t_offset_accuracy.tsv`; source: `Data_S2.zip` (`offset_delta_per_seed.tsv`); built by `benchmarks/mvp_ms_supp_rank_tables.R`.", "",
             md_row(c("Candidate set", "Replicates", "Median accuracy", "Lowest accuracy", "Highest accuracy",
                      "Median difference", "Minimum", "Maximum", "Replicates with difference ≥ 0, %")),
             md_row(rep("---", 9)),
             vapply(seq_len(nrow(DB)), function(i) with(DB[i], md_row(c(set_lab(panel), n, f_dif(median_accuracy),
                    f_dif(min_accuracy), f_dif(max_accuracy), f_dif(median), f_dif(min), f_dif(max), f_pct(ge0_pct)))),
                    character(1)),
             with(CAUSAL, md_row(c("Causal loci", n, f_dif(median_accuracy), f_dif(min_accuracy), f_dif(max_accuracy),
                                   "—", "—", "—", "—")))),   # [audit3b] causal loci row
           file.path(out_dir, "t_offset_accuracy.md"))
writeLines(c("# t_rank_tests", "",
             "Kruskal–Wallis tests; Holm p over the 12 tests. Values: `t_rank_tests.tsv`; source: `Data_S2.zip` (`rank_long.tsv`, `detection_per_seed.tsv`); built by `benchmarks/mvp_ms_supp_rank_tables.R`.", "",
             md_row(c("Test", "Candidate set", "H", "df", "p", "Holm p", "n", "Groups")), md_row(rep("---", 8)),
             vapply(seq_len(nrow(KWT)), function(i) with(KWT[i], md_row(c(KW_LAB[[factor]], set_lab(panel),
                    formatC(H, format = "f", digits = 2), df, f_p(p), f_p(p_holm), n, groups))), character(1))),
           file.path(out_dir, "t_rank_tests.md"))
message("wrote t_ranks, t_ranks_strata, t_offset_accuracy, t_rank_tests {tsv,md} to ", out_dir)
