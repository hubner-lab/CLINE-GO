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
             "LFMM" = "LFMM", "RDA" = "RDA", "EMMAX" = "EMMAX")
set_lab <- function(x) unname(SET_LAB[as.character(x)])
N_PANEL <- length(RULES)
NREP <- 600L
S1_SEED <- 20261005L
S1_B <- 10000L
cell_of <- function(a, p, r) paste(a, p, r, sep = ", ")                      # SIM:862

# ------------------------------------------------------------------ inputs (Data S2 members)
read_member <- function(member, ...) fread(cmd = paste("unzip -p", shQuote(zip), shQuote(member)), sep = "\t", ...)
RL  <- read_member("figures_ssclines_offset13/rank/rank_long.tsv", colClasses = c(seed = "character"))
DET <- read_member("figures_ssclines_gallery13/detection_per_seed.tsv", colClasses = c(seed = "character"))[label %in% RULES]
OFF <- read_member("figures_ssclines_gallery13/offset_delta_per_seed.tsv", colClasses = c(seed = "character"))[label %in% RULES]
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
stopifnot(identical(PANEL_ORDER, c("2/3 methods", "1/3 methods", "LFMM", "3/3 methods", "RDA", "EMMAX")))   # SIM:332
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
chk["ex_23_least"]  <- near3(ex("2/3 methods", "r6"), 1.83, 0.83, 3.00)
chk["ex_emmax_least"] <- near3(ex("EMMAX", "r6"), 32.33, 28.67, 36.17)
chk["ex_23_last2"] <- near3(ex("2/3 methods", "r56"), 11.08, 8.58, 13.67)
if (!all(chk)) stop("t_ranks assertion failed: ", paste(names(chk)[!chk], collapse = ", "))
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
                   want_lo = c(30, 5, 6, 10), want_hi = c(50, 14, 10, 18))
RNG <- merge(RNG, WANT, by = c("regime", "panel"))
RNG[, ok := half_up(lo, 0) == want_lo & half_up(hi, 0) == want_hi]
for (i in seq_len(nrow(RNG)))
    message(sprintf("two least accurate by design cell, %s, %s: %.4f-%.4f %% (stated %g-%g), n = 50 per cell",
                    RNG$regime[i], set_lab(RNG$panel[i]), RNG$lo[i], RNG$hi[i], RNG$want_lo[i], RNG$want_hi[i]))
chk2["design_cell_ranges"] <- nrow(RNG) == 4L && all(RNG$ok)
if (!all(chk2)) stop("t_ranks_strata assertion failed: ", paste(names(chk2)[!chk2], collapse = ", "))
message("t_ranks_strata asserts passed: ", paste(names(chk2), collapse = ", "))

# ==================================================================== t_offset_accuracy
DB <- OFF[, .(n = .N, median = median(delta), min = min(delta), max = max(delta), ge0_pct = 100 * mean(delta >= 0)),
          by = .(panel = as.character(label))]
B_MET <- c(median = "median delta accuracy", min = "min delta accuracy", max = "max delta accuracy",
           ge0_pct = "% replicates delta >= 0")
xb <- merge(melt(DB, id.vars = c("panel", "n"), variable.name = "stat", variable.factor = FALSE)[
                , .(set = set_lab(panel), metric = unname(B_MET[stat]), n, value)],
            TS1[section == "b" & set %in% SET_LAB, .(set, metric, n1 = n, v1 = value)], by = c("set", "metric"))
want_b <- data.table(set = c("2/3 set", "1/3 set", "3/3 set", "LFMM", "RDA", "EMMAX"),
                     med = c(-0.032, -0.035, -0.040, -0.036, -0.026, -0.039),
                     ge0 = c(27.8, 26.3, 31.1, 25.8, 37.7, 25.5))
xw <- merge(DB[, .(set = set_lab(panel), median, ge0_pct)], want_b, by = "set")
chk3 <- c(table_s1_b_rows = nrow(xb) == 4L * N_PANEL, table_s1_b_values = max(abs(xb$value - xb$v1)) < 1e-9,
          table_s1_b_n = all(xb$n == xb$n1), stated = nrow(xw) == N_PANEL &&
              all(abs(xw$median - xw$med) < 0.0005 & abs(xw$ge0_pct - xw$ge0) < 0.05))
message(sprintf("t_offset_accuracy vs Table_S1.tsv (b): %d values, max |difference| %.2e", nrow(xb), max(abs(xb$value - xb$v1))))
if (!all(chk3)) stop("t_offset_accuracy assertion failed: ", paste(names(chk3)[!chk3], collapse = ", "))
message("t_offset_accuracy asserts passed: ", paste(names(chk3), collapse = ", "))

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
          ex_23_cells = with(e23("design cell"), abs(H - 9.82) < 0.005 && df == 11 && abs(p - 0.546) < 0.0005),
          ex_23_scen  = with(e23("demographic scenario"), abs(H - 30.71) < 0.005 && df == 4 &&
                             abs(p - 3.5e-06) < 0.05e-06 && abs(p_holm - 1.4e-05) < 0.05e-05))
message(sprintf("t_rank_tests vs Table_S1.tsv (e): 12 tests, max relative |difference| H %.2e, p %.2e, Holm p %.2e",
                max(rel(xe$H, xe$H1)), max(rel(xe$p, xe$p1)), max(rel(xe$p_holm, xe$ph1))))
if (!all(chk4)) stop("t_rank_tests assertion failed: ", paste(names(chk4)[!chk4], collapse = ", "))
message("t_rank_tests asserts passed: ", paste(names(chk4), collapse = ", "))

# ==================================================================== TSVs (full precision, display labels)
ord_rows <- function(D) D[, `:=`(set_o = match(panel, PANEL_ORDER), stat_o = match(stat, names(STATS)))]
TR <- ord_rows(TR)[, eng_o := match(engine, ENGS)][order(eng_o, set_o, stat_o)]
fwrite(TR[, .(engine = unname(ENG_LAB[engine]), set = set_lab(panel), statistic = unname(STATS[stat]), value, lo, hi, n)],
       file.path(out_dir, "t_ranks.tsv"), sep = "\t")
TRS[, str_o := match(paste(kind, stratum), paste(STRATA$kind, STRATA$stratum))]
TRS <- ord_rows(TRS)[order(str_o, set_o, stat_o)]
fwrite(TRS[, .(stratum_kind = kind, stratum, set = set_lab(panel), statistic = unname(STATS[stat]), value, lo, hi, n)],
       file.path(out_dir, "t_ranks_strata.tsv"), sep = "\t")
DB <- DB[match(PANEL_ORDER, panel)]
fwrite(DB[, .(set = set_lab(panel), n, median_difference = median, min_difference = min, max_difference = max,
              pct_replicates_difference_ge_0 = ge0_pct)], file.path(out_dir, "t_offset_accuracy.tsv"), sep = "\t")
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
writeLines(c("# t_ranks", "",
             "Cells: % of replicates [95% bootstrap CI]. Values: `t_ranks.tsv`; source: `Data_S2.zip` (`rank_long.tsv`); built by `benchmarks/mvp_ms_supp_rank_tables.R`.", "",
             unlist(lapply(ENGS, function(e) c(sprintf("## %s", ENG_LAB[[e]]), "",
                                               md_row(c("Candidate set", RANK_COLS)), md_row(rep("---", 9)),
                                               rank_rows(TR[engine == e], NULL), "")))),
           file.path(out_dir, "t_ranks.md"))
KIND_TITLE <- c(`genetic architecture` = "Genetic architecture", `demographic scenario` = "Demographic scenario",
                `design cell` = "Design cell")
writeLines(c("# t_ranks_strata", "",
             "Median of the three offset engines. Cells: % of the stratum's replicates [95% bootstrap CI]; n = replicates in the stratum. Values: `t_ranks_strata.tsv`; source: `Data_S2.zip` (`rank_long.tsv`, `detection_per_seed.tsv`); built by `benchmarks/mvp_ms_supp_rank_tables.R`.", "",
             unlist(lapply(names(KIND_TITLE), function(k) {
                 D <- TRS[kind == k]
                 c(sprintf("## %s", KIND_TITLE[[k]]), "", md_row(c(KIND_TITLE[[k]], "Candidate set", RANK_COLS)), md_row(rep("---", 10)),
                   unlist(lapply(STRATA[kind == k, stratum], function(s) {
                       Ds <- D[stratum == s]
                       rank_rows(Ds, function(i) if (i == 1L) sprintf("%s (n = %d)", s, Ds$n[1]) else "")
                   })), "")
             }))),
           file.path(out_dir, "t_ranks_strata.md"))
writeLines(c("# t_offset_accuracy", "",
             "Difference = offset accuracy of the set − offset accuracy of the causal loci of the same replicate, median of the three offset engines. Values: `t_offset_accuracy.tsv`; source: `Data_S2.zip` (`offset_delta_per_seed.tsv`); built by `benchmarks/mvp_ms_supp_rank_tables.R`.", "",
             md_row(c("Candidate set", "Replicates", "Median difference", "Minimum", "Maximum", "Replicates with difference ≥ 0, %")),
             md_row(rep("---", 6)),
             vapply(seq_len(nrow(DB)), function(i) with(DB[i], md_row(c(set_lab(panel), n, f_dif(median), f_dif(min), f_dif(max), f_pct(ge0_pct)))),
                    character(1))),
           file.path(out_dir, "t_offset_accuracy.md"))
writeLines(c("# t_rank_tests", "",
             "Kruskal–Wallis tests; Holm p over the 12 tests. Values: `t_rank_tests.tsv`; source: `Data_S2.zip` (`rank_long.tsv`, `detection_per_seed.tsv`); built by `benchmarks/mvp_ms_supp_rank_tables.R`.", "",
             md_row(c("Test", "Candidate set", "H", "df", "p", "Holm p", "n", "Groups")), md_row(rep("---", 8)),
             vapply(seq_len(nrow(KWT)), function(i) with(KWT[i], md_row(c(KW_LAB[[factor]], set_lab(panel),
                    formatC(H, format = "f", digits = 2), df, f_p(p), f_p(p_holm), n, groups))), character(1))),
           file.path(out_dir, "t_rank_tests.md"))
message("wrote t_ranks, t_ranks_strata, t_offset_accuracy, t_rank_tests {tsv,md} to ", out_dir)
