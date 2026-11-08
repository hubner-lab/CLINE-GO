#!/usr/bin/env Rscript
# =============================================================================
# mvp_land_tables.R -- the simulation tables of the manuscript, per landscape and pooled (Job B of the
# landscapes pack, plan let-s-run-the-lanpscape-dynamic-pony.md).
#
# The published SS-Clines tables came from the closeout copies (work/closeout_lab/scratch/new/
# supp_rank_tables_s3.R, table_s1_detection_ld.R), not from the committed mvp_ms_*.R, which still hold the
# pre-LD-window numbers. This script extends the closeout logic; each computation below names its source.
#
# Every table has a `landscape` column: SS-Clines, SS-Mtn, Est-Clines, all (the 1,800 replicates).
#   t_detection, t_detection_arch   median (min-max) per candidate set, as published
#   t_replicates                    as published + K, its floor of 2, LD window, declared 619 bp, max VIF
#   t_offset_accuracy               per set and offset method: accuracy, difference from the causal loci,
#                                   accuracy lost against the replicate's best set (median, 90th percentile)
#   t_ranks, t_ranks_strata         % of replicates at each rank + ranks 1-2, 1-3, 4-6, 5-6, bootstrap CI
#   t_ranks_omnibus                 the chi-square of the set x rank table against permutations of the set labels
#                                   within each replicate (the only test of the ranks: the per-set tests of each
#                                   place were dropped 2026-11-08, the user; the bootstrap CIs carry each place)
#   t_rank_tests                    Kruskal-Wallis of each set's rank across design cells / scenarios
#   t_ranks_conditions              landscape x selection (equal, unequal): rank profile, loss
#   t_set_pairs                     3/3 - 2/3, 3/3 - RDA, 2/3 - RDA: ranks 4-6 and accuracy
#   fig_s1_numbers, fig_s3_numbers  the numbers of Figures S1 and S3
#   numbers.tsv                     the columns of figures_ld_ms/numbers.tsv + landscape + set
#   gate_ssclines.tsv               SS-Clines rows against the published tables (REF_DIR), cell by cell
#
# RULES OF THE PUBLISHED TABLES, kept: ties split evenly (mvp_rank_ties.R: 1/m at each tied rank), a set of
# fewer than 3 SNPs ranked last, percentile bootstrap over replicates with 10,000 resamples, 10,000
# permutations of the set labels within each replicate.
# STREAMS: bootstrap and permutation draw only on n and the seed, so extra statistics leave them alone.
#   SS-Clines keeps 20261005 and the published draws (one sample.int(n, n, TRUE) per resample; one
#   runif(n * 6) per permutation, replicates in radix order, sets in RULES order). SS-Mtn 20261106,
#   Est-Clines 20261107, all 20261108, each set before every stratum as published. t_ranks_conditions and
#   t_set_pairs add an offset to the scope's seed. The bootstrap of `all` is plain over its 1,800 replicates,
#   not stratified by landscape, the rule of the published tables (which do not stratify by design cell).
# HOLM: t_ranks_omnibus p_holm = over the four rankings of one scope. The per-set p and their Holm family of 30
#   are gone (2026-11-08); the permutation draws are unchanged, so the omnibus k and p reproduce the 2026-11-07 run.
# ORDER: rows and figures use one order, fewest replicates in ranks 4-6 over the 1,800, median of the three
#   offset methods; computed here and asserted against the laptop's. Permutation work keeps RULES order.
#
# Usage: OUT_DIR=<dir> REF_DIR=<published t_*.tsv> Rscript benchmarks/mvp_land_tables.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OUT  <- Sys.getenv("OUT_DIR", file.path(EVAL, "landscape_ms/tables"))
REF  <- Sys.getenv("REF_DIR", "/reference")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "benchmarks/mvp_rank_ties.R"))

# ------------------------------------------------------------------ vocabulary
LANDS <- data.table(
    land  = c("SS-Clines", "SS-Mtn", "Est-Clines"),
    lkey  = c("ssclines", "ssmtn", "estclines"),
    gal   = c("figures_ld_gallery14", "figures_land_ssmtn_gallery", "figures_land_estclines_gallery"),
    rank  = c("figures_ld_offset14/rank", "figures_land_ssmtn_offset15_ssmtn/rank",
              "figures_land_estclines_offset15_estclines/rank"),
    rem   = c("remeasure600", "remeasure600_ssmtn", "remeasure600_estclines"),
    seeds = c("work/ld_phase3_2026-10-26/scripts/seeds600.txt", "work/mvp_land/seeds_ssmtn.txt",
              "work/mvp_land/seeds_estclines.txt"))
SCOPES  <- c(LANDS$land, "all")
SEED_OF <- c("SS-Clines" = 20261005L, "SS-Mtn" = 20261106L, "Est-Clines" = 20261107L, all = 20261108L)
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
BLOCK_LAB <- c("N variable, m variable", "N cline north-south", "N equal, m constant",
               "N cline centre-edge", "N equal, m breaks")
block_show <- function(x) sub("centre-edge", "center-edge", x, fixed = TRUE)
RULES <- c("2/3 methods", "1/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")
ENG_LAB <- c(gradientForest = "Gradient Forest", LFMM2 = "Geometric offset", RDA = "RDA-based offset",
             `median of 3` = "Median of the three")
ENGS <- names(ENG_LAB)
SET_LAB <- c("2/3 methods" = "2/3 set", "1/3 methods" = "1/3 set", "3/3 methods" = "3/3 set",
             "LFMM" = "LFMM2", "RDA" = "RDA", "EMMAX" = "EMMAX")
set_lab <- function(x) unname(SET_LAB[as.character(x)])
N_PANEL <- length(RULES)
S1_B <- as.integer(Sys.getenv("LAND_B", "10000"))   # resamples and permutations; a smaller LAND_B is for a dry run only
cell_of <- function(a, p, r) paste(a, p, r, sep = ", ")
in_scope <- function(D, s) if (s == "all") copy(D) else D[landscape == s]
SEL_LAB <- c("equal-S" = "equal selection", "unequal-S" = "unequal selection")

# ------------------------------------------------------------------ inputs
rd <- function(...) fread(file.path(...), colClasses = c(seed = "character"))
man <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
LDW <- rd(ROOT, "work/mvp_land/ld_window_by_replicate.tsv")
fac <- function(D) {
    D[, arch := factor(arch, levels = ARCH_LABELS)]
    if ("regime" %in% names(D)) D[, regime := factor(regime, levels = c("equal-S", "unequal-S"))]
    if ("pleio"  %in% names(D)) D[, pleio  := factor(pleio,  levels = c("no pleiotropy", "pleiotropy"))]
    if ("block"  %in% names(D)) D[, block  := factor(block,  levels = BLOCK_LAB)]
    D[, label := factor(label, levels = RULES)]
    stopifnot(!anyNA(D$arch), !anyNA(D$label))
    D[]
}
SEEDS_OF <- list(); RL <- list(); DET <- list(); OFF <- list(); TESTED <- list()
for (i in seq_len(nrow(LANDS))) {
    L <- LANDS$land[i]
    s <- trimws(readLines(file.path(ROOT, LANDS$seeds[i]))); s <- s[nzchar(s)]
    stopifnot(length(s) == 600L, !anyDuplicated(s))
    SEEDS_OF[[L]] <- s
    rl <- rd(EVAL, LANDS$rank[i], "rank_long.tsv")
    stopifnot(setequal(rl$seed, s), setequal(unique(rl$engine), ENGS), setequal(unique(rl$panel), RULES),
              nrow(rl) == 600L * length(ENGS) * N_PANEL, !anyDuplicated(rl[, .(seed, engine, panel)]))
    RL[[L]] <- rl[, landscape := L]
    d <- fac(rd(EVAL, LANDS$gal[i], "detection_per_seed.tsv")[label %in% RULES])
    stopifnot(setequal(d$seed, s), nrow(d) == 600L * N_PANEL, !anyNA(d$block), !anyNA(d$regime), !anyNA(d$pleio))
    DET[[L]] <- d[, landscape := L]
    o <- fac(rd(EVAL, LANDS$gal[i], "offset_delta_per_seed.tsv")[label %in% RULES])
    stopifnot(all(o$seed %in% s))
    OFF[[L]] <- o[, landscape := L]
    TESTED[[L]] <- rd(EVAL, LANDS$rem[i], "rdaunc_rda2x_ldb/rank_metrics.tsv")[method == "LFMM", .(seed, n_snps)]
    stopifnot(setequal(TESTED[[L]]$seed, s))
}
RL <- rbindlist(RL); DET <- rbindlist(DET); OFF <- rbindlist(OFF); TESTED <- rbindlist(TESTED)
stopifnot(!anyDuplicated(unlist(SEEDS_OF)), uniqueN(RL$seed) == 1800L)
# tie groups found from the ranks are the ties of the raw accuracies (closeout, SIM:552, 562-564)
RL[, key := fifelse(empty, -Inf, accuracy)]
RL[, g_key := .GRP, by = .(seed, engine, key)][, g_rank := .GRP, by = .(seed, engine, rank)]
stopifnot(uniqueN(RL[, .(g_key, g_rank)]) == uniqueN(RL$g_key), uniqueN(RL$g_key) == uniqueN(RL$g_rank))
RL[, c("g_key", "g_rank") := NULL]
RX <- rank_split(RL, c("seed", "engine"))
stopifnot(RX[tie_n == 1L, all(r == rank & w == 1)])
message(sprintf("ties: %d tie groups of two sets, %d of three or more",
                uniqueN(RX[tie_n == 2L, .(seed, engine, rank)]), uniqueN(RX[tie_n > 2L, .(seed, engine, rank)])))

STATS <- c(r1 = "rank 1", r2 = "rank 2", r3 = "rank 3", r4 = "rank 4", r5 = "rank 5", r6 = "rank 6",
           r12 = "ranks 1–2", r123 = "ranks 1–3", r456 = "ranks 4–6", r56 = "ranks 5–6")
ST <- RX[, .(r1 = sum(w[r == 1L]), r2 = sum(w[r == 2L]), r3 = sum(w[r == 3L]), r4 = sum(w[r == 4L]),
             r5 = sum(w[r == 5L]), r6 = sum(w[r == 6L]), r12 = sum(w[r <= 2L]), r123 = sum(w[r <= 3L]),
             r456 = sum(w[r >= 4L]), r56 = sum(w[r >= 5L])), by = .(seed, engine, panel)]
stopifnot(nrow(ST) == 1800L * length(ENGS) * N_PANEL, max(abs(ST$r123 + ST$r456 - 1)) < 1e-12)
LAND_OF <- unique(RL[, .(seed, landscape)])
ST <- merge(ST, LAND_OF, by = "seed")

# the one order: fewest replicates in ranks 4-6 over the 1,800, median of the three
ORD <- ST[engine == "median of 3", .(pct = 100 * mean(r456)), by = panel][order(pct)]
PANEL_ORDER <- ORD$panel
message("set order (ranks 4-6, all 1,800, median of the three): ",
        paste(sprintf("%s %.2f", set_lab(ORD$panel), ORD$pct), collapse = " < "))
stopifnot(identical(PANEL_ORDER, c("RDA", "3/3 methods", "EMMAX", "2/3 methods", "1/3 methods", "LFMM")))
fwrite(ORD[, .(set = set_lab(panel), key = panel, pct_ranks_4_6 = pct, order = seq_len(.N))],
       file.path(OUT, "set_order.tsv"), sep = "\t")

COV <- unique(DET[, .(seed, landscape, arch, pleio, regime, block)])
stopifnot(nrow(COV) == 1800L, !anyDuplicated(COV$seed), !anyNA(COV))
# selection: the regime column of detection_per_seed.tsv must be the deposit's architecture string
SEL <- man[seed %in% COV$seed, .(seed, from_arch = fifelse(grepl("unequal-S$", architecture), "unequal-S",
                                                           fifelse(grepl("equal-S$", architecture), "equal-S", NA_character_)))]
SEL <- merge(SEL, COV[, .(seed, regime = as.character(regime))], by = "seed")
stopifnot(nrow(SEL) == 1800L, !anyNA(SEL$from_arch), all(SEL$from_arch == SEL$regime))
message("selection regime = suffix of the architecture string in all 1,800 replicates")
COV[, cell := cell_of(arch, pleio, regime)]
CELLS <- CJ(arch = ARCH_LABELS, pleio = c("no pleiotropy", "pleiotropy"), regime = c("equal-S", "unequal-S"),
            sorted = FALSE)[, cell_of(arch, pleio, regime)]

# ------------------------------------------------------------------ bootstrap (closeout build_big, boot_ci)
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
boot_ci <- function(M, n, seed) {
    set.seed(seed)
    B <- matrix(NA_real_, S1_B, ncol(M))
    for (b in seq_len(S1_B)) {
        i <- sample.int(n, n, replace = TRUE)
        B[b, ] <- colMeans(M[i, , drop = FALSE])
    }
    data.table(value = 100 * colMeans(M),
               lo = apply(100 * B, 2, quantile, 0.025, names = FALSE),
               hi = apply(100 * B, 2, quantile, 0.975, names = FALSE))
}

# ================================================================= t_ranks, omnibus, per-set p
TR <- list(); OMNI <- list()
for (sc in SCOPES) {
    S <- in_scope(ST, sc)
    SEEDS <- sort(unique(S$seed), method = "radix")
    NREP <- length(SEEDS)
    message(sprintf("=== %s: %d replicates, seed %d", sc, NREP, SEED_OF[[sc]]))
    BG <- build_big(S, ENGS, SEEDS)
    tr <- cbind(BG$idx, boot_ci(BG$M, NREP, SEED_OF[[sc]]))
    tr[, `:=`(landscape = sc, n = NREP)]
    TR[[sc]] <- tr
    # whole-table test against permutations of the set labels within each replicate (closeout: Table S4's test)
    X <- RX[seed %in% SEEDS]
    VM <- lapply(setNames(ENGS, ENGS), function(e) {
        Y <- X[engine == e]
        V <- array(0, c(NREP, N_PANEL, N_PANEL))                                  # replicate x set x rank
        si <- match(Y$seed, SEEDS); pj <- match(Y$panel, RULES)
        stopifnot(!anyNA(si), !anyNA(pj), all(Y$r %in% seq_len(N_PANEL)))
        for (j in seq_len(nrow(Y))) V[si[j], pj[j], Y$r[j]] <- V[si[j], pj[j], Y$r[j]] + Y$w[j]
        stopifnot(max(abs(apply(V, c(1, 2), sum) - 1)) < 1e-12, max(abs(apply(V, c(1, 3), sum) - 1)) < 1e-12)
        matrix(V, nrow = NREP * N_PANEL)                                          # row (s - 1) * NREP + i = V[i, s, ]
    })
    GRP <- rep(seq_len(N_PANEL), each = NREP)
    E_CELL <- NREP / N_PANEL
    chi_of <- function(O) sum((O - E_CELL)^2 / E_CELL)
    OBS_T <- lapply(VM, function(M) rowsum(M, GRP))                               # set x rank counts
    for (e in ENGS) {   # tie: the counts are the table's rank 1..6 percentages
        tr6 <- tr[engine == e & stat %in% sprintf("r%d", 1:6)]
        stopifnot(nrow(tr6) == 36L, max(abs(OBS_T[[e]][cbind(match(tr6$panel, RULES), as.integer(sub("r", "", tr6$stat)))] -
                                        tr6$value * NREP / 100)) < 1e-9)
    }
    CHI_OBS <- vapply(OBS_T, chi_of, 0)
    BLK <- rep(seq_len(NREP), each = N_PANEL)
    PERM_CHI <- matrix(NA_real_, S1_B, length(ENGS), dimnames = list(NULL, ENGS))
    set.seed(SEED_OF[[sc]])
    for (b in seq_len(S1_B)) {
        o <- order(BLK, runif(NREP * N_PANEL))
        P <- matrix((o - 1L) %% N_PANEL + 1L, nrow = NREP, byrow = TRUE)           # row i: a random order of the six sets
        idx <- as.vector((P - 1L) * NREP + seq_len(NREP))                          # set s of replicate i takes set P[i, s]'s row
        for (e in ENGS) {
            PERM_CHI[b, e] <- chi_of(rowsum(VM[[e]][idx, , drop = FALSE], GRP))
        }
    }
    om <- data.table(landscape = sc, engine = ENGS, chi2 = unname(CHI_OBS[ENGS]), permutations = S1_B,
                     k = vapply(ENGS, function(e) sum(PERM_CHI[, e] >= CHI_OBS[[e]] - 1e-9), 0L),
                     max_permuted_chi2 = apply(PERM_CHI, 2, max))
    om[, p := (1 + k) / (S1_B + 1)][, p_holm := p.adjust(p, method = "holm")][, seed := SEED_OF[[sc]]]
    OMNI[[sc]] <- om
    print(om[, .(engine, chi2, k, p, p_holm, max_permuted_chi2)])
}
TR <- rbindlist(TR); OMNI <- rbindlist(OMNI)
ord_rows <- function(D) D[, `:=`(set_o = match(panel, PANEL_ORDER), stat_o = match(stat, names(STATS)))]
TR <- ord_rows(TR)[, `:=`(eng_o = match(engine, ENGS), sc_o = match(landscape, SCOPES))][order(sc_o, eng_o, set_o, stat_o)]
fwrite(TR[, .(landscape, offset_method = unname(ENG_LAB[engine]), set = set_lab(panel), statistic = unname(STATS[stat]),
              value, lo, hi, n)], file.path(OUT, "t_ranks.tsv"), sep = "\t")
fwrite(OMNI[, .(landscape, offset_method = unname(ENG_LAB[engine]), chi2, permutations, k, p, p_holm, seed, max_permuted_chi2)],
       file.path(OUT, "t_ranks_omnibus.tsv"), sep = "\t")

# ================================================================= t_ranks_strata
S3 <- merge(ST[engine == "median of 3"], COV[, .(seed, arch, block, cell)], by = "seed")
STRATA <- rbind(data.table(kind = "genetic architecture", stratum = ARCH_LABELS, col = "arch"),
                data.table(kind = "demographic scenario", stratum = BLOCK_LAB, col = "block"),
                data.table(kind = "design cell", stratum = CELLS, col = "cell"))
TRS <- rbindlist(lapply(SCOPES, function(sc) {
    message("strata: ", sc)
    S <- in_scope(S3, sc)
    rbindlist(lapply(seq_len(nrow(STRATA)), function(k) {
        s <- STRATA[k]
        D <- S[as.character(get(s$col)) == s$stratum]
        sd <- sort(unique(D$seed), method = "radix")
        G <- build_big(D, "median of 3", sd)
        cbind(landscape = sc, kind = s$kind, stratum = s$stratum, G$idx[, .(stat, panel)],
              boot_ci(G$M, length(sd), SEED_OF[[sc]]), n = length(sd))
    }))
}))
stopifnot(nrow(TRS) == 4L * nrow(STRATA) * length(STATS) * N_PANEL)
TRS[, `:=`(str_o = match(paste(kind, stratum), paste(STRATA$kind, STRATA$stratum)), sc_o = match(landscape, SCOPES))]
TRS <- ord_rows(TRS)[order(sc_o, str_o, set_o, stat_o)]
fwrite(TRS[, .(landscape, stratum_kind = kind, stratum = block_show(stratum), set = set_lab(panel),
               statistic = unname(STATS[stat]), value, lo, hi, n)], file.path(OUT, "t_ranks_strata.tsv"), sep = "\t")

# ================================================================= t_offset_accuracy
# Median of the three: from offset_delta_per_seed.tsv, the published source; rank_long.tsv must hold the same values.
chk <- merge(RL[engine == "median of 3" & empty == FALSE, .(seed, label = panel, a = accuracy, o = oracle, d = delta)],
             OFF[, .(seed, label = as.character(label), accuracy, oracle, delta)], by = c("seed", "label"), all = TRUE)
stopifnot(!anyNA(chk), max(abs(chk$a - chk$accuracy), abs(chk$o - chk$oracle), abs(chk$d - chk$delta)) < 1e-12)
ACC_RL <- rbind(RL[engine != "median of 3", .(landscape, seed, engine, panel, empty, accuracy, oracle, delta)],
                OFF[, .(landscape, seed, engine = "median of 3", panel = as.character(label), empty = FALSE, accuracy, oracle, delta)])
# the rows of a set without an offset: no accuracy, no loss, n says so
ACC_RL <- ACC_RL[empty == FALSE & !is.na(accuracy)]
ACC_RL[, best := max(accuracy), by = .(seed, engine)][, loss := best - accuracy]
stopifnot(all(ACC_RL$loss >= 0))
TOA <- rbindlist(lapply(SCOPES, function(sc) {
    A <- in_scope(ACC_RL, sc)
    sets <- A[, .(n = .N, median_accuracy = median(accuracy), min_accuracy = min(accuracy), max_accuracy = max(accuracy),
                  median_difference = median(delta), min_difference = min(delta), max_difference = max(delta),
                  pct_replicates_difference_ge_0 = 100 * mean(delta >= 0),
                  median_loss = median(loss), p90_loss = quantile(loss, 0.9, names = FALSE)), by = .(engine, panel)]
    orc <- unique(A[, .(seed, engine, oracle)])
    stopifnot(!anyDuplicated(orc[, .(seed, engine)]))
    cau <- orc[, .(n = .N, median_accuracy = median(oracle), min_accuracy = min(oracle), max_accuracy = max(oracle)), by = engine]
    rbind(sets[order(match(engine, ENGS), match(panel, PANEL_ORDER))][, set := set_lab(panel)],
          cau[, `:=`(panel = "causal", set = "Causal loci")], fill = TRUE)[, landscape := sc]
}))
TOA[, offset_method := unname(ENG_LAB[engine])]
TOA <- TOA[order(match(landscape, SCOPES), match(engine, ENGS), match(panel, c(PANEL_ORDER, "causal")))]
fwrite(TOA[, .(landscape, offset_method, set, n, median_accuracy, min_accuracy, max_accuracy, median_difference, min_difference,
               max_difference, pct_replicates_difference_ge_0, median_loss, p90_loss)],
       file.path(OUT, "t_offset_accuracy.tsv"), sep = "\t")

# ================================================================= t_rank_tests (closeout KWT, Holm per scope)
M3 <- merge(RL[engine == "median of 3", .(seed, panel, rank)], COV, by = "seed")
KWT <- rbindlist(lapply(SCOPES, function(sc) {
    M <- in_scope(M3, sc)
    k <- rbindlist(lapply(c("design cell", "demographic scenario"), function(f) rbindlist(lapply(PANEL_ORDER, function(s) {
        D <- M[panel == s]
        g <- if (f == "design cell") factor(cell_of(D$arch, D$pleio, D$regime)) else factor(D$block)
        kt <- kruskal.test(D$rank, g)
        data.table(factor = f, panel = s, groups = nlevels(g), n = nrow(D), H = unname(kt$statistic),
                   df = unname(kt$parameter), p = kt$p.value)
    }))))
    k[, `:=`(p_holm = p.adjust(p, method = "holm"), landscape = sc)]
}))
stopifnot(nrow(KWT) == 48L)
KW_LAB <- c(`design cell` = "across the 12 design cells", `demographic scenario` = "across the 5 demographic scenarios")
fwrite(KWT[, .(landscape, test = unname(KW_LAB[factor]), set = set_lab(panel), H, df, p, p_holm, n, groups)],
       file.path(OUT, "t_rank_tests.tsv"), sep = "\t")

# ================================================================= t_ranks_conditions
LOSS_M3 <- merge(ACC_RL[engine == "median of 3", .(seed, panel, loss)], COV[, .(seed, regime)], by = "seed")
COND_STATS <- c("r1", "r2", "r3", "r4", "r5", "r6", "r456")
TC <- rbindlist(lapply(SCOPES, function(sc) rbindlist(lapply(seq_along(SEL_LAB), function(k) {
    reg <- names(SEL_LAB)[k]
    keep <- COV[regime == reg & (sc == "all" | landscape == sc), seed]
    sd <- sort(keep, method = "radix")
    G <- build_big(S3[seed %in% sd], "median of 3", sd)
    seed_k <- SEED_OF[[sc]] + 100L + k
    out <- cbind(G$idx[, .(stat, panel)], boot_ci(G$M, length(sd), seed_k), n = length(sd))[stat %in% COND_STATS]
    # the median loss, with its CI over the same resamples of the condition's replicates
    L <- dcast(LOSS_M3[seed %in% sd], seed ~ panel, value.var = "loss")
    L <- merge(data.table(seed = sd), L, by = "seed", all.x = TRUE)
    stopifnot(identical(L$seed, sd))
    LM <- as.matrix(L[, RULES, with = FALSE])
    set.seed(seed_k)
    B <- matrix(NA_real_, S1_B, N_PANEL)
    for (b in seq_len(S1_B)) {
        i <- sample.int(length(sd), length(sd), replace = TRUE)
        B[b, ] <- vapply(seq_len(N_PANEL), function(j) median(LM[i, j], na.rm = TRUE), 0)
    }
    loss <- data.table(stat = "loss", panel = RULES, value = apply(LM, 2, median, na.rm = TRUE),
                       lo = apply(B, 2, quantile, 0.025, names = FALSE), hi = apply(B, 2, quantile, 0.975, names = FALSE),
                       n = colSums(!is.na(LM)))
    rbind(out, loss)[, `:=`(landscape = sc, selection = unname(SEL_LAB[reg]), seed = seed_k)]
}))))
STAT_C <- c(STATS[COND_STATS], loss = "median accuracy lost against the best set")
TC[, `:=`(sc_o = match(landscape, SCOPES), sel_o = match(selection, SEL_LAB), set_o = match(panel, PANEL_ORDER), st_o = match(stat, names(STAT_C)))]
TC <- TC[order(sc_o, sel_o, set_o, st_o)]
fwrite(TC[, .(landscape, selection, set = set_lab(panel), statistic = unname(STAT_C[stat]), value, lo, hi, n, seed)],
       file.path(OUT, "t_ranks_conditions.tsv"), sep = "\t")

# ================================================================= t_set_pairs
PAIRS <- list(c("3/3 methods", "2/3 methods"), c("3/3 methods", "RDA"), c("2/3 methods", "RDA"))
pair_lab <- function(p) sprintf("%s minus %s", set_lab(p[1]), set_lab(p[2]))
q7 <- function(x, p) {                    # type-7 quantile of sorted x without names, as quantile()
    n <- length(x); h <- (n - 1) * p + 1; f <- floor(h); c <- min(f + 1, n)
    x[f] + (h - f) * (x[c] - x[f])
}
TP <- rbindlist(lapply(SCOPES, function(sc) {
    message("pairs: ", sc)
    S <- in_scope(ST, sc)[, .(seed, engine, panel, r456)]
    sd <- sort(unique(S$seed), method = "radix"); NREP <- length(sd)
    A <- in_scope(RL, sc)[, .(seed, engine, panel, accuracy = fifelse(empty, NA_real_, accuracy))]
    res <- rbindlist(lapply(seq_along(ENGS), function(ei) {
        e <- ENGS[ei]
        WR <- dcast(S[engine == e], seed ~ panel, value.var = "r456")
        WA <- dcast(A[engine == e], seed ~ panel, value.var = "accuracy")
        stopifnot(identical(WR$seed, sd), identical(WA$seed, sd))
        rbindlist(lapply(seq_along(PAIRS), function(k) {
            p <- PAIRS[[k]]
            dr <- 100 * (WR[[p[1]]] - WR[[p[2]]])
            da <- WA[[p[1]]] - WA[[p[2]]]
            ok <- which(!is.na(da))
            set.seed(SEED_OF[[sc]] + 1000L * ei + 10L * k)
            B <- matrix(NA_real_, S1_B, 5L)
            for (b in seq_len(S1_B)) {
                i <- sample.int(NREP, NREP, replace = TRUE)
                x <- da[i]; x <- x[!is.na(x)]
                sx <- sort(x)
                B[b, ] <- c(mean(dr[i]), q7(sx, 0.5), q7(sx, 0.1), q7(sx, 0.9), 100 * mean(x > 0))
            }
            sa <- sort(da[ok])
            val <- c(mean(dr), q7(sa, 0.5), q7(sa, 0.1), q7(sa, 0.9), 100 * mean(da[ok] > 0))
            data.table(landscape = sc, offset_method = unname(ENG_LAB[e]), pair = pair_lab(p),
                       measure = c("ranks 4–6, percentage points", "accuracy, median difference",
                                   "accuracy, 10th percentile of the difference", "accuracy, 90th percentile of the difference",
                                   "% replicates in which the first is more accurate"),
                       value = val, lo = apply(B, 2, quantile, 0.025, names = FALSE), hi = apply(B, 2, quantile, 0.975, names = FALSE),
                       n = c(NREP, rep(length(ok), 4L)), n_tied = c(NA, rep(sum(da[ok] == 0), 4L)),
                       seed = SEED_OF[[sc]] + 1000L * ei + 10L * k)
        }))
    }))
    res
}))
fwrite(TP, file.path(OUT, "t_set_pairs.tsv"), sep = "\t")

# ================================================================= t_detection, t_detection_arch, t_replicates
metrics_a <- c(n = "SNPs identified", TP = "True positives (causal or linked SNPs)",
               FP = "False positives (neutral SNPs)", prec_pct = "Precision, %",
               n_causal = "Causal loci identified", rec_pct = "Causal loci identified, % of the replicate's",
               n_linked = "Linked SNPs identified")
D2 <- copy(DET)
D2[, `:=`(prec_pct = 100 * prec_cl, rec_pct = 100 * rec_c)]
D2[, (names(metrics_a)) := lapply(.SD, as.numeric), .SDcols = names(metrics_a)]
D2[, set := set_lab(label)]
LONG <- melt(D2, id.vars = c("seed", "landscape", "set", "arch"), measure.vars = names(metrics_a),
             variable.name = "metric", variable.factor = FALSE)[!is.na(value)]
sum_s <- function(D, by) D[, .(median = median(value), min = min(value), max = max(value), n = .N), by = by]
SET_ORDER <- set_lab(PANEL_ORDER)
TD <- rbindlist(lapply(SCOPES, function(sc) sum_s(in_scope(LONG, sc), c("set", "metric"))[, landscape := sc]))
TD[, `:=`(sc_o = match(landscape, SCOPES), set_o = match(set, SET_ORDER), m_o = match(metric, names(metrics_a)),
          label = unname(metrics_a[metric]))]
TD <- TD[order(sc_o, set_o, m_o)]
fwrite(TD[, .(landscape, set, label, median, min, max, n)], file.path(OUT, "t_detection.tsv"), sep = "\t")
TDA <- rbindlist(lapply(SCOPES, function(sc) sum_s(in_scope(LONG, sc), c("arch", "set", "metric"))[, landscape := sc]))
TDA[, `:=`(sc_o = match(landscape, SCOPES), a_o = match(arch, ARCH_LABELS), set_o = match(set, SET_ORDER),
           m_o = match(metric, names(metrics_a)), label = unname(metrics_a[metric]))]
TDA <- TDA[order(sc_o, a_o, set_o, m_o)]
fwrite(TDA[, .(landscape, architecture = arch, set, label, median, min, max, n)], file.path(OUT, "t_detection_arch.tsv"), sep = "\t")

# the replicates
REPS <- unique(DET[, .(seed, landscape, arch, causal_total)])
stopifnot(nrow(REPS) == 1800L, !anyDuplicated(REPS$seed))
REPS <- merge(REPS, TESTED, by = "seed")
REPS[, causal_pct := 100 * causal_total / n_snps]
KT <- man[seed %in% REPS$seed, .(seed, k_best, K_authors)]
REPS <- merge(REPS, KT, by = "seed")
stopifnot(all(REPS$k_best >= 2L), all(REPS[K_authors >= 2L, k_best == K_authors]), all(REPS[K_authors < 2L, k_best == 2L]))
REPS <- merge(REPS, LDW[, .(seed, ld_bp = clumping_distance_bp, distance_source)], by = "seed")
vif_hdrc <- function(s) {
    f <- Sys.glob(file.path(ROOT, sprintf("MVP%shdrc_results/GEA/tables/methods/RDA/RDA_diagnostics_K*.tsv", s)))
    stopifnot(length(f) == 1L)
    x <- readLines(f, warn = FALSE)
    val <- function(k) sub(paste0("^", k, "\t"), "", grep(paste0("^", k, "\t"), x, value = TRUE))
    data.table(seed = s, vif = as.numeric(val("max_vif")), cond = as.integer(val("condition_pcs")))
}
VIF <- rbindlist(lapply(c(SEEDS_OF[["SS-Mtn"]], SEEDS_OF[["Est-Clines"]]), vif_hdrc))
VC <- rd(EVAL, "remeasure600/onefit_diag/rda_diag.tsv")
VC <- dcast(VC[key %in% c("max_vif", "condition_pcs")], seed ~ key, value.var = "value")
VIF <- rbind(VIF, VC[, .(seed, vif = as.numeric(max_vif), cond = as.integer(condition_pcs))])
REPS <- merge(REPS, VIF, by = "seed")
stopifnot(nrow(REPS) == 1800L, !anyNA(REPS$vif), all(REPS$cond == REPS$k_best))
summ <- function(x) list(median = median(x), min = min(x), max = max(x), n = length(x))
TRP <- rbindlist(lapply(SCOPES, function(sc) {
    R <- in_scope(REPS, sc)
    rbind(R[, c(stratum = "all", label = "Tested SNPs", summ(n_snps))],
          R[, c(stratum = "all", label = "Causal loci", summ(causal_total))],
          R[, c(label = "Causal loci", summ(causal_total)), by = .(stratum = arch)][order(-match(stratum, ARCH_LABELS))],
          R[, c(stratum = "all", label = "Causal loci, % of tested SNPs", summ(causal_pct))],
          R[, c(stratum = "all", label = "K", summ(k_best))],
          R[, .(count = .N), by = k_best][order(k_best)][, .(stratum = "all", label = sprintf("Replicates at K = %d", k_best), n = nrow(R), count)],
          data.table(stratum = "all", label = "Replicates raised to the floor of K = 2 (K of the deposit is 1)",
                     n = nrow(R), count = R[K_authors < 2L, .N]),
          R[, c(stratum = "all", label = "LD window, bp", summ(ld_bp))],
          data.table(stratum = "all", label = "Replicates at a declared 619 bp", n = nrow(R),
                     count = R[distance_source == "declared_median_619", .N]),
          R[, c(stratum = "all", label = "Max VIF of the conditioned RDA", summ(vif))],
          data.table(stratum = "all", label = "Replicates with max VIF >= 10", n = nrow(R), count = R[vif >= 10, .N]),
          use.names = TRUE, fill = TRUE)[, landscape := sc]
}))
fwrite(TRP[, .(landscape, stratum, label, median, min, max, n, count)], file.path(OUT, "t_replicates.tsv"), sep = "\t")

# ================================================================= numbers of Figures S1 and S3
AR <- rbindlist(lapply(seq_len(nrow(LANDS)), function(i) {
    u <- rd(EVAL, LANDS$rem[i], "rdaunc_rda2x_ldb/rank_metrics.tsv")
    c <- rd(EVAL, LANDS$rem[i], "onefit_rda2x_ldb/rank_metrics.tsv")
    chk <- merge(u[method != "RDA", .(seed, method, a = aucpr)], c[method != "RDA", .(seed, method, b = aucpr)], by = c("seed", "method"))
    stopifnot(nrow(chk) == 1200L, isTRUE(all.equal(chk$a, chk$b, tolerance = 0)))
    dcast(rbind(u[, .(seed, method, aucpr)][method != "RDA"], u[method == "RDA", .(seed, method = "RDA uncorrected", aucpr)],
                c[method == "RDA", .(seed, method = "RDA corrected", aucpr)]), seed ~ method, value.var = "aucpr")
}))
AR <- merge(AR, COV[, .(seed, landscape, arch)], by = "seed")
stopifnot(nrow(AR) == 1800L, !anyNA(AR))
s1_row <- function(D) D[, .(n = .N, wins_uncorrected = sum(`RDA uncorrected` > `RDA corrected`),
                            median_aucpr_rda_uncorrected = median(`RDA uncorrected`), median_aucpr_rda_corrected = median(`RDA corrected`),
                            median_aucpr_lfmm2 = median(LFMM), median_aucpr_emmax = median(EMMAX))]
FS1 <- rbindlist(lapply(SCOPES, function(sc) {
    D <- in_scope(AR, sc)
    rbind(s1_row(D)[, architecture := "pooled"], D[, s1_row(.SD), by = .(architecture = arch)][order(match(architecture, ARCH_LABELS))],
          use.names = TRUE)[, landscape := sc]
}))
setcolorder(FS1, c("landscape", "architecture"))
fwrite(FS1, file.path(OUT, "fig_s1_numbers.tsv"), sep = "\t")

CUTS <- c(top_0.001 = 0.1, top_0.0025 = 0.25, top_0.005 = 0.5, top_0.01 = 1, top_0.02 = 2)
CUT <- rbindlist(lapply(seq_len(nrow(LANDS)), function(i) {
    arm <- file.path(LANDS$rem[i], "rdaunc_rda2x_ldb")
    cpm <- rd(EVAL, arm, "calls_per_method.tsv")
    ccb <- rd(EVAL, arm, "calls_combine.tsv")
    stopifnot(ccb[combine == "union" & window_kb > 0, uniqueN(n_called), by = .(seed, rung)][, all(V1 == 1L)],
              cpm[, all(n_called == tp_any_causal + expected_linked + fp_background)],
              ccb[, all(n_called == tp_any_causal + expected_linked + fp_background)])
    cnt <- function(D, lab) D[, .(seed, rung, label = lab, n = n_called, TP = tp_any_causal + expected_linked, FP = fp_background)]
    rbind(cnt(ccb[combine == "union" & window_kb > 0], "1/3 methods"), cnt(ccb[combine == "ge2" & window_kb > 0], "2/3 methods"),
          cnt(ccb[combine == "all3" & window_kb > 0], "3/3 methods"),
          cpm[, .(seed, rung, label = method, n = n_called, TP = tp_any_causal + expected_linked, FP = fp_background)])[, landscape := LANDS$land[i]]
}))
CUT[, cutoff := unname(CUTS[rung])]
stopifnot(!anyNA(CUT$cutoff), setequal(CUT$label, RULES), CUT[, .N, by = .(landscape, cutoff, label)][, all(N == 600L)])
FS3 <- rbindlist(lapply(SCOPES, function(sc)
    in_scope(CUT, sc)[, .(n = .N, TP_median = as.numeric(median(TP)), FP_median = as.numeric(median(FP)), n_empty = sum(n == 0L)),
                      by = .(cutoff, label)][, landscape := sc]))
FS3[, set := set_lab(label)]
FS3 <- FS3[order(match(landscape, SCOPES), cutoff, match(label, PANEL_ORDER))]
fwrite(FS3[, .(landscape, cutoff_pct = cutoff, set, n, TP_median, FP_median, n_empty)], file.path(OUT, "fig_s3_numbers.tsv"), sep = "\t")

# ================================================================= numbers.tsv: the published columns + landscape + set
num <- list()
add_num <- function(D, figure, source) {
    D <- as.data.table(D)
    D[, `:=`(figure = figure, source = source)]
    num[[length(num) + 1L]] <<- D
}
KEY_OF <- setNames(names(SET_LAB), unname(SET_LAB))
for (sc in SCOPES) {
    D <- in_scope(DET, sc)
    st <- function(x) x[, .(n = .N, n_empty = as.numeric(sum(n == 0)), TP_median = median(TP), TP_q25 = quantile(TP, .25, names = FALSE),
                            TP_q75 = quantile(TP, .75, names = FALSE), FP_median = median(FP), FP_q25 = quantile(FP, .25, names = FALSE),
                            FP_q75 = quantile(FP, .75, names = FALSE), prec_cl_median = median(prec_cl, na.rm = TRUE), F1_median = median(F1, na.rm = TRUE)),
                         by = .(rule = as.character(label))]
    x <- rbind(st(D)[, stratum := "pooled"], D[, st(.SD), by = .(stratum = as.character(arch))], use.names = TRUE)
    add_num(melt(x, id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value")[, landscape := sc],
            "Figure 2a / S4 / S5 detection", "detection_per_seed.tsv")
    A <- melt(in_scope(AR, sc)[, .(seed, arch, LFMM, EMMAX, RDA = `RDA uncorrected`)], id.vars = c("seed", "arch"), variable.name = "rule", value.name = "aucpr")
    sa <- function(x) x[, .(n = .N, median = median(aucpr), q25 = quantile(aucpr, .25, names = FALSE), q75 = quantile(aucpr, .75, names = FALSE)), by = .(rule = as.character(rule))]
    y <- rbind(sa(A)[, stratum := "pooled"], A[, sa(.SD), by = .(stratum = as.character(arch))], use.names = TRUE)
    add_num(melt(y, id.vars = c("rule", "stratum", "n"), variable.name = "metric", value.name = "value")[, landscape := sc],
            "Figure S2 AUC-PR per method (RDA without structure correction)", "rdaunc_rda2x_ldb/rank_metrics.tsv")
}
add_num(TR[, .(landscape, rule = panel, stratum = "pooled", engine = unname(ENG_LAB[engine]), metric = unname(STATS[stat]), value, lo, hi,
               expected = fcase(stat %in% c("r1", "r2", "r3", "r4", "r5", "r6"), 100 / N_PANEL, stat %in% c("r12", "r56"), 200 / N_PANEL, default = 300 / N_PANEL),
               n)],
    "Figure 2b / t_ranks within-replicate rank: % (95% bootstrap CI, 10,000 resamples)", "rank_long.tsv")
add_num(TRS[, .(landscape, rule = panel, stratum = paste0(kind, ": ", block_show(stratum)), engine = "median of 3", metric = unname(STATS[stat]), value, lo, hi, n)],
    "t_ranks_strata / s_rank_strata / s_rank_architecture", "rank_long.tsv + detection_per_seed.tsv")
add_num(ORD[, .(landscape = "all", rule = panel, stratum = "pooled", engine = "median of 3", metric = "set_order", value = seq_len(.N))],
    "candidate-set order (every per-set figure and table)", "rank_long.tsv (ranks 4-6, the 1,800 replicates)")
add_num(melt(TOA[panel != "causal", .(landscape, rule = panel, stratum = "pooled", engine = offset_method, n, median_accuracy, median_difference, median_loss, p90_loss,
                                      pct_replicates_difference_ge_0)],
             id.vars = c("landscape", "rule", "stratum", "engine", "n"), variable.name = "metric", value.name = "value"),
    "t_offset_accuracy: accuracy, difference from the causal loci, accuracy lost against the best set", "rank_long.tsv + offset_delta_per_seed.tsv")
add_num(TOA[panel == "causal", .(landscape, rule = "causal loci", stratum = "pooled", engine = offset_method, n, metric = "median_accuracy", value = median_accuracy)],
    "t_offset_accuracy: causal loci (reference)", "rank_long.tsv + offset_delta_per_seed.tsv")
add_num(FS3[, .(landscape, rule = label, stratum = sprintf("cutoff %s%%", cutoff), n, TP_median, FP_median, n_empty = as.numeric(n_empty))][
            , melt(.SD, id.vars = c("landscape", "rule", "stratum", "n"), variable.name = "metric", value.name = "value")],
    "Figure S3 detection at five cutoffs (LD-decay distance of each replicate; RDA at twice each cutoff)", "rdaunc_rda2x_ldb/{calls_per_method,calls_combine}.tsv")
add_num(melt(copy(FS1)[, wins_uncorrected := as.numeric(wins_uncorrected)], id.vars = c("landscape", "architecture", "n"), variable.name = "metric", value.name = "value")[
            , .(landscape, rule = "RDA uncorrected vs corrected", stratum = architecture, n, metric, value)],
    "Figure S1 AUC-PR of RDA with and without structure correction", "{rdaunc,onefit}_rda2x_ldb/rank_metrics.tsv")
add_num(TRP[!is.na(median), .(landscape, rule = "replicates", stratum, n, metric = label, value = median, lo = min, hi = max)],
    "t_replicates: median, lo = minimum, hi = maximum", "detection_per_seed.tsv + rank_metrics.tsv + ld_window_by_replicate.tsv + RDA diagnostics")
add_num(TRP[!is.na(count), .(landscape, rule = "replicates", stratum, n, metric = label, value = as.numeric(count))],
    "t_replicates: counts", "mvp_seeds.tsv + ld_window_by_replicate.tsv + RDA diagnostics")
NUMS <- rbindlist(num, use.names = TRUE, fill = TRUE)
NUMS[, set := fifelse(rule %in% names(SET_LAB), unname(SET_LAB[rule]), rule)]
NUMS[, metric := as.character(metric)]
# metric_key: the keys of the published numbers.tsv for the places it had, the metric itself otherwise
KEY_STAT <- c("rank 1" = "top1", "ranks 1–2" = "top2", "ranks 1–3" = "top3", "rank 6" = "bottom1", "ranks 5–6" = "bottom2", "ranks 4–6" = "bottom3")
NUMS[, metric_key := fifelse(metric %in% names(KEY_STAT), unname(KEY_STAT[metric]), metric)]
setcolorder(NUMS, c("landscape", "figure", "rule", "set", "stratum", "engine", "metric", "metric_key", "value", "lo", "hi", "expected", "n", "source"))
fwrite(NUMS[, c("landscape", "figure", "rule", "set", "stratum", "engine", "metric", "metric_key", "value", "lo", "hi", "expected", "n", "source"), with = FALSE],
       file.path(OUT, "numbers.tsv"), sep = "\t")
message(sprintf("numbers.tsv: %d rows", nrow(NUMS)))

# ================================================================= markdown (the published layouts, one section per scope)
minus <- function(s) sub("^-", "−", s)
f_pct <- function(x) minus(formatC(x, format = "f", digits = 1))
f_dif <- function(x) minus(formatC(x, format = "f", digits = 3))
f_ci  <- function(v, l, h) sprintf("%s [%s–%s]", f_pct(v), f_pct(l), f_pct(h))
f_ph  <- function(p) trimws(formatC(signif(p, 2), format = "fg", digits = 2))
SUP <- c(`-` = "⁻", `0` = "⁰", `1` = "¹", `2` = "²", `3` = "³", `4` = "⁴", `5` = "⁵", `6` = "⁶", `7` = "⁷", `8` = "⁸", `9` = "⁹")
f_p <- function(p) vapply(p, function(x) {
    if (x >= 0.001) return(formatC(signif(x, 3), format = "fg", digits = 3, flag = "#"))
    e <- floor(log10(x)); m <- x / 10^e
    if (round(m, 2) >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("%s × 10%s", formatC(m, format = "f", digits = 2), paste(SUP[strsplit(as.character(e), "")[[1]]], collapse = ""))
}, character(1))
md_row <- function(x) paste0("| ", paste(x, collapse = " | "), " |")
RANK_COLS <- unname(c("Rank 1", "Rank 2", "Rank 3", "Rank 4", "Rank 5", "Rank 6", "Ranks 1–2", "Ranks 1–3", "Ranks 4–6", "Ranks 5–6"))
fmt_count <- function(x) ifelse(x == round(x), format(round(x), big.mark = ",", trim = TRUE), formatC(x, format = "f", digits = 1, big.mark = ","))
fmt_cell <- function(med, lo, hi, kind) {
    f <- switch(kind, count = fmt_count, pct1 = function(x) formatC(x, format = "f", digits = 1), pct2 = function(x) formatC(x, format = "f", digits = 2),
                int = function(x) format(round(x), big.mark = ",", trim = TRUE), dec2 = function(x) formatC(x, format = "f", digits = 2))
    sprintf("%s (%s–%s)", f(med), f(lo), f(hi))
}
kind_a <- c(n = "count", TP = "count", FP = "count", prec_pct = "pct1", n_causal = "count", rec_pct = "pct1", n_linked = "count")
md_table_a <- function(S) {
    S <- copy(S)
    S[, cell := fmt_cell(median, min, max, kind_a[metric]), by = seq_len(nrow(S))]
    W <- dcast(S, set_o + set ~ m_o, value.var = "cell")
    c(md_row(c("Candidate set", unname(metrics_a))), md_row(rep("---", length(metrics_a) + 1)),
      vapply(seq_len(nrow(W)), function(i) md_row(c(W$set[i], unlist(W[i, -(1:2)]))), character(1)))
}
prec_note <- function(S) {
    e <- S[metric == "prec_pct" & n < max(S$n)]
    if (!nrow(e)) return(sprintf("Every value is taken over %s replicates.", max(S$n)))
    sprintf("Every value is taken over %s replicates, except the precision of %s: empty there, every other column counts 0.",
            max(S$n), paste(sprintf("the %s (%s replicates)", e$set, e$n), collapse = ", "))
}
writeLines(c("# t_detection", "", "Cells: median (minimum–maximum) over the replicates. Values: `t_detection.tsv`; built by `benchmarks/mvp_land_tables.R`.", "",
             unlist(lapply(SCOPES, function(sc) { S <- TD[landscape == sc]; c(sprintf("## %s (%d replicates)", sc, max(S$n)), "", md_table_a(S), "", prec_note(S), "") }))),
           file.path(OUT, "t_detection.md"))
writeLines(c("# t_detection_arch", "", "Cells: median (minimum–maximum) over the replicates. Values: `t_detection_arch.tsv`; built by `benchmarks/mvp_land_tables.R`.", "",
             unlist(lapply(SCOPES, function(sc) unlist(lapply(ARCH_LABELS, function(a) {
                 S <- TDA[landscape == sc & arch == a]
                 c(sprintf("## %s, %s (%d replicates)", sc, a, max(S$n)), "", md_table_a(S), "", prec_note(S), "")
             }))))), file.path(OUT, "t_detection_arch.md"))
writeLines(c("# t_replicates", "", "Cells: median (minimum–maximum), or a count of replicates. Values: `t_replicates.tsv`; built by `benchmarks/mvp_land_tables.R`.", "",
             unlist(lapply(SCOPES, function(sc) {
                 S <- TRP[landscape == sc]
                 c(sprintf("## %s", sc), "", md_row(c("Quantity", "Stratum", "Replicates", "Median (minimum–maximum) or count")), md_row(rep("---", 4)),
                   vapply(seq_len(nrow(S)), function(i) with(S[i], md_row(c(label, stratum, format(n, big.mark = ","),
                       if (!is.na(count)) format(count, big.mark = ",") else fmt_cell(median, min, max,
                           if (grepl("%", label)) "pct2" else if (grepl("VIF", label)) "dec2" else "int")))), character(1)), "")
             }))), file.path(OUT, "t_replicates.md"))
rank_rows_ci <- function(D) {
    D <- copy(D)[, cell := f_ci(value, lo, hi)]
    W <- dcast(D, set_o + panel ~ stat_o, value.var = "cell")
    vapply(seq_len(nrow(W)), function(i) md_row(c(set_lab(W$panel[i]), unlist(W[i, -(1:2)]))), character(1))
}
omni_line <- function(sc, e) with(OMNI[landscape == sc & engine == e], sprintf(
    "Whole table (6 × 6 counts of candidate set by rank): χ² = %s; permutation p = %s (k = %d of %s permutations), Holm-adjusted over the four rankings p = %s.",
    formatC(chi2, format = "f", digits = 1, big.mark = ","), f_ph(p), k, format(permutations, big.mark = ","), f_ph(p_holm)))
writeLines(c("# t_ranks", "",
             paste0("Cells: % of replicates [95% bootstrap CI]; expected if the six candidate sets were interchangeable: 16.7% for one rank, 33.3% for two, 50% for three. ",
                    "Last line of each ranking: the χ² of its table of candidate set by rank against 10,000 permutations of the candidate-set labels within each replicate ",
                    "(P = (k + 1) / 10,001, Holm-adjusted over the four rankings); values in `t_ranks_omnibus.tsv`. Built by `benchmarks/mvp_land_tables.R`."), "",
             unlist(lapply(SCOPES, function(sc) unlist(lapply(ENGS, function(e) {
                 c(sprintf("## %s, %s", sc, ENG_LAB[[e]]), "", md_row(c("Candidate set", RANK_COLS)), md_row(rep("---", 11)),
                   rank_rows_ci(TR[landscape == sc & engine == e]), "", omni_line(sc, e), "")
             }))))), file.path(OUT, "t_ranks.md"))
KIND_TITLE <- c(`genetic architecture` = "Genetic architecture", `demographic scenario` = "Demographic scenario", `design cell` = "Design cell")
rank_rows <- function(D, lead) {
    D <- copy(D)[, cell := f_ci(value, lo, hi)]
    W <- dcast(D, set_o + panel ~ stat_o, value.var = "cell")
    vapply(seq_len(nrow(W)), function(i) md_row(c(lead(i), set_lab(W$panel[i]), unlist(W[i, -(1:2)]))), character(1))
}
writeLines(c("# t_ranks_strata", "",
             "Median of the three offset methods. Cells: % of the stratum's replicates [95% bootstrap CI]; n = replicates in the stratum. Values: `t_ranks_strata.tsv`; built by `benchmarks/mvp_land_tables.R`.", "",
             unlist(lapply(SCOPES, function(sc) unlist(lapply(names(KIND_TITLE), function(k) {
                 D <- TRS[landscape == sc & kind == k]
                 c(sprintf("## %s, %s", sc, KIND_TITLE[[k]]), "", md_row(c(KIND_TITLE[[k]], "Candidate set", RANK_COLS)), md_row(rep("---", 12)),
                   unlist(lapply(STRATA[kind == k, stratum], function(s) {
                       Ds <- D[stratum == s]
                       rank_rows(Ds, function(i) if (i == 1L) sprintf("%s (n = %d)", block_show(s), Ds$n[1]) else "")
                   })), "")
             }))))), file.path(OUT, "t_ranks_strata.md"))
writeLines(c("# t_offset_accuracy", "",
             "Accuracy = offset accuracy of the set, or of the causal loci in the last row; difference = the set's accuracy − the causal loci's accuracy in the same replicate; loss = the accuracy of the replicate's most accurate set − the set's accuracy. Each under the ranking named in the heading. A set without an offset (fewer than 3 SNPs) has no accuracy and no loss; n says so. Values: `t_offset_accuracy.tsv`; built by `benchmarks/mvp_land_tables.R`.", "",
             unlist(lapply(SCOPES, function(sc) unlist(lapply(c("median of 3", "gradientForest", "LFMM2", "RDA"), function(e) {
                 D <- TOA[landscape == sc & engine == e]
                 c(sprintf("## %s, %s", sc, ENG_LAB[[e]]), "",
                   md_row(c("Candidate set", "Replicates", "Median accuracy", "Lowest accuracy", "Highest accuracy", "Median difference", "Minimum", "Maximum",
                            "Replicates with difference ≥ 0, %", "Median loss", "90th percentile of loss")), md_row(rep("---", 11)),
                   vapply(seq_len(nrow(D)), function(i) with(D[i], if (is.na(median_difference))
                       md_row(c(set, n, f_dif(median_accuracy), f_dif(min_accuracy), f_dif(max_accuracy), "—", "—", "—", "—", "—", "—")) else
                       md_row(c(set, n, f_dif(median_accuracy), f_dif(min_accuracy), f_dif(max_accuracy), f_dif(median_difference), f_dif(min_difference),
                                f_dif(max_difference), f_pct(pct_replicates_difference_ge_0), f_dif(median_loss), f_dif(p90_loss)))), character(1)), "")
             }))))), file.path(OUT, "t_offset_accuracy.md"))
writeLines(c("# t_rank_tests", "", "Kruskal–Wallis tests of each set's within-replicate rank (median of the three offset methods); Holm p over the 12 tests of each landscape. Values: `t_rank_tests.tsv`; built by `benchmarks/mvp_land_tables.R`.", "",
             unlist(lapply(SCOPES, function(sc) {
                 D <- KWT[landscape == sc]
                 c(sprintf("## %s", sc), "", md_row(c("Test", "Candidate set", "H", "df", "p", "Holm p", "n", "Groups")), md_row(rep("---", 8)),
                   vapply(seq_len(nrow(D)), function(i) with(D[i], md_row(c(KW_LAB[[factor]], set_lab(panel), formatC(H, format = "f", digits = 2), df, f_p(p), f_p(p_holm), n, groups))), character(1)), "")
             }))), file.path(OUT, "t_rank_tests.md"))
writeLines(c("# t_ranks_conditions", "", "Median of the three offset methods. Rank columns: % of the condition's replicates [95% bootstrap CI]. Loss: median accuracy lost against the replicate's best set [95% bootstrap CI]; n = replicates, fewer where a set has no offset. Values: `t_ranks_conditions.tsv`; built by `benchmarks/mvp_land_tables.R`.", "",
             unlist(lapply(SCOPES, function(sc) unlist(lapply(SEL_LAB, function(se) {
                 D <- TC[landscape == sc & selection == se]
                 W <- dcast(copy(D)[, cell := fifelse(stat == "loss", sprintf("%s [%s–%s], n = %d", f_dif(value), f_dif(lo), f_dif(hi), n), f_ci(value, lo, hi))],
                            set_o + panel ~ st_o, value.var = "cell")
                 c(sprintf("## %s, %s", sc, se), "", md_row(c("Candidate set", unname(STAT_C))), md_row(rep("---", 1 + length(STAT_C))),
                   vapply(seq_len(nrow(W)), function(i) md_row(c(set_lab(W$panel[i]), unlist(W[i, -(1:2)]))), character(1)), "")
             }))))), file.path(OUT, "t_ranks_conditions.md"))
writeLines(c("# t_set_pairs", "", "Paired differences within replicates, first set minus second, with 95% bootstrap CI over replicates. Ranks 4–6: percentage points of replicates. Accuracy: the difference in offset accuracy; n excludes replicates in which a set has no offset. A positive accuracy difference favours the first set; a positive ranks 4–6 difference counts against it. Values: `t_set_pairs.tsv`; built by `benchmarks/mvp_land_tables.R`.", "",
             unlist(lapply(SCOPES, function(sc) unlist(lapply(unname(ENG_LAB[c(4, 1, 2, 3)]), function(e) {
                 D <- TP[landscape == sc & offset_method == e]
                 c(sprintf("## %s, %s", sc, e), "", md_row(c("Pair", "Measure", "Value [95% CI]", "n")), md_row(rep("---", 4)),
                   vapply(seq_len(nrow(D)), function(i) with(D[i], md_row(c(pair, measure,
                       if (grepl("accuracy,", measure)) sprintf("%s [%s–%s]", f_dif(value), f_dif(lo), f_dif(hi)) else f_ci(value, lo, hi), n))), character(1)), "")
             }))))), file.path(OUT, "t_set_pairs.md"))
message("tables written to ", OUT)

# ================================================================= the gate: SS-Clines rows against the published tables
if (dir.exists(REF)) {
    GATE <- list()
    rs <- function(f) fread(f, sep = "\t", colClasses = "character", quote = "")
    cmp <- function(name, new, ref, keys, cols) {
        stopifnot(!anyDuplicated(ref[, ..keys]), !anyDuplicated(new[, ..keys]))
        x <- merge(ref, new, by = keys, suffixes = c(".ref", ".new"), all = TRUE)
        miss <- x[rowSums(is.na(x)) > 0 & (is.na(get(paste0(cols[1], ".ref"))) | is.na(get(paste0(cols[1], ".new"))))]
        nd <- 0L; mx <- 0
        for (c in cols) {
            a <- x[[paste0(c, ".ref")]]; b <- x[[paste0(c, ".new")]]
            bad <- !(a == b | (is.na(a) & is.na(b)) | (a == "" & is.na(b)) | (is.na(a) & b == ""))
            bad[is.na(bad)] <- TRUE
            nd <- nd + sum(bad)
            na <- suppressWarnings(as.numeric(a)); nb <- suppressWarnings(as.numeric(b))
            if (any(bad)) mx <- max(mx, abs(na - nb)[bad], na.rm = TRUE)
        }
        GATE[[length(GATE) + 1L]] <<- data.table(table = name, rows_published = nrow(ref), rows_new = nrow(new), cells = nrow(x) * length(cols),
                                                 cells_differ = nd, max_abs_diff = mx, rows_unmatched = nrow(miss))
        message(sprintf("gate %-22s %5d rows, %6d cells, %d differ (max %.3g), %d unmatched", name, nrow(ref), nrow(x) * length(cols), nd, mx, nrow(miss)))
    }
    rn <- function(f) { D <- rs(file.path(OUT, f)); D[landscape == "SS-Clines"][, landscape := NULL][] }
    cmp("t_detection", rn("t_detection.tsv"), rs(file.path(REF, "t_detection.tsv")), c("set", "label"), c("median", "min", "max", "n"))
    cmp("t_detection_arch", rn("t_detection_arch.tsv"), rs(file.path(REF, "t_detection_arch.tsv")), c("architecture", "set", "label"), c("median", "min", "max", "n"))
    ref_rp <- rs(file.path(REF, "t_replicates.tsv"))
    cmp("t_replicates", rn("t_replicates.tsv")[paste(stratum, label) %in% ref_rp[, paste(stratum, label)]], ref_rp, c("stratum", "label"), c("median", "min", "max", "n"))
    cmp("t_offset_accuracy", rn("t_offset_accuracy.tsv")[offset_method == "Median of the three"], rs(file.path(REF, "t_offset_accuracy.tsv")), "set",
        c("n", "median_accuracy", "min_accuracy", "max_accuracy", "median_difference", "min_difference", "max_difference", "pct_replicates_difference_ge_0"))
    ref_tr <- rs(file.path(REF, "t_ranks.tsv"))
    cmp("t_ranks", rn("t_ranks.tsv")[statistic %in% ref_tr$statistic], ref_tr, c("offset_method", "set", "statistic"), c("value", "lo", "hi", "n"))
    ref_st <- rs(file.path(REF, "t_ranks_strata.tsv"))
    cmp("t_ranks_strata", rn("t_ranks_strata.tsv")[statistic %in% ref_st$statistic], ref_st, c("stratum_kind", "stratum", "set", "statistic"), c("value", "lo", "hi", "n"))
    cmp("t_ranks_omnibus", rn("t_ranks_omnibus.tsv"), rs(file.path(REF, "t_ranks_omnibus.tsv")), "offset_method", c("chi2", "permutations", "k", "p", "p_holm", "seed"))
    cmp("t_rank_tests", rn("t_rank_tests.tsv"), rs(file.path(REF, "t_rank_tests.tsv")), c("test", "set"), c("H", "df", "p", "p_holm", "n", "groups"))
    # Figure S3 and S1 numbers against the published numbers
    PN <- fread(file.path(EVAL, "figures_ld_ms/numbers.tsv"), sep = "\t")
    p3 <- PN[startsWith(figure, "Figure S3") & metric %in% c("TP_median", "FP_median", "n_empty"), .(label = rule, cutoff = as.numeric(sub("cutoff ([0-9.]+)%", "\\1", stratum)), metric, v = value)]
    n3 <- melt(FS3[landscape == "SS-Clines", .(label, cutoff, TP_median, FP_median, n_empty = as.numeric(n_empty))], id.vars = c("label", "cutoff"), variable.name = "metric", value.name = "w", variable.factor = FALSE)
    x3 <- merge(p3, n3, by = c("label", "cutoff", "metric"))
    stopifnot(nrow(x3) == 90L)
    GATE[[length(GATE) + 1L]] <- data.table(table = "Figure S3 numbers", rows_published = 90L, rows_new = nrow(x3), cells = 90L, cells_differ = x3[abs(v - w) > 1e-9, .N],
                                            max_abs_diff = x3[, max(abs(v - w))], rows_unmatched = 90L - nrow(x3))
    tie_num <- function(fig, name) {
        a <- PN[figure == fig, .(rule, stratum, metric, v = value)]
        b <- NUMS[landscape == "SS-Clines" & figure == fig, .(rule, stratum, metric, w = value)]
        x <- merge(a, b, by = c("rule", "stratum", "metric"), all = TRUE)
        bad <- x[is.na(v) != is.na(w) | (!is.na(v) & !is.na(w) & abs(v - w) > 1e-9)]
        GATE[[length(GATE) + 1L]] <<- data.table(table = name, rows_published = nrow(a), rows_new = nrow(b), cells = nrow(x), cells_differ = nrow(bad),
                                                 max_abs_diff = x[!is.na(v) & !is.na(w), max(abs(v - w))], rows_unmatched = x[is.na(v) | is.na(w), .N] - x[is.na(v) & is.na(w), .N])
    }
    tie_num("Figure 2a / S4 / S5 detection", "numbers.tsv: detection rows")
    tie_num("Figure S2 AUC-PR per method (RDA without structure correction)", "numbers.tsv: AUC-PR rows")
    RN <- fread(file.path(EVAL, "figures_ld_rda/numbers.tsv"), sep = "\t", colClasses = "character")
    rv <- function(k) as.numeric(RN[key == k, value])
    s1 <- FS1[landscape == "SS-Clines" & architecture == "pooled"]
    vc <- REPS[landscape == "SS-Clines"]
    fx <- c(wins = s1$wins_uncorrected - rv("wins_unc_pooled"), unc = s1$median_aucpr_rda_uncorrected - rv("aucpr_RDA_uncorrected_pooled_median"),
            cor = s1$median_aucpr_rda_corrected - rv("aucpr_RDA_corrected_pooled_median"), lfmm = s1$median_aucpr_lfmm2 - rv("aucpr_LFMM_pooled_median"),
            vif_median = median(vc$vif) - rv("vif_median"), vif_ge10 = sum(vc$vif >= 10) - rv("vif_ge10"))
    GATE[[length(GATE) + 1L]] <- data.table(table = "Figure S1 and VIF numbers", rows_published = length(fx), rows_new = length(fx), cells = length(fx),
                                            cells_differ = sum(abs(fx) > 1e-12), max_abs_diff = max(abs(fx)), rows_unmatched = 0L)
    GATE <- rbindlist(GATE)
    fwrite(GATE, file.path(OUT, "gate_ssclines.tsv"), sep = "\t")
    print(GATE)
    if (sum(GATE$cells_differ) + sum(GATE$rows_unmatched) > 0) warning("GATE: SS-Clines rows differ from the published tables, see gate_ssclines.tsv")
    else message("GATE PASSED: every SS-Clines cell equals the published tables")
}
