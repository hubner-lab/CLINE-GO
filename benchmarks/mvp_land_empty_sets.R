#!/usr/bin/env Rscript
# =============================================================================
# mvp_land_empty_sets.R -- Job C of the landscapes pack: how often the 2/3 and 3/3 sets are empty or hold fewer than
# 3 SNPs when each method keeps its SNPs at Bonferroni 0.05 or q 0.05 per predictor (lists pooled, RDA on its one p)
# instead of at the fixed 0.25 % cutoff, at each replicate's own LD window (Est-Clines: 39 replicates at a declared
# 619 bp). Reads landscape_ms/empty_sets/<land>/calls_combine.tsv (mvp_land_empty_sets.sh).
# Controls: SS-Clines at window 0 equals remeasure600/rdaunc/calls_combine.tsv cell by cell (same p-tables, same
# harness), and the older 5 kb figures are recomputed from that file.
# Writes landscape_ms/empty_sets.tsv and landscape_ms/empty_sets_control.tsv.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OUT  <- file.path(EVAL, "landscape_ms")
LANDS <- c("SS-Clines" = "ssclines", "SS-Mtn" = "ssmtn", "Est-Clines" = "estclines")
# q 0.10 is extra: with Bonferroni 0.05 it spans the older quoted ranges (5 kb: 2/3 set 45.7-52.7 %, 3/3 set 84.5-93.7 % empty)
RUNGS <- c(bonf_0.05 = "Bonferroni 0.05", qval_0.05 = "q 0.05", qval_0.1 = "q 0.10")
SETS  <- c(ge2 = "2/3 set", all3 = "3/3 set")
rd <- function(...) fread(file.path(...), colClasses = c(seed = "character"))
CC <- rbindlist(lapply(names(LANDS), function(L) rd(OUT, "empty_sets", LANDS[[L]], "calls_combine.tsv")[, landscape := L]))
stopifnot(all(names(RUNGS) %in% CC$rung), all(names(SETS) %in% CC$combine))
SEL <- CC[rung %in% names(RUNGS) & combine %in% names(SETS) & window_kb > 0]
stopifnot(SEL[, .N, by = .(landscape, rung, combine)][, all(N == 600L)])
summ <- function(D) D[, .(n = .N, n_empty = sum(n_called == 0L), pct_empty = 100 * mean(n_called == 0L),
                          n_under_3 = sum(n_called < 3L), pct_under_3 = 100 * mean(n_called < 3L),
                          median_size = median(n_called)), by = .(rung, combine)]
R <- rbind(SEL[, summ(.SD), by = landscape], summ(SEL)[, landscape := "all"], use.names = TRUE)
R[, `:=`(threshold = unname(RUNGS[rung]), set = unname(SETS[combine]))]
R <- R[order(match(landscape, c(names(LANDS), "all")), match(rung, names(RUNGS)), match(combine, names(SETS)))]
fwrite(R[, .(landscape, threshold, set, n, n_empty, pct_empty, n_under_3, pct_under_3, median_size)], file.path(OUT, "empty_sets.tsv"), sep = "\t")
print(R[, .(landscape, threshold, set, n, pct_empty = round(pct_empty, 1), pct_under_3 = round(pct_under_3, 1), median_size)])
# controls
RE <- rd(EVAL, "remeasure600/rdaunc/calls_combine.tsv")[rung %in% names(RUNGS)]
A <- CC[landscape == "SS-Clines" & rung %in% names(RUNGS) & window_kb == 0, .(seed, rung, combine, n_new = n_called)]
B <- RE[window_kb == 0, .(seed, rung, combine, n_old = n_called)]
X <- merge(A, B, by = c("seed", "rung", "combine"))
stopifnot(nrow(X) == nrow(A), nrow(X) == nrow(B))
old5 <- RE[window_kb == 5 & combine %in% names(SETS), .(n = .N, pct_empty = 100 * mean(n_called == 0L), pct_under_3 = 100 * mean(n_called < 3L)), by = .(rung, combine)]
print(old5)
CTRL <- rbind(data.table(check = "SS-Clines window 0 = remeasure600/rdaunc", rung = "bonf_0.05, qval_0.05, qval_0.1", combine = "all", n = nrow(X), n_agree = X[n_new == n_old, .N],
                         pct_empty = NA_real_, pct_under_3 = NA_real_),
              old5[, .(check = "older figures: 5 kb window, remeasure600/rdaunc", rung, combine, n, n_agree = NA_integer_, pct_empty, pct_under_3)])
fwrite(CTRL, file.path(OUT, "empty_sets_control.tsv"), sep = "\t")
message(sprintf("control: %d of %d SS-Clines window-0 rows equal", X[n_new == n_old, .N], nrow(X)))
if (X[n_new != n_old, .N] > 0) stop("SS-Clines window-0 rows differ from remeasure600/rdaunc")
