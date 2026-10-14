# =============================================================================
# mvp_rank_ties.R -- the one tie rule for the within-replicate rank of the six candidate sets.
# source() it (needs data.table). Used by mvp_panel_rank.R (rank_summary.tsv, V4) and
# mvp_ms_sim_figures.R (Fig 2b / S7 / S8 bars, Table S1).
#
# RULE (user, 2026-10-14): split a tie evenly. Ranks come from frank(ties.method = "average"), so
# the m sets of one tie share the average k + (m - 1) / 2 of the ranks k .. k + m - 1 they
# occupy. Each of them counts 1/m at each of those ranks: a set tied for ranks k and k + 1 counts
# 0.5 at each. It replaces three conventions that disagreed on the tied units:
#   - mvp_panel_rank.R counted last = rank >= 6, last two = rank >= 5, so a tie at 5.5 was never
#     "last" and a tie at 4.5 never "last two";
#   - the figures drew round(rank), and R rounds half to even (4.5 -> 4, 5.5 -> 6);
#   - the ad hoc script behind the first Table S1 counted ties fully on the inclusive side.
#
# A tie group is found from the ranks alone: the sets of a unit that share one average rank. Two
# different ties occupy disjoint rank ranges, so they never share an average; a three-way tie,
# whose average is a whole number, is still found, because m is counted, not read from a ".5".
# No float comparison of accuracies is involved, so the rule gives the same rows on in-memory
# ranks and on ranks read back from rank_long.tsv.
# =============================================================================

#' One row per (unit, set, tied rank): the candidate set's rank turned into (rank, weight) rows.
#'   D     data.table, one row per candidate set and unit, columns `unit`, `set` and rank
#'         (the average rank from frank(ties.method = "average"))
#'   unit  key columns of one ranking, e.g. c("seed", "engine")
#' Returns the unit and set columns, rank, tie_n (sets sharing the rank), r (whole rank), w (1 / tie_n).
#' Untied sets get one row with r = rank and w = 1.
rank_split <- function(D, unit, set = "panel") {
    X <- D[, c(unit, set, "rank"), with = FALSE]
    if (anyDuplicated(X[, c(unit, set), with = FALSE])) stop("rank_split: duplicated unit x set rows")
    X[, tie_n := .N, by = c(unit, "rank")]
    X <- X[rep(seq_len(.N), tie_n)]
    X[, r := rank - (tie_n - 1) / 2 + seq_len(.N) - 1, by = c(unit, set)]
    if (any(abs(X$r - round(X$r)) > 1e-9)) stop("rank_split: a tied rank is not a whole number")
    X[, `:=`(r = as.integer(round(r)), w = 1 / tie_n)]
    n_set <- uniqueN(D[[set]])
    ok <- c(weight_per_set  = X[, abs(sum(w) - 1) < 1e-9, by = c(unit, set)][, all(V1)],
            weight_per_rank = X[, abs(sum(w) - 1) < 1e-9, by = c(unit, "r")][, all(V1)],
            every_rank_once = X[, setequal(r, seq_len(n_set)), by = unit][, all(V1)])
    if (!all(ok)) stop("rank_split: invariant failed: ", paste(names(ok)[!ok], collapse = ", "))
    X[]
}

#' One row per (unit, set) from rank_split() rows: the set's weight at the top k / bottom k ranks
#' (k = 1, 2, 3; each in {0, 1/m, ..., 1}) and its mean rank (sum of w * r = the average rank).
rank_indicators <- function(X, unit, set = "panel", n_set) {
    X[, .(top1 = sum(w[r <= 1L]), top2 = sum(w[r <= 2L]), top3 = sum(w[r <= 3L]),
          bottom1 = sum(w[r >= n_set]), bottom2 = sum(w[r >= n_set - 1L]),
          bottom3 = sum(w[r >= n_set - 2L]), mean_rank = sum(w * r)), by = c(unit, set)]
}
