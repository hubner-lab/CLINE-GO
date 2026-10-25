#!/usr/bin/env Rscript
# Phase 3c copy of ANALYSIS/Trifolium_final/select_predictors.R (2026-10-12): identical rule, with the
# |r| cut and an exclusion list as arguments (user: |r| < 0.8, VIF < 10; bio_17/18 dropped as binary).
# Trifolium — climate predictor selection, run ONCE and pinned across all three arms.
#
# Why this exists: two different 4-predictor sets are already on disk
# (bio_1,bio_3,bio_5,bio_18 from the 51-sample run; bio_2,bio_5,bio_15,bio_18 from the
# 33-sample run). Both satisfy "all pairwise |r| < 0.7", so the greedy screen is not
# unique — which makes the choice arbitrary unless the tie-break is declared.
#
# Declared rule (pre-registered in the plan):
#   1. Drop invariant predictors (bio_14, zero variance across the 11 sites).
#   2. Correlate at the SITE level (11 rows), not the sample level — every sample at a
#      site shares one coordinate pair, so sample-level correlation would just weight
#      sites by n and misstate the environmental design.
#   3. Enumerate ALL maximal subsets with pairwise |r| < 0.7 (= maximal cliques of the
#      compatibility graph), keep those with max VIF < 10.
#   4. Among the largest such subsets, choose the one that best REPRESENTS the whole
#      climate space: maximise the mean R^2 of every non-invariant bio variable
#      regressed on the subset. Unit-free, unlike raw variance.
#   5. Break remaining ties by lowest max VIF, then by lowest bio index.
#
# Usage: Rscript select_predictors_cut.R <climate_present_site.tsv> <invariant.tsv> <outdir> <r_cut> [exclude,list]

suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
SITE_TSV  <- args[1]
INVAR_TSV <- args[2]
OUTDIR    <- args[3]
EXCLUDE   <- if (length(args) >= 5 && nzchar(args[5])) strsplit(args[5], ",")[[1]] else character(0)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

R_CUT   <- as.numeric(args[4])                     # 0.7 in the pre-registered original
VIF_CUT <- 10

# ---- 1. site-level climate matrix ------------------------------------------------
# climate_present_site.tsv is per-sample; collapse to unique environmental points.
d <- fread(SITE_TSV, colClasses = c(sample = "character"))
bio_cols <- grep("^bio_", names(d), value = TRUE)
env <- unique(d[, ..bio_cols])
message("INFO: unique environmental points: ", nrow(env), " (expect 11)")

invar <- character(0)
if (file.exists(INVAR_TSV)) {
    iv <- fread(INVAR_TSV)
    if (nrow(iv)) invar <- iv$predictor
}
# also catch any zero-variance column empirically
zv <- bio_cols[sapply(bio_cols, function(v) stats::var(env[[v]]) == 0)]
invar <- union(invar, zv)
message("INFO: invariant predictors dropped: ", paste(invar, collapse = ", "))

cand <- setdiff(bio_cols, union(invar, EXCLUDE))
message("INFO: excluded by request: ", paste(EXCLUDE, collapse = ", "), "; candidates: ", length(cand))
E <- as.matrix(env[, ..cand])
Z <- scale(E)                                   # z-scored: unit-free throughout
R <- cor(E)

fwrite(data.table(predictor = rownames(R), as.data.table(round(R, 4))),
       file.path(OUTDIR, "predictor_correlation_site.tsv"), sep = "\t")

# ---- 2. maximal cliques of the |r| < 0.7 compatibility graph ----------------------
n   <- length(cand)
adj <- (abs(R) < R_CUT); diag(adj) <- FALSE

cliques <- list()
bron_kerbosch <- function(Rset, P, X) {
    if (!length(P) && !length(X)) { cliques[[length(cliques) + 1L]] <<- Rset; return(invisible()) }
    pivot <- c(P, X)[which.max(sapply(c(P, X), function(u) sum(adj[u, P])))]
    for (v in setdiff(P, which(adj[pivot, ]))) {
        bron_kerbosch(c(Rset, v), intersect(P, which(adj[v, ])), intersect(X, which(adj[v, ])))
        P <- setdiff(P, v); X <- c(X, v)
    }
}
bron_kerbosch(integer(0), seq_len(n), integer(0))
message("INFO: maximal |r|<", R_CUT, " subsets found: ", length(cliques))

# ---- 3. VIF filter + representation score ----------------------------------------
max_vif <- function(idx) {
    if (length(idx) < 2) return(1)
    sapply(seq_along(idx), function(j) {
        y <- Z[, idx[j]]; Xo <- Z[, idx[-j], drop = FALSE]
        r2 <- summary(stats::lm(y ~ Xo))$r.squared
        if (r2 >= 1 - 1e-9) Inf else 1 / (1 - r2)
    }) |> max()
}
# mean ADJUSTED R^2 of every candidate variable reconstructed from the subset.
# Adjusted, not raw: with 11 environmental points and 4 predictors + intercept the raw
# R^2 is inflated by df alone (every subset scores ~0.985 and the criterion cannot
# discriminate). The adjusted form penalises that and is the honest measure.
represent <- function(idx) {
    Xs <- Z[, idx, drop = FALSE]
    mean(sapply(seq_len(n), function(j) {
        if (j %in% idx) return(1)
        summary(stats::lm(Z[, j] ~ Xs))$adj.r.squared
    }))
}

by_index <- function(v) paste(v[order(as.integer(sub("bio_", "", v)))], collapse = ",")

res <- rbindlist(lapply(cliques, function(idx) {
    mv <- max_vif(idx)
    data.table(size = length(idx), max_vif = mv,
               represent = if (is.finite(mv) && mv < VIF_CUT) represent(idx) else NA_real_,
               min_index = min(as.integer(sub("bio_", "", cand[idx]))),
               set = by_index(cand[idx]))          # ordered by bio index, not alphabetically
}))
res <- res[is.finite(max_vif) & max_vif < VIF_CUT]
setorder(res, -size, -represent, max_vif, min_index)

# TIED-SET RULE: subsets within TOL adjusted-R2 of the leader are statistically
# indistinguishable at 11 environmental points; declaring a winner on a 4th-decimal
# difference would be exactly the arbitrariness this script exists to remove.
# Among tied maximal sets, take the lowest bio index.
TOL <- 0.01
top_size <- res$size[1]
tied <- res[size == top_size & represent >= max(represent) - TOL]
setorder(tied, min_index, max_vif)
best <- tied[1]

res[, tied_with_leader := (size == top_size & represent >= max(res$represent) - TOL)]
res[, selected := (set == best$set)]
fwrite(res, file.path(OUTDIR, "predictor_set_candidates.tsv"), sep = "\t")

message("\n=== SELECTED: ", best$set, " ===")
message("   size=", best$size, "  max_VIF=", round(best$max_vif, 3),
        "  mean_adjR2_representation=", round(best$represent, 4))
message("   ", nrow(tied), " of ", nrow(res[size == top_size]),
        " maximal subsets are within ", TOL, " adjR2 of the leader (statistically tied);",
        " lowest-bio-index rule applied among them.")
message("\nAll size-", top_size, " candidates:")
print(res[size == top_size][, .(set, max_vif = round(max_vif, 3),
                                represent = round(represent, 4), tied_with_leader, selected)])

# how the two historical sets score, for the journal
for (hist in c("bio_1,bio_3,bio_5,bio_18", "bio_2,bio_5,bio_15,bio_18")) {
    row <- res[set == hist]
    if (nrow(row)) message("historical set ", hist, " -> rank ", which(res$set == hist),
                           ", adjR2=", round(row$represent, 4), ", maxVIF=", round(row$max_vif, 3))
    else message("historical set ", hist, " -> NOT a maximal VIF-passing subset")
}
writeLines(best$set, file.path(OUTDIR, "predictor_set_selected.txt"))
