#!/usr/bin/env Rscript
# =============================================================================
# mvp_method_redundancy.R -- how much INDEPENDENT evidence does each extra method add?
#
# THE QUESTION THIS EXISTS TO ANSWER
# ----------------------------------
# GAPIT ships eight models and the pipeline can run all of them. Four of them --
# GLM, MLM, CMLM, ECMLM -- are a nested family: the same mixed model with progressively
# more structure in the kinship/grouping. Three more -- MLMM, FarmCPU, BLINK -- are the
# multi-locus family, each an iterative variant of the previous. Running all eight and
# taking their UNION inflates the call set with no new evidence; taking their
# INTERSECTION converges on whatever any one of them would have said. Either way the
# combine rule reports a number that looks like corroboration and is not.
#
# So before any subset search, measure the redundancy directly and pick one
# representative per family. Everything here is post-hoc on p-value tables that already
# exist -- no refitting.
#
# TWO MEASURES, because they disagree in an informative way:
#   corr    Spearman of the per-SNP ranking score (-log10 of the min p across the
#           method's traits). Genome-wide, threshold-free. This is the right comparison
#           for RDA, which has ONE pseudo-trait column (climate_multivariate) and no
#           per-axis p-value -- a per-trait correlation cannot include it at all.
#   jaccard Overlap of the CALLED sets at a fixed operating point. Two methods can rank
#           near-identically and still call disjoint sets if their p-value calibration
#           differs, which is exactly RDA's situation (MVP_README.md:283-306: pmax gives
#           correct order, not a correct p-value). Ranking similarity without call
#           similarity is the signature of a calibration problem, not a redundancy one.
#
# Usage:
#   Rscript mvp_method_redundancy.R --methods="EMMAX=/p/a.tsv,..." --truth=FILE \
#       --tag=STR [--outdir=DIR] [--adjust=bonf_0.05] [--k=NULL] [--cut-height=0.3]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(qvalue)
    library(ggplot2)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/pval_threshold.R"))
source(file.path(PIPELINE_ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args <- parse_kv_args(commandArgs(trailingOnly = TRUE))
req  <- function(k) {
    if (is.null(args[[k]]) || !nzchar(args[[k]])) stop("Missing required --", gsub("_", "-", k))
    args[[k]]
}
opt  <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]

TAG      <- req("tag")
TRUTH_F  <- req("truth")
OUTDIR   <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/redundancy"))
ADJUST   <- opt("adjust", "bonf_0.05")
CUT_H    <- as.numeric(opt("cut_height", "0.3"))

spec <- strsplit(req("methods"), ",", fixed = TRUE)[[1]]
methods <- list()
for (s in spec) {
    kv <- strsplit(s, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Bad --methods entry (expected NAME=PATH): ", s)
    if (!file.exists(kv[2])) stop("Method p-value table not found: ", kv[2])
    methods[[kv[1]]] <- kv[2]
}
if (length(methods) < 2) stop("--methods needs at least two methods")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
truth <- load_truth(TRUTH_F)

ap  <- strsplit(ADJUST, "_", fixed = TRUE)[[1]]
adj <- paste(head(ap, -1), collapse = "_"); val <- tail(ap, 1)

score <- list(); called <- list(); aucs <- list()
for (m in names(methods)) {
    lp <- load_pvalues(methods[[m]], "all")
    rk <- rank_keys(lp$pv, lp$trait_cols)
    tc <- testable_causal(lp$pv, lp$trait_cols, truth)
    # Ranking score per SNP: -log10 of the min p across the method's traits. Same
    # quantity rank_keys() orders on, kept as a value so Spearman can use it.
    minp <- do.call(pmin, c(lapply(lp$trait_cols, function(x) lp$pv[[x]]), list(na.rm = TRUE)))
    score[[m]]  <- setNames(-log10(pmax(minp, .Machine$double.xmin)), lp$pv$key)
    called[[m]] <- call_by_threshold(lp$pv, lp$trait_cols, adj, val, quiet = TRUE)$called
    aucs[[m]]   <- data.table(method = m, auc_pr = auc_pr_from_rank(rk$keys, truth, length(tc)),
                              n_testable = length(tc), n_called = length(called[[m]]))
    message(sprintf("INFO: %-8s AUC-PR=%.4f  calls@%s=%d", m, aucs[[m]]$auc_pr, ADJUST,
                    length(called[[m]])))
}
AUC <- rbindlist(aucs)

# Common SNP set. Methods can drop SNPs independently (GAPIT's own filtering,
# monomorphic-in-subset), and correlating over a union with NAs would let coverage
# differences masquerade as ranking differences.
common <- Reduce(intersect, lapply(score, names))
message("INFO: ", length(common), " SNPs scored by all ", length(methods), " methods")
if (length(common) < 100) stop("Fewer than 100 SNPs are common to all methods -- refusing to correlate")

M  <- do.call(cbind, lapply(score, function(v) v[common]))
CS <- suppressWarnings(cor(M, method = "spearman", use = "pairwise.complete.obs"))

jac <- function(a, b) {
    u <- length(union(a, b))
    if (u == 0) return(NA_real_)
    length(intersect(a, b)) / u
}
JA <- outer(names(called), names(called), Vectorize(function(i, j) jac(called[[i]], called[[j]])))
dimnames(JA) <- list(names(called), names(called))

long <- rbindlist(list(
    as.data.table(as.table(CS))[, .(method_a = V1, method_b = V2, measure = "spearman", value = N)],
    as.data.table(as.table(JA))[, .(method_a = V1, method_b = V2, measure = "jaccard",  value = N)]
))
long <- long[method_a != method_b]
fwrite(long, file.path(OUTDIR, paste0(TAG, "_method_similarity.tsv")), sep = "\t")

# --------------------------------------------------------------------- families
# Clustering on the RANKING correlation, not on Jaccard: Jaccard at a fixed threshold is
# confounded by call-set size (a method calling 3 SNPs cannot overlap much with one calling
# 300, however similar their rankings), and the redundancy question is about the evidence,
# not about where each method happens to put its cutoff.
hc  <- hclust(as.dist(1 - CS), method = "average")
grp <- cutree(hc, h = CUT_H)
if (!is.null(args$k) && nzchar(args$k)) grp <- cutree(hc, k = as.integer(args$k))

fam <- merge(data.table(method = names(grp), cluster = as.integer(grp)), AUC, by = "method")
# Representative = highest AUC-PR in the family. AUC-PR rather than F1 because it is
# threshold-free, so the choice of representative does not depend on the operating point
# the subset search has not chosen yet.
fam[, representative := auc_pr == max(auc_pr, na.rm = TRUE), by = cluster]
setorder(fam, cluster, -auc_pr)
fwrite(fam, file.path(OUTDIR, paste0(TAG, "_method_clusters.tsv")), sep = "\t")

message("\n=== families at 1 - rho <= ", CUT_H, " ===")
print(fam[, .(method, cluster, auc_pr = round(auc_pr, 4), n_called, representative)])
message("\nRepresentatives: ", paste(fam[representative == TRUE, method], collapse = ", "))

# ------------------------------------------------------------------- figures
ord <- hc$labels[hc$order]
hm  <- long[measure == "spearman"][, `:=`(method_a = factor(method_a, ord),
                                          method_b = factor(method_b, ord))]
p <- ggplot(hm, aes(method_a, method_b, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), size = 2.6, colour = "white") +
    scale_fill_gradient(low = CLINEGO_THRESHOLD, high = CLINEGO_RETAINED, limits = c(-1, 1)) +
    labs(title = paste0("Method ranking redundancy, Spearman of -log10 min-p (", TAG, ")"),
         x = NULL, y = NULL, fill = "rho") +
    theme_clinego() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTDIR, paste0(TAG, "_redundancy_heatmap.png")), p, width = 9, height = 7.5, dpi = 150)
ggsave(file.path(OUTDIR, paste0(TAG, "_redundancy_heatmap.svg")), p, width = 9, height = 7.5)

png(file.path(OUTDIR, paste0(TAG, "_redundancy_dendrogram.png")), width = 1400, height = 900, res = 150)
plot(hc, main = paste0("Method families (1 - Spearman rho), ", TAG), xlab = "", sub = "")
abline(h = CUT_H, col = CLINEGO_REMOVED, lty = 2)
dev.off()

message("INFO: redundancy outputs in ", OUTDIR)
