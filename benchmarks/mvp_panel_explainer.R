#!/usr/bin/env Rscript
# mvp_panel_explainer.R -- plain-language figure set for the marker-panel result.
#
# Three panels, each answering one question a non-specialist would actually ask:
#   A  Which marker set predicts best?          (offset accuracy, higher = better)
#   B  Why? Is it about finding MORE loci?      (recall vs accuracy)
#   C  What is in each set?                     (composition: causal / linked / background)
#
# The offset metric is flipped to "accuracy" = -tau so that HIGHER IS BETTER on
# every axis in this figure. The raw tau is negative by construction (a good
# offset predicts LOW fitness), which is why a ratio of two taus can exceed 100%
# and read as nonsense -- that framing is deliberately not used here.

suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/theme_clinego.R"))
OUT <- Sys.getenv("FIG_OUT", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/figures"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Recomputed per panel from each replicate's own truth table (chr:pos key), NOT
# snp_sets_summary.tsv -- that file is overwritten by whichever --sets run wrote it
# last and therefore covers different replicate counts per panel.
SS <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset09/panel_pr_recomputed.tsv"),
            colClasses = c(seed = "character"))
SS[, `:=`(n_snps = n, written = TRUE,
          n_background = n - n_causal - n_linked)]
S  <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset09/phase1_seed_medians_solo.tsv"),
            colClasses = c(seed = "character"))
map <- c(truth = "true QTNs (oracle)", intersect3 = "all 3 methods agree",
         best = "≥2 methods agree", solo_emmax = "EMMAX only",
         solo_rda = "RDA only", solo_lfmm = "LFMM only", union = "any method (union)",
         all = "ALL SNPs (no selection)", neutral_all = "neutral SNPs",
         rand_best1 = "random, size-matched")
rmap <- c(adaptive = "true QTNs (oracle)", gea_strict = "all 3 methods agree",
          gea_best = "≥2 methods agree", gea_emmax_only = "EMMAX only",
          gea_rda_only = "RDA only", gea_lfmm_only = "LFMM only", gea_union = "any method (union)",
          all = "ALL SNPs (no selection)", neutral = "neutral SNPs",
          random_matched = "random, size-matched")
SS[, panel := map[set]]; SS <- SS[!is.na(panel) & written == TRUE]
S[,  panel := rmap[marker_set]]; S <- S[!is.na(panel)]


A <- S[, .(accuracy = -median(tau)), by = .(seed, panel)]
D <- merge(SS[, .(seed, panel, n_snps, precision, recall, n_causal, n_linked, n_background)],
           A, by = c("seed", "panel"))
M <- D[, .(n_snps = median(n_snps), precision = median(precision), recall = median(recall),
           accuracy = median(accuracy)), by = panel][order(-accuracy)]
M[, panel := factor(panel, levels = M$panel)]
fwrite(M, file.path(OUT, "panel_summary.tsv"), sep = "\t")

# ---- A: which set predicts best -------------------------------------------
oracle <- M[panel == "true QTNs (oracle)", accuracy]
pA <- ggplot(M, aes(panel, accuracy,
                    fill = ifelse(panel == "true QTNs (oracle)", "oracle", "obtainable"))) +
    geom_col(width = 0.68) +
    geom_hline(yintercept = oracle, linetype = "dashed", colour = CLINEGO_THRESHOLD, linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", accuracy)), vjust = -0.4, size = 3, colour = CLINEGO_COL$fg) +
    scale_fill_manual(values = c(oracle = CLINEGO_THRESHOLD, obtainable = CLINEGO_RETAINED)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = "A. How well does each marker set predict where plants grow badly?",
         subtitle = "Higher is better. Dashed line = the true causal loci (an oracle you cannot have on real data).",
         x = NULL, y = "prediction accuracy", fill = NULL) +
    theme_clinego() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(OUT, "explain_A_accuracy.png"), pA, width = 9, height = 5, dpi = 300, bg = "white")

# ---- B: more loci is worse, not better -------------------------------------
pB <- ggplot(M, aes(recall, accuracy)) +
    geom_smooth(method = "lm", se = TRUE, colour = CLINEGO_THRESHOLD, linewidth = 0.5, alpha = 0.12) +
    geom_point(size = 3.2, colour = CLINEGO_RETAINED) +
    ggrepel::geom_text_repel(aes(label = panel), size = 3, colour = CLINEGO_COL$fg, seed = 1) +
    labs(title = "B. Catching MORE of the true loci makes the prediction WORSE",
         subtitle = "Each dot is one marker set. Right = finds more true loci. Up = predicts better.",
         x = "fraction of true causal loci captured (recall)", y = "prediction accuracy") +
    theme_clinego()
ggsave(file.path(OUT, "explain_B_recall.png"), pB, width = 8, height = 5, dpi = 300, bg = "white")

# ---- C: what is actually in each set ---------------------------------------
CM <- D[, .(causal = median(n_causal), linked = median(n_linked),
            background = median(n_background)), by = panel]
CL <- melt(CM, id.vars = "panel", variable.name = "category", value.name = "n")
CL[, panel := factor(panel, levels = rev(M$panel))]
CL[, category := factor(category, levels = c("causal", "linked", "background"),
                        labels = c("causal (a real adaptive locus)",
                                   "linked (near a real one)",
                                   "background (noise)"))]
pC <- ggplot(CL, aes(panel, n, fill = category)) +
    geom_col(position = "fill", width = 0.7) + coord_flip() +
    scale_fill_manual(values = c(CLINEGO_RETAINED, CLINEGO_CATEGORICAL[1], CLINEGO_REMOVED)) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "C. What each marker set is made of",
         subtitle = "Sets that reach for more loci pay for it in background noise.",
         x = NULL, y = "share of the marker set", fill = NULL) +
    theme_clinego() + guides(fill = guide_legend(nrow = 1))
ggsave(file.path(OUT, "explain_C_composition.png"), pC, width = 9, height = 5, dpi = 300, bg = "white")

print(M)
message("\nwritten to ", OUT)
