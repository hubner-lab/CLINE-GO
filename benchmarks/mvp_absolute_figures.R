#!/usr/bin/env Rscript
# mvp_absolute_figures.R -- absolute-count figures for the marker-panel comparison.
#
# WHY THESE EXIST. The composition figure (explain_C) is normalised to 100%, which
# makes `EMMAX only` look excellent: 20% of its markers are causal, the highest of
# any obtainable panel. But that is 20% OF 12 MARKERS. Normalising hides the thing
# that actually decides whether a panel works -- how many usable markers it puts in
# front of the offset model, and whether that number moves the right way as the
# trait gets more polygenic.
#
#   D  absolute make-up of each panel (counts, not shares)
#   E  the failure mode: panel size vs architecture, next to accuracy vs architecture
#   F  usable markers vs accuracy, one point per panel x architecture
#
# All three are per-architecture, because the pooled view is what created the
# misleading impression in the first place.

suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/mvp_arm.R"))
# Which scored offset root to read. Defaults to offset09, the 32-replicate run, so an
# unparameterised call reproduces the frozen figure set; a new cohort passes OFFSET_DIR.
OFF <- Sys.getenv("OFFSET_DIR", "offset09")
OUT <- Sys.getenv("FIG_OUT", file.path(PIPELINE_ROOT,
        "benchmarks/mvp_eval/figures_FINAL_32seeds_112gardens"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

D  <- fread(file.path(PIPELINE_ROOT, file.path("benchmarks/mvp_eval", OFF, "panel_pr_recomputed.tsv")),
            colClasses = c(seed = "character"))
D[, n_background := n - n_causal - n_linked]
S  <- fread(file.path(PIPELINE_ROOT, file.path("benchmarks/mvp_eval", OFF, "phase1_seed_medians_solo.tsv")),
            colClasses = c(seed = "character"))
man <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))

RL <- c(adaptive = "truth", gea_strict = "intersect3", gea_best = "best",
        gea_emmax_only = "solo_emmax", gea_rda_only = "solo_rda",
        gea_lfmm_only = "solo_lfmm", gea_union = "union",
        all = "all", neutral = "neutral_all", random_matched = "rand_best1")
S[, set := RL[marker_set]]
A <- S[!is.na(set), .(accuracy = -median(tau)), by = .(seed, set)]
M <- merge(D, A, by = c("seed", "set"))
M <- merge(M, mvp_prim(man)[, .(seed, arch_level, arm)], by = "seed")

LAB <- c(truth = "true QTNs\n(oracle)", best = "≥2 methods\nagree",
         intersect3 = "all 3\nagree", union = "union\n(any method)",
         solo_lfmm = "LFMM\nonly", solo_rda = "RDA\nonly", solo_emmax = "EMMAX\nonly",
         rand_best1 = "random\n(size-matched)")
M <- M[set %in% names(LAB)]
M[, panel := factor(LAB[set], levels = LAB[c("truth","best","intersect3","union",
                                             "solo_rda","solo_lfmm","solo_emmax","rand_best1")])]
M[, arch := factor(arch_level, levels = c("oliogenic","mod-polygenic","highly-polygenic"),
                   labels = c("oligogenic\n(~7 causal loci)","mod-polygenic\n(~42 causal loci)",
                              "highly-polygenic\n(~410 causal loci)"))]

# ---- D: absolute composition, per architecture -----------------------------
# Only the OBTAINABLE panels: the oracle and the random floor are references, not
# candidates, and including them compresses the axis and invites the wrong comparison.
# Ordered by usable markers (causal + linked) pooled over architectures, so the same
# order holds in every facet and panels stay comparable across them.
DM <- M[!panel %in% c("true QTNs\n(oracle)", "random\n(size-matched)")]
DM[, panel := droplevels(panel)]
ordD <- DM[, .(u = as.numeric(median(n_causal + n_linked))), by = panel][order(u)]$panel
DM[, panel := factor(as.character(panel), levels = as.character(ordD))]
CM <- DM[, .(causal = as.numeric(median(n_causal)), linked = as.numeric(median(n_linked)),
            background = as.numeric(median(n_background)), total = as.numeric(median(n))), by = .(panel, arch)]
CL <- melt(CM, id.vars = c("panel","arch"), measure.vars = c("causal","linked","background"),
           variable.name = "category", value.name = "n")
CL[, category := factor(category, levels = c("causal","linked","background"),
                        labels = c("causal (a real adaptive locus)",
                                   "linked (next to a real one)",
                                   "background (noise)"))]
pD <- ggplot(CL, aes(panel, n, fill = category)) +
    geom_col(width = 0.7) +
    geom_text(data = CM, aes(panel, total, label = as.integer(causal + linked)),
              inherit.aes = FALSE, vjust = -0.35, size = 2.7, colour = CLINEGO_COL$fg) +
    facet_wrap(~ arch, nrow = 1) +
    coord_flip() +
    scale_fill_manual(values = c(CLINEGO_RETAINED, CLINEGO_CATEGORICAL[1], CLINEGO_REMOVED)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "D. How many markers each selection rule actually delivers",
         subtitle = paste("ABSOLUTE counts, obtainable panels only. Number = usable markers (causal + linked).",
                          "\nEMMAX's high purity is 20% of ~12 markers; >=2-agree is 6% of ~150."),
         x = NULL, y = "markers (median per replicate)", fill = NULL) +
    theme_clinego() +
    theme(axis.text.y = element_text(size = 7.5), legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1))
ggsave(file.path(OUT, "abs_D_composition_counts.png"), pD, width = 12, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "abs_D_composition_counts.svg"), pD, width = 12, height = 5.5, bg = "white")
# Duplicated under a MAIN_ name: this is the headline panel for GEA marker-selection
# efficiency, and the filename should say so without needing the README.
ggsave(file.path(OUT, "MAIN_FIGURE_gea_marker_efficiency.png"), pD, width = 12, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "MAIN_FIGURE_gea_marker_efficiency.svg"), pD, width = 12, height = 5.5, bg = "white")
fwrite(CM, file.path(OUT, "abs_D_composition_counts.tsv"), sep = "\t")

# ---- E: the failure mode ---------------------------------------------------
# Left: how panel size responds to polygenicity. Right: accuracy.
# A usable rule must NOT shrink as the trait gets more polygenic.
KEEP <- c("≥2 methods\nagree","all 3\nagree","EMMAX\nonly","RDA\nonly",
          "LFMM\nonly","union\n(any method)")
E <- M[panel %in% KEEP, .(markers = median(n), accuracy = median(accuracy)), by = .(panel, arch)]
EL <- melt(E, id.vars = c("panel","arch"), variable.name = "what", value.name = "v")
EL[, what := factor(what, levels = c("markers","accuracy"),
                    labels = c("markers in the set (median)", "prediction accuracy"))]
pE <- ggplot(EL, aes(arch, v, colour = panel, group = panel)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
    facet_wrap(~ what, scales = "free_y") +
    scale_colour_manual(values = CLINEGO_CATEGORICAL[1:6]) +
    labs(title = "E. Why ≥2-agree wins: it is the only rule that grows with the trait",
         subtitle = "EMMAX shrinks 22→10→5 markers exactly as more are needed, and its accuracy collapses. ≥2-agree grows 140→150→158.",
         x = NULL, y = NULL, colour = NULL) +
    theme_clinego() +
    theme(axis.text.x = element_text(size = 7.5)) +
    guides(colour = guide_legend(nrow = 1))
ggsave(file.path(OUT, "abs_E_size_vs_architecture.png"), pE, width = 12, height = 5.5, dpi = 300, bg = "white")
fwrite(E, file.path(OUT, "abs_E_size_vs_architecture.tsv"), sep = "\t")

# ---- F: usable markers vs accuracy -----------------------------------------
# "Usable" = causal + linked. Too few and the panel cannot span the gradient;
# too many and background noise takes over.
FF <- M[, .(usable = as.numeric(median(n_causal + n_linked)), noise = as.numeric(median(n_background)),
            accuracy = median(accuracy)), by = .(panel, arch)]
pF <- ggplot(FF, aes(usable, accuracy, colour = panel, shape = arch)) +
    geom_point(size = 3.4, alpha = 0.9) +
    scale_x_log10() +
    scale_colour_manual(values = CLINEGO_CATEGORICAL[1:8]) +
    labs(title = "F. Too few usable markers is as bad as too many noisy ones",
         subtitle = "x = causal + linked markers in the set (log scale). Each point is one panel in one architecture.",
         x = "usable markers (causal + linked)", y = "prediction accuracy",
         colour = NULL, shape = NULL) +
    theme_clinego() + guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(file.path(OUT, "abs_F_usable_vs_accuracy.png"), pF, width = 10, height = 6, dpi = 300, bg = "white")
fwrite(FF, file.path(OUT, "abs_F_usable_vs_accuracy.tsv"), sep = "\t")

# ---- the table the figures rest on -----------------------------------------
TB <- M[, .(markers = as.integer(median(n)),
            causal = as.integer(median(n_causal)),
            linked = as.integer(median(n_linked)),
            noise = as.integer(median(n_background)),
            usable = as.integer(median(n_causal + n_linked)),
            accuracy = round(median(accuracy), 3)), by = .(arch, panel)][order(arch, -accuracy)]
fwrite(TB, file.path(OUT, "abs_table_by_architecture.tsv"), sep = "\t")
print(TB)
message("\nwritten to ", OUT)
