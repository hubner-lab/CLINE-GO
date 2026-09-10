#!/usr/bin/env Rscript
# =============================================================================
# mvp_beta_sweep_figs.R -- "which selection rule wins depends on how you weight
# precision against recall". Four established ways of showing that, so one can
# be chosen on how it reads rather than on what was drawn first.
#
# THE CLAIM. Rules differ in their true-positive / false-positive balance. No
# single detection score identifies the rule that predicts best downstream --
# precision-weighted scores pick EMMAX, recall-weighted scores pick 1/3-methods,
# and the rule that wins on genomic offset (2/3) is intermediate on all of them.
# That is why every rule was carried forward into the maladaptation modelling
# rather than pre-selected on a detection statistic.
#
# Positives = causal + linked (a marker in LD with a causal locus carries real
# signal); the causal-only variant is emitted too, since that is where the
# precision spread between rules is widest.
#
#   B1  F-beta curves + winner strip   (quantitative, shows the crossovers)
#   B2  rank bump chart               (cleanest read of "who leads when")
#   B3  precision-recall space + iso-F contours  (classic; shows WHY it flips)
#   B4  rule x beta heatmap of rank   (most compact, no lines to trace)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(cowplot)})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_main"))
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

MINOU <- c(teal="#00798c", red="#d1495b", amber="#edae49", sage="#66a182",
           navy="#2e4057", grey="#8d96a3")
sh <- function(b,n) c(b, grDevices::colorRampPalette(c(b,"white"))(6)[2:(1+(n-1)%/%2+(n-1)%%2)],
                        grDevices::colorRampPalette(c(b,"black"))(6)[2:(1+(n-1)%/%2)])[seq_len(n)]
g3 <- sh(MINOU[["sage"]],3); a3 <- sh(MINOU[["amber"]],3)
RULE_COL <- c("2/3 methods"=g3[1], "3/3 methods"=g3[2], "1/3 methods"=g3[3],
              "RDA"=a3[1], "EMMAX"=a3[2], "LFMM"=a3[3])
RULES <- names(RULE_COL)
ARCHL <- c(oliogenic="oligogenic", `mod-polygenic`="moderately polygenic",
           `highly-polygenic`="highly polygenic")
S <- c(best="2/3 methods", intersect3="3/3 methods", union="1/3 methods",
       solo_rda="RDA", solo_emmax="EMMAX", solo_lfmm="LFMM")

man <- fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character"))[arm=="primary"]
D   <- fread(file.path(EVAL,OFF,"panel_pr_recomputed.tsv"), colClasses=c(seed="character"))
D   <- merge(D, man[, .(seed, arch = factor(ARCHL[arch_level], levels=ARCHL))], by="seed")
den <- D[set=="all", .(seed, C=n_causal, L=n_linked)]
M   <- merge(D[set %in% names(S)], den, by="seed")
M[, rule := factor(S[set], levels=RULES)]

pr_for <- function(mode) {
    x <- copy(M)
    if (mode=="usable") { x[, `:=`(tp = n_causal + n_linked, POS = C + L)] }
    else                { x[, `:=`(tp = n_causal,            POS = C)] }
    x <- x[n > 0 & POS > 0]
    x[, .(precision = median(tp/n), recall = median(tp/POS)), by=.(rule, arch)]
}
BETAS <- 2^seq(log2(0.1), log2(16), length.out = 40)
fbeta <- function(p,r,b) { d <- b*b*p + r; fifelse(d>0, (1+b*b)*p*r/d, 0) }

build <- function(mode, tag_suffix) {
  PR <- pr_for(mode)
  fwrite(PR, file.path(OUT, sprintf("B_precision_recall_%s.tsv", mode)), sep="\t")
  SW <- PR[, .(beta = BETAS, F = fbeta(precision, recall, BETAS)), by=.(rule, arch)]
  SW[, rank := frank(-F, ties.method="min"), by=.(arch, beta)]
  win <- SW[rank == 1L]
  fwrite(SW, file.path(OUT, sprintf("B_beta_sweep_%s.tsv", mode)), sep="\t")

  # ---- B1: F-beta curves + winner strip ----------------------------------
  p1 <- ggplot(SW, aes(beta, F, colour = rule)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ arch, nrow = 1, scales = "free_y") +
      scale_x_log10(breaks = c(0.1,0.25,0.5,1,2,4,8,16)) +
      scale_colour_manual(values = RULE_COL, name = NULL) +
      labs(tag = paste0("B1", tag_suffix), x = NULL, y = "F-beta score") +
      theme_clinego() +
      theme(plot.tag = element_text(face="bold", size=13), legend.position="none",
            strip.text = element_text(size=8.5), panel.border = element_blank(),
            axis.text.x = element_blank(), axis.ticks.x = element_blank())
  strip <- ggplot(win, aes(beta, y = 1, fill = rule)) +
      geom_tile(height = 1) +
      facet_wrap(~ arch, nrow = 1) +
      scale_x_log10(breaks = c(0.1,0.25,0.5,1,2,4,8,16)) +
      scale_fill_manual(values = RULE_COL, name = NULL) +
      labs(x = "β        (← precision weighted            recall weighted →)", y = NULL) +
      theme_clinego() +
      theme(legend.position = "bottom", strip.text = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.border = element_blank(), strip.background = element_blank()) +
      guides(fill = guide_legend(nrow = 1))
  clinego_save_both(file.path(OUT, sprintf("B1_fbeta_curves_%s", mode)),
                  plot_grid(p1, strip, ncol=1, rel_heights=c(1, 0.42), align="v", axis="lr"),
                  w = 11, h = 5)

  # ---- B2: rank bump chart -----------------------------------------------
  p2 <- ggplot(SW, aes(beta, rank, colour = rule)) +
      geom_line(linewidth = 1.1) +
      geom_point(data = SW[beta %in% range(BETAS)], size = 2) +
      facet_wrap(~ arch, nrow = 1) +
      scale_x_log10(breaks = c(0.1,0.5,1,4,16)) +
      scale_y_reverse(breaks = 1:6) +
      scale_colour_manual(values = RULE_COL, name = NULL) +
      labs(tag = paste0("B2", tag_suffix),
           x = "β        (← precision weighted            recall weighted →)",
           y = "rank (1 = best)") +
      theme_clinego() +
      theme(plot.tag = element_text(face="bold", size=13), legend.position="bottom",
            strip.text = element_text(size=8.5), panel.border = element_blank()) +
      guides(colour = guide_legend(nrow = 1))
  clinego_save_both(file.path(OUT, sprintf("B2_rank_bump_%s", mode)), p2, w = 11, h = 4.6)

  # ---- B3: precision-recall space with iso-F1 contours --------------------
  gr <- CJ(recall = seq(max(1e-4, min(PR$recall)*0.7), max(PR$recall)*1.25, length.out = 220),
           precision = seq(max(0, min(PR$precision)-0.05), 1, length.out = 220))
  gr[, F1 := fbeta(precision, recall, 1)]
  p3 <- ggplot() +
      geom_contour(data = gr, aes(recall, precision, z = F1),
                   colour = MINOU[["grey"]], linewidth = 0.3, bins = 12) +
      geom_point(data = PR, aes(recall, precision, colour = rule), size = 3.2) +
      ggrepel::geom_text_repel(data = PR, aes(recall, precision, label = rule, colour = rule),
                               size = 2.7, show.legend = FALSE, seed = 1, max.overlaps = 20) +
      facet_wrap(~ arch, nrow = 1, scales = "free") +
      scale_x_log10() +
      scale_colour_manual(values = RULE_COL, guide = "none") +
      labs(tag = paste0("B3", tag_suffix),
           x = "recall (log scale)", y = "precision") +
      theme_clinego() +
      theme(plot.tag = element_text(face="bold", size=13),
            strip.text = element_text(size=8.5), panel.border = element_blank())
  clinego_save_both(file.path(OUT, sprintf("B3_pr_space_%s", mode)), p3, w = 11.5, h = 4.6)

  # ---- B4: rank heatmap ---------------------------------------------------
  HB <- c(0.1,0.15,0.25,0.4,0.6,1,1.6,2.5,4,6.5,10,16)
  HM <- PR[, .(beta = HB, F = fbeta(precision, recall, HB)), by=.(rule, arch)]
  HM[, rank := frank(-F, ties.method="min"), by=.(arch, beta)]
  p4 <- ggplot(HM, aes(factor(round(beta,2)), rule, fill = factor(rank))) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = rank), size = 2.6, colour = "white") +
      facet_wrap(~ arch, nrow = 1) +
      scale_fill_manual(values = colorRampPalette(c(MINOU[["sage"]], MINOU[["grey"]],
                                                    MINOU[["red"]]))(6), name = "rank") +
      labs(tag = paste0("B4", tag_suffix),
           x = "β        (← precision weighted            recall weighted →)", y = NULL) +
      theme_clinego() +
      theme(plot.tag = element_text(face="bold", size=13),
            axis.text.x = element_text(size=6, angle=90, vjust=0.5),
            strip.text = element_text(size=8.5), panel.border = element_blank(),
            legend.position = "none")
  clinego_save_both(file.path(OUT, sprintf("B4_rank_heatmap_%s", mode)), p4, w = 11, h = 4)

  message(sprintf("  OK B1-B4 (%s)", mode))
  print(dcast(win[, .(n = .N), by=.(arch, rule)], rule ~ arch, value.var="n", fill=0L))
}
suppressPackageStartupMessages(library(ggrepel))
build("usable", "")        # positives = causal + linked  (as requested)
build("causal", "′")       # positives = causal only      (widest precision spread)
message("\nwritten to ", OUT)
