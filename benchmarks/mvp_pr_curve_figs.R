#!/usr/bin/env Rscript
# =============================================================================
# mvp_pr_curve_figs.R -- the precision/recall trade-off as a LINE plot.
#   x = recall, y = precision. No rank strip, no heatmap.
#
# Each selection rule is one fixed operating point, so the line does not come
# from a threshold sweep -- it joins the rules in recall order, tracing the
# trade-off the rule set spans. Three layouts:
#   B5a  faceted by architecture, rules labelled on the curve
#   B5b  one panel, one line per architecture
#   B5c  as B5a with the 30 per-replicate points behind the median curve
#
# Positives = causal only. With causal+linked every rule is 77-100% precise
# (the linked class is huge), the precision axis collapses, and the trade-off
# is invisible -- see B3_pr_space_usable. Both are emitted; causal is the
# informative one.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})

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
ARCH_COL <- c("oligogenic"=MINOU[["teal"]], "moderately polygenic"=MINOU[["amber"]],
              "highly polygenic"=MINOU[["navy"]])
ARCHL <- c(oliogenic="oligogenic", `mod-polygenic`="moderately polygenic",
           `highly-polygenic`="highly polygenic")
S <- c(best="2/3 methods", intersect3="3/3 methods", union="1/3 methods",
       solo_rda="RDA", solo_emmax="EMMAX", solo_lfmm="LFMM")

man <- fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character"))[arm=="primary"]
D   <- merge(fread(file.path(EVAL,OFF,"panel_pr_recomputed.tsv"), colClasses=c(seed="character")),
             man[, .(seed, arch = factor(ARCHL[arch_level], levels=ARCHL))], by="seed")
den <- D[set=="all", .(seed, C=n_causal, L=n_linked)]
M   <- merge(D[set %in% names(S)], den, by="seed")
M[, rule := factor(S[set], levels=names(RULE_COL))]

go <- function(mode, sfx) {
  X <- copy(M)
  if (mode=="usable") X[, `:=`(tp=n_causal+n_linked, POS=C+L)] else X[, `:=`(tp=n_causal, POS=C)]
  X <- X[n>0 & POS>0][, `:=`(precision = tp/n, recall = tp/POS)]
  PR <- X[, .(precision = median(precision), recall = median(recall)), by=.(rule, arch)]
  setorder(PR, arch, recall)
  fwrite(PR, file.path(OUT, sprintf("B5_pr_curve_%s.tsv", mode)), sep="\t")

  base <- function(p) p +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "recall  (of the causal loci present, how many were found)",
           y = "precision  (of the markers returned, how many are causal)") +
      theme_clinego() +
      theme(plot.tag = element_text(face="bold", size=13),
            strip.text = element_text(size=8.5), panel.border = element_blank())

  # ---- B5a: faceted, rules labelled along the curve -----------------------
  p_a <- base(ggplot(PR, aes(recall, precision)) +
      geom_path(colour = MINOU[["grey"]], linewidth = 0.7) +
      geom_point(aes(colour = rule), size = 3.4) +
      geom_text_repel(aes(label = rule, colour = rule), size = 2.8, seed = 1,
                      show.legend = FALSE, max.overlaps = 20, box.padding = 0.45) +
      facet_wrap(~ arch, nrow = 1, scales = "free") +
      scale_colour_manual(values = RULE_COL, guide = "none") +
      labs(tag = paste0("B5a", sfx)))
  clinego_save_both(file.path(OUT, sprintf("B5a_pr_curve_faceted_%s", mode)), p_a, w = 11.5, h = 4.4)

  # ---- B5b: one panel, one line per architecture --------------------------
  p_b <- base(ggplot(PR, aes(recall, precision, colour = arch)) +
      geom_path(linewidth = 0.9) +
      geom_point(size = 3.2) +
      geom_text_repel(aes(label = rule), size = 2.7, seed = 1, show.legend = FALSE,
                      max.overlaps = 24, box.padding = 0.4) +
      scale_colour_manual(values = ARCH_COL, name = NULL) +
      labs(tag = paste0("B5b", sfx))) +
      theme(legend.position = "bottom")
  clinego_save_both(file.path(OUT, sprintf("B5b_pr_curve_single_%s", mode)), p_b, w = 7.6, h = 5.6)

  # ---- B5c: per-replicate spread behind the median curve ------------------
  p_c <- base(ggplot() +
      geom_point(data = X, aes(recall, precision, colour = rule), alpha = 0.16, size = 1) +
      geom_path(data = PR, aes(recall, precision), colour = MINOU[["grey"]], linewidth = 0.7) +
      geom_point(data = PR, aes(recall, precision, colour = rule), size = 3.2) +
      geom_text_repel(data = PR, aes(recall, precision, label = rule, colour = rule),
                      size = 2.7, seed = 1, show.legend = FALSE, max.overlaps = 20,
                      box.padding = 0.45) +
      facet_wrap(~ arch, nrow = 1, scales = "free") +
      scale_colour_manual(values = RULE_COL, guide = "none") +
      labs(tag = paste0("B5c", sfx)))
  clinego_save_both(file.path(OUT, sprintf("B5c_pr_curve_spread_%s", mode)), p_c, w = 11.5, h = 4.4)

  message(sprintf("  OK B5a/b/c (%s)", mode)); print(PR)
}
go("causal", "")
go("usable", "′")
