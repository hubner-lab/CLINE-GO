#!/usr/bin/env Rscript
# =============================================================================
# mvp_regime_panel.R -- G3 (transposed distribution grid) with every replicate
# point coloured by SELECTION REGIME, to expose the bimodality directly.
#
# Regime is read off the corpus's own `architecture` label:
#   1-trait     only one trait under selection -- no axis asymmetry possible
#   equal-S     two traits, equal selection strength on both axes
#   unequal-S   two traits, one axis under stronger selection than the other
#
# The low cloud in the highly-polygenic panels is exactly the unequal-S set
# (9/9), which is why the points are coloured rather than the boxes.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- Sys.getenv("PIPELINE_ROOT","/pipeline"); EVAL <- file.path(ROOT,"benchmarks/mvp_eval")
OFF <- Sys.getenv("OFFSET_DIR","offset11"); METHOD <- Sys.getenv("METHOD","GFoffset")
OUT <- Sys.getenv("FIG_OUT", file.path(EVAL,"figures_main"))
source(file.path(ROOT,"scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))

MINOU <- c(teal="#00798c", red="#d1495b", amber="#edae49", sage="#66a182",
           navy="#2e4057", grey="#8d96a3")
REG_COL <- c("equal-S"=MINOU[["teal"]], "unequal-S"=MINOU[["red"]],
             "1-trait"=MINOU[["amber"]])
ORD <- c("2/3 methods","3/3 methods","1/3 methods","RDA","EMMAX","LFMM",
         "causal loci","all loci")
SETLAB <- c(gea_best="2/3 methods", gea_strict="3/3 methods", gea_union="1/3 methods",
            gea_rda_only="RDA", gea_emmax_only="EMMAX", gea_lfmm_only="LFMM",
            adaptive="causal loci", all="all loci")
ARCHL <- c(oliogenic="oligogenic",`mod-polygenic`="moderately polygenic",
           `highly-polygenic`="highly polygenic")

man <- mvp_prim(fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character")))
man[, regime := fifelse(grepl("1-trait", architecture), "1-trait",
                fifelse(grepl("unequal-S", architecture), "unequal-S", "equal-S"))]
man[, regime := factor(regime, levels=names(REG_COL))]
S <- fread(file.path(EVAL,OFF,"phase1_seed_medians_solo.tsv"), colClasses=c(seed="character"))
S <- S[seed %in% man$seed & marker_set %in% names(SETLAB) & method_label==METHOD]
S[, `:=`(accuracy=-tau, label=factor(SETLAB[marker_set], levels=ORD))]
S <- merge(S, man[, .(seed, regime, arch=factor(ARCHL[arch_level], levels=ARCHL))], by="seed")
fwrite(S[, .(seed, label, arch, regime, accuracy)],
       file.path(OUT,"G3b_regime_points.tsv"), sep="\t")
print(S[, .(n=uniqueN(seed)), by=.(arch, regime)][order(arch, regime)])

p <- ggplot(S, aes(accuracy, label)) +
    geom_boxplot(outlier.shape=NA, fill=NA, colour=MINOU[["grey"]],
                 linewidth=0.35, width=0.62) +
    geom_jitter(aes(colour=regime), height=0.16, size=1.1, alpha=0.8) +
    facet_wrap(~ arch, nrow=1) +
    scale_y_discrete(limits=rev(ORD)) +
    scale_colour_manual(values=REG_COL, name="selection regime") +
    labs(tag="S", x="Accuracy", y=NULL) +
    theme_clinego() +
    theme(plot.tag=element_text(face="bold",size=13), panel.border=element_blank(),
          strip.text=element_text(size=8.5), legend.position="bottom",
          axis.text.y=element_text(size=8)) +
    guides(colour=guide_legend(nrow=1, override.aes=list(size=2.4, alpha=1)))
clinego_save_both(file.path(OUT,"G3b_dist_by_regime"), p, w=11, h=4.8)
message("  OK G3b_dist_by_regime")
