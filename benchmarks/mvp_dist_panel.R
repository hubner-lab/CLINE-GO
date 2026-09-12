#!/usr/bin/env Rscript
# =============================================================================
# mvp_dist_panel.R -- accuracy across architecture, showing the DISTRIBUTION.
#
# WHY. A slopegraph of medians cannot support the robustness claim (the claim is
# about spread), but whiskers layered on crossing lines are unreadable. The fix
# is to give the distribution its own grid: one small multiple per selection
# rule, three architectures inside each, all 30 replicates visible.
#
# Every panel carries the causal-loci reference as a dashed green line at that
# architecture's value, so each rule is read against the ceiling in place.
#
#   G1  boxplot + all replicate points     (standard, shows every observation)
#   G2  ggdist half-eye                    (density slab + median/interval)
#   G3  transposed: architecture as facet, rules on y
#   G4  sina                               (points positioned by density)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)
                                library(ggdist); library(ggforce)})
ROOT <- Sys.getenv("PIPELINE_ROOT","/pipeline"); EVAL <- file.path(ROOT,"benchmarks/mvp_eval")
OFF <- Sys.getenv("OFFSET_DIR","offset11"); METHOD <- Sys.getenv("METHOD","GFoffset")
OUT <- Sys.getenv("FIG_OUT", file.path(EVAL,"figures_main"))
source(file.path(ROOT,"scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))

MINOU <- c(teal="#00798c", red="#d1495b", amber="#edae49", sage="#66a182",
           navy="#2e4057", grey="#8d96a3")
sh <- function(b,n) c(b, grDevices::colorRampPalette(c(b,"white"))(6)[2:(1+(n-1)%/%2+(n-1)%%2)],
                        grDevices::colorRampPalette(c(b,"black"))(6)[2:(1+(n-1)%/%2)])[seq_len(n)]
g3 <- sh(MINOU[["sage"]],3); a3 <- sh(MINOU[["amber"]],3)
REF_GOOD <- MINOU[["sage"]]; REF_BAD <- MINOU[["red"]]
RULE_COL <- c("2/3 methods"=g3[1],"3/3 methods"=g3[2],"1/3 methods"=g3[3],
              "RDA"=a3[1],"EMMAX"=a3[2],"LFMM"=a3[3],
              "causal loci"=REF_GOOD,"all loci"=REF_BAD)
ORD <- c("2/3 methods","3/3 methods","1/3 methods","RDA","EMMAX","LFMM",
         "causal loci","all loci")
SETLAB <- c(gea_best="2/3 methods", gea_strict="3/3 methods", gea_union="1/3 methods",
            gea_rda_only="RDA", gea_emmax_only="EMMAX", gea_lfmm_only="LFMM",
            adaptive="causal loci", all="all loci")
ARCHL <- c(oliogenic="oligo",`mod-polygenic`="moderate",`highly-polygenic`="high")

man <- mvp_prim(fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character")))
S <- fread(file.path(EVAL,OFF,"phase1_seed_medians_solo.tsv"), colClasses=c(seed="character"))
S <- S[seed %in% man$seed & marker_set %in% names(SETLAB) & method_label==METHOD]
S[, `:=`(accuracy=-tau, label=factor(SETLAB[marker_set], levels=ORD))]
S <- merge(S, man[, .(seed, arch=factor(ARCHL[arch_level], levels=ARCHL))], by="seed")
fwrite(S[, .(seed, label, arch, accuracy)], file.path(OUT,"G_distributions.tsv"), sep="\t")

# causal-loci ceiling per architecture, replicated into every facet
CEIL <- S[label=="causal loci", .(ceil=median(accuracy)), by=arch]
GRID <- CJ(label=factor(ORD, levels=ORD), arch=factor(levels(S$arch), levels=levels(S$arch)))
CEIL <- merge(GRID, CEIL, by="arch")

base <- function(p,tag) p + theme_clinego() +
    theme(plot.tag=element_text(face="bold",size=13), panel.border=element_blank(),
          strip.text=element_text(size=8), legend.position="none",
          axis.text.x=element_text(size=7.5)) + labs(tag=tag)

# ---- G1: boxplot + every replicate -----------------------------------------
p1 <- base(ggplot(S, aes(arch, accuracy)) +
    geom_hline(data=CEIL, aes(yintercept=ceil), linetype="22",
               colour=REF_GOOD, linewidth=0.45) +
    geom_boxplot(aes(fill=label), outlier.shape=NA, alpha=0.55, linewidth=0.3, width=0.62) +
    geom_jitter(width=0.14, size=0.5, alpha=0.4, colour=CLINEGO_COL$fg) +
    facet_wrap(~ label, nrow=2) +
    scale_fill_manual(values=RULE_COL) +
    labs(x=NULL, y="Accuracy"), "C")
clinego_save_both(file.path(OUT,"G1_dist_box_points"), p1, w=10.5, h=5.6)

# ---- G2: ggdist half-eye ----------------------------------------------------
p2 <- base(ggplot(S, aes(arch, accuracy, fill=label)) +
    geom_hline(data=CEIL, aes(yintercept=ceil), linetype="22",
               colour=REF_GOOD, linewidth=0.45) +
    stat_halfeye(adjust=0.9, width=0.65, .width=c(0.5,0.95), point_size=1.4,
                 slab_alpha=0.6, interval_size_range=c(0.3,0.8)) +
    facet_wrap(~ label, nrow=2) +
    scale_fill_manual(values=RULE_COL) +
    labs(x=NULL, y="Accuracy"), "C")
clinego_save_both(file.path(OUT,"G2_dist_halfeye"), p2, w=10.5, h=5.6)

# ---- G3: transposed -- architecture as facet, rules on y -------------------
CE3 <- S[label=="causal loci", .(ceil=median(accuracy)), by=arch]
p3 <- base(ggplot(S, aes(accuracy, label, fill=label)) +
    geom_vline(data=CE3, aes(xintercept=ceil), linetype="22",
               colour=REF_GOOD, linewidth=0.45) +
    geom_boxplot(outlier.shape=NA, alpha=0.55, linewidth=0.3, width=0.6) +
    geom_jitter(height=0.14, size=0.45, alpha=0.35, colour=CLINEGO_COL$fg) +
    facet_wrap(~ arch, nrow=1) +
    scale_y_discrete(limits=rev(ORD)) +
    scale_fill_manual(values=RULE_COL) +
    labs(x="Accuracy", y=NULL), "C") +
    theme(axis.text.y=element_text(size=8))
clinego_save_both(file.path(OUT,"G3_dist_transposed"), p3, w=11, h=4.6)

# ---- G4: sina ---------------------------------------------------------------
p4 <- base(ggplot(S, aes(arch, accuracy, colour=label)) +
    geom_hline(data=CEIL, aes(yintercept=ceil), linetype="22",
               colour=REF_GOOD, linewidth=0.45) +
    geom_sina(size=0.8, alpha=0.55, maxwidth=0.7) +
    stat_summary(fun=median, geom="crossbar", width=0.5, linewidth=0.3,
                 colour=CLINEGO_COL$fg) +
    facet_wrap(~ label, nrow=2) +
    scale_colour_manual(values=RULE_COL) +
    labs(x=NULL, y="Accuracy"), "C")
clinego_save_both(file.path(OUT,"G4_dist_sina"), p4, w=10.5, h=5.6)
message("  OK G1-G4")
