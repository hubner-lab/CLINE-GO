#!/usr/bin/env Rscript
# =============================================================================
# mvp_arch_panel.R -- marker set x genetic architecture, GFoffset only.
#
# THE MESSAGE. Selection rules are not uniformly good or bad: they trade places
# as the trait becomes polygenic. EMMAX is competitive on oligogenic traits and
# collapses on polygenic ones; 3/3-agreement does the same; 2/3-agreement barely
# moves. Architecture therefore has to be a first-class axis, not a facet buried
# behind a marker-set axis.
#
# Five structures, same numbers:
#   V1 dot plot     rules on y, one dot per architecture -- spread WITHIN a row
#                   is instability, which is the quantity being claimed
#   V2 small multiples  one ranked whisker panel per architecture
#   V3 slopegraph   architecture on x, one line per rule (rank crossings)
#   V4 heatmap      rule x architecture matrix with values
#   V5 deviation    gap to the causal-loci reference, diverging from zero
#
# Architecture uses a light->dark sequential ramp so "more polygenic" reads as a
# gradient rather than as three unrelated categories.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_main"))
METHOD <- Sys.getenv("METHOD", "GFoffset")
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

MINOU <- c(teal="#00798c", red="#d1495b", amber="#edae49", sage="#66a182",
           navy="#2e4057", grey="#8d96a3")
sh <- function(b,n) c(b, grDevices::colorRampPalette(c(b,"white"))(6)[2:(1+(n-1)%/%2+(n-1)%%2)],
                        grDevices::colorRampPalette(c(b,"black"))(6)[2:(1+(n-1)%/%2)])[seq_len(n)]
g3 <- sh(MINOU[["sage"]],3); a3 <- sh(MINOU[["amber"]],3); gy <- sh(MINOU[["grey"]],3)
QTN <- grDevices::colorRampPalette(c(MINOU[["navy"]],"white"))(10)[5]
REF_GOOD <- MINOU[["sage"]]; REF_BAD <- MINOU[["red"]]
RULE_COL <- c("2/3 methods"=g3[1],"3/3 methods"=g3[2],"1/3 methods"=g3[3],
              "RDA"=a3[1],"EMMAX"=a3[2],"LFMM"=a3[3],
              "causal loci"=QTN,"all loci"=gy[3])
# sequential: pale = few causal loci, dark = many
ARCH_COL <- setNames(grDevices::colorRampPalette(c("#9ecfd6", MINOU[["teal"]],
                                                   MINOU[["navy"]]))(3),
                     c("oligogenic","moderately polygenic","highly polygenic"))
SETLAB <- c(gea_best="2/3 methods", gea_strict="3/3 methods", gea_union="1/3 methods",
            gea_rda_only="RDA", gea_emmax_only="EMMAX", gea_lfmm_only="LFMM",
            adaptive="causal loci", all="all loci")
OBT <- c("2/3 methods","3/3 methods","1/3 methods","RDA","EMMAX","LFMM")
ARCHL <- c(oliogenic="oligogenic",`mod-polygenic`="moderately polygenic",
           `highly-polygenic`="highly polygenic")

man <- fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character"))[arm=="primary"]
S <- fread(file.path(EVAL,OFF,"phase1_seed_medians_solo.tsv"), colClasses=c(seed="character"))
S <- S[seed %in% man$seed & marker_set %in% names(SETLAB) & method_label==METHOD]
S[, `:=`(accuracy=-tau, label=SETLAB[marker_set])]
S <- merge(S, man[, .(seed, arch=factor(ARCHL[arch_level], levels=ARCHL))], by="seed")

M <- S[, .(accuracy=median(accuracy), q25=quantile(accuracy,.25),
           q75=quantile(accuracy,.75), n=.N), by=.(label, arch)]
# order rules by their WORST architecture: a rule is only as good as where it fails
ord <- M[label %in% OBT, .(worst=min(accuracy)), by=label][order(worst)]$label
ordall <- c(setdiff(unique(M$label), OBT), rev(ord))
M[, label := factor(label, levels=ordall)]
fwrite(M, file.path(OUT,"V_arch_by_set.tsv"), sep="\t")
print(dcast(M, label ~ arch, value.var="accuracy"))

base <- function(p,tag) p + theme_clinego() +
    theme(plot.tag=element_text(face="bold",size=13), panel.border=element_blank(),
          strip.text=element_text(size=8.5), legend.position="bottom") + labs(tag=tag)

# ---- V1: dot plot, spread within a row = instability -----------------------
p1 <- base(ggplot(M, aes(accuracy, label)) +
    geom_line(aes(group=label), colour=MINOU[["grey"]], linewidth=1.6, alpha=0.35) +
    geom_point(aes(colour=arch), size=3.1) +
    scale_colour_manual(values=ARCH_COL, name=NULL) +
    labs(x="Accuracy", y=NULL) +
    guides(colour=guide_legend(nrow=1)), "V1")
clinego_save_both(file.path(OUT,"V1_arch_dotplot"), p1, w=8.2, h=4.6)

# ---- V2: small multiples, ranked within each architecture ------------------
M2 <- copy(M)[, lab2 := reorder(label, accuracy), by=arch]
p2 <- base(ggplot(M2, aes(accuracy, tidytext_reorder <- reorder(paste(label, arch), accuracy))) +
    geom_errorbarh(aes(xmin=q25, xmax=q75, colour=label), height=0.25, linewidth=0.45) +
    geom_point(aes(colour=label), size=2.6) +
    facet_wrap(~ arch, nrow=1, scales="free_y") +
    scale_y_discrete(labels=function(x) sub(" (oligogenic|moderately polygenic|highly polygenic)$","",x)) +
    scale_colour_manual(values=RULE_COL, guide="none") +
    labs(x="Accuracy", y=NULL), "V2")
clinego_save_both(file.path(OUT,"V2_arch_small_multiples"), p2, w=11.5, h=4.4)

# ---- V3: slopegraph --------------------------------------------------------
lab3 <- M[arch==levels(arch)[1L]]
.w <- max(strwidth(as.character(lab3$label), units="inches", cex=2.9*.pt/par("ps")))
PAD <- max(0.16,(.w/(8.0-1.0))*2*1.25+0.05)
p3 <- base(ggplot(M, aes(as.integer(arch), accuracy, colour=label, group=label)) +
    geom_line(linewidth=0.95) + geom_point(size=2.2) +
    geom_text(data=lab3, aes(x=1-0.03, label=label), hjust=1, size=2.9) +
    scale_x_continuous(breaks=1:3, labels=levels(M$arch),
                       expand=expansion(add=c(PAD,0.06))) +
    scale_colour_manual(values=RULE_COL, guide="none") +
    labs(x=NULL, y="Accuracy"), "V3")
clinego_save_both(file.path(OUT,"V3_arch_slopegraph"), p3, w=8.0, h=5.0)

# ---- V4: heatmap -----------------------------------------------------------
p4 <- base(ggplot(M, aes(arch, label, fill=accuracy)) +
    geom_tile(colour="white", linewidth=0.6) +
    geom_text(aes(label=sprintf("%.3f",accuracy)), size=2.7, colour="white") +
    scale_fill_gradientn(colours=c(MINOU[["red"]],"grey92",MINOU[["teal"]]), name=NULL) +
    labs(x=NULL, y=NULL), "V4") + theme(legend.position="right")
clinego_save_both(file.path(OUT,"V4_arch_heatmap"), p4, w=7.6, h=4.2)

# ---- V5: deviation from the causal-loci reference --------------------------
ref <- M[label=="causal loci", .(arch, r=accuracy)]
D <- merge(M[label!="causal loci"], ref, by="arch")[, gap := accuracy - r]
p5 <- base(ggplot(D, aes(gap, label, fill=arch)) +
    geom_vline(xintercept=0, colour=MINOU[["navy"]], linetype="22", linewidth=0.6) +
    geom_col(position=position_dodge(width=0.75), width=0.68) +
    scale_fill_manual(values=ARCH_COL, name=NULL) +
    labs(x="accuracy relative to the causal loci", y=NULL) +
    guides(fill=guide_legend(nrow=1)), "V5")
clinego_save_both(file.path(OUT,"V5_arch_deviation"), p5, w=8.2, h=4.8)
message("  OK V1-V5")
