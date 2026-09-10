#!/usr/bin/env Rscript
# =============================================================================
# mvp_panel_try.R -- trial panels before committing the main figure.
#   B  method comparison, lines + IQR whiskers  (the reason GFoffset is chosen)
#   C  accuracy per marker set, GFoffset only, points + whiskers instead of bars
#      (two variants: with and without the marker sets kept on the x axis)
# Written to *_try_* names so nothing already agreed is overwritten.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_main"))
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

MINOU <- c(teal="#00798c", red="#d1495b", amber="#edae49", sage="#66a182",
           navy="#2e4057", grey="#8d96a3")
sh <- function(b,n) c(b, grDevices::colorRampPalette(c(b,"white"))(6)[2:(1+(n-1)%/%2+(n-1)%%2)],
                        grDevices::colorRampPalette(c(b,"black"))(6)[2:(1+(n-1)%/%2)])[seq_len(n)]
g3 <- sh(MINOU[["sage"]],3); a3 <- sh(MINOU[["amber"]],3); gy <- sh(MINOU[["grey"]],3)
QTN_BAR <- grDevices::colorRampPalette(c(MINOU[["navy"]],"white"))(10)[5]
REF_GOOD <- MINOU[["sage"]]; REF_BAD <- MINOU[["red"]]
METH_COL <- c("GFoffset"=MINOU[["teal"]], "LFMM2offset"=MINOU[["navy"]],
              "RDA-uncorrected"=MINOU[["grey"]], "RDA-corrected"=MINOU[["red"]])
RULE_COL <- c("2/3 methods"=g3[1],"3/3 methods"=g3[2],"1/3 methods"=g3[3],
              "RDA"=a3[1],"EMMAX"=a3[2],"LFMM"=a3[3],
              "causal loci"=QTN_BAR,"all loci"=gy[3])
OBT <- c("2/3 methods","3/3 methods","1/3 methods","RDA","EMMAX","LFMM")
ROLE <- c(setNames(rep("combined",3), OBT[1:3]), setNames(rep("single",3), OBT[4:6]),
          "causal loci"="reference","all loci"="reference")
SETLAB <- c(gea_best="2/3 methods", gea_strict="3/3 methods", gea_union="1/3 methods",
            gea_rda_only="RDA", gea_emmax_only="EMMAX", gea_lfmm_only="LFMM",
            adaptive="causal loci", all="all loci")
ARCHL <- c(oliogenic="oligogenic",`mod-polygenic`="moderately polygenic",
           `highly-polygenic`="highly polygenic")

man <- fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character"))[arm=="primary"]
S <- fread(file.path(EVAL,OFF,"phase1_seed_medians_solo.tsv"), colClasses=c(seed="character"))
S <- S[seed %in% man$seed & marker_set %in% names(SETLAB)]
S[, `:=`(accuracy = -tau, label = SETLAB[marker_set])]
S <- merge(S, man[, .(seed, arch = factor(ARCHL[arch_level], levels=ARCHL))], by="seed")

q <- function(d) list(accuracy=median(d), q25=quantile(d,.25), q75=quantile(d,.75))

# ---- B: method comparison, lines + IQR whiskers ----------------------------
# RDA-corrected is dropped here (S-RDA shows it anti-predicts, ignores the marker
# set on 6% of replicates, and inverts -- it is not a peer of the other three).
S <- S[method_label != "RDA-corrected"]
MB <- S[, q(accuracy), by=.(label, method_label)]
ordset <- MB[method_label!="RDA-corrected", .(a=median(accuracy)), by=label][order(-a)]$label
MB[, `:=`(label = factor(label, levels=ordset),
          method_label = factor(method_label,
                                levels=setdiff(names(METH_COL), "RDA-corrected")))]
# Direct labels at the left edge instead of a legend, matching panels C/D. Margin
# sized to the widest method name so no dead space opens up.
labB <- MB[label == ordset[1]]
.wB <- max(strwidth(as.character(labB$method_label), units="inches",
                    cex = 2.9 * .pt / par("ps")))
PADB <- max(0.20, (.wB / (8.6 - 1.0)) * (length(ordset) - 1) * 1.25 + 0.08)
fwrite(MB, file.path(OUT,"B_try_method_comparison.tsv"), sep="\t")
pd <- position_dodge(width=0.42)
pB <- ggplot(MB, aes(as.integer(label), accuracy, colour=method_label,
                     group=method_label)) +
    geom_line(position=pd, linewidth=0.9) +
    geom_errorbar(aes(ymin=q25, ymax=q75), position=pd, width=0.22, linewidth=0.4) +
    geom_point(position=pd, size=2.2) +
    geom_text(data=labB, aes(x=1-0.06, y=accuracy, label=method_label, colour=method_label),
              hjust=1, size=2.9, show.legend=FALSE) +
    scale_x_continuous(breaks=seq_along(ordset), labels=ordset,
                       expand=expansion(add=c(PADB, 0.12))) +
    scale_colour_manual(values=METH_COL, guide="none") +
    labs(tag="B", x=NULL, y="Accuracy") +
    theme_clinego() +
    theme(plot.tag=element_text(face="bold",size=13), panel.border=element_blank(),
          axis.text.x=element_text(angle=25,hjust=1,size=8))
clinego_save_both(file.path(OUT,"B_try_method_comparison"), pB, w=8.6, h=5.0)

# ---- C: GFoffset only, points + whiskers, dashed references ---------------
G <- S[method_label=="GFoffset"]
MC <- G[, q(accuracy), by=label]
MC[, grp := factor(c(combined="combining rules", single="single method",
                     reference="reference")[ROLE[label]],
                   levels=c("reference","single method","combining rules"))]
setorder(MC, grp, accuracy)
MC[, label := factor(as.character(label), levels=as.character(label))]
fwrite(MC, file.path(OUT,"C_try_accuracy_GFoffset.tsv"), sep="\t")
rq <- MC[label=="causal loci", accuracy]; ra <- MC[label=="all loci", accuracy]

mk <- function(dt, tag, faceted) {
    p <- ggplot(dt, aes(label, accuracy, colour=label)) +
        geom_hline(yintercept=rq, linetype="22", colour=REF_GOOD, linewidth=0.6) +
        geom_hline(yintercept=ra, linetype="22", colour=REF_BAD,  linewidth=0.6) +
        geom_errorbar(aes(ymin=q25, ymax=q75), width=0.18, linewidth=0.5) +
        geom_point(size=3) +
        geom_text(aes(y=q75, label=sprintf("%.3f", accuracy)), hjust=-0.3,
                  size=2.9, colour=CLINEGO_COL$fg) +
        coord_flip(ylim=c(min(dt$q25)-0.02, max(dt$q75)+0.05)) +
        scale_colour_manual(values=RULE_COL, guide="none") +
        labs(tag=tag, x=NULL, y="Accuracy") +
        theme_clinego() +
        theme(plot.tag=element_text(face="bold",size=13), panel.border=element_blank(),
              axis.text.y=element_text(size=8.5))
    if (faceted) p <- p + facet_grid(grp ~ ., scales="free_y", space="free_y", switch="y") +
        theme(strip.background=element_blank(), strip.placement="outside",
              strip.text.y.left=element_text(angle=90,size=8.5,face="bold"),
              panel.spacing.y=unit(0.5,"lines"))
    p
}
clinego_save_both(file.path(OUT,"C_try_whisker_grouped"), mk(MC,"C",TRUE),  w=7.8, h=5.0)
clinego_save_both(file.path(OUT,"C_try_whisker_flat"),
                mk(MC[order(accuracy)][, label:=factor(as.character(label),
                    levels=as.character(label))],"C",FALSE), w=7.6, h=4.2)
message("  OK trial panels"); print(MC[order(-accuracy), .(label, accuracy=round(accuracy,4),
        q25=round(q25,3), q75=round(q75,3))])
