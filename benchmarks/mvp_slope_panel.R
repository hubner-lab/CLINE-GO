#!/usr/bin/env Rscript
# =============================================================================
# mvp_slope_panel.R -- accuracy across genetic architecture as a slopegraph.
# Three refinements of the same structure, GFoffset only.
#
#   S1  plain      dashed green/red references, labels at BOTH ends
#   S2  focus      as S1 but the recommended rule emphasised, others muted
#   S3  whiskers   as S1 plus per-point IQR
#
# WHY LABELS AT BOTH ENDS. The lines cross -- that crossing IS the result -- so a
# reader who only has left-hand labels loses track of which line ends where.
#
# WHY THE AXIS CARRIES LOCUS COUNTS. "Polygenic" is a word; 7 -> 42 -> 410 causal
# loci is the thing that actually changes. Put the number on the axis.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11"); METHOD <- Sys.getenv("METHOD","GFoffset")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL,"figures_main"))
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
OBT <- c("2/3 methods","3/3 methods","1/3 methods","RDA","EMMAX","LFMM")
REFS <- c("causal loci","all loci")
SETLAB <- c(gea_best="2/3 methods", gea_strict="3/3 methods", gea_union="1/3 methods",
            gea_rda_only="RDA", gea_emmax_only="EMMAX", gea_lfmm_only="LFMM",
            adaptive="causal loci", all="all loci")
ARCHL <- c(oliogenic="oligogenic",`mod-polygenic`="moderately polygenic",
           `highly-polygenic`="highly polygenic")

man <- mvp_prim(fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character")))
S <- fread(file.path(EVAL,OFF,"phase1_seed_medians_solo.tsv"), colClasses=c(seed="character"))
S <- S[seed %in% man$seed & marker_set %in% names(SETLAB) & method_label==METHOD]
S[, `:=`(accuracy=-tau, label=SETLAB[marker_set])]
S <- merge(S, man[, .(seed, arch=factor(ARCHL[arch_level], levels=ARCHL))], by="seed")
M <- S[, .(accuracy=median(accuracy), q25=quantile(accuracy,.25),
           q75=quantile(accuracy,.75)), by=.(label, arch)]
M[, `:=`(x=as.integer(arch), is_ref = label %in% REFS)]
fwrite(M, file.path(OUT,"S_slope_values.tsv"), sep="\t")

# axis ticks carry the median causal-locus count for that architecture
PR <- fread(file.path(EVAL,OFF,"panel_pr_recomputed.tsv"), colClasses=c(seed="character"))
NC <- merge(PR[set=="all", .(seed, nc=n_causal)],
            man[, .(seed, arch=factor(ARCHL[arch_level], levels=ARCHL))], by="seed")[
            , .(nc=round(median(nc))), by=arch][order(arch)]
xlabs <- sprintf("%s\n(~%d causal loci)", NC$arch, NC$nc)

L <- M[x==1]; R <- M[x==3]
cx <- 2.9*.pt/par("ps")
padL <- max(0.18, (max(strwidth(as.character(L$label),units="inches",cex=cx))/(8.6-1.0))*2*1.3+0.06)
padR <- max(0.18, (max(strwidth(as.character(R$label),units="inches",cex=cx))/(8.6-1.0))*2*1.3+0.06)

mk <- function(tag, focus=FALSE, whisk=FALSE) {
    M2 <- copy(M)
    M2[, `:=`(lw = if (focus) fifelse(label=="2/3 methods",1.5,0.7) else 0.95,
              al = if (focus) fifelse(label=="2/3 methods",1,0.45) else 1)]
    p <- ggplot(M2, aes(x, accuracy, colour=label, group=label)) +
        # subtle bands so the three architectures read as separate regimes
        annotate("rect", xmin=c(0.5,2.5), xmax=c(1.5,3.5), ymin=-Inf, ymax=Inf,
                 fill=MINOU[["grey"]], alpha=0.06) +
        geom_line(aes(linetype=is_ref, linewidth=lw, alpha=al)) +
        geom_point(aes(alpha=al), size=2.1)
    if (whisk) p <- p + geom_errorbar(aes(ymin=q25, ymax=q75, alpha=al),
                                      width=0.055, linewidth=0.4)
    p + geom_text(data=M2[x==1], aes(x=1-0.04, label=label, alpha=al,
                                     fontface=fifelse(is_ref,"italic","plain")),
                  hjust=1, size=2.9, show.legend=FALSE) +
        geom_text(data=M2[x==3], aes(x=3+0.04, label=label, alpha=al,
                                     fontface=fifelse(is_ref,"italic","plain")),
                  hjust=0, size=2.9, show.legend=FALSE) +
        scale_linetype_manual(values=c(`FALSE`="solid",`TRUE`="22"), guide="none") +
        scale_linewidth_identity() + scale_alpha_identity() +
        scale_colour_manual(values=RULE_COL, guide="none") +
        scale_x_continuous(breaks=1:3, labels=xlabs,
                           expand=expansion(add=c(padL,padR))) +
        labs(tag=tag, x=NULL, y="Accuracy") +
        theme_clinego() +
        theme(plot.tag=element_text(face="bold",size=13), panel.border=element_blank(),
              axis.text.x=element_text(size=8))
}
clinego_save_both(file.path(OUT,"S1_slope_plain"),    mk("C"),                    w=8.6,h=5.0)
clinego_save_both(file.path(OUT,"S2_slope_focus"),    mk("C",focus=TRUE),         w=8.6,h=5.0)
clinego_save_both(file.path(OUT,"S3_slope_whiskers"), mk("C",whisk=TRUE),         w=8.6,h=5.4)
message("  OK S1-S3"); print(dcast(M, label ~ arch, value.var="accuracy"))
