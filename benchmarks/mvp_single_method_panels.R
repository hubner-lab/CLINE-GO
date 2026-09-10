#!/usr/bin/env Rscript
# =============================================================================
# mvp_single_method_panels.R -- panels B and C rebuilt on ONE offset method
# instead of the median across methods.
#
# WHY. mvp_main_figure.R collapses the three working offset methods by taking a
# per-replicate median (step 3 of the tau chain). That hides method performance
# and invites the question of whether the marker-set result is carried by one
# method. It is not -- `2/3 methods` ranks first under all three -- but the
# single-method version is the honest way to present it once a method is chosen.
#
# Writes to *_<METHOD>.png/svg/tsv, so the median-across-methods originals are
# never overwritten. Choose with METHOD=<label>.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_main"))
METHOD <- Sys.getenv("METHOD", "LFMM2offset")
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

MINOU <- c(teal="#00798c", red="#d1495b", amber="#edae49", sage="#66a182",
           navy="#2e4057", grey="#8d96a3")
sh <- function(b,n) c(b, grDevices::colorRampPalette(c(b,"white"))(6)[2:(1+(n-1)%/%2+(n-1)%%2)],
                        grDevices::colorRampPalette(c(b,"black"))(6)[2:(1+(n-1)%/%2)])[seq_len(n)]
g3 <- sh(MINOU[["sage"]],3); a3 <- sh(MINOU[["amber"]],3); gy <- sh(MINOU[["grey"]],3)
QTN_BAR  <- grDevices::colorRampPalette(c(MINOU[["navy"]],"white"))(10)[5]
REF_GOOD <- MINOU[["sage"]]; REF_BAD <- MINOU[["red"]]
RULE_COL <- c("2/3 methods"=g3[1], "3/3 methods"=g3[2], "1/3 methods"=g3[3],
              "RDA"=a3[1], "EMMAX"=a3[2], "LFMM"=a3[3],
              "causal loci"=QTN_BAR, "all loci"=gy[3])
OBT   <- c("2/3 methods","3/3 methods","1/3 methods","RDA","EMMAX","LFMM")
ROLE  <- c("2/3 methods"="combined","3/3 methods"="combined","1/3 methods"="combined",
           "RDA"="single","EMMAX"="single","LFMM"="single",
           "causal loci"="reference","all loci"="reference")
SETLAB <- c(gea_best="2/3 methods", gea_strict="3/3 methods", gea_union="1/3 methods",
            gea_rda_only="RDA", gea_emmax_only="EMMAX", gea_lfmm_only="LFMM",
            adaptive="causal loci", all="all loci")
ARCHL <- c(oliogenic="oligogenic", `mod-polygenic`="moderately polygenic",
           `highly-polygenic`="highly polygenic")

man <- fread(file.path(ROOT,"benchmarks/mvp_seeds.tsv"), colClasses=c(seed="character"))[arm=="primary"]
S <- fread(file.path(EVAL,OFF,"phase1_seed_medians_solo.tsv"), colClasses=c(seed="character"))
S <- S[seed %in% man$seed & marker_set %in% names(SETLAB) & method_label == METHOD]
stopifnot(nrow(S) > 0)
S[, `:=`(accuracy = -tau, label = SETLAB[marker_set])]
S <- merge(S, man[, .(seed, arch = factor(ARCHL[arch_level], levels = ARCHL))], by="seed")
message(sprintf("%s: %d replicates x %d marker sets", METHOD, uniqueN(S$seed), uniqueN(S$label)))

# ---------------- B: accuracy per marker set --------------------------------
B <- S[, .(accuracy = median(accuracy), q25 = quantile(accuracy,.25),
           q75 = quantile(accuracy,.75), n_seeds = .N), by = label]
B[, grp := factor(c(combined="combining rules", single="single method",
                    reference="reference")[ROLE[label]],
                  levels = c("reference","single method","combining rules"))]
setorder(B, grp, accuracy)
B[, label := factor(as.character(label), levels = as.character(label))]
fwrite(B, file.path(OUT, sprintf("B_accuracy_by_panel_%s.tsv", METHOD)), sep="\t")
print(B[order(-accuracy), .(label, accuracy = round(accuracy,4), n_seeds)])

rq <- B[label=="causal loci", accuracy]; ra <- B[label=="all loci", accuracy]
pB <- ggplot(B, aes(label, accuracy, fill = label)) +
    geom_hline(yintercept = rq, linetype="22", colour=REF_GOOD, linewidth=0.6) +
    geom_hline(yintercept = ra, linetype="22", colour=REF_BAD,  linewidth=0.6) +
    geom_col(width=0.68) +
    geom_errorbar(aes(ymin=q25, ymax=q75), width=0.2, colour=CLINEGO_COL$fg, linewidth=0.4) +
    geom_text(aes(y=q75, label=sprintf("%.3f", accuracy)), hjust=-0.22, size=2.9,
              colour=CLINEGO_COL$fg) +
    coord_flip(ylim = c(min(B$q25)-0.03, max(B$q75)+0.055)) +
    facet_grid(grp ~ ., scales="free_y", space="free_y", switch="y") +
    scale_fill_manual(values=RULE_COL, guide="none") +
    labs(tag="B", x=NULL, y="Accuracy") +
    theme_clinego() +
    theme(plot.tag=element_text(face="bold", size=13), strip.background=element_blank(),
          strip.placement="outside",
          strip.text.y.left=element_text(angle=90, size=8.5, face="bold"),
          panel.border=element_blank(), panel.spacing.y=unit(0.5,"lines"),
          axis.text.y=element_text(size=8.5))
clinego_save_both(file.path(OUT, sprintf("B_accuracy_by_panel_%s", METHOD)), pB, w=7.8, h=5.2)

# ---------------- C: across architecture ------------------------------------
C <- S[, .(accuracy = median(accuracy)), by = .(label, arch)]
fwrite(C, file.path(OUT, sprintf("C_robustness_%s.tsv", METHOD)), sep="\t")
c_rule <- C[label %in% OBT]; c_ref <- C[!label %in% OBT]
lab_all <- rbind(c_rule[arch==levels(arch)[1L]][, face:="plain"],
                 c_ref [arch==levels(arch)[1L]][, face:="italic"])
.w <- max(strwidth(lab_all$label, units="inches", cex = 2.9 * .pt / par("ps")))
LAB_PAD <- max(0.14, (.w / (7.8-1.0)) * (nlevels(C$arch)-1) * 1.25 + 0.05)

pC <- ggplot(mapping = aes(as.integer(arch), accuracy, group = label)) +
    geom_line(data=c_ref,  aes(colour=label), linetype="22", linewidth=0.75) +
    geom_line(data=c_rule, aes(colour=label), linewidth=0.95) +
    geom_point(data=c_rule, aes(colour=label), size=2.2) +
    geom_text(data=lab_all, aes(x=1-0.03, y=accuracy, label=label, colour=label,
                                fontface=face), hjust=1, size=2.9, show.legend=FALSE) +
    scale_x_continuous(breaks=seq_len(nlevels(C$arch)), labels=levels(C$arch),
                       expand=expansion(add=c(LAB_PAD, 0.06))) +
    scale_colour_manual(values=c(RULE_COL[OBT], "causal loci"=REF_GOOD,
                                 "all loci"=REF_BAD), guide="none") +
    labs(tag="C", x=NULL, y="Accuracy") +
    theme_clinego() +
    theme(plot.tag=element_text(face="bold", size=13), panel.border=element_blank(),
          axis.text.x=element_text(size=8.5))
clinego_save_both(file.path(OUT, sprintf("C_robustness_%s", METHOD)), pC, w=7.8, h=5.0)
message(sprintf("  OK B/C for %s", METHOD))
