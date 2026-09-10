#!/usr/bin/env Rscript
# =============================================================================
# mvp_varpart_screen.R -- can anything measurable BEFORE an offset analysis
# predict whether that analysis will work?
#
# Motivation. Panel E showed that realised local adaptation predicts offset
# accuracy (rho +0.53) while mean Fst does not (+0.05). But final_LA requires a
# reciprocal transplant -- the very fitness data offset exists to substitute for
# -- so it is not a usable screening statistic. Variance partitioning IS
# computable from genotypes + coordinates + climate alone, which makes it the
# natural candidate: if a dataset carries more climate signal INDEPENDENT of
# geography and structure, offset ought to work better on it.
#
# This script tests that, per replicate, over all 90 primary MVP replicates
# (mode=climate run for each; see benchmarks/mvp_eval/figures_main/README.md).
#
# Targets tested: offset accuracy (2/3 rule, all-SNPs, oracle) and GEA causal
# recovery. Multiplicity is handled explicitly -- 10 predictors x 4 targets.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_main"))
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))
MINOU <- c(teal="#00798c", red="#d1495b", amber="#edae49", sage="#66a182",
           navy="#2e4057", grey="#8d96a3")

man <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed="character"))
man <- man[arm == "primary"]

# ---- collect varpart per replicate ----------------------------------------
vp <- rbindlist(lapply(man$seed, function(s) {
    f <- file.path(ROOT, sprintf("MVP%s_results/climate/tables/varpart", s))
    v <- file.path(f, "variance_partition.tsv"); if (!file.exists(v)) return(NULL)
    d <- fread(v); out <- as.list(setNames(d$variance_pct, d$component))
    cf <- file.path(f, "climate_confounding.tsv")
    if (file.exists(cf)) { c1 <- fread(cf); out$shared_pct <- c1$shared_pct[1] }
    px <- file.path(f, "px_per_variable.tsv")
    if (file.exists(px)) { p <- fread(px)
        out$Px_climate_partial_geo <- p[variable=="bio_1" & model=="partial_geo", Px_pct][1]
        out$Px_climate_uncond      <- p[variable=="bio_1" & model=="unconditional", Px_pct][1] }
    c(list(seed = s), out)
}), fill = TRUE)
setnames(vp, c("Climate","Structure","Geography"),
         c("climate_unique","structure_unique","geography_unique"), skip_absent = TRUE)
vp <- merge(vp, man[, .(seed, final_LA, meanFst, arch_level)], by = "seed")
fwrite(vp, file.path(OUT, "E2_varpart_per_replicate.tsv"), sep = "\t")
message(sprintf("varpart parsed for %d replicates", nrow(vp)))

# ---- targets ---------------------------------------------------------------
acc <- fread(file.path(EVAL, OFF, "phase1_seed_medians_solo.tsv"), colClasses=c(seed="character"))
acc <- acc[method_label != "RDA-corrected" & seed %in% man$seed]
acc <- acc[, .(v = median(-tau)), by = .(seed, marker_set)]
tg  <- dcast(acc[marker_set %in% c("gea_best","all","adaptive")], seed ~ marker_set, value.var="v")
setnames(tg, c("gea_best","all","adaptive"),
         c("offset: 2/3 methods","offset: all SNPs","offset: oracle"), skip_absent=TRUE)

pr <- fread(file.path(EVAL, OFF, "panel_pr_recomputed.tsv"), colClasses=c(seed="character"))
den <- pr[set=="all", .(seed, C = n_causal)]
rec <- merge(pr[set=="best", .(seed, c = n_causal)], den, by="seed")[C>0,
        .(seed, `GEA causal recall` = c / C)]
TG <- merge(merge(tg, rec, by="seed"), vp, by="seed")

PRED <- c("climate_unique","structure_unique","geography_unique","shared_pct",
          "Px_climate_uncond","Px_climate_partial_geo","final_LA","meanFst")
TARG <- c("offset: 2/3 methods","offset: all SNPs","offset: oracle","GEA causal recall")

res <- rbindlist(lapply(PRED, function(p) rbindlist(lapply(TARG, function(t) {
    x <- TG[[p]]; y <- TG[[t]]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman"))
    data.table(predictor=p, target=t, n=sum(ok),
               rho=unname(ct$estimate), p=ct$p.value)
}))))
# Multiplicity: the varpart screen is 6 predictors x 4 targets = 24 tests.
# final_LA / meanFst are the pre-registered comparisons, corrected separately.
res[, family := fifelse(predictor %in% c("final_LA","meanFst"), "reference", "varpart screen")]
res[, p_holm := p.adjust(p, method="holm"), by = family]
res[, sig := fifelse(p_holm < 0.05, "survives Holm", "not significant")]
res[, absrho := abs(rho)]
setorder(res, family, -absrho)
res[, absrho := NULL]
fwrite(res, file.path(OUT, "E2_varpart_screen.tsv"), sep="\t")
print(res[, .(predictor, target, n, rho=round(rho,3), p=signif(p,2),
              p_holm=signif(p_holm,2), sig)])

# ---- figure: every candidate screen, against every target -------------------
res[, predictor := factor(predictor, levels = rev(PRED))]
pE2 <- ggplot(res, aes(rho, predictor, fill = sig)) +
    geom_vline(xintercept = 0, colour = CLINEGO_COL$fg, linewidth = 0.3) +
    geom_col(width = 0.66) +
    facet_wrap(~ target, nrow = 1) +
    scale_fill_manual(values = c("survives Holm" = MINOU[["sage"]],
                                 "not significant" = MINOU[["grey"]]), name = NULL) +
    scale_x_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3)) +
    labs(tag = "E2", x = "Spearman ρ with the outcome", y = NULL) +
    theme_clinego() +
    theme(plot.tag = element_text(face="bold", size=13), legend.position = "bottom",
          strip.text = element_text(size = 8), axis.text.y = element_text(size = 8),
          panel.border = element_blank())
clinego_save_both(file.path(OUT, "E2_what_predicts_success"), pE2, w = 11, h = 4.4)
message("  OK E2_what_predicts_success")

# ---- scatter: the one that works vs the best varpart candidate --------------
sc <- melt(TG, id.vars = "offset: 2/3 methods",
           measure.vars = c("final_LA","climate_unique","meanFst"),
           variable.name = "screen", value.name = "value")
setnames(sc, "offset: 2/3 methods", "accuracy")
lb <- c(final_LA = "local adaptation\n(needs a transplant experiment)",
        climate_unique = "climate-unique variance\n(computable before analysis)",
        meanFst = "mean F_ST\n(computable before analysis)")
sc[, screen := factor(lb[as.character(screen)], levels = lb)]
rho <- sc[, .(rho = cor(value, accuracy, method="spearman")), by = screen]
pE3 <- ggplot(sc, aes(value, accuracy)) +
    geom_point(alpha = 0.5, size = 1.4, colour = MINOU[["teal"]]) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = MINOU[["navy"]], linewidth = 0.9) +
    geom_text(data = rho, aes(x = Inf, y = -Inf, label = sprintf("ρ = %+.2f", rho)),
              hjust = 1.15, vjust = -1.0, size = 3.3, colour = CLINEGO_COL$fg) +
    facet_wrap(~ screen, scales = "free_x") +
    labs(tag = "E3", x = NULL, y = "Accuracy") +
    theme_clinego() +
    theme(plot.tag = element_text(face="bold", size=13),
          strip.text = element_text(size = 8.5))
clinego_save_both(file.path(OUT, "E3_screen_scatter"), pE3, w = 10, h = 4)
message("  OK E3_screen_scatter")

# =============================================================================
# E4 -- the one varpart relationship that survives inspection, shown honestly.
#
# Pooling all 90 replicates is unsafe here: architectures differ in BOTH the
# varpart components and the detection rate, so a pooled correlation can be a
# between-group artifact. `Climate` is exactly that -- pooled rho +0.32, but
# +0.004 / +0.030 / +0.341 within the three strata, i.e. essentially absent in
# two of them. `shared_pct` (the climate/structure/geography confounded
# fraction) is the opposite: same sign in all three strata and strongest where
# detection is hardest. That is the one worth plotting, faceted, never pooled.
# =============================================================================
det <- merge(pr[set=="best", .(seed, c = n_causal)], den, by="seed")[C>0,
        .(seed, recall = c / C)]
E4 <- merge(det, vp[, .(seed, shared_pct, climate_unique)], by="seed")
E4 <- merge(E4, man[, .(seed, arch_level)], by="seed")
ARCHL <- c(oliogenic="oligogenic", `mod-polygenic`="moderately polygenic",
           `highly-polygenic`="highly polygenic")
E4[, arch := factor(ARCHL[arch_level], levels = ARCHL)]

e4 <- melt(E4, id.vars = c("seed","arch","recall"),
           measure.vars = c("shared_pct","climate_unique"),
           variable.name = "component", value.name = "value")
e4[, component := factor(component, levels = c("shared_pct","climate_unique"),
     labels = c("confounded fraction\n(climate ∩ structure ∩ geography)",
                "climate-unique variance"))]
r4 <- e4[, {ct <- suppressWarnings(cor.test(value, recall, method="spearman"))
            .(rho = unname(ct$estimate), p = ct$p.value, n = .N)},
         by = .(component, arch)]
fwrite(r4, file.path(OUT, "E4_within_stratum_correlations.tsv"), sep="\t")
print(r4)

pE4 <- ggplot(e4, aes(value, recall)) +
    geom_point(alpha = 0.6, size = 1.5, colour = MINOU[["teal"]]) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = MINOU[["navy"]], fill = MINOU[["grey"]],
                alpha = 0.18, linewidth = 0.85) +
    geom_text(data = r4, aes(x = Inf, y = Inf, label = sprintf("ρ = %+.2f\np = %.3f", rho, p)),
              hjust = 1.1, vjust = 1.2, size = 2.9, colour = CLINEGO_COL$fg) +
    facet_grid(component ~ arch, scales = "free") +
    scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) +
    labs(tag = "E4", x = "% of genetic variance", y = "causal loci recovered by the 2/3 rule") +
    theme_clinego() +
    theme(plot.tag = element_text(face="bold", size=13),
          strip.text = element_text(size = 7.5), panel.border = element_blank())
clinego_save_both(file.path(OUT, "E4_varpart_vs_detection"), pE4, w = 10.5, h = 6)
message("  OK E4_varpart_vs_detection")

# =============================================================================
# F -- climate-associated genetic variance vs the QUALITY of the GEA scan.
#
# x is computable from genotypes + coordinates + climate alone (no fitness), so
# unlike final_LA it is a usable pre-flight statistic. It is the sum of every
# varpart fraction that involves climate:
#     Climate(unique) + Climate∩Geography + All three
# NOTE this is what `shared_pct` in climate_confounding.tsv actually tracks
# (rho +0.993 with the sum above); its name suggests "confounding", which
# misreads the sign. The named sum is used here instead.
#
# y is the user's signal-vs-noise question: of the markers a rule returned, what
# fraction is causal or linked to causal, as opposed to neutral noise.
#
# Never pooled across architecture -- architectures differ in both axes, and
# pooling manufactures a between-group correlation (demonstrated for `Climate`
# vs recall, +0.32 pooled but +0.004/+0.030 within two of three strata).
# =============================================================================
CLIM <- vp[, .(seed,
    climate_var = rowSums(cbind(climate_unique,
                                get("Climate ∩ Geography"),
                                get("All three")), na.rm = TRUE))]
shadesF <- function(base, n) c(base,
    grDevices::colorRampPalette(c(base,"white"))(6)[2:(1+(n-1)%/%2+(n-1)%%2)],
    grDevices::colorRampPalette(c(base,"black"))(6)[2:(1+(n-1)%/%2)])[seq_len(n)]
green3 <- shadesF(MINOU[["sage"]], 3); amber3 <- shadesF(MINOU[["amber"]], 3)
LABF <- c(best="2/3 methods", intersect3="3/3 methods", union="1/3 methods",
          solo_rda="RDA", solo_emmax="EMMAX", solo_lfmm="LFMM")
RULE_COL_F <- c("2/3 methods"=green3[1], "3/3 methods"=green3[2], "1/3 methods"=green3[3],
                "RDA"=amber3[1], "EMMAX"=amber3[2], "LFMM"=amber3[3])
prq <- merge(pr[set %in% names(LABF)], CLIM, by = "seed")
prq <- merge(prq, man[, .(seed, arch_level)], by = "seed")
prq[, `:=`(rule = factor(LABF[set], levels = names(RULE_COL_F)),
           arch = factor(ARCHL[arch_level], levels = ARCHL),
           signal_frac = (n_causal + n_linked) / n)]
prq <- prq[is.finite(signal_frac) & n > 0]

rF <- prq[, {ct <- suppressWarnings(cor.test(climate_var, signal_frac, method = "spearman"))
             .(rho = unname(ct$estimate), p = ct$p.value, n = .N)}, by = .(rule, arch)]
setorder(rF, rule, arch)
fwrite(rF, file.path(OUT, "F_climate_vs_gea_quality_correlations.tsv"), sep = "\t")
fwrite(prq[, .(seed, arch, rule, climate_var, n, n_causal, n_linked, signal_frac)],
       file.path(OUT, "F_climate_vs_gea_quality.tsv"), sep = "\t")
print(dcast(rF, rule ~ arch, value.var = "rho"))

pF <- ggplot(prq, aes(climate_var, signal_frac, colour = rule)) +
    geom_point(alpha = 0.42, size = 1.3) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.85) +
    facet_wrap(~ arch, nrow = 1) +
    scale_colour_manual(values = RULE_COL_F, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(tag = "F",
         x = "climate-associated genetic variance (% of total)",
         y = "of markers returned, % causal or linked") +
    theme_clinego() +
    theme(plot.tag = element_text(face = "bold", size = 13),
          legend.position = "bottom", strip.text = element_text(size = 8.5),
          panel.border = element_blank()) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 2)))
clinego_save_both(file.path(OUT, "F_climate_vs_gea_quality"), pF, w = 11, h = 4.6)
message("  OK F_climate_vs_gea_quality")
