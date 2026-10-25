#!/usr/bin/env Rscript
# Supplementary figure + table for the choice of K = 5 (user 2026-10-16; plan ~/.claude/plans/peaceful-kindling-micali.md,
# Phase 5). Cross-entropy falls monotonically over K 2-8 (ANALYSIS/Trifolium86MainKsearch/README.md) and is not used.
#   a  % variance of PC1-PC15 (LD-pruned set, LEA pca), filled = Tracy-Widom p < 0.05; dashed line between PC5 and PC6,
#      the largest eigenvalue drop of the first ten axes. TW assumes independent individuals, which selfing violates,
#      so the figure leads with the break, not with the count of significant axes.
#   b  number of sNMF clusters (dominant cluster per site) that hold more than one site vs only one site, K 2-8.
# Inputs: <base>.pca/<base>.tracywidom of the project (PCA is upstream of GEA, identical in K5ER and K5ERpc5),
#   ANALYSIS/Trifolium86MainKsearch/ksearch/site_by_K.tsv (same filters, K 2-8, 10 repeats).
# Outputs (paper_pack/): fig_s_k_choice.{png,svg}, k_choice_eigen.tsv, k_choice_sites.tsv (site x K dominant cluster).
# Run inside the image from /pipeline: Rscript ANALYSIS/Trifolium86MainK5ERpc5/maladaptation/fig_k_choice.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
B   <- file.path("/pipeline", Sys.getenv("AN", "ANALYSIS/Trifolium86MainK5ERpc5/maladaptation"))   # env AN: minou preview dir (2026-10-16)
TW  <- "/pipeline/Trifolium86MainK5ERpc5_results/_work/maf0.1_miss0.3_smiss0.8_rel0.999/ld0.2_win100_step10/Trifolium86_chr1-6.pca/Trifolium86_chr1-6.tracywidom"
SBK <- "/pipeline/ANALYSIS/Trifolium86MainKsearch/ksearch/site_by_K.tsv"
OUT <- file.path(B, "paper_pack"); dir.create(OUT, showWarnings = FALSE)
COL <- c(sig = "#0072B2", brk = "#35274A", multi = "#0072B2", single = "#D55E00")
# PALETTE=minou (2026-10-16, plan scalable-doodling-allen.md): the whole paper in minou; default = the zipped v2 colours.
if (Sys.getenv("PALETTE") == "minou") COL <- c(sig = "#00798c", brk = "#2e4057", multi = "#00798c", single = "#d1495b")

# LEA writes double tabs, so fread sees empty columns: keep the non-empty ones (N, eigenvalues, twstats, pvalues,
# effectn, percentage)
tw <- fread(TW)
tw <- tw[, which(vapply(tw, function(x) !all(is.na(x)), TRUE)), with = FALSE]
stopifnot(ncol(tw) == 6)
setnames(tw, c("pc", "eigenvalue", "tw_stat", "tw_p", "effect_n", "fraction"))
tw <- tw[pc <= 15][, `:=`(pct = 100 * fraction, sig = tw_p < 0.05,
                          drop_pct = 100 * (1 - eigenvalue / shift(eigenvalue)))]
brk <- tw[pc <= 10][which.max(drop_pct), pc]   # first axis after the break
stopifnot(brk == 6)
fwrite(tw[, .(pc, eigenvalue, pct = round(pct, 2), drop_from_previous_pct = round(drop_pct, 1), tw_p = signif(tw_p, 3))],
       file.path(OUT, "k_choice_eigen.tsv"), sep = "\t")

pa <- ggplot(tw, aes(pc, pct)) +
    geom_line(colour = "grey60", linewidth = 0.4) +
    geom_point(aes(fill = sig), shape = 21, size = 2.2, colour = "grey20", stroke = 0.4) +
    geom_vline(xintercept = brk - 0.5, linetype = 2, colour = COL[["brk"]], linewidth = 0.5) +
    scale_fill_manual(values = c(`TRUE` = COL[["sig"]], `FALSE` = "white"),
                      labels = c(`TRUE` = "Tracy-Widom p < 0.05", `FALSE` = "p >= 0.05"), name = NULL) +
    scale_x_continuous(breaks = 1:15) +
    labs(title = "a  PCA eigenvalues", x = "Principal component", y = "Variance explained (%)") +
    theme_classic(base_size = 9) +
    theme(legend.position = c(0.75, 0.85), plot.title = element_text(face = "bold"))

sbk <- fread(SBK, colClasses = c(site = "character"))
cl  <- sbk[, .(n_sites = .N), by = .(K, dominant)]
cnt <- cl[, .(`more than one site` = sum(n_sites > 1), `one site` = sum(n_sites == 1)), by = K]
cntl <- melt(cnt, id.vars = "K", variable.name = "kind", value.name = "clusters")
pb <- ggplot(cntl, aes(factor(K), clusters, fill = kind)) +
    geom_col(width = 0.7, colour = "grey20", linewidth = 0.3) +
    scale_fill_manual(values = c(`more than one site` = COL[["multi"]], `one site` = COL[["single"]]), name = "sNMF clusters holding") +
    labs(title = "b  sNMF clusters per K", x = "K", y = "Clusters (dominant per site)") +
    theme_classic(base_size = 9) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))
stopifnot(cnt[K <= 5, all(`one site` == 0)], cnt[K == 6, `one site` > 0])

# KCHOICE_PANELS=a (2026-10-16, plan fancy-scribbling-garden.md Phase 5; user: "keep just (a)"): the eigenvalue panel
# alone, no panel letter, at its own width. Unset = both panels, unchanged.
fig <- pa + pb + plot_layout(widths = c(1.4, 1)); fig_w <- 8
if (Sys.getenv("KCHOICE_PANELS") == "a") { fig <- pa + labs(title = "PCA eigenvalues"); fig_w <- 4.7 }
ggsave(file.path(OUT, "fig_s_k_choice.png"), fig, width = fig_w, height = 3.4, dpi = 250, bg = "white")
ggsave(file.path(OUT, "fig_s_k_choice.svg"), fig, width = fig_w, height = 3.4, bg = "white")
# SUPP_MS=1 (2026-10-17, plan idempotent-conjuring-noodle.md Phase 5 follow-up; user: "redraw S15-S19 at final size"):
# panel a alone at 88 x 60 mm in the manuscript style (benchmarks/ms_style.R ms_save, 7/6 pt Arial), no title.
# Writes fig_s_k_choice_ms.{svg,png} beside the unchanged figure.
if (Sys.getenv("SUPP_MS") == "1") {
    MSE <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MSE)
    MSE$ms_save(file.path(OUT, "fig_s_k_choice_ms"),
                pa + labs(title = NULL) + MSE$theme_ms() + theme(legend.position = c(0.75, 0.85)), 88, 60)
}
fwrite(dcast(sbk, site ~ K, value.var = "dominant"), file.path(OUT, "k_choice_sites.tsv"), sep = "\t")
print(tw[pc <= 8]); print(cnt)
# SUPP2_MS=1 (2026-10-25, plan let-s-follow-the-plan-lively-hellman.md, supplement build B2; every line of this block is
# added): s_k_choice = a, the pipeline's own cross-entropy plot object of the archived K sweep (per K the minimum over the
# 10 repeats, as points), re-themed; b, panel a of fig_s_k_choice_ms with the manuscript's legend labels. Last block on
# purpose: every output above is unchanged. Builds its own ms_style env (the SUPP_MS one exists only inside that block).
# Writes paper_pack/s_k_choice.{svg,png}.
if (Sys.getenv("SUPP2_MS") == "1") {
    MS2 <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MS2)
    stopifnot(identical(unname(COL[["sig"]]), unname(MS2$MINOU[["teal"]])))   # minou env: the significant-point teal
    KS  <- "/pipeline/_archive/trifolium_pre_2026-10-17"
    ce  <- qs::qread(file.path(KS, "Trifolium86MainKsearch_results/_intermediate/cross_entropy_K2-8.qs"))
    # The object's x mapping reads K_START:K_END, globals of scripts/plot_cross_entropy.R that the saved quosure
    # environment does not hold (it is empty): put back the sweep's own values from its config. Nothing is recomputed.
    kcfg <- yaml::read_yaml(file.path(KS, "config_Trifolium86MainKsearch.yaml"))$sNMF
    stopifnot(kcfg$k_start == 2, kcfg$k_end == 8, identical(as.integer(ce$data$name), 2:8))
    qenv <- rlang::quo_get_env(ce$mapping$x)
    assign("K_START", kcfg$k_start, envir = qenv); assign("K_END", kcfg$k_end, envir = qenv)
    ces <- fread("/pipeline/ANALYSIS/Trifolium86MainKsearch/ksearch/cross_entropy.tsv")
    d_ce <- max(abs(ce$data$value - ces$min))
    stopifnot(identical(as.integer(ces$K), 2:8), d_ce < 1e-12,
              all(abs(ce$data$value - c(0.6302, 0.5924, 0.5729, 0.5452, 0.5235, 0.5078, 0.4968)) < 5e-5))
    message("s_k_choice a: cross-entropy (min over repeats) = cross_entropy.tsv min, max |diff| ", signif(d_ce, 3))
    stopifnot(length(ce$layers) == 1, inherits(ce$layers[[1]]$geom, "GeomPoint"))
    l1 <- ce$layers[[1]]; l1$aes_params$colour <- COL[["sig"]]   # was CLINEGO_NEUTRAL #0B775E
    pa_ce <- ce + MS2$theme_ms()
    bce <- ggplot_build(pa_ce)
    stopifnot(identical(bce$layout$panel_params[[1]]$x$get_labels(), paste0("K", 2:8)),
              all(bce$data[[1]]$colour == COL[["sig"]]))
    # b: PC6 is both the first axis after the largest drop among PC1-10 and the first axis with Tracy-Widom p >= 0.05
    stopifnot(brk == 6, tw[tw_p >= 0.05, min(pc)] == 6)
    message("s_k_choice b: significant axes after PC6 (PC1-15): ", paste0("PC", tw[pc > 6 & sig == TRUE, pc], collapse = ", "))
    pb_tw <- pa + scale_fill_manual(values = c(`TRUE` = COL[["sig"]], `FALSE` = "white"),
                                    labels = c(`TRUE` = "Tracy–Widom p < 0.05", `FALSE` = "p ≥ 0.05"), name = NULL) +
        labs(title = NULL) + MS2$theme_ms() + theme(legend.position = c(0.75, 0.85))
    fig_s2 <- pa_ce + pb_tw + plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(face = "bold", size = 8, colour = MS2$MS_INK))
    MS2$ms_save(file.path(OUT, "s_k_choice"), fig_s2, 180, 60)
}
