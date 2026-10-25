#!/usr/bin/env Rscript
# Draft candidates for Fig 3 and S9-S13 (Trifolium, EMMAX + RDA at sNMF K 5 with 5 PCs), all in minou
# (plan ~/.claude/plans/scalable-doodling-allen.md, Phase 2). Deliberately redundant: the user prunes on the review page.
# Every panel is drawn at final size in mm by ms_save() (editable SVG + 600 dpi PNG); no titles, captions explain.
# Inputs are pipeline tables of Trifolium86MainK5ERpc5 (GEA), Trifolium86MainK5PC5 / Trifolium86MainImpMAF (LFMM
# supplement, K 5 vs K 3) and the Phase 7 region tables. Genome-wide point clouds are rasterized with scattermore
# (CLAUDE.md rule 6); only the 35 significant region SNPs are vector points.
# Writes manuscript_fig3/draft/<id>.{svg,png} and draft/figures.tsv (id, file, width_mm, height_mm, source tables).
# Run inside the image from /pipeline: Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/build_fig3.R [only=<id prefix>]
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(patchwork); library(scattermore); library(ggrepel)
})
source("/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/fig3_style.R")
ONLY <- sub("^only=", "", grep("^only=", commandArgs(TRUE), value = TRUE))
AN  <- "/pipeline/ANALYSIS/Trifolium86MainK5ERpc5"
OUT <- file.path(AN, "manuscript_fig3/draft"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P5  <- "/pipeline/Trifolium86MainK5ERpc5_results"
T5  <- file.path(P5, "GEA/tables")
RT  <- file.path(AN, "regions/tables")
DUP <- c("bio_11", "bio_19")   # identical to bio_8 / bio_16 at every plant
reg_items <- list()
want <- function(id) !length(ONLY) || startsWith(id, ONLY)
save_item <- function(id, p, w, h, src) {
    ms_save(file.path(OUT, id), p, w, h)
    reg_items[[id]] <<- data.table(id = id, file = paste0(id, ".svg"), width_mm = w, height_mm = h, sources = src)
    message("wrote ", id)
}
nl10 <- function(p) -log10(pmax(p, 1e-300))

# ---- shared data ------------------------------------------------------------------------------------------------------
em  <- fread(file.path(T5, "methods/EMMAX/EMMAX_pvalues_K5.tsv"), colClasses = list(character = "SNPID"))
rda <- fread(file.path(T5, "methods/RDA/RDA_pvalues_K5.tsv"), colClasses = list(character = "SNPID"))
stopifnot(identical(em$SNPID, rda$SNPID))
PRED <- setdiff(grep("^bio_", names(em), value = TRUE), DUP)
stopifnot(length(PRED) == 14)
em[, minp := do.call(pmin, .SD), .SDcols = PRED]
em[, rda := rda$climate_multivariate]
chr_len <- em[, .(len = max(pos)), by = chr][order(chr)]
chr_len[, off := cumsum(c(0, head(len, -1)))]
em[chr_len, x := pos + i.off, on = "chr"]
em[, chr_tint := FIG3_CHR[(chr %% 2) + 1]]
chr_mid <- chr_len[, .(chr, mid = off + len / 2)]
thr_e <- fread(file.path(T5, "methods/EMMAX/EMMAX_pvalues_K5_sig_snps_bonf_0.05.tsv"))$pval_threshold[1]
thr_r <- fread(file.path(T5, "methods/RDA/RDA_pvalues_K5_sig_snps_qval_0.05.tsv"))$pval_threshold[1]
sel <- fread(file.path(T5, "selected_snps.tsv"), colClasses = list(character = c("SNPID", "EMMAX", "RDA")))
sel[, status := fifelse(EMMAX != "" & RDA != "", "both methods", fifelse(EMMAX != "", "EMMAX only", "RDA only"))]
stopifnot(nrow(sel) == 30)
sel <- merge(sel[, .(SNPID, status)], em[, .(SNPID, chr, pos, x, minp, rda)], by = "SNPID")
regs <- fread(file.path(T5, "regions_combined.tsv"))
regs[chr_len, `:=`(x0 = start + i.off, x1 = end + i.off), on = "chr"]
x_chr <- scale_x_continuous(breaks = chr_mid$mid, labels = chr_mid$chr, expand = expansion(mult = 0.005))
status_scale <- scale_colour_manual(values = FIG3_STATUS, name = "Agreement SNPs", drop = FALSE)

# ---- structure ----------------------------------------------------------------------------------------------------------
LDW <- file.path(P5, "_work/maf0.1_miss0.3_smiss0.8_rel0.999/ld0.2_win100_step10")
cl  <- fread(file.path(P5, "PreStructure/tables/K5/clusters_K5.tsv"), colClasses = list(character = c("sample", "site")))
fam <- fread(file.path(LDW, "Trifolium86_chr1-6.fam"), header = FALSE, colClasses = "character")$V2
stopifnot(identical(fam, cl$sample))
QC  <- paste0("C", 1:5)
cl[, cluster := factor(QC[max.col(as.matrix(.SD), ties.method = "first")], QC), .SDcols = QC]
CL_COL <- setNames(FIG3_CLUSTERS, QC)
eig <- fread(file.path(AN, "maladaptation/paper_pack/k_choice_eigen.tsv"))
pr  <- fread(file.path(LDW, "Trifolium86_chr1-6.pca/Trifolium86_chr1-6.projections"), header = FALSE)
pca <- cbind(cl[, .(sample, site, cluster)], pr[, 1:4])
setnames(pca, paste0("V", 1:4), paste0("PC", 1:4))
pc_lab <- function(i) sprintf("PC%d (%.1f%%)", i, eig$pct[eig$pc == i])
site_ctr <- function(a, b) pca[, .(x = mean(get(a)), y = mean(get(b))), by = site]
pca_plot <- function(a, b) {
    ggplot(pca, aes(.data[[a]], .data[[b]])) +
        geom_point(aes(fill = cluster), shape = 21, size = 1.6, stroke = 0.2, colour = MS_INK) +
        geom_text_repel(data = site_ctr(a, b), aes(x, y, label = site), size = MS_LAB, colour = MS_INK,
                        segment.size = 0.2, min.segment.length = 0, box.padding = 0.3, seed = 1) +
        scale_fill_manual(values = CL_COL, name = "sNMF cluster (K = 5)") +
        labs(x = pc_lab(as.integer(sub("PC", "", a))), y = pc_lab(as.integer(sub("PC", "", b)))) + theme_ms()
}
if (want("st_pca")) {
    save_item("st_pca12", pca_plot("PC1", "PC2"), 88, 80, "LEA .projections; k_choice_eigen.tsv; clusters_K5.tsv")
    save_item("st_pca34", pca_plot("PC3", "PC4"), 88, 80, "LEA .projections; k_choice_eigen.tsv; clusters_K5.tsv")
}
if (want("st_admix")) {
    q <- melt(cl, id.vars = c("sample", "site", "cluster"), measure.vars = QC, variable.name = "k", value.name = "q")
    site_ord <- cl[, .N, by = .(site, cluster)][order(cluster, -N)][!duplicated(site), site]
    q[, site := factor(site, site_ord)]
    cl[, qmax := do.call(pmax, .SD), .SDcols = QC]
    smp_ord <- cl[order(factor(site, site_ord), cluster, -qmax), sample]
    q[, sample := factor(sample, smp_ord)]
    p <- ggplot(q, aes(sample, q, fill = k)) + geom_col(width = 1) +
        facet_grid(~site, scales = "free_x", space = "free_x", switch = "x") +
        scale_fill_manual(values = CL_COL, name = "sNMF cluster (K = 5)") +
        scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) + labs(x = NULL, y = "Ancestry") + theme_ms() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
              panel.spacing = unit(0.6, "mm"), strip.text = element_text(size = MS_TICK, angle = 90, hjust = 1),
              strip.placement = "outside")
    save_item("st_admix", p, 180, 45, "clusters_K5.tsv")
}

# ---- association: Miami, stacked tracks, per-predictor, QQ --------------------------------------------------------------
bg_layer <- function(d, y) geom_scattermore(data = d, aes(x, .data[[y]], colour = I(chr_tint)), pointsize = 1.6,
                                             pixels = c(2400, 900))
reg_rect <- function(ymin, ymax) geom_rect(data = regs, aes(xmin = x0, xmax = x1), ymin = ymin, ymax = ymax,
                                           inherit.aes = FALSE, fill = FIG3_REGION, colour = NA)
if (want("as_miami") || want("as_tracks")) {
    em[, `:=`(ye = nl10(minp), yr = -nl10(rda))]
    sel[, `:=`(ye = nl10(minp), yr = -nl10(rda), status = factor(status, names(FIG3_STATUS)))]
    lim_e <- ceiling(max(em$ye)); lim_r <- ceiling(max(-em$yr))
    if (want("as_miami")) {
        p <- ggplot() + reg_rect(-Inf, Inf) + bg_layer(em, "ye") + bg_layer(em, "yr") +
            geom_hline(yintercept = nl10(thr_e), colour = FIG3_THR[["EMMAX"]], linewidth = 0.3, linetype = "22") +
            geom_hline(yintercept = -nl10(thr_r), colour = FIG3_THR[["RDA"]], linewidth = 0.3, linetype = "22") +
            geom_hline(yintercept = 0, colour = MS_INK, linewidth = 0.3) +
            geom_point(data = sel, aes(x, ye, colour = status), size = 1.1) +
            geom_point(data = sel, aes(x, yr, colour = status), size = 1.1) +
            annotate("text", x = -Inf, y = lim_e, label = "EMMAX", hjust = -0.2, vjust = 1, size = MS_LAB) +
            annotate("text", x = -Inf, y = -lim_r, label = "RDA", hjust = -0.2, vjust = 0, size = MS_LAB) +
            x_chr + scale_y_continuous(labels = function(v) abs(v), breaks = scales::pretty_breaks(8)) +
            status_scale + labs(x = "Chromosome", y = "-log10 P") + theme_ms()
        save_item("as_miami", p, 180, 90, "EMMAX_pvalues_K5.tsv (min over 14 distinct predictors); RDA_pvalues_K5.tsv; selected_snps.tsv; regions_combined.tsv")
    }
    if (want("as_tracks")) {
        track <- function(y, thr, col, lab, sgn) {
            ggplot() + reg_rect(-Inf, Inf) + bg_layer(em[, .(x, chr_tint, y = sgn * get(y))], "y") +
                geom_hline(yintercept = nl10(thr), colour = col, linewidth = 0.3, linetype = "22") +
                geom_point(data = sel[, .(x, status, y = sgn * get(y))], aes(x, y, colour = status), size = 1.1) +
                x_chr + status_scale + labs(x = NULL, y = paste(lab, "-log10 P")) + theme_ms()
        }
        p <- track("ye", thr_e, FIG3_THR[["EMMAX"]], "EMMAX", 1) / track("yr", thr_r, FIG3_THR[["RDA"]], "RDA", -1) +
            plot_layout(guides = "collect") & theme(legend.position = "bottom")
        save_item("as_tracks", p, 180, 100, "EMMAX_pvalues_K5.tsv; RDA_pvalues_K5.tsv; selected_snps.tsv; regions_combined.tsv")
    }
}
if (want("as_pred14")) {
    long <- melt(em[, c("x", "chr", PRED), with = FALSE], id.vars = c("x", "chr"), variable.name = "pred", value.name = "p")
    long[, `:=`(y = nl10(p), pred = factor(pred, PRED))]
    tint <- vapply(FIG3_PRED_COL, function(c) mix("#FFFFFF", c, 0.45), "")   # once per predictor, not per row
    long[, col := fifelse(chr %% 2 == 1, FIG3_PRED_COL[as.character(pred)], tint[as.character(pred)])]
    p <- ggplot(long, aes(x, y, colour = I(col))) + geom_scattermore(pointsize = 1.4, pixels = c(1200, 300)) +
        geom_hline(yintercept = nl10(thr_e), colour = FIG3_THR[["EMMAX"]], linewidth = 0.25, linetype = "22") +
        facet_wrap(~pred, ncol = 2, scales = "free_y", labeller = as_labeller(function(v)
            ifelse(v == "bio_8", "bio_8 (= bio_11)", ifelse(v == "bio_16", "bio_16 (= bio_19)", v)))) +
        x_chr + labs(x = "Chromosome", y = "EMMAX -log10 P") + theme_ms() +
        theme(strip.text = element_text(size = MS_TICK))
    save_item("as_pred14", p, 180, 200, "EMMAX_pvalues_K5.tsv (14 distinct predictors)")
}
if (want("as_qq")) {
    qq <- function(p, lab) {
        p <- sort(p); n <- length(p); e <- (seq_len(n) - 0.5) / n
        keep <- unique(c(which(p < 1e-3), round(seq(1, n, length.out = 3000))))
        data.table(series = lab, exp = nl10(e[keep]), obs = nl10(p[keep]))
    }
    lam <- function(p) qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1)
    lab_e <- sapply(PRED, function(v) sprintf("%s (λ %.2f)", v, lam(em[[v]])))
    d <- rbind(rbindlist(lapply(PRED, function(v) qq(em[[v]], lab_e[[v]]))), qq(em$rda, sprintf("RDA (λ %.2f after GIF)", lam(em$rda))))
    cols <- c(setNames(FIG3_PRED_COL[PRED], lab_e[PRED]), setNames(FIG3_THR[["RDA"]], tail(unique(d$series), 1)))
    d[, series := factor(series, names(cols))]
    p <- ggplot(d, aes(exp, obs, colour = series)) + geom_abline(colour = MS_REF, linewidth = 0.3) +
        geom_point(size = 0.35) + scale_colour_manual(values = cols, name = NULL) +
        labs(x = "Expected -log10 P", y = "Observed -log10 P") +
        theme_ms() + theme(legend.position = "right") + guides(colour = guide_legend(override.aes = list(size = 1.5)))
    save_item("as_qq", p, 120, 90, "EMMAX_pvalues_K5.tsv; RDA_pvalues_K5.tsv")
}
# [added 2026-10-25, supplement build B4] s_qq = as_qq with the manuscript's labels (deviation D28): "BIO1 (λ 1.03)"
# in numeric order BIO1 ... BIO16 (the 14 distinct EMMAX scans; DUP left out as in as_qq) and "RDA (λ 1.00 after
# correction)". Same points, colours (FIG3_PRED_COL, RDA sage) and size as as_qq; as_qq itself is unchanged.
if (want("s_qq")) {
    qq <- function(p, lab) {
        p <- sort(p); n <- length(p); e <- (seq_len(n) - 0.5) / n
        keep <- unique(c(which(p < 1e-3), round(seq(1, n, length.out = 3000))))
        data.table(series = lab, exp = nl10(e[keep]), obs = nl10(p[keep]))
    }
    lam <- function(p) qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1)
    pred_n <- PRED[order(as.integer(sub("^bio_", "", PRED)))]
    stopifnot(length(pred_n) == 14, identical(sub("^bio_", "", pred_n), as.character(c(1:10, 12, 13, 15, 16))))
    lam_e <- vapply(pred_n, function(v) lam(em[[v]]), 0)
    lab_e <- sprintf("BIO%s (λ %.2f)", sub("^bio_", "", pred_n), lam_e); names(lab_e) <- pred_n
    lam_r <- lam(em$rda)
    lab_r <- sprintf("RDA (λ %.2f after correction)", lam_r)
    stopifnot(lab_r == "RDA (λ 1.00 after correction)")
    d <- rbind(rbindlist(lapply(pred_n, function(v) qq(em[[v]], lab_e[[v]]))), qq(em$rda, lab_r))
    cols <- c(setNames(FIG3_PRED_COL[pred_n], lab_e[pred_n]), setNames(FIG3_THR[["RDA"]], lab_r))
    d[, series := factor(series, names(cols))]
    fwrite(data.table(series = c(pred_n, "RDA"), label = names(cols), lambda = c(lam_e, lam_r), colour = unname(cols),
                      class = c(ifelse(pred_n %in% FIG3_TEMP, "temperature", "precipitation"), "RDA")),
           file.path(OUT, "s_qq_lambda.tsv"), sep = "\t")
    p <- ggplot(d, aes(exp, obs, colour = series)) + geom_abline(colour = MS_REF, linewidth = 0.3) +
        geom_point(size = 0.35) + scale_colour_manual(values = cols, name = NULL) +
        labs(x = "Expected -log10 P", y = "Observed -log10 P") +
        theme_ms() + theme(legend.position = "right") + guides(colour = guide_legend(override.aes = list(size = 1.5)))
    save_item("s_qq", p, 120, 90, "EMMAX_pvalues_K5.tsv; RDA_pvalues_K5.tsv")
}

# ---- association v2 (plan ~/.claude/plans/fancy-scribbling-garden.md, Phase 2; user review 2026-10-16) -----------------
# Miami: bigger background; significant SNPs in three classes, each drawn on the side(s) where it is significant:
#   same SNP significant in both methods; linked = significant in one method with a significant SNP of the other within
#   the agreement distance (the rest of selected_snps.tsv); one method, no partner = significant, outside the set.
sig_e <- fread(file.path(T5, "methods/EMMAX/EMMAX_pvalues_K5_sig_snps_bonf_0.05.tsv"), colClasses = list(character = "SNPID"))
sig_r <- fread(file.path(T5, "methods/RDA/RDA_pvalues_K5_sig_snps_qval_0.05.tsv"), colClasses = list(character = "SNPID"))
AGREE_D <- unique(c(sig_e$overlap_distance, sig_r$overlap_distance))
stopifnot(length(AGREE_D) == 1, AGREE_D == 130811)
if (want("as_miami_v2")) {
    se <- unique(sig_e[trait %in% PRED, SNPID]); sr <- unique(sig_r$SNPID)
    stopifnot(setequal(se, unique(sig_e$SNPID)))   # the dropped duplicates add no SNP
    cls <- rbind(data.table(SNPID = se, side = "EMMAX"), data.table(SNPID = sr, side = "RDA"))
    cls[, class := fifelse(SNPID %in% intersect(se, sr), names(FIG3_AGREE)[1],
                   fifelse(SNPID %in% sel$SNPID, names(FIG3_AGREE)[2], names(FIG3_AGREE)[3]))]
    # the split must reproduce selected_snps.tsv: 15 identical + 15 linked = 30, identical = sel's "both methods"
    stopifnot(setequal(cls[class == names(FIG3_AGREE)[1], SNPID], sel[status == "both methods", SNPID]),
              setequal(cls[class != names(FIG3_AGREE)[3], SNPID], sel$SNPID))
    # v3 (user review v2): two classes only. Agreement = identical + linked = selected_snps.tsv; the rest = one method.
    cls[, class := fifelse(class == names(FIG3_AGREE)[3], names(FIG3_AGREE2)[2], names(FIG3_AGREE2)[1])]
    n_cls <- cls[, .(n_snps = uniqueN(SNPID)), by = .(side, class)][order(side, class)]
    print(n_cls)
    stopifnot(n_cls[side == "EMMAX" & class == names(FIG3_AGREE2)[1], n_snps] == 23,
              n_cls[side == "EMMAX" & class == names(FIG3_AGREE2)[2], n_snps] == 314,
              n_cls[side == "RDA" & class == names(FIG3_AGREE2)[1], n_snps] == 22,
              n_cls[side == "RDA" & class == names(FIG3_AGREE2)[2], n_snps] == 17)
    fwrite(n_cls, file.path(OUT, "as_miami_v2_classes.tsv"), sep = "\t")
    cls <- merge(cls, em[, .(SNPID, x, ye = nl10(minp), yr = -nl10(rda))], by = "SNPID")
    cls[, `:=`(y = fifelse(side == "EMMAX", ye, yr), class = factor(class, names(FIG3_AGREE2)))]
    setorder(cls, -class)   # agreement drawn last, on top
    lim_e <- ceiling(max(nl10(em$minp))); lim_r <- ceiling(max(nl10(em$rda)))
    # thresholds in ink: red now marks significant SNPs; the side of the axis names the method
    # 2026-10-18 (D15): no region shading behind the Miami (reg_rect stays in the other plots)
    p <- ggplot() +
        geom_scattermore(data = em, aes(x, nl10(minp), colour = I(chr_tint)), pointsize = 3.2, pixels = c(2400, 1200)) +
        geom_scattermore(data = em, aes(x, -nl10(rda), colour = I(chr_tint)), pointsize = 3.2, pixels = c(2400, 1200)) +
        geom_hline(yintercept = nl10(thr_e), colour = MS_INK, linewidth = 0.3, linetype = "22") +
        geom_hline(yintercept = -nl10(thr_r), colour = MS_INK, linewidth = 0.3, linetype = "22") +
        geom_hline(yintercept = 0, colour = MS_INK, linewidth = 0.3) +
        geom_point(data = cls, aes(x, y, fill = class, size = class), shape = 21, colour = MS_INK, stroke = 0.12) +
        scale_fill_manual(values = FIG3_AGREE2, name = "Significant SNPs", drop = FALSE) +
        scale_size_manual(values = setNames(c(1.6, 0.8), names(FIG3_AGREE2)), name = "Significant SNPs", drop = FALSE) +
        annotate("text", x = -Inf, y = lim_e, label = "EMMAX", hjust = -0.2, vjust = 1, size = MS_LAB) +
        annotate("text", x = -Inf, y = -lim_r, label = "RDA", hjust = -0.2, vjust = 0, size = MS_LAB) +
        x_chr + scale_y_continuous(labels = function(v) abs(v), breaks = scales::pretty_breaks(8)) +
        labs(x = "Chromosome", y = "-log10 P") + theme_ms()
    save_item("as_miami_v2", p, 180, 90, paste("EMMAX_pvalues_K5.tsv (min over 14 distinct predictors); RDA_pvalues_K5.tsv;",
              "EMMAX/RDA sig_snps tables; selected_snps.tsv; regions_combined.tsv; class counts as_miami_v2_classes.tsv"))
}
# Per-predictor EMMAX: bigger background dots; each predictor's significant SNPs as larger red points on top.
if (want("as_pred14_v2")) {
    # all 16 predictors (user 2026-10-17: draw bio_11 and bio_19 too, no "bio_8 (= bio_11)" labels)
    PRED16 <- grep("^bio_", names(em), value = TRUE)
    long <- melt(em[, c("SNPID", "x", "chr", PRED16), with = FALSE], id.vars = c("SNPID", "x", "chr"),
                 variable.name = "pred", value.name = "p")
    # numeric bioclim order (user review v2). Both layers must carry the same factor: a character `pred` in the
    # significant layer made facet_wrap sort the panels as text (bio_1, bio_10, bio_12, ...) in v2.
    PRED_NUM <- PRED16[order(as.integer(sub("bio_", "", PRED16)))]
    long[, `:=`(y = nl10(p), pred = factor(pred, PRED_NUM))]
    tint <- vapply(FIG3_PRED_COL, function(c) mix("#FFFFFF", c, 0.45), "")
    long[, col := fifelse(chr %% 2 == 1, FIG3_PRED_COL[as.character(pred)], tint[as.character(pred)])]
    sig14 <- merge(long, unique(sig_e[trait %in% PRED16, .(SNPID, pred = factor(trait, PRED_NUM))]), by = c("SNPID", "pred"))
    stopifnot(is.factor(sig14$pred), identical(levels(sig14$pred), PRED_NUM))
    stopifnot(all(sig14$p <= thr_e))
    message("as_pred14_v2: ", nrow(sig14), " significant SNP x predictor points")
    p <- ggplot(long, aes(x, y, colour = I(col))) + geom_scattermore(pointsize = 2.8, pixels = c(1200, 300)) +
        geom_hline(yintercept = nl10(thr_e), colour = MS_INK, linewidth = 0.25, linetype = "22") +
        geom_point(data = sig14, aes(x, y), inherit.aes = FALSE, shape = 21, size = 0.9, stroke = 0.1,
                   fill = MINOU[["red"]], colour = MS_INK) +
        facet_wrap(~pred, ncol = 2, scales = "free_y") +
        x_chr + labs(x = "Chromosome", y = "EMMAX -log10 P") + theme_ms() +
        theme(strip.text = element_text(size = MS_TICK))
    lay <- ggplot_build(p)$layout$layout
    stopifnot(identical(as.character(lay$pred[order(lay$PANEL)]), PRED_NUM))
    save_item("as_pred14_v2", p, 180, 225, "EMMAX_pvalues_K5.tsv (all 16 predictors; bio_11 = bio_8, bio_19 = bio_16 at every plant); EMMAX sig_snps table")
}

# ---- regions --------------------------------------------------------------------------------------------------------------
calls <- fread(file.path(AN, "regions/figures/input/regionplot_calls.tsv"), colClasses = "character")
calls[, order := as.integer(order)]
genes <- fread(file.path(AN, "regions/figures/input/topr_labelled.tsv"))
rsig  <- fread(file.path(RT, "region_sig_snps.tsv"), colClasses = list(character = "snp_id"))
region_plot <- function(i, legend = TRUE, base_h = 46) {
    cr <- calls[i]
    w  <- as.numeric(strsplit(sub("^\\d+:", "", cr$custom_region), "-")[[1]]); ch <- as.integer(sub(":.*", "", cr$custom_region))
    edges <- as.numeric(strsplit(cr$vlines, ",")[[1]])
    series <- rbindlist(lapply(strsplit(cr$traits_arg, ",")[[1]], function(s) {
        tm <- strsplit(s, ":")[[1]]; data.table(trait = tm[1], method = tm[2])
    }))
    series <- series[!trait %in% DUP]
    win <- em[chr == ch & pos >= w[1] & pos <= w[2]]
    pts <- rbindlist(lapply(seq_len(nrow(series)), function(j) {
        s <- series[j]
        v <- if (s$method == "RDA") win$rda else win[[s$trait]]
        data.table(pos = win$pos, y = nl10(v), series = if (s$method == "RDA") "RDA" else s$trait, SNPID = win$SNPID)
    }))
    lev <- c(intersect(PRED, unique(pts$series)), "RDA")
    pts[, series := factor(series, lev)]
    cols <- c(FIG3_PRED_COL, RDA = FIG3_THR[["RDA"]])[lev]
    hi <- pts[SNPID %in% rsig[region_id == cr$region_id, snp_id]]
    main <- ggplot(pts, aes(pos / 1e6, y, colour = series, shape = series == "RDA")) +
        annotate("rect", xmin = edges[1] / 1e6, xmax = edges[2] / 1e6, ymin = -Inf, ymax = Inf, fill = FIG3_REGION, alpha = 0.5) +
        geom_hline(yintercept = nl10(thr_e), colour = FIG3_THR[["EMMAX"]], linewidth = 0.3, linetype = "22") +
        geom_hline(yintercept = nl10(thr_r), colour = FIG3_THR[["RDA"]], linewidth = 0.3, linetype = "22") +
        geom_point(size = 0.45, alpha = 0.7) + geom_point(data = hi, size = 1.3, stroke = 0.3) +
        scale_colour_manual(values = cols, name = NULL, labels = function(v) sub("^bio_8$", "bio_8 (= bio_11)", v)) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
        scale_x_continuous(limits = w / 1e6, expand = c(0, 0)) +
        labs(x = NULL, y = "-log10 P") + theme_ms() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
              legend.position = if (legend) "right" else "none")
    g <- genes[chrom == ch & gene_end >= w[1] & gene_start <= w[2]][order(gene_start)]
    g[, row := 0L]; ends <- numeric(0)
    for (k in seq_len(nrow(g))) {       # greedy rows: a gene goes to the first row whose last gene ended 2 kb earlier
        r <- which(ends < g$gene_start[k] - 2000)[1]
        if (is.na(r)) { ends <- c(ends, g$gene_end[k]); r <- length(ends) } else ends[r] <- g$gene_end[k]
        g$row[k] <- r
    }
    lab <- g[gene_symbol != ""]
    gp <- ggplot(g) +
        annotate("rect", xmin = edges[1] / 1e6, xmax = edges[2] / 1e6, ymin = -Inf, ymax = Inf, fill = FIG3_REGION, alpha = 0.5) +
        geom_rect(aes(xmin = pmax(gene_start, w[1]) / 1e6, xmax = pmin(gene_end, w[2]) / 1e6, ymin = -row - 0.3, ymax = -row + 0.3),
                  fill = FIG3_GENE, colour = NA) +
        geom_text_repel(data = lab, aes(x = (gene_start + gene_end) / 2e6, y = -row, label = gene_symbol), size = MS_LAB * 0.9,
                        colour = MS_INK, segment.size = 0.15, min.segment.length = 0, direction = "both",
                        max.overlaps = Inf, box.padding = 0.15, seed = 1) +
        scale_x_continuous(limits = w / 1e6, expand = c(0, 0)) +
        labs(x = sprintf("Chromosome %d (Mb)", ch), y = NULL) + theme_ms() +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
    main / gp + plot_layout(heights = c(base_h, 8 + 3.5 * max(g$row, 1)))
}
if (want("rg_")) {
    for (i in seq_len(nrow(calls))) {
        id <- sprintf("rg_%02d_%s", calls$order[i], calls$region_id[i])
        if (!want(id)) next
        save_item(id, region_plot(i), 180, 46 + 12 + 3.5 * 6,
                  "EMMAX_pvalues_K5.tsv; RDA_pvalues_K5.tsv; regions/figures/input/{regionplot_calls,topr_labelled}.tsv; region_sig_snps.tsv")
    }
}
if (want("rg_carriers")) {
    pat <- fread(file.path(RT, "sig_snp_patterns.tsv"), colClasses = list(character = "snp_id"))
    sizes <- cl[, .N, by = site]
    car <- rbindlist(lapply(seq_len(nrow(pat)), function(j) {
        s <- strsplit(pat$carrier_sites[j], ";")[[1]]
        m <- regmatches(s, regexec("^(\\d+)\\((\\d+)/(\\d+)\\)$", s))
        d <- rbindlist(lapply(m, function(x) data.table(site = x[2], share = as.numeric(x[3]) / as.numeric(x[4]))))
        merge(sizes[, .(site)], d, by = "site", all.x = TRUE)[, `:=`(snp_id = pat$snp_id[j], region_id = pat$region_id[j])]
    }))
    car[is.na(share), share := 0]
    car <- merge(car, rsig[, .(snp_id, agreed, sig_methods)], by = "snp_id")
    car[, lab := paste0(snp_id, fifelse(agreed, "", " *"))]
    ord <- unique(car[order(match(region_id, calls$region_id), as.numeric(sub(".*:", "", snp_id))), lab])
    car[, `:=`(lab = factor(lab, rev(ord)), region_id = factor(region_id, calls$region_id))]
    p <- ggplot(car, aes(site, lab, fill = share)) + geom_tile(colour = "white", linewidth = 0.2) +
        facet_grid(region_id ~ ., scales = "free_y", space = "free_y", switch = "y") +
        scale_fill_gradientn(colours = MS_SEQUENTIAL, limits = c(0, 1), labels = scales::percent,
                             name = "Plants carrying the minor allele") +
        labs(x = "Site", y = NULL) + theme_ms() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), strip.text.y.left = element_text(angle = 0, size = MS_TICK * 0.8),
              strip.placement = "outside", panel.spacing = unit(0.4, "mm"), axis.text.y = element_text(size = MS_TICK * 0.8))
    save_item("rg_carriers", p, 150, 190, "regions/tables/sig_snp_patterns.tsv; region_sig_snps.tsv (* = significant, not in the agreement set)")
}
if (want("rg_grid")) {
    pl <- lapply(seq_len(nrow(calls)), function(i) region_plot(i, legend = FALSE, base_h = 30) +
                     plot_annotation(title = calls$region_id[i]) & theme(plot.title = element_text(size = MS_TICK)))
    pl <- lapply(pl, wrap_elements)
    p <- wrap_plots(pl, ncol = 3)
    save_item("rg_grid", p, 180, 240, "as rg_*")
}

# ---- stress screen ------------------------------------------------------------------------------------------------------
if (want("sx_")) {
    st <- fread(file.path(AN, "regions/set_stress.tsv"))
    d <- rbind(st[, .(set, line = "GO stress terms", observed = go_obs, expected = go_exp, ratio = go_ratio, p = go_p)],
               st[, .(set, line = "Stress keywords", observed = kw_strict_obs, expected = kw_strict_exp, ratio = kw_strict_ratio, p = kw_strict_p)])
    d[, set := factor(fifelse(set == "PC5", "5 PCs (final)", "3 PCs"), c("3 PCs", "5 PCs (final)"))]
    long <- melt(d, measure.vars = c("observed", "expected"), variable.name = "what", value.name = "arrays")
    p <- ggplot(long, aes(set, arrays, fill = what)) + geom_col(position = position_dodge(0.75), width = 0.7) +
        geom_text(data = d, aes(set, pmax(observed, expected), label = sprintf("%.2f×\nP = %s", ratio, format(p, digits = 2))),
                  inherit.aes = FALSE, vjust = -0.2, size = MS_LAB, lineheight = 0.9) +
        facet_wrap(~line) + scale_fill_manual(values = c(observed = MINOU[["navy"]], expected = MS_REF), name = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + labs(x = NULL, y = "Stress-annotated gene arrays") + theme_ms()
    save_item("sx_stress", p, 88, 70, "regions/set_stress.tsv")
}

# ---- LFMM supplement: K 3 (Trifolium86MainImpMAF) vs K 5 (Trifolium86MainK5PC5, 5 PCs) ---------------------------------
if (want("lf_")) {
    IC <- "/pipeline/ANALYSIS/Trifolium86Main/imputation_check"
    TREES <- data.table(K = c("K = 3", "K = 5"), project = c("Trifolium86MainImpMAF", "Trifolium86MainK5PC5"),
                        stripes = c("tables_ImpMAF", "tables_Trifolium86MainK5PC5"), k = c(3, 5), tag = c("ImpMAF", "K5PC5"))
    METH <- c("EMMAX", "LFMM", "RDA")
    b <- rbindlist(lapply(seq_len(2), function(i) {
        s <- fread(file.path(IC, TREES$stripes[i], "snp_stripe_imputation.tsv"))
        rbindlist(lapply(METH, function(m) { x <- s[get(paste0("sig_", m)) == 1]
            data.table(K = TREES$K[i], method = m, class = names(FIG3_STRIPE), n = c(sum(x$pattern_n >= 10), sum(x$pattern_n < 10))) }))
    }))
    b[, class := factor(class, names(FIG3_STRIPE))]
    pb <- ggplot(b, aes(K, n, fill = class)) + geom_col(width = 0.7) + facet_wrap(~method, scales = "free_y") +
        scale_fill_manual(values = FIG3_STRIPE, name = NULL) + labs(x = NULL, y = "Significant SNPs") + theme_ms()
    gif <- rbindlist(lapply(seq_len(2), function(i) {
        f <- Sys.glob(sprintf("/pipeline/%s_logs/gea/assoc_lfmm_bio_*.log", TREES$project[i]))
        rbindlist(lapply(f, function(p) { l <- readLines(p)
            x <- regmatches(l, regexpr("\\$ gif\\s*:\\s*num\\s+[0-9.eE+-]+", l))
            if (!length(x)) return(NULL)
            data.table(K = TREES$K[i], predictor = sub("^assoc_lfmm_(.*)\\.log$", "\\1", basename(p)), gif = as.numeric(sub(".*num\\s+", "", x[1]))) }))
    }))
    stopifnot(nrow(gif) == 32)
    gif <- gif[!predictor %in% DUP]
    lev_k3 <- gif[K == "K = 3"][order(gif), predictor]
    gif[, predictor := factor(predictor, lev_k3)]
    pc <- ggplot(gif, aes(gif, predictor, colour = K)) + geom_line(aes(group = predictor), colour = MS_REF, linewidth = 0.3) +
        geom_point(size = 1.2) + geom_vline(xintercept = 1, linetype = "22", colour = MS_REF, linewidth = 0.3) +
        scale_colour_manual(values = FIG3_K, name = NULL) + labs(x = "LFMM inflation factor before genomic control", y = NULL) + theme_ms()
    mk <- fread("/pipeline/ANALYSIS/Trifolium86MainK5/compare/marker_stripes.tsv")
    d <- rbindlist(lapply(seq_len(2), function(i) rbindlist(lapply(METH, function(m)
        mk[, .(K = TREES$K[i], method = m, stripe, share = get(paste0(TREES$tag[i], "_", m)) / n_snps)]))))
    pd <- ggplot(d, aes(method, stripe, fill = share)) + geom_tile(colour = "white") +
        geom_text(aes(label = sprintf("%.0f%%", 100 * share), colour = I(ms_text_on(grDevices::colorRampPalette(MS_SEQUENTIAL)(101)[round(share * 100) + 1]))), size = MS_LAB) +
        facet_wrap(~K) + scale_fill_gradientn(colours = MS_SEQUENTIAL, limits = c(0, 1), labels = scales::percent, name = "Called significant") +
        labs(x = NULL, y = NULL) + theme_ms()
    man <- function(file, col, lab) {
        x <- fread(file, select = c("chr", "pos", "bio_8"))
        x[chr_len, xx := pos + i.off, on = "chr"]
        x[, tint := fifelse(chr %% 2 == 1, col, mix("#FFFFFF", col, 0.45))]
        thr <- 0.05 / nrow(x)
        ggplot(x, aes(xx, nl10(bio_8), colour = I(tint))) + geom_scattermore(pointsize = 1.4, pixels = c(1600, 400)) +
            geom_hline(yintercept = nl10(thr), colour = MS_REF, linewidth = 0.3, linetype = "22") +
            x_chr + labs(x = NULL, y = paste(lab, "-log10 P")) + theme_ms()
    }
    pa <- man(sprintf("/pipeline/%s_results/GEA/tables/methods/LFMM/LFMM_pvalues_K3.tsv", TREES$project[1]), FIG3_K[["K = 3"]], "LFMM K3") /
          man(sprintf("/pipeline/%s_results/GEA/tables/methods/LFMM/LFMM_pvalues_K5.tsv", TREES$project[2]), FIG3_K[["K = 5"]], "LFMM K5") /
          man(file.path(T5, "methods/EMMAX/EMMAX_pvalues_K5.tsv"), MINOU[["red"]], "EMMAX")
    p <- (wrap_elements(pa) | (pb / pc / pd)) + plot_layout(widths = c(1.3, 1)) + plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(size = MS_BASE, face = "bold"))
    save_item("lf_supp", p, 180, 170, "snp_stripe_imputation.tsv (ImpMAF, K5PC5); LFMM logs (gif); compare/marker_stripes.tsv; LFMM/EMMAX p tables (bio_8)")
    fwrite(b, file.path(OUT, "lf_supp_panel_b.tsv"), sep = "\t"); fwrite(gif, file.path(OUT, "lf_supp_panel_c.tsv"), sep = "\t")
    fwrite(d, file.path(OUT, "lf_supp_panel_d.tsv"), sep = "\t")
}

# ---- register -----------------------------------------------------------------------------------------------------------
reg <- rbindlist(reg_items)
idx <- file.path(OUT, "figures.tsv")
if (file.exists(idx)) { old <- fread(idx); reg <- rbind(old[!id %in% reg$id], reg) }
fwrite(reg[order(id)], idx, sep = "\t")
message("figures.tsv: ", nrow(reg), " items")
