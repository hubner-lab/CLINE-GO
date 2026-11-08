#!/usr/bin/env Rscript
# =============================================================================
# mvp_land_data_s1_trifolium.R -- Job F3 and F4 of docs/manuscript/lab/data_s1_landscapes_LAB_PROMPT.md (2026-11-08):
# the two Trifolium merges of the story restructure, drawn by the scripts that drew the published members.
#   F3 s_structure_details   a = s_pca34 (plot_pca_structure_ms.R, PC3 x PC4, 88 x 80 mm)
#                            b = s_ancestry_by_site (plot_structure_ms.R, split layout, 180 x 45 mm)
#                            arguments as manuscript/trifolium/Trifolium86MainK5ERpc5/fig3/run_structure.sh
#   F4 s_random_sets_climate a = s_random_sets (closeout copy build_paper_pack_d150.R, SUPP2_MS block, 88 x 110 mm)
#                            b = s_random_climate (closeout copy plot_random_climate_ms_d150.R, 180 x 82 mm)
#                            environment as work/closeout_lab/results/scripts/tools/run_closeout.sh steps s17 / s18
#                            (the MINOU block: AN = the 150 kb offset run's offset_minou, its tables)
# Each script runs unchanged through sys.source() in its own environment (commandArgs() answered from that
# environment), writes its member files into CHECK_OUT, where their SVGs are compared with the published ones, and
# leaves its plot object (gPCA, gStructure, p7, fig), which is composed here: a above b, a at its own width at the
# left, panel letters a and b only, the members at their published sizes. Nothing is refitted or recomputed.
# Usage (inside the image, the MINOU variables of run_closeout.sh set):
#   FIG_OUT=<dir> CHECK_OUT=<dir> Rscript benchmarks/mvp_land_data_s1_trifolium.R
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
ROOT  <- "/pipeline"
OUT   <- Sys.getenv("FIG_OUT");   stopifnot(nzchar(OUT))
CHECK <- Sys.getenv("CHECK_OUT"); stopifnot(nzchar(CHECK))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(CHECK, recursive = TRUE, showWarnings = FALSE)
MS <- new.env(); sys.source(file.path(ROOT, "benchmarks/ms_style.R"), envir = MS)
FAM <- MS$theme_ms()$text$family

run_script <- function(file, args = character()) {
    e <- new.env(parent = globalenv())
    e$commandArgs <- function(trailingOnly = FALSE) if (trailingOnly) args else c("R", "--args", args)
    message("== ", file, " ", paste(args, collapse = " "))
    sys.source(file, envir = e)
    e
}
HEAD_MM <- 4                                              # a strip above each panel for its letter, clear of the member's text
compose <- function(a, b, wa, ha, hb) {                  # a above b, a at its own width at the left of a 180 mm page
    top <- plot_grid(a, NULL, nrow = 1, rel_widths = c(wa, 180 - wa))
    with_head <- function(p, h) plot_grid(NULL, p, ncol = 1, rel_heights = c(HEAD_MM, h))
    plot_grid(with_head(top, ha), with_head(b, hb), ncol = 1, rel_heights = c(HEAD_MM + ha, HEAD_MM + hb),
              labels = c("a", "b"), label_size = 8, label_fontface = "bold", label_fontfamily = FAM)
}

# ---- F3 -------------------------------------------------------------------------------------------------------------
F3  <- file.path(ROOT, "manuscript/trifolium/Trifolium86MainK5ERpc5/fig3")
RES <- file.path(ROOT, "Trifolium86MainK5ERpc5_results")
PCA <- file.path(RES, "_work/maf0.1_miss0.3_smiss0.8_rel0.999/ld0.2_win100_step10/Trifolium86_chr1-6.pca/Trifolium86_chr1-6")
CL  <- file.path(ROOT, "ANALYSIS/Trifolium86MainK5ERpc5/fig3/inputs/clusters_K5.short.tsv")
e3a <- run_script(file.path(F3, "plot_pca_structure_ms.R"),
                  c(CL, paste0(PCA, ".projections"), paste0(PCA, ".eigenvalues"), "5", file.path(CHECK, "s_pca34"), "88x80", "3", "4", "F"))
e3b <- run_script(file.path(F3, "plot_structure_ms.R"), c(CL, "5", file.path(CHECK, "s_ancestry_by_site"), "180x45", "split"))
MS$ms_save(file.path(OUT, "s_structure_details"), compose(e3a$gPCA, e3b$gStructure, 88, 80, 45), 180, 133)
message("s_structure_details 180 x 133 mm")

# ---- F4 -------------------------------------------------------------------------------------------------------------
CO <- file.path(ROOT, "work/closeout_lab/results/scripts/patched")
Sys.setenv(CLOSEOUT_OUT = file.path(CHECK, "s17"), SOFT_LOG = file.path(CHECK, "s17_soft_asserts.tsv"))
e4a <- run_script(file.path(CO, "build_paper_pack_d150.R"))
Sys.setenv(CLOSEOUT_OUT = file.path(CHECK, "s18"), SOFT_LOG = file.path(CHECK, "s18_soft_asserts.tsv"))
e4b <- run_script(file.path(CO, "plot_random_climate_ms_d150.R"))
MS$ms_save(file.path(OUT, "s_random_sets_climate"), compose(e4a$p7, e4b$fig, 88, 110, 82), 180, 200)
message("s_random_sets_climate 180 x 200 mm")
message("done")
