#!/usr/bin/env Rscript
# =============================================================================
# mvp_design_figure.R -- a schematic of the MVP simulation benchmark for a
# general biology / ecology audience.
#
# NOT a reproduction of Lotterhos 2023 Fig 1. That figure is a complete parameter
# inventory (and is CC BY-NC-ND, so it may be shown verbatim but not adapted).
# This one answers the four questions a listener actually asks:
#
#   A  What is one simulated world?      lattice, two environments, sampling
#   B  What is a seed / population /     the vocabulary, as a flow
#      garden?
#   C  What was varied?                  the four knobs, as icons, not as text
#   D  What did WE run?                  our 34 of 225 cells, and why this slice
#
# Every panel writes the values it plots to a same-stem .tsv, each carrying a
# `source` column:
#   source = "deposit"    verbatim from the MVP deposit / our cohort manifest.
#                         Panels A, D and the C1 landscape curves. Traceable.
#   source = "schematic"  hand-authored ICONS that illustrate a design axis
#                         without reproducing SLiM's actual values. Panels C2,
#                         C3, C4. Do NOT cite a number off these.
# Reads only tables already on disk; fits nothing, runs no pipeline mode.
#
# Usage:
#   Rscript /pipeline/benchmarks/mvp_design_figure.R
#   FIG_OUT=/pipeline/benchmarks/mvp_eval/figures_design Rscript ...
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
OUT  <- Sys.getenv("FIG_OUT", file.path(ROOT, "benchmarks/mvp_eval/figures_design"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
# Occupied cells expected in panel D. A literal here fired the moment the
# manifest grew past one arm; keep it an env so a block can be inventoried.
N_CELLS <- as.integer(Sys.getenv("MVP_N_CELLS", "34"))

# House palette -- same MINOU hexes as mvp_main_figure.R / mvp_arch_panel.R so
# the design figure sits in the same deck as the results figures.
MINOU <- c(teal = "#00798c", red = "#d1495b", amber = "#edae49",
           sage = "#66a182", navy = "#2e4057", grey = "#8d96a3")
INK   <- MINOU[["navy"]]
MUTED <- MINOU[["grey"]]
# One diverging ramp for BOTH environments, so "cold/low" and "warm/high" mean
# the same thing in every map on the slide.
ENV_LO <- "#2c7fb8"; ENV_MID <- "#f7f7f7"; ENV_HI <- MINOU[["red"]]

SEED_SSMTN <- "1231288"   # any SS-Mtn seed: the env surface is a landscape
SEED_EST   <- "1232568"   # property, verified identical across seeds (md5)

wtsv <- function(dt, stem) fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
save3 <- function(stem, p, w, h) {
    clinego_save_both(file.path(OUT, stem), p, w = w, h = h, dpi = 300)
    ggsave(file.path(OUT, paste0(stem, ".pdf")), p, width = w, height = h, device = cairo_pdf)
}

# Shared look for the schematic panels: no axes, no grid, nothing but the drawing.
theme_schematic <- function() {
    theme_void(base_size = 11) +
        theme(plot.title    = element_text(face = "bold", size = 13, colour = INK,
                                           hjust = 0, margin = margin(b = 2)),
              plot.subtitle = element_text(size = 9.5, colour = MUTED, hjust = 0,
                                           margin = margin(b = 6)),
              legend.position = "none",
              plot.margin = margin(6, 8, 6, 8))
}

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
read_env <- function(seed) {
    f <- file.path(ROOT, "data/mvp", paste0("MVP", seed),
                   paste0("MVP", seed, "_env_present.tsv"))
    e <- fread(f)
    e[, idx := as.integer(sub("deme_", "", site))]
    e[, `:=`(x = ((idx - 1L) %% 10L) + 1L, y = ((idx - 1L) %/% 10L) + 1L)]
    e[]
}
env_mtn <- read_env(SEED_SSMTN)
env_est <- read_env(SEED_EST)

seeds  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"))
params <- fread(file.path(ROOT, "data/mvp/selection/0b-final_params-20220428.txt"))
setnames(params, 1, "row_idx")
summ <- fread(file.path(ROOT, "data/mvp/selection/summary_20220428_20220726.csv"),
              select = c("seed", "meanFst", "cor_PC1_temp", "cor_PC1_sal"))

# Inventory ONE arm at a time, plus the degenerate controls -- they are panel
# D's third status colour, so dropping them would make that branch dead code.
# Default MVP_ARM=primary reproduces the legacy slice exactly: 92 rows over 34
# occupied cells. Unfiltered, the grown manifest (452 rows, 70 cells) would
# silently redraw panel D as a different experiment.
seeds <- seeds[arm %in% c(mvp_arm(), "control_degenerate")]
message(sprintf("design inventory: %s + controls -- %d worlds",
                mvp_arm_label(), nrow(seeds)))
cohort <- merge(seeds, params[, .(seed, demog_name, arch, level)], by = "seed", all.x = TRUE)
stopifnot(nrow(cohort) == nrow(seeds), !anyNA(cohort$level))

N_GARDEN <- length(list.files(file.path(ROOT, "data/mvp", paste0("MVP", SEED_SSMTN), "gardens")))

# =============================================================================
# A -- one simulated world
# =============================================================================
panelA_data <- rbind(
    env_mtn[, .(x, y, value = bio_1, env = "Temperature\nsmooth north-south gradient")],
    env_mtn[, .(x, y, value = bio_2, env = "Second environment (Env2)\nV-shaped ridge: cool edges, warm centre")]
)
panelA_data[, env := factor(env, levels = unique(env))]
panelA_data[, source := "deposit"]
wtsv(panelA_data, "panelA_world")

pA <- ggplot(panelA_data, aes(x, y, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    facet_wrap(~env, nrow = 1) +
    scale_fill_gradient2(low = ENV_LO, mid = ENV_MID, high = ENV_HI, midpoint = 0) +
    coord_equal(expand = FALSE) +
    labs(title = "A   One simulated world",
         subtitle = paste0("100 populations on a 10 x 10 grid. Each square has its own climate. ",
                           "10 individuals sampled per square = 1,000 genomes.")) +
    theme_schematic() +
    theme(strip.text = element_text(size = 9.5, colour = INK, lineheight = 1.15,
                                    margin = margin(b = 4)),
          panel.spacing = unit(18, "pt"))

save3("panelA_world", pA, 8.2, 3.4)

# =============================================================================
# B -- the vocabulary, as a flow
# =============================================================================
steps <- data.table(
    i     = 1:4,
    head  = c("1 SEED", "100 POPULATIONS", paste0(N_GARDEN, " GARDENS"), "GROUND TRUTH"),
    body  = c("one whole simulated\nworld, run from scratch",
              "one grid square each\n10 individuals sampled",
              paste0("\"move this population to\nclimate X\" -- ", N_GARDEN - 12L,
                     " real climates\n+ 12 new ones"),
              paste0(format(100L * N_GARDEN, big.mark = ","),
                     " fitness values\nper world")),
    fill  = c(MINOU[["navy"]], MINOU[["teal"]], MINOU[["amber"]], MINOU[["sage"]])
)
steps[, `:=`(xmin = (i - 1) * 2.55, xmax = (i - 1) * 2.55 + 2.15)]
wtsv(steps[, .(i, head, body = gsub("\n", " ", body), source = "deposit")], "panelB_vocabulary")

pB <- ggplot(steps) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1.55, fill = fill),
              alpha = 0.13, colour = NA) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1.30, ymax = 1.55, fill = fill),
              colour = NA) +
    geom_text(aes(x = (xmin + xmax) / 2, y = 1.425, label = head),
              colour = "white", fontface = "bold", size = 3.5) +
    geom_text(aes(x = (xmin + xmax) / 2, y = 0.62, label = body),
              colour = INK, size = 3.1, lineheight = 1.2) +
    geom_segment(data = steps[i < 4],
                 aes(x = xmax + 0.06, xend = xmax + 0.34, y = 0.775, yend = 0.775),
                 arrow = arrow(length = unit(6, "pt"), type = "closed"),
                 colour = MUTED, linewidth = 0.7) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(-0.1, 9.3), ylim = c(-0.05, 1.62), expand = FALSE) +
    labs(title = "B   What the words mean",
         subtitle = paste0("A garden is a question asked AFTER the simulation, not part of its history. ",
                           "Every population is tested in every garden.")) +
    theme_schematic()

save3("panelB_vocabulary", pB, 8.2, 2.4)

# =============================================================================
# C -- the four knobs
# =============================================================================
# C1 landscape: the second environment across one row of the grid, per landscape
# The two clinal landscapes share an identical environmental surface and differ
# only in migration (m_x 0.49 vs 0.03, from the deposit's params table), so the
# migration rate goes in the facet label -- otherwise the first two panels read
# as an accidental duplicate.
lsX <- rbind(
    data.table(opt = "Estuary\nstrong E-W flow  (m = 0.49)",
               x = 1:10, v = env_est[y == 1][order(x)]$bio_2),
    data.table(opt = "Stepping-stone\nweak, even flow  (m = 0.03)",
               x = 1:10, v = env_est[y == 1][order(x)]$bio_2),
    data.table(opt = "Stepping-stone MOUNTAIN\nweak, even flow  (m = 0.03)",
               x = 1:10, v = env_mtn[y == 1][order(x)]$bio_2)
)
lsX[, opt := factor(opt, levels = unique(opt))]
lsX[, source := "deposit"]
wtsv(lsX, "panelC1_landscape")

pC1 <- ggplot(lsX, aes(x, v)) +
    geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.4) +
    geom_line(colour = INK, linewidth = 0.9) +
    geom_point(aes(fill = v), shape = 21, colour = "white", size = 2.6, stroke = 0.5) +
    facet_wrap(~opt, nrow = 1) +
    scale_fill_gradient2(low = ENV_LO, mid = ENV_MID, high = ENV_HI, midpoint = 0) +
    coord_cartesian(ylim = c(-1.35, 1.35)) +
    labs(title = "1.  Landscape  (3 options)",
         subtitle = "How the second environment changes across space. In the mountain, two distant places share a climate.") +
    theme_schematic() +
    theme(strip.text = element_text(size = 8.5, colour = INK, lineheight = 1.1))

# C2 demography: population size across the grid
# ICONS, not SLiM output: these surfaces resemble the five `Nequal`/`isVariableM`/
# `MIG_breaks` codes without reproducing them. All five levels are drawn so the
# panel title "(5 options)" matches what is on screen.
lat <- function(lab) data.table(opt = lab, x = rep(1:10, 10), y = rep(1:10, each = 10))
demo <- rbind(
    lat("equal size,\neven migration")[, n := 1],
    lat("equal size,\nmigration BARRIERS")[, n := 1],
    lat("size declines\nnorth to south")[, n := y / 10],
    lat("big centre,\nsmall edges")[, n := 1 - (abs(x - 5.5) + abs(y - 5.5)) / 11],
    lat("irregular size AND\nirregular migration")[, n := 0.25 + 0.75 * abs(sin(x * 1.7) * cos(y * 1.3))]
)
demo[, opt := factor(opt, levels = unique(opt))]
# the barrier icon: two hard cuts through the lattice
barrier <- data.table(opt = factor("equal size,\nmigration BARRIERS", levels = levels(demo$opt)),
                      x = c(0.4, 0.4), xend = c(10.6, 10.6), y = c(3.5, 7.5), yend = c(3.5, 7.5))
demo[, source := "schematic"]
wtsv(demo, "panelC2_demography")

pC2 <- ggplot(demo, aes(x, y, size = n)) +
    geom_point(colour = MINOU[["teal"]], alpha = 0.85) +
    geom_segment(data = barrier, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, colour = MINOU[["red"]], linewidth = 0.8) +
    facet_wrap(~opt, nrow = 1) +
    scale_size_continuous(range = c(0.35, 2.2)) +
    coord_equal() +
    labs(title = "2.  Demography  (5 options)",
         subtitle = "How many individuals live where, and whether migration is even. Dot size = population size.") +
    theme_schematic() +
    theme(strip.text = element_text(size = 8.5, colour = INK))

# C3 genic level: same genome, different numbers and sizes of causal loci
set.seed(11)
mk_gene <- function(lab, n, h) data.table(opt = lab, pos = sort(runif(n, 0.02, 0.98)), h = h)
gen <- rbind(
    mk_gene("a few big genes\n(oligogenic)", 5, 1.00),
    mk_gene("dozens, medium\n(moderately polygenic)", 45, 0.45),
    mk_gene("hundreds, tiny\n(highly polygenic)", 320, 0.16)
)
gen[, opt := factor(opt, levels = unique(opt))]
wtsv(gen[, .(opt, pos, h, source = "schematic")], "panelC3_genic")

pC3 <- ggplot(gen) +
    annotate("rect", xmin = 0, xmax = 1, ymin = -0.06, ymax = 0, fill = "grey88") +
    geom_segment(aes(x = pos, xend = pos, y = 0, yend = h),
                 colour = MINOU[["navy"]], linewidth = 0.45, alpha = 0.85) +
    facet_wrap(~opt, nrow = 1) +
    coord_cartesian(ylim = c(-0.1, 1.15), xlim = c(-0.02, 1.02), expand = FALSE) +
    labs(title = "3.  Genetic architecture  (3 options)",
         subtitle = "Bar = the genome. Tick = a gene affecting the trait; tick height = effect size.") +
    theme_schematic() +
    theme(strip.text = element_text(size = 8.5, colour = INK, lineheight = 1.1))

# C4 pleiotropy: gene -> trait wiring
wir <- rbind(
    data.table(opt = "one trait only",       gx = c(1, 2, 3), tx = c(2, 2, 2), tname = "T1"),
    data.table(opt = "two traits,\nseparate genes", gx = c(1, 2, 3), tx = c(1, 1, 3), tname = c("T1", "T1", "T2")),
    data.table(opt = "two traits,\nSHARED genes",   gx = c(1, 2, 2, 3), tx = c(1, 1, 3, 3), tname = c("T1", "T1", "T2", "T2"))
)
wir[, opt := factor(opt, levels = unique(opt))]
trt <- unique(wir[, .(opt, tx, tname)])
gns <- unique(wir[, .(opt, gx)])
wtsv(copy(wir)[, source := "schematic"], "panelC4_pleiotropy")

pC4 <- ggplot() +
    geom_segment(data = wir, aes(x = gx, xend = tx, y = 0, yend = 1),
                 colour = MUTED, linewidth = 0.55) +
    geom_point(data = gns, aes(gx, 0), size = 4, colour = MINOU[["amber"]]) +
    geom_point(data = trt, aes(tx, 1), size = 6.5, shape = 21,
               fill = MINOU[["sage"]], colour = "white", stroke = 0.7) +
    geom_text(data = trt, aes(tx, 1, label = tname), colour = "white",
              size = 2.6, fontface = "bold") +
    facet_wrap(~opt, nrow = 1) +
    coord_cartesian(xlim = c(0.4, 3.6), ylim = c(-0.35, 1.35), expand = FALSE) +
    labs(title = "4.  Pleiotropy  (5 options: these 3 wirings, the two 2-trait ones run under equal or unequal selection)",
         subtitle = "Dot below = a gene, circle above = a trait. Pleiotropy means one gene does two jobs at once.") +
    theme_schematic() +
    theme(strip.text = element_text(size = 8.5, colour = INK, lineheight = 1.1))

pC <- plot_grid(pC1, pC2, pC3, pC4, ncol = 1, align = "v", axis = "lr",
                rel_heights = c(1, 1, 1, 1))
pC <- plot_grid(
    ggdraw() + draw_label("C   The four things the study varied", fontface = "bold",
                          size = 13, colour = INK, x = 0.005, hjust = 0) +
        draw_label(paste0("3 landscapes  x  5 demographies  x  3 architectures  x  5 pleiotropy levels",
                          "  =  225 combinations,  each run 10 times  =  2,250 worlds"),
                   size = 9.5, colour = MUTED, x = 0.005, hjust = 0, y = 0.18),
    pC, ncol = 1, rel_heights = c(0.10, 1))

save3("panelC_knobs", pC, 8.2, 8.6)
for (nm in c("C1", "C2", "C3", "C4"))
    save3(paste0("panelC_", nm), get(paste0("p", nm)), 8.2, 2.2)

# =============================================================================
# D -- our slice, and why
# =============================================================================
grid_all <- unique(params[, .(demog_name, arch, demog_level, demog_level_sub,
                              arch_level, arch_level_sub)])
LAND_ORD  <- c("Est-Clines", "SS-Clines", "SS-Mtn")
GENIC_ORD <- c("oliogenic", "mod-polygenic", "highly-polygenic")
grid_all[, `:=`(demog_level = factor(demog_level, levels = LAND_ORD),
                arch_level  = factor(arch_level,  levels = GENIC_ORD))]
setorder(grid_all, demog_level, demog_level_sub, arch_level, arch_level_sub)
grid_all[, row := .GRP, by = .(demog_level, demog_level_sub)]
grid_all[, col := .GRP, by = .(arch_level, arch_level_sub)]

used <- cohort[, .(n = .N, primary = sum(arm == mvp_arm())), by = .(demog_name, arch)]
grid_all <- merge(grid_all, used, by = c("demog_name", "arch"), all.x = TRUE)
grid_all[is.na(n), `:=`(n = 0L, primary = 0L)]
grid_all[, status := fifelse(n == 0L, "not run",
                     fifelse(primary > 0L, "our replicates", "degenerate control"))]
wtsv(grid_all[, .(demog_level, demog_level_sub, arch_level, arch_level_sub,
                  row, col, n, primary, status, source = "deposit")], "panelD_slice")

stopifnot(sum(grid_all$n) == nrow(cohort), nrow(grid_all[n > 0]) == N_CELLS)

band_rows <- grid_all[, .(ymin = min(row) - 0.5, ymax = max(row) + 0.5), by = demog_level]
band_cols <- grid_all[, .(xmin = min(col) - 0.5, xmax = max(col) + 0.5), by = arch_level]

pD1 <- ggplot(grid_all) +
    geom_tile(aes(col, row, fill = status), colour = "white", linewidth = 0.9) +
    geom_rect(data = band_rows, aes(xmin = 0.5, xmax = 15.5, ymin = ymin, ymax = ymax),
              fill = NA, colour = INK, linewidth = 0.55) +
    geom_rect(data = band_cols, aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = 15.5),
              fill = NA, colour = INK, linewidth = 0.55) +
    scale_fill_manual(values = c("not run" = "grey92",
                                 "our replicates" = MINOU[["teal"]],
                                 "degenerate control" = MINOU[["amber"]])) +
    scale_y_reverse(breaks = band_rows[, (ymin + ymax) / 2], labels = band_rows$demog_level) +
    scale_x_continuous(position = "top",
                       breaks = band_cols[, (xmin + xmax) / 2], labels = band_cols$arch_level) +
    coord_equal(expand = FALSE) +
    labs(title = "D   What we ran",
         subtitle = paste0("Each small square is one combination, run 10 times. We used ",
                           nrow(grid_all[n > 0]), " of 225 combinations = ",
                           nrow(cohort), " of 2,250 worlds."),
         x = NULL, y = NULL) +
    theme_schematic() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 9, colour = INK),
          axis.text.x = element_text(size = 8.5, colour = INK, margin = margin(b = 3)),
          axis.text.y = element_text(size = 8.5, colour = INK, angle = 90, hjust = 0.5,
                                     margin = margin(r = 3)))

# D2 -- why the mountain landscape
conf <- merge(params[, .(seed, demog_level)], summ, by = "seed")
conf <- conf[, .(r2 = median(cor_PC1_temp^2, na.rm = TRUE)), by = demog_level]
conf[, demog_level := factor(demog_level, levels = LAND_ORD)]
conf[, lab := c("Estuary", "Stepping-stone", "Stepping-stone MOUNTAIN")[as.integer(demog_level)]]
conf[, chosen := demog_level == "SS-Mtn"]
wtsv(copy(conf)[, source := "deposit"], "panelD2_confounding")

pD2 <- ggplot(conf, aes(reorder(lab, -r2), r2, fill = chosen)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = sprintf("%.2f", r2)), hjust = -0.18, size = 3.2, colour = INK) +
    scale_fill_manual(values = c(`TRUE` = MINOU[["teal"]], `FALSE` = "grey85")) +
    coord_flip(ylim = c(0, 1.13), expand = FALSE) +
    labs(title = "Why the mountain?",
         subtitle = paste0("How much of the main genetic axis is just temperature.\n",
                           "Near 1 = geography and climate are the same variable, so there is\n",
                           "nothing to disentangle. We chose the only landscape where they differ."),
         x = NULL, y = NULL) +
    theme_schematic() +
    theme(axis.text.y = element_text(size = 9, colour = INK, hjust = 1,
                                     margin = margin(r = 4)))

save3("panelD_slice", pD1, 6.2, 6.4)
save3("panelD2_confounding", pD2, 5.4, 2.6)

# =============================================================================
# Composites
# =============================================================================
talk <- plot_grid(pA, pB, ncol = 1, rel_heights = c(1.32, 1))
talk <- plot_grid(talk, plot_grid(pD1, pD2, ncol = 1, rel_heights = c(1.9, 1)),
                  nrow = 1, rel_widths = c(1.25, 1))
save3("composite_talk", talk, 15.5, 8.0)

full <- plot_grid(
    plot_grid(pA, pB, ncol = 1, rel_heights = c(1.35, 1)),
    pC,
    plot_grid(pD1, pD2, ncol = 1, rel_heights = c(1.9, 1)),
    nrow = 1, rel_widths = c(1.05, 1.05, 0.95))
save3("composite_full", full, 20.5, 9.2)

message("wrote ", length(list.files(OUT)), " files to ", OUT)
