#!/usr/bin/env Rscript
# =============================================================================
# mvp_slide_figure.R -- ONE conference slide. Not a manuscript figure.
#
# The manuscript version (four panels, full parameter inventory) is
# mvp_design_figure.R. This is the opposite brief: three blocks, one row,
# ~25 words on screen, readable from the back of a room. Everything the
# speaker can say out loud has been taken OFF the slide.
#
# The three blocks are the three things an ecologist needs to accept the
# benchmark, in order:
#   1  it is a real landscape with real populations
#   2  the right answer is known by construction   <- the whole point
#   3  it was repeated across the range of genetic architectures
#
# Deliberately NOT on the slide: demography, pleiotropy, the 225-cell grid,
# migration rates, Fst. Those are spoken, or live in the manuscript figure.
#
# Usage:  Rscript /pipeline/benchmarks/mvp_slide_figure.R
# =============================================================================

suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(cowplot)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
OUT  <- Sys.getenv("FIG_OUT", file.path(ROOT, "benchmarks/mvp_eval/figures_design"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))

MINOU  <- c(teal = "#00798c", red = "#d1495b", amber = "#edae49",
            sage = "#66a182", navy = "#2e4057", grey = "#8d96a3")
INK    <- MINOU[["navy"]]; MUTED <- MINOU[["grey"]]; ACC <- MINOU[["teal"]]
ENV_LO <- "#2c7fb8"; ENV_MID <- "#f7f7f7"; ENV_HI <- MINOU[["red"]]

SEED_MAP   <- "1231288"   # env surface is a landscape property (md5-identical across seeds)
SEED_TRUTH <- "1232353"   # a mid-range moderately polygenic seed: 30 causal loci, legible as ticks

# Slide type sizes. Everything is large; nothing below 9 pt survives a projector.
TS <- list(head = 15, num = 26, sub = 11, tick = 9.5)

blank <- theme_void() + theme(legend.position = "none", plot.margin = margin(2, 6, 2, 6))

# -----------------------------------------------------------------------------
# 1  A WORLD -- two climates over one grid of populations
# -----------------------------------------------------------------------------
e <- fread(file.path(ROOT, "data/mvp", paste0("MVP", SEED_MAP),
                     paste0("MVP", SEED_MAP, "_env_present.tsv")))
e[, idx := as.integer(sub("deme_", "", site))]
e[, `:=`(x = ((idx - 1L) %% 10L) + 1L, y = ((idx - 1L) %/% 10L) + 1L)]

world <- rbind(e[, .(x, y, v = bio_1, panel = "climate 1")],
               e[, .(x, y, v = bio_2, panel = "climate 2")])
world[, panel := factor(panel, levels = c("climate 1", "climate 2"))]

p1 <- ggplot(world, aes(x, y, fill = v)) +
    geom_tile(colour = "white", linewidth = 1.1) +
    facet_wrap(~panel, nrow = 1) +
    scale_fill_gradient2(low = ENV_LO, mid = ENV_MID, high = ENV_HI, midpoint = 0) +
    coord_equal(expand = FALSE) +
    blank +
    theme(strip.text = element_text(size = TS$tick, colour = MUTED,
                                    margin = margin(b = 3)),
          panel.spacing = unit(14, "pt"))

# -----------------------------------------------------------------------------
# 2  KNOWN TRUTH -- the genome, with the real adaptive loci marked
# -----------------------------------------------------------------------------
tru <- fread(file.path(ROOT, "data/mvp", paste0("MVP", SEED_TRUTH), "truth_any.tsv"))
tru[, gpos := (chr - 1) * 50000 + pos]          # 20 LGs x 50 kb = 1 Mb
caus <- tru[category == "causal"]
lgs  <- data.table(lg = 1:20, xmin = (0:19) * 50000, xmax = (1:20) * 50000)
lgs[, fill := fifelse(lg <= 10, "#dfe6ea", "#f2f4f5")]   # QTNs only on LG 1-10

p2 <- ggplot() +
    geom_rect(data = lgs, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = fill),
              colour = "white", linewidth = 0.8) +
    geom_segment(data = caus, aes(x = gpos, xend = gpos, y = 0, yend = 1.85),
                 colour = ACC, linewidth = 1.0) +
    geom_point(data = caus, aes(gpos, 1.85), colour = ACC, size = 1.5) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, 1e6), ylim = c(-0.35, 2.3), expand = FALSE) +
    blank

# -----------------------------------------------------------------------------
# 3  92 WORLDS -- one dot each, in three architecture columns
# -----------------------------------------------------------------------------
s <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"))
LAB <- c(oliogenic = "a few", `mod-polygenic` = "dozens", `highly-polygenic` = "hundreds")
s[, grp := factor(LAB[arch_level], levels = LAB)]
setorder(s, grp, -arm, seed)
dots <- mvp_prim(s)[, .(grp, i = seq_len(.N)), by = grp][, .(grp, i)]
# The two control dots below are hand-placed for a 90-dot grid (dy = 8). A
# larger arm overflows into them, so fail rather than draw a wrong slide.
stopifnot(nrow(dots) == mvp_n_expect())
dots[, `:=`(dx = ((i - 1L) %% 5L) + 1L, dy = ((i - 1L) %/% 5L) + 1L)]
ctrl <- data.table(grp = factor("dozens", levels = LAB), dx = c(2, 4), dy = c(8, 8))

p3 <- ggplot() +
    geom_point(data = dots, aes(dx, dy), colour = ACC, size = 3.1) +
    geom_point(data = ctrl, aes(dx, dy), colour = MINOU[["amber"]], size = 3.1) +
    facet_wrap(~grp, nrow = 1) +
    scale_y_reverse() +
    coord_cartesian(xlim = c(0.3, 5.7), ylim = c(8.9, 0.3), expand = FALSE) +
    blank +
    theme(strip.text = element_text(size = TS$tick, colour = MUTED,
                                    margin = margin(b = 4)),
          panel.spacing = unit(12, "pt"))

# -----------------------------------------------------------------------------
# Assembly: headline / drawing / one number, three times across
# -----------------------------------------------------------------------------
blk <- function(head, art, num, sub, art_h = 1) {
    plot_grid(
        ggdraw() + draw_label(head, fontface = "bold", size = TS$head, colour = INK),
        art,
        ggdraw() + draw_label(num, fontface = "bold", size = TS$num, colour = ACC,
                              y = 0.68) +
                   draw_label(sub, size = TS$sub, colour = MUTED, y = 0.16),
        ncol = 1, rel_heights = c(0.16, art_h, 0.26))
}

arrow_col <- ggdraw() +
    draw_line(x = c(0.28, 0.72), y = c(0.52, 0.52), colour = "grey78", linewidth = 1.1,
              arrow = arrow(length = unit(7, "pt"), type = "closed"))

slide <- plot_grid(
    blk("A simulated world",   p1, "100", "populations, 1,000 individuals", 1.00),
    arrow_col,
    blk("The answer is known", p2, "4 – 548", "truly adaptive genes per world", 0.62),
    arrow_col,
    blk("Repeated",            p3, "92", "worlds, easy to hard genetics", 1.00),
    nrow = 1, rel_widths = c(1, 0.13, 1.05, 0.13, 0.82))

for (ext in c("png", "svg", "pdf")) {
    f <- file.path(OUT, paste0("slide_benchmark.", ext))
    if (ext == "svg") ggsave(f, slide, width = 13.33, height = 5.0, device = svglite::svglite,
                             bg = "white")
    else if (ext == "pdf") ggsave(f, slide, width = 13.33, height = 5.0, device = cairo_pdf)
    else ggsave(f, slide, width = 13.33, height = 5.0, dpi = 300, bg = "white")
}

# Values behind the three numbers, so nothing on the slide is unsourced.
fwrite(data.table(
    block  = c("A simulated world", "The answer is known", "Repeated"),
    number = c("100", "4 - 548", "92"),
    meaning = c("demes per world (10 x 10 lattice), 10 individuals sampled each",
                paste0("range of MAF-0.01 causal loci over the 90 primary replicates; ",
                       "genome drawn from seed ", SEED_TRUTH, " (", nrow(caus), " causal)"),
                "replicates: 30 per architecture level + 2 degenerate controls (amber)"),
    source = "deposit"), file.path(OUT, "slide_benchmark.tsv"), sep = "\t")

message("wrote slide_benchmark.{png,svg,pdf,tsv} to ", OUT)
