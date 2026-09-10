library(LEA)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(qs)
library(scattermore)

source("/pipeline/scripts/R/utils/theme_clinego.R")

args = commandArgs(trailingOnly=TRUE)
#################################
SNMF_PROJECT = args[1]
K = args[2] %>% as.numeric
PLOIDY = args[3] %>% as.numeric
OUT_PNG = args[4]
INTER_DIR = args[5]
#################################

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Load SNMF project
snmf_res <- load.snmfProject(SNMF_PROJECT)

# Population differentiation tests
p <- snmf.pvalues(snmf_res,
                  entropy = TRUE,
                  ploidy = PLOIDY,
                  K = K)
pvalues.df <- data.frame(pvalues = p$pvalues)
pvalues.df$idx <- 1:nrow(pvalues.df)
pvalues.df$log10p <- -log10(pvalues.df$pvalues)

n_snps <- nrow(pvalues.df)
message(paste0('INFO: Using scattermore for rendering (', n_snps, ' SNPs)'))

# Histogram plot
gHist <-
  ggplot(pvalues.df, aes(x = pvalues)) +
    geom_histogram(aes(y = after_stat(density)), color = 'black', fill = CLINEGO_NEUTRAL) +
    geom_density(alpha = 0.2, color = CLINEGO_THRESHOLD, fill = CLINEGO_THRESHOLD) +
    theme_clinego(base_size = 16)

# Manhattan-style plot
gPval <- ggplot(pvalues.df, aes(x = idx, y = log10p)) +
    # pixels should match output: 12.8in x 4.8in (half height) @ 300dpi = ~3840x1440
    geom_scattermore(color = CLINEGO_NEUTRAL, pointsize = 12, pixels = c(3840, 1440),
                     interpolate = FALSE)

gPval <- gPval +
    labs(x = 'Index', y = expression(-log[10](p-value))) +
    theme_clinego(base_size = 16)

# Combine plots
gGrid <- ggarrange(gHist, gPval, nrow = 2)

# Save
ggsave(OUT_PNG, gGrid, width = 2 * 6.4, height = 9.6)
ggsave(OUT_SVG, gGrid, width = 2 * 6.4, height = 9.6,
       device = svglite::svglite, bg = "transparent")
qsave(gGrid, OUT_QS)
