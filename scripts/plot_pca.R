library(LEA)
library(ggplot2)
library(dplyr)
library(data.table)
library(qs)

source("/pipeline/scripts/R/utils/theme_clinego.R")

args = commandArgs(trailingOnly=TRUE)
#################################
LFMM = args[1]
SAMPLES = args[2]
OUT_PCA_PNG = args[3]
OUT_TRACY_PNG = args[4]
INTER_DIR = args[5]
#################################

# Derive SVG and QS paths from PNG paths
OUT_PCA_SVG <- sub("\\.png$", ".svg", OUT_PCA_PNG)
OUT_PCA_QS  <- file.path(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PCA_PNG)))
OUT_TRACY_SVG <- sub("\\.png$", ".svg", OUT_TRACY_PNG)
OUT_TRACY_QS  <- file.path(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_TRACY_PNG)))

# Load samples info
samples.df <- fread(SAMPLES,
                    colClasses = c("site" = "character", 'sample' = 'character')) %>%
              dplyr::select(sample, site)

# Run PCA
pc <- LEA::pca(LFMM, scale = TRUE)

# Preprocess for plotting
			message(pc$projections %>% str)
PCA <- as.data.frame(pc$projections) %>%
       setNames(paste0('PC', 1:ncol(.)))
isrPCA <- cbind(samples.df, PCA)
isrPCA$site <- factor(isrPCA$site)

# Each axis's share of the summed eigenvalues, in %: the variance_fraction column of
# tracy_widom.tsv (576d53: the labels said only PC1/PC2, which hid whether PC1 holds 40% or 2.7%).
var.explained <- c((pc$eigenvalues[1,] / sum(pc$eigenvalues)) * 100,
                   (pc$eigenvalues[2,] / sum(pc$eigenvalues)) * 100)
pct_label <- function(x) if (x >= 1) sprintf("%.1f%%", x) else sprintf("%.2g%%", x)

# PCA plot
gPCA <- ggplot(isrPCA, aes(x = PC1, y = PC2,
                           label = site,
                           color = site)) +
        geom_label() +
        scale_color_manual(name = "Site",
                           values = clinego_cluster_palette(nlevels(isrPCA$site))) +
        xlab(paste0("PC1 (", pct_label(var.explained[1]), ")")) +
        ylab(paste0("PC2 (", pct_label(var.explained[2]), ")")) +
        guides(color = guide_legend(ncol=2)) +
        theme_clinego()

ggsave(OUT_PCA_PNG, gPCA)
ggsave(OUT_PCA_SVG, gPCA, device = svglite::svglite, bg = 'transparent')
qsave(gPCA, OUT_PCA_QS)

# Tracy-Widom test: the variance each axis explains, the axes LEA's test calls significant
# (p < 0.05) filled (fff847: the plot showed the fraction under a "percentage" label and marked
# no axis). The numbers are tracy_widom.tsv's (6bf499). The last axis's statistic is undefined
# (LEA writes -NAN), so a missing p counts as not significant.
tw <- tracy.widom(pc)
tw_df <- data.frame(axis        = tw$N,
                    variance    = tw$percentage * 100,
                    significant = factor(ifelse(!is.na(tw$pvalues) & tw$pvalues < 0.05,
                                                "p < 0.05", "p ≥ 0.05"),
                                         levels = c("p < 0.05", "p ≥ 0.05")))
tw_df <- tw_df[tw_df$axis <= 25, ]

gTW <- ggplot(data = tw_df, aes(x = axis, y = variance, shape = significant)) +
       geom_point(color = CLINEGO_NEUTRAL, size = 2.2, stroke = 0.8) +
       scale_shape_manual(name = "Tracy–Widom", values = c(16, 1), drop = FALSE) +
       scale_x_continuous(breaks = scales::breaks_pretty()) +
       labs(title = "Tracy–Widom test", x = "PCA axis", y = "Variance explained (%)") +
       theme_clinego() +
       theme(plot.title = element_text(hjust = 0.5))

ggsave(OUT_TRACY_PNG, gTW)
ggsave(OUT_TRACY_SVG, gTW, device = svglite::svglite, bg = 'transparent')
qsave(gTW, OUT_TRACY_QS)
