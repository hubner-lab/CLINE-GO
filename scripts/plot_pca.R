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
OUT_PCA_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PCA_PNG)))
OUT_TRACY_SVG <- sub("\\.png$", ".svg", OUT_TRACY_PNG)
OUT_TRACY_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_TRACY_PNG)))

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

var.explained <- c((pc$eigenvalues[1,] / sum(pc$eigenvalues)) * 100,
                   (pc$eigenvalues[2,] / sum(pc$eigenvalues)) * 100)

# PCA plot
gPCA <- ggplot(isrPCA, aes(x = PC1, y = PC2,
                           label = site,
                           color = site)) +
        geom_label() +
        scale_color_manual(name = "Site",
                           values = clinego_cluster_palette(nlevels(isrPCA$site))) +
        xlab("PC1") +
        ylab("PC2") +
        guides(color = guide_legend(ncol=2)) +
        theme_clinego()

ggsave(OUT_PCA_PNG, gPCA)
ggsave(OUT_PCA_SVG, gPCA, device = svglite::svglite, bg = 'transparent')
qsave(gPCA, OUT_PCA_QS)

# Tracy-Widom test
tw <- tracy.widom(pc)

gTW <- ggplot(data = tw, aes(x = N, y = percentage)) +
       geom_point(color = CLINEGO_NEUTRAL) +
       xlim(1, 25) +
       ggtitle("Tracy Widom") +
       theme_clinego() +
       theme(plot.title = element_text(hjust = 0.5))

ggsave(OUT_TRACY_PNG, gTW)
ggsave(OUT_TRACY_SVG, gTW, device = svglite::svglite, bg = 'transparent')
qsave(gTW, OUT_TRACY_QS)
