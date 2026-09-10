library(ggplot2)
library(dplyr)
library(data.table)
library(qs)

source("/pipeline/scripts/R/utils/theme_clinego.R")

args = commandArgs(trailingOnly=TRUE)
#################################
CLUSTERS = args[1]
K = args[2] %>% as.numeric
OUT_PNG = args[3]
INTER_DIR = args[4]
#################################

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Load clusters
clusters <- fread(CLUSTERS)

# Make wide table
IsrPopKs <- clusters %>%
            dplyr::select(-site) %>%
            melt.data.table(id.vars = 'sample')

# Plot
gStructure <-
  ggplot(data = IsrPopKs, aes(y = value, x = sample, fill = variable)) +
    geom_bar(show.legend = TRUE, stat = "identity", position = "fill") +
    scale_fill_manual(name = "Clusters", values = clinego_cluster_palette(K)) +
    ylab("Proportion of assignment") +
    xlab("Accessions") +
    theme_clinego(base_size = 16) +
    theme(axis.text.x = element_blank())

# Save
ggsave(OUT_PNG, gStructure, width = 3 * 6.4, height = 4.8)
ggsave(OUT_SVG, gStructure, width = 3 * 6.4, height = 4.8,
       device = svglite::svglite, bg = "transparent")
qsave(gStructure, OUT_QS)
