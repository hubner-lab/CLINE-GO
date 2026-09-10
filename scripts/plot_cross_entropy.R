library(LEA)
library(ggplot2)
library(dplyr)
library(tibble)
library(qs)

source("/pipeline/scripts/R/utils/theme_clinego.R")

args = commandArgs(trailingOnly=TRUE)
#################################
SNMF_PROJECT = args[1]
K_START = args[2] %>% as.numeric
K_END = args[3] %>% as.numeric
OUT_PNG = args[4]
INTER_DIR = args[5]
#################################

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Load SNMF project
snmf_res <- load.snmfProject(SNMF_PROJECT)

# Calculate cross-entropy for each K
message(paste0("INFO: Run Cross Entropy test on K", K_START, '-', K_END))

gCrossEntropy <-
  sapply(K_START:K_END, function(K){ min(cross.entropy(snmf_res, K)) }) %>%
  setNames(K_START:K_END) %>%
  enframe %>%
  dplyr::arrange(as.numeric(name)) %>%
  ggplot(aes(x = factor(paste0('K', name),
                         levels = paste0('K', K_START:K_END)), y = value)) +
    geom_point(color = CLINEGO_NEUTRAL) +
    theme_clinego() +
    labs(x = 'Number of ancestral populations',
         y = 'Cross-entropy')

# Save outputs
ggsave(OUT_PNG, gCrossEntropy)
ggsave(OUT_SVG, gCrossEntropy, device = svglite::svglite, bg = 'transparent')
qsave(gCrossEntropy, OUT_QS)
