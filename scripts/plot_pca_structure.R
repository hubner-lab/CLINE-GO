library(ggplot2)
library(dplyr)
library(data.table)
library(scatterpie)
library(qs)

source("/pipeline/scripts/R/utils/theme_clinego.R")

args = commandArgs(trailingOnly=TRUE)
#################################
CLUSTERS = args[1]
PROJECTIONS = args[2]  # LEA PCA projections file
EIGENVALUES = args[3]  # LEA PCA eigenvalues file
K = args[4] %>% as.numeric
OUT_PNG = args[5]
INTER_DIR = args[6]
#################################

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Load clusters
clusters <- fread(CLUSTERS)
message(clusters %>% str)

# Load LEA PCA projections (space-separated, no header, rows = samples, cols = PCs)
projections <- fread(PROJECTIONS, header = FALSE)
# Rename columns to PC1, PC2, etc.
colnames(projections) <- paste0('PC', 1:ncol(projections))
message(projections %>% str)

# Load LEA eigenvalues (one value per line)
eigenvalues <- fread(EIGENVALUES, header = FALSE)$V1
message(eigenvalues %>% str)

# Validate dimensions
if (nrow(clusters) != nrow(projections)) {
    stop(paste0("ERROR: Row count mismatch! Clusters has ", nrow(clusters),
                " rows but PCA projections has ", nrow(projections), " rows."))
}

# Combine clusters with PCA projections
df <- cbind(clusters, projections[, 1:min(4, ncol(projections))])

# Calculate variance explained
var.explained <- c((eigenvalues[1] / sum(eigenvalues)) * 100,
                   (eigenvalues[2] / sum(eigenvalues)) * 100)

# Plot PCA with SNMF clusters
gPCA <- ggplot(df, aes(PC1, PC2)) +
  geom_scatterpie(data = df,
                  aes(x = PC1, y = PC2),
                  cols = paste0('C', 1:K),
                  color = 'black', alpha = 0.8) +
  scale_fill_manual(name = "Clusters", values = clinego_cluster_palette(K)) +
  xlab(paste0('PC1 (', round(var.explained[1], 1), '%)')) +
  ylab(paste0('PC2 (', round(var.explained[2], 1), '%)')) +
  theme_clinego()

# Save
ggsave(OUT_PNG, gPCA)
ggsave(OUT_SVG, gPCA, device = svglite::svglite, bg = "transparent")
qsave(gPCA, OUT_QS)

message('INFO: PCA structure plot complete')
