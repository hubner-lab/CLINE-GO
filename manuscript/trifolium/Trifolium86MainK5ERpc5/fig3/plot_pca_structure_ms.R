#!/usr/bin/env Rscript
# DERIVED COPY of scripts/plot_pca_structure.R (pipeline, 2026-10-16 state) for the Trifolium paper figures (plan
# ~/.claude/plans/fancy-scribbling-garden.md, Phase 1). Lines marked `# DERIVED` differ; each is an entry in
# manuscript_fig3/pipeline_deviations.md awaiting the user's approval:
#   D1 pie colours = minou FIG3_CLUSTERS instead of clinego_cluster_palette(K)
#   D2 any PC pair (args 7-8), variance explained computed for that pair with the pipeline's own formula
#   D3 optional site labels at the per-site centroid of the pies (arg 9), site = short name from site_names.tsv
#      (the CLUSTERS input is fig3/inputs/clusters_K5.short.tsv)
#   D4 theme_ms() + ms_save() at a fixed mm size (7/6 pt, editable SVG) instead of theme_clinego() + ggsave defaults;
#      no .qs (args 5-6 are an output stem and a width/height in mm)
#   D12 axis labels with one decimal, sprintf("%.1f") instead of round(x, 1) ("PC2 (8.0%)", not "8%"; applies N2)
#   D20 optional arg 10 = legend position; absent = theme_ms() default (bottom). "none" for Fig 3a, whose cluster
#      colours the piemap legend beside it already keys (user 2026-10-18)
library(ggplot2)
library(dplyr)
library(data.table)
library(scatterpie)
library(ggrepel)  # DERIVED (D3)

source("/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/fig3_style.R")  # DERIVED (D1, D4)

args = commandArgs(trailingOnly=TRUE)
#################################
CLUSTERS = args[1]
PROJECTIONS = args[2]  # LEA PCA projections file
EIGENVALUES = args[3]  # LEA PCA eigenvalues file
K = args[4] %>% as.numeric
OUT_STEM = args[5]                                                    # DERIVED (D4)
SIZE_MM = as.numeric(strsplit(args[6], "x")[[1]])                     # DERIVED (D4) "88x80"
PC_X = as.integer(args[7]); PC_Y = as.integer(args[8])                # DERIVED (D2)
SITE_LABELS = toupper(substr(args[9], 1, 1)) == "T"                   # DERIVED (D3)
LEGEND_POS = if (length(args) >= 10) args[10] else NA                 # DERIVED (D20)
#################################

# Load clusters
clusters <- fread(CLUSTERS, colClasses = list(character = c("sample", "site")))  # DERIVED (colClasses, CLAUDE.md)
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
df <- cbind(clusters, projections[, 1:max(4, PC_X, PC_Y)])            # DERIVED (D2)
PX <- paste0('PC', PC_X); PY <- paste0('PC', PC_Y)                     # DERIVED (D2)

# Calculate variance explained
var.explained <- c((eigenvalues[PC_X] / sum(eigenvalues)) * 100,       # DERIVED (D2)
                   (eigenvalues[PC_Y] / sum(eigenvalues)) * 100)       # DERIVED (D2)

# Plot PCA with SNMF clusters
gPCA <- ggplot(df, aes(.data[[PX]], .data[[PY]])) +                    # DERIVED (D2)
  geom_scatterpie(data = df,
                  aes(x = .data[[PX]], y = .data[[PY]]),               # DERIVED (D2)
                  cols = paste0('C', 1:K),
                  color = 'black', alpha = 0.8,
                  linewidth = 0.1) +                                   # DERIVED (D4) hairline pie borders at mm size
  scale_fill_manual(name = "Clusters", values = FIG3_CLUSTERS[1:K]) +  # DERIVED (D1)
  xlab(paste0(PX, ' (', sprintf('%.1f', var.explained[1]), '%)')) +   # DERIVED (D12)
  ylab(paste0(PY, ' (', sprintf('%.1f', var.explained[2]), '%)')) +   # DERIVED (D12)
  coord_equal() +                                                      # DERIVED (D4) pies stay round at any panel size
  theme_ms()                                                           # DERIVED (D4)

if (SITE_LABELS) {                                                     # DERIVED (D3)
  ctr <- df[, .(x = mean(get(PX)), y = mean(get(PY))), by = site]
  gPCA <- gPCA + geom_label_repel(data = ctr, aes(x, y, label = site), inherit.aes = FALSE,
                                  size = MS_LAB, colour = MS_INK, fill = "white", label.size = 0.15,
                                  label.padding = unit(0.6, "mm"), segment.size = 0.2,
                                  min.segment.length = 0, box.padding = 0.4, seed = 1, max.overlaps = Inf)
}

if (!is.na(LEGEND_POS)) gPCA <- gPCA + theme(legend.position = LEGEND_POS)  # DERIVED (D20)

# Save
ms_save(OUT_STEM, gPCA, SIZE_MM[1], SIZE_MM[2])                        # DERIVED (D4)

message('INFO: PCA structure plot complete')
