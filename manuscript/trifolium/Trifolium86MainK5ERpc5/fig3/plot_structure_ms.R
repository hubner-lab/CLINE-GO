#!/usr/bin/env Rscript
# DERIVED COPY of scripts/plot_structure.R (pipeline, 2026-10-16 state) for the Trifolium paper figures (plan
# ~/.claude/plans/fancy-scribbling-garden.md, Phase 1). Lines marked `# DERIVED` differ; each is an entry in
# manuscript_fig3/pipeline_deviations.md awaiting the user's approval:
#   D1 bar colours = minou FIG3_CLUSTERS instead of clinego_cluster_palette(K)
#   D4 theme_ms() + ms_save() at a fixed mm size instead of theme_clinego(base_size = 16) + ggsave 19.2 x 4.8 in; no .qs
#   D5 LAYOUT = "split" adds the per-site split (user 2026-10-16: wanted as a pipeline option, both drawn, app toggle):
#      one facet per site, plants inside a site ordered by their dominant cluster then its share, short site name
#      under each facet. "split": sites grouped by the cluster most of their plants belong to (the draft st_admix
#      order the user kept); "split_site": sites in short-name order 1..11 (south to north), offered as a variant.
#      LAYOUT = "pipeline" keeps the pipeline's single bar in input order.
library(ggplot2)
library(dplyr)
library(data.table)

source("/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/fig3_style.R")  # DERIVED (D1, D4)

args = commandArgs(trailingOnly=TRUE)
#################################
CLUSTERS = args[1]
K = args[2] %>% as.numeric
OUT_STEM = args[3]                                                    # DERIVED (D4)
SIZE_MM = as.numeric(strsplit(args[4], "x")[[1]])                     # DERIVED (D4) "180x45"
LAYOUT = args[5]                                                      # DERIVED (D5) pipeline | split | split_site
#################################

# Load clusters
clusters <- fread(CLUSTERS, colClasses = list(character = c("sample", "site")))  # DERIVED (colClasses, CLAUDE.md)

if (LAYOUT %in% c("split", "split_site")) {                         # DERIVED (D5)
  QC <- paste0("C", 1:K)
  clusters[, dom := max.col(as.matrix(.SD), ties.method = "first"), .SDcols = QC]
  clusters[, qmax := do.call(pmax, .SD), .SDcols = QC]
  site_lev <- if (LAYOUT == "split") {
    # sites grouped by the cluster most of their plants belong to (the draft st_admix rule the user kept)
    clusters[, .N, by = .(site, dom)][order(dom, -N)][!duplicated(site), site]
  } else as.character(sort(as.integer(unique(clusters$site))))   # split_site: short-name order 1..11
  stopifnot(length(site_lev) == uniqueN(clusters$site))
  smp_ord <- clusters[order(factor(site, site_lev), dom, -qmax), sample]
  IsrPopKs <- melt(clusters[, c("sample", "site", QC), with = FALSE], id.vars = c("sample", "site"))
  IsrPopKs[, `:=`(sample = factor(sample, smp_ord), site = factor(site, site_lev))]
} else {
# Make wide table
IsrPopKs <- clusters %>%
            dplyr::select(-site) %>%
            melt.data.table(id.vars = 'sample')
}

# Plot
gStructure <-
  ggplot(data = IsrPopKs, aes(y = value, x = sample, fill = variable)) +
    geom_bar(show.legend = TRUE, stat = "identity", position = "fill", width = 1) +  # DERIVED (D4) no gaps at mm size
    scale_fill_manual(name = "Clusters", values = FIG3_CLUSTERS[1:K]) +              # DERIVED (D1)
    scale_y_continuous(expand = c(0, 0)) +                                           # DERIVED (D4)
    ylab("Proportion of assignment") +
    xlab("Accessions") +
    theme_ms() +                                                                      # DERIVED (D4)
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())             # DERIVED (D4) ticks off with labels

if (LAYOUT %in% c("split", "split_site")) {                         # DERIVED (D5)
  gStructure <- gStructure +
    facet_grid(~site, scales = "free_x", space = "free_x", switch = "x") +
    xlab("Site") +
    theme(panel.spacing = unit(0.6, "mm"), strip.placement = "outside", axis.line.x = element_blank())
}

# Save
ms_save(OUT_STEM, gStructure, SIZE_MM[1], SIZE_MM[2])                 # DERIVED (D4)
