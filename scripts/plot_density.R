#!/usr/bin/env Rscript
# Density plots for climate predictors - combined into single figure
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(qs)

source("/pipeline/scripts/R/utils/factor_order.R")
source("/pipeline/scripts/R/utils/theme_clinego.R")

args = commandArgs(trailingOnly=TRUE)
#################################
CLIMATE = args[1]       # climate site values TSV
PREDICTORS = args[2]    # comma-separated list of bio variables (e.g., "bio_1,bio_2,bio_3")
OUT_PNG = args[3]
INTER_DIR = args[4]
#################################

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- file.path(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Load data
climate <- fread(CLIMATE)
message(paste0('INFO: Loaded data with ', nrow(climate), ' rows and ', ncol(climate), ' columns'))

# Parse predictors — "all" means auto-detect numeric columns (exclude metadata cols)
if (PREDICTORS == "all" || PREDICTORS == "") {
  exclude_cols <- c("site", "sample", "latitude", "longitude", "lat", "lon")
  all_cols <- colnames(climate)
  bio_vars <- all_cols[!tolower(all_cols) %in% tolower(exclude_cols)]
  bio_vars <- bio_vars[vapply(bio_vars, function(x) is.numeric(climate[[x]]), logical(1))]
} else {
  bio_vars <- strsplit(PREDICTORS, ",")[[1]]
}

# Numbered names (bio_N, trait_N) in numeric order, the others after them as given: the order
# every figure uses (factor_order.R)
bio_vars <- bio_vars[numbered_order(bio_vars)]

message(paste0('INFO: Creating density plots for ', length(bio_vars), ' predictors'))

# Plot theme
plot_theme <-
  theme_clinego(base_size = 10) +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        axis.text = element_text(size = 9),
        plot.margin = margin(5, 10, 5, 10))

# Create individual density plots
plot_list <- list()
for (bio in bio_vars) {
  if (!bio %in% colnames(climate)) {
    warning(paste0("Variable '", bio, "' not found in climate data, skipping"))
    next
  }

  gDens <- climate %>%
    ggplot(aes(x = .data[[bio]])) +
    geom_density(fill = CLINEGO_NEUTRAL, alpha = 0.6, color = CLINEGO_NEUTRAL) +
    labs(title = bio, x = "Value", y = "Density") +
    plot_theme

  plot_list[[bio]] <- gDens

  # Save individual plot to intermediate for potential reuse
  qsave(gDens, file.path(INTER_DIR, paste0('density_plot_', bio, '.qs')))
}

message(paste0('INFO: Created ', length(plot_list), ' density plots'))

# Determine grid layout based on number of plots
n_plots <- length(plot_list)
if (n_plots == 0) {
  stop("ERROR: No valid predictors found in climate data")
}

# Calculate optimal grid: prefer wider layouts
if (n_plots <= 3) {
  ncol <- n_plots
  nrow <- 1
} else if (n_plots <= 6) {
  ncol <- 3
  nrow <- ceiling(n_plots / 3)
} else if (n_plots <= 12) {
  ncol <- 4
  nrow <- ceiling(n_plots / 4)
} else {
  ncol <- 5
  nrow <- ceiling(n_plots / 5)
}

# Combine plots with ggarrange
combined_plot <- ggarrange(
  plotlist = plot_list,
  ncol = ncol,
  nrow = nrow,
  align = "hv"
)

# Calculate figure dimensions based on grid
fig_width <- ncol * 3.5
fig_height <- nrow * 3

# Save combined plot
ggsave(OUT_PNG, combined_plot, width = fig_width, height = fig_height, dpi = 300)
ggsave(OUT_SVG, combined_plot, device = svglite::svglite, bg = 'transparent',
       width = fig_width, height = fig_height)
qsave(combined_plot, OUT_QS)

message(paste0('INFO: Saved combined density plot (', ncol, 'x', nrow, ' grid)'))
message('INFO: Complete')
