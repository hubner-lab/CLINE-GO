library(qs)
library(ggplot2)
library(ggpubr)
library(tibble)
library(forcats)
library(gradientForest)
library(stringr)
library(svglite)
args = commandArgs(trailingOnly=TRUE)
###############
GF_ADAPTIVE_PATH = args[1]
GF_RANDOM_PATH = args[2]    # path or "NULL"
OUT_PNG = args[3]
INTER_DIR = args[4]
###############

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

message('INFO: Plotting Overall Importance')

plot_theme <- theme_classic(base_size = 12, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 17, face = 'bold'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 14))

gf <- qread(GF_ADAPTIVE_PATH)
has_random <- GF_RANDOM_PATH != 'NULL'

# Adaptive importance
imp <- importance(gf, type = 'Weighted') %>%
  enframe() %>%
  dplyr::mutate(name = as.factor(name))

max_val <- max(imp$value) + 0.005

gAdapt <- ggplot(imp, aes(y = fct_reorder(name, value), x = value)) +
  geom_bar(stat = 'identity', fill = 'red', color = 'black') +
  plot_theme +
  labs(x = expression(paste("R"^2, " weighted importance")),
       title = 'Adaptive') +
  xlim(c(0, max_val))

if (has_random) {
  gf_random <- qread(GF_RANDOM_PATH)
  # gradient_forest_model.R writes a list sentinel (status empty_forest) when no random SNP
  # had a positive R^2 — the null explaining nothing. Show that instead of erroring.
  if (!inherits(gf_random, 'gradientForest')) {
    message('INFO: random model is a sentinel (', gf_random$status %||% 'unknown', '): ', gf_random$reason %||% '')
    has_random <- FALSE
    gNeutral_note <- ggplot() + theme_void() +
      annotate('text', x = 0.5, y = 0.5, size = 3.5,
               label = paste0('Neutral model empty:\nno random SNP with R\u00b2 > 0')) +
      labs(title = 'Neutral')
    gImp <- ggarrange(gAdapt, gNeutral_note, ncol = 2)
  }
}
if (has_random) {
  imp_random <- importance(gf_random, type = 'Weighted') %>%
    enframe() %>%
    dplyr::mutate(name = as.factor(name))

  gNeutral <- ggplot(imp_random, aes(y = fct_reorder(name, value), x = value)) +
    geom_bar(stat = 'identity', fill = 'blue', color = 'black') +
    plot_theme +
    labs(x = expression(paste("R"^2, " weighted importance")),
         title = 'Neutral') +
    xlim(c(0, max_val))

  gImp <- ggarrange(gAdapt, gNeutral, ncol = 2)
} else if (!exists('gImp')) {
  gImp <- gAdapt
}

ggsave(OUT_PNG, gImp)
ggsave(OUT_SVG, gImp, device = svglite::svglite, bg = "transparent", fix_text_size = FALSE)
qsave(gImp, OUT_QS)

message('INFO: Overall importance plot complete')
