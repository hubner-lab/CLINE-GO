#!/usr/bin/env Rscript
# rho_scale.R: CLINE-GO's offsets on one rho scale with GrENE-net's own predictors, and how closely
# each of our runs follows Table S12 garden by garden and plot by plot.
#
#   Rscript manuscript/arabidopsis/rho_scale.R <project>
#
# Reads ANALYSIS/Arabidopsis1135/offset/<project>/{summary,garden_means,replicate_rho}.tsv, written by
# offset_eval.R, and GrENE-net's Table S12. Writes into the same folder:
#   rho_scale.tsv         one row per predictor: our 18 runs (3 methods x 6 sets), our random draws of the
#                         union set (20 per method where they exist), our climate distance, GrENE-net's three
#                         predictors as its text states them, and its replicate prediction, which only
#                         Table S12 gives.
#   garden_agreement.tsv  Spearman between our rho and Table S12's per-plot rho for the matching predictor,
#                         over the 30 gardens (garden means) and over the 325 plots, each with its two-sided
#                         p (cor.test, asymptotic for ties). The garden is the independent unit; the plot-level
#                         p ignores that plots share a garden and is written for completeness only.
#
# Rules:
#   statistic  The paper's: the mean over the 325 plots of the per-plot Spearman, with R's type-7
#              quartiles. Our runs' values are summary.tsv's own cells; the draws' are computed the same
#              way from replicate_rho.tsv, and the computation is checked against summary.tsv on the runs.
#   paper      GrENE-net's values come from its text: the literal table below names the line of
#              ANALYSIS/Arabidopsis1135/literature/grenenet/10.1126_science.adz0777.txt each was read
#              from (user decision 2026-11-09: the text's numbers are the references). Table S12's
#              genomic-offset column (mean 0.340) disagrees with the text's 0.263, so it never enters
#              rho_scale.tsv; its one use is garden_agreement.tsv, whose rows name it.
#   experiment uses_experiment_data is "yes" for a predictor trained on GrENE-net's evolved populations:
#              leave-one-out over the other 29 gardens, or the other plots of the same garden.

suppressPackageStartupMessages(library(data.table))

usage <- "usage: Rscript manuscript/arabidopsis/rho_scale.R <project>"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || startsWith(args[1], "-")) stop(usage)
project <- args[1]

offset_dir <- file.path("ANALYSIS", "Arabidopsis1135", "offset", project)
s12_path <- file.path("data", "arabidopsis1135", "sources", "dl_grene", "Supplemental_tables",
                      "Table_S12_LOO_accuracy.csv")
text_path <- "ANALYSIS/Arabidopsis1135/literature/grenenet/10.1126_science.adz0777.txt"
for (f in c(file.path(offset_dir, c("summary.tsv", "garden_means.tsv", "replicate_rho.tsv")), s12_path)) {
  if (!file.exists(f)) stop("missing input: ", f)
}

method_label <- c(gradient_forest = "Gradient Forest", geometric_offset = "Geometric offset",
                  rda_offset = "RDA offset")
set_label <- c(gea.union = "union set", gea.emmax = "EMMAX set", gea.lfmm = "LFMM set",
               gea.rda = "RDA set", gea.two_of_three = "2-of-3 set", gea.three_of_three = "3-of-3 set")
draw_pattern <- "^random\\.(.+)\\.([0-9]+)$"
our_climate <- "WorldClim 2.1 30s normals (homes and gardens)"
their_climate <- "ERA5-Land BIOCLIM (accessions 1995-2005, gardens 2018-2022)"

quartiles <- function(x) quantile(x, c(0.25, 0.75), names = FALSE, type = 7)
spearman <- function(x, y) cor(x, y, method = "spearman")
spearman_p <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))$p.value

# GrENE-net's predictors, as its main text states them (mean over 325 plots, with the text's IQR).
grene_text <- data.table(
  predictor = c("grene.climate_distance", "grene.genomic_offset", "grene.genomic_offset_gaussian"),
  label = c("Climate distance (BIO1)",
            "Genomic offset (BIO1, trained on the other 29 gardens)",
            "Genomic offset with Gaussian local adaptation parameters (Vs-1, Wmax)"),
  uses_experiment_data = c("no", "yes", "yes"),
  mean = c(0.181, 0.263, 0.415), q1 = c(0.020, 0.160, 0.311), q3 = c(0.332, 0.368, 0.544),
  citation = paste(text_path, c("line 1088", "lines 1086-1087", "lines 1094-1095")))

# Table S12 -------------------------------------------------------------------------------------------

s12 <- fread(s12_path)
s12_spread <- s12[, .(mean = mean(sp_r), q1 = quartiles(sp_r)[1], q3 = quartiles(sp_r)[2], n = .N),
                  by = source]
check_text <- function(s12_source, text_key) {
  s <- s12_spread[source == s12_source]
  t <- grene_text[predictor == text_key]
  same <- all(round(c(s$mean, s$q1, s$q3), 3) == c(t$mean, t$q1, t$q3))
  if (!same) stop(sprintf("Table S12 '%s' (%.4f, %.4f-%.4f) no longer equals the text (%s)",
                          s12_source, s$mean, s$q1, s$q3, t$citation))
}
check_text("climate_distance", "grene.climate_distance")
check_text("stabilizing_era5", "grene.genomic_offset_gaussian")
s12_offset_mean <- s12_spread[source == "genomic offset", mean]
message(sprintf("Table S12 equals the text on climate distance and the Gaussian model; its genomic offset is %.3f against the text's %.3f",
                s12_offset_mean, grene_text[predictor == "grene.genomic_offset", mean]))

replicate_row <- s12_spread[source == "plot_mean"]
grene_rows <- rbind(
  grene_text[, .(source = "GrENE-net text", predictor, label, method = NA_character_, set = NA_character_,
                 uses_experiment_data, mean, q1, q3, n_plots = 325L, n_founders = 231L,
                 climate_source = their_climate, citation)],
  data.table(source = "GrENE-net Table S12", predictor = "grene.replicate_prediction",
             label = "Replicate prediction (mean of the other plots of the garden)",
             method = NA_character_, set = NA_character_, uses_experiment_data = "yes",
             mean = replicate_row$mean, q1 = replicate_row$q1, q3 = replicate_row$q3,
             n_plots = replicate_row$n, n_founders = 231L, climate_source = NA_character_,
             citation = paste(s12_path, "rows source == plot_mean, column sp_r")))

# Our runs, draws and climate distance -----------------------------------------------------------------

summary_runs <- fread(file.path(offset_dir, "summary.tsv"))
per_replicate <- fread(file.path(offset_dir, "replicate_rho.tsv"))
garden_means <- fread(file.path(offset_dir, "garden_means.tsv"))
if (any(grepl(draw_pattern, summary_runs$set))) stop("summary.tsv unexpectedly holds random draws")

recomputed <- per_replicate[!grepl(draw_pattern, set),
                            .(mean = mean(rho_offset), q1 = quartiles(rho_offset)[1],
                              q3 = quartiles(rho_offset)[2], n = .N), by = .(method, set)]
check <- summary_runs[recomputed, on = .(method, set)]
if (nrow(check) != nrow(summary_runs) ||
    check[, max(abs(c(rho_mean - mean, rho_q1 - q1, rho_q3 - q3)))] > 1e-12 ||
    any(check$n_replicates != check$n)) {
  stop("replicate_rho.tsv does not reproduce summary.tsv's rho_mean / rho_q1 / rho_q3")
}

run_rows <- summary_runs[, .(source = "CLINE-GO run",
                             predictor = sprintf("ours.%s.%s", method, set),
                             label = sprintf("%s, %s", method_label[method], set_label[set]),
                             method, set, uses_experiment_data = "no",
                             mean = rho_mean, q1 = rho_q1, q3 = rho_q3,
                             n_plots = n_replicates, n_founders = n_founders_min,
                             climate_source = our_climate,
                             citation = sprintf("%s/summary.tsv: rho_mean, rho_q1, rho_q3 (method %s, set %s)",
                                                offset_dir, method, set))]
if (anyNA(run_rows$label)) stop("a run has a method or set without a label")

draws <- per_replicate[grepl(draw_pattern, set),
                       .(mean = mean(rho_offset), q1 = quartiles(rho_offset)[1],
                         q3 = quartiles(rho_offset)[2], n_plots = .N, n_founders = max(n_founders)),
                       by = .(method, set)]
draw_rows <- draws[, .(source = "CLINE-GO random draw",
                       predictor = sprintf("ours.%s.%s", method, set),
                       label = sprintf("%s, random draw %s of the %s", method_label[method],
                                       sub(draw_pattern, "\\2", set), set_label[sub(draw_pattern, "\\1", set)]),
                       method, set, uses_experiment_data = "no", mean, q1, q3, n_plots, n_founders,
                       climate_source = our_climate,
                       citation = sprintf("%s/replicate_rho.tsv: rho_offset over the plots (method %s, set %s)",
                                          offset_dir, method, set))]

climate_values <- unique(summary_runs[, .(climate_mean, climate_q1, climate_q3, n_founders_min, n_replicates)])
if (nrow(climate_values) != 1L) stop("the runs do not share one climate-distance baseline")
climate_row <- climate_values[, .(source = "CLINE-GO baseline", predictor = "ours.climate_distance",
                                  label = "Climate distance (BIO1)", method = NA_character_,
                                  set = NA_character_, uses_experiment_data = "no",
                                  mean = climate_mean, q1 = climate_q1, q3 = climate_q3,
                                  n_plots = n_replicates, n_founders = n_founders_min,
                                  climate_source = our_climate,
                                  citation = sprintf("%s/summary.tsv: climate_mean, climate_q1, climate_q3 (every run)",
                                                     offset_dir))]

rho_scale <- rbind(grene_rows, climate_row, run_rows, draw_rows)
rho_scale[, project := project]
setcolorder(rho_scale, "project")
fwrite(rho_scale, file.path(offset_dir, "rho_scale.tsv"), sep = "\t")

# Agreement with Table S12, garden by garden and plot by plot ------------------------------------------

climate_per_garden <- garden_means[!grepl(draw_pattern, set), .(n_values = uniqueN(rho_climate)), by = garden]
if (any(climate_per_garden$n_values != 1L)) stop("rho_climate differs between runs")
one_run <- garden_means[method == "gradient_forest" & set == "gea.union"]
one_run_plots <- per_replicate[method == "gradient_forest" & set == "gea.union"]
if (nrow(one_run) != 30L || nrow(one_run_plots) != 325L) stop("expected 30 gardens and 325 plots")

s12_offset_note <- sprintf("Table S12 genomic offset, per plot (its mean %.3f differs from the text's %.3f)",
                           s12_offset_mean, grene_text[predictor == "grene.genomic_offset", mean])
agree_climate <- data.table(
  source = "CLINE-GO baseline", predictor = "ours.climate_distance", method = NA_character_,
  set = NA_character_, against = "Table S12 climate distance, per plot",
  gardens_n = nrow(one_run), gardens_spearman = spearman(one_run$rho_climate, one_run$s12_climate),
  gardens_p = spearman_p(one_run$rho_climate, one_run$s12_climate),
  plots_n = nrow(one_run_plots),
  plots_spearman = spearman(one_run_plots$rho_climate, one_run_plots$s12_climate),
  plots_p = spearman_p(one_run_plots$rho_climate, one_run_plots$s12_climate))

by_garden <- garden_means[, .(gardens_n = .N, gardens_spearman = spearman(rho_offset, s12_offset),
                              gardens_p = spearman_p(rho_offset, s12_offset)), by = .(method, set)]
by_plot <- per_replicate[, .(plots_n = .N, plots_spearman = spearman(rho_offset, s12_offset),
                             plots_p = spearman_p(rho_offset, s12_offset)), by = .(method, set)]
agree_offsets <- by_garden[by_plot, on = .(method, set)]
agree_offsets[, `:=`(source = fifelse(grepl(draw_pattern, set), "CLINE-GO random draw", "CLINE-GO run"),
                     predictor = sprintf("ours.%s.%s", method, set), against = s12_offset_note)]
agreement <- rbind(agree_climate, agree_offsets, use.names = TRUE)
agreement[, project := project]
setcolorder(agreement, c("project", "source", "predictor", "method", "set", "against",
                         "gardens_n", "gardens_spearman", "gardens_p", "plots_n", "plots_spearman", "plots_p"))
fwrite(agreement, file.path(offset_dir, "garden_agreement.tsv"), sep = "\t")

# Report -----------------------------------------------------------------------------------------------

message(sprintf("%s: rho_scale.tsv %d rows (%d runs, %d draws), garden_agreement.tsv %d rows",
                project, nrow(rho_scale), nrow(run_rows), nrow(draw_rows), nrow(agreement)))
message(sprintf("climate distance: ours %.3f, theirs (text) %.3f; agreement with Table S12 %.2f over gardens, %.2f over plots",
                climate_row$mean, grene_text[predictor == "grene.climate_distance", mean],
                agree_climate$gardens_spearman, agree_climate$plots_spearman))
for (m in names(method_label)) {
  u <- agree_offsets[method == m & set == "gea.union"]
  d <- agree_offsets[method == m & grepl(draw_pattern, set)]
  if (nrow(u) == 1L) {
    message(sprintf("%s union: agreement %.2f gardens, %.2f plots%s", method_label[m],
                    u$gardens_spearman, u$plots_spearman,
                    if (nrow(d)) sprintf("; its %d draws %.2f-%.2f gardens, %.2f-%.2f plots", nrow(d),
                                         min(d$gardens_spearman), max(d$gardens_spearman),
                                         min(d$plots_spearman), max(d$plots_spearman)) else ""))
  }
}
