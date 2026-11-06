#!/usr/bin/env Rscript
# pc_climate_r2.R: how much of each leading genomic PC the RDA predictors explain.
#
# When a PC tracks climate, using it as a covariate removes climate signal, and lambda does not show
# it. No pipeline table holds this number, so the choice of n_pcs, K and condition_pcs reads it here.
#
# Reads the LD-pruned PCA (LEA projections and eigenvalues), samples_order.list, the scaled climate
# table, metadata.tsv, clusters_K<k_best>.tsv and the config. Writes
# ANALYSIS/Arabidopsis1135/pregea/pc_climate_r2.tsv, one row per PC:
#   pc pct_variance cum_pct_variance r2_samples r2_adj_samples r2_sites r2_adj_sites r2_main_cluster
#
# r2_samples is lm(PC ~ predictors) over the samples. r2_sites is the same over one row per site,
# with the PC averaged within the site, as varpart counts sites. r2_main_cluster is the PC against
# each sample's main sNMF cluster: it reads the PCs through samples_order.list, so a wrong sample
# order would show as values near zero, and the script stops if PC1 falls below min_cluster_r2.
#
# Run from the repo root, in the image:
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp -e OPENBLAS_NUM_THREADS=1 \
#     -v $PWD:/data -w /data cline-go:6cdedcc Rscript manuscript/arabidopsis/pc_climate_r2.R

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(yaml)
})

project <- "Arabidopsis1135"
n_pcs <- 20L
min_cluster_r2 <- 0.2

results <- here(paste0(project, "_results"))
config <- read_yaml(here(paste0("config_", project, ".yaml")))
predictors <- strsplit(config$PreGEA$predictors, ",", fixed = TRUE)[[1]]
out_file <- here("ANALYSIS", project, "pregea", "pc_climate_r2.tsv")

#-------------------------------------------------------------------------- inputs
pca_dir <- Sys.glob(file.path(results, "_work", "*", "ld*", "*.pca"))
stopifnot(length(pca_dir) == 1L)
stem <- sub("\\.pca$", "", basename(pca_dir))
scores <- as.matrix(fread(file.path(pca_dir, paste0(stem, ".projections")), sep = " ",
                          header = FALSE, select = seq_len(n_pcs)))
eigenvalues <- scan(file.path(pca_dir, paste0(stem, ".eigenvalues")), quiet = TRUE)

samples <- readLines(file.path(results, "_intermediate", "samples", "samples_order.list"))
stopifnot(nrow(scores) == length(samples), !anyDuplicated(samples))

climate <- fread(file.path(results, "climate", "tables", "present", "climate_present_site_scaled.tsv"),
                 colClasses = c(sample = "character"))
metadata <- fread(file.path(results, "Processing", "tables", "metadata.tsv"),
                  colClasses = c(site = "character", sample = "character"))
clusters <- fread(file.path(results, "PreStructure", "tables", paste0("K", config$sNMF$k_best),
                            paste0("clusters_K", config$sNMF$k_best, ".tsv")),
                  colClasses = c(sample = "character", site = "character"))
q_cols <- grep("^C[0-9]+$", names(clusters), value = TRUE)
clusters[, main := q_cols[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = q_cols]

stopifnot(setequal(samples, climate$sample), setequal(samples, metadata$sample),
          setequal(samples, clusters$sample), all(predictors %in% names(climate)))
X <- as.matrix(climate[match(samples, sample), ..predictors])
site <- metadata[match(samples, sample), site]
main <- clusters[match(samples, sample), main]

#-------------------------------------------------------------------------- R2 per PC
fit_r2 <- function(y, x) {
  s <- summary(lm(y ~ x))
  c(s$r.squared, s$adj.r.squared)
}

site_ids <- unique(site)
X_sites <- X[match(site_ids, site), , drop = FALSE]
stopifnot(all.equal(X, X_sites[match(site, site_ids), , drop = FALSE], check.attributes = FALSE))

rows <- lapply(seq_len(n_pcs), function(i) {
  on_samples <- fit_r2(scores[, i], X)
  site_mean <- tapply(scores[, i], factor(site, levels = site_ids), mean)
  on_sites <- fit_r2(as.numeric(site_mean), X_sites)
  data.table(pc = i, r2_samples = on_samples[1], r2_adj_samples = on_samples[2],
             r2_sites = on_sites[1], r2_adj_sites = on_sites[2],
             r2_main_cluster = summary(lm(scores[, i] ~ main))$r.squared)
})
result <- rbindlist(rows)
result[, `:=`(pct_variance = 100 * eigenvalues[pc] / sum(eigenvalues),
              cum_pct_variance = 100 * cumsum(eigenvalues)[pc] / sum(eigenvalues))]
setcolorder(result, c("pc", "pct_variance", "cum_pct_variance"))

if (result[pc == 1L, r2_main_cluster] < min_cluster_r2) {
  stop("PC1 explains under ", min_cluster_r2, " of the main cluster: check the sample order")
}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
fwrite(result[, lapply(.SD, signif, 4)], out_file, sep = "\t")
message("Predictors: ", paste(predictors, collapse = ", "))
message("Wrote ", out_file, ": ", nrow(result), " PCs, ", length(samples), " samples, ",
        length(site_ids), " sites")
