#!/usr/bin/env Rscript
# K search readout (2026-10-14, Trifolium86MainKsearch). The pipeline renders cross-entropy as a PNG only (filed), so the
# values are read back from the .snmfProject here. Writes, under ANALYSIS/Trifolium86MainKsearch/ksearch/:
#   cross_entropy.tsv  K, mean, sd, min, max over the repeats
#   site_by_K.tsv      per K: each site's dominant cluster (best run) and its mean max-Q, to see which sites split off
# Run inside the image from /pipeline: Rscript ANALYSIS/Trifolium86MainKsearch/ksearch_tables.R
suppressPackageStartupMessages({ library(LEA); library(data.table) })
W   <- "Trifolium86MainKsearch_results/_work/maf0.1_miss0.3_smiss0.8_rel0.999/ld0.2_win100_step10"
OUT <- "ANALYSIS/Trifolium86MainKsearch/ksearch"
dir.create(OUT, showWarnings = FALSE)
proj <- load.snmfProject(file.path(W, "Trifolium86_chr1-6.snmfProject"))
meta <- fread("Trifolium86MainKsearch_results/Processing/tables/metadata.tsv",
              colClasses = c(site = "character", sample = "character"))
fam  <- fread(file.path(W, "Trifolium86_chr1-6.fam"), header = FALSE, colClasses = "character")
samples <- sub("^0_", "", fam$V2)
stopifnot(all(samples %in% meta$sample))
site <- meta$site[match(samples, meta$sample)]

ce <- rbindlist(lapply(2:8, function(k) {
  x <- cross.entropy(proj, K = k)
  data.table(K = k, mean = mean(x), sd = sd(x), min = min(x), max = max(x))
}))
fwrite(ce, file.path(OUT, "cross_entropy.tsv"), sep = "\t")
print(ce)

sk <- rbindlist(lapply(2:8, function(k) {
  q <- Q(proj, K = k, run = which.min(cross.entropy(proj, K = k)))
  stopifnot(nrow(q) == length(samples))
  dom <- max.col(q, ties.method = "first")
  data.table(site = site, dom = dom, maxq = apply(q, 1, max))[
    , .(K = k, n = .N, clusters = paste(sort(unique(dom)), collapse = "+"),
        dominant = as.integer(names(which.max(table(dom)))), mean_maxq = round(mean(maxq), 3)), by = site]
}))
setorder(sk, K, dominant, site)
fwrite(sk, file.path(OUT, "site_by_K.tsv"), sep = "\t")
print(dcast(sk, site + n ~ K, value.var = "clusters"))
