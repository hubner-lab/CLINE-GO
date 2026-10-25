#!/usr/bin/env Rscript
# predictor_screen_check (supplement build 2026-10-25, C5): the RDA predictor screen at site level for this run.
# The selection itself is recorded: ANALYSIS/Trifolium_rerun_v3/compare/select_predictors_cut.R, run on 2026-10-12 with
# `0.8 bio_17,bio_18` into ANALYSIS/Trifolium_rerun_v3/predictor_selection_r08/ (it picked bio_1,bio_3,bio_5), and
# re-run unchanged on this run's climate_present_site.tsv on 2026-10-25 with byte-identical outputs. This check only
# prints the numbers of the chosen three on this run's data: every pairwise |r| and each VIF at site level (11 sites,
# one row per site, as the selection script collapses them), and asserts the rule |r| < 0.8, VIF < 10.
# Writes manuscript_fig3/tables/predictor_screen.tsv.
# Run inside the image from /pipeline: Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/predictor_screen_check.R
suppressPackageStartupMessages(library(data.table))
AN  <- "/pipeline/ANALYSIS/Trifolium86MainK5ERpc5"
P5  <- "/pipeline/Trifolium86MainK5ERpc5_results"
OUT <- file.path(AN, "manuscript_fig3/tables")
R08 <- "/pipeline/ANALYSIS/Trifolium_rerun_v3/predictor_selection_r08"
SEL <- c("bio_1", "bio_3", "bio_5")

cp   <- fread(file.path(P5, "climate/tables/present/climate_present_site.tsv"), colClasses = c(sample = "character"))
meta <- fread(file.path(P5, "Processing/tables/metadata.tsv"), colClasses = "character")
stopifnot(nrow(cp) == 57, setequal(cp$sample, meta$sample))
cp[, site := meta$site[match(sample, meta$sample)]]
bio <- grep("^bio_", names(cp), value = TRUE)
stopifnot(cp[, lapply(.SD, uniqueN), by = site, .SDcols = bio][, all(as.matrix(.SD) == 1), .SDcols = bio])
env <- unique(cp[, c("site", bio), with = FALSE])
stopifnot(nrow(env) == 11, !anyDuplicated(env$site))
stopifnot(identical(readLines(file.path(R08, "predictor_set_selected.txt")), paste(SEL, collapse = ",")))

E <- as.matrix(env[, ..SEL])
R <- cor(E)
vif <- sapply(seq_along(SEL), function(j) 1 / (1 - summary(lm(E[, j] ~ E[, -j]))$r.squared))
names(vif) <- SEL
pairs <- data.table(t(combn(SEL, 2)))
setnames(pairs, c("predictor_1", "predictor_2"))
pairs[, r := mapply(function(a, b) R[a, b], predictor_1, predictor_2)]
pairs[, abs_r := abs(r)]

# the selection script's own site-level correlation table (4 decimals) agrees
r08 <- fread(file.path(R08, "predictor_correlation_site.tsv"))
for (i in seq_len(nrow(pairs))) stopifnot(abs(r08[predictor == pairs$predictor_1[i]][[pairs$predictor_2[i]]] - pairs$r[i]) < 5e-5)
cand <- fread(file.path(R08, "predictor_set_candidates.tsv"))[selected == TRUE]
stopifnot(nrow(cand) == 1, abs(cand$max_vif - max(vif)) < 1e-9)

stopifnot(all(pairs$abs_r < 0.8), all(vif < 10))
message(sprintf("max pairwise |r| = %.4f (%s x %s); max VIF = %.4f (%s)",
                max(pairs$abs_r), pairs[which.max(abs_r), predictor_1], pairs[which.max(abs_r), predictor_2],
                max(vif), names(vif)[which.max(vif)]))
print(pairs); print(vif)
out <- rbind(pairs[, .(quantity = "pairwise r", predictor_1, predictor_2, value = r)],
             data.table(quantity = "VIF", predictor_1 = SEL, predictor_2 = "", value = unname(vif)))
out[, n_sites := nrow(env)]
fwrite(out, file.path(OUT, "predictor_screen.tsv"), sep = "\t")
message("done")
