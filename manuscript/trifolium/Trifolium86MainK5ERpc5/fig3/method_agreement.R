#!/usr/bin/env Rscript
# method_agreement (supplement build 2026-10-25, C4): Kendall's tau-b and its exact two-sided p between the site offsets
# of each pair of offset methods (Gradient Forest x geometric, Gradient Forest x RDA, geometric x RDA), for the
# associated SNPs (set gea.ER_overlap) of fig3/offset_minou/site_maps/site_offsets.tsv. These p-values replace an
# ad hoc recompute. Exact p: R cor.test(method = "kendall", exact = TRUE), which needs no ties (asserted); with no
# ties tau-b = tau-a. The two-sided p is 2 x the smaller tail, capped at 1.
# Writes manuscript_fig3/tables/method_agreement.tsv. Reads offset_minou, writes nothing there (outside the B9 gate).
# Run inside the image from /pipeline: Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/method_agreement.R
suppressPackageStartupMessages(library(data.table))
AN  <- "/pipeline/ANALYSIS/Trifolium86MainK5ERpc5"
SM  <- file.path(AN, "fig3/offset_minou/site_maps")
OUT <- file.path(AN, "manuscript_fig3/tables")

so <- fread(file.path(SM, "site_offsets.tsv"), colClasses = c(site = "character"))[set == "gea.ER_overlap"]
w  <- dcast(so[method %in% c("gf", "geometric", "rda")], site ~ method, value.var = "offset")
stopifnot(nrow(w) == 11, !anyNA(w))
for (m in c("gf", "geometric", "rda")) stopifnot(!anyDuplicated(w[[m]]))        # exact p needs no ties

PAIRS <- list(c("gf", "geometric"), c("gf", "rda"), c("geometric", "rda"))
LAB   <- c(gf = "Gradient Forest", geometric = "Geometric offset", rda = "RDA-based offset")
res <- rbindlist(lapply(PAIRS, function(p) {
    ct <- cor.test(w[[p[1]]], w[[p[2]]], method = "kendall", exact = TRUE)
    stopifnot(names(ct$statistic) == "T")          # T = the exact statistic; "z" would mean the normal approximation
    data.table(method_1 = LAB[[p[1]]], method_2 = LAB[[p[2]]], pair = paste0(p[1], "~", p[2]),
               n_sites = nrow(w), kendall_tau_b = unname(ct$estimate),
               concordant_minus_discordant = round(unname(ct$estimate) * choose(nrow(w), 2)),
               p_exact_two_sided = ct$p.value, test = ct$method)
}))
print(res)
stopifnot(all(res$concordant_minus_discordant == c(9, -1, 37)),
          all(abs(res$kendall_tau_b * 55 - res$concordant_minus_discordant) < 1e-9))

# ---- asserts: tau against kendall_check.tsv (method_agreement, kendall); p against the text -------------------------
kc <- fread(file.path(SM, "kendall_check.tsv"))[test == "method_agreement" & stat == "kendall", .(pair = method, kc = adaptive)]
chk <- merge(res, kc, by = "pair")
stopifnot(nrow(chk) == 3, chk[, all(abs(kendall_tau_b - kc) < 5e-4)])
stopifnot(abs(res[pair == "gf~geometric", p_exact_two_sided] - 0.54) < 0.005,
          sprintf("%.1f", res[pair == "gf~rda", p_exact_two_sided]) == "1.0",
          abs(res[pair == "geometric~rda", p_exact_two_sided] - 0.0031) < 5e-5)
fwrite(res[, .(method_1, method_2, n_sites, kendall_tau_b, p_exact_two_sided, test)],
       file.path(OUT, "method_agreement.tsv"), sep = "\t")
message("done")
