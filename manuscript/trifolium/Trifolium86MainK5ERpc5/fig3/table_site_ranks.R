#!/usr/bin/env Rscript
# t_site_ranks (supplement build 2026-10-25, C3): the sites in consensus-score order, as the lab table
# manuscript_fig3/tables/offset_Table1_consensus_ranks.tsv orders them, with each method's offset rank, the consensus
# score and rank, a display category, and BIO1/BIO3/BIO5 present, future and signed change (future - present).
# Inputs: fig3/offset_minou/paper_pack/{consensus_sites,site_climate_shift}.tsv (the minou pack of this run; the V1
# pack under Trifolium86MainK5/ belongs to the 3-PC project) and fig3/site_names.tsv.
# Row order: a stable sort on consensus_score, which keeps consensus_sites.tsv's row order within a tied score.
# Writes manuscript_fig3/tables/t_site_ranks.{tsv,md}. The raw category values stay in sources/consensus_sites.tsv.
# Run inside the image from /pipeline: Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/table_site_ranks.R
suppressPackageStartupMessages(library(data.table))
AN  <- "/pipeline/ANALYSIS/Trifolium86MainK5ERpc5"
PP  <- file.path(AN, "fig3/offset_minou/paper_pack")
OUT <- file.path(AN, "manuscript_fig3/tables")
P5  <- "/pipeline/Trifolium86MainK5ERpc5_results"
MINUS <- "−"

SN <- fread(file.path(AN, "fig3/site_names.tsv"), colClasses = "character")
cs <- fread(file.path(PP, "consensus_sites.tsv"), colClasses = c(site = "character"))
sh <- fread(file.path(PP, "site_climate_shift.tsv"), colClasses = c(site = "character"))
lab <- fread(file.path(OUT, "offset_Table1_consensus_ranks.tsv"), colClasses = c(site = "character", site_id = "character"))
stopifnot(nrow(cs) == 11, nrow(sh) == 11, setequal(cs$site, SN$site_id), setequal(sh$site, cs$site))

# ---- ranks, score, category ------------------------------------------------------------------------------------------
stopifnot(identical(cs$rank_gf, cs$gf_rank))
for (r in c("rank_gf", "rank_geometric", "rank_rda")) stopifnot(setequal(cs[[r]], 1:11))   # no tied offset ranks
stopifnot(all(abs(cs[, (rank_gf + rank_geometric + rank_rda) / 3] - cs$consensus_score) < 1e-12),
          identical(as.integer(cs$consensus_rank), as.integer(frank(cs$consensus_score, ties.method = "min"))))
# category rule of build_paper_pack.R: upper half (rank < 6) under all three = higher, lower half (> 6) = lower
rule <- cs[, fifelse(rank_gf < 6 & rank_geometric < 6 & rank_rda < 6, "higher risk",
                     fifelse(rank_gf > 6 & rank_geometric > 6 & rank_rda > 6, "lower risk", "methods disagree"))]
stopifnot(identical(rule, cs$category))
CAT_LAB <- c("higher risk" = "higher offset in all three methods", "lower risk" = "lower offset in all three methods",
             "methods disagree" = "methods disagree")
stopifnot(all(cs$category %in% names(CAT_LAB)))
cs[, short := as.integer(SN$short[match(site, SN$site_id)])]
stopifnot(setequal(cs[category == "higher risk", short], c(1L, 3L)),
          setequal(cs[category == "lower risk", short], c(8L, 9L)),
          sum(cs$category == "methods disagree") == 7)
ties <- cs[, .(sites = paste(sort(short), collapse = " and "), n = .N), by = consensus_rank][n > 1]
message("consensus-rank ties: ", paste(sprintf("rank %d: sites %s", ties$consensus_rank, ties$sites), collapse = "; "))
stopifnot(nrow(ties) == 2, all(ties$n == 2))

# ---- order: stable sort on the score, equal to the lab table's order ---------------------------------------------------
ord <- cs[order(consensus_score)]                 # data.table order() is stable
stopifnot(identical(ord$site, cs$site), identical(ord$site, lab$site_id),
          identical(as.integer(lab$site), ord$short),
          identical(lab$rank_gf, ord$rank_gf), identical(lab$rank_geometric, ord$rank_geometric),
          identical(lab$rank_rda, ord$rank_rda), identical(lab$consensus_rank, ord$consensus_rank),
          identical(lab$category, ord$category))

# ---- climate: present, future, signed change; site means over the retained plants -----------------------------------
B <- c("bio_1", "bio_3", "bio_5")
for (b in B) {
    d <- sh[[paste0(b, "_future")]] - sh[[paste0(b, "_present")]]
    stopifnot(max(abs(d - sh[[paste0("d_", b)]])) < 1e-9)
}
stopifnot(all(sh$d_bio_3 < 0), all(sh$d_bio_1 > 0), all(sh$d_bio_5 > 0))
m <- merge(cs[, c("site", paste0(B, "_present"), paste0("d_", B)), with = FALSE], sh, by = "site", suffixes = c(".cs", ""))
for (b in B) stopifnot(identical(m[[paste0(b, "_present.cs")]], m[[paste0(b, "_present")]]),
                       identical(m[[paste0("d_", b, ".cs")]], m[[paste0("d_", b)]]))
# the per-plant present climate of the 57 retained plants: one value per site, equal to the site mean
cp <- fread(file.path(P5, "climate/tables/present/climate_present_site.tsv"), colClasses = c(sample = "character"))
meta <- fread(file.path(P5, "Processing/tables/metadata.tsv"), colClasses = "character")
stopifnot(nrow(cp) == 57, setequal(cp$sample, meta$sample))
cp[, site := meta$site[match(sample, meta$sample)]]
pm <- cp[, c(lapply(.SD, mean), lapply(.SD, function(v) diff(range(v)))), by = site, .SDcols = B]
setnames(pm, c("site", paste0(B, "_mean"), paste0(B, "_range")))
stopifnot(all(as.matrix(pm[, paste0(B, "_range"), with = FALSE]) == 0))
pm <- merge(pm, sh, by = "site")
for (b in B) stopifnot(max(abs(pm[[paste0(b, "_mean")]] - pm[[paste0(b, "_present")]])) < 1e-9)
stopifnot(identical(as.integer(cs$n), as.integer(meta[, .N, by = site][match(cs$site, site), N])))

# ---- write ---------------------------------------------------------------------------------------------------------
t <- merge(ord[, .(site, short, rank_gf, rank_geometric, rank_rda, consensus_score, consensus_rank, category)],
           sh, by = "site", sort = FALSE)
stopifnot(identical(t$site, ord$site))
out <- t[, .(site = short, rank_gradient_forest = rank_gf, rank_geometric_offset = rank_geometric,
             rank_rda_based_offset = rank_rda, consensus_score, consensus_rank,
             category = unname(CAT_LAB[category]),
             BIO1_present_C = bio_1_present, BIO1_future_C = bio_1_future, BIO1_change_C = d_bio_1,
             BIO3_present_pct = bio_3_present, BIO3_future_pct = bio_3_future, BIO3_change_pct = d_bio_3,
             BIO5_present_C = bio_5_present, BIO5_future_C = bio_5_future, BIO5_change_C = d_bio_5)]
fwrite(out, file.path(OUT, "t_site_ranks.tsv"), sep = "\t")

num <- function(x, d) sub("^-", MINUS, sprintf(paste0("%.", d, "f"), x))   # formatted once, true minus sign
md <- c("# t_site_ranks", "",
        paste0("Sites in consensus-score order (rows with an equal score keep the order of `consensus_sites.tsv`). ",
               "Offset rank under each method, 1 = highest offset of the 11 sites. Consensus score = mean of the three ",
               "ranks, 2 decimals; consensus rank = rank of that score, ties sharing the lower rank. Climate: site means ",
               "over the retained plants, present and future (2061–2080, SSP585), change = future − present; ",
               "BIO1 and BIO5 in °C, BIO3 in %, 2 decimals. Values: `t_site_ranks.tsv`; sources: ",
               "`paper_pack/consensus_sites.tsv`, `paper_pack/site_climate_shift.tsv`; built by `fig3/table_site_ranks.R`."), "",
        paste0("| Site | Gradient Forest rank | Geometric offset rank | RDA-based offset rank | Consensus score | Consensus rank | Category | ",
               "BIO1 present (°C) | BIO1 future (°C) | BIO1 change (°C) | BIO3 present (%) | BIO3 future (%) | BIO3 change (%) | ",
               "BIO5 present (°C) | BIO5 future (°C) | BIO5 change (°C) |"),
        paste(c("|", rep(" --- |", 16)), collapse = ""),
        out[, sprintf("| %d | %d | %d | %d | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                      site, rank_gradient_forest, rank_geometric_offset, rank_rda_based_offset, num(consensus_score, 2),
                      consensus_rank, category, num(BIO1_present_C, 2), num(BIO1_future_C, 2), num(BIO1_change_C, 2),
                      num(BIO3_present_pct, 2), num(BIO3_future_pct, 2), num(BIO3_change_pct, 2),
                      num(BIO5_present_C, 2), num(BIO5_future_C, 2), num(BIO5_change_C, 2))])
stopifnot(!any(grepl("-[0-9]", md)))
writeLines(md, file.path(OUT, "t_site_ranks.md"))
print(out)
message("done")
