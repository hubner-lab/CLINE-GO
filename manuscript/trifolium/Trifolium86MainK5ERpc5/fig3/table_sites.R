#!/usr/bin/env Rscript
# t_sites (supplement build 2026-10-25, C1): one row per sampling site, sites 1-11 south to north.
# Inputs: the run's input metadata (86 plants, by its `site` column), Processing/tables/metadata.tsv (the 57 plants
# retained after quality control, coordinates as the pipeline wrote them), PreStructure/tables/K5/clusters_K5.tsv
# (sNMF Q at K = 5) and fig3/site_names.tsv (short names 1-11).
# Most frequent cluster: each plant assigned to its largest proportion, as build_fig3.R assigns `cluster`
# (max.col, ties.method = "first"); per site the cluster with the most plants. Ties are reported, never broken silently.
# Writes manuscript_fig3/tables/t_sites.{tsv,md}.
# Run inside the image from /pipeline: Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/table_sites.R
suppressPackageStartupMessages(library(data.table))
AN  <- "/pipeline/ANALYSIS/Trifolium86MainK5ERpc5"
P5  <- "/pipeline/Trifolium86MainK5ERpc5_results"
OUT <- file.path(AN, "manuscript_fig3/tables")

SN   <- fread(file.path(AN, "fig3/site_names.tsv"), colClasses = "character")
inp  <- fread("/pipeline/data/Trifolium_location_86.tsv", colClasses = "character")
meta <- fread(file.path(P5, "Processing/tables/metadata.tsv"), colClasses = "character")
cl   <- fread(file.path(P5, "PreStructure/tables/K5/clusters_K5.tsv"),
              colClasses = list(character = c("sample", "site")))
QC <- paste0("C", 1:5)

stopifnot(nrow(inp) == 86, !anyDuplicated(inp$sample), nrow(SN) == 11,
          nrow(meta) == 57, !anyDuplicated(meta$sample),
          setequal(meta$sample, cl$sample), nrow(cl) == 57,
          all(meta$sample %in% inp$sample),
          identical(cl[order(sample), site], meta[order(sample), site]),
          setequal(unique(inp$site), SN$site_id))
# the plant-name prefix is not the site: 29107-6 belongs to site 29017 (input metadata)
stopifnot(inp[sample == "29107-6", site] == "29017")
message("plants whose name prefix differs from their site: ",
        paste(inp[sub("-.*", "", sample) != site, paste0(sample, " (site ", site, ")")], collapse = ", "))

# ---- coordinates: metadata.tsv form printed, input form agrees within 1e-9 --------------------------------------------
crd <- unique(meta[, .(site, latitude, longitude)])
stopifnot(!anyDuplicated(crd$site), nrow(crd) == 11)
crd_in <- unique(inp[, .(site, lat_in = as.numeric(latitude), lon_in = as.numeric(longitude))])
stopifnot(!anyDuplicated(crd_in$site))
crd <- merge(crd, crd_in, by = "site")
dlat <- abs(as.numeric(crd$latitude) - crd$lat_in); dlon <- abs(as.numeric(crd$longitude) - crd$lon_in)
message(sprintf("coordinates, metadata.tsv vs input: max |dlat| %.3g, max |dlon| %.3g", max(dlat), max(dlon)))
stopifnot(max(dlat) < 1e-9, max(dlon) < 1e-9)
crd <- merge(crd, SN[, .(site = site_id, short, lat_sn = latitude, lon_sn = longitude)], by = "site")
stopifnot(identical(crd$latitude, crd$lat_sn), identical(crd$longitude, crd$lon_sn))

# ---- counts --------------------------------------------------------------------------------------------------------
n_gen <- inp[, .(plants_genotyped = .N), by = site]
n_ret <- meta[, .(plants_retained = .N), by = site]

# ---- clusters ------------------------------------------------------------------------------------------------------
cl[, cluster := QC[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = QC]
qmax_ties <- cl[, sum(apply(as.matrix(.SD), 1, function(q) sum(q == max(q)) > 1)), .SDcols = QC]
message("plants with a tie for their largest proportion: ", qmax_ties)
cnt <- cl[, .N, by = .(site, cluster)]
top <- cnt[, .(most_frequent_cluster = paste(sort(cluster[N == max(N)]), collapse = "/"),
               n_top = max(N), tie = sum(N == max(N)) > 1,
               plants_per_cluster = paste(sprintf("%s %d", cluster[order(cluster)], N[order(cluster)]), collapse = ", ")),
           by = site]
if (any(top$tie)) message("TIE for the most frequent cluster at site(s): ", paste(top[tie == TRUE, site], collapse = ", "))
mq <- cl[, lapply(.SD, mean), by = site, .SDcols = QC]
setnames(mq, QC, paste0("mean_", QC))

t <- Reduce(function(a, b) merge(a, b, by = "site"), list(crd, n_gen, n_ret, top, mq))
t[, short := as.integer(short)]
setorder(t, short)
stopifnot(identical(t$short, 1:11), !is.unsorted(as.numeric(t$latitude)))   # south to north
stopifnot(t$plants_genotyped == c(9, 9, 5, 10, 7, 10, 10, 10, 5, 2, 9),
          t$plants_retained  == c(4, 5, 3, 7, 5, 5, 8, 8, 4, 2, 6),
          sum(t$plants_genotyped) == 86, sum(t$plants_retained) == 57)
stopifnot(t[, all(abs(mean_C1 + mean_C2 + mean_C3 + mean_C4 + mean_C5 - 1) < 1e-4)])
cs <- fread(file.path(AN, "fig3/offset_minou/paper_pack/consensus_sites.tsv"), colClasses = c(site = "character"))
stopifnot(identical(as.integer(cs[match(t$site, site), n]), t$plants_retained))

# ---- s_ancestry_by_site's site order (deviation D5): 4, 5, 9, 11, 10, 6, 3, 2, 1, 8, 7 -------------------------------
# build_fig3.R orders sites by cnt[order(cluster, -N)][!duplicated(site)]: a site's FIRST row there is its
# lowest-numbered cluster that holds any plant, which equals its most frequent cluster only where they coincide.
fig_ord <- cnt[order(cluster, -N)][!duplicated(site)]
fig_ord[, short := as.integer(SN$short[match(site, SN$site_id)])]
FIG_ORDER <- c(4L, 5L, 9L, 11L, 10L, 6L, 3L, 2L, 1L, 8L, 7L)
stopifnot(identical(fig_ord$short, FIG_ORDER))
fig_ord <- merge(fig_ord, t[, .(site, most_frequent_cluster, n_top)], by = "site", sort = FALSE)
fig_ord[, agrees := cluster == most_frequent_cluster]
message("s_ancestry_by_site groups (figure rule vs most frequent cluster):")
print(fig_ord[, .(short, site, figure_group = cluster, plants_in_figure_group = N, most_frequent_cluster, n_top, agrees)])
# The figure's rule (asserted above to reproduce its order exactly) is NOT "most frequent cluster" (user, 2026-10-25:
# deliver the true most frequent cluster, assert the figure's actual rule, report the difference). The two rules
# differ at exactly one site: site 9 (29072), 1 plant C1 / 3 plants C2, grouped with C1 in the figure.
mf <- t[, .(short, cluster = most_frequent_cluster, n_top)][order(cluster, -n_top)]
message("order by most frequent cluster (for comparison only): ", paste(mf$short, collapse = ", "))
message("sites whose figure group differs from their most frequent cluster: ",
        paste(fig_ord[agrees == FALSE, sprintf("site %d (%s): figure %s, most frequent %s", short, site, cluster,
                                               most_frequent_cluster)], collapse = "; "))
stopifnot(!any(t$tie), identical(fig_ord[agrees == FALSE, short], 9L),
          fig_ord[short == 9L, cluster == "C1" & most_frequent_cluster == "C2"])

# ---- write ---------------------------------------------------------------------------------------------------------
out <- t[, .(site = short, site_id = site, latitude, longitude, plants_genotyped, plants_retained,
             most_frequent_cluster, plants_per_cluster, mean_C1, mean_C2, mean_C3, mean_C4, mean_C5)]
fwrite(out, file.path(OUT, "t_sites.tsv"), sep = "\t")

f3 <- function(x) sprintf("%.3f", x)
md <- c("# t_sites", "",
        paste0("One row per sampling site, south to north. Coordinates as written in the run's `Processing/tables/metadata.tsv`. ",
               "Plants genotyped: input metadata (86 plants, by its site column); plants retained: the 57 that passed ",
               "quality control. Most frequent cluster: each retained plant assigned to its largest sNMF proportion at ",
               "K = 5. Mean proportions over the site's retained plants, 3 decimals. ",
               "Values: `t_sites.tsv`; sources: `Processing/tables/metadata.tsv`, `PreStructure/tables/K5/clusters_K5.tsv`, ",
               "`fig3/site_names.tsv`, the input metadata; built by `fig3/table_sites.R`."), "",
        "| Site | Sampling identifier | Latitude | Longitude | Plants genotyped | Plants retained | Most frequent cluster | Mean C1 | Mean C2 | Mean C3 | Mean C4 | Mean C5 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        out[, sprintf("| %d | %s | %s | %s | %d | %d | %s | %s | %s | %s | %s | %s |", site, site_id, latitude, longitude,
                      plants_genotyped, plants_retained, most_frequent_cluster,
                      f3(mean_C1), f3(mean_C2), f3(mean_C3), f3(mean_C4), f3(mean_C5))])
stopifnot(!any(grepl("-[0-9]", md[-(1:3)])))   # no negative numbers expected; never an ASCII minus before a digit
writeLines(md, file.path(OUT, "t_sites.md"))
print(out)
message("done")
