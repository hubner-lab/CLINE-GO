#!/usr/bin/env Rscript
# relatedness_thresholds.R: related pairs and removals at each Filter.relatedness threshold.
#
#   Rscript manuscript/arabidopsis/relatedness_thresholds.R <project> [--compare <other project>]
#
# The pipeline writes relatedness_removed.tsv for the configured threshold only, and its histogram
# bins IBS at 0.02, too coarse to read a clone tail (finding of 2026-11-06). This replays rule
# relatedness_filter (workflow/rules/processing.smk) on the project's relatedness_pairs.tsv at every
# threshold from 0.900 to 0.999 in steps of 0.001. The replay runs the rule's own commands (tail, awk,
# sort -k3,3gr and its greedy awk program: pairs in decreasing IBS, each pair above the threshold
# with both members still kept loses the member with more missing calls), in the image's locale, so
# ties at equal IBS break as they do in the pipeline. When the project's config removes relatives and
# its tree holds relatedness_removed.tsv, the replay at that threshold must reproduce the file byte
# for byte, or the script stops.
#
# Reads <project>_results/Processing/tables/relatedness_pairs.tsv, _intermediate/samples/
# {samples_missing_stats.tsv, samples_filtered.list} (the rule's inputs when no het filter is set; a
# config with Filter.het_outlier_sd stops), <project>_logs/processing/compute_relatedness.log (the
# SNP counts) and data/arabidopsis1135/founder_sites.tsv. Writes under ANALYSIS/<project>/relatedness/:
#   thresholds.tsv  threshold, n_pairs_above, n_removed, n_kept, founders_removed, founders_kept, and
#                   with --compare: n_removed_other (the other tree's relatedness_removed.tsv) and
#                   shared_removed (samples both remove)
#   density.tsv     IBS bin (width 0.001, lower edge, from 0.900 up), pairs in it
#   summary.tsv     samples, pairs, SNPs passing the filters and after pruning, IBS quantiles
#   ibs_tail.png    log10(pairs + 1) per 0.001 IBS from 0.900 up, and samples removed per threshold
#
# Run from the repo root, in the image:
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp -v $PWD:/data -w /data \
#     cline-go:6cdedcc Rscript manuscript/arabidopsis/relatedness_thresholds.R <project> [--compare P]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(here)
  library(yaml)
})
source(here("scripts", "R", "utils", "theme_clinego.R"))

usage <- "usage: relatedness_thresholds.R <project> [--compare <other project>]"
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 3L) || (length(args) == 3L && args[2] != "--compare")) stop(usage)
project <- args[1]
other_removed <- NULL
if (length(args) == 3L) {
  other_file <- here(paste0(args[3], "_results"), "Processing", "tables", "relatedness_removed.tsv")
  if (!file.exists(other_file)) stop("no ", other_file)
  other_removed <- fread(other_file, colClasses = c(IID = "character"))$IID
}
config <- read_yaml(here(paste0("config_", project, ".yaml")))
if (!is.null(config$Filter$het_outlier_sd))
  stop("Filter.het_outlier_sd is set: relatedness_filter would read the het-filtered list, not replayed here")

results <- here(paste0(project, "_results"))
pairs_file <- file.path(results, "Processing", "tables", "relatedness_pairs.tsv")
imiss_file <- file.path(results, "_intermediate", "samples", "samples_missing_stats.tsv")
samples_file <- file.path(results, "_intermediate", "samples", "samples_filtered.list")
removed_file <- file.path(results, "Processing", "tables", "relatedness_removed.tsv")
log_file <- here(paste0(project, "_logs"), "processing", "compute_relatedness.log")
out_dir <- here("ANALYSIS", project, "relatedness")
stopifnot(file.exists(pairs_file), file.exists(imiss_file), file.exists(samples_file))
thresholds <- round(seq(0.900, 0.999, by = 0.001), 3)

#-------------------------------------------------------------------------- the rule, replayed
# The shell of rule relatedness_filter with Snakemake's doubled braces undone and its paths as $1-$5:
# $1 pairs, $2 imiss, $3 samples, $4 threshold, $5 output folder (removed.tsv, kept.list).
replay_sh <- r"---(
tail -n +2 "$2" > "$5/imiss.tsv"
tail -n +2 "$1" | awk -v thresh="$4" '$3 > thresh' | sort -k3,3gr > "$5/sorted_pairs.tsv"
awk -v thresh="$4" -v removed_out="$5/removed.tsv" '
BEGIN { OFS="\t"; file_num = 0; print "IID\tIBS\tF_MISS" > removed_out }
FNR==1 { file_num++ }
file_num==1 { fmiss[$2]=$6; next }
file_num==2 { keep[$2]=1; fid[$2]=$1; order[++n]=$2; next }
file_num==3 {
    iid1=$1; iid2=$2; ibs=$3
    if (keep[iid1]==1 && keep[iid2]==1) {
        if (fmiss[iid1]+0 >= fmiss[iid2]+0) { drop=iid1 } else { drop=iid2 }
        keep[drop]=0
        print drop, ibs, fmiss[drop] >> removed_out
    }
}
END {
    for (i=1;i<=n;i++) { iid=order[i]; if (keep[iid]==1) print fid[iid], iid }
}
' "$5/imiss.tsv" "$3" "$5/sorted_pairs.tsv" > "$5/kept.list"
)---"
replay <- function(threshold) {
  dir <- tempfile("relatedness_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  status <- system2("bash", shQuote(c("-c", replay_sh, "replay", pairs_file, imiss_file,
                                       samples_file, format(threshold, nsmall = 3), dir)))
  if (status != 0L) stop("the replay at ", threshold, " exited ", status)
  list(removed = readLines(file.path(dir, "removed.tsv")),
       kept = fread(file.path(dir, "kept.list"), header = FALSE, colClasses = "character")$V2)
}

#-------------------------------------------------------------------------- check against the tree
configured <- config$Filter$relatedness
if (!is.null(configured) && identical(tolower(config$Filter$relatedness_action), "remove") &&
    file.exists(removed_file)) {
  if (!identical(replay(configured)$removed, readLines(removed_file)))
    stop("the replay at ", configured, " does not reproduce ", removed_file)
  message("The replay at ", configured, " reproduces ", basename(removed_file), " byte for byte")
}

#-------------------------------------------------------------------------- tables
pairs <- fread(pairs_file, colClasses = c(IID1 = "character", IID2 = "character"))
samples <- fread(samples_file, header = FALSE, colClasses = "character")$V2
founders <- fread(here("data", "arabidopsis1135", "founder_sites.tsv"),
                  colClasses = c(sample = "character"))[sample %in% samples, sample]

by_threshold <- rbindlist(lapply(thresholds, function(t) {
  r <- replay(t)
  removed <- fread(text = r$removed, colClasses = c(IID = "character"))$IID
  if (length(removed) + length(r$kept) != length(samples)) stop("samples do not add up at ", t)
  row <- data.table(threshold = t, n_pairs_above = sum(pairs$IBS > t), n_removed = length(removed),
                    n_kept = length(r$kept), founders_removed = sum(founders %in% removed),
                    founders_kept = sum(founders %in% r$kept))
  if (!is.null(other_removed))
    row[, `:=`(n_removed_other = length(other_removed),
               shared_removed = sum(removed %in% other_removed))]
  row
}))

# Bins on integer thousandths, so the join never compares doubles.
density <- pairs[IBS >= 0.9, .(pairs = .N), keyby = .(milli = as.integer(floor(IBS * 1000 + 1e-9)))]
density <- density[CJ(milli = 900:1000), on = "milli"][is.na(pairs), pairs := 0L]
density[, bin := milli / 1000]
density <- density[, .(bin, pairs)]

log_lines <- if (file.exists(log_file)) readLines(log_file) else character(0)
snp_count <- function(pattern) {
  hit <- regmatches(log_lines, regexpr(pattern, log_lines))
  if (length(hit) == 0L) NA_integer_ else as.integer(sub(" .*", "", sub("^--extract: ", "", hit[1])))
}
quantiles <- quantile(pairs$IBS, c(0.5, 0.9, 0.99, 0.999))
summary <- data.table(
  samples = length(samples), pairs = nrow(pairs),
  snps_passing_filters = snp_count("^[0-9]+ variants and [0-9]+ (samples|people) pass filters"),
  snps_pruned = snp_count("^--extract: [0-9]+ variants remaining"),
  ibs_median = quantiles[[1]], ibs_q90 = quantiles[[2]], ibs_q99 = quantiles[[3]],
  ibs_q999 = quantiles[[4]], pairs_at_1 = sum(pairs$IBS >= 1))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(by_threshold, file.path(out_dir, "thresholds.tsv"), sep = "\t")
fwrite(density, file.path(out_dir, "density.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "summary.tsv"), sep = "\t")

#-------------------------------------------------------------------------- figure
panels <- c("log10(pairs + 1) per 0.001 IBS", "Samples removed at threshold")
plot_data <- rbind(
  density[, .(ibs = bin + 0.0005, value = log10(pairs + 1), panel = panels[1])],
  by_threshold[, .(ibs = threshold, value = n_removed, panel = panels[2])])
plot_data[, panel := factor(panel, levels = panels)]
plot <- ggplot(plot_data, aes(ibs, value)) +
  geom_col(data = plot_data[panel == panels[1]], width = 0.001, fill = CLINEGO_NEUTRAL) +
  geom_step(data = plot_data[panel == panels[2]], colour = CLINEGO_REMOVED, direction = "mid") +
  geom_point(data = plot_data[panel == panels[2]], colour = CLINEGO_REMOVED, size = 0.8) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  labs(title = paste0(project, ": relatedness tail"), x = "IBS (proportion of shared alleles)",
       y = NULL) +
  theme_clinego()
ggsave(file.path(out_dir, "ibs_tail.png"), plot, width = 8, height = 6, dpi = 150)

print(summary)
print(by_threshold[threshold %in% c(0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.995, 0.996, 0.997, 0.998, 0.999)])
message("Wrote ", out_dir)
