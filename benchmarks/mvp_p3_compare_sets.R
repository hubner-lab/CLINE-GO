#!/usr/bin/env Rscript
# mvp_p3_compare_sets.R -- LD-window plan, Phase 3: compare the framework's declared sets
# (scripts/build_snp_set.R on the OLD p-value trio, staged by work/ld_window_lab/p3_one.sh) with
# the sets benchmarks/mvp_build_snp_sets.R wrote for offset13 into
# MVP<seed>_results/_intermediate/snp_sets/<name>/selected_snps.tsv.
#
# Writes into OUTDIR:
#   set_equality.tsv  seed, set, n_old, n_new, n_offset13 (summary cross-check), identical,
#                     only_old / only_new (SNPIDs), tie_explained (every differing SNP has p exactly
#                     equal to a cell's top-N cutoff)
#   ld_distances.tsv  seed, resolved_bp (build log of best_ld), out600 decay_critical_bp, equal
#   set_sizes.tsv     seed, distance_bp, best / intersect3 at 5 kb and at the critical distance
#
# Usage: Rscript benchmarks/mvp_p3_compare_sets.R <seeds_csv> <p3_dir> <out600_dir> <outdir>
#   (repo root as cwd; old trees are only read)

suppressPackageStartupMessages(library(data.table))
args   <- commandArgs(trailingOnly = TRUE)
seeds  <- strsplit(args[1], ",")[[1]]
P3     <- args[2]; OUT600 <- args[3]; OUTDIR <- args[4]
stopifnot(length(args) == 4, dir.exists(P3), dir.exists(OUT600))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OLD_SETS <- c("solo_lfmm", "solo_emmax", "solo_rda", "union", "best", "intersect3")
summ13 <- fread("benchmarks/mvp_eval/offset13/snp_sets_summary.tsv", colClasses = c(seed = "character"))

read_ids <- function(f) {
    if (!file.exists(f)) return(character(0))       # the old builder writes no file for an empty set
    fread(f, colClasses = c(SNPID = "character", chr = "character"))$SNPID
}
# Cell cutoffs as mvp_build_snp_sets.R defines them: top N = ceiling(share x rows) per column.
cutoffs <- function(trio_dir) {
    out <- list()
    for (m in c("LFMM", "EMMAX", "RDA")) {
        d <- fread(file.path(trio_dir, m, paste0(m, "_pvalues_K3.tsv")),
                   colClasses = c(SNPID = "character", chr = "character"))
        share <- if (m == "RDA") 0.005 else 0.0025
        N <- ceiling(share * nrow(d) - 1e-9)
        for (col in setdiff(names(d), c("SNPID", "chr", "pos"))) {
            cut <- sort(d[[col]])[N]
            out[[paste(col, m, sep = "::")]] <- d[get(col) == cut, SNPID]
        }
    }
    unique(unlist(out))
}

eq <- list(); ld <- list(); sz <- list()
for (s in seeds) {
    pdir <- file.path(P3, paste0("MVP", s))
    ties <- cutoffs(file.path(pdir, "trio"))
    for (n in OLD_SETS) {
        old <- read_ids(sprintf("MVP%s_results/_intermediate/snp_sets/%s/selected_snps.tsv", s, n))
        new <- read_ids(file.path(pdir, "sets", n, "selected_snps.tsv"))
        oo <- setdiff(old, new); on <- setdiff(new, old)
        n13 <- summ13[seed == s & set == n, n_snps]
        eq[[length(eq) + 1]] <- data.table(
            seed = s, set = n, n_old = length(old), n_new = length(new),
            n_offset13 = if (length(n13)) n13[1] else NA_integer_,
            identical = setequal(old, new),
            only_old = paste(oo, collapse = ","), only_new = paste(on, collapse = ","),
            tie_explained = if (setequal(old, new)) NA else all(c(oo, on) %in% ties))
    }
    lg <- readLines(file.path(pdir, "sets", "best_ld", "build.log"))
    bp <- as.integer(sub(".*\\(critical\\): ([0-9]+) bp.*", "\\1",
                         grep("AUTO genome-wide distance \\(critical\\)", lg, value = TRUE)[1]))
    t6 <- fread(file.path(OUT600, paste0("MVP", s), "table.tsv"))[group == "All" & scope == "genome_wide"]
    ld[[length(ld) + 1]] <- data.table(seed = s, resolved_bp = bp,
                                       out600_decay_critical_bp = t6$decay_critical_bp,
                                       equal = isTRUE(bp == round(t6$decay_critical_bp)))
    nn <- function(n) length(read_ids(file.path(pdir, "sets", n, "selected_snps.tsv")))
    sz[[length(sz) + 1]] <- data.table(seed = s, distance_bp = bp,
        best_5kb = nn("best"), best_ld = nn("best_ld"),
        intersect3_5kb = nn("intersect3"), intersect3_ld = nn("intersect3_ld"))
}
eq <- rbindlist(eq); ld <- rbindlist(ld); sz <- rbindlist(sz)
fwrite(eq, file.path(OUTDIR, "set_equality.tsv"), sep = "\t")
fwrite(ld, file.path(OUTDIR, "ld_distances.tsv"), sep = "\t")
fwrite(sz, file.path(OUTDIR, "set_sizes.tsv"), sep = "\t")
message(sprintf("sets identical: %d / %d; offset13 summary agrees with old files: %d / %d; distances equal: %d / %d",
                sum(eq$identical), nrow(eq), sum(eq$n_old == eq$n_offset13, na.rm = TRUE), nrow(eq),
                sum(ld$equal), nrow(ld)))
