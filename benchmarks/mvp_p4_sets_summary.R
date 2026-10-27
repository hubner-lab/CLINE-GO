#!/usr/bin/env Rscript
# =============================================================================
# mvp_p4_sets_summary.R -- snp_sets_summary.tsv for an offset arm whose candidate sets are the
# framework's declared sets (MVP<seed>hd trees, LD-window plan Phase 4), in the schema
# mvp_build_snp_sets.R wrote for offset13. Read by mvp_panel_rank.R (usable_for_offset of an
# empty panel) and mvp_ms_sim_figures.R (set "all": tested SNPs and causal loci per replicate).
#
#   sets     every row of <outdir>/panel_pr_recomputed.tsv (mvp_panel_tables.R, MVP_PANELS=hd):
#            n_linked is option B there when MVP_LINKED_BP was set; n_background = n - causal -
#            linked; status "framework"; untracked 0.
#   all      copied from <ref>/snp_sets_summary.tsv: the hd trees carry the old p-value tables
#            (copied, and identical on the 5 refitted pilot seeds), so the tested SNPs are the
#            same; status "unchanged".
#
#   Rscript mvp_p4_sets_summary.R <outdir> <ref_dir>
# =============================================================================
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 2L) stop("usage: mvp_p4_sets_summary.R <outdir> <ref_dir>")
OUT <- a[1]; REF <- a[2]
PR  <- fread(file.path(OUT, "panel_pr_recomputed.tsv"), colClasses = c(seed = "character"))
OLD <- fread(file.path(REF, "snp_sets_summary.tsv"), colClasses = c(seed = "character"))
co  <- unique(OLD[, .(seed, cohort)])
stopifnot(!anyDuplicated(co$seed), all(PR$seed %in% co$seed))
S <- PR[, .(seed, set, n_snps = n, n_causal, n_linked, n_background = n - n_causal - n_linked,
            status = "framework", usable_for_offset, untracked = 0L)]
stopifnot(all(S$n_background >= 0))
ALL <- OLD[set == "all" & seed %in% PR$seed]
stopifnot(nrow(ALL) == uniqueN(PR$seed))
S <- rbind(merge(S, co, by = "seed")[, names(OLD), with = FALSE], ALL)
setorder(S, seed, set)
fwrite(S, file.path(OUT, "snp_sets_summary.tsv"), sep = "\t")
message(sprintf("wrote %s/snp_sets_summary.tsv -- %d rows, %d seeds, sets: %s", OUT, nrow(S),
                uniqueN(S$seed), paste(sort(unique(S$set)), collapse = ",")))
print(S[, .(n = .N, usable = sum(usable_for_offset), empty = sum(n_snps == 0)), by = set])
