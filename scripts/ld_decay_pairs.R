#!/usr/bin/env Rscript
# LD-decay pairs per group: PLINK 1.9 --r2 (squared allele-count correlation, missing
# calls excluded pairwise) between SNPs of one chromosome within max_distance, summed
# per bp distance into PopLDdecay-layout .stat.gz files that ld_decay_analyze.R reads:
# stat_gw/<stem>.stat.gz and, when the scope includes per_chromosome,
# stat_chr/<stem>_<chr>.stat.gz. <stem> is the group's sample-list name in
# manifest.tsv. When all pairs within the window fit the pair budget the run is
# exhaustive; otherwise PLINK pairs a seeded sample of anchor SNPs with every SNP of
# their window (--ld-snp-list), which thins all distances evenly. SNPs on the contigs
# in LDdecay.exclude_chromosomes are left out (--exclude) and never anchor.
#
# Chromosomes are read with --chr-set 95 --allow-extra-chr: in PLINK's default human
# mode a numeric chromosome 23-26 is X/Y/XY/MT (Y calls of non-males are dropped) and a
# code above 26 is an error. Here 1-95 are autosomes, X/Y/XY/MT take codes 96-99 and
# other names (Pt, scaffolds) pass through. X, MT and Pt give the r^2 of an autosome; Y
# gives no pairs, since a VCF carries no sex and PLINK treats non-males' Y calls as missing.
# (`no-x no-y no-xy no-mt` made PLINK refuse any VCF with such a contig.)
library(data.table)
source("/pipeline/scripts/R/lib/ld_decay.R")

args <- commandArgs(trailingOnly = TRUE)
################
VCF_PATH         <- args[1]  # filtered VCF (observed calls)
SAMPLE_LISTS_DIR <- args[2]  # ld_decay_prepare.R output: manifest.tsv + <stem>.txt
MAX_DIST_KB      <- as.numeric(args[3])
MAX_PAIRS        <- as.numeric(args[4])
SEED             <- as.integer(args[5])
SCOPE            <- args[6]  # genome_wide | per_chromosome | both
STAT_GW_DIR      <- args[7]
STAT_CHR_DIR     <- args[8]
OUT_SUMMARY      <- args[9]  # one row per group: SNPs, anchors, pairs
EXCLUDE          <- args[10] # LDdecay.exclude_chromosomes, comma-separated, or "NULL"
################
stopifnot(length(args) == 10, MAX_DIST_KB > 0, MAX_PAIRS > 0,
          SCOPE %in% c("genome_wide", "per_chromosome", "both"))
PLINK_CHR <- c("--chr-set", "95", "--allow-extra-chr")
max_dist <- MAX_DIST_KB * 1000

manifest <- fread(file.path(SAMPLE_LISTS_DIR, "manifest.tsv"), colClasses = "character")
if (!"list_file" %in% names(manifest)) manifest[, list_file := group]
dir.create(STAT_GW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_CHR_DIR, recursive = TRUE, showWarnings = FALSE)
# ld_decay_analyze.R reads whatever file a group or chromosome has, so every one must
# come from this run: a chromosome without pairs now would keep an earlier run's.
unlink(list.files(c(STAT_GW_DIR, STAT_CHR_DIR), pattern = "\\.stat\\.gz$", full.names = TRUE))
# PLINK's pair table is gigabytes at the default budget: it stays beside the output.
tmp <- ld_scratch_dir(OUT_SUMMARY, "pairs_tmp")

plink <- function(...) {
    a <- c(...)
    out <- suppressWarnings(system2("plink", a, stdout = TRUE, stderr = TRUE))
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) {
        stop("plink failed (status ", status, "): plink ", paste(a, collapse = " "),
             "\n", paste(tail(out, 15), collapse = "\n"))
    }
    invisible(out)
}

all_set <- file.path(tmp, "all")
plink("--vcf", VCF_PATH, "--const-fid", PLINK_CHR, "--vcf-half-call", "missing",
      "--make-bed", "--threads", "1", "--out", all_set)
bim <- fread(paste0(all_set, ".bim"), header = FALSE,
             colClasses = c(V1 = "character", V2 = "character"))
setnames(bim, c("chr", "id", "cm", "pos", "a1", "a2"))
# PLINK renames chromosomes ("chr1" -> "1", X -> 96) and writes them in its own order
# (1, 2, 10, then X/Y/XY/MT, then other names), whatever the VCF's order, so a .bim row
# is never paired with a VCF row. The per-chromosome files must carry the VCF's names,
# which ld_decay_prepare.R lists in chromosomes.txt: PLINK itself maps them, reading a
# one-sample VCF with one record per VCF chromosome, each record's ID the chromosome's
# index, under the same chromosome flags.
vcf_counts <- system(paste("bcftools query -f '%CHROM\\n'", shQuote(VCF_PATH), "| uniq -c"),
                     intern = TRUE)
if (!length(vcf_counts)) stop("bcftools read no variant from ", VCF_PATH)
vcf_counts <- data.table(n    = as.numeric(sub("^\\s*([0-9]+) .*$", "\\1", vcf_counts)),
                         name = sub("^\\s*[0-9]+ ", "", vcf_counts))
vcf_counts <- vcf_counts[, .(n = sum(n)), by = name]
vcf_chr <- vcf_counts$name
probe <- file.path(tmp, "chr_probe")
writeLines(c("##fileformat=VCFv4.2",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                     "p"), collapse = "\t"),
             paste(vcf_chr, "1", paste0("c", seq_along(vcf_chr)), "A", "G", ".", ".", ".",
                   "GT", "0/0", sep = "\t")),
           paste0(probe, ".vcf"))
plink("--vcf", paste0(probe, ".vcf"), "--const-fid", PLINK_CHR, "--make-just-bim",
      "--threads", "1", "--out", probe)
probe_bim <- fread(paste0(probe, ".bim"), header = FALSE, select = 1:2,
                   colClasses = c(V1 = "character", V2 = "character"))
code <- setNames(probe_bim$V1, vcf_chr[as.integer(sub("^c", "", probe_bim$V2))])
if (length(code) != length(vcf_chr)) {
    stop("PLINK read ", length(code), " of the VCF's ", length(vcf_chr), " chromosome names")
}
code <- code[vcf_chr]
merged <- code[duplicated(code) | duplicated(code, fromLast = TRUE)]
if (length(merged)) {
    stop("PLINK merged distinct VCF chromosome names: ",
         paste(names(merged), "->", merged, collapse = ", "))
}
chr_name <- setNames(names(code), code)
bim_n <- table(factor(bim$chr, levels = code))
if (anyNA(match(bim$chr, code)) || any(as.numeric(bim_n) != vcf_counts$n)) {
    stop("PLINK kept ", nrow(bim), " of the VCF's ", sum(vcf_counts$n), " variants")
}
# Internal variant IDs: VCF IDs may be "." or repeated, and --ld-snp-list needs them unique.
bim[, id := paste0("v", .I)]
fwrite(bim, paste0(all_set, ".bim"), sep = "\t", col.names = FALSE)

exclude_args <- character(0)
excluded <- ld_excluded_contigs(vcf_chr, EXCLUDE)
if (length(excluded)) {
    drop <- chr_name[bim$chr] %in% excluded
    if (all(drop)) stop("LDdecay.exclude_chromosomes leaves no SNP")
    writeLines(bim$id[drop], file.path(tmp, "excluded.txt"))
    exclude_args <- c("--exclude", file.path(tmp, "excluded.txt"))
    message("INFO: Left out (LDdecay.exclude_chromosomes): ", sum(drop), " SNPs on ",
            paste(excluded, collapse = ", "))
    bim <- bim[!drop]
}

plan <- ld_anchor_plan(bim$chr, bim$pos, max_dist, MAX_PAIRS)
anchor_args <- character(0)
if (!plan$exhaustive) {
    set.seed(SEED)
    anchors <- bim$id[sort(sample.int(nrow(bim), plan$n_anchors))]
    anchor_file <- file.path(tmp, "anchors.txt")
    writeLines(anchors, anchor_file)
    anchor_args <- c("--ld-snp-list", anchor_file)
}
message("INFO: ", nrow(bim), " SNPs, ", round(plan$total_pairs), " pairs within ",
        MAX_DIST_KB, " kb; ",
        if (plan$exhaustive) "all pairs used" else
            paste0(plan$n_anchors, " anchors (seed ", SEED, "), about ",
                   round(plan$n_anchors * plan$mean_partners), " pairs"))

per_chr <- SCOPE %in% c("per_chromosome", "both")
rows <- list()
for (i in seq_len(nrow(manifest))) {
    grp  <- manifest$group[i]
    stem <- manifest$list_file[i]
    keep_args <- character(0)
    if (grp != "All") {
        smp  <- readLines(file.path(SAMPLE_LISTS_DIR, paste0(stem, ".txt")))
        keep <- file.path(tmp, paste0(stem, ".keep"))
        fwrite(data.table("0", smp), keep, sep = "\t", col.names = FALSE)
        keep_args <- c("--keep", keep)
    }
    out <- file.path(tmp, stem)
    plink("--bfile", all_set, PLINK_CHR, keep_args, exclude_args, "--r2", anchor_args,
          "--ld-window", "999999", "--ld-window-kb", MAX_DIST_KB, "--ld-window-r2", "0",
          "--threads", "1", "--out", out)
    pairs <- fread(paste0(out, ".ld"), select = c("CHR_A", "BP_A", "SNP_A", "BP_B", "SNP_B", "R2"),
                   colClasses = c(CHR_A = "character", SNP_A = "character", SNP_B = "character"))
    unlink(paste0(out, ".ld"))
    st <- ld_pairs_to_stat(pairs, max_dist, dedupe = !plan$exhaustive)
    ld_write_stat(st$gw, file.path(STAT_GW_DIR, paste0(stem, ".stat.gz")))
    if (per_chr) {
        for (ch in unique(st$chr$chr)) {
            ld_write_stat(st$chr[chr == ch],
                          file.path(STAT_CHR_DIR, paste0(stem, "_", chr_name[[ch]], ".stat.gz")))
        }
    }
    rows[[grp]] <- data.table(group = grp, n_snps = nrow(bim), exhaustive = plan$exhaustive,
                              n_anchors = plan$n_anchors, n_pairs = sum(st$gw$n_pairs),
                              max_distance_kb = MAX_DIST_KB, seed = SEED)
    message("INFO: ", grp, ": ", sum(st$gw$n_pairs), " pairs")
}

unlink(tmp, recursive = TRUE)
fwrite(rbindlist(rows), OUT_SUMMARY, sep = "\t")
message("INFO: written ", OUT_SUMMARY)
