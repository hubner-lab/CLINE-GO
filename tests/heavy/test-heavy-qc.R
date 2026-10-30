# Processing QC against the real tools (release Phase 4c): the shell recipes of
# compute_depth_stats (7ca55e) and filter_vcf_steps (7e11b6) in workflow/rules/processing.smk,
# copied verbatim, on fixtures whose answers are known by construction.

VCFTOOLS <- "/app/vcftools-0.1.16/bin/vcftools"

test_that("compute_depth_stats reads each genotype's DP from the raw VCF for the keep-list's samples", {
    # Until 2026-10-30 the rule read W['vcf_filt'], which plink's --recode vcf writes with GT
    # only: vcftools then found no DP (N_SITES 0, MEAN_DEPTH -nan) and the rule exited 0.
    # Known depths: A has a DP at one site (10); C at three (30, 6, 8); B (20, 0, 4) is not in
    # the keep-list. Site 300 carries no DP; a '.' DP is skipped.
    require_tools(c("plink", VCFTOOLS))
    d <- withr::local_tempdir()
    vcf <- file.path(d, "raw.vcf")
    writeLines(c(
        "##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000>",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC",
        "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:10\t0/0:20\t1/1:30",
        "chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:.\t./.:0\t1/1:6",
        "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1",
        "chr1\t400\t.\tA\tG\t.\tPASS\t.\tGT:DP\t./.:.\t0/1:4\t0/0:8"), vcf)
    samples <- file.path(d, "samples.list")            # extract_samples' "0 IID" lines
    writeLines(c("0 A", "0 C"), samples)
    prefix <- file.path(d, "depth_tmp")

    # Verbatim shape of rule compute_depth_stats.
    system2("sh", c("-c", shQuote(sprintf("awk '{print $2}' %s > %s.keep", samples, prefix))))
    for (opt in c("--depth", "--site-mean-depth"))
        system2(VCFTOOLS, c("--vcf", vcf, "--keep", paste0(prefix, ".keep"), opt, "--out", prefix),
                stdout = FALSE, stderr = FALSE)
    idepth <- data.table::fread(paste0(prefix, ".idepth"))
    expect_identical(idepth$INDV, c("A", "C"))
    expect_identical(idepth$N_SITES, c(1L, 3L))
    expect_equal(idepth$MEAN_DEPTH, c(10, 44 / 3), tolerance = 1e-4)
    site <- readLines(paste0(prefix, ".ldepth.mean"))
    expect_identical(vapply(strsplit(site[-1], "\t"), `[[`, "", 3L), c("20", "6", "-nan", "8"))

    # The plink list fed as is: vcftools takes a line's FIRST field as the ID, keeps nobody
    # and exits 1 (system2 warns about the status, which is the point here).
    as_is <- suppressWarnings(system2(VCFTOOLS, c("--vcf", vcf, "--keep", samples, "--depth",
                                                  "--out", file.path(d, "as_is")),
                                      stdout = TRUE, stderr = TRUE))
    expect_identical(attr(as_is, "status"), 1L)
    expect_true(any(grepl("kept 0 out of 3 Individuals", as_is)))

    # What the rule used to read: filter_vcf's plink output has GT only, so no depth at all.
    writeLines(paste(0, c("A", "B", "C")), file.path(d, "all.list"))
    system2("plink", c("--vcf", vcf, "--const-fid", "--chr-set", "95", "--allow-extra-chr",
                       "--output-chr", "MT", "--set-missing-var-ids", "@:#",
                       "--keep", file.path(d, "all.list"), "--maf", "0.01", "--geno", "0.5",
                       "--recode", "vcf", "--out", file.path(d, "filt")),
            stdout = FALSE, stderr = FALSE)
    expect_false(any(grepl("ID=DP", readLines(file.path(d, "filt.vcf")))))
    system2(VCFTOOLS, c("--vcf", file.path(d, "filt.vcf"), "--depth", "--out", file.path(d, "filt")),
            stdout = FALSE, stderr = FALSE)
    expect_true(all(data.table::fread(file.path(d, "filt.idepth"))$N_SITES == 0L))
})

test_that("filter_vcf_steps' plink report gives filter_vcf's own step counts, --geno before --maf", {
    # 7e11b6: filtering_summary.tsv's "After MAF filter" was the raw SNP count minus every raw
    # SNP below the MAF threshold, a step plink never takes. Seven SNPs on 10 samples, at
    # --geno 0.2 --maf 0.1:
    #   s1 s5 s6  pass (s6 misses 1 of 10 calls)
    #   s2 s7     missing 4 and 3 of 10: removed by --geno
    #   s3        MAF 0.05: removed by --maf
    #   s4        missing 4 of 10 AND MAF 1/12: removed by --geno, so --maf never sees it
    # plink: 7 loaded, 3 removed by --geno, 1 by --maf, 3 pass. The old row read 7 - 2 = 5
    # (s3 and s4 below 0.1 on the raw table); the real --geno step leaves 4.
    require_tools("plink")
    d <- withr::local_tempdir()
    gt <- function(...) paste(c(...), collapse = "\t")
    rep_gt <- function(...) { x <- c(...); rep(names(x), x) }
    snps <- list(
        s1 = rep_gt("0/0" = 3, "0/1" = 4, "1/1" = 3),
        s2 = rep_gt("./." = 4, "0/1" = 3, "0/0" = 3),
        s3 = rep_gt("0/1" = 1, "0/0" = 9),
        s4 = rep_gt("./." = 4, "0/1" = 1, "0/0" = 5),
        s5 = rep_gt("0/1" = 6, "0/0" = 4),
        s6 = rep_gt("./." = 1, "0/1" = 4, "1/1" = 2, "0/0" = 3),
        s7 = rep_gt("./." = 3, "0/1" = 5, "0/0" = 2))
    ids <- sprintf("S%02d", 1:10)
    vcf <- file.path(d, "raw.vcf")
    writeLines(c("##fileformat=VCFv4.2",
                 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                 gt("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ids),
                 vapply(seq_along(snps), function(i)
                     gt("chr1", i * 100, names(snps)[i], "A", "G", ".", "PASS", ".", "GT", snps[[i]]),
                     "")), vcf)
    keep <- file.path(d, "keep.list")
    writeLines(paste(0, ids), keep)
    flags <- c("--vcf", vcf, "--const-fid", "--chr-set", "95", "--allow-extra-chr",
               "--output-chr", "MT", "--set-missing-var-ids", "@:#", "--keep", keep,
               "--maf", "0.1", "--geno", "0.2")
    # Verbatim shape of rule filter_vcf_steps, and of rule filter_vcf on the same flags.
    system2("plink", c(flags, "--write-snplist", "--out", file.path(d, "steps")),
            stdout = FALSE, stderr = FALSE)
    system2("plink", c(flags, "--recode", "vcf", "--out", file.path(d, "filt")),
            stdout = FALSE, stderr = FALSE)

    p <- parse_plink_filter_log(readLines(file.path(d, "steps.log")), geno = 0.2)
    expect_identical(p, list(n_loaded = 7, n_removed_geno = 3, n_removed_maf = 1,
                             n_pass = 3, n_samples_pass = 10))
    filt <- readLines(file.path(d, "filt.vcf"))
    body <- filt[!startsWith(filt, "#")]
    expect_identical(length(body), 3L)
    expect_identical(vapply(strsplit(body, "\t"), `[[`, "", 3L), c("s1", "s5", "s6"))
    expect_identical(readLines(file.path(d, "steps.snplist")), c("s1", "s5", "s6"))

    a <- filtering_attrition(10, 7, list("After sample missingness filter" = 10L), p, 3)
    expect_identical(a[stage == "After SNP missingness filter", n_snps], 4)
    expect_identical(a[stage == "After MAF + SNP missingness filter", n_snps], 3)
})
