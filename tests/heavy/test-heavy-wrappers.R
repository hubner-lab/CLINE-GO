# HEAVY tier: wrapper smoke tests that need a genomics toolchain.
#
# Run with:  tests/run_all.sh --heavy
# Or alone:  Rscript /pipeline/tests/run_heavy.R
#
# Expected GREEN. See tests/heavy/heavy_wrappers.R for the contract and for why a
# missing tool fails rather than skips.
#
# Known host limitation, not a test bug: the emmax.R row fails fast on
# amd64-under-emulation (Docker Desktop / colima on Apple silicon), because
# scripts/emmax_run.sh refuses to launch the Intel-MKL binary there — the
# alternative being an unkillable hang that orphans the container. On a native
# x86-64 Linux host it runs.

test_that("the repo and the tracked test dataset are both mounted", {
    expect_true(dir.exists(SCRIPTS),
                info = paste0("scripts/ not found at ", SCRIPTS))
    expect_true(file.exists(file.path(HEAVY_TESTDATA, "testdata.vcf")),
                info = "test_data/testdata.vcf missing — it is tracked, so this means no mount")
})

for (spec in HEAVY_WRAPPERS) {
    local({
        s <- spec
        test_that(paste0("heavy wrapper runs and writes its outputs: ", s$label), {
            skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
            require_tools(s$tools)
            expect_wrapper_ok(s)
        })
    })
}

test_that("vcf2lfmm.R writes beside its INPUT, which is why the fixture is a copy", {
    # Pins the reason test_data/ must never be passed directly: the output path is
    # derived from the input path (vcf2lfmm.R:13-14), so pointing this script at
    # the tracked fixture would write into the repo.
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    vcf <- heavy_copy_vcf(d, "geno.vcf")
    res <- run_wrapper("vcf2lfmm.R", c(vcf, "TRUE"), wd = d)

    expect_identical(res$status, 0L, info = res$output)
    expect_true(file.exists(file.path(d, "geno.lfmm")))
    # And nothing appeared next to the tracked original.
    expect_false(file.exists(file.path(HEAVY_TESTDATA, "testdata.lfmm")))
})

test_that("filter_vcf's plink + sed leave the same contig names normalize_gff.py produces", {
    # The chromosome-name contract has two halves in two languages (audit 2026-09-13
    # B1): plink re-emits the VCF through its own parser (Chr1/CHR1 -> 1, and X/Y/MT as
    # LETTERS only because of --output-chr MT — the default is 23/24/26), the rule's sed
    # then strips any-case `chr` from what plink passed through verbatim (chr2H, Chr5H),
    # and scripts/normalize_gff.py applies the same strip to the GFF. SIMDATA is
    # lowercase `chr1..chr5`, the one case that always worked, so this is the only place
    # the mixed-case / sex-chromosome path is exercised.
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    require_tools(c("plink", "python3"))
    d <- withr::local_tempdir()
    vcf <- file.path(d, "raw.vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##contig=<ID=Chr1>", "##contig=<ID=CHR2>", "##contig=<ID=chr3>",
        "##contig=<ID=Chr5H>", "##contig=<ID=chr2H>", "##contig=<ID=X>", "##contig=<ID=MT>",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
        paste(c("Chr1", "CHR2", "chr3", "Chr5H", "chr2H", "X", "MT"),
              seq(100, 700, by = 100), ".", "A", "G", ".", "PASS", ".", "GT",
              "0/0", "0/1", "1/1", "0/1", sep = "\t")), vcf)
    keep <- file.path(d, "keep.txt")
    writeLines(paste(0, c("S1", "S2", "S3", "S4")), keep)
    prefix <- file.path(d, "filt")
    # Verbatim shape of rule filter_vcf (workflow/rules/processing.smk).
    system2("plink", c("--vcf", vcf, "--const-fid", "--allow-extra-chr", "--output-chr", "MT",
                       "--set-missing-var-ids", "@:#", "--keep", keep, "--maf", "0.01",
                       "--geno", "0.5", "--recode", "vcf", "--out", prefix),
            stdout = FALSE, stderr = FALSE)
    out <- paste0(prefix, ".vcf")
    expect_true(file.exists(out))
    system2("sed", c("-E", "-i", shQuote("/^#/! s/^[Cc][Hh][Rr]//"), out))
    system2("sed", c("-E", "-i", shQuote("s/^(##contig=<ID=)[Cc][Hh][Rr]/\\1/"), out))

    lines <- readLines(out)
    body  <- lines[!startsWith(lines, "#")]
    expect_identical(vapply(strsplit(body, "\t"), `[[`, "", 1L),
                     c("1", "2", "3", "5H", "2H", "X", "MT"))
    hdr_ids <- sub("^##contig=<ID=([^,>]+).*$", "\\1", grep("^##contig", lines, value = TRUE))
    expect_setequal(hdr_ids, c("1", "2", "3", "5H", "2H", "X", "MT"))

    # The GFF side, on the same raw names, must land on the same set.
    gff <- file.path(d, "raw.gff3")
    writeLines(c("##gff-version 3",
                 paste(c("Chr1", "CHR2", "chr3", "Chr5H", "chr2H", "X", "M"),
                       "src", "gene", 1, 50, ".", "+", ".", "ID=g", sep = "\t")), gff)
    norm <- file.path(d, "normalized.gff3")
    res <- system2("python3", c(file.path(SCRIPTS, "normalize_gff.py"), gff, out, norm),
                   stdout = TRUE, stderr = TRUE)
    status <- attr(res, "status"); if (is.null(status)) status <- 0L
    expect_identical(status, 0L, info = paste(res, collapse = "\n"))
    gl <- readLines(norm); gl <- gl[!startsWith(gl, "#")]
    expect_identical(vapply(strsplit(gl, "\t"), `[[`, "", 1L),
                     c("1", "2", "3", "5H", "2H", "X", "MT"))
    expect_true(any(grepl("7 chromosome name\\(s\\) shared", res)))
})
