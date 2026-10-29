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

test_that("snmf.R is reproducible run to run, without collapsing its repetitions", {
    # LEA::snmf() draws its OWN seed at random unless given seed=; R's set.seed()
    # never reaches its C RNG. Without the pass-through, two runs of identical inputs
    # produced different Q-matrices, a different cross-entropy curve (hence a
    # different k_best), different sNMF-driven imputation and therefore different
    # structure-corrected p-values everywhere downstream (finding 8d07b1).
    #
    # Both halves matter: a single seed must make the run reproducible AND must not
    # collapse the `repetitions` into one answer, because best-run selection by
    # cross-entropy depends on them differing.
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    require_tools("Rscript")
    d <- withr::local_tempdir()

    # A .geno is one line per SNP, one character per individual (0/1/2, 9 = NA).
    set.seed(7)
    n_ind <- 24L; n_snp <- 60L
    lines <- vapply(seq_len(n_snp), function(i)
        paste(sample(c(0L, 1L, 2L), n_ind, replace = TRUE), collapse = ""),
        character(1))

    run_once <- function(tag) {
        sub <- file.path(d, tag)
        dir.create(sub, showWarnings = FALSE)
        geno <- file.path(sub, "t.geno")
        writeLines(lines, geno)
        res <- run_wrapper("snmf.R", c(geno, "2", "3", "2", "2", "1", "new"))
        expect_identical(res$status, 0L, info = res$output)
        # The seed actually used is logged, so a run is auditable from its log.
        expect_match(res$output, "sNMF seed = 42")
        sub
    }

    a <- run_once("a")
    b <- run_once("b")

    q_of <- function(sub) {
        f <- sort(list.files(sub, pattern = "\\.Q$", recursive = TRUE,
                             full.names = TRUE))
        expect_gt(length(f), 1L)
        lapply(f, function(x) as.matrix(utils::read.table(x)))
    }
    qa <- q_of(a); qb <- q_of(b)
    expect_identical(length(qa), length(qb))
    for (i in seq_along(qa)) expect_equal(qa[[i]], qb[[i]], tolerance = 0)

    # Not collapsed: at a given K the repetitions must not all be the same matrix.
    k3 <- sort(list.files(a, pattern = "\\.3\\.Q$", recursive = TRUE,
                          full.names = TRUE))
    skip_if(length(k3) < 2, "fewer than 2 repetitions written")
    m <- lapply(k3, function(x) as.matrix(utils::read.table(x)))
    expect_false(isTRUE(all.equal(m[[1]], m[[2]], tolerance = 1e-10)))
})

test_that("snmf.R gives the same Q matrices and cross-entropies at 1, 4 and 8 threads", {
    # snmf.R hands LEA::snmf() CPU = {threads}, and cpu defaults to the host's cores - 2,
    # so the thread count differs between machines and between runs of one project. The
    # seed must make a run independent of it: the Q matrices become clusters_K*.tsv and
    # drive the sNMF imputation, and each run's cross-entropy picks the best run and the K
    # (c0f228). Measured 2026-10-29 in image 005aa13b65f0: identical at 1, 4 and 8 threads
    # and run to run at 4, with or without OPENBLAS_NUM_THREADS=1.
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    require_tools("Rscript")
    d <- withr::local_tempdir()

    # Three source populations, so the runs have structure to converge on.
    set.seed(11)
    n_ind <- 60L; n_snp <- 400L
    pop <- rep(1:3, each = n_ind / 3)
    freqs <- matrix(stats::runif(3 * n_snp, 0.05, 0.95), 3)
    lines <- vapply(seq_len(n_snp), function(j)
        paste(stats::rbinom(n_ind, 2, freqs[pop, j]), collapse = ""), character(1))

    run_at <- function(tag, cpu) {
        sub <- file.path(d, tag)
        dir.create(sub)
        geno <- file.path(sub, "t.geno")
        writeLines(lines, geno)
        res <- run_wrapper("snmf.R", c(geno, "2", "4", "2", "3", cpu, "new"))
        expect_identical(res$status, 0L, info = res$output)
        runs <- sort(list.files(sub, "\\.snmfClass$", recursive = TRUE, full.names = TRUE))
        expect_length(runs, 9L)   # K 2..4 x 3 repetitions
        ce <- vapply(runs, function(f) {
            txt <- paste(readLines(f), collapse = " ")
            as.numeric(sub(".*crossEntropy = ([-0-9.eE+]+).*", "\\1", txt))
        }, numeric(1), USE.NAMES = FALSE)
        expect_false(anyNA(ce))
        list(run = substring(runs, nchar(sub) + 2L), ce = ce,
             q = lapply(sub("\\.snmfClass$", ".Q", runs),
                        function(f) as.matrix(utils::read.table(f))))
    }

    one   <- run_at("cpu1", "1")
    four  <- run_at("cpu4", "4")
    again <- run_at("cpu4_again", "4")
    eight <- run_at("cpu8", "8")
    for (other in list(again, one, eight)) {
        expect_identical(other$run, four$run)
        expect_identical(other$ce, four$ce)
        expect_identical(other$q, four$q)
    }
})

test_that("the EMMAX TPED rules keep the VCF's full sample names as FID and IID", {
    # plink's --vcf importer splits a sample ID on '_': WBDC_001 becomes
    # FID=WBDC / IID=001, and an ID with TWO underscores is a hard error
    # ("Multiple instances of '_' in sample ID"). emmax.R:65 reads the TFAM's
    # first two columns, but emmax_phenotypes.R:91,96 keys its phen/covar rows on
    # the FULL VCF sample name, so the split made the two sides unjoinable for any
    # underscore-bearing scheme (finding d8f924). The five tped rules
    # (gea.smk:69,113 gwas.smk:33,188 pregea.smk:71) therefore pass --double-id,
    # which is what the removed TFAM-fixup awk was trying to reproduce.
    #
    # SIMDATA's own IDs carry no underscore, which is the one case that always
    # worked and the reason this went unnoticed — so this is the only place the
    # WBDC/PG-style schemes are exercised.
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    require_tools("plink")
    d <- withr::local_tempdir()

    write_vcf <- function(path, samples) {
        writeLines(c(
            "##fileformat=VCFv4.2", "##contig=<ID=1>",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                    "INFO", "FORMAT", samples), collapse = "\t"),
            paste(c("1", 100, ".", "A", "G", ".", "PASS", ".", "GT",
                    rep(c("0/0", "0/1", "1/1"), length.out = length(samples))),
                  collapse = "\t"),
            paste(c("1", 200, ".", "C", "T", ".", "PASS", ".", "GT",
                    rep(c("0/1", "1/1", "0/0"), length.out = length(samples))),
                  collapse = "\t")), path)
        path
    }

    # Verbatim shape of the five tped rules.
    run_tped <- function(vcf, prefix) {
        system2("plink", c("--vcf", vcf, "--double-id", "--chr-set", "95", "--allow-extra-chr",
                           "--output-chr", "MT", "--recode12", "transpose",
                           "--output-missing-genotype", "0", "--out", prefix),
                stdout = FALSE, stderr = FALSE)
        tfam <- paste0(prefix, ".tfam")
        expect_true(file.exists(tfam))
        f <- data.table::fread(tfam, header = FALSE,
                               select = list(character = 1:2))
        data.table::setnames(f, c("FID", "IID"))
        f
    }

    for (nm in c("no_underscore", "one_underscore", "two_underscores")) {
        samples <- switch(nm,
            no_underscore   = c("NEG01", "NEG02", "NEG03"),
            one_underscore  = c("WBDC_001", "WBDC_002", "WBDC_003"),
            two_underscores = c("PG_B_001", "PG_B_002", "PG_B_003"))
        vcf <- write_vcf(file.path(d, paste0(nm, ".vcf")), samples)
        f <- run_tped(vcf, file.path(d, nm))
        # FID == IID == the VCF sample name, verbatim, for every scheme.
        expect_identical(f$IID, samples, info = nm)
        expect_identical(f$FID, samples, info = nm)
    }

    # And the half that actually broke: emmax_phenotypes.R writes FID/IID from the
    # VCF sample names, so those rows must join the TFAM on both keys.
    f <- run_tped(file.path(d, "one_underscore.vcf"), file.path(d, "join"))
    phen <- data.table::data.table(FID = c("WBDC_001", "WBDC_002", "WBDC_003"),
                                   IID = c("WBDC_001", "WBDC_002", "WBDC_003"),
                                   value = c(1.0, 2.0, 3.0))
    joined <- merge(f, phen, by = c("FID", "IID"))
    expect_identical(nrow(joined), 3L)
})

test_that("filter_vcf's plink + sed leave the same contig names normalize_gff.py produces", {
    # The chromosome-name contract has two halves in two languages (audit 2026-09-13
    # B1): plink re-emits the VCF through its own parser (Chr1/CHR1 -> 1, X/Y/MT as
    # LETTERS only because of --output-chr MT, and a contig numbered 23-26 as a NUMBER only
    # because of --chr-set 95 — the human default reads 23/26 as X/MT, 12f364), the rule's
    # sed then strips any-case `chr` from what plink passed through verbatim (chr2H, Chr5H),
    # and scripts/normalize_gff.py applies the same strip to the GFF. SIMDATA is
    # lowercase `chr1..chr5`, the one case that always worked, so this is the only place
    # the mixed-case / sex-chromosome / >22-chromosome path is exercised.
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    require_tools(c("plink", "python3"))
    d <- withr::local_tempdir()
    vcf <- file.path(d, "raw.vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##contig=<ID=Chr1>", "##contig=<ID=CHR2>", "##contig=<ID=chr3>",
        "##contig=<ID=Chr5H>", "##contig=<ID=chr2H>", "##contig=<ID=23>", "##contig=<ID=chr26>",
        "##contig=<ID=X>", "##contig=<ID=MT>",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
        paste(c("Chr1", "CHR2", "chr3", "Chr5H", "chr2H", "23", "chr26", "X", "MT"),
              seq(100, 900, by = 100), ".", "A", "G", ".", "PASS", ".", "GT",
              "0/0", "0/1", "1/1", "0/1", sep = "\t")), vcf)
    keep <- file.path(d, "keep.txt")
    writeLines(paste(0, c("S1", "S2", "S3", "S4")), keep)
    prefix <- file.path(d, "filt")
    # Verbatim shape of rule filter_vcf (workflow/rules/processing.smk).
    system2("plink", c("--vcf", vcf, "--const-fid", "--chr-set", "95", "--allow-extra-chr",
                       "--output-chr", "MT",
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
                     c("1", "2", "3", "5H", "2H", "23", "26", "X", "MT"))
    hdr_ids <- sub("^##contig=<ID=([^,>]+).*$", "\\1", grep("^##contig", lines, value = TRUE))
    expect_setequal(hdr_ids, c("1", "2", "3", "5H", "2H", "23", "26", "X", "MT"))

    # The GFF side, on the same raw names, must land on the same set.
    gff <- file.path(d, "raw.gff3")
    writeLines(c("##gff-version 3",
                 paste(c("Chr1", "CHR2", "chr3", "Chr5H", "chr2H", "23", "Chr26", "X", "M"),
                       "src", "gene", 1, 50, ".", "+", ".", "ID=g", sep = "\t")), gff)
    norm <- file.path(d, "normalized.gff3")
    res <- system2("python3", c(file.path(SCRIPTS, "normalize_gff.py"), gff, out, norm),
                   stdout = TRUE, stderr = TRUE)
    status <- attr(res, "status"); if (is.null(status)) status <- 0L
    expect_identical(status, 0L, info = paste(res, collapse = "\n"))
    gl <- readLines(norm); gl <- gl[!startsWith(gl, "#")]
    expect_identical(vapply(strsplit(gl, "\t"), `[[`, "", 1L),
                     c("1", "2", "3", "5H", "2H", "23", "26", "X", "MT"))
    expect_true(any(grepl("9 chromosome name\\(s\\) shared", res)))
})

test_that("lfmm.R stops when the vcfsnp does not label every tested SNP column", {
    # Since sNMF.imputed_maf, the matrix LFMM tests and the vcfsnp that names its columns
    # come from filter_imputed_maf.R as a pair. A rule pairing the filtered matrix with the
    # UNFILTERED vcfsnp used to run clean: data.table recycles the shorter p-value vector
    # into the SNP table when the lengths divide, so every label after the first removed
    # SNP was silently wrong. The guard fires before LEA fits anything, so this is cheap;
    # it lives here only because lfmm.R attaches LEA and WGCNA at the top.
    skip_if_not(dir.exists(SCRIPTS), "scripts/ not mounted")
    d <- withr::local_tempdir()
    f <- fx_imputed_triple(d)                      # 25 samples x 6 SNP columns
    writeLines(readLines(f$vcfsnp)[-6], f$vcfsnp)   # 5 labels
    clim <- write_tsv(data.table::data.table(bio_1 = seq_len(f$n_ind) / 10),
                      file.path(d, "climate.tsv"))
    out <- file.path(d, "lfmm_pvalues_bio_1.tsv")
    res <- run_wrapper("lfmm.R", c(f$lfmm, f$lfmm, clim, "2", "bio_1", f$vcfsnp, out), wd = d)

    expect_gt(res$status, 0L)
    expect_match(res$output, "has 6 SNP columns but .* has 5 rows")
    expect_false(file.exists(out))
})
