# HEAVY-tier wrapper specs — the ones that need a genomics toolchain.
#
# Kept in a side-effect-free file of its own so that BOTH the heavy test
# (tests/heavy/test-heavy-wrappers.R) and the quick-tier registry canary
# (tests/testthat/test-heavy-tier-registry.R) can read the table. The canary is
# the reason this file must stay free of tool calls at source time: it runs in the
# default gate, where none of these tools is exercised.
#
# CONTRACT: the heavy tier is EXPECTED GREEN. It is separated from the quick gate
# by RUNTIME and tool dependency, not by correctness — unlike --invariants, which
# is expected red by design. A missing tool is therefore a FAILURE here, never a
# skip: inside cline-go:latest plink, LEA, vcftools and the EMMAX binaries are all
# present, so their absence means the image is wrong.
#
# Every row that consumes test_data/ must COPY the fixture into its sandbox
# first. test_data/ is tracked (5 files on origin/main) and vcf2lfmm.R derives its
# outputs from the INPUT path — pointed at the tracked copy it would write
# testdata.lfmm straight into the repo.

HEAVY_TESTDATA <- file.path(REPO, "test_data")

# Copy the tracked VCF into d and return the copy's path.
heavy_copy_vcf <- function(d, name = "testdata.vcf") {
    src <- file.path(HEAVY_TESTDATA, "testdata.vcf")
    dst <- file.path(d, name)
    stopifnot(file.exists(src))
    file.copy(src, dst, overwrite = TRUE)
    dst
}

heavy_copy_metadata <- function(d) {
    src <- file.path(HEAVY_TESTDATA, "testdata_metadata.tsv")
    dst <- file.path(d, "metadata.tsv")
    file.copy(src, dst, overwrite = TRUE)
    dst
}

# Number of samples in a VCF, read off the #CHROM header.
heavy_vcf_n_samples <- function(vcf) {
    hdr <- system2("grep", c("-m1", shQuote("^#CHROM"), shQuote(vcf)), stdout = TRUE)
    length(strsplit(hdr, "\t")[[1]]) - 9L
}

# Build everything emmax.R needs: TPED/TFAM via plink, BN kinship via the
# preflight-guarded emmax-kin, a trait table and a LEA-shaped covariates file.
# This is the most expensive fixture in the suite and the reason emmax.R is heavy
# rather than quick — but it is also the only test anywhere that exercises the
# plink -> emmax-kin -> emmax chain end to end.
heavy_build_emmax_inputs <- function(d) {
    vcf    <- heavy_copy_vcf(d)
    prefix <- file.path(d, "emmax_in")

    plink_log <- system2("plink", c("--vcf", shQuote(vcf), "--allow-extra-chr",
                                    "--recode12", "transpose",
                                    "--output-missing-genotype", "0",
                                    "--out", shQuote(prefix)),
                         stdout = TRUE, stderr = TRUE)
    tfam <- paste0(prefix, ".tfam")
    if (!file.exists(tfam)) {
        stop("plink did not write ", tfam, "\n", paste(plink_log, collapse = "\n"))
    }
    # Mirror gea.smk:71-72's FID/IID collapse, so the TFAM matches production.
    awk <- paste0("awk '{split($1,a,\"_\"); split($2,b,\"_\"); ",
                  "if(a[1]==b[1]){$1=a[1];$2=a[1]} print}' ",
                  shQuote(tfam), " > ", shQuote(paste0(prefix, "_tmp.tfam")))
    system(awk)
    file.rename(paste0(prefix, "_tmp.tfam"), tfam)

    kin_log <- system2(file.path(SCRIPTS, "emmax_run.sh"),
                       c(shQuote(file.path(SCRIPTS, "emmax-kin-intel64")),
                         "-v", "-d", "10", "-x", shQuote(prefix)),
                       stdout = TRUE, stderr = TRUE)
    kinship <- paste0(prefix, ".aBN.kinf")
    if (!file.exists(kinship)) {
        stop("emmax-kin did not write ", kinship, "\n", paste(kin_log, collapse = "\n"))
    }

    n <- heavy_vcf_n_samples(vcf)
    trait_file <- file.path(d, "trait.tsv")
    data.table::fwrite(
        data.table::data.table(bio_1 = as.numeric(seq_len(n)) + 0.5),
        trait_file, sep = "\t")

    # LEA-shaped projections: plain space-separated numerics, one row per sample.
    covar <- file.path(d, "projections.txt")
    m <- matrix(round(stats::rnorm(n * 3), 4), nrow = n)
    utils::write.table(m, covar, sep = " ", row.names = FALSE, col.names = FALSE)

    samples_file <- file.path(d, "samples.list")
    writeLines(sprintf("S%03d", seq_len(n)), samples_file)

    inter <- file.path(d, "inter")
    dir.create(inter, showWarnings = FALSE)

    list(vcf = vcf, trait = trait_file, covar = covar, samples = samples_file,
         inter = paste0(inter, "/"), prefix = prefix, kinship = kinship)
}

HEAVY_WRAPPERS <- list(
    list(
        label   = "vcf_to_gapit_numeric.R (vcfR)",
        script  = "vcf_to_gapit_numeric.R",
        tools   = character(0),          # pure R: vcfR, no external binary
        build   = function(d) heavy_copy_vcf(d),
        args    = function(d, vcf) c(vcf, file.path(d, "GD.txt"),
                                     file.path(d, "GM.txt")),
        outputs = function(d) file.path(d, c("GD.txt", "GM.txt"))
    ),
    list(
        # The hermeticity demo of the set: vcf2lfmm() derives <base>.lfmm and
        # <base>.lfmm_nmissing from the INPUT path, so `outputs` must name what is
        # written, not what is passed, and the fixture must be a copy.
        label   = "vcf2lfmm.R (LEA)",
        script  = "vcf2lfmm.R",
        tools   = character(0),          # LEA is an R package
        build   = function(d) heavy_copy_vcf(d, "geno.vcf"),
        args    = function(d, vcf) c(vcf, "TRUE"),
        outputs = function(d) file.path(d, c("geno.lfmm", "geno.lfmm_nmissing"))
    ),
    list(
        label   = "ld_decay_prepare.R (plink + shell)",
        script  = "ld_decay_prepare.R",
        tools   = "plink",
        build   = function(d) {
            vcf  <- heavy_copy_vcf(d)
            meta <- heavy_copy_metadata(d)
            dir.create(file.path(d, "sample_lists"), showWarnings = FALSE)
            dir.create(file.path(d, "chr_vcfs"), showWarnings = FALSE)
            c(vcf, meta)
        },
        args    = function(d, f) c(f[1], f[2], "site", "3", "genome_wide",
                                   file.path(d, "sample_lists"),
                                   file.path(d, "chr_vcfs"), "NULL"),
        outputs = function(d) character(0),
        out_min = function(d) list(dir = file.path(d, "sample_lists"),
                                   pattern = "\\.", n = 1)
    ),
    list(
        # The full plink -> emmax-kin -> emmax chain. emmax_run.sh's preflight
        # (ELF e_machine vs uname -m, then an emulation probe on /proc/cpuinfo)
        # passes on a native x86-64 Linux host and correctly REFUSES on
        # amd64-under-emulation (mac-studio), where the Intel MKL topology probe
        # hangs forever. So this row is expected to fail fast, with a readable
        # reason, on an Apple-silicon host — that is the wrapper working.
        label   = "emmax.R (plink + emmax-kin + emmax)",
        script  = "emmax.R",
        tools   = c("plink", file.path(SCRIPTS, "emmax-intel64"),
                    file.path(SCRIPTS, "emmax-kin-intel64"),
                    file.path(SCRIPTS, "emmax_run.sh")),
        build   = function(d) heavy_build_emmax_inputs(d),
        args    = function(d, f) c(f$vcf, "3", f$trait, f$covar, "bio_1",
                                   f$inter, f$samples,
                                   file.path(d, "EMMAX_pvalues.tsv"),
                                   f$prefix, f$kinship),
        outputs = function(d) file.path(d, "EMMAX_pvalues.tsv")
    )
)
