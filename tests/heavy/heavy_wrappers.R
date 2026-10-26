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

# A VCF for the LD-decay scripts that PLINK's default human chromosome mode would
# break: 26 numeric chromosomes (23-26 read as X/Y/XY/MT), heterozygous calls and ~15%
# missing calls. `chr_names` replaces the numeric names (e.g. X, MT and Pt contigs).
# Returns the paths and the dosage matrix (samples x SNPs), so a test can recompute
# every r^2 in R. Seeded; written into d.
heavy_ld_fixture <- function(d, n = 30L, n_chr = 26L, per_chr = 15L, seed = 11L,
                             chr_names = as.character(seq_len(n_chr))) {
    set.seed(seed)
    info <- data.table::rbindlist(lapply(chr_names, function(ch) data.table::data.table(
        chr = ch, pos = sort(sample.int(20000L, per_chr)))))
    info[, id := paste0("s", .I)]
    X <- matrix(NA_real_, n, nrow(info))
    for (ch in unique(info$chr)) {
        base <- stats::rbinom(n, 2, 0.4)
        for (k in which(info$chr == ch)) {
            g <- ifelse(stats::runif(n) < 0.7, base, stats::rbinom(n, 2, stats::runif(1, 0.1, 0.5)))
            g[stats::runif(n) < 0.15] <- NA
            X[, k] <- g
        }
    }
    gt <- ifelse(is.na(X), "./.", c("0/0", "0/1", "1/1")[X + 1])
    smp <- sprintf("p%02d", seq_len(n))
    vcf <- file.path(d, "ld.vcf")
    con <- file(vcf, "w")
    writeLines(c("##fileformat=VCFv4.2",
                 paste0("##contig=<ID=", unique(info$chr), ">"),
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                 paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                         "FORMAT", smp), collapse = "\t")), con)
    for (i in seq_len(nrow(info))) {
        writeLines(paste(c(info$chr[i], info$pos[i], ".", "A", "G", ".", "PASS", ".", "GT",
                           gt[, i]), collapse = "\t"), con)
    }
    close(con)
    sl <- file.path(d, "sample_lists")
    dir.create(sl, showWarnings = FALSE)
    writeLines(smp, file.path(sl, "All.txt"))
    writeLines(smp[1:15], file.path(sl, "G_1.txt"))
    data.table::fwrite(data.table::data.table(group = c("All", "G 1"), n_samples = c(n, 15L),
                                              list_file = c("All", "G_1")),
                       file.path(sl, "manifest.tsv"), sep = "\t")
    list(vcf = vcf, sl = sl, info = info, X = X, samples = smp)
}

# r^2 of every within-chromosome pair at 0 < distance <= max_dist, recomputed in R as
# the squared correlation with missing calls excluded pairwise, summed per distance
# (genome-wide). What ld_decay_pairs.R must reproduce from PLINK.
heavy_ld_expected <- function(fx, max_dist, rows = seq_len(nrow(fx$X))) {
    out <- list()
    for (ch in unique(fx$info$chr)) {
        k <- which(fx$info$chr == ch)
        for (a in seq_along(k)) for (b in seq_along(k)) if (a < b) {
            dd <- abs(fx$info$pos[k[b]] - fx$info$pos[k[a]])
            if (dd == 0 || dd > max_dist) next
            r <- suppressWarnings(stats::cor(fx$X[rows, k[a]], fx$X[rows, k[b]],
                                             use = "pairwise.complete.obs"))
            if (is.finite(r)) out[[length(out) + 1]] <- data.table::data.table(dist_bp = dd, r2 = r^2)
        }
    }
    data.table::rbindlist(out)[, .(sum_r2 = sum(r2), n_pairs = .N), keyby = dist_bp]
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
                                   file.path(d, "chr_vcfs"), "NULL", "NULL"),
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
    ),
    list(
        # Background LD of unlinked markers: plink --r2 inter-chr on a seeded sample
        # of markers per chromosome, for All and one sub-group. testdata.vcf has five
        # chromosomes, so both rows get unlinked pairs and a finite critical r2.
        label   = "ld_decay_background.R (plink inter-chr r2)",
        script  = "ld_decay_background.R",
        tools   = c("plink", "bcftools"),
        build   = function(d) {
            vcf <- heavy_copy_vcf(d)
            sl  <- file.path(d, "sample_lists")
            dir.create(sl, showWarnings = FALSE)
            hdr <- system2("grep", c("-m1", shQuote("^#CHROM"), shQuote(vcf)), stdout = TRUE)
            smp <- strsplit(hdr, "\t")[[1]][-(1:9)]
            writeLines(smp, file.path(sl, "All.txt"))
            writeLines(smp[seq(1, length(smp), 2)], file.path(sl, "G1.txt"))
            data.table::fwrite(data.table::data.table(
                group = c("All", "G1"),
                n_samples = c(length(smp), length(seq(1, length(smp), 2)))),
                file.path(sl, "manifest.tsv"), sep = "\t")
            list(vcf = vcf, sl = sl)
        },
        args    = function(d, f) c(f$vcf, f$sl, "20", "0.95", "1",
                                   file.path(d, "ld_decay_background.tsv"), "NULL"),
        outputs = function(d) file.path(d, "ld_decay_background.tsv")
    ),
    list(
        # LD-decay pairs: PLINK --r2 within 20 kb on a 26-chromosome fixture with
        # heterozygous and missing calls, summed per distance for All and one group
        # whose name is not file-safe. Content is asserted in test-heavy-ld-decay.R.
        label   = "ld_decay_pairs.R (plink r2 per distance, 26 chromosomes)",
        script  = "ld_decay_pairs.R",
        tools   = c("plink", "bcftools"),
        build   = function(d) heavy_ld_fixture(d),
        args    = function(d, f) c(f$vcf, f$sl, "20", "1e6", "1", "both",
                                   file.path(d, "stat_gw"), file.path(d, "stat_chr"),
                                   file.path(d, "pairs_summary.tsv"), "NULL"),
        outputs = function(d) file.path(d, c("pairs_summary.tsv", "stat_gw/All.stat.gz",
                                             "stat_gw/G_1.stat.gz", "stat_chr/All_26.stat.gz"))
    )
)
