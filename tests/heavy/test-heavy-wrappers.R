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
