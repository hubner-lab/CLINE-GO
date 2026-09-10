# Shared pipeline libraries, loaded into the package namespace at LOAD time.
#
# ── Why this file exists ──────────────────────────────────────────────────────
# The app reuses two pipeline libraries so Shiny and Snakemake apply identical
# logic (same clustering, same p-value thresholds). They live outside the
# package, under /pipeline/scripts/R/, and used to be pulled in with a plain
# top-level source() in fct_regions.R / fct_data_loading.R.
#
# That cannot work in an installed package. Top-level code in R/ is evaluated
# ONCE, during `R CMD INSTALL`'s lazy-load step — inside `docker build`, where
# no /pipeline bind mount exists. Unguarded it threw "cannot open file", R
# reported "lazy loading failed" and deleted the half-installed package, and
# remotes::install_local() downgraded that to a warning, so `docker build`
# exited 0 while shipping an image with no clinego.app in it. Guarded with
# file.exists() it installs, but the guard is false at install time, so the
# functions are simply absent from the installed namespace — silently.
#
# .onLoad() fires at library() time instead, when the mount IS present, so the
# functions land in the namespace where package code can find them.
#
# ── The two run paths ─────────────────────────────────────────────────────────
# PACKAGE PATH (clinego.app::run_app) — needs the package installed in the
#   image. This file is what makes it work. It is the intended end state: once
#   the pipeline's scripts/ is COPYed into the image rather than bind-mounted,
#   the paths below resolve at build time too and the app becomes fully
#   self-contained. See the note at the app-install block in the Dockerfile.
#
# FILE PATH (scripts/clinego.app/dev.R) — sources R/*.R off the bind mount and
#   never calls library(clinego.app), so .onLoad() NEVER FIRES HERE. dev.R
#   therefore sources CLINEGO_SHARED_LIBS itself. If you add a library to the
#   vector below, dev.R picks it up automatically — but if you replace this
#   mechanism, dev.R must be updated in the same commit or dev mode breaks.
#
# Keep both paths working. Changing one and not the other is exactly how the
# package path rotted unnoticed for five weeks: nobody runs it day to day.

#' Pipeline libraries the app reuses. Absolute container paths.
#' @noRd
CLINEGO_SHARED_LIBS <- c(
    "/pipeline/scripts/R/lib/regions.R",         # cluster_snps_to_regions() etc.
    "/pipeline/scripts/R/utils/pval_threshold.R" # compute_pval_threshold() etc.
)

#' Source the shared libraries into an environment.
#'
#' Returns the paths that were NOT found, so callers can report a partial load
#' instead of failing later with "could not find function".
#' @noRd
load_shared_libs <- function(envir, paths = CLINEGO_SHARED_LIBS) {
    missing <- character(0)
    for (f in paths) {
        if (file.exists(f)) sys.source(f, envir = envir) else missing <- c(missing, f)
    }
    missing
}

.onLoad <- function(libname, pkgname) {
    # asNamespace() is writable here: R seals the namespace AFTER .onLoad returns.
    load_shared_libs(asNamespace(pkgname))
}

.onAttach <- function(libname, pkgname) {
    missing <- CLINEGO_SHARED_LIBS[!file.exists(CLINEGO_SHARED_LIBS)]
    if (length(missing)) {
        packageStartupMessage(
            "clinego.app: shared pipeline libraries not found:\n  ",
            paste(missing, collapse = "\n  "),
            "\nRegion computation and interactive thresholds will fail. ",
            "Mount the pipeline root at /pipeline (-v $PWD:/pipeline)."
        )
    }
}
