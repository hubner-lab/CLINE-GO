#!/usr/bin/env bash
# Install scripts/clinego.app FROM THE MOUNT into a throwaway library, put that
# library first on R's search path, then exec the command given as arguments.
#
# WHY THIS EXISTS
# ---------------
# scripts/clinego.app/tests/testthat.R is `library(clinego.app)` + test_check(),
# and Tier 5 (tests/testthat/test-equivalence-app-pipeline.R) reaches the app
# through `clinego.app:::`. Both resolve against /usr/local/lib/R/site-library,
# i.e. the copy baked into the image at Dockerfile:176 — NOT against the
# -v $PWD:/pipeline mount. So an edit under scripts/clinego.app/R/ was invisible
# to the gate until the next `docker build`: the suite went green while testing
# code that is no longer in the repository. That is the one failure mode a merge
# gate must not have.
#
# The throwaway library goes first on .libPaths() via R_LIBS, so no test file
# needs to change; the assertion below is what makes the redirection non-silent.
# NOT R_LIBS_USER: rocker's Renviron.site sets R_LIBS=site-library:library, and
# R orders .libPaths() as R_LIBS, R_LIBS_USER, R_LIBS_SITE — so an R_LIBS_USER
# lib lands AFTER site-library and the image copy still wins. Renviron.site uses
# ${R_LIBS-...}, so an R_LIBS set here replaces its value; site-library stays
# reachable through R_LIBS_SITE and .Library. Measured 2026-09-23: the
# R_LIBS_USER version tripped the assertion below on every run in the image.
#
# NOT remotes::install_local(): it short-circuits when the local package SHA1 is
# unchanged, leaving the target library EMPTY while find.package() quietly falls
# through to site-library — which is exactly the failure being fixed.
set -euo pipefail

APP_SRC="${CLINEGO_APP_SRC:-/pipeline/scripts/clinego.app}"
APP_LIB="${CLINEGO_APP_LIB:-/tmp/clinego_applib}"

mkdir -p "$APP_LIB"
log="$(mktemp)"
if ! R CMD INSTALL --no-docs --no-help --no-byte-compile \
        -l "$APP_LIB" "$APP_SRC" >"$log" 2>&1; then
    echo "with_app_from_mount: R CMD INSTALL of ${APP_SRC} FAILED" >&2
    cat "$log" >&2
    exit 1
fi
rm -f "$log"

export R_LIBS="$APP_LIB"

# Fail loudly if the redirection did not take. Without this the fallback is
# silent and indistinguishable from the bug.
Rscript -e '
    lib <- normalizePath(Sys.getenv("R_LIBS"))
    p   <- normalizePath(find.package("clinego.app"))
    if (!startsWith(p, lib))
        stop("clinego.app resolves to ", p, ", not the mounted source install in ",
             lib, " — the suite would be testing the image copy")
    cat("clinego.app under test:", p, "\n")
'

exec "$@"
