# dev.R — Development runner for CLINE-GO Shiny app
# Sources R files from mounted volume instead of installed package.
# Avoids Docker rebuilds: edit files and shiny.autoreload restarts the app.
#
# Usage:
#   docker run --user $(id -u):$(id -g) --rm -p 3838:3838 -v $PWD:/pipeline cline-go:latest \
#     Rscript /pipeline/scripts/clinego.app/dev.R

options(shiny.autoreload = TRUE)

# Auto-install any packages from DESCRIPTION that are not yet in the image.
# Uses a writable user library inside the pipeline bind-mount so no rebuild needed.
.dev_user_lib <- "/pipeline/.R_libs_dev"
dir.create(.dev_user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(.dev_user_lib, .libPaths()))

.dev_pkgs <- c("processx", "shinyjs", "shinyFiles")
for (.p in .dev_pkgs) {
    if (!requireNamespace(.p, quietly = TRUE)) {
        message("dev.R: installing missing package: ", .p)
        install.packages(.p,
                         lib   = .dev_user_lib,
                         repos = "https://cloud.r-project.org",
                         quiet = TRUE)
    }
}
rm(.dev_user_lib, .dev_pkgs, .p)

# Load all dependencies (from DESCRIPTION Imports, minus base R packages)
library(shiny)
library(bslib)
library(bsicons)
library(golem)
library(htmltools)
library(data.table)
library(plotly)
library(DT)
library(yaml)
library(jsonlite)
library(base64enc)
library(cachem)
library(config)
library(processx)
library(shinyjs)
library(shinyFiles)

# Source all app R files from mounted volume
app_r_dir <- "/pipeline/scripts/clinego.app/R"
for (f in sort(list.files(app_r_dir, pattern = "\\.R$", full.names = TRUE))) {
    source(f, local = FALSE)
}

# Shared pipeline libraries (lib/regions.R, utils/pval_threshold.R).
#
# On the package path these are loaded by .onLoad() in zzz.R. This path never
# calls library(clinego.app) -- it source()s R/*.R directly -- so .onLoad() does
# NOT fire here and the libraries must be loaded explicitly. Without this,
# compute_all_regions() dies with "could not find function
# cluster_snps_to_regions", and interactive thresholds silently return NA for
# every method (the compute_pval_threshold() call sits inside a tryCatch that
# maps errors to status = "error"), i.e. no significant SNPs anywhere.
#
# CLINEGO_SHARED_LIBS and load_shared_libs() are defined in zzz.R, already
# sourced by the loop above -- add a library there and this picks it up.
.missing_libs <- load_shared_libs(globalenv())
if (length(.missing_libs)) {
    stop("dev.R: shared pipeline libraries not found:\n  ",
         paste(.missing_libs, collapse = "\n  "),
         "\nIs the pipeline root mounted at /pipeline?")
}
rm(.missing_libs)

# Serve static assets (CSS/JS/SCSS from inst/app/www/)
shiny::addResourcePath("www", "/pipeline/scripts/clinego.app/inst/app/www")

# Run app
shiny::shinyApp(
    ui     = app_ui,
    server = app_server,
    options = list(host = "0.0.0.0", port = 3838)
)
