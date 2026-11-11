#!/usr/bin/env Rscript
# dump_layers.R -- the value check of the final redraw (lab, 2026-11-11). Runs one drawing script through sys.source()
# in its own environment (commandArgs() answered from DUMP_ARGS, tab-separated), then writes the built data of every
# layer of every ggplot left in that environment (patchwork members unpacked) as <DUMP_OUT>/<object>[_<member>]_L<i>.tsv,
# plus objects.tsv (object, member, layer, rows, geom). The data are what ggplot_build() hands the geoms: every plotted
# value, read from the plot, not from pixels. Compare two dumps with compare_dumps.py.
# Usage: DUMP_SCRIPT=<file> DUMP_OUT=<dir> [DUMP_ARGS=<a1\ta2...>] [DUMP_OBJECTS=<o1,o2>] Rscript dump_layers.R
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
SCRIPT <- Sys.getenv("DUMP_SCRIPT"); OUTD <- Sys.getenv("DUMP_OUT"); stopifnot(nzchar(SCRIPT), nzchar(OUTD))
ARGS <- Sys.getenv("DUMP_ARGS"); ARGS <- if (nzchar(ARGS)) strsplit(ARGS, "\t", fixed = TRUE)[[1]] else character()
ONLY <- Sys.getenv("DUMP_OBJECTS"); ONLY <- if (nzchar(ONLY)) strsplit(ONLY, ",", fixed = TRUE)[[1]] else NULL
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
e <- new.env(parent = globalenv())
e$commandArgs <- function(trailingOnly = FALSE) if (trailingOnly) ARGS else c("R", "--args", ARGS)
message("== dump_layers: ", SCRIPT, " ", paste(ARGS, collapse = " "))
sys.source(SCRIPT, envir = e)

members <- function(p, name) {                      # list of (name, ggplot) for a ggplot or a patchwork
    if (inherits(p, "patchwork")) {
        kids <- p$patches$plots
        out <- list()
        for (i in seq_along(kids)) out <- c(out, members(kids[[i]], paste0(name, "_m", i)))
        self <- p; self$patches <- NULL; class(self) <- setdiff(class(self), "patchwork")
        if (length(self$layers) > 0) out <- c(out, setNames(list(self), paste0(name, "_self")))
        return(out)
    }
    if (inherits(p, "ggplot")) return(setNames(list(p), name))
    list()
}
objs <- ls(e, all.names = TRUE)
if (!is.null(ONLY)) { stopifnot(all(ONLY %in% objs)); objs <- ONLY }
index <- list()
for (o in objs) {
    v <- get(o, envir = e)
    if (!inherits(v, c("ggplot", "patchwork"))) next
    for (m in names(ms <- members(v, o))) {
        b <- ggplot_build(ms[[m]])
        for (i in seq_along(b$data)) {
            d <- as.data.table(b$data[[i]])
            keep <- names(d)[vapply(d, function(x) is.atomic(x) && !is.list(x), logical(1))]
            d <- d[, ..keep]
            f <- sprintf("%s_L%d.tsv", m, i)
            fwrite(d, file.path(OUTD, f), sep = "\t")
            index[[length(index) + 1]] <- data.table(object = o, member = m, layer = i, rows = nrow(d),
                                                     geom = class(ms[[m]]$layers[[i]]$geom)[1], file = f)
        }
    }
}
fwrite(rbindlist(index), file.path(OUTD, "objects.tsv"), sep = "\t")
message("== dump_layers: ", length(index), " layers -> ", OUTD)
