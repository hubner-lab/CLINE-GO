#!/usr/bin/env Rscript
# =============================================================================
# mvp_geo_offset_equiv_compare.R -- compare two geometric_offset output trees file by file.
#
# Why: the 2026-10-01 single-fit scripts/geometric_offset.R must write the same offsets as the
# per-scenario LEA::genetic.gap() version. This pairs every file under
#   <root>/Maladaptation/{tables,plots}/geometric_offset/<panel>_nospatial/
# (+ the offset rasters under <root>/_intermediate/geometric_offset/<panel>_nospatial/)
# of tree A with the same relative path in tree B and reports, per file:
#   .tsv  every numeric column / cell: max |A - B|, max relative diff, NA pattern equal,
#         non-numeric cells identical, byte-identical
#   .tif  cell values (terra): max |A - B|, max relative diff, NA pattern equal
#   other byte-identical only (e.g. the importance PNG)
# A file present in only one tree is reported, never skipped. Nothing is written but --out.
#
# Usage: Rscript mvp_geo_offset_equiv_compare.R --a=DIR --b=DIR --panel=truth --label=X --out=FILE
#   DIR = the directory that CONTAINS Maladaptation/ and _intermediate/ (a results tree or a replay dir)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(terra) })
source(file.path(Sys.getenv("PIPELINE_ROOT", "/pipeline"), "benchmarks/lib_detection.R"))   # parse_kv_args
a <- parse_kv_args(commandArgs(trailingOnly = TRUE))
for (k in c("a", "b", "panel", "label", "out")) if (is.null(a[[k]])) stop("missing --", k)

rel_files <- function(root) {
    dirs <- c(file.path(root, "Maladaptation", c("tables", "plots"), "geometric_offset",
                        paste0(a$panel, "_nospatial")),
              file.path(root, "_intermediate", "geometric_offset", paste0(a$panel, "_nospatial")))   # rasters
    f <- unlist(lapply(dirs[dir.exists(dirs)], list.files, recursive = TRUE, full.names = TRUE))
    sub(paste0("^", root, "/?"), "", f)
}
fa <- rel_files(a$a); fb <- rel_files(a$b)
if (!length(fa)) stop("no geometric_offset outputs for panel ", a$panel, " under ", a$a)
both <- intersect(fa, fb)

num_cmp <- function(x, y) {
    ok <- !is.na(x) & !is.na(y)
    d  <- abs(x[ok] - y[ok])
    s  <- pmax(abs(x[ok]), abs(y[ok]))
    list(n = length(x), max_abs = if (length(d)) max(d) else 0,
         max_rel = if (length(d)) max(ifelse(s > 0, d / s, 0)) else 0,
         na_equal = identical(is.na(x), is.na(y)))
}
md5 <- function(f) unname(tools::md5sum(f))

rows <- lapply(both, function(r) {
    pa <- file.path(a$a, r); pb <- file.path(a$b, r)
    bytes <- md5(pa) == md5(pb)
    ext <- tolower(tools::file_ext(r))
    res <- list(n = NA_integer_, max_abs = NA_real_, max_rel = NA_real_, na_equal = NA, text_equal = NA)
    if (ext == "tsv") {
        A <- fread(pa, header = FALSE, colClasses = "character")
        B <- fread(pb, header = FALSE, colClasses = "character")
        if (!identical(dim(A), dim(B))) {
            res$text_equal <- FALSE
        } else {
            va <- unlist(A, use.names = FALSE); vb <- unlist(B, use.names = FALSE)
            na <- suppressWarnings(as.numeric(va)); nb <- suppressWarnings(as.numeric(vb))
            isnum <- !is.na(na) | !is.na(nb)
            res$text_equal <- identical(va[!isnum], vb[!isnum])
            res[c("n", "max_abs", "max_rel", "na_equal")] <- num_cmp(na[isnum], nb[isnum])
        }
    } else if (ext == "tif") {
        res[c("n", "max_abs", "max_rel", "na_equal")] <- num_cmp(values(rast(pa))[, 1], values(rast(pb))[, 1])
    }
    data.table(label = a$label, file = r, type = ext, bytes_identical = bytes, n = res$n,
               max_abs = res$max_abs, max_rel = res$max_rel, na_equal = res$na_equal,
               text_equal = res$text_equal)
})
OUT <- rbindlist(rows)
only <- c(setdiff(fa, fb), setdiff(fb, fa))
if (length(only)) OUT <- rbind(OUT, data.table(label = a$label, file = only, type = "ONLY_IN_ONE_TREE",
                                               bytes_identical = FALSE), fill = TRUE)
fwrite(OUT, a$out, sep = "\t")
cat(sprintf("%s: %d files compared, %d in one tree only; byte-identical %d; max_abs %.3g; max_rel %.3g; NA pattern equal %s; text equal %s\n",
            a$label, length(both), length(only), sum(OUT$bytes_identical, na.rm = TRUE),
            max(OUT$max_abs, na.rm = TRUE), max(OUT$max_rel, na.rm = TRUE),
            all(OUT$na_equal, na.rm = TRUE), all(OUT$text_equal, na.rm = TRUE)))
print(OUT[, .(files = .N, bytes_identical = sum(bytes_identical), max_abs = max(max_abs, na.rm = TRUE),
              max_rel = max(max_rel, na.rm = TRUE)), by = type])
