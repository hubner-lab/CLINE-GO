#!/usr/bin/env Rscript
# Compare two results trees and report every difference.
#
# Usage:
#   Rscript /pipeline/scripts/compare_results.R <tree_A> <tree_B>
#       [--expect FILE] [--tol T] [--name-a NAME] [--name-b NAME]
#       [--out REPORT.tsv] [--exclude GLOB]...
#
# A is the reference (a baseline), B the tree under test: the same config run
# twice under two project names, or once before and once after a code change.
# Nothing in the pipeline compares two {PROJECT}_results/ directories; this does.
#
# WHAT IS COMPARED
#   * The file sets: every file below each tree except hidden path components
#     (.snakemake/) and --exclude globs. In B only -> new_file, in A only ->
#     removed_file.
#   * Tables (*.tsv, *.csv, also gzipped): byte-identical files are skipped, the
#     rest are read as LINES and split on the separator. Never through fread's type
#     and header guessing, which can turn "NA" into a missing value or a header-less
#     grid's first row into column names, and so hide or invent a difference. A
#     rectangular file is compared as a table: line 1 is its header unless every
#     field of it is a number (genetic_offset_map.tsv is a bare grid, whose columns
#     are then V1..Vn). Kinds: new_column, removed_column, changed_rows (the row
#     count differs, or the same rows come in another order), changed_values (one
#     row per column: how many rows, the largest difference, the first one). A
#     ragged file (bcftools stats output) is compared line by line.
#   * Every other file: md5. A text file that differs is re-read with the project
#     names normalised; if it still differs, and for binary files and renders (png,
#     svg, pdf, qs, qs2, tif, log), the kind is changed_bytes. changed_bytes NEVER
#     fails the run, because renders and serialized objects carry bytes that differ
#     between identical runs. It is there to find where a table difference starts
#     (an sNMF .Q or a .kinf that moved upstream of every TSV).
#
# Numbers: two fields that both parse as finite numbers are equal when
# |a - b| <= tol * max(|a|, |b|). The default --tol 0 is exact numeric equality
# (1e-05 equals 1.0e-05). Everything else is compared as text.
#
# Project names: each tree's name (default: the directory name minus "_results")
# becomes {PROJECT} in paths and in every field before anything is compared.
# "<name>_results" and "<name>_logs" go first, then the bare name as a whole token,
# with no letter, digit or "_" on either side: base1 never matches inside base10,
# and a project named after its VCF keeps the VCF's own prefixes (SIMDATA_K3imp).
#
# Expected differences (--expect FILE): a TSV with columns path, kind, column,
# finding; column and finding may be blank. path and column are globs whose *
# matches anything, "/" included, and a blank column matches any column. kind is
# one of new_file, removed_file, new_column, removed_column, changed_rows,
# changed_values, changed_bytes. A difference matched by a row is "expected". A row
# that matched nothing is reported as "expected_not_seen" and does not fail the run:
# a fix may leave SIMDATA untouched. Lines starting with # are comments.
#
# Report: one row per difference, columns status (unexpected / expected /
# expected_not_seen / info), kind, path, column, detail, finding. Written to --out,
# or to stdout without it; a summary goes to stderr. Nothing is written into the
# working directory.
#
# Exit status: 0 when nothing unexpected differs, 1 when something does, 2 on a
# usage or read error.

suppressPackageStartupMessages({
    library(data.table)
})

TABLE_RE    <- "\\.(tsv|csv)(\\.gz)?$"
VOLATILE_RE <- "\\.(png|svg|pdf|qs|qs2|tif|tiff|log)$"
KINDS       <- c("new_file", "removed_file", "new_column", "removed_column",
                 "changed_rows", "changed_values", "changed_bytes")
STATUSES    <- c("unexpected", "expected", "expected_not_seen", "info")
TEXT_MAX_BYTES <- 50e6   # a non-table file above this is compared by md5 alone

USAGE <- paste("usage: compare_results.R <tree_A> <tree_B> [--expect FILE] [--tol T]",
               "[--name-a NAME] [--name-b NAME] [--out FILE] [--exclude GLOB]...")

# ── arguments ─────────────────────────────────────────────────────────────────

parse_args <- function(args) {
    o <- list(expect = NULL, tol = "0", name_a = NULL, name_b = NULL, out = NULL,
              exclude = character(0))
    trees <- character(0)
    i <- 1L
    while (i <= length(args)) {
        a <- args[[i]]
        if (!startsWith(a, "--")) {
            trees <- c(trees, a)
            i <- i + 1L
            next
        }
        if (i == length(args)) stop("missing value for ", a, "\n", USAGE)
        v <- args[[i + 1L]]
        switch(a,
               "--expect"  = o$expect <- v,
               "--tol"     = o$tol <- v,
               "--name-a"  = o$name_a <- v,
               "--name-b"  = o$name_b <- v,
               "--out"     = o$out <- v,
               "--exclude" = o$exclude <- c(o$exclude, v),
               stop("unknown option ", a, "\n", USAGE))
        i <- i + 2L
    }
    if (length(trees) != 2L) stop("need exactly two trees, got ", length(trees), "\n", USAGE)
    for (t in trees) if (!dir.exists(t)) stop("not a directory: ", t)
    o$a <- normalizePath(trees[[1L]])
    o$b <- normalizePath(trees[[2L]])
    tol <- suppressWarnings(as.numeric(o$tol))
    if (length(tol) != 1L || !is.finite(tol) || tol < 0) stop("--tol must be a number >= 0, got ", o$tol)
    o$tol <- tol
    if (!is.null(o$expect) && !file.exists(o$expect)) stop("--expect file not found: ", o$expect)
    default_name <- function(dir) {
        b <- basename(dir)
        if (grepl("_results$", b)) sub("_results$", "", b) else ""
    }
    if (is.null(o$name_a)) o$name_a <- default_name(o$a)
    if (is.null(o$name_b)) o$name_b <- default_name(o$b)
    o
}

# ── normalisation ─────────────────────────────────────────────────────────────

normaliser <- function(name) {
    if (!nzchar(name)) return(identity)
    q <- paste0("\\Q", name, "\\E")
    re_dir  <- paste0("(?<![A-Za-z0-9_])", q, "_(results|logs)(?![A-Za-z0-9_])")
    re_bare <- paste0("(?<![A-Za-z0-9_])", q, "(?![A-Za-z0-9_])")
    function(x) gsub(re_bare, "{PROJECT}",
                     gsub(re_dir, "{PROJECT}_\\1", x, perl = TRUE), perl = TRUE)
}

# ── reading ───────────────────────────────────────────────────────────────────

# gzfile() reads plain files as they are, so one reader serves .tsv and .tsv.gz.
read_lines <- function(path) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    readLines(con, warn = FALSE)
}

is_binary <- function(path) {
    con <- gzfile(path, "rb")
    on.exit(close(con))
    any(readBin(con, "raw", n = 8000L) == as.raw(0L))
}

# strsplit() drops one trailing empty field; the appended separator gives it back,
# so "a\t" is two fields and "" is one.
split_fields <- function(lines, sep) strsplit(paste0(lines, sep), sep, fixed = TRUE)

num_or_na <- function(x) suppressWarnings(as.numeric(x))

looks_like_header <- function(fields) {
    numberish <- fields %in% c("", "NA", "NaN", "Inf", "-Inf") | !is.na(num_or_na(fields))
    any(!numberish)
}

# ── comparing ─────────────────────────────────────────────────────────────────

# Equal as text, or as two finite numbers within the relative tolerance.
fields_equal <- function(a, b, tol) {
    eq <- a == b
    i  <- which(!eq)
    if (length(i)) {
        x <- num_or_na(a[i])
        y <- num_or_na(b[i])
        num <- is.finite(x) & is.finite(y)
        eq[i[num]] <- abs(x[num] - y[num]) <= tol * pmax(abs(x[num]), abs(y[num]))
    }
    eq
}

clip <- function(s, n = 60L) {
    s <- gsub("\t", " | ", s, fixed = TRUE)
    ifelse(nchar(s) > n, paste0(substr(s, 1L, n - 3L), "..."), s)
}

# How many elements of x have no partner in y, counting duplicates.
only_in <- function(x, y) {
    tx <- table(x)
    ty <- table(y)
    m  <- as.integer(ty[match(names(tx), names(ty))])
    m[is.na(m)] <- 0L
    sum(pmax(as.integer(tx) - m, 0L))
}

diff_row <- function(kind, column, detail) data.table(kind = kind, column = column, detail = detail)

describe_values <- function(a, b, bad, line0) {
    k <- which(bad)
    x <- num_or_na(a[k])
    y <- num_or_na(b[k])
    num <- is.finite(x) & is.finite(y)
    s <- sprintf("%d of %d rows differ", length(k), length(a))
    if (any(num)) {
        d <- abs(x[num] - y[num])
        s <- sprintf("%s; max |diff| %.3g, max rel %.3g", s, max(d),
                     max(d / pmax(abs(x[num]), abs(y[num]))))
    }
    sprintf("%s; first at line %d: '%s' vs '%s'", s, line0 + k[1L] - 1L,
            clip(a[k[1L]]), clip(b[k[1L]]))
}

# A ragged file has no columns to name, so the unit is the line.
compare_lines <- function(la, lb, sa, sb, tol) {
    if (length(la) != length(lb)) {
        return(diff_row("changed_rows", "", sprintf(
            "%d lines in A, %d in B (ragged file, compared as lines); %d only in A, %d only in B",
            length(la), length(lb), only_in(la, lb), only_in(lb, la))))
    }
    bad <- vapply(seq_along(la), function(i) {
        if (la[[i]] == lb[[i]]) return(FALSE)
        if (length(sa[[i]]) != length(sb[[i]])) return(TRUE)
        !all(fields_equal(sa[[i]], sb[[i]], tol))
    }, logical(1))
    if (!any(bad)) return(NULL)
    if (identical(sort(la, method = "radix"), sort(lb, method = "radix"))) {
        return(diff_row("changed_rows", "",
                        sprintf("same %d lines in a different order (ragged file)", length(la))))
    }
    k <- which(bad)[1L]
    diff_row("changed_values", "", sprintf(
        "%d of %d lines differ (ragged file, compared as lines); first at line %d: '%s' vs '%s'",
        sum(bad), length(la), k, clip(la[[k]]), clip(lb[[k]])))
}

to_matrix <- function(rows, header) {
    matrix(as.character(unlist(rows, use.names = FALSE)), ncol = length(header),
           byrow = TRUE, dimnames = list(NULL, header))
}

row_keys <- function(m, cols) {
    if (!length(cols)) return(character(nrow(m)))
    do.call(paste, c(lapply(cols, function(col) m[, col]), sep = "\x1f"))
}

compare_table <- function(la, lb, sep, tol) {
    sa <- split_fields(la, sep)
    sb <- split_fields(lb, sep)
    wa <- unique(lengths(sa))
    wb <- unique(lengths(sb))
    if (length(wa) != 1L || length(wb) != 1L) return(compare_lines(la, lb, sa, sb, tol))

    header <- looks_like_header(sa[[1L]]) || looks_like_header(sb[[1L]])
    if (header) {
        ha <- make.unique(sa[[1L]]); hb <- make.unique(sb[[1L]])
        sa <- sa[-1L]; sb <- sb[-1L]
        line0 <- 2L
    } else {
        ha <- paste0("V", seq_len(wa)); hb <- paste0("V", seq_len(wb))
        line0 <- 1L
    }
    out <- c(lapply(setdiff(hb, ha), function(col) diff_row("new_column", col, "")),
             lapply(setdiff(ha, hb), function(col) diff_row("removed_column", col, "")))
    common <- intersect(ha, hb)
    ma <- to_matrix(sa, ha)
    mb <- to_matrix(sb, hb)
    ka <- row_keys(ma, common)
    kb <- row_keys(mb, common)

    if (nrow(ma) != nrow(mb)) {
        out <- c(out, list(diff_row("changed_rows", "", sprintf(
            "%d rows in A, %d in B; %d only in A, %d only in B (over the shared columns)",
            nrow(ma), nrow(mb), only_in(ka, kb), only_in(kb, ka)))))
        return(rbindlist(out))
    }
    bad <- setNames(lapply(common, function(col) !fields_equal(ma[, col], mb[, col], tol)), common)
    hit <- common[vapply(bad, any, logical(1))]
    if (length(hit)) {
        if (identical(sort(ka, method = "radix"), sort(kb, method = "radix"))) {
            out <- c(out, list(diff_row("changed_rows", "",
                                        sprintf("same %d rows in a different order", nrow(ma)))))
        } else {
            out <- c(out, lapply(hit, function(col) diff_row(
                "changed_values", col, describe_values(ma[, col], mb[, col], bad[[col]], line0))))
        }
    }
    rbindlist(out)
}

compare_other <- function(pa, pb, rel, norm_a, norm_b) {
    if (grepl(VOLATILE_RE, rel, ignore.case = TRUE)) {
        return(diff_row("changed_bytes", "", "render or serialized object; bytes differ"))
    }
    if (max(file.size(pa), file.size(pb)) > TEXT_MAX_BYTES || is_binary(pa) || is_binary(pb)) {
        return(diff_row("changed_bytes", "", "binary or large file; bytes differ"))
    }
    la <- norm_a(read_lines(pa))
    lb <- norm_b(read_lines(pb))
    if (identical(la, lb)) return(NULL)
    diff_row("changed_bytes", "", sprintf("text: %d lines in A, %d in B; %d only in A, %d only in B",
                                          length(la), length(lb), only_in(la, lb), only_in(lb, la)))
}

# ── expectations ──────────────────────────────────────────────────────────────

read_expect <- function(path) {
    empty <- data.table(path = character(0), kind = character(0),
                        column = character(0), finding = character(0))
    if (is.null(path)) return(empty)
    lines <- readLines(path, warn = FALSE)
    lines <- lines[!grepl("^\\s*(#|$)", lines)]
    if (!length(lines)) return(empty)
    hdr <- split_fields(lines[[1L]], "\t")[[1L]]
    if (!all(c("path", "kind") %in% hdr)) stop("--expect header must name path and kind: ", path)
    rows <- split_fields(lines[-1L], "\t")
    if (any(lengths(rows) > length(hdr))) stop("--expect row with more fields than its header: ", path)
    rows <- lapply(rows, function(r) c(r, rep("", length(hdr) - length(r))))
    ex <- as.data.table(to_matrix(rows, hdr))
    for (col in c("column", "finding")) if (!col %in% names(ex)) ex[, (col) := ""]
    bad <- setdiff(ex$kind, KINDS)
    if (length(bad)) stop("--expect: unknown kind(s) ", paste(bad, collapse = ", "),
                          "; allowed: ", paste(KINDS, collapse = ", "))
    ex[, .(path, kind, column, finding)]
}

apply_expect <- function(report, ex) {
    report[, `:=`(status = fifelse(kind == "changed_bytes", "info", "unexpected"), finding = "")]
    seen <- logical(nrow(ex))
    for (j in seq_len(nrow(ex))) {
        hit <- report$kind == ex$kind[[j]] &
               grepl(glob2rx(ex$path[[j]]), report$path) &
               (!nzchar(ex$column[[j]]) | grepl(glob2rx(ex$column[[j]]), report$column))
        seen[[j]] <- any(hit)
        first <- hit & !report$status %in% "expected"
        report[first, `:=`(status = "expected", finding = ex$finding[[j]])]
    }
    unseen <- ex[!seen]
    if (nrow(unseen)) {
        report <- rbind(report, unseen[, .(kind, path, column,
                                           detail = "listed in --expect, not observed",
                                           status = "expected_not_seen", finding)])
    }
    report
}

# ── main ──────────────────────────────────────────────────────────────────────

main <- function(args) {
    o <- parse_args(args)
    norm_a <- normaliser(o$name_a)
    norm_b <- normaliser(o$name_b)
    ex <- read_expect(o$expect)

    list_tree <- function(root, norm) {
        f <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
        f <- f[!grepl("(^|/)\\.", f)]
        key <- norm(f)
        keep <- !Reduce(`|`, lapply(o$exclude, function(g) grepl(glob2rx(g), key)), FALSE)
        if (anyDuplicated(key[keep])) stop("two files normalise to one path in ", root)
        setNames(f[keep], key[keep])
    }
    fa <- list_tree(o$a, norm_a)
    fb <- list_tree(o$b, norm_b)

    out <- c(lapply(setdiff(names(fb), names(fa)), function(k) diff_row("new_file", "", "")[, path := k]),
             lapply(setdiff(names(fa), names(fb)), function(k) diff_row("removed_file", "", "")[, path := k]))
    both <- intersect(names(fa), names(fb))
    pa <- file.path(o$a, fa[both])
    pb <- file.path(o$b, fb[both])
    differ <- which(unname(tools::md5sum(pa)) != unname(tools::md5sum(pb)))
    for (i in differ) {
        rel <- both[[i]]
        d <- if (grepl(TABLE_RE, rel, ignore.case = TRUE)) {
            sep <- if (grepl("\\.csv(\\.gz)?$", rel, ignore.case = TRUE)) "," else "\t"
            la <- norm_a(read_lines(pa[[i]]))
            lb <- norm_b(read_lines(pb[[i]]))
            if (identical(la, lb)) NULL else compare_table(la, lb, sep, o$tol)
        } else {
            compare_other(pa[[i]], pb[[i]], rel, norm_a, norm_b)
        }
        if (!is.null(d) && nrow(d)) out <- c(out, list(d[, path := rel]))
    }

    report <- rbindlist(c(list(data.table(kind = character(0), column = character(0),
                                          detail = character(0), path = character(0))), out),
                        use.names = TRUE)
    report <- apply_expect(report, ex)
    report[, ord := match(status, STATUSES)]
    setorderv(report, c("ord", "path", "kind", "column"))
    report <- report[, .(status, kind, path, column, detail, finding)]

    fwrite(report, if (is.null(o$out)) "" else o$out, sep = "\t", quote = FALSE)

    n <- function(s) sum(report$status == s)
    message(sprintf("INFO: A = %s (project '%s'), B = %s (project '%s'), tol %g",
                    o$a, o$name_a, o$b, o$name_b, o$tol))
    message(sprintf("INFO: %d files in A, %d in B, %d in both, %d byte-identical",
                    length(fa), length(fb), length(both), length(both) - length(differ)))
    message(sprintf("INFO: %d unexpected, %d expected, %d expected_not_seen, %d info (changed_bytes)",
                    n("unexpected"), n("expected"), n("expected_not_seen"), n("info")))
    u <- report[status == "unexpected"]
    for (r in seq_len(min(nrow(u), 50L))) {
        message(sprintf("UNEXPECTED: %s %s%s%s", u$kind[[r]], u$path[[r]],
                        if (nzchar(u$column[[r]])) paste0(" [", u$column[[r]], "]") else "",
                        if (nzchar(u$detail[[r]])) paste0(": ", u$detail[[r]]) else ""))
    }
    if (nrow(u) > 50L) message("UNEXPECTED: ... ", nrow(u) - 50L, " more in the report")
    if (nrow(u)) 1L else 0L
}

status <- tryCatch(main(commandArgs(trailingOnly = TRUE)),
                   error = function(e) {
                       message("ERROR: ", conditionMessage(e))
                       2L
                   })
quit(save = "no", status = status)
