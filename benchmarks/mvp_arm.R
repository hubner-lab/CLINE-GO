# =============================================================================
# mvp_arm.R -- which slice of benchmarks/mvp_seeds.tsv a figure script reports.
#
# WHY THIS EXISTS. Every figure script used to open with `seeds[arm == "primary"]`
# and, in one case, `stopifnot(nrow(PRIM) == 90L)`. That was deliberate: it kept a
# second landscape from being silently absorbed into the legacy 90-replicate
# analysis. It also made a second arm unreportable without editing 18 scripts.
#
# This file replaces the literal with a selector whose DEFAULTS ARE THE OLD
# LITERALS, so a script sourcing it and calling mvp_prim() with no environment set
# behaves byte-identically to before. The isolation is preserved -- it now has to
# be opted out of explicitly, per run, instead of edited in.
#
#   MVP_ARM         arm to report          default "primary"
#   MVP_ADDED       cohort tag to report   default ""  (no filter)
#   MVP_N_EXPECT    asserted replicates    default 90
#   MVP_N_PER_ARCH  asserted per stratum   default 30
#
# MVP_ADDED IS NOT OPTIONAL FOR THE SS-CLINES ARM. mvp_seeds.tsv holds 120 rows
# each for ssclines b1, b2 and b3 under arm == "primary_ssclines" (360 total), so
# filtering on arm alone pools three demographies into one marginal. Any script
# reporting that arm must set both.
#
# Usage (from a script that already has data.table loaded):
#   source(file.path(ROOT, "benchmarks/mvp_arm.R"))
#   PRIM <- mvp_prim(seeds)
#   stopifnot(nrow(PRIM) == mvp_n_expect())
# =============================================================================

mvp_arm        <- function() Sys.getenv("MVP_ARM", "primary")
mvp_added      <- function() Sys.getenv("MVP_ADDED", "")
mvp_n_expect   <- function() as.integer(Sys.getenv("MVP_N_EXPECT", "90"))
mvp_n_per_arch <- function() as.integer(Sys.getenv("MVP_N_PER_ARCH", "30"))
# Trait-count restriction. Default "" = no filter, which is the legacy behaviour.
# Exists because the two arms are not trait-matched: the legacy 90 include 10
# 1-trait replicates and the SS-Clines blocks are 2-trait only by design. The
# 1-trait replicates REVERSE the headline (the >=2-of-3 panel loses to the causal
# loci on 10/10 of them), so a legacy-vs-block comparison that leaves them in is
# not the comparable number -- set MVP_N_TRAITS=2 to get it.
mvp_n_traits   <- function() Sys.getenv("MVP_N_TRAITS", "")

# Human-readable label for messages, figure subtitles and output provenance.
mvp_arm_label <- function() {
    a <- mvp_arm(); g <- mvp_added(); nt <- mvp_n_traits()
    if (nzchar(g)) {
        tags <- trimws(strsplit(g, ",", fixed = TRUE)[[1]])
        tags <- tags[nzchar(tags)]
        a <- if (length(tags) > 1L) sprintf("%s / %d blocks pooled (%s)",
                                            a, length(tags), paste(tags, collapse = " + "))
             else sprintf("%s / %s", a, g)
    }
    if (nzchar(nt)) a <- sprintf("%s / %s-trait", a, nt)
    a
}

# The one filter. Returns the manifest rows this run reports, and nothing else.
# `seeds` is benchmarks/mvp_seeds.tsv already read as a data.table.
mvp_prim <- function(seeds) {
    stopifnot(is.data.frame(seeds), "arm" %in% names(seeds))
    out <- seeds[seeds$arm == mvp_arm(), ]
    g <- mvp_added()
    if (nzchar(g)) {
        if (!"added" %in% names(seeds))
            stop("MVP_ADDED is set but mvp_seeds.tsv has no `added` column")
        # Comma-separated list, so pooling blocks is an EXPLICIT act. Reporting the
        # whole SS-Clines arm is MVP_ADDED="ssclines_nvar_mvar,ssclines_ncline_ns,..."
        # -- not "leave MVP_ADDED unset", which is indistinguishable from forgetting it.
        tags <- trimws(strsplit(g, ",", fixed = TRUE)[[1]])
        tags <- tags[nzchar(tags)]
        missing <- setdiff(tags, unique(seeds$added))
        if (length(missing))
            stop("MVP_ADDED names cohort tag(s) not in mvp_seeds.tsv: ",
                 paste(missing, collapse = ", "))
        out <- out[out$added %in% tags, ]
    }
    nt <- mvp_n_traits()
    if (nzchar(nt)) {
        if (!"n_traits" %in% names(seeds))
            stop("MVP_N_TRAITS is set but mvp_seeds.tsv has no `n_traits` column")
        out <- out[as.character(out$n_traits) == nt, ]
    }
    if (nrow(out) == 0L)
        stop(sprintf("no manifest rows for arm=%s added=%s n_traits=%s -- check MVP_ARM/MVP_ADDED/MVP_N_TRAITS",
                     mvp_arm(), if (nzchar(g)) g else "<unset>",
                     if (nzchar(nt)) nt else "<unset>"))
    # Guard the exact trap this file exists for: a multi-cohort arm reported
    # without a cohort tag silently pools demographies.
    if (!nzchar(g) && "added" %in% names(seeds) && length(unique(out$added)) > 1L)
        message(sprintf("!! arm %s spans %d cohort tags (%s) and MVP_ADDED is unset -- ",
                        mvp_arm(), length(unique(out$added)),
                        paste(sort(unique(out$added)), collapse = ", ")),
                "this POOLS them. Set MVP_ADDED to report one block.")
    message(sprintf("manifest slice: %s -- %d replicates", mvp_arm_label(), nrow(out)))
    out
}
