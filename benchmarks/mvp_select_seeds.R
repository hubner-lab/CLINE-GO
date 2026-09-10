#!/usr/bin/env Rscript
# mvp_select_seeds.R -- build the MVP benchmark seed manifest.
#
# Inputs (both persisted under data/mvp/selection/):
#   summary_20220428_20220726.csv  -- BCO-DMO doi:10.26008/1912/bco-dmo.889769.1, one row per seed
#   0b-final_params-20220428.txt   -- ModelValidationProgram/MVP-NonClinalAF src/, seed -> design cell
#
# Output:
#   benchmarks/mvp_seeds.tsv  -- 14 rows: 12 primary + 2 degenerate controls
#
# Selection rule (docs/gea_simulation_benchmarks.md Sec 8.2):
#   primary  : meanFst >= 0.05 AND cor_PC1_temp^2 in [0.30, 0.60] AND cor_PC1_sal^2 <= 0.20
#              -> 55 seeds, all SS-Mtn. The 12 used here are the stratified subset in Sec 8.3.
#   controls : the INVERSE on the alignment axis (cor_PC1_temp^2 > 0.9, i.e. degenerate = the
#              Laruson failure mode) while holding demography and architecture constant against
#              two of the primaries, so LANDSCAPE is the only thing that differs.
#
# The 12 primary seeds are transcribed from the doc, then re-verified here against the raw
# summary table. A mismatch means the doc and the deposit disagree -- fail loudly, do not proceed.

suppressPackageStartupMessages({
  library(data.table)
})

ROOT      <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
SEL_DIR   <- file.path(ROOT, "data/mvp/selection")
OUT_TSV   <- file.path(ROOT, "benchmarks/mvp_seeds.tsv")

# Seeds proposed in docs/gea_simulation_benchmarks.md Sec 8.3 (stratified over architecture).
PRIMARY_SEEDS <- c(1231288, 1232218, 1231993, 1232398,   # oligogenic
                   1232578, 1232548, 1231858, 1233208,   # moderately polygenic
                   1232923, 1231798, 1232908, 1233133)   # highly polygenic

message("== reading selection inputs ==")
summ <- fread(file.path(SEL_DIR, "summary_20220428_20220726.csv"))
# The params file was written by R with row.names -- header carries one fewer field than the
# data rows, so read.table's row-name inference is what we want here (fread would misalign).
par  <- as.data.table(read.table(file.path(SEL_DIR, "0b-final_params-20220428.txt"),
                                 header = TRUE, stringsAsFactors = FALSE))

message(sprintf("  summary: %d seeds x %d cols", nrow(summ), ncol(summ)))
message(sprintf("  params : %d seeds x %d cols", nrow(par),  ncol(par)))

keep_par <- c("seed", "arch", "arch_level", "arch_level_sub", "demog_level", "demog_level_sub",
              "N_traits", "ispleiotropy", "MIG_x", "MIG_y")
stopifnot("K" %in% names(summ))
d <- merge(summ, par[, ..keep_par], by = "seed")
stopifnot(nrow(d) == nrow(summ))
message(sprintf("  joined : %d / %d seeds", nrow(d), nrow(summ)))

# The deposit ships Kendall correlations between structure and environment; the screening
# criteria are stated on R^2, so square them here rather than anywhere downstream.
d[, r2_pc1_temp := cor_PC1_temp^2]
d[, r2_pc1_sal  := cor_PC1_sal^2]

# ---------------------------------------------------------------- primary set ----
message("== primary set ==")
missing_seeds <- setdiff(PRIMARY_SEEDS, d$seed)
if (length(missing_seeds)) {
  stop("Seeds from docs Sec 8.3 not present in the deposit summary: ",
       paste(missing_seeds, collapse = ", "))
}
primary <- d[seed %in% PRIMARY_SEEDS]

# Re-verify the doc's own selection rule holds for every transcribed seed. This is the check
# that the doc was not mis-transcribed, and it costs nothing.
bad <- primary[!(meanFst >= 0.05 & r2_pc1_temp >= 0.30 & r2_pc1_temp <= 0.60 & r2_pc1_sal <= 0.20)]
if (nrow(bad)) {
  print(bad[, .(seed, meanFst, r2_pc1_temp, r2_pc1_sal)])
  stop("Primary seeds violate the Sec 8.2 selection rule -- doc and deposit disagree.")
}
if (length(unique(primary$demog_level)) != 1L || primary$demog_level[1] != "SS-Mtn") {
  stop("Expected all primary seeds in the SS-Mtn landscape; got: ",
       paste(unique(primary$demog_level), collapse = ", "))
}
message(sprintf("  %d seeds verified against the Sec 8.2 rule, all %s / %s",
                nrow(primary), primary$demog_level[1],
                paste(unique(primary$demog_level_sub), collapse = ",")))

# ---------------------------------------------------- previously selected panel ----
# APPEND-ONLY SELECTION. Every seed already in the manifest has been fetched, converted,
# swept and scored -- dropping one silently invalidates hours of compute and every summary
# built from it. The naive path does exactly that: pick_stratum() takes evenly spaced ranks
# out of the pool, so raising PER_STRATUM 10 -> 16 RE-PICKS rather than extends (measured:
# 4 highly-polygenic, 2 mod-polygenic and 2 oligogenic seeds of the group-1 panel fall out).
#
# So the seeds already on disk are the fixed point: they are carried through verbatim with
# their own arm/group/added tags, and the expansion only fills the gap between what exists
# and PER_STRATUM. Re-running with an unchanged PER_STRATUM is a no-op.
PREV <- if (file.exists(OUT_TSV)) fread(OUT_TSV) else NULL
if (!is.null(PREV)) {
  message(sprintf("== existing manifest: %d seeds (%d primary, %d control) ==",
                  nrow(PREV), sum(PREV$arm == "primary"), sum(PREV$arm != "primary")))
  gone <- setdiff(PREV$seed, d$seed)
  if (length(gone)) stop("Manifest seeds absent from the deposit summary: ",
                         paste(gone, collapse = ", "))
} else {
  message("== no existing manifest -- building from scratch ==")
}
# Primaries already committed to. Union with the transcribed set so a fresh build behaves
# exactly as before this change.
HELD_PRIMARY <- if (is.null(PREV)) {
  PRIMARY_SEEDS
} else {
  union(PRIMARY_SEEDS, PREV[arm == "primary"]$seed)
}
held <- d[seed %in% HELD_PRIMARY]

# ---------------------------------------------------------------- expansion -----
# Extend the SS-Mtn panel to PER_STRATUM seeds in each architecture level.
#
# Why the band widens from [0.30, 0.60] to [0.20, 0.60]: the 12 transcribed seeds were
# picked under the narrower rule and all 12 remain inside the wider one, so nothing is
# dropped -- the widening only enlarges the pool to draw from (55 -> 82 seeds).
#
# BAND_HI/BAND_LO are overridable because 30 per stratum does not fit [0.20, 0.60]: the
# in-band pool is 22 oligogenic / 27 mod-polygenic / 33 highly-polygenic. meanFst and
# R^2(PC1~sal) are non-binding inside SS-Mtn -- only the temp band moves the count.
#
# Why stratify on ARCHITECTURE and not on R^2: measured across the 14-seed panel, within
# the band R^2 explains little of the offset accuracy (rho -0.06..-0.24) while architecture
# swings Kendall's tau by 0.3-0.4. Architecture is the axis the claim rests on, so it is
# the axis that gets balanced.
#
# Selection is DETERMINISTIC -- no sample(). Within each stratum the in-band pool is sorted
# by R^2(PC1~temp) and the new seeds are taken at evenly spaced rank positions, so the
# masking band stays evenly covered and re-running reproduces the manifest exactly.
PER_STRATUM <- as.integer(Sys.getenv("MVP_PER_STRATUM", "10"))
BAND_LO <- as.numeric(Sys.getenv("MVP_BAND_LO", "0.20"))
BAND_HI <- as.numeric(Sys.getenv("MVP_BAND_HI", "0.60"))
# Cohort tags for the rows this run ADDS. Existing rows keep whatever they were tagged with.
COHORT     <- Sys.getenv("MVP_COHORT",     "group1")
COHORT_TAG <- Sys.getenv("MVP_COHORT_TAG", "group1_expansion")

# ---------------------------------------------------------------- complete block ----
# MVP_BLOCK = "landscape:demography", e.g. "SS-Clines:N-variable_m-variable".
#
# When set, the expansion is EVERY 2-trait seed of that design block -- 3 genic levels x
# 4 architecture sub-levels x 10 replicates = 120 -- with no threshold on any realized
# statistic. The band/stratum machinery below is bypassed entirely: a complete block has no
# selection rule inside it, which is the point. Confounding (the PC1~env correlation) then
# becomes a reported covariate spanning whatever the design produces, not a criterion, so
# there is nothing left for "you picked the easy replicates" to attach to.
#
# 1-trait is excluded on identifiability grounds, not convenience: with a single selected
# axis, structure is near-collinear with the environment by construction (those replicates
# sit at |tau|(PC1,temp) > 0.95), which is the same failure that retired the Laruson set.
#
# New rows get their OWN arm value, not "primary". Every existing figure script filters
# arm == "primary", so a new landscape stays invisible until a script opts in explicitly.
# The alternative silently mixes two landscapes into every marginal.
BLOCK     <- Sys.getenv("MVP_BLOCK", "")
BLOCK_ARM <- Sys.getenv("MVP_BLOCK_ARM", "primary_ssclines")

if (nzchar(BLOCK)) {
  bp <- strsplit(BLOCK, ":", fixed = TRUE)[[1]]
  if (length(bp) != 2L) stop("MVP_BLOCK must be 'landscape:demography', got: ", BLOCK)
  message(sprintf("== complete block %s / %s, cohort %s/%s, arm %s ==",
                  bp[1], bp[2], COHORT, COHORT_TAG, BLOCK_ARM))
  block_pool <- d[demog_level == bp[1] & demog_level_sub == bp[2] &
                  arch_level_sub != "1-trait"]
  if (!nrow(block_pool)) stop("MVP_BLOCK matched no seeds in the deposit: ", BLOCK)
  # The deposit is a full factorial at 10 reps per cell. Assert it rather than assume it --
  # a partial block is not a complete block and must not be reported as one.
  .cells <- block_pool[, .N, by = .(arch_level, arch_level_sub)]
  if (nrow(.cells) != 12L || any(.cells$N != 10L)) {
    print(.cells[order(arch_level, arch_level_sub)])
    stop("block is not the expected 3 genic x 4 sub-level x 10 replicate design")
  }
  expansion <- block_pool[!seed %in% PREV$seed]
  message(sprintf("  block holds %d seeds, %d already in the manifest, adding %d",
                  nrow(block_pool), nrow(block_pool) - nrow(expansion), nrow(expansion)))
} else {
  message(sprintf("== expansion to %d per stratum, band [%.2f, %.2f], cohort %s/%s ==",
                  PER_STRATUM, BAND_LO, BAND_HI, COHORT, COHORT_TAG))
  pool <- d[demog_level == "SS-Mtn" & meanFst >= 0.05 &
            r2_pc1_temp >= BAND_LO & r2_pc1_temp <= BAND_HI & r2_pc1_sal <= 0.20]
  message(sprintf("  in-band SS-Mtn pool [%.2f, %.2f]: %d seeds", BAND_LO, BAND_HI, nrow(pool)))
  print(pool[, .N, by = arch_level][order(arch_level)])

  if (!all(PRIMARY_SEEDS %in% pool$seed)) {
    stop("Widening the band dropped a transcribed primary seed -- refusing to proceed: ",
         paste(setdiff(PRIMARY_SEEDS, pool$seed), collapse = ", "))
  }
  # A held seed outside the current band is kept anyway (it is already computed), but it is
  # not silently ignored either -- narrowing the band after a cohort has run is a mistake.
  out_of_band <- setdiff(HELD_PRIMARY, pool$seed)
  if (length(out_of_band)) {
    message("  NOTE: ", length(out_of_band), " already-selected primary seed(s) fall outside ",
            "the current band and are kept regardless: ", paste(out_of_band, collapse = ", "))
  }

  pick_stratum <- function(lvl) {
    have  <- held[arch_level == lvl]$seed
    avail <- pool[arch_level == lvl & !seed %in% have][order(r2_pc1_temp)]
    need  <- PER_STRATUM - length(have)
    if (need <= 0L) {
      message(sprintf("  %-18s have %d, adding 0 (already at or above target)", lvl, length(have)))
      return(avail[0L])
    }
    if (nrow(avail) < need) {
      stop(sprintf("stratum %s: need %d more seeds but only %d available in band [%.2f, %.2f]",
                   lvl, need, nrow(avail), BAND_LO, BAND_HI))
    }
    idx <- unique(round(seq(1, nrow(avail), length.out = need)))
    message(sprintf("  %-18s have %d, adding %d of %d available",
                    lvl, length(have), need, nrow(avail)))
    avail[idx]
  }

  strata <- sort(unique(held$arch_level))
  expansion <- rbindlist(lapply(strata, pick_stratum))
  message(sprintf("  expansion: %d seeds", nrow(expansion)))
}

# ---------------------------------------------------------------- controls ------
# Degenerate arm: same demography sub-level and architecture as a primary seed, but the
# Est-Clines landscape, where R^2(PC1~temp) ~ 0.96. Expectation when these are run: structure
# correction destroys the signal, reproducing the Laruson collapse on purpose.
#
# Derived from the 12 TRANSCRIBED primaries only, never from the expanded panel: the Fst
# window below is a function of that set, so letting it grow with the panel would re-pick
# controls that have already been run. Once a manifest exists its controls are carried
# through unchanged and this block only re-verifies them.
message("== degenerate controls ==")
fst_lo <- min(primary$meanFst); fst_hi <- max(primary$meanFst)
cand <- d[demog_level == "Est-Clines" &
          demog_level_sub == unique(primary$demog_level_sub)[1] &
          r2_pc1_temp > 0.90 &
          meanFst >= fst_lo & meanFst <= fst_hi]
message(sprintf("  %d Est-Clines candidates in the primary set's Fst range [%.3f, %.3f]",
                nrow(cand), fst_lo, fst_hi))

# One control per architecture level that carries the most primaries, highest Fst within each.
target_arch <- c("moderately-polygenic", "highly-polygenic")
target_arch <- intersect(target_arch, unique(cand$arch_level))
if (length(target_arch) < 2L) {
  # Fall back to whatever architecture levels are available, still 2 controls.
  target_arch <- head(unique(cand$arch_level), 2L)
}
controls <- rbindlist(lapply(target_arch, function(a) {
  cand[arch_level == a][order(-meanFst)][1L]
}))
if (nrow(controls) != 2L || anyNA(controls$seed)) {
  print(cand[, .(seed, arch_level, meanFst, r2_pc1_temp, r2_pc1_sal)])
  stop("Could not pick 2 degenerate controls -- inspect the candidate table above.")
}
message(sprintf("  picked %s", paste(controls$seed, collapse = ", ")))

# Controls are append-only too. If the manifest already carries some, they win, and a
# re-derivation that disagrees is a loud failure rather than a silent swap.
if (!is.null(PREV) && any(PREV$arm != "primary")) {
  prev_ctrl <- PREV[arm != "primary"]$seed
  if (!setequal(prev_ctrl, controls$seed)) {
    stop("Re-derived controls (", paste(sort(controls$seed), collapse = ", "),
         ") disagree with the manifest's (", paste(sort(prev_ctrl), collapse = ", "),
         ") -- refusing to swap a control that has already been run.")
  }
  message("  match the existing manifest")
}

# ---------------------------------------------------------------- manifest ------
mk <- function(x, arm, added, group = COHORT) {
  data.table(
    seed            = as.integer(x$seed),
    arm             = arm,
    project         = paste0("MVP", x$seed),
    architecture    = x$arch,
    arch_level      = x$arch_level,
    demog_level     = x$demog_level,
    demog_level_sub = x$demog_level_sub,
    n_traits        = x$N_traits,
    ispleiotropy    = x$ispleiotropy,
    # "number of populations used in analyses" -- the authors' own per-seed K, shipped as a
    # column. It is the defensible source for sNMF.k_best here because sNMF cross-entropy has
    # no interior minimum on this landscape (100 demes on a continuous stepping-stone grid;
    # verified on seed 1231288, where it decreases monotonically from K=2 to K=10). Floored at
    # 2 for the pipeline, which needs >= 2 clusters for structure plots and pop_diff_test.
    K_authors       = as.integer(x$K),
    k_best          = pmax(2L, as.integer(x$K)),
    meanFst         = round(x$meanFst, 4),
    r2_pc1_temp     = round(x$r2_pc1_temp, 4),
    r2_pc1_sal      = round(x$r2_pc1_sal, 4),
    n_snps          = x$nSNPs,
    n_causal_pre    = x$num_causal_prefilter,
    n_causal_maf01  = x$num_causal_postfilter,
    n_causal_temp   = x$num_causal_temp,
    n_causal_sal    = x$num_causal_sal,
    final_LA        = round(x$final_LA, 4),
    # Cohort bookkeeping: `group` names the panel this seed belongs to, `added` records
    # whether it came from the original 14-seed transcription or the Group 1 expansion.
    #
    # APPENDED AT THE END ON PURPOSE. Several consumers read this manifest by COLUMN
    # POSITION, not by name -- mvp_run_sweep.sh's kbest_of() takes `$11`. Inserting these
    # after `arm` shifted every later column right by two, so k_best silently became
    # ispleiotropy and the harvest would have looked for {method}_pvalues_K1.tsv. New
    # columns go last; never in the middle.
    group           = group,
    added           = added
  )
}

if (is.null(PREV)) {
  # Fresh build: unchanged behaviour.
  manifest <- rbind(
    mk(primary,   "primary",            "initial"),
    mk(expansion, "primary",            "group1_expansion"),
    mk(controls,  "control_degenerate", "initial")
  )
} else {
  # Append: every existing row is REBUILT from the same deposit inputs (so a column added
  # later is filled in consistently) but keeps its own arm/group/added tags. Rebuilding and
  # then verifying is stronger than copying -- if the deposit or the derivation ever moved,
  # the check below fails loudly instead of leaving two eras of rows in one file.
  kept <- rbindlist(lapply(seq_len(nrow(PREV)), function(i) {
    r <- PREV[i]
    mk(d[seed == r$seed], r$arm, r$added, group = r$group)
  }))
  cmp <- intersect(names(PREV), names(kept))
  same <- vapply(cmp, function(cl) isTRUE(all.equal(kept[[cl]], PREV[[cl]],
                                                    check.attributes = FALSE)), logical(1))
  if (!all(same)) {
    stop("Rebuilt rows differ from the existing manifest in column(s): ",
         paste(cmp[!same], collapse = ", "),
         " -- the deposit or the derivation moved; refusing to rewrite.")
  }
  message(sprintf("== %d existing rows reproduced exactly (%d columns checked) ==",
                  nrow(kept), length(cmp)))
  manifest <- rbind(kept, mk(expansion,
                             if (nzchar(BLOCK)) BLOCK_ARM else "primary",
                             COHORT_TAG, group = COHORT))
}
setorder(manifest, arm, -meanFst)

# The panel must be balanced by construction, not by hope -- check it here rather than
# discovering an unbalanced stratum after 18 seeds have been fetched and run.
# Block mode does not go through pick_stratum at all -- its rows carry BLOCK_ARM, and the
# SS-Mtn primaries it leaves untouched sit at whatever PER_STRATUM built them (30), not at
# this run's PER_STRATUM default. Checking them here would fail on rows this run did not
# create. The block's own completeness is asserted where it is built (3 x 4 x 10).
.bal <- manifest[arm == "primary", .N, by = arch_level]
if (!nzchar(BLOCK) && any(.bal$N != PER_STRATUM)) {
  print(.bal)
  stop(sprintf("architecture strata are not balanced at %d each", PER_STRATUM))
}
if (nzchar(BLOCK)) {
  .blk <- manifest[arm == BLOCK_ARM, .N, by = .(demog_level, demog_level_sub)]
  message("== block arm rows ==")
  print(.blk)
}
if (anyDuplicated(manifest$seed)) stop("duplicate seed in manifest")

# Append-only, enforced: nothing that was in the manifest may leave it.
if (!is.null(PREV)) {
  lost <- setdiff(PREV$seed, manifest$seed)
  if (length(lost)) stop("Selection dropped already-computed seed(s): ",
                         paste(lost, collapse = ", "))
  message(sprintf("== append-only check passed: %d kept + %d new = %d ==",
                  nrow(PREV), nrow(manifest) - nrow(PREV), nrow(manifest)))
}

dir.create(dirname(OUT_TSV), showWarnings = FALSE, recursive = TRUE)
fwrite(manifest, OUT_TSV, sep = "\t")
message(sprintf("== wrote %s (%d seeds) ==", OUT_TSV, nrow(manifest)))
print(manifest[, .(seed, arm, arch_level, meanFst, r2_pc1_temp, r2_pc1_sal,
                   n_snps, n_causal_maf01, K_authors, k_best)])
