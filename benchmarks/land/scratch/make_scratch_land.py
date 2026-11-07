#!/usr/bin/env python3
"""Landscape copies of the SS-Clines scoring scripts (plan nifty-prancing-grove, Phase 8).

Committed benchmarks/ stay untouched. Each copy differs from its source only by the exact-string
replacements below, each asserted to apply exactly once (count given per replacement):
  * demography: a landscape's replicates carry cohort tags that do not name their demography
    (SS-Mtn's 80 legacy seeds are group1_expansion / group2_phase1 / group3_phase2 / initial), so
    right after mvp_prim() the copies set `added` from the manifest's demog_level_sub, in the
    SS-Clines tag spelling (ssclines_nvar_mvar ...), so every downstream label, level and the
    figure scripts' BLOCK_LAB keep working unchanged. LAND_TAG_CHECK asserts the 5 x 120 balance.
  * detection option B (mvp_p4_detection_b.R): d is read from MVP_LINKED_BP (the per-seed window the
    GEA sets used, work/mvp_land/results/seed_ld_bp_<L>.tsv) instead of round(decay_critical_bp):
    39 Est-Clines seeds use a declared 619 bp (user, 2026-11-03), for which the LD table holds NA
    (7) or the unreliable tiny value (32).
  * A3 reference (mvp_remeasure_600.R): A3_REF names the landscape's own AUC-PR table, written by the
    detection_600 copy from the same params dir; detection600/ is SS-Clines'. A missing row still
    fails A3 (anyNA on the merge), so it cannot pass vacuously.
  * detection_600 copy: PARAMS_DIR replaces the frozen params/ arm; the legacy gate is kept as is.
  * panel tables (mvp_panel_tables.R): eval_offset_lind.R:191 flags EVERY Est-Clines seed `control`
    (written when the only Est-Clines seeds were the 2 degenerate controls), and line 65 drops control
    rows, which empties Est-Clines. The copy slices by the seed list alone (`seed %in% SEEDS`); for
    SS-Mtn (control all FALSE) it is a no-op, checked byte for byte.
  * rank (mvp_panel_rank.R:191): its old-vs-new tie check assumed two-way ties (moves 0 or 0.5);
    mvp_rank_ties.R splits an m-way tie in 1/m steps, and Est-Clines has one three-way tie, so the
    copy accepts multiples of 1/m below 1. SS-Clines rank tables + SVGs unchanged (10/10).
Usage: python3 make_scratch_land.py   (writes next to itself)
"""
import pathlib

ROOT = pathlib.Path("/mnt/data/eugene/ADAPTOGENE")
OUT = pathlib.Path(__file__).resolve().parent

TAG_R = r'''
# ---- [landscape scratch copy] demography tag from demog_level_sub (SS-Clines tag spelling)
LAND_DEMOG_TAG <- c("N-variable_m-variable"             = "ssclines_nvar_mvar",
                    "N-cline-N-to-S_m-constant"         = "ssclines_ncline_ns",
                    "N-equal_m-constant"                = "ssclines_nequal_mconst",
                    "N-cline-center-to-edge_m-constant" = "ssclines_ncline_ctredge",
                    "N-equal_m_breaks"                  = "ssclines_nequal_mbreaks")
stopifnot(all(PRIM$demog_level_sub %in% names(LAND_DEMOG_TAG)))
PRIM[, added := unname(LAND_DEMOG_TAG[demog_level_sub])]
stopifnot(uniqueN(PRIM$added) == 5L, all(PRIM[, .N, by = added]$N == nrow(PRIM) / 5L))
message("[landscape scratch] added <- demog_level_sub: ", paste(PRIM[, .N, by = added][, paste(added, N)], collapse = ", "))
'''

def patch(src, dst, reps):
    t = (ROOT / src).read_text()
    for old, new, want in reps:
        n = t.count(old)
        if n != want:
            raise SystemExit(f"{src}: expected {want} x {old[:70]!r}, found {n}")
        t = t.replace(old, new)
    hdr = f"# [landscape scratch copy of {src}, made by work/mvp_land/scratch/make_scratch_land.py]\n"
    (OUT / dst).write_text(t.replace("\n", "\n" + hdr, 1) if t.startswith("#!") else hdr + t)
    print(f"wrote {dst} ({len(reps)} patches)")

def tag_after(line):
    return (line, line + TAG_R, 1)

patch("benchmarks/mvp_p4_detection_b.R", "detection_b_land.R", [
    ("    d <- round(lt$decay_critical_bp)\n",
     "    d <- LAND_BP$d_bp[LAND_BP$seed == s]   # [landscape scratch] the window the sets used\n"
     "    stopifnot(length(d) == 1L, is.finite(d), d > 0)\n", 1),
    ("SETS <- c(",
     "LAND_BP <- fread(Sys.getenv(\"MVP_LINKED_BP\"), colClasses = c(seed = \"character\"))  # [landscape scratch]\n"
     "stopifnot(nrow(LAND_BP) > 0, all(seeds %in% LAND_BP$seed))\n"
     "SETS <- c(", 1),
])
patch("benchmarks/mvp_detection_600.R", "detection_600_land.R", [
    tag_after("PRIM <- mvp_prim(man)\n"),
    ('score_seed(s, file.path(EVAL, "params", paste0("MVP", s), CELL),',
     'score_seed(s, file.path(Sys.getenv("PARAMS_DIR"), paste0("MVP", s), CELL),', 1),
])
patch("benchmarks/mvp_remeasure_600.R", "remeasure_600_land.R", [
    tag_after("PRIM <- mvp_prim(man)\n"),
    ('ref <- fread(file.path(EVAL, "detection600", "aucpr_per_seed.tsv"), colClasses = c(seed = "character"))',
     'ref <- fread(Sys.getenv("A3_REF"), colClasses = c(seed = "character"))  # [landscape scratch]', 1),
    ('AUC-PR does not reproduce detection600', 'AUC-PR does not reproduce A3_REF', 1),
    ('AUC-PR reproduces detection600 on', 'AUC-PR reproduces A3_REF on', 1),
])
patch("benchmarks/mvp_j16_gallery.R", "j16_gallery_land.R", [
    tag_after("PRIM <- mvp_prim(man)\n"),
])
patch("benchmarks/mvp_panel_rank.R", "panel_rank_land.R", [
    tag_after("PRIM  <- mvp_prim(seeds)\n"),
    # the old-vs-new count check assumed two-way ties (moves of 0 or 0.5); mvp_rank_ties.R splits an
    # m-way tie in 1/m steps (its header: "a three-way tie ... is still found"). Est-Clines has one
    # three-way tie (SS-Clines none), so the check allows multiples of 1/m below 1 (m = the panel's tie size).
    ('CMP <- merge(CMP, RX[, .(tied = any(tie_n > 1L)), by = .(seed, engine)], by = c("seed", "engine"))\n',
     'CMP <- merge(CMP, RX[, .(tied = any(tie_n > 1L)), by = .(seed, engine)], by = c("seed", "engine"))\n'
     'CMP <- merge(CMP, RX[, .(m_tie = max(tie_n)), by = .(seed, engine, panel)], by = c("seed", "engine", "panel"))  # [landscape scratch]\n', 1),
    ('    if (!all(abs(d[CMP$tied]) %in% c(0, 0.5))) stop(v, ": a tied panel moved by other than 0 or 0.5")\n',
     '    dm <- abs(d[CMP$tied]) * CMP$m_tie[CMP$tied]   # [landscape scratch] moves in 1/m steps\n'
     '    if (any(abs(dm - round(dm)) > 1e-9 | abs(d[CMP$tied]) >= 1)) stop(v, ": a tied panel moved by other than k / m (m = its tie size)")\n', 1),
])
patch("benchmarks/mvp_panel_tables.R", "panel_tables_land.R", [
    ('GP[garden_type == "landscape" & control == FALSE & seed %in% SEEDS & !is.na(tau),',
     'GP[garden_type == "landscape" & seed %in% SEEDS & !is.na(tau),   # [landscape scratch] no control filter', 1),
])
