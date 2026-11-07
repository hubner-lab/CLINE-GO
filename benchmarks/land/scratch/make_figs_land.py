#!/usr/bin/env python3
"""Figure 2 + Figure S2 per landscape (plan nifty-prancing-grove Phase 8, scope chosen by the user
2026-11-06: "Figure 2 + S2 only"; the rest of S1-S8 / Table S1 only if the landscapes enter the paper).

Source: work/ld_window_lab/scratch_figs/sim_figures_ld.R (SS-Clines' LD-lab figure copy, which already
turns fixed-number checks into SOFT ASSERT warnings logged to SOFT_LOG). The copy keeps, verbatim:
  * the head (vocabulary, helpers, inputs, panel order, the detection plane A1 / D5 / D6),
  * the within-replicate rank section (R1, R2, Figure 2 = A1 + R1, reliability numbers),
  * the Figure S2 section (AUC-PR per method),
and drops the offset-accuracy, size-curve, S1, operating-point, S3-S8, Table S1, Data S2, captions
and zip sections. Each cut is located by its section title, asserted once. Patches (each once):
  * SIZE_DIR is not needed (its section is dropped); RDA_FIG_DIR must hold detection_table.tsv,
    written per landscape by detection_table_land.R (the head's harness tie reads it);
  * RNUM (S1's RDA brief numbers, SS-Clines only) is an empty table, so S2's tie to it is logged as
    a soft assert instead of compared;
  * ARCH_WRAP (defined in the dropped offset section, used by S2) is copied in front of S2;
  * a tail writes numbers.tsv and figures.tsv.
Usage: python3 make_figs_land.py
"""
import pathlib

SRC = pathlib.Path("/mnt/data/eugene/ADAPTOGENE/work/ld_window_lab/scratch_figs/sim_figures_ld.R")
OUT = pathlib.Path(__file__).resolve().parent / "fig2_s2_land.R"
BAR = "# =============================================================================\n"

t = SRC.read_text()

def section_start(title):
    key = BAR + "# " + title
    n = t.count(key)
    if n != 1:
        raise SystemExit(f"section {title!r}: found {n}")
    return t.index(key)

i_off  = section_start("Offset accuracy (C2, D5-offset)")
i_rank = section_start("Within-replicate rank (R1 main, R2 supp)")
i_size = section_start("Size curve (block b1)")
i_s2   = section_start("Figure S2: AUC-PR of the three single methods")
i_s3   = section_start("Figure S3: detection at five cutoffs")
assert i_off < i_rank < i_size < i_s2 < i_s3

head, rank, s2 = t[:i_off], t[i_rank:i_size], t[i_s2:i_s3]

def rep(s, old, new):
    n = s.count(old)
    if n != 1:
        raise SystemExit(f"expected 1 x {old[:60]!r}, found {n}")
    return s.replace(old, new)

head = rep(head, 'SIZ  <- need("SIZE_DIR")', 'SIZ  <- NA_character_   # [landscape figs] size section dropped')
# ARCH_WRAP is defined in the dropped offset-accuracy section; S2 uses it: copy its definition (asserted once)
AW = 'ARCH_WRAP <- c("oligogenic" = "oligogenic", "moderately polygenic" = "moderately\\npolygenic",\n               "highly polygenic" = "highly\\npolygenic")\n'
assert t[i_off:i_rank].count(AW) == 1, "ARCH_WRAP definition moved"
s2 = AW + s2
s2 = rep(s2, 'message("=== S2 AUC-PR per method")',
         'message("=== S2 AUC-PR per method")\n'
         'RNUM <- setnames(data.table(k_ = character(), value = character()), "k_", "key")  # [landscape figs] no S1 brief: tie logged as soft assert')

tail = r'''
# =============================================================================
# [landscape figs] numbers + figure list (the full report's captions / Data S2 / zip are not built)
# =============================================================================
NUMS <- rbindlist(NUM, use.names = TRUE, fill = TRUE)
setcolorder(NUMS, intersect(c("figure", "rule", "stratum", "engine", "metric", "value", "lo", "hi",
                              "expected", "n", "p_perm", "source"), names(NUMS)))
fwrite(NUMS, file.path(OUT, "numbers.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(FIGS), function(k) data.table(stem = k, role = FIGS[[k]]$role,
       w_mm = FIGS[[k]]$w, h_mm = FIGS[[k]]$h, source = FIGS[[k]]$source))), file.path(OUT, "figures.tsv"), sep = "\t")
message(sprintf("=== landscape figs done: %d figures, numbers.tsv %d rows -> %s", length(FIGS), nrow(NUMS), OUT))
'''
hdr = "# [landscape figs copy of work/ld_window_lab/scratch_figs/sim_figures_ld.R, made by work/mvp_land/scratch/make_figs_land.py]\n"
OUT.write_text(hdr + head + rank + s2 + tail)
print(f"wrote {OUT.name}: {len((head + rank + s2).splitlines())} lines kept of {len(t.splitlines())}")
