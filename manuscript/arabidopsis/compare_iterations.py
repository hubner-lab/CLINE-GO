#!/usr/bin/env python3
"""Iteration 1 (Arabidopsis1135) against iteration 2 (Arabidopsis1135maf01), plan shiny-seeking-scott Phase 13.

Reads only tables that earlier phases wrote; nothing is refitted. Writes into ANALYSIS/Arabidopsis1135/compare/:

  ld_single_draws.tsv    the genome-wide critical distance each of iteration 2's 100 background draws would give on
                         its own (the rule before 6e7efbb took one draw, seed 1): the draw's critical r2 crossed with
                         the project's genome-wide curve by ld_first_crossing() of scripts/R/lib/ld_decay.R. Checked
                         first: the median and the 2.5 / 97.5 % r2 must give the pipeline's own decay_critical_bp,
                         _lo and _hi.
  level3_draws.tsv       per project and offset method, the union set's mean rho and each of its 20 MAF-matched draws'
                         (mean over replicates, as offset_eval.R), with the union's rank among the 21. Checked first
                         against summary.tsv (rho_mean, random_mean, random_q95).
  level3_holm.tsv        Holm over the three union level-3 tests of each project: a sensitivity, since offset_eval.R
                         has no multiple-testing step at level 3.
  iteration_compare.tsv  one row per metric: iteration 1, iteration 2, the file it comes from.

Usage: python3 -I manuscript/arabidopsis/compare_iterations.py [repo root]
"""
import csv
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "ANALYSIS" / "Arabidopsis1135"
OUT = ANALYSIS / "compare"
P1, P2 = "Arabidopsis1135", "Arabidopsis1135maf01"
K = {P1: 6, P2: 11}
ARM = {P1: "bonf_79kb", P2: "bonf_56kb"}
SETS = ["emmax", "lfmm", "rda", "union", "two_of_three", "three_of_three"]  # the first project's order, everywhere
METHODS = ["gradient_forest", "geometric_offset", "rda_offset"]
CLASSIC = ["FRI", "FLC", "FT", "SVP", "DOG1"]


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {path.relative_to(ROOT)} ({len(rows)} rows)")


def check(ok, what):
    if not ok:
        sys.exit(f"check failed: {what}")


def quantile7(xs, p):
    """R's default quantile (type 7)."""
    s = sorted(xs)
    h = (len(s) - 1) * p
    lo = math.floor(h)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (h - lo) * (s[hi] - s[lo])


def first_crossing(curve, target):
    """ld_first_crossing() of scripts/R/lib/ld_decay.R: first distance where r2_mono reaches target, log10 interpolation."""
    y = [r2 for _, r2 in curve]
    if y[0] <= target:
        return float("nan")
    for i, v in enumerate(y):
        if v <= target:
            lx0, lx1 = math.log10(curve[i - 1][0]), math.log10(curve[i][0])
            f = (y[i - 1] - target) / (y[i - 1] - y[i])
            return 10 ** (lx0 + f * (lx1 - lx0))
    return float("nan")


def holm(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, running = [0.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(ps) - rank) * ps[i]))
        adj[i] = running
    return adj


def results(p):
    return ROOT / f"{p}_results"


def ld_genome_wide(p):
    rows = read(results(p) / "Structure/tables/ld_decay_half_distances.tsv")
    gw = [r for r in rows if r["group"] == "All" and r["scope"] == "genome_wide"]
    check(len(gw) == 1, f"{p}: one genome-wide All row")
    return gw[0]


def ld_single_draws():
    tables = results(P2) / "Structure/tables"
    curve = [(float(r["dist_bp"]), float(r["r2_mono"])) for r in read(tables / "ld_decay_curve.tsv")
             if r["group"] == "All" and r["scope"] == "genome_wide"]
    gw = ld_genome_wide(P2)
    # The table prints r2 to 5 significant digits, so a crossing can move by a few bp; 0.1 % is far below any claim.
    for r2_col, bp_col in (("background_r2", "decay_critical_bp"), ("background_r2_hi", "decay_critical_bp_lo"),
                           ("background_r2_lo", "decay_critical_bp_hi")):
        got, want = first_crossing(curve, float(gw[r2_col])), float(gw[bp_col])
        check(abs(got - want) / want < 1e-3, f"crossing of {r2_col} gives {got:.0f}, the pipeline {want:.0f}")
    draws = [r for r in read(tables / "ld_decay_background_draws.tsv") if r["group"] == "All"]
    check(len(draws) == int(gw["n_background_draws"]), "one row per background draw")
    rows = []
    for r in sorted(draws, key=lambda r: int(r["draw"])):
        bp = first_crossing(curve, float(r["critical_r2"]))
        rows.append([int(r["draw"]), int(r["seed"]), float(r["critical_r2"]), round(bp)])
    order = sorted(rows, key=lambda r: r[3])
    rank = {r[0]: i + 1 for i, r in enumerate(order)}
    out = [r + [rank[r[0]]] for r in rows]
    write(OUT / "ld_single_draws.tsv", ["draw", "seed", "critical_r2", "critical_bp", "rank_shortest_first"], out)
    bps = [r[3] for r in rows]
    return {"single_min": min(bps), "single_max": max(bps), "draw1_bp": rows[0][3], "draw1_rank": rank[1],
            "draw1_r2": rows[0][2], "n": len(rows)}


def level3():
    draws_out, holm_out, summary = [], [], {}
    for p in (P1, P2):
        offset = ANALYSIS / "offset" / p
        summ = {(r["method"], r["set"]): r for r in read(offset / "summary.tsv")}
        sums, counts = {}, {}
        for r in read(offset / "replicate_rho.tsv"):
            if r["set"] == "gea.union" or r["set"].startswith("random.gea.union."):
                key = (r["method"], r["set"])
                sums[key] = sums.get(key, 0.0) + float(r["rho_offset"])
                counts[key] = counts.get(key, 0) + 1
        means = {k: sums[k] / counts[k] for k in sums}
        ps = []
        for m in METHODS:
            s = summ[(m, "gea.union")]
            union = means[(m, "gea.union")]
            draw = sorted((int(k[1].rsplit(".", 1)[1]), v) for k, v in means.items()
                          if k[0] == m and k[1].startswith("random."))
            vals = [v for _, v in draw]
            check(len(vals) == int(s["random_draws"]) == 20, f"{p} {m}: 20 draws")
            check(abs(union - float(s["rho_mean"])) < 1e-12, f"{p} {m}: union rho_mean")
            check(abs(statistics.fmean(vals) - float(s["random_mean"])) < 1e-12, f"{p} {m}: random_mean")
            check(abs(quantile7(vals, 0.95) - float(s["random_q95"])) < 1e-12, f"{p} {m}: random_q95")
            rank = 1 + sum(v > union for v in vals)
            draws_out.append([p, m, "gea.union", "", f"{union:.6f}", rank])
            draws_out += [[p, m, f"random.gea.union.{d}", d, f"{v:.6f}", ""] for d, v in draw]
            ps.append(float(s["random_p"]))
            summary[(p, m)] = {"union": union, "draw_min": min(vals), "draw_max": max(vals),
                               "draw_mean": statistics.fmean(vals), "rank": rank, "p": float(s["random_p"]),
                               "gardens_above": int(s["gardens_above_random"]), "level_3": s["level_3"],
                               "climate_mean": float(s["climate_mean"])}
        for m, p_raw, p_adj in zip(METHODS, ps, holm(ps)):
            holm_out.append([p, m, f"{p_raw:.5f}", f"{p_adj:.5f}", str(p_adj < 0.05).upper()])
            summary[(p, m)]["holm"] = p_adj
    write(OUT / "level3_draws.tsv", ["project", "method", "run", "draw", "rho_mean", "union_rank_of_21"], draws_out)
    write(OUT / "level3_holm.tsv", ["project", "method", "random_p", "random_p_holm3", "met_after_holm"], holm_out)
    return summary


def merged_loci(rows):
    """Distinct stretches over the union's GEA x GWAS rows of both traits (overlapping rows are one locus)."""
    iv = sorted((r["chr"], int(r["start"]), int(r["end"])) for r in rows)
    loci = []
    for c, s, e in iv:
        if loci and loci[-1][0] == c and s <= loci[-1][2]:
            loci[-1][2] = max(loci[-1][2], e)
        else:
            loci.append([c, s, e])
    return loci


def metrics(p, ld, l3):
    res, reg = results(p), ANALYSIS / "regions" / p
    m = {}
    filt = {r["stage"]: r for r in read(res / "Processing/tables/filtering_summary.tsv")}
    src_filt = f"{p}_results/Processing/tables/filtering_summary.tsv"
    m["samples kept"] = (filt["After relatedness filter"]["n_samples"], src_filt)
    m["SNPs after MAF and SNP missingness"] = (filt["After MAF + SNP missingness filter"]["n_snps"], src_filt)
    m["SNPs after LD pruning"] = (filt["After LD pruning"]["n_snps"], src_filt)
    thr = read(res / f"GEA/tables/methods/EMMAX/EMMAX_thresholds_K{K[p]}.tsv")[0]
    src_thr = f"{p}_results/GEA/tables/methods/EMMAX/EMMAX_thresholds_K{K[p]}.tsv"
    m["SNPs tested in GEA and GWAS"] = (thr["n_tested"], src_thr)
    m["Bonferroni threshold per predictor"] = (f"{float(thr['threshold']):.3g}", src_thr)
    gw = ld_genome_wide(p)
    src_ld = f"{p}_results/Structure/tables/ld_decay_half_distances.tsv"
    draws = gw.get("n_background_draws") or "1 (before the key)"
    m["background draws for the critical r2"] = (draws, src_ld)
    lo, hi = gw.get("decay_critical_bp_lo") or "", gw.get("decay_critical_bp_hi") or ""
    rng = f" ({lo}-{hi})" if lo else ""
    m["critical distance, genome-wide (bp)"] = (gw["decay_critical_bp"] + rng, src_ld)
    m["half-decay distance (bp)"] = (gw["decay_half_bp"], src_ld)
    m["short-range r2"] = (gw["short_r2"], src_ld)
    sets = {r["set"]: r for r in read(res / "Maladaptation/tables/snp_sets_summary.tsv")}
    src_sets = f"{p}_results/Maladaptation/tables/snp_sets_summary.tsv"
    for s in SETS:
        m[f"set {s} (SNPs)"] = (sets[f"gea.{s}"]["n_snps"], src_sets)
    chk = {r["scan"]: r for r in read(reg / "check/sets.tsv") if r["arm"] == ARM[p]}
    src_chk = f"ANALYSIS/Arabidopsis1135/regions/{p}/check/sets.tsv"
    u = chk["gea.union"]
    m["union regions"] = (u["n_regions"], src_chk)
    m["union regions span (Mb, share of genome)"] = (f"{float(u['span_mb']):.1f} ({100 * float(u['genome_share']):.1f} %)",
                                                     src_chk)
    m["union regions genes"] = (u["n_genes"], src_chk)
    scores = read(reg / "scores.tsv")
    src_sc = f"ANALYSIS/Arabidopsis1135/regions/{p}/scores.tsv"
    gea = [r for r in scores if r["scan"].startswith("gea.")]
    check(len(gea) == 18, f"{p}: 18 GEA rows in scores.tsv")
    best = min(gea, key=lambda r: float(r["empirical_p"]))
    m["region test: smallest p of 18 (line 0.0028)"] = (
        f"{float(best['empirical_p']):.3f} ({best['scan'][4:]}, {best['list']})", src_sc)
    m["region test: rows below 0.05 / below 0.0028"] = (
        f"{sum(float(r['empirical_p']) < 0.05 for r in gea)} / {sum(float(r['empirical_p']) < 0.05 / 18 for r in gea)}",
        src_sc)
    ft10 = [r for r in scores if r["scan"] == "gwas.FT10" and r["list"] == "aragwas_ft"][0]
    m["GWAS FT10 regions touching AraGWAS (p)"] = (
        f"{ft10['n_with_item']} of {ft10['n_regions']} ({float(ft10['empirical_p']):.1g})", src_sc)
    gw_snps = read(res / "GWAS/tables/selected_snps.tsv")
    gw_reg = read(res / "GWAS/tables/regions_combined.tsv")
    m["GWAS distinct SNP calls"] = (str(len(gw_snps)), f"{p}_results/GWAS/tables/selected_snps.tsv")
    m["GWAS combined regions"] = (str(len(gw_reg)), f"{p}_results/GWAS/tables/regions_combined.tsv")
    shared = [r for r in read(reg / "shared_regions.tsv") if r["scan"].startswith("gxg.union.")]
    loci = merged_loci(shared)
    m["GEA x GWAS shared loci (union, both traits)"] = (
        f"{len(loci)}: " + "; ".join(f"chr{c} {s / 1e6:.2f}-{e / 1e6:.2f} Mb" for c, s, e in loci),
        f"ANALYSIS/Arabidopsis1135/regions/{p}/shared_regions.tsv")
    near = read(reg / "check/classic_nearest.tsv")
    src_near = f"ANALYSIS/Arabidopsis1135/regions/{p}/check/classic_nearest.tsv"
    for g in CLASSIC:
        gu = [r for r in near if r["arm"] == ARM[p] and r["scan"] == "gea.union" and r["gene"] == g][0]
        gwas = [r for r in near if r["arm"] == ARM[p] and r["scan"].startswith("gwas.") and r["gene"] == g
                and r["distance_bp"] != ""]
        wza = [r for r in near if r["arm"] == "wza_bonf" and r["scan"].startswith("gwas.") and r["gene"] == g
               and r["in_region"] == "TRUE"]
        gd = min((int(r["distance_bp"]) for r in gwas), default=None)
        gtxt = "none" if gd is None else ("inside" if gd == 0 else f"{gd / 1000:.1f} kb")
        m[f"{g}: nearest GEA union call"] = ("inside" if gu["distance_bp"] == "0" else
                                             f"{int(gu['distance_bp']) / 1000:.1f} kb", src_near)
        m[f"{g}: nearest GWAS call; WZA window holds it"] = (f"{gtxt}; {'yes' if wza else 'no'}", src_near)
    summ = read(ANALYSIS / "offset" / p / "summary.tsv")
    src_off = f"ANALYSIS/Arabidopsis1135/offset/{p}/summary.tsv"
    real = [r for r in summ if r["set"].startswith("gea.")]
    check(len(real) == 18, f"{p}: 18 real offset runs")
    m["offsets: level 1 / level 2 (of 18)"] = (
        f"{sum(r['level_1'] == 'TRUE' for r in real)} / {sum(r['level_2'] == 'TRUE' for r in real)}", src_off)
    m["offsets: climate-distance baseline (mean rho)"] = (f"{float(real[0]['climate_mean']):.3f}", src_off)
    for meth in METHODS:
        s = l3[(p, meth)]
        m[f"{meth}: union mean rho vs its 20 draws"] = (
            f"{s['union']:.3f} vs {s['draw_mean']:.3f} ({s['draw_min']:.3f}-{s['draw_max']:.3f})", src_off)
        m[f"{meth}: level 3 p (Holm over 3); gardens above; rank of 21"] = (
            f"{s['p']:.3f} ({s['holm']:.3f}); {s['gardens_above']}/30; {s['rank']}", "compare/level3_holm.tsv")
    tot = {r["project"]: r for r in read(ANALYSIS / "time_memory/totals.tsv")}[p]
    src_tot = "ANALYSIS/Arabidopsis1135/time_memory/totals.tsv"
    m["chain wall time (each mode once)"] = (tot["chain_wall"], src_tot)
    m["chain peak memory (GiB anon, mode)"] = (f"{tot['chain_max_peak_anon_gib']} ({tot['chain_max_peak_run']})",
                                               src_tot)
    return m


def main():
    OUT.mkdir(exist_ok=True)
    ld = ld_single_draws()
    l3 = level3()
    a, b = metrics(P1, ld, l3), metrics(P2, ld, l3)
    check(list(a) == list(b), "both projects give the same metrics")
    write(OUT / "iteration_compare.tsv", ["metric", P1, P2, "source_iteration_1", "source_iteration_2"],
          [[k, a[k][0], b[k][0], a[k][1], b[k][1]] for k in a])
    print(f"single-draw critical distance on {P2}: {ld['single_min']}-{ld['single_max']} bp over {ld['n']} draws; "
          f"draw 1 (seed 1, the old default) {ld['draw1_bp']} bp, rank {ld['draw1_rank']} of {ld['n']}")


if __name__ == "__main__":
    main()
