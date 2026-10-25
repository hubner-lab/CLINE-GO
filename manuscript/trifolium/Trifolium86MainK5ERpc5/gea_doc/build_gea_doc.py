#!/usr/bin/env python3
"""Builds GEA_RESULTS_METHODS_v2.md, the GEA Methods/Results drafting document of the PC5 set (plan
~/.claude/plans/peaceful-kindling-micali.md, Phase 8), from the hand-written GEA_RESULTS_METHODS_v2.src.md.

Every number in the text is a token, filled here from the tables of BOTH runs, so the change log is checked the
same way as the v2 text:
  {{v5.key}} / {{v3.key}}   a value computed below for the PC5 (Trifolium86MainK5ERpc5) or PC3 (Trifolium86MainK5ER) run
  {{c.key}}                 a config value (config_Trifolium86MainK5ERpc5.yaml) or a run-independent value
  {{dist:gNNN}}             distance of a v2 region gene to its nearest significant SNP (region_gene_evidence.tsv)
  {{M7}}                    the approved "Choice of K" paragraph, copied verbatim from manuscript_offset_v2
  {{REGION_TABLE}}          one row per v2 region: span, genes, SNPs, best p, predictors, status, support (the support
                            cell is copied from the overview table of regions/REGIONS_v2.md)
Checks (exit 1 on any failure):
  - every token resolves and every computed v5 key is used at least once (an unused key = a number left out);
  - no digit outside a token, except inside identifiers (SNP ids, gene ids, region ids, bio_N, PCn, PMIDs, Author Year,
    file names) — so no hand-typed number can enter the text;
  - outside the "## Change log" section no PC3 value appears: every {{v3.*}} token sits inside that section (except
    {{v3.n_pcs}}, the PC3 setting that names the run);
  - the M7 paragraph equals the source paragraph byte for byte;
  - gene ids (gNNN) and SNP ids (chr:pos) in a region bullet belong to that region (genes: region_gene_evidence.tsv or
    neighbourhood.tsv; SNPs: region_sig_snps.tsv); elsewhere they belong to some v2 region.
Writes gea_doc_numbers.tsv (every key, value, source) beside this script.
Run from the repo root: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/gea_doc/build_gea_doc.py
"""
import collections
import csv
import math
import os
import re
import statistics
import sys
from statistics import NormalDist

ROOT = os.getcwd()
B = os.path.dirname(os.path.abspath(__file__))
AN5 = os.path.dirname(B)
AN3 = os.path.join(ROOT, "ANALYSIS/Trifolium86MainK5")
SRC = os.path.join(B, "GEA_RESULTS_METHODS_v2.src.md")
OUT = os.path.join(AN5, "GEA_RESULTS_METHODS_v2.md")
NUM = os.path.join(B, "gea_doc_numbers.tsv")
RUNS = {
    "v5": dict(proj="Trifolium86MainK5ERpc5", reg=os.path.join(AN5, "regions/tables")),
    "v3": dict(proj="Trifolium86MainK5ER", reg=os.path.join(AN3, "regionsER/tables")),
}
CONFIG = os.path.join(ROOT, "config_Trifolium86MainK5ERpc5.yaml")
fail = []


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ints(n):
    return f"{n:,}"


def pfmt(x, sig=2):
    m, e = f"{float(x):.{sig - 1}e}".split("e")
    return f"{m}e{int(e)}"


def kb(bp):
    bp = int(bp)
    return f"{bp / 1000:.1f} kb" if bp < 10000 else f"{round(bp / 1000)} kb"


N, SRCS = {}, {}


def put(key, val, src):
    N[key], SRCS[key] = str(val), src


# ---- config (flat grep: the YAML is simple and PyYAML is not in the stdlib) --------------------------------------
cfg = open(CONFIG).read()


def cget(pat):
    m = re.search(pat, cfg, re.M)
    if not m:
        sys.exit(f"config: {pat} not found")
    return m.group(1).strip().strip("'\"")


put("c.maf", cget(r"^  maf:\s*(\S+)"), "Filter.maf")
put("c.snp_miss", cget(r"^  snp_miss:\s*(\S+)"), "Filter.snp_miss")
put("c.sample_miss", cget(r"^  sample_miss:\s*(\S+)"), "Filter.sample_miss")
put("c.relatedness", cget(r"^  relatedness:\s*(\S+)"), "Filter.relatedness")
put("c.ld_window", cget(r"^LD:\n  window:\s*(\S+)"), "LD.window")
put("c.ld_step", cget(r"^LD:\n  window:.*\n  step:\s*(\S+)"), "LD.step")
put("c.ld_r2", cget(r"^LD:\n(?:  .*\n){2}  r2:\s*(\S+)"), "LD.r2")
put("c.K", cget(r"^  k_best:\s*(\S+)"), "sNMF.k_best")
put("c.snmf_repeats", cget(r"^  repeats:\s*(\S+)"), "sNMF.repeats")
put("c.n_pcs", cget(r"^\s+n_pcs:\s*(\d+)"), "GEA EMMAX n_pcs")
put("c.condition_pcs", cget(r"condition_pcs:\s*(\d+)\n\s*predictor_set"), "GEA RDA condition_pcs")
put("c.alpha", cget(r"method: EMMAX\n\s*adjust: bonf\n\s*threshold:\s*(\S+)"), "GEA EMMAX threshold")
put("c.rda_q", cget(r"method: RDA\n\s*adjust: qval\n\s*threshold:\s*(\S+)"), "GEA RDA threshold")
put("c.rda_perm", cget(r"axis_alpha: 0.05\n\s*condition_pcs: \d+\n.*\n\s*permutations:\s*(\d+)"), "GEA RDA permutations")
put("c.clump_bp", ints(int(cget(r"snp_clumping_distance:\s*(\d+)"))), "GEA.snp_clumping_distance")
put("c.promoter_kb", int(cget(r"^GEA:(?:\n.*)*?\n  promoter_length:\s*(\d+)")) // 1000, "GEA.promoter_length")
pred16 = cget(r"^Climate:\n  enabled:.*\n  predictors:\s*(\S+)").split(",")
put("c.n_pred", len(pred16), "Climate.predictors")
rda_pred = cget(r"predictor_set:\s*(\S+)").split(",")
put("c.rda_pred", ", ".join(rda_pred), "GEA RDA predictor_set")
put("c.n_rda_pred", len(rda_pred), "GEA RDA predictor_set")

# ---- climate facts behind the predictor list (plant-level climate of the 57 analysed plants) --------------------
clim = rows(os.path.join(ROOT, "Trifolium86MainK5ERpc5_results/climate/tables/present/climate_present_site.tsv"))
if len({r["bio_17"] for r in clim}) != len({r["bio_18"] for r in clim}):
    fail.append("bio_17 and bio_18 differ in their number of values")
for v in ("bio_14", "bio_17", "bio_18"):
    put(f"c.nvals_{v}", len({r[v] for r in clim}), f"climate_present_site.tsv distinct {v}")
for a, b in (("bio_8", "bio_11"), ("bio_16", "bio_19")):
    if not all(r[a] == r[b] for r in clim):
        fail.append(f"{a} and {b} differ in climate_present_site.tsv; the 'identical predictors' text is wrong")


# ---- per run ------------------------------------------------------------------------------------------------------
def chi2_lambda(ps):
    m = statistics.median(ps)
    z = NormalDist().inv_cdf(1 - m / 2)
    return z * z / 0.454936


def run(v, proj, reg):
    R = os.path.join(ROOT, f"{proj}_results")
    T = os.path.join(R, "GEA/tables")
    summ = {(r["step"], r["metric"]): r["value"] for r in rows(os.path.join(R, "pipeline_summary.tsv"))}
    if v == "v5":
        for k, m in (("plants_raw", "samples_total"), ("plants", "samples_after_filtering"),
                     ("plants_rm_rel", "samples_removed_relatedness"), ("snps_raw", "snps_raw"),
                     ("snps_filt", "snps_after_filtering"), ("snps_ld", "snps_after_ld_pruning")):
            put(f"v5.{k}", ints(int(summ[("processing", m)])), f"pipeline_summary processing {m}")
        put("v5.plants_rm_miss", int(summ[("processing", "samples_total")]) - int(summ[("processing", "samples_after_filtering")])
            - int(summ[("processing", "samples_removed_relatedness")]), "pipeline_summary: total - after - relatedness")
        put("v5.snps_ld_imp", ints(int(summ[("structure", "snps_after_imputed_maf")])), "pipeline_summary structure")
        put("v5.ld_bp", ints(int(summ[("structure", "ld_decay_r2_02_genome_wide")])), "pipeline_summary structure")
        if N["v5.ld_bp"] != N["c.clump_bp"]:
            fail.append("snp_clumping_distance != genome-wide r2 0.2 distance")
        kept = {r["sample"] for r in clim}
        put("v5.sites", len({r["site"] for r in rows(os.path.join(R, "Processing/tables/metadata.tsv")) if r["sample"] in kept}),
            "metadata.tsv sites of the analysed plants")
    m_imp = int(summ[("gea", "snps_after_imputed_maf")])
    put(f"{v}.snps_imp", ints(m_imp), "pipeline_summary gea snps_after_imputed_maf")
    put(f"{v}.snps_imp_pre", ints(int(summ[("gea", "snps_before_imputed_maf")])), "pipeline_summary gea")
    em = rows(os.path.join(T, "methods/EMMAX/EMMAX_pvalues_K5_sig_snps_bonf_0.05.tsv"))
    rd = rows(os.path.join(T, "methods/RDA/RDA_pvalues_K5_sig_snps_qval_0.05.tsv"))
    thr = float(em[0]["pval_threshold"])
    if abs(thr - 0.05 / m_imp) / thr > 1e-6:
        fail.append(f"{v}: EMMAX threshold {thr} != 0.05 / {m_imp}")
    put(f"{v}.emmax_thr", pfmt(thr, 3), "EMMAX sig table pval_threshold")
    put(f"{v}.rda_thr", pfmt(rd[0]["pval_threshold"], 3), "RDA sig table pval_threshold")
    put(f"{v}.emmax_snps", len({r["SNPID"] for r in em}), "EMMAX sig table, unique SNPID")
    put(f"{v}.emmax_pairs", len(em), "EMMAX sig table rows (SNP x predictor)")
    put(f"{v}.rda_snps", len(rd), "RDA sig table rows")
    cnt = collections.Counter(r["trait"] for r in em).most_common()
    for i, (t, n) in enumerate(cnt[:3], 1):
        put(f"{v}.emmax_top{i}", f"{t} ({n})", "EMMAX sig table, rows per predictor")
    dup = {"bio_11", "bio_19"}  # identical to bio_8 / bio_16 at every plant (checked below)
    put(f"{v}.emmax_pairs_distinct", sum(1 for r in em if r["trait"] not in dup), "EMMAX sig rows minus bio_11 / bio_19 rows")
    put(f"{v}.emmax_n_pred_hit_distinct", len([t for t, _ in cnt if t not in dup]), "predictors hit minus bio_11 / bio_19")
    # lambda per predictor, recomputed from the p-value table (no pipeline table holds it)
    with open(os.path.join(T, "methods/EMMAX/EMMAX_pvalues_K5.tsv")) as fh:
        h = next(fh).rstrip("\n").split("\t")
        cols = h[3:]
        ps = {c: [] for c in cols}
        ident = {("bio_8", "bio_11"): True, ("bio_16", "bio_19"): True}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for c, x in zip(cols, f[3:]):
                ps[c].append(float(x))
            for a, b in ident:
                if f[h.index(a)] != f[h.index(b)]:
                    ident[(a, b)] = False
    if not all(ident.values()):
        fail.append(f"{v}: EMMAX columns bio_8/bio_11 or bio_16/bio_19 are not identical")
    if v == "v5":
        put("v5.n_tests", len(cols) - sum(ident.values()), "EMMAX p table: predictors minus identical column pairs")
    lam = {c: chi2_lambda(p) for c, p in ps.items()}
    lo, hi = min(lam, key=lam.get), max(lam, key=lam.get)
    put(f"{v}.lam_lo", f"{lam[lo]:.2f} ({lo})", "median chi2 of EMMAX p, per predictor")
    put(f"{v}.lam_hi", f"{lam[hi]:.2f} ({hi})", "median chi2 of EMMAX p, per predictor")
    put(f"{v}.lam_med", f"{statistics.median(lam.values()):.2f}", "median over predictors")
    # RDA diagnostics / anova
    dg = {r["key"]: r["value"] for r in rows(os.path.join(T, "methods/RDA/RDA_diagnostics_K5.tsv"))}
    put(f"{v}.rda_axes", dg["rda_axes"], "RDA_diagnostics rda_axes")
    put(f"{v}.rda_adjr2", f"{float(dg['adj_r_squared']):.3f}", "RDA_diagnostics adj_r_squared")
    put(f"{v}.rda_gif", f"{float(dg['gif_lambda']):.2f}", "RDA_diagnostics gif_lambda")
    put(f"{v}.rda_vif", f"{float(dg['max_vif']):.2f}", "RDA_diagnostics max_vif (predictors + conditioning PCs)")
    put(f"{v}.rda_n", dg["n_samples"], "RDA_diagnostics n_samples")
    an = rows(os.path.join(T, "methods/RDA/RDA_anova_K5.tsv"))
    put(f"{v}.rda_model_p", next(r["p"] for r in an if r["term"] == "Model"), "RDA_anova Model p")
    put(f"{v}.rda_axis_p", ", ".join(r["p"] for r in an if r["test"] == "axis" and r["p"]), "RDA_anova axis p")
    # agreement set
    sel = rows(os.path.join(T, "selected_snps.tsv"))
    put(f"{v}.agreed", len(sel), "selected_snps.tsv rows")
    put(f"{v}.agreed_both", sum(1 for r in sel if r["EMMAX"] and r["RDA"]), "selected_snps: EMMAX and RDA")
    put(f"{v}.agreed_e", sum(1 for r in sel if r["EMMAX"] and not r["RDA"]), "selected_snps: EMMAX only")
    put(f"{v}.agreed_r", sum(1 for r in sel if r["RDA"] and not r["EMMAX"]), "selected_snps: RDA only")
    selc = collections.Counter(t for r in sel for t in r["EMMAX"].split(",") if t).most_common()
    put(f"{v}.sel_top2", f"{selc[0][0]} ({selc[0][1]} SNPs) and {selc[1][0]} ({selc[1][1]} SNPs)", "selected_snps EMMAX column")
    reg_rows = rows(os.path.join(T, "regions_combined.tsv"))
    put(f"{v}.regions", len(reg_rows), "regions_combined.tsv rows")
    lens = sorted(int(r["length"]) for r in reg_rows)
    put(f"{v}.reg_kb_range", f"{round(lens[0] / 1000)}-{round(lens[-1] / 1000)} kb", "regions_combined length")
    put(f"{v}.reg_mb_total", f"{sum(lens) / 1e6:.1f} Mb", "regions_combined length sum")
    put(f"{v}.reg_single", sum(1 for r in reg_rows if r["snp_count"] == "1"), "regions_combined snp_count == 1")
    genes = rows(os.path.join(T, "genes_combined.tsv"))
    put(f"{v}.genes", len({g["gene_id"] for g in genes}), "genes_combined.tsv unique gene_id")
    put(f"{v}.regions_per_trait", summ[("gea", "regions_per_trait")], "pipeline_summary gea regions_per_trait")
    put(f"{v}.min_p", pfmt(min(float(r["min_pvalue"]) for r in reg_rows)), "regions_combined min_pvalue")
    # significant SNPs in the regions and their genotype patterns (snp_patterns.py output)
    rs = rows(os.path.join(reg, "region_sig_snps.tsv"))
    pat = {r["snp_id"]: r for r in rows(os.path.join(reg, "sig_snp_patterns.tsv"))}
    put(f"{v}.reg_sig", len(rs), "region_sig_snps.tsv rows")
    put(f"{v}.c29072", sum(1 for r in rs if "29072(4/4)" in pat[r["snp_id"]]["carrier_sites"]),
        "sig_snp_patterns: carrier_sites has 29072(4/4)")
    put(f"{v}.locus_unique", sum(1 for r in rs if pat[r["snp_id"]]["locus_unique"] == "True"), "sig_snp_patterns locus_unique")
    put(f"{v}.bio2_best", sum(1 for r in rs if r["EMMAX_pred"] == "bio_2"), "region_sig_snps EMMAX_pred == bio_2")
    regs_29072 = {r["region_id"] for r in rs if "29072" in pat[r["snp_id"]]["carrier_signature"].split("+")}
    put(f"{v}.reg_29072", len(regs_29072), "regions with a SNP whose carrier_signature includes 29072")
    agreed = [r for r in rs if r["agreed"] == "TRUE"]
    stripe = sum(1 for r in agreed if int(pat[r["snp_id"]]["pattern_n"]) >= 10)
    if len(agreed) != len(sel):
        fail.append(f"{v}: agreed SNPs in region_sig_snps ({len(agreed)}) != selected_snps ({len(sel)})")
    put(f"{v}.stripe_n", stripe, "agreed SNPs with pattern_n >= 10 (compare_sets.py stripe rule)")
    put(f"{v}.stripe_share", f"{stripe / len(agreed):.2f}", "stripe_n / agreed")
    if v == "v5":
        big = max(pat.values(), key=lambda r: int(r["pattern_n"]))
        put("v5.coast_n", big["pattern_n"], f"sig_snp_patterns pattern_n of {big['snp_id']}")
        put("v5.coast_chr", big["pattern_n_chr"], "sig_snp_patterns pattern_n_chr")
        put("v5.coast_sig", big["pattern_n_sig"], "sig_snp_patterns pattern_n_sig")
        put("v5.coast_in_set", sum(1 for r in agreed if pat[r["snp_id"]]["pattern_id"] == big["pattern_id"]),
            "agreed SNPs carrying the largest vector")
        sig29 = {r["signature_n"] for r in pat.values() if r["carrier_signature"] == "29010+29072"}
        put("v5.sig_29010_29072", ints(int(sig29.pop())), "signature_n of 29010+29072")
    return {r["SNPID"] for r in sel}


s5 = run("v5", **RUNS["v5"])
s3 = run("v3", **RUNS["v3"])
cfg3 = open(os.path.join(ROOT, "config_Trifolium86MainK5ER.yaml")).read()
put("v3.n_pcs", re.search(r"^\s+n_pcs:\s*(\d+)", cfg3, re.M).group(1), "config_Trifolium86MainK5ER.yaml EMMAX n_pcs")
put("v3.condition_pcs", re.search(r"condition_pcs:\s*(\d+)\n\s*predictor_set", cfg3).group(1),
    "config_Trifolium86MainK5ER.yaml RDA condition_pcs")
put("v5.shared", len(s5 & s3), "selected_snps v5 & v3")
put("v5.only", len(s5 - s3), "selected_snps v5 - v3")
put("v3.only", len(s3 - s5), "selected_snps v3 - v5")

# ---- regions old/new/lost, stress screen, curation and literature (Phase 2 / 7 tables) ----------------------------
rv = rows(os.path.join(AN5, "regions/regions_v2.tsv"))
put("v5.reg_old", sum(1 for r in rv if r["status"] == "old"), "regions_v2.tsv status old")
put("v5.reg_new", sum(1 for r in rv if r["status"] == "new"), "regions_v2.tsv status new")
put("v5.genes_new_regions", sum(int(r["n_genes"]) for r in rv if r["status"] == "new"), "regions_v2.tsv n_genes, new")
put("v5.reg_indep", len(rv) - 1, "regions minus one: 4_100.9M and 4_101.3M tag one haplotype (REGIONS_v2.md)")
shifts = []
for r in rv:
    if r["status"] == "old" and r["pc3_region"] != r["region_id"]:
        c5, a5, b5 = re.match(r"(\d)_(\d+)-(\d+)", r["region_id"]).groups()
        _, a3, b3 = re.match(r"(\d)_(\d+)-(\d+)", r["pc3_region"]).groups()
        d = []
        if a5 != a3:
            d.append(f"left edge {kb(abs(int(a5) - int(a3)))} {'further left' if int(a5) < int(a3) else 'further right'}")
        if b5 != b3:
            d.append(f"right edge {kb(abs(int(b5) - int(b3)))} {'further right' if int(b5) > int(b3) else 'further left'}")
        shifts.append(f"{r['region_id']} ({'; '.join(d)})")
put("v5.span_changes", ", ".join(shifts), "regions_v2.tsv region_id vs pc3_region")
lost = rows(os.path.join(AN5, "regions/regions_lost.tsv"))
put("v5.reg_lost", len(lost), "regions_lost.tsv rows")
for i, r in enumerate(sorted(lost, key=lambda r: r["region_id"]), 1):
    put(f"v5.lost{i}_genes", r["n_genes"], f"regions_lost.tsv n_genes {r['region_id']}")
st = {r["set"]: r for r in rows(os.path.join(AN5, "regions/set_stress.tsv"))}
for v, s in (("v3", "PC3"), ("v5", "PC5")):
    r = st[s]
    put(f"{v}.go_obs", r["go_obs"], "set_stress.tsv go_obs")
    put(f"{v}.go_exp", r["go_exp"], "set_stress.tsv go_exp")
    put(f"{v}.go_ratio", r["go_ratio"], "set_stress.tsv go_ratio")
    put(f"{v}.go_p", r["go_p"], "set_stress.tsv go_p")
    put(f"{v}.kw_obs", r["kw_strict_obs"], "set_stress.tsv kw_strict_obs")
    put(f"{v}.kw_exp", r["kw_strict_exp"], "set_stress.tsv kw_strict_exp")
    put(f"{v}.kw_ratio", r["kw_strict_ratio"], "set_stress.tsv kw_strict_ratio")
    put(f"{v}.kw_p", r["kw_strict_p"], "set_stress.tsv kw_strict_p")
ev = rows(os.path.join(AN5, "regions/tables/region_gene_evidence.tsv"))
put("v5.genes_carried", sum(1 for r in ev if r["source"] == "carried"), "region_gene_evidence source carried")
put("v5.genes_called_new", sum(1 for r in ev if r["source"] != "carried"), "region_gene_evidence source != carried")
lit = rows(os.path.join(AN5, "regions/tables/literature.tsv"))
put("v5.lit_rows", len(lit), "literature.tsv rows")
put("v5.lit_pmids", len({r["pmid"] for r in lit}), "literature.tsv unique pmid")
put("v5.lit_cited", len({r["pmid"] for r in lit if r["role"] != "rejected"}), "literature.tsv unique pmid, role != rejected")
nbh = rows(os.path.join(AN5, "regions/tables/neighbourhood_guess.tsv"))
put("v5.nbh_genes", len(nbh), "neighbourhood_guess.tsv rows")
dist = {r["gene_id"].replace(".t1", ""): r["dist_to_sig_snp"] for r in ev}

# LFMM (supplementary figure draft, K5PC3 track; LFMM uses no PCs, so it is identical in the PC5 track)
FIG = os.path.join(AN3, "figures")
pb = [r for r in rows(os.path.join(FIG, "figS_lfmm_panel_b.tsv")) if r["K"] == "K = 5"]
lf = {r["class"]: int(r["n"]) for r in pb if r["method"] == "LFMM"}
put("c.lfmm_sig", ints(sum(lf.values())), "figS_lfmm_panel_b K5 LFMM")
for m in ("LFMM", "EMMAX", "RDA"):
    d = {r["class"]: int(r["n"]) for r in pb if r["method"] == m}
    put(f"c.{m.lower()}_k5pc3_stripe", f"{ints(d['in stripe'])} of {ints(sum(d.values()))}", f"figS_lfmm_panel_b K5 {m} in stripe")
gif = [float(r["gif"]) for r in rows(os.path.join(FIG, "figS_lfmm_panel_c.tsv")) if r["K"] == "K = 5"]
put("c.lfmm_gif", f"median {statistics.median(gif):.1f}, range {min(gif):.1f}-{max(gif):.1f}", "figS_lfmm_panel_c K5 median (range)")

# ---- the approved Choice-of-K paragraph -------------------------------------------------------------------------------
mo = open(os.path.join(AN5, "manuscript_offset_v2/manuscript_offset.md")).read()
M7 = re.search(r"### Choice of K\n\n(.+?)\n", mo).group(1)

# ---- the region table ---------------------------------------------------------------------------------------------------
rmd = open(os.path.join(AN5, "regions/REGIONS_v2.md")).read()
support = {}
for line in rmd.splitlines():
    m = re.match(r"\| (\d+) \| (\d_\d+-\d+) \|.*\| ([^|]+) \|$", line)
    if m:
        support[m.group(2)] = m.group(3).strip().replace("as region 8", "as 4_100931231-101192853")
reg5 = rows(os.path.join(ROOT, "Trifolium86MainK5ERpc5_results/GEA/tables/regions_combined.tsv"))
rs5 = rows(os.path.join(AN5, "regions/tables/region_sig_snps.tsv"))
status = {r["region_id"]: r["status"] for r in rv}
tab = ["| Region | kb | Genes | Sig. SNPs | Best p EMMAX / RDA | Best EMMAX predictor of its EMMAX SNPs | At 5 PCs | Candidate gene support |",
       "|---|---|---|---|---|---|---|---|"]
for r in reg5:
    rid = r["region_id"]
    s = [x for x in rs5 if x["region_id"] == rid]
    pe = min(float(x["EMMAX_minp"]) for x in s if x["EMMAX_minp"] not in ("", "NA"))
    pr = min(float(x["RDA_p"]) for x in s if x["RDA_p"] not in ("", "NA"))
    preds = sorted({x["EMMAX_pred"] for x in s if "EMMAX" in x["sig_methods"].split(",")}, key=lambda t: int(t.split("_")[1]))
    ng = sum(1 for g in ev if g["region_id"] == rid)
    tab.append(f"| {rid} | {round(int(r['length']) / 1000)} | {ng} | {len(s)} | {pfmt(pe)} / {pfmt(pr)} | "
               f"{', '.join(preds)} | {status[rid]} | {support[rid]} |")
TABLE = "\n".join(tab)

# ---- render with token markers, check, strip -----------------------------------------------------------------------
src = open(SRC).read()
used = set()
O, C = "\x02", "\x03"


def sub(m):
    k = m.group(1)
    if k == "M7":
        return O + M7 + C
    if k == "REGION_TABLE":
        return O + TABLE + C
    if k.startswith("dist:"):
        g = k[5:]
        if g not in dist:
            fail.append(f"{{{{{k}}}}}: not a v2 region gene")
            return O + "??" + C
        return O + kb(dist[g]) + C
    if k not in N:
        fail.append(f"unknown token {{{{{k}}}}}")
        return O + "??" + C
    used.add(k)
    return O + N[k] + C


out = re.sub(r"\{\{([^}]+)\}\}", sub, src)
for k in sorted(set(N) - used):
    if k.startswith("v5."):
        fail.append(f"v5 key computed but never used: {k} = {N[k]}")

# PC3 values only inside the change log
cl = src.find("\n## Change log")
end = src.find("\n## ", cl + 1)
for m in re.finditer(r"\{\{v3\.(?!n_pcs\}\})", src):  # the PC3 setting itself may be named anywhere
    if not (cl < m.start() < (end if end > 0 else len(src))):
        line = src.count("\n", 0, m.start()) + 1
        fail.append(f"src line {line}: PC3 token outside the change log")

# digits outside tokens
IDENT = re.compile(r"`[^`]*`|\b\d+:\d+\b|\bg\d+\b|\b\d_\d+(?:\.\d+)?M\b|\b\d_\d+-\d+\b|\bbio_\d+\b|\bPC\d+\b|PMID \d+"
                   r"|[A-Z][a-zé]+(?: (?:et al\.|& [A-Z][a-z]+))? (?:19|20)\d\d|\b(?:19|20)\d\d-\d\d-\d\d\b|\bv[35]\b|\bPC[35]\b"
                   r"|\bK\d\b|\bFig(?:ure)?\. S\d+|\bTable S\d+|\bResults \d|\bPhase \d|\b\d\d{3,}\b(?=\)|,? \()"
                   r"|\b2[89]\d{3}\b|\b[A-Za-z]+_\d+\b|\bv\d\b|\bchromosome[ -]\d\b|\b[A-Z]+-\d\b|^#+ .*$|\(\d[a-z]?\)|\b[A-Za-z]+\d+[A-Za-z0-9]*\b", re.M)
plain = re.sub(O + r".*?" + C, " ", out, flags=re.S)
for i, line in enumerate(plain.splitlines(), 1):
    rest = IDENT.sub(" ", line)
    for m in re.finditer(r"\d", rest):
        fail.append(f"out line {i}: untracked digit in: {line.strip()[:140]}")
        break

out = out.replace(O, "").replace(C, "")
# gene and SNP ids: a region bullet names only genes of that region (or its neighbourhood) and significant SNPs of
# that region; anywhere else a gene id must be a v2 region gene and a SNP id a significant SNP of some v2 region
genes_of = collections.defaultdict(set)
for r in ev:
    genes_of[r["region_id"]].add(r["gene_id"].replace(".t1", ""))
for r in rows(os.path.join(AN5, "regions/tables/neighbourhood.tsv")):
    for rid in r["region_ids"].split(","):
        genes_of[rid].add(r["gene_id"].replace(".t1", ""))
snps_of = collections.defaultdict(set)
for r in rs5:
    snps_of[r["region_id"]].add(r["snp_id"])
all_genes, all_snps = set().union(*genes_of.values()), set().union(*snps_of.values())
cur = None
for i, line in enumerate(out.splitlines(), 1):
    if line.startswith("- "):
        m = re.match(r"- \*\*(\d_\d+-\d+)", line)
        cur = m.group(1) if m else None
    elif not line.startswith("  "):
        cur = None
    gs, ss = (genes_of[cur], snps_of[cur]) if cur else (all_genes, all_snps)
    for g in re.findall(r"\bg\d+\b", line):
        if g not in gs:
            fail.append(f"out line {i}: gene {g} not a gene of {cur or 'any v2 region'}")
    for snp in re.findall(r"\b\d:\d{5,}\b", line):
        if snp not in ss:
            fail.append(f"out line {i}: SNP {snp} not a significant SNP of {cur or 'any v2 region'}")
if M7 not in out:
    fail.append("M7 paragraph missing")
with open(OUT, "w") as fh:
    fh.write(out)
with open(NUM, "w") as fh:
    fh.write("key\tvalue\tsource\n")
    for k in sorted(N):
        fh.write(f"{k}\t{N[k]}\t{SRCS[k]}\n")
print(f"wrote {OUT} ({len(out.split())} words), {NUM} ({len(N)} keys)")
if fail:
    print(f"{len(fail)} FAILURES:")
    print("\n".join("  " + f for f in fail))
    sys.exit(1)
print("checks: 0 failures")
