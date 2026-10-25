#!/usr/bin/env python3
"""Express check (2026-10-14): how much of the combined-Manhattan "stripes" of Trifolium86Main comes from sNMF
imputation?

A stripe = many SNPs sharing one imputed genotype vector (folded 0<->2), so they get one p value wherever they sit.
LEA mode-fill never changes an observed call, so SNPs of one imputed pattern already agree on every observed cell;
the question is how much of their identity comes from filled cells.

Per tested SNP (imputed VCF = the matrix EMMAX / LFMM / RDA test; observed = the filtered VCF before imputation):
  miss          missing fraction of observed calls (plink --geno 0.3 caps it at 0.3)
  empty_sites   sites with zero observed calls at this SNP (filled entirely from Q)
  pattern_n     tested SNPs genome-wide sharing the folded imputed vector
  mac           minor-allele carriers (samples with the minor homozygote, imputed vector; no hets exist)
  carriers_obs  how many of those carriers are observed calls (rest are fills)
  noncar_filled share of non-carrier samples whose major call is a fill
  one_site      the carriers are exactly all samples of one site
  sig_<method>  significant in that method (any predictor), from the pipeline's sig tables

The missingness sweep SUBSETS the existing results; a real rerun at a stricter snp_miss also changes Q, kinship,
the Bonferroni count and the imputation itself. It is an approximation, labelled as such in the output.

Run from the repo root: python3 -I ANALYSIS/Trifolium86Main/imputation_check/stripe_imputation.py
  [PROJECT [OUT_SUBDIR [IMP_VCF_BASENAME]]]   defaults: Trifolium86Main tables Trifolium86_chr1-6_K3imp.vcf
(2026-10-14: parameterized for the Trifolium86MainImpMAF rerun; observed SNPs absent from the imputed VCF --
removed by sNMF.imputed_maf -- are skipped and counted.)
"""
import collections
import csv
import os
import sys

ROOT = os.getcwd()
B = os.path.dirname(os.path.abspath(__file__))
ARGS = sys.argv[1:] + [None] * 3
PROJECT = ARGS[0] or "Trifolium86Main"
OUT = os.path.join(B, ARGS[1] or "tables")
os.makedirs(OUT, exist_ok=True)
R = os.path.join(ROOT, f"{PROJECT}_results")
W = os.path.join(R, "_work/maf0.1_miss0.3_smiss0.8_rel0.999")
VCF_IMP = os.path.join(W, "climate", ARGS[2] or "Trifolium86_chr1-6_K3imp.vcf")
VCF_OBS = os.path.join(W, "Trifolium86_chr1-6.vcf")
META = os.path.join(ROOT, "data/Trifolium_location_86.tsv")
M = os.path.join(R, "GEA/tables/methods")
with open(os.path.join(R, "pipeline_summary.tsv")) as fh:  # sig tables are named by K_best
    K = next(l.split("\t")[2].strip() for l in fh if l.startswith("structure\tK_best\t"))
SIG = {}  # every method with exactly one SNP sig table (WZA excluded); 2026-10-14: was hardcoded EMMAX/LFMM/RDA
for m in sorted(os.listdir(M)):
    hits = [f for f in os.listdir(os.path.join(M, m)) if f.startswith(f"{m}_pvalues_K{K}_sig_snps_") and "_wza_" not in f]
    if hits:
        assert len(hits) == 1, (m, hits)
        SIG[m] = f"{m}/{hits[0]}"
CODE = {"0/0": "0", "1/1": "2", "0/1": "1", "1/0": "1", "./.": ".", "0|0": "0", "1|1": "2", ".": "."}
FLIP = str.maketrans("02", "20")
STRIPE = 10  # pattern_n >= STRIPE counts as a stripe


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def gt_vectors(path):
    samples = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            yield samples, f[0], int(f[1]), f[2], f[3], "".join(CODE[x.split(":")[0]] for x in f[9:])


site_of = {r["sample"]: r["site"] for r in rows(META)}

# ---- significance per method, and per (method, predictor) ------------------------------------------------
sig = collections.defaultdict(set)
sig_pred = collections.defaultdict(set)
for m, f in SIG.items():
    for r in rows(os.path.join(M, f)):
        sig[m].add(r["SNPID"])
        sig_pred[(m, r["trait"])].add(r["SNPID"])
sig_any = set().union(*sig.values())

# ---- pass 1: imputed vectors ------------------------------------------------------------------------------
imp, pos = {}, {}
ref_imp = {}
for smp, c, p, sid, ref, v in gt_vectors(VCF_IMP):
    assert "." not in v, f"imputed VCF still has a missing call at {sid}"
    imp[sid] = v
    ref_imp[sid] = ref
    pos[sid] = (c, p)
    samples = smp
sites = [site_of[s] for s in samples]
n_site = collections.Counter(sites)
site_members = {st: frozenset(i for i, s in enumerate(sites) if s == st) for st in n_site}
fold = {sid: min(v, v.translate(FLIP)) for sid, v in imp.items()}
pat_n = collections.Counter(fold.values())
chr_len = collections.defaultdict(int)
for c, p in pos.values():
    chr_len[c] = max(chr_len[c], p)

# ---- pass 2: observed calls -------------------------------------------------------------------------------
recs = []
n_flip = 0
n_absent = 0
for smp, c, p, sid, ref, v in gt_vectors(VCF_OBS):
    assert smp == samples
    if sid not in imp:  # dropped by the imputed-MAF filter
        n_absent += 1
        continue
    vi = imp[sid]
    if ref != ref_imp[sid]:  # plink rewrote the imputed VCF with REF/ALT swapped at this SNP
        vi = vi.translate(FLIP)
        imp[sid] = vi
        n_flip += 1
    assert all(o == "." or o == x for o, x in zip(v, vi)), f"imputation changed an observed call at {sid}"
    n = len(v)
    miss = v.count(".") / n
    empty = sum(1 for st, idx in site_members.items() if all(v[i] == "." for i in idx))
    n2 = vi.count("2")
    minor = "2" if n2 <= n - n2 else "0"
    car = [i for i, x in enumerate(vi) if x == minor]
    non = [i for i, x in enumerate(vi) if x != minor]
    car_obs = sum(1 for i in car if v[i] != ".")
    non_fill = sum(1 for i in non if v[i] == ".") / len(non) if non else 0.0
    cset = frozenset(car)
    one_site = next((st for st, idx in site_members.items() if idx == cset), "")
    recs.append(dict(SNPID=sid, chr=c, pos=p, rel=p / chr_len[c], miss=miss, empty_sites=empty,
                     pattern_n=pat_n[fold[sid]], mac=len(car), carriers_obs=car_obs,
                     noncar_filled=non_fill, one_site=one_site,
                     **{f"sig_{m}": int(sid in sig[m]) for m in SIG}, sig_any=int(sid in sig_any)))
assert len(recs) == len(imp)

with open(os.path.join(OUT, "snp_stripe_imputation.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(recs[0]), delimiter="\t")
    w.writeheader()
    for r in recs:
        w.writerow({k: (f"{x:.4f}" if isinstance(x, float) else x) for k, x in r.items()})


def summ(sel):
    n = len(sel)
    if not n:
        return dict(n=0)
    f = lambda key: sum(r[key] for r in sel) / n
    return dict(n=n,
                miss_mean=f"{f('miss'):.3f}",
                share_miss0=f"{sum(r['miss'] == 0 for r in sel) / n:.3f}",
                share_empty_site=f"{sum(r['empty_sites'] > 0 for r in sel) / n:.3f}",
                noncar_filled_mean=f"{f('noncar_filled'):.3f}",
                carriers_filled_share=f"{1 - sum(r['carriers_obs'] for r in sel) / max(1, sum(r['mac'] for r in sel)):.3f}",
                share_one_site=f"{sum(bool(r['one_site']) for r in sel) / n:.3f}",
                share_stripe=f"{sum(r['pattern_n'] >= STRIPE for r in sel) / n:.3f}",
                share_sig=f"{f('sig_any'):.4f}")


def write(name, table):
    keys = list(dict.fromkeys(k for row in table for k in row))
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(table)
    print(f"\n== {name}")
    print("\t".join(keys))
    for row in table:
        print("\t".join(str(row.get(k, "")) for k in keys))


stripe = lambda r: r["pattern_n"] >= STRIPE

# 1. groups: stripe sig / non-stripe sig / background (non-sig), per method and pooled
groups = []
for m in list(SIG) + ["any"]:
    key = f"sig_{m}"
    for g, sel in (("sig_stripe", [r for r in recs if r[key] and stripe(r)]),
                   ("sig_nonstripe", [r for r in recs if r[key] and not stripe(r)])):
        groups.append(dict(method=m, group=g, **summ(sel)))
groups.append(dict(method="-", group="background_stripe", **summ([r for r in recs if not r["sig_any"] and stripe(r)])))
groups.append(dict(method="-", group="background_nonstripe", **summ([r for r in recs if not r["sig_any"] and not stripe(r)])))
groups.append(dict(method="-", group="all_tested", **summ(recs)))
write("groups.tsv", groups)

# 2. dose-response by missing fraction
bins = [(0, 0), (1e-9, 0.05), (0.05, 0.1), (0.1, 0.15), (0.15, 0.2), (0.2, 0.25), (0.25, 0.3001)]
dose = []
for lo, hi in bins:
    sel = [r for r in recs if (r["miss"] == 0 if hi == 0 else lo < r["miss"] <= hi)]
    n = len(sel)
    dose.append(dict(miss_bin="0" if hi == 0 else f"({lo:.2f},{min(hi, 0.3):.2f}]", n=n,
                     share_stripe=f"{sum(map(stripe, sel)) / n:.3f}" if n else "",
                     share_sig=f"{sum(r['sig_any'] for r in sel) / n:.4f}" if n else "",
                     share_sig_stripe=f"{sum(r['sig_any'] and stripe(r) for r in sel) / n:.4f}" if n else "",
                     share_one_site=f"{sum(bool(r['one_site']) for r in sel) / n:.3f}" if n else ""))
write("dose_by_missingness.tsv", dose)

# 3. approximate snp_miss sweep (subset of existing results)
sweep = []
for cut in (0.3, 0.2, 0.1, 0.05, 0.0):
    keep = [r for r in recs if r["miss"] <= cut + 1e-12]
    pn = collections.Counter(fold[r["SNPID"]] for r in keep)
    big = sorted(pn.values(), reverse=True)
    row = dict(snp_miss_max=cut, tested=len(keep),
               sig_any=sum(r["sig_any"] for r in keep),
               sig_in_stripe_after=sum(r["sig_any"] and pn[fold[r["SNPID"]]] >= STRIPE for r in keep),
               largest_pattern=big[0] if big else 0,
               patterns_ge10=sum(1 for x in big if x >= STRIPE))
    for m in SIG:
        row[f"sig_{m}"] = sum(r[f"sig_{m}"] for r in keep)
    sweep.append(row)
write("approx_snp_miss_sweep.tsv", sweep)

# 4. the largest stripes: are their members observed or filled?
by_pat = collections.defaultdict(list)
for r in recs:
    by_pat[fold[r["SNPID"]]].append(r)
top = sorted(by_pat.values(), key=len, reverse=True)[:15]
tops = []
for mem in top:
    s = summ(mem)
    tops.append(dict(pattern_n=len(mem), n_chr=len({r["chr"] for r in mem}), mac=mem[0]["mac"],
                     one_site=mem[0]["one_site"] or "multi", sig_any=sum(r["sig_any"] for r in mem),
                     miss_mean=s["miss_mean"], share_miss0=s["share_miss0"],
                     share_empty_site=s["share_empty_site"], noncar_filled_mean=s["noncar_filled_mean"],
                     carriers_filled_share=s["carriers_filled_share"],
                     share_chr_end=f"{sum(r['rel'] < 0.1 or r['rel'] > 0.9 for r in mem) / len(mem):.3f}"))
write("top_stripes.tsv", tops)

# 5. relative chromosome position deciles
pos_t = []
for d in range(10):
    sel = [r for r in recs if d / 10 <= r["rel"] < (d + 1) / 10 or (d == 9 and r["rel"] == 1)]
    s = summ(sel)
    pos_t.append(dict(decile=d + 1, n=s["n"], miss_mean=s["miss_mean"], share_stripe=s["share_stripe"],
                      share_sig=s["share_sig"],
                      share_sig_stripe=f"{sum(r['sig_any'] and stripe(r) for r in sel) / len(sel):.4f}"))
write("by_chr_position_decile.tsv", pos_t)

# 6. per (method, predictor): stripe share of sig SNPs, and which single site carries them
per = []
rec_of = {r["SNPID"]: r for r in recs}
for (m, t), ids in sorted(sig_pred.items()):
    sel = [rec_of[i] for i in ids if i in rec_of]
    sites_c = collections.Counter(r["one_site"] for r in sel if r["one_site"])
    per.append(dict(method=m, predictor=t, n_sig=len(sel),
                    share_stripe=f"{sum(map(stripe, sel)) / len(sel):.3f}" if sel else "",
                    miss_mean=summ(sel).get("miss_mean", ""),
                    share_one_site=f"{sum(bool(r['one_site']) for r in sel) / len(sel):.3f}" if sel else "",
                    top_one_site=";".join(f"{k}:{v}" for k, v in sites_c.most_common(2))))
write("per_method_predictor.tsv", per)

print(f"\nINFO: {len(recs)} tested SNPs, {len(samples)} samples, {len(n_site)} sites; "
      f"site sizes {dict(sorted(n_site.items()))}; overall missing {sum(r['miss'] for r in recs) / len(recs):.3f}; REF/ALT swapped in imputed VCF at {n_flip} SNPs; observed SNPs absent from the imputed VCF {n_absent}")
