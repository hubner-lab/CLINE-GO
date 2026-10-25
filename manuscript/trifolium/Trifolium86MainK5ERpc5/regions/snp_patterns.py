#!/usr/bin/env python3
"""Genotype-pattern sharing of the significant SNPs of the 15 EMMAX+RDA regions of Trifolium86MainK5ERpc5 (5 PCs)
(plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 7). COPY of ANALYSIS/Trifolium86MainK5/regionsER/snp_patterns.py
(PC3 run); inputs switched to the PC5 project (same imputed matrix, 177,259 SNPs).

Per significant SNP (tables/region_sig_snps.tsv, run evidence_inputs.py first):
  pattern_id / pattern_n  folded imputed genotype vector (vector and its 0<->2 complement are one); pattern_n = tested
                          SNPs genome-wide with that vector. Identical vectors get identical p in every method.
  pattern_n_chr / pattern_regions  chromosomes / which of the 15 regions hold a SNP of the pattern
  locus_unique            pattern_n SNPs all inside this region (a local haplotype, not a genome-wide site indicator)
  carrier_signature       sites where the minor-allele frequency among OBSERVED calls is >= 0.5
  signature_n / _chr / _sig  tested SNPs genome-wide with that signature, their chromosomes, how many are significant
                          (EMMAX or RDA)
  carrier_sites           minor-allele carriers per site among imputed calls, site(n carriers/n samples)
Reported only, never used to prune or re-rank (user rule: pipeline-only confounding control).

Outputs: tables/sig_snp_patterns.tsv
Run: python3 -I ANALYSIS/Trifolium86MainK5ERpc5/regions/snp_patterns.py   (from the repo root)
"""
import collections
import csv
import hashlib
import os

ROOT = os.getcwd()
B = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(ROOT, "Trifolium86MainK5ERpc5_results")
W = os.path.join(P, "_work/maf0.1_miss0.3_smiss0.8_rel0.999")
VCF_IMP = os.path.join(W, "Trifolium86_chr1-6_K5imp_maf0.1.vcf")
VCF_OBS = os.path.join(W, "Trifolium86_chr1-6.vcf")
META = os.path.join(ROOT, "data/Trifolium_location_86.tsv")
M = os.path.join(P, "GEA/tables/methods")
CODE = {"0/0": "0", "1/1": "2", "0/1": "1", "1/0": "1", "./.": "."}
FLIP = str.maketrans("02", "20")


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
            yield samples, f[0], int(f[1]), f[2], "".join(CODE[x.split(":")[0]] for x in f[9:])


site_of = {r["sample"]: r["site"] for r in rows(META)}
reg = {r["region_id"]: (r["chr"], int(r["start"]), int(r["end"])) for r in rows(os.path.join(P, "GEA/tables/regions_combined.tsv"))}
rsig = rows(os.path.join(B, "tables/region_sig_snps.tsv"))
want = {r["snp_id"] for r in rsig}
sig_any = set()
for f in ("EMMAX/EMMAX_pvalues_K5_sig_snps_bonf_0.05.tsv", "RDA/RDA_pvalues_K5_sig_snps_qval_0.05.tsv"):
    sig_any |= {r["SNPID"] for r in rows(os.path.join(M, f))}


def region_of(c, p):
    return [k for k, (rc, s, e) in reg.items() if rc == c and s <= p <= e]


# ---- pass 1: imputed vectors, genome-wide pattern index -------------------------------------------------
pat_of, members, vec_imp, samples = {}, collections.defaultdict(list), {}, None
for samples, c, p, sid, v in gt_vectors(VCF_IMP):
    assert "." not in v
    key = hashlib.md5(min(v, v.translate(FLIP)).encode()).hexdigest()[:10]
    pat_of[sid] = key
    members[key].append((c, p, sid))
    if sid in want:
        vec_imp[sid] = v
assert len(pat_of) == 177259, len(pat_of)
assert set(vec_imp) == want
assert all(s in site_of for s in samples)
sites = [site_of[s] for s in samples]
n_site = collections.Counter(sites)
site_idx = {st: [i for i, s in enumerate(sites) if s == st] for st in sorted(set(sites))}
print(f"INFO: {len(samples)} samples, {len(n_site)} sites, {len(members)} distinct folded patterns")


def minor_of(v):
    called = [x for x in v if x != "."]
    alt = sum(1 for x in called if x == "2") + 0.5 * sum(1 for x in called if x == "1")
    return "2" if alt <= len(called) / 2 else "0"


def signature(v):
    m = minor_of(v)
    out = []
    for st, idx in site_idx.items():
        cs = [v[i] for i in idx if v[i] != "."]
        if cs and sum(1.0 if x == m else 0.5 if x == "1" else 0.0 for x in cs) / len(cs) >= 0.5:
            out.append(st)
    return "+".join(out) if out else "none"


def carrier_sites(v):
    m = minor_of(v)
    c = collections.Counter(s for s, x in zip(sites, v) if x in (m, "1"))
    return ";".join(f"{k}({c[k]}/{n_site[k]})" for k in sorted(c, key=lambda k: (-c[k], k)))


# ---- pass 2: observed calls of the TESTED SNPs; signatures genome-wide -----------------------------------
sg_n, sg_sig, sg_chr, vec_obs = collections.Counter(), collections.Counter(), collections.defaultdict(set), {}
for smp, c, p, sid, v in gt_vectors(VCF_OBS):
    assert smp == samples
    if sid not in pat_of:
        continue
    g = signature(v)
    sg_n[g] += 1
    sg_chr[g].add(c)
    if sid in sig_any:
        sg_sig[g] += 1
    if sid in want:
        vec_obs[sid] = v
assert sum(sg_n.values()) == 177259 and set(vec_obs) == want

with open(os.path.join(B, "tables/sig_snp_patterns.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["region_id", "snp_id", "pattern_id", "pattern_n", "pattern_n_chr", "pattern_regions", "locus_unique",
                "pattern_n_sig", "carrier_signature", "signature_n", "signature_n_chr", "signature_n_sig", "carrier_sites"])
    for r in rsig:
        sid = r["snp_id"]
        key = pat_of[sid]
        mem = members[key]
        prs = sorted({x for c, p, _ in mem for x in region_of(c, p)})
        rc, s, e = reg[r["region_id"]]
        lu = all(c == rc and s <= p <= e for c, p, _ in mem)
        sg = signature(vec_obs[sid])
        w.writerow([r["region_id"], sid, key, len(mem), len({c for c, _, _ in mem}), ",".join(prs), lu,
                    sum(1 for _, _, x in mem if x in sig_any), sg, sg_n[sg], len(sg_chr[sg]), sg_sig[sg],
                    carrier_sites(vec_imp[sid])])
print(f"INFO: sig_snp_patterns.tsv written for {len(rsig)} significant SNPs")
