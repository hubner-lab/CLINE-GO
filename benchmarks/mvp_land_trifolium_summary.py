#!/usr/bin/env python3
"""Job C of the landscapes pack: uncorrected RDA on Trifolium beside the conditioned fit.

Reads the published fit (Trifolium86MainK5d150_results, condition_pcs 5) and the two runs of
mvp_land_trifolium_rda.sh under landscape_ms/trifolium_rda_uncorrected/{control,uncorrected}/, compares the control's
p-values with the published ones, and writes summary.tsv (one row per fit) and control.tsv.
"""
import filecmp
import pandas as pd

ROOT = "/mnt/data/eugene/ADAPTOGENE"
PUB = f"{ROOT}/Trifolium86MainK5d150_results/GEA/tables/methods/RDA"
OUT = f"{ROOT}/benchmarks/mvp_eval/landscape_ms/trifolium_rda_uncorrected"
FITS = {"conditioned, published (condition_pcs 5)": PUB, "conditioned, control of this run (condition_pcs 5)": f"{OUT}/control",
        "uncorrected (condition_pcs 0)": f"{OUT}/uncorrected"}


def diag(d):
    out = {}
    with open(f"{d}/RDA_diagnostics_K5.tsv", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                out[p[0]] = p[1]
    return out


rows = []
for name, d in FITS.items():
    g = diag(d)
    pv = pd.read_csv(f"{d}/RDA_pvalues_K5.tsv", sep="\t", float_precision="round_trip", dtype={"chr": str, "SNPID": str})
    p = pv["climate_multivariate"]
    rows.append({"fit": name, "condition_pcs": int(float(g["condition_pcs"])), "predictors": g["predictors"].replace("bio_", "BIO"), "n_tested": int(p.notna().sum()),
                 "n_sig_q_0.05": int(g["n_candidates_q_0.05"]), "n_sig_q_0.10": int(g["n_candidates_q_0.10"]), "n_sig_bonf_0.05": int(g["n_candidates_bonf_0.05"]),
                 "candidate_threshold_p": g.get("candidate_threshold", ""), "candidate_threshold_status": g.get("candidate_threshold_status", ""),
                 "smallest_p": float(p.min()), "gif": float(g["gif_lambda"]), "adj_r_squared": float(g["adj_r_squared"]), "axes_kept": int(float(g["rda_axes"])),
                 "axes_max": int(float(g["rda_axes_max"])), "K_selection": g.get("K_selection", ""), "K_floored": g.get("K_floored", ""),
                 "max_vif": g.get("max_vif", ""), "threads": g.get("cpu", "")})
S = pd.DataFrame(rows)
S.to_csv(f"{OUT}/summary.tsv", sep="\t", index=False)
print(S.T.to_string())
a = pd.read_csv(f"{PUB}/RDA_pvalues_K5.tsv", sep="\t", float_precision="round_trip", dtype={"chr": str, "SNPID": str})
b = pd.read_csv(f"{OUT}/control/RDA_pvalues_K5.tsv", sep="\t", float_precision="round_trip", dtype={"chr": str, "SNPID": str})
same = filecmp.cmp(f"{PUB}/RDA_pvalues_K5.tsv", f"{OUT}/control/RDA_pvalues_K5.tsv", shallow=False)
ids = bool((a.SNPID == b.SNPID).all()) if len(a) == len(b) else False
d = (a.climate_multivariate - b.climate_multivariate).abs()
rel = (d / a.climate_multivariate.clip(lower=1e-300)).max()
C = pd.DataFrame([{"check": "control p-values = published", "byte_identical": same, "snp_ids_identical": ids, "n": len(a), "max_abs_diff": float(d.max()), "max_rel_diff": float(rel),
                   "cand_published": FITS and rows[0]["n_sig_q_0.05"], "cand_control": rows[1]["n_sig_q_0.05"], "axes_published": rows[0]["axes_kept"], "axes_control": rows[1]["axes_kept"]}])
C.to_csv(f"{OUT}/control.tsv", sep="\t", index=False)
print(C.T.to_string())
