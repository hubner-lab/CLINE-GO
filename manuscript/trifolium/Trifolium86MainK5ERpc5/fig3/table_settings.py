#!/usr/bin/env python3
"""t_settings (supplement build 2026-10-25, C2): the settings the Trifolium86MainK5ERpc5 run used, one row per setting.

The configuration is parsed with PyYAML (never regex): yaml.compose() for the values as the file writes them (`yes`,
`30000.0`, `bio_1,bio_3,bio_5`), and Snakemake's own load_configfile() for the data, which is checked against the
`Config MD5` Snakemake logged for every run that named the file (md5 of json.dumps(config + mode, sort_keys=True),
snakemake/workflow.py). List entries are keyed by `method` / `name`, never by position.

Steps follow the framework's sections. Groups read only by modes the run never ran are left out (modes read from
Trifolium86MainK5ERpc5_logs/: a mode ran if its log folder holds files). Keys of a group that ran are all kept; a key
another key switches off keeps its row with a note. Rows outside this configuration (K sweep, LFMM sibling, choice of
K, RDA predictor screen, random-set test) are asserted against the files that record them.

Writes manuscript_fig3/tables/t_settings.{tsv,md} and t_settings_lfmm_sibling_diff.tsv; prints the report facts.
Run inside the image from /pipeline with the Snakemake venv's Python (the image's python3 has no PyYAML):
  /usr/local/pipx/venvs/snakemake/bin/python -I ANALYSIS/Trifolium86MainK5ERpc5/fig3/table_settings.py
"""
import csv
import glob
import hashlib
import json
import os
import re
import sys

import yaml
from snakemake.common.configfile import load_configfile

ROOT = "/pipeline"
AN = f"{ROOT}/ANALYSIS/Trifolium86MainK5ERpc5"
OUT = f"{AN}/manuscript_fig3/tables"
RES = f"{ROOT}/Trifolium86MainK5ERpc5_results"
LOGS = f"{ROOT}/Trifolium86MainK5ERpc5_logs"
CFG = "config_Trifolium86MainK5ERpc5.yaml"
CFG_UP = "config_Trifolium86MainK5ER.yaml"          # processing / prestructure / structure ran under this file
CFG_LFMM = "config_Trifolium86MainK5PC5.yaml"
CFG_KS = "_archive/trifolium_pre_2026-10-17/config_Trifolium86MainKsearch.yaml"
MINUS = "−"


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---- flatten a composed YAML tree to (path, value as written) ------------------------------------------------------
def list_key(node):
    keys = {k.value: v.value for k, v in node.value if isinstance(v, yaml.ScalarNode)}
    for k in ("method", "name"):
        if k in keys:
            return k, keys[k]
    raise ValueError("list entry without method/name")


def flatten(node, path=""):
    if isinstance(node, yaml.MappingNode):
        for k, v in node.value:
            yield from flatten(v, f"{path}.{k.value}" if path else k.value)
    elif isinstance(node, yaml.SequenceNode):
        if not node.value:
            yield path, "[]"
        elif all(isinstance(v, yaml.MappingNode) for v in node.value):
            for v in node.value:
                k, name = list_key(v)
                p = f"{path}[{name}]"
                for kk, vv in v.value:
                    if kk.value == k:
                        continue
                    yield from flatten(vv, f"{p}.{kk.value}")
        else:
            yield path, ", ".join(v.value for v in node.value)
    else:
        yield path, node.value if node.value != "" else "''"


def flat_data(d, path=""):
    """The same flattening over parsed data, for the sibling diff (lists keyed by method/name)."""
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat_data(v, f"{path}.{k}" if path else str(k)))
    elif isinstance(d, list) and d and all(isinstance(v, dict) for v in d):
        for v in d:
            k = "method" if "method" in v else "name"
            for kk, vv in v.items():
                if kk != k:
                    out.update(flat_data(vv, f"{path}[{v[k]}].{kk}"))
    else:
        out[path] = d
    return out


def config_md5(cfg, mode):
    c = load_configfile(cfg)
    c["mode"] = mode
    return hashlib.md5(json.dumps(c, sort_keys=True).encode("utf-8")).hexdigest()


os.chdir(ROOT)
report = []
def rep(s):
    report.append(s)
    print("REPORT:", s)


# ---- 1. the runs: Snakemake logs naming the three configs, their command lines and Config MD5 ---------------------
runs = []
for f in sorted(glob.glob(".snakemake/log/*.snakemake.log")):
    with open(f, errors="replace") as fh:
        head = [next(fh, "") for _ in range(60)]
    cmd = next((l.split("Command:", 1)[1].strip() for l in head if "Command:" in l), "")
    m = re.search(r"--configfile (\S+)", cmd)
    if not m or m.group(1) not in (CFG, CFG_UP, CFG_LFMM):
        continue
    md5 = next((l.split("Config MD5:", 1)[1].strip() for l in head if "Config MD5:" in l), "")
    mode = re.search(r"mode=(\S+)", cmd).group(1)
    runs.append(dict(log=os.path.basename(f), cfg=m.group(1), mode=mode, cmd=cmd, md5=md5,
                     md5_today=config_md5(m.group(1), mode)))
for r in runs:
    extra = re.sub(r"^/usr/local/bin/snakemake ", "", r["cmd"])
    rep(f"run {r['log'][:17]} {r['cfg']} mode={r['mode']}: logged MD5 {r['md5'][:8]} "
        f"{'= ' if r['md5'] == r['md5_today'] else '!= '}today's file | {extra}")
this = [r for r in runs if r["cfg"] == CFG]
assert {r["mode"] for r in this} == {"gea", "maladaptation"}
assert all(r["md5"] == r["md5_today"] for r in this), "today's config does not parse to what the ERpc5 runs used"
up = [r for r in runs if r["cfg"] == CFG_UP]
assert {"processing", "prestructure", "structure"} <= {r["mode"] for r in up}
# command-line overrides other than mode=
for r in this + up:
    cfgs = re.findall(r"--config (.*?)(?= --|$)", r["cmd"])
    keys = [kv.split("=")[0] for c in cfgs for kv in c.split()]
    assert keys == ["mode"], f"override other than mode= in {r['log']}: {keys}"
flags = sorted({re.sub(r"--config \S+|--configfile \S+|-s Snakefile|--scheduler greedy", "", r["cmd"]
                       .replace("/usr/local/bin/snakemake", "")).strip() for r in this + up})
rep("no --config key other than mode= in any run of the two files; Snakemake flags used: " + " | ".join(flags))

# ---- 2. modes the run ran (log folders with files) ------------------------------------------------------------------
ALL_MODES = ["processing", "prestructure", "structure", "climate", "traits", "pregea", "gea", "gwas", "gea_x_gwas",
             "maladaptation"]
ran = set()
for d in sorted(os.listdir(LOGS)):
    if os.listdir(f"{LOGS}/{d}"):
        ran.add(d.lower())
ran_modes = [m for m in ALL_MODES if m in ran]
not_run = [m for m in ALL_MODES if m not in ran]
rep(f"log folders with files: {sorted(ran)}; modes run: {ran_modes}; never run: {not_run}")
assert ran_modes == ["processing", "prestructure", "structure", "gea", "maladaptation"]
assert {r["mode"] for r in this + up} == set(ran_modes)

# ---- 3. rows of this configuration ------------------------------------------------------------------------------------
with open(CFG) as f:
    root = yaml.compose(f)
rows_cfg = list(flatten(root))
data = load_configfile(CFG)
assert len(rows_cfg) == len(flat_data(data)), "compose and load flatten to different key sets"

STEP = [  # (path prefix, step); first match wins
    ("Input", "data processing"), ("project_name", "data processing"), ("cpu", "data processing"),
    ("Regime", "data processing"), ("Filter", "data processing"), ("LD.", "data processing"),
    ("sNMF", "population structure"), ("Population", "population structure"), ("Piemap", "population structure"),
    ("LDdecay", "population structure"),
    ("Map", "environmental factors"), ("Climate", "environmental factors"), ("Future", "environmental factors"),
    ("GEA.configs", "GEA"), ("GEA.wza", "GEA"),
    ("GEA.combine_method", "candidate regions"), ("GEA.snp_clumping_distance", "candidate regions"),
    ("GEA.promoter_length", "candidate regions"), ("GFF", "candidate regions"),
    ("GEA.snp_sets", "genetic offset"), ("Maladaptation", "genetic offset"),
]
LEFT_OUT = [  # (path prefix, reason)
    ("Climate.Varpart", "read by mode=climate, never run"), ("Climate.dbMEM", "read by mode=climate, never run"),
    ("PreGEA", "read by mode=pregea, never run"), ("GWAS", "read by mode=gwas, never run"),
    ("GEAxGWAS", "read by mode=gea_x_gwas, never run"),
    ("Enrichment", "read by no Snakemake rule (common.smk:408-412 sets ENRICHMENT_* and nothing uses them; "
                   "GO enrichment runs in the Shiny app)"),
]
NOTE = {
    "Population.window_size": "inactive: Population.calc_stats is no",
    "Population.custom_trait_file": "inactive: Population.calc_stats is no",
    "cpu": "processing, prestructure and structure ran under config_Trifolium86MainK5ER.yaml (cpu 90, Snakemake -c80); "
           "gea and maladaptation with -c170",
    "Maladaptation.methods.gradient_forest.spatial_correction": "without: no dbMEM input is read",
    "GEA.snp_clumping_distance": "entered as a fixed value: the genome-wide LD-decay distance (r2 = 0.2, 57 plants)",
}
step_of = lambda p: next((s for pre, s in STEP if p == pre or p.startswith(pre)), None)
left_of = lambda p: next((r for pre, r in LEFT_OUT if p == pre or p.startswith(pre + ".") or p.startswith(pre + "[")
                          or p == pre), None)
out, left = [], []
for p, v in rows_cfg:
    lr = left_of(p)
    if lr:
        left.append((p, lr))
        continue
    s = step_of(p)
    assert s, f"no step for {p}"
    out.append(dict(step=s, setting=p, value=v, scope="this configuration", note=NOTE.get(p, "")))
groups_left = sorted({next(pre for pre, _ in LEFT_OUT if p == pre or p.startswith(pre)) for p, _ in left})
rep(f"left out ({len(left)} keys): " + "; ".join(f"{g} ({next(r for pre, r in LEFT_OUT if pre == g)})" for g in groups_left))
assert all(k in NOTE for k in ("Population.window_size", "Population.custom_trait_file"))
assert data["Population"]["calc_stats"] is False and data["Maladaptation"]["methods"]["gradient_forest"]["spatial_correction"] == "without"

# ---- 4. the upstream modes ran under config_Trifolium86MainK5ER.yaml: their settings from records on disk -------------
up_data = load_configfile(CFG_UP)
diff_up = sorted(k for k in set(flat_data(up_data)) | set(flat_data(data))
                 if flat_data(up_data).get(k) != flat_data(data).get(k))
rep(f"{CFG_UP} (today) differs from {CFG} at: {diff_up}")
assert diff_up == ["GEA.configs[EMMAX].params.n_pcs", "GEA.configs[RDA].params.condition_pcs", "cpu", "project_name"]
# The ER file was edited on 10-15, after the 10-14 upstream runs. Undoing its three dated 10-15 edits (GEA.snp_sets,
# Maladaptation.predictors, rda_offset.axes 3 <- auto) must reproduce the Config MD5 those runs logged: then the
# upstream modes read exactly that data, which differs from this file only at keys no upstream mode reads.
recon = load_configfile(CFG_UP)
recon["GEA"].pop("snp_sets"); recon["Maladaptation"].pop("predictors")
recon["Maladaptation"]["methods"]["rda_offset"]["axes"] = "auto"
def md5_of(c, mode):
    c = dict(c, mode=mode)
    return hashlib.md5(json.dumps(c, sort_keys=True).encode("utf-8")).hexdigest()
up_1014 = [r for r in up if r["log"].startswith("2026-10-14")]
assert {r["mode"] for r in up_1014} == {"processing", "prestructure", "structure", "gea"}
assert all(md5_of(recon, r["mode"]) == r["md5"] for r in up_1014), "10-14 upstream config not reproduced"
fr = flat_data(recon)
fa_ = flat_data(data)
diff_recon = sorted(k for k in set(fr) | set(fa_) if fr.get(k, "<absent>") != fa_.get(k, "<absent>"))
rep("the 10-14 upstream runs' config reproduced exactly (Config MD5 of processing, prestructure, structure, gea): "
    f"{CFG_UP} without its 10-15 edits; it differs from {CFG} at {len(diff_recon)} keys: {diff_recon}")
assert all(k in ("cpu", "project_name", "GEA.configs[EMMAX].params.n_pcs", "GEA.configs[RDA].params.condition_pcs",
                 "Maladaptation.predictors", "Maladaptation.methods.rda_offset.axes") or k.startswith("GEA.snp_sets")
           for k in diff_recon)
summ = {(r["step"], r["metric"]): r["value"] for r in read_tsv(f"{RES}/pipeline_summary.tsv")}
F, L = data["Filter"], data["LD"]
work = f"{RES}/_work/maf{F['maf']}_miss{F['snp_miss']}_smiss{F['sample_miss']}_rel{F['relatedness']}/" \
       f"ld{L['r2']}_win{L['window']}_step{L['step']}"
assert os.path.isdir(work), work
snmf_runs = sorted(os.listdir(f"{work}/Trifolium86_chr1-6.snmf/K5"))
assert [d for d in snmf_runs if d.startswith("run")] and len([d for d in snmf_runs if d.startswith("run")]) == data["sNMF"]["repeats"]
assert summ[("prestructure", "K_range")] == f"{data['sNMF']['k_start']}-{data['sNMF']['k_end']}"
assert summ[("structure", "K_best")] == str(data["sNMF"]["k_best"])
assert float(summ[("structure", "imputed_maf_threshold")]) == F["maf"] == float(summ[("gea", "imputed_maf_threshold")])
assert summ[("structure", "ld_decay_group_by")] == data["LDdecay"]["group_by"]
assert summ[("structure", "ld_decay_scope")] == data["LDdecay"]["scope"]
assert int(summ[("structure", "n_climate_variables")]) == len(data["Climate"]["predictors"].split(","))
assert int(summ[("processing", "samples_total")]) == 86 and int(summ[("processing", "samples_after_filtering")]) == 57
rep(f"upstream settings confirmed on disk: _work path {os.path.relpath(work, RES)}; sNMF {sum(d.startswith('run') for d in snmf_runs)} runs "
    f"under K5; pipeline_summary K_range {summ[('prestructure', 'K_range')]}, K_best {summ[('structure', 'K_best')]}, "
    f"imputed MAF {summ[('structure', 'imputed_maf_threshold')]}, LD decay {data['LDdecay']['group_by']}/{data['LDdecay']['scope']}, "
    f"{summ[('structure', 'n_climate_variables')]} climate variables")

# ---- 5. the clumping distance = the genome-wide LD-decay distance ------------------------------------------------------
ld = [r for r in read_tsv(f"{RES}/Structure/tables/ld_decay_half_distances.tsv")
      if r["group"] == "All" and r["scope"] == "genome_wide"]
assert len(ld) == 1 and int(ld[0]["n_samples"]) == 57
assert data["GEA"]["snp_clumping_distance"] == int(ld[0]["r2_02_bp"]) == 130811
assert summ[("structure", "ld_decay_r2_02_genome_wide")] == "130811" == summ[("gea", "region_distance_spec")]
rep("GEA.snp_clumping_distance 130811 = ld_decay_half_distances.tsv All/genome_wide r2_02_bp (57 plants, method "
    f"{ld[0]['method']}); the table stores round() of the fitted distance to whole bp (scripts/ld_decay_analyze.R:200, "
    "unchanged since the run), so the config value is that rounded number; the unrounded value is not on disk")

# ---- 6. rows outside this configuration --------------------------------------------------------------------------------
ks = load_configfile(CFG_KS)
ce = read_tsv(f"{ROOT}/ANALYSIS/Trifolium86MainKsearch/ksearch/cross_entropy.tsv")
assert (ks["sNMF"]["k_start"], ks["sNMF"]["k_end"], ks["sNMF"]["repeats"]) == (2, 8, 10)
assert [int(r["K"]) for r in ce] == list(range(2, 9))
out.append(dict(step="population structure", setting="K sweep (sNMF.k_start, sNMF.k_end, sNMF.repeats)",
                value="2, 8, 10", scope="outside this configuration",
                note=f"its own configuration ({os.path.basename(CFG_KS)}): K 2–8, 10 repeats each; per K the minimum cross-entropy"))
eig = read_tsv(f"{AN}/fig3/offset_minou/paper_pack/k_choice_eigen.tsv")
twp = {int(r["pc"]): float(r["tw_p"]) for r in eig}
first_ns = min(pc for pc, p in twp.items() if p >= 0.05)
assert first_ns == 6 and all(twp[pc] < 0.05 for pc in range(1, 6)) and data["sNMF"]["k_best"] == 5
out.append(dict(step="population structure", setting="choice of K (sNMF.k_best)", value="5", scope="outside this configuration",
                note="from the Tracy–Widom test of the run's PCA: PC1–PC5 p < 0.05, PC6 the first with p ≥ 0.05"))

lf_root = yaml.compose(open(CFG_LFMM))
lf_rows = [(p, v) for p, v in flatten(lf_root) if p.startswith("GEA.configs[LFMM]")]
assert lf_rows, "no LFMM entry in the sibling config"
for p, v in lf_rows:
    out.append(dict(step="GEA", setting=p, value=v, scope="outside this configuration",
                    note="LFMM ran in the sibling project Trifolium86MainK5PC5 (config_Trifolium86MainK5PC5.yaml)"))
fa, fb = flat_data(data), flat_data(load_configfile(CFG_LFMM))
sib = [dict(setting=k, Trifolium86MainK5ERpc5=json.dumps(fa.get(k, "<absent>")), Trifolium86MainK5PC5=json.dumps(fb.get(k, "<absent>")))
       for k in sorted(set(fa) | set(fb)) if fa.get(k, "<absent>") != fb.get(k, "<absent>")]
with open(f"{OUT}/t_settings_lfmm_sibling_diff.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["setting", "Trifolium86MainK5ERpc5", "Trifolium86MainK5PC5"], delimiter="\t")
    w.writeheader(); w.writerows(sib)
rep(f"keys that differ between {CFG} and {CFG_LFMM} ({len(sib)}): " + "; ".join(
    f"{d['setting']}: {d['Trifolium86MainK5ERpc5']} vs {d['Trifolium86MainK5PC5']}" for d in sib))
lf_runs = [r for r in runs if r["cfg"] == CFG_LFMM]
rep("LFMM sibling runs: " + "; ".join(f"{r['log'][:17]} mode={r['mode']} MD5 {'=' if r['md5'] == r['md5_today'] else '!='} today's file"
                                       for r in lf_runs))

scr = {(r["quantity"], r["predictor_1"], r["predictor_2"]): float(r["value"]) for r in read_tsv(f"{OUT}/predictor_screen.tsv")}
rda_cfg = next(c for c in data["GEA"]["configs"] if c["method"] == "RDA")
assert rda_cfg["params"]["predictor_set"] == "bio_1,bio_3,bio_5" == data["Maladaptation"]["predictors"]
assert all(abs(v) < 0.8 for (q, _, _), v in scr.items() if q == "pairwise r") and \
       all(v < 10 for (q, _, _), v in scr.items() if q == "VIF")
with open(f"{ROOT}/ANALYSIS/Trifolium_rerun_v3/predictor_selection_r08/predictor_set_selected.txt") as f:
    assert f.read().strip() == "bio_1,bio_3,bio_5"
out.append(dict(step="GEA", setting="RDA predictor screen (GEA.configs[RDA].params.predictor_set, Maladaptation.predictors)",
                value="bio_1,bio_3,bio_5", scope="outside this configuration",
                note="site-level |r| < 0.8 and VIF < 10 over 16 candidates (BIO14 invariant, BIO17 and BIO18 excluded); "
                     "the largest passing subsets have 3 predictors, 6 of them tied within 0.01 adjusted R², "
                     "lowest-index rule (ANALYSIS/Trifolium_rerun_v3/compare/select_predictors_cut.R)"))

n_sel = sum(1 for _ in open(f"{RES}/_intermediate/snp_sets/gea.ER_overlap/selected_snps.tsv")) - 1
pool = sum(1 for _ in open(f"{RES}/_work/maf0.1_miss0.3_smiss0.8_rel0.999/Trifolium86_chr1-6_K5imp_maf0.1.vcfsnp"))
assert n_sel == 30 and pool == 177259 == int(summ[("gea", "snps_after_imputed_maf")])
SM = f"{AN}/fig3/offset_minou/site_maps"
gfd = [r for r in read_tsv(f"{SM}/gf_random_null/site_offsets_draws.tsv") if r["kind"] == "random"]
gseeds = sorted({int(r["seed"]) for r in gfd})
assert gseeds == list(range(100001, 100200))
grd = read_tsv(f"{SM}/random_null/site_offsets_draws.tsv")
for e in ("geometric", "rda"):
    assert sorted({int(r["draw"]) for r in grd if r["engine"] == e}) == list(range(1, 200))
out.append(dict(step="genetic offset", setting="random-set offset test", value="199 sets of 30 SNPs",
                scope="outside this configuration",
                note="sets drawn from the 177,259 tested SNPs, one seed per set (set.seed(100000 + i), i = 1–199), "
                     "the same sets for the three methods (the random-set scripts of manuscript/INDEX.tsv)"))

# ---- 7. defaults: keys the app's template sets and the file does not -------------------------------------------------
DEF_COMMIT = "309691b"   # first commit holding the code the ERpc5 runs used (Maladaptation.predictors, uncommitted at the run)
# extracted on the host: git show 309691b:scripts/clinego.app/inst/config_default.yaml > <this file>
with open(f"{ROOT}/work/manuscript_supp_2026-10-25_work/scratch/config_default_{DEF_COMMIT}.yaml") as f:   # Phase 7: scratch moved out of the deliverable
    dt = flat_data(yaml.safe_load(f))
# a method entry the file does not list is not a default (the method list is the setting); Traits = mode=traits
methods = {c["method"] for c in data["GEA"]["configs"]}
miss = sorted(k for k in dt if k not in fa and not any(k.startswith(pre) for pre, _ in LEFT_OUT + [("Traits", "")])
              and not (k.startswith("GEA.configs[") and k[len("GEA.configs["):].split("]")[0] not in methods))
rep(f"defaults ({DEF_COMMIT}:scripts/clinego.app/inst/config_default.yaml, the app's template; the pipeline's own "
    f"fallbacks are the config.get() defaults in workflow/rules/common.smk) - keys of run groups the file does not set "
    f"({len(miss)}): " + "; ".join(f"{k} = {json.dumps(dt[k])}" for k in miss))

# ---- 8. write ---------------------------------------------------------------------------------------------------------
ORDER = ["data processing", "population structure", "environmental factors", "GEA", "candidate regions", "genetic offset"]
out.sort(key=lambda r: ORDER.index(r["step"]))            # stable: file order within a step, outside rows last
cols = ["step", "setting", "value", "scope", "note"]
with open(f"{OUT}/t_settings.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(out)
esc = lambda s: s.replace("|", "\\|")
md = ["# t_settings", "",
      "One row per setting the run used, by step. Settings and values as `config_Trifolium86MainK5ERpc5.yaml` writes "
      "them; rows marked “outside this configuration” were set elsewhere (note). Settings of modes the run "
      "never ran (climate, traits, pregea, gwas, gea_x_gwas) are left out, and so is `Enrichment`, which no pipeline rule reads. Processing, population structure and the environmental-factor extraction ran on 2026-10-14 under `config_Trifolium86MainK5ER.yaml`, (proven by Snakemake's logged config hash), which differs from this file only in `project_name`, `cpu` (90) and GEA and offset keys those steps do not read. Values: `t_settings.tsv`; built by "
      "`fig3/table_settings.py`.", "",
      "| Step | Setting | Value | Scope | Note |", "| --- | --- | --- | --- | --- |"]
md += [f"| {r['step']} | `{esc(r['setting'])}` | `{esc(r['value'])}` | {r['scope']} | {esc(r['note'])} |" for r in out]
with open(f"{OUT}/t_settings.md", "w") as f:
    f.write("\n".join(md) + "\n")
print(f"rows: {len(out)} ({sum(r['scope'] == 'this configuration' for r in out)} from the file, "
      f"{sum(r['scope'] != 'this configuration' for r in out)} outside)")
for s in ORDER:
    print(f"  {s}: {sum(r['step'] == s for r in out)}")
print("done")
