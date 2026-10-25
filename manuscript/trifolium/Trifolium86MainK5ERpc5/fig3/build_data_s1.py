#!/usr/bin/env python3
"""data_s1 (supplement build 2026-10-25, D1): the configuration files of the reported analyses, with no comments.

Every shipped YAML is a comment-free copy of the file a reported run read. Each copy is checked three ways against its
original: yaml.safe_load gives the same data, yaml.scan gives the same token stream (type, value, style; PyYAML's
scanner drops comments, so equal streams mean only comments and blank lines were removed), and no '#' is left outside
a quoted scalar. The originals are only read.

Which file a run read is not taken on trust: Snakemake logs `Config MD5` = md5(json.dumps(config + {"mode": m},
sort_keys=True)) for every run (snakemake/workflow.py), and every file shipped here reproduces the MD5 of its runs.

Families of the simulation (600 replicates = the `primary_ssclines` seeds of benchmarks/mvp_seeds.tsv):
  gea_scan         config_MVP{seed}_c1.yaml: processing, prestructure, structure, gea (LFMM, EMMAX, RDA with
                   structure correction). Varies only in per-replicate fields -> one template + replicates table.
  rda_uncorrected  the RDA scan without structure correction never read a YAML: benchmarks/mvp_rda_uncorrected_run.sh
                   re-issued each replicate's recorded scripts/rda.R command with condition_pcs 0. Shipped as the
                   gea_scan template with that one key changed.
  offset           config_MVP{seed}_sweep.yaml: maladaptation over the garden scenarios. Its SNP-set roster differs
                   between replicates, which is not a per-replicate field of the template rule -> all 600 files.
The base family config_MVP{seed}.yaml is not shipped: no reported run read it (upstream ran under _c1).

Writes work/manuscript_supp_2026-10-25/data_s1/ (refuses to write into an existing folder). Run from the repo root:
  python3 -I ANALYSIS/Trifolium86MainK5ERpc5/fig3/build_data_s1.py
"""
import collections
import csv
import glob
import hashlib
import json
import os
import re
import sys

import yaml

ROOT = os.getcwd()
assert os.path.isfile(os.path.join(ROOT, "Snakefile")), "run from the repo root"
OUT = "work/manuscript_supp_2026-10-25/data_s1"
TRI_CFG = {
    "config_Trifolium86MainK5ERpc5.yaml": "config_Trifolium86MainK5ERpc5.yaml",
    "config_Trifolium86MainK5PC5.yaml": "config_Trifolium86MainK5PC5.yaml",
    "config_Trifolium86MainKsearch.yaml": "_archive/trifolium_pre_2026-10-17/config_Trifolium86MainKsearch.yaml",
}
K_EXPECTED = {2: 200, 3: 272, 4: 74, 5: 50, 7: 3, 8: 1}
K_LINES = ("k_best", "n_pcs", "K", "condition_pcs")
MD5_RE = re.compile(r"--config mode=(\S+) --configfile (\S+).*?\n(?:.*\n){0,8}?\s*Config MD5: ([0-9a-f]{32})")


def rep(s):
    print("REPORT:", s, flush=True)


def load(path):
    with open(path) as f:
        return yaml.safe_load(f)


def config_md5(data, mode):
    c = dict(data)
    c["mode"] = mode
    return hashlib.md5(json.dumps(c, sort_keys=True).encode("utf-8")).hexdigest()


def flat(d, p=""):
    """Flatten parsed data; a list of mappings is keyed by position."""
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, f"{p}.{k}" if p else str(k)))
    elif isinstance(d, list) and d and all(isinstance(x, dict) for x in d):
        for i, x in enumerate(d):
            out.update(flat(x, f"{p}[{i}]"))
    else:
        out[p] = d
    return out


# ---- comment stripping -----------------------------------------------------------------------------------------------
def strip_line(line):
    """Cut a line at its comment. A quote opens a quoted scalar only where a token can start."""
    q = None
    last = ""                      # last non-blank character outside a quoted scalar
    i = 0
    while i < len(line):
        ch = line[i]
        if q is None:
            if ch == "#" and (i == 0 or line[i - 1] in " \t"):
                break
            if ch in "\"'" and last in ("", ":", "-", "[", "{", ",", "?"):
                q = ch
            elif not ch.isspace():
                last = ch
        elif q == "'":
            if ch == "'":
                if i + 1 < len(line) and line[i + 1] == "'":
                    i += 1
                else:
                    q, last = None, "'"
        else:
            if ch == "\\":
                i += 1
            elif ch == '"':
                q, last = None, '"'
        i += 1
    assert q is None, f"quoted scalar runs past the end of the line: {line!r}"
    return line[:i].rstrip()


def strip_comments(text):
    assert not re.search(r":\s*[|>][-+0-9]*\s*(#.*)?$", text, re.M), "block scalar present; stripper is line-based"
    lines = [s for s in (strip_line(l) for l in text.splitlines()) if s.strip()]
    out = []
    for s in lines:
        if out and not s[0].isspace() and s[0] != "-":
            out.append("")           # one blank line before every top-level key, for reading only
        out.append(s)
    return "\n".join(out) + "\n"


def tokens(text):
    return [(type(t).__name__, getattr(t, "value", None), getattr(t, "style", None)) for t in yaml.scan(text)]


def check_copy(orig_text, copy_text, what):
    assert yaml.safe_load(copy_text) == yaml.safe_load(orig_text), f"{what}: parsed data differ"
    assert tokens(copy_text) == tokens(orig_text), f"{what}: token streams differ"
    for l in copy_text.splitlines():
        assert strip_line(l) == l, f"{what}: a comment is left: {l!r}"


os.makedirs(os.path.dirname(OUT), exist_ok=True)
assert not os.path.exists(OUT), f"{OUT} exists; never write into an existing output folder"
os.makedirs(f"{OUT}/trifolium")
os.makedirs(f"{OUT}/simulation/offset")
written = []


def write(rel, text):
    with open(f"{OUT}/{rel}", "w", newline="\n") as f:
        f.write(text)
    written.append(rel)


# ---- 1. Trifolium: the three configs and the runs that read them ------------------------------------------------------
runs = collections.defaultdict(list)          # config file name -> [(log, mode, logged md5)]
for f in sorted(glob.glob(".snakemake/log/*.snakemake.log")):
    with open(f, errors="replace") as fh:
        head = "".join(next(fh, "") for _ in range(60))
    for mode, cfg, md5 in MD5_RE.findall(head):
        if cfg in TRI_CFG:
            runs[cfg].append((os.path.basename(f), mode, md5))
# the K sweep's driver also kept each run's stdout; those are the same runs, logged twice -> must match, not count
for f in sorted(glob.glob("ANALYSIS/Trifolium86MainKsearch/logs/*.log")):
    with open(f, errors="replace") as fh:
        for mode, cfg, md5 in MD5_RE.findall(fh.read()):
            assert (mode, md5) in {(m, h) for _, m, h in runs[cfg]}, f"{f}: a run missing from .snakemake/log"

for name, src in TRI_CFG.items():
    with open(src) as f:
        orig = f.read()
    data = yaml.safe_load(orig)
    rs = runs[name]
    assert rs, f"no logged run read {name}"
    for log, mode, md5 in rs:
        assert config_md5(data, mode) == md5, f"{name}: {log} mode={mode} logged MD5 is not today's file"
    rep(f"{name} ({src}): Config MD5 of {len(rs)} logged runs reproduced, modes "
        + ", ".join(sorted({m for _, m, _ in rs})))
    copy = strip_comments(orig)
    check_copy(orig, copy, name)
    write(f"trifolium/{name}", copy)
    rep(f"trifolium/{name}: {len(orig.splitlines())} -> {len(copy.splitlines())} lines, data and tokens equal")
assert {m for _, m, _ in runs["config_Trifolium86MainK5ERpc5.yaml"]} == {"gea", "maladaptation"}
assert {m for _, m, _ in runs["config_Trifolium86MainKsearch.yaml"]} == {"processing", "prestructure"}
ks = load(TRI_CFG["config_Trifolium86MainKsearch.yaml"])["sNMF"]
assert (ks["k_start"], ks["k_end"], ks["repeats"]) == (2, 8, 10), ks

# ---- 2. simulation: the 600 replicates and their K --------------------------------------------------------------------
with open("benchmarks/mvp_seeds.tsv", newline="") as f:
    MAN = [r for r in csv.DictReader(f, delimiter="\t") if r["arm"] == "primary_ssclines"]
SEEDS = sorted(r["seed"] for r in MAN)
assert len(SEEDS) == 600 == len(set(SEEDS))
KB = {r["seed"]: int(r["k_best"]) for r in MAN}
for r in MAN:
    assert KB[r["seed"]] == max(int(r["K_authors"]), 2), r["seed"]

C1 = {s: f"config_MVP{s}_c1.yaml" for s in SEEDS}
SW = {s: f"config_MVP{s}_sweep.yaml" for s in SEEDS}
TXT_C1 = {s: open(C1[s]).read() for s in SEEDS}
D_C1 = {s: yaml.safe_load(TXT_C1[s]) for s in SEEDS}
D_SW = {s: load(SW[s]) for s in SEEDS}

k_cfg = collections.Counter()
for s in SEEDS:
    d = D_C1[s]
    by = {c["method"]: c for c in d["GEA"]["configs"]}
    assert sorted(by) == ["EMMAX", "LFMM", "RDA"], (s, sorted(by))
    k = d["sNMF"]["k_best"]
    assert k == KB[s] == by["EMMAX"]["params"]["n_pcs"] == by["LFMM"]["params"]["K"] == \
        by["RDA"]["params"]["condition_pcs"] == D_SW[s]["sNMF"]["k_best"], s
    mm = D_SW[s]["Maladaptation"]["methods"]
    assert mm["rda_offset"]["condition_pcs"] == 0 and mm["rda_offset_corrected"]["condition_pcs"] == k, s
    k_cfg[k] += 1
assert dict(k_cfg) == K_EXPECTED, dict(k_cfg)
rep("K (sNMF.k_best = EMMAX n_pcs = LFMM K = RDA condition_pcs, = max(K_authors, 2) of mvp_seeds.tsv): "
    + ", ".join(f"K = {k} in {n}" for k, n in sorted(k_cfg.items())))

# ---- 3. which runs read them: Config MD5 of every logged run ----------------------------------------------------------
def md5_runs(logdirs, cfg_of, data_of):
    seen = collections.Counter()
    for s in SEEDS:
        hits = []
        for d in logdirs:
            p = f"{d}/MVP{s}.log"
            if os.path.exists(p):
                with open(p, errors="replace") as fh:
                    hits += MD5_RE.findall(fh.read())
        assert hits, f"no logged run for {s}"
        for mode, cfg, md5 in hits:
            assert cfg == cfg_of[s], (s, cfg)
            assert config_md5(data_of[s], mode) == md5, f"{s} mode={mode}: logged MD5 is not today's {cfg}"
            seen[mode] += 1
    return seen

seen = md5_runs(["benchmarks/mvp_eval/runlogs_onefit", "benchmarks/mvp_eval/onefit_pilot/runlogs"], C1, D_C1)
assert dict(seen) == {m: 600 for m in ("processing", "prestructure", "structure", "gea")}, dict(seen)
rep("gea_scan: today's 600 _c1 files reproduce the Config MD5 of every logged run "
    "(benchmarks/mvp_eval/runlogs_onefit 597 + onefit_pilot/runlogs 3; " + ", ".join(f"{m} {n}" for m, n in seen.items())
    + "); mode=pregea was not run")
seen = md5_runs(["benchmarks/mvp_eval/offset13/sweeplogs"], SW, D_SW)
assert dict(seen) == {"maladaptation": 600}, dict(seen)
rep("offset: today's 600 _sweep files reproduce the Config MD5 of benchmarks/mvp_eval/offset13/sweeplogs "
    "(mode=maladaptation, 600)")

# _sweep was written from the BASE config, not from _c1: list where the two differ
d_keys = collections.Counter()
for s in SEEDS:
    a, b = flat(D_C1[s]), flat(D_SW[s])
    for k in set(a) | set(b):
        if a.get(k, "<absent>") != b.get(k, "<absent>"):
            d_keys[re.sub(r"\.\w+$", "", k) if k.startswith("Maladaptation.") else k] += 1
allowed = re.compile(r"^(GEA\.configs\[\d\]\.params\.\w+|LDdecay\.scope|PreGEA\.Advanced\.permutations|"
                     r"Climate\.custom\.future_tables|cpu|Maladaptation(\..*)?)$")
assert all(allowed.match(k) for k in d_keys), sorted(d_keys)
rep("offset vs gea_scan, the keys that differ (replicates): " + "; ".join(f"{k} {n}" for k, n in sorted(d_keys.items())))

# ---- 4. gea_scan: one template, a replicates table, regenerated and compared --------------------------------------------
def templ(s, text):
    t = text.replace(s, "{seed}")
    k = str(KB[s])
    n = 0
    lines = []
    for l in t.splitlines():
        m = re.match(r"^(\s*(?:-\s+)?(" + "|".join(K_LINES) + r"):\s+)(\S+)$", l)
        if m and m.group(3) == k:
            l = m.group(1) + "{k_best}"
            n += 1
        lines.append(l)
    assert n == 4, (s, n)
    return "\n".join(lines) + "\n"

def render(t, s):
    return t.replace("{seed}", s).replace("{k_best}", str(KB[s]))

STRIP_C1 = {s: strip_comments(TXT_C1[s]) for s in SEEDS}
for s in SEEDS:
    check_copy(TXT_C1[s], STRIP_C1[s], C1[s])
TEMPLATES = {templ(s, STRIP_C1[s]) for s in SEEDS}
assert len(TEMPLATES) == 1, f"{len(TEMPLATES)} distinct gea_scan templates"
T_C1 = TEMPLATES.pop()
assert "{" not in T_C1.replace("{seed}", "").replace("{k_best}", ""), "a brace other than the placeholders"
for s in SEEDS:
    assert render(T_C1, s) == STRIP_C1[s], s
    assert yaml.safe_load(render(T_C1, s)) == D_C1[s], s
rep("gea_scan: the 600 comment-free _c1 texts are byte-identical once {seed} and {k_best} are substituted; "
    "the rendered template parses to each file's data (600/600)")

def write_replicates(rel):
    with open(f"{OUT}/{rel}", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["seed", "k_best"])
        for s in SEEDS:
            w.writerow([s, KB[s]])
    written.append(rel)

write("simulation/gea_scan_template.yaml", T_C1)
write_replicates("simulation/gea_scan_replicates.tsv")

# ---- 5. rda_uncorrected: the gea_scan template with RDA condition_pcs 0 -------------------------------------------------
cp = [i for i, l in enumerate(T_C1.splitlines()) if re.match(r"^\s*condition_pcs:\s+\{k_best\}$", l)]
assert len(cp) == 1, cp
L = T_C1.splitlines()
L[cp[0]] = L[cp[0]].replace("{k_best}", "0")
T_UNC = "\n".join(L) + "\n"
for s in SEEDS:
    a, b = flat(D_C1[s]), flat(yaml.safe_load(render(T_UNC, s)))
    diff = [k for k in set(a) | set(b) if a.get(k, "<absent>") != b.get(k, "<absent>")]
    rda = [i for i, c in enumerate(D_C1[s]["GEA"]["configs"]) if c["method"] == "RDA"][0]
    assert diff == [f"GEA.configs[{rda}].params.condition_pcs"], (s, diff)
    assert b[diff[0]] == 0
    with open(f"benchmarks/mvp_eval/params_rdaunc/MVP{s}/c1/DONE") as f:
        done = dict(l.rstrip("\n").split("\t", 1) for l in f if "\t" in l)
    assert done["seed"] == s and int(done["k_best"]) == KB[s] and done["condition_pcs"] == "0", (s, done)
rep("rda_uncorrected: differs from gea_scan only at GEA.configs[RDA].params.condition_pcs (k_best -> 0) in 600/600; "
    "all 600 benchmarks/mvp_eval/params_rdaunc/MVP*/c1/DONE record condition_pcs 0 and the manifest's k_best")
write("simulation/rda_uncorrected_template.yaml", T_UNC)
write_replicates("simulation/rda_uncorrected_replicates.tsv")

# ---- 6. offset: all 600 files (the SNP-set roster is not a per-replicate field of the template rule) ---------------------
rosters = collections.Counter()
for s in SEEDS:
    with open(SW[s]) as f:
        orig = f.read()
    copy = strip_comments(orig)
    check_copy(orig, copy, SW[s])
    write(f"simulation/offset/{SW[s]}", copy)
    rosters[tuple(D_SW[s]["Maladaptation"]["snp_sets"])] += 1
for r, n in rosters.most_common():
    rep(f"offset roster in {n} replicates: {', '.join(r)}")
assert sorted(rosters.values()) == [2, 118, 480], rosters
ft = {json.dumps(D_SW[s]["Climate"]["custom"]["future_tables"]) for s in SEEDS}
assert len(ft) == 1
n_gardens = len(json.loads(ft.pop()))
mm = D_SW[SEEDS[0]]["Maladaptation"]["methods"]
rep(f"offset: {n_gardens} garden scenarios (the same list in every replicate); methods {', '.join(mm)}")

# ---- 7. manifest ----------------------------------------------------------------------------------------------------
MANIFEST = [
    ("trifolium/config_Trifolium86MainK5ERpc5.yaml", "Trifolium: GEA (EMMAX, RDA) and genetic offset",
     "the run Trifolium86MainK5ERpc5, modes gea and maladaptation (Gradient Forest, geometric offset, RDA-based "
     "offset). Its processing, prestructure and structure outputs were made 2026-10-14 under "
     "config_Trifolium86MainK5ER.yaml, which as those runs read it differs from this file only at project_name, cpu, "
     "EMMAX n_pcs, RDA condition_pcs, GEA.snp_sets, Maladaptation.predictors and "
     "Maladaptation.methods.rda_offset.axes"),
    ("trifolium/config_Trifolium86MainK5PC5.yaml", "Trifolium: LFMM",
     "the sibling run Trifolium86MainK5PC5; the reported analysis is its LFMM scan (K 5)"),
    ("trifolium/config_Trifolium86MainKsearch.yaml", "Trifolium: choice of K",
     "the sNMF sweep over K 2–8 with 10 repeats per K (modes processing and prestructure; cross-entropy)"),
    ("simulation/gea_scan_template.yaml", "simulation: GEA scan",
     "modes processing, prestructure, structure and gea (LFMM, EMMAX, RDA with structure correction) of each of the "
     "600 replicates; replace {seed} and {k_best} with a row of gea_scan_replicates.tsv before parsing. It configures "
     "the p-value scans only: the reported candidate sets were called afterwards from the harvested p-value tables "
     "by benchmarks/mvp_build_snp_sets.R (top 0.25% of the tested SNPs per method and predictor, uncorrected RDA "
     "at 0.5%, agreement window 5 kb; benchmarks/mvp_eval/offset13/provenance.tsv), not by this file's adjust, "
     "threshold or snp_clumping_distance"),
    ("simulation/gea_scan_replicates.tsv", "simulation: GEA scan",
     "one row per replicate: seed and k_best (the depositors' K, at least 2) for the template"),
    ("simulation/rda_uncorrected_template.yaml", "simulation: RDA without structure correction",
     "settings of the direct `scripts/rda.R` re-run by `benchmarks/mvp_rda_uncorrected_run.sh`; this arm never read "
     "a YAML. The gea_scan template with RDA condition_pcs 0; fill it from rda_uncorrected_replicates.tsv"),
    ("simulation/rda_uncorrected_replicates.tsv", "simulation: RDA without structure correction",
     "one row per replicate: seed and k_best for the template (equal to gea_scan_replicates.tsv)"),
    ("simulation/offset/config_MVP{seed}_sweep.yaml (600 files)", "simulation: genetic offset",
     f"mode maladaptation of each replicate over {n_gardens} garden scenarios: Gradient Forest, geometric offset, "
     "RDA-based offset (rda_offset, without structure correction) and rda_offset_corrected (run, not reported); "
     "snp_sets lists the SNP sets of that replicate, built by benchmarks/mvp_build_snp_sets.R (truth = its causal "
     "loci; union_top* = size-curve sets, run, not reported)"),
]
with open(f"{OUT}/manifest.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["file", "analysis", "what it configures"])
    w.writerows(MANIFEST)
written.append("manifest.tsv")

# ---- 8. the folder holds exactly what was written --------------------------------------------------------------------
on_disk = sorted(os.path.relpath(os.path.join(d, f), OUT) for d, _, fs in os.walk(OUT) for f in fs)
assert on_disk == sorted(written), set(on_disk) ^ set(written)
assert len(written) == 3 + 4 + 600 + 1
for row in MANIFEST:
    p = row[0].split(" (")[0]
    assert glob.glob(f"{OUT}/" + p.replace("{seed}", "*")), p
banned = re.compile(r"risk|adaptive|null|Climate-only|>=", re.I)
with open(f"{OUT}/manifest.tsv") as f:
    assert not banned.search(f.read()), "banned word in manifest.tsv"
rep(f"data_s1: {len(written)} files in {OUT}")
print("done")
