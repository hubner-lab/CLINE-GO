#!/usr/bin/env python3
"""mvp_land_data_s1.py -- Job A of docs/manuscript/lab/data_s1_landscapes_LAB_PROMPT.md (2026-11-08): Data S1's
configuration files of the SS-Mtn and Est-Clines simulation runs, built like the SS-Clines part
(manuscript/trifolium/Trifolium86MainK5ERpc5/fig3/build_data_s1.py, lab copy work/ld_window_lab/scratch_supp/
build_data_s1_d150.py), whose comment stripping and copy checks are copied here unchanged.

The two landscapes were not run in SS-Clines' three families (user, 2026-11-08: ship the families that ran). Every
replicate ran three, each read by an executed run in .snakemake/log (the authority; -n runs are listed, never counted):
  sets           config_MVP{seed}hd.yaml: processing, prestructure, structure and gea from the deposit's files
                 (LFMM, EMMAX, RDA without structure correction, the LD decay, the six declared candidate sets)
  rda_corrected  config_MVP{seed}hdrc.yaml: gea only, RDA conditioned on K PCs, on the hd tree's upstream outputs
                 (benchmarks/mvp_clone_project.sh hardlinks them; its candidate sets and offsets were not run)
  offset         config_MVP{seed}hd_sweep.yaml: maladaptation over the garden scenarios, all 600 files
The base config_MVP{seed}.yaml and _c1 files are not shipped: no reported run read them (generator parents).

Per-replicate fields of the templates: {seed}, {k_best} (sNMF k_best = EMMAX n_pcs = LFMM K; hdrc RDA condition_pcs),
{k_end} (10, or 11 where the depositors' K is 11) and {clumping_distance} (auto_genome_wide, or 619 bp declared for the
39 Est-Clines replicates whose LD curve was unreliable). The replicates tables carry those plus the LD-decay distance as
used (the clumping_distance of the six snp_set.yaml files the hd run wrote, cross-checked against
work/mvp_land/ld_window_by_replicate.tsv and, for an auto distance, the tree's own ld_decay_half_distances.tsv).

SS-Clines is checked, not rewritten: the published simulation/ files (work/ld_window_lab/results/phase8/data_s1, the
bundle the laptop hash-checked into its pack) must match their sha256 list, and today's SS-Clines configurations
(_archive/mvp_ssclines/configs/) must render to them byte for byte.

Writes OUT/simulation_ssmtn/, OUT/simulation_estclines/ and OUT/manifest.tsv (refuses existing folders) and
OUT/checks/data_s1_check.tsv. Run from the repo root:  python3 -I benchmarks/mvp_land_data_s1.py
"""
import collections
import csv
import glob
import hashlib
import json
import os
import re

import yaml

ROOT = os.getcwd()
assert os.path.isfile(os.path.join(ROOT, "Snakefile")), "run from the repo root"
OUT = "benchmarks/mvp_eval/landscape_ms/data_s1_2026-11"
PUB = "work/ld_window_lab/results/phase8/data_s1"
ARCH = "_archive/mvp_ssclines/configs"
LANDS = {"ssmtn": "SS-Mtn", "estclines": "Est-Clines"}
K_LINES = ("k_best", "n_pcs", "K", "condition_pcs")
CMD_RE = re.compile(r"Command: (.*?)\n(?:.*\n){0,8}?\s*Config MD5: ([0-9a-f]{32})")
CHECKS = []


def check(item, what, ok, detail=""):
    CHECKS.append((item, what, "OK" if ok else "FAIL", detail))
    print(f"{'OK  ' if ok else 'FAIL'} {item}: {what} {detail}", flush=True)
    assert ok, f"{item}: {what} {detail}"


def load(path):
    with open(path) as f:
        return yaml.safe_load(f)


def config_md5(data, mode):
    c = dict(data)
    c["mode"] = mode
    return hashlib.md5(json.dumps(c, sort_keys=True).encode("utf-8")).hexdigest()


def flat(d, p=""):
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


def diff_keys(a, b):
    a, b = flat(a), flat(b)
    return sorted(k for k in set(a) | set(b) if a.get(k, "<absent>") != b.get(k, "<absent>"))


# ---- comment stripping and copy checks: verbatim from build_data_s1.py -------------------------------------------------
def strip_line(line):
    q = None
    last = ""
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
            out.append("")
        out.append(s)
    return "\n".join(out) + "\n"


def tokens(text):
    return [(type(t).__name__, getattr(t, "value", None), getattr(t, "style", None)) for t in yaml.scan(text)]


def check_copy(orig_text, copy_text, what):
    assert yaml.safe_load(copy_text) == yaml.safe_load(orig_text), f"{what}: parsed data differ"
    assert tokens(copy_text) == tokens(orig_text), f"{what}: token streams differ"
    for l in copy_text.splitlines():
        assert strip_line(l) == l, f"{what}: a comment is left: {l!r}"


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# ---- templates --------------------------------------------------------------------------------------------------------
def templ(s, text, rep, n_k, n_cd):
    """Replace the replicate's fields by placeholders; n_k K lines, one k_end line, n_cd clumping-distance lines."""
    t = text.replace(s, "{seed}")
    k, kend, cd = str(rep["k_best"]), str(rep["k_end"]), str(rep["clumping_distance"])
    nk = ne = nc = 0
    lines = []
    for l in t.splitlines():
        m = re.match(r"^(\s*(?:-\s+)?(" + "|".join(K_LINES) + r"):\s+)(\S+)$", l)
        e = re.match(r"^(\s*k_end:\s+)(\S+)$", l)
        c = re.match(r"^(\s*(?:-\s+)?(?:snp_)?clumping_distance:\s+)(\S+)$", l)
        if m and m.group(3) == k:
            l = m.group(1) + "{k_best}"
            nk += 1
        elif e and e.group(2) == kend:
            l = e.group(1) + "{k_end}"
            ne += 1
        elif c and c.group(2) == cd:
            l = c.group(1) + "{clumping_distance}"
            nc += 1
        lines.append(l)
    assert (nk, ne, nc) == (n_k, 1, n_cd), (s, nk, ne, nc)
    return "\n".join(lines) + "\n"


def render(t, s, rep):
    return (t.replace("{seed}", s).replace("{k_best}", str(rep["k_best"])).replace("{k_end}", str(rep["k_end"]))
            .replace("{clumping_distance}", str(rep["clumping_distance"])))


# ---- 0. folders ---------------------------------------------------------------------------------------------------------
for land in LANDS:
    assert not os.path.exists(f"{OUT}/simulation_{land}"), f"{OUT}/simulation_{land} exists; never write into it"
assert not os.path.exists(f"{OUT}/manifest.tsv"), "manifest.tsv exists"
os.makedirs(f"{OUT}/checks", exist_ok=True)
written = []


def write(rel, text):
    os.makedirs(os.path.dirname(f"{OUT}/{rel}"), exist_ok=True)
    with open(f"{OUT}/{rel}", "w", newline="\n") as f:
        f.write(text)
    written.append(rel)


# ---- 1. seeds, manifest, LD table ---------------------------------------------------------------------------------------
with open("benchmarks/mvp_seeds.tsv", newline="") as f:
    MAN = {r["seed"]: r for r in csv.DictReader(f, delimiter="\t")}
SEEDS = {l: open(f"work/mvp_land/seeds_{l}.txt").read().split() for l in LANDS}
for l, s in SEEDS.items():
    check(LANDS[l], "600 distinct seeds, all in benchmarks/mvp_seeds.tsv", len(s) == 600 == len(set(s)) and all(x in MAN for x in s))
check("seeds", "no seed in both landscapes", not set(SEEDS["ssmtn"]) & set(SEEDS["estclines"]))
with open("work/mvp_land/ld_window_by_replicate.tsv", newline="") as f:
    LDW = {(r["landscape"], r["seed"]): r for r in csv.DictReader(f, delimiter="\t")}

# ---- 2. every logged run of the three families --------------------------------------------------------------------------
ALL = {s: l for l, ss in SEEDS.items() for s in ss}
RUNS = collections.defaultdict(list)          # (seed, family) -> [(log, mode, kind, md5)]
for f in sorted(glob.glob(".snakemake/log/*.snakemake.log")):
    with open(f, errors="replace") as fh:
        head = "".join(next(fh, "") for _ in range(60))
    for cmd, h in CMD_RE.findall(head):
        m = re.search(r"--configfile config_MVP(\d+)(hd|hdrc|hd_sweep)\.yaml", cmd)
        if not m or m.group(1) not in ALL:
            continue
        kind = "dry" if re.search(r"(^| )-n( |$)", cmd) else "touch" if "--touch" in cmd else \
            "omit-from" if "--omit-from" in cmd else "run"
        RUNS[(m.group(1), m.group(2))].append((os.path.basename(f), re.search(r"mode=(\S+)", cmd).group(1), kind, h))

FAM = {"hd": "sets", "hdrc": "rda_corrected", "hd_sweep": "offset"}
MODES = {"hd": ["gea", "prestructure", "processing", "structure"], "hdrc": ["gea"], "hd_sweep": ["maladaptation"]}
REPS, TXT, DATA = {}, {}, {}
for land, name in LANDS.items():
    REPS[land] = {}
    n_backup = 0
    for s in SEEDS[land]:
        for sfx in FAM:
            p = f"config_MVP{s}{sfx}.yaml"
            TXT[(s, sfx)] = open(p).read()
            DATA[(s, sfx)] = yaml.safe_load(TXT[(s, sfx)])
        hd, rc, sw = DATA[(s, "hd")], DATA[(s, "hdrc")], DATA[(s, "hd_sweep")]
        # ---- K, k_end, clumping distance, from the files
        by = {c["method"]: c for c in hd["GEA"]["configs"]}
        k = hd["sNMF"]["k_best"]
        assert sorted(by) == ["EMMAX", "LFMM", "RDA"], s
        assert k == by["EMMAX"]["params"]["n_pcs"] == by["LFMM"]["params"]["K"] and by["RDA"]["params"]["condition_pcs"] == 0, s
        rb = {c["method"]: c for c in rc["GEA"]["configs"]}
        assert sorted(rb) == ["RDA"] and rc["sNMF"]["k_best"] == k == rb["RDA"]["params"]["condition_pcs"], s
        mm = sw["Maladaptation"]["methods"]
        assert sw["sNMF"]["k_best"] == k and mm["rda_offset"]["condition_pcs"] == 0 and mm["rda_offset_corrected"]["condition_pcs"] == k, s
        ka = int(MAN[s]["K_authors"])
        assert k == max(ka, 2) == int(MAN[s]["k_best"]), (s, k, ka)
        kend = hd["sNMF"]["k_end"]
        assert kend == rc["sNMF"]["k_end"] == sw["sNMF"]["k_end"] and (kend == 11) == (k > 10) and kend in (10, 11), s
        cd = hd["GEA"]["snp_clumping_distance"]
        assert cd in ("auto_genome_wide", 619), (s, cd)
        assert all(x["clumping_distance"] == cd for x in hd["GEA"]["snp_sets"]) and rc["GEA"]["snp_clumping_distance"] == cd \
            and sw["GEA"]["snp_clumping_distance"] == cd, s
        declared = cd == 619
        # ---- the distance as used
        used = set()
        sets = sorted(glob.glob(f"MVP{s}hd_results/_intermediate/snp_sets/gea.*/snp_set.yaml"))
        assert len(sets) == 6, (s, sets)
        for y in sets:
            d = load(y)
            used.add(int(d["clumping_distance"]))
            assert d["k_best"] == k, (s, y)
        assert len(used) == 1, (s, used)
        used = used.pop()
        w = LDW[(land, s)]
        assert int(w["clumping_distance_bp"]) == used, (s, w["clumping_distance_bp"], used)
        assert w["distance_source"].startswith("declared") == declared, (s, w["distance_source"])
        if declared:
            assert used == 619
        else:
            with open(f"MVP{s}hd_results/Structure/tables/ld_decay_half_distances.tsv", newline="") as f:
                ga = [r for r in csv.DictReader(f, delimiter="\t") if r["group"] == "All" and r["scope"] == "genome_wide"]
            assert len(ga) == 1 and int(float(ga[0]["decay_critical_bp"])) == used, (s, ga, used)
        with open(f"MVP{s}hd_results/Structure/tables/ld_decay_half_distances.tsv", "rb") as f1, \
                open(f"MVP{s}hdrc_results/Structure/tables/ld_decay_half_distances.tsv", "rb") as f2:
            assert f1.read() == f2.read(), f"{s}: hdrc's LD-decay table differs from hd's"   # hdrc's gea ran its own
        REPS[land][s] = dict(k_best=k, k_authors=ka, k_end=kend, clumping_distance=cd, ld_distance_bp=used,
                             distance_source=w["distance_source"], declared_619="yes" if declared else "no")
        # ---- which runs read them
        for sfx in FAM:
            rs = [r for r in RUNS[(s, sfx)] if r[2] != "dry"]
            assert sorted(r[1] for r in rs) == MODES[sfx] and all(r[2] == "run" for r in rs), (s, sfx, rs)
            for log, mode, kind, h in rs:
                if config_md5(DATA[(s, sfx)], mode) == h:
                    continue
                bk = f"work/mvp_land/config_backup_ld619/config_MVP{s}{sfx}.yaml"
                assert declared and sfx == "hd" and mode in ("processing", "prestructure", "structure") and \
                    os.path.exists(bk) and config_md5(load(bk), mode) == h, f"{s} {sfx} {mode} {log}: logged MD5 matches no file"
                dk = diff_keys(load(bk), hd)
                assert all(re.match(r"^GEA\.(snp_clumping_distance|snp_sets\[\d\]\.clumping_distance)$", x) for x in dk), (s, dk)
                n_backup += 1
    nd = sum(r["declared_619"] == "yes" for r in REPS[land].values())
    check(name, "K = sNMF k_best = EMMAX n_pcs = LFMM K = hdrc RDA condition_pcs = rda_offset_corrected condition_pcs "
          "= max(K_authors, 2) of mvp_seeds.tsv; hd RDA and rda_offset condition_pcs 0", True,
          "K " + ", ".join(f"{k}: {n}" for k, n in sorted(collections.Counter(r["k_best"] for r in REPS[land].values()).items()))
          + f"; floor of 2 applied in {sum(r['k_authors'] < 2 for r in REPS[land].values())}")
    check(name, "k_end 11 exactly where K > 10", True,
          f"{sum(r['k_end'] == 11 for r in REPS[land].values())} replicates at k_end 11")
    check(name, "distance as used = the six snp_set.yaml = ld_window_by_replicate.tsv; auto = the tree's critical-r2 "
          "distance; declared = 619", True, f"declared 619 bp: {nd}")
    check(name, "hdrc's gea run wrote its own LD-decay table (ld_decay_half_distances.tsv), byte-identical to hd's, so "
          "the rda_corrected replicates table repeats the sets row", True, "600/600")
    check(name, "every executed run of hd (4 modes), hdrc (gea) and hd_sweep (maladaptation) reproduces today's file's "
          "Config MD5; the exceptions are the declared replicates' upstream runs, which reproduce the pre-declaration "
          "copy (config_backup_ld619), differing only at the clumping-distance keys", True,
          f"600 x 6 runs; {n_backup} upstream runs on the pre-declaration copy; dry runs listed "
          f"{sum(1 for s in SEEDS[land] for sfx in FAM for r in RUNS[(s, sfx)] if r[2] == 'dry')}")
    if land == "estclines":
        check(name, "39 replicates declared at 619 bp", nd == 39)
    else:
        check(name, "no declared distance", nd == 0)

# ---- 3. templates, replicates tables, offset files ----------------------------------------------------------------------
COLS = ["seed", "k_best", "k_authors", "k_end", "clumping_distance", "ld_distance_bp", "distance_source", "declared_619"]
ROSTER = {}
for land, name in LANDS.items():
    D = f"simulation_{land}"
    for sfx, fam, n_k, n_cd in (("hd", "sets", 3, 7), ("hdrc", "rda_corrected", 2, 1)):
        strip = {}
        for s in SEEDS[land]:
            strip[s] = strip_comments(TXT[(s, sfx)])
            check_copy(TXT[(s, sfx)], strip[s], f"config_MVP{s}{sfx}.yaml")
        T = {templ(s, strip[s], REPS[land][s], n_k, n_cd) for s in SEEDS[land]}
        assert len(T) == 1, f"{name} {fam}: {len(T)} distinct templates"
        T = T.pop()
        assert "{" not in re.sub(r"\{(seed|k_best|k_end|clumping_distance)\}", "", T), "a brace other than the placeholders"
        for s in SEEDS[land]:
            assert render(T, s, REPS[land][s]) == strip[s], s
            assert yaml.safe_load(render(T, s, REPS[land][s])) == DATA[(s, sfx)], s
        write(f"{D}/{fam}_template.yaml", T)
        with open(f"{OUT}/{D}/{fam}_replicates.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(COLS)
            for s in SEEDS[land]:
                w.writerow([s] + [REPS[land][s][c] for c in COLS[1:]])
        written.append(f"{D}/{fam}_replicates.tsv")
        check(name, f"{fam}: one template; the 600 comment-free config_MVP{{seed}}{sfx}.yaml texts are byte-identical "
              f"once {{seed}}, {{k_best}} ({n_k} lines), {{k_end}} and {{clumping_distance}} ({n_cd} lines) are filled "
              "from the replicates table, and parse to each file's data", True, "600/600")
    rost = collections.Counter()
    for s in SEEDS[land]:
        c = strip_comments(TXT[(s, "hd_sweep")])
        check_copy(TXT[(s, "hd_sweep")], c, f"config_MVP{s}hd_sweep.yaml")
        write(f"{D}/offset/config_MVP{s}hd_sweep.yaml", c)
        rost[tuple(DATA[(s, "hd_sweep")]["Maladaptation"]["snp_sets"])] += 1
    ROSTER[land] = rost
    n_off = len(glob.glob(f"{OUT}/{D}/offset/config_MVP*hd_sweep.yaml"))
    have = sorted(re.search(r"MVP(\d+)hd_sweep", p).group(1) for p in glob.glob(f"{OUT}/{D}/offset/*.yaml"))
    check(name, "exactly one offset configuration per seed of work/mvp_land/seeds_*.txt", n_off == 600 and have == sorted(SEEDS[land]),
          "; ".join(f"roster of {n}: {len(r)} sets" for r, n in rost.most_common()))
    ft = {json.dumps(DATA[(s, "hd_sweep")]["Climate"]["custom"]["future_tables"]) for s in SEEDS[land]}
    check(name, "one garden list in every offset configuration", len(ft) == 1, f"{len(json.loads(ft.pop()))} garden scenarios")
    dk = collections.Counter(x for s in SEEDS[land] for x in diff_keys(DATA[(s, "hd")], DATA[(s, "hdrc")]))
    check(name, "rda_corrected differs from sets only at project_name, the GEA configs and the declared sets",
          all(re.match(r"^(project_name|GEA\.configs\[\d\]\..*|GEA\.snp_sets(\[\d+\]\..*)?)$", x) for x in dk), "; ".join(sorted(dk)))
    dk = collections.Counter(re.sub(r"^(Maladaptation)\..*", r"\1", x) for s in SEEDS[land]
                             for x in diff_keys(DATA[(s, "hd")], DATA[(s, "hd_sweep")]))
    check(name, "offset differs from sets only at cpu, Climate.custom.future_tables and Maladaptation",
          all(re.match(r"^(cpu|Climate\.custom\.future_tables|Maladaptation)$", x) for x in dk), "; ".join(sorted(dk)))

# ---- 4. SS-Clines: the published files, unchanged? ----------------------------------------------------------------------
with open(os.path.join(os.path.dirname(PUB), "data_s1_sha256.txt")) as f:
    SHA = {l.split("  ", 1)[1].strip().replace("data_s1/", "", 1): l.split("  ", 1)[0] for l in f if "/simulation/" in l}
bad = [p for p, h in SHA.items() if sha(f"{PUB}/{p}") != h]
check("SS-Clines", "published simulation/ files equal their sha256 list (phase8 data_s1_sha256.txt)", not bad and len(SHA) == 606,
      f"{len(SHA)} files")
SC = sorted(r["seed"] for r in MAN.values() if r["arm"] == "primary_ssclines")
KB = {s: int(MAN[s]["k_best"]) for s in SC}
def rd(p):
    return open(p).read()
T_C1 = rd(f"{PUB}/simulation/gea_scan_template.yaml")
T_HD = rd(f"{PUB}/simulation/sets_template.yaml")
T_UN = rd(f"{PUB}/simulation/rda_uncorrected_template.yaml")
r1 = lambda t, s: t.replace("{seed}", s).replace("{k_best}", str(KB[s]))
n_ok = 0
for s in SC:
    c1 = strip_comments(rd(f"{ARCH}/config_MVP{s}_c1.yaml"))
    hd = strip_comments(rd(f"{ARCH}/config_MVP{s}hd.yaml"))
    assert r1(T_C1, s) == c1, f"SS-Clines {s}: gea_scan template no longer renders today's _c1 file"
    assert r1(T_HD, s) == hd, f"SS-Clines {s}: sets template no longer renders today's hd file"
    assert strip_comments(rd(f"{ARCH}/config_MVP{s}hd_sweep.yaml")) == rd(f"{PUB}/simulation/offset/config_MVP{s}hd_sweep.yaml"), s
    a, b = flat(yaml.safe_load(r1(T_C1, s))), flat(yaml.safe_load(r1(T_UN, s)))
    assert [x for x in set(a) | set(b) if a.get(x) != b.get(x)] == [x for x in a if x.endswith("params.condition_pcs")
                                                                     and a[x] == KB[s] and x.startswith("GEA.configs")], s
    n_ok += 1
check("SS-Clines", "today's configurations (_archive/mvp_ssclines/configs/) render to the published templates and offset "
      "files byte for byte, so the SS-Clines folder is unchanged and not rewritten", n_ok == 600, f"{n_ok}/600 replicates")

# ---- 5. manifest ------------------------------------------------------------------------------------------------------
def nk(land, key, val):
    return sum(r[key] == val for r in REPS[land].values())

MANI = []
for land, name in LANDS.items():
    D = f"simulation_{land}"
    rost = ROSTER[land].most_common()
    extra = f" Upstream modes of the {nk(land, 'declared_619', 'yes')} replicates declared at 619 bp ran before the " \
            "declaration, under the same file with snp_clumping_distance and the sets' clumping_distance at " \
            "auto_genome_wide; those modes do not read them." if land == "estclines" else ""
    MANI += [
        (f"{D}/sets_template.yaml", f"simulation, {name}: GEA scans, LD decay and candidate sets",
         f"modes processing, prestructure, structure and gea of each of the 600 replicates (project MVP{{seed}}hd), "
         "pipeline commit c80710e: LFMM and EMMAX with K latent factors / PCs, RDA without structure correction, the "
         "framework's LD decay (all plants as one group, genome-wide) and the declared candidate sets, built by the "
         "framework: top 0.25% of the tested SNPs per method and predictor, RDA at 0.5%; each method alone, their "
         "union, and agreement of at least two and of all three methods within the replicate's LD-decay distance "
         "(critical r² rule). Replace {seed}, {k_best}, {k_end} and {clumping_distance} with a row of "
         "sets_replicates.tsv before parsing." + extra),
        (f"{D}/sets_replicates.tsv", f"simulation, {name}: GEA scans, LD decay and candidate sets",
         "one row per replicate: seed; k_best (the depositors' K, raised to 2 where it is 1: "
         f"{sum(r['k_authors'] < 2 for r in REPS[land].values())} replicates); k_end (11 where K is 11: "
         f"{nk(land, 'k_end', 11)} replicates); clumping_distance as written in the file; ld_distance_bp, the "
         "distance the candidate sets used; distance_source; declared_619 (yes for a distance declared at 619 bp "
         f"because the replicate's LD curve was unreliable: {nk(land, 'declared_619', 'yes')} replicates)"),
        (f"{D}/rda_corrected_template.yaml", f"simulation, {name}: RDA with structure correction",
         "mode gea of each replicate (project MVP{seed}hdrc), pipeline commit c80710e: RDA conditioned on K PCs, on "
         "the processing, prestructure and structure outputs of the sets run (linked by "
         "benchmarks/mvp_clone_project.sh, not recomputed). Replace {seed}, {k_best}, {k_end} and {clumping_distance} "
         "with a row of rda_corrected_replicates.tsv before parsing."),
        (f"{D}/rda_corrected_replicates.tsv", f"simulation, {name}: RDA with structure correction",
         "one row per replicate, equal to sets_replicates.tsv"),
        (f"{D}/offset/config_MVP{{seed}}hd_sweep.yaml (600 files)", f"simulation, {name}: genetic offset",
         "mode maladaptation of each replicate over 112 garden scenarios, pipeline commit c80710e: Gradient Forest, "
         "geometric offset, RDA-based offset (rda_offset, without structure correction) and rda_offset_corrected (run, "
         "not reported); snp_sets lists the replicate's declared sets of at least 3 SNPs (built in the sets run) and "
         "its causal loci (truth)" + ("" if len(rost) == 1 else
                                      f"; {sum(n for r, n in rost[1:])} replicates list fewer sets")),
    ]
with open(f"{OUT}/manifest.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["file", "analysis", "what it configures"])
    w.writerows(MANI)
written.append("manifest.tsv")
txt = open(f"{OUT}/manifest.tsv").read()
check("manifest.tsv", "no banned word (risk, adaptive, null, Climate-only, >=, salinity, elevation)",
      not re.search(r"risk|adaptive|null|Climate-only|>=|salinity|elevation", txt, re.I))
for row in MANI:
    p = row[0].split(" (")[0]
    assert glob.glob(f"{OUT}/" + p.replace("{seed}", "*")), p
on_disk = sorted(os.path.relpath(os.path.join(d, f), OUT) for land in LANDS
                 for d, _, fs in os.walk(f"{OUT}/simulation_{land}") for f in fs) + ["manifest.tsv"]
check("Job A", "the folders hold exactly what was written", sorted(on_disk) == sorted(written), f"{len(written)} files")

with open(f"{OUT}/checks/data_s1_check.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["item", "check", "result", "detail"])
    w.writerows(CHECKS)
print("done")
