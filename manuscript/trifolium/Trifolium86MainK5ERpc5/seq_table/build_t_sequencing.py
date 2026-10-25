#!/usr/bin/env python3
"""Sequencing of the Trifolium israeliticum plants (supplement slug t_sequencing, provisional Table S12):
read pairs, mapping, duplicates, mean depth, missing calls and filter outcome of each called plant.

Inputs (all copied read-only from the lab server into <in_dir>):
  per_plant_seq_stats.tsv    per_plant_seq_stats.sh output (read pairs, mapped %, duplicates %,
                             mean depth over Chr_1-Chr_6, missing calls %), one row per called plant;
                             its site column is the plant-name prefix and is not used here
  Trifolium_location_86.tsv  the framework's input metadata, the site of every called plant
                             ([AG]/data/Trifolium_location_86.tsv)
  site_names.tsv             site_id -> short name 1-11, south to north
                             (ANALYSIS/Trifolium86MainK5ERpc5/fig3/site_names.tsv)
  samples_missing_stats.tsv  plink --missing over the framework's input VCF
                             (<run>_results/_intermediate/samples/samples_missing_stats.tsv)
  relatedness_removed.tsv    plants removed as near-identical (<run>_results/Processing/tables/)
  metadata.tsv               plants retained after all filters (<run>_results/Processing/tables/)
The framework's order is sample missingness (F_MISS > 0.8 removed), then relatedness (IBS > 0.999,
the member with more missing calls removed), so a plant has exactly one outcome.
Each plant's site comes from the metadata, never from its name: plant 29107-6 belongs to site 29017.
The script prints every plant whose name prefix differs from its site; the table has no note column (user,
2026-10-25: not needed).

Usage: python3 -I build_t_sequencing.py <in_dir> <out_prefix>
Writes <out_prefix>.tsv and <out_prefix>.md, and prints the median and range of mean depth over all
called plants and over the retained ones.
"""
import csv
import statistics
import sys


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def summary(rows):
    d = [float(o["mean_depth"]) for o in rows]
    p = [int(o["read_pairs"]) for o in rows]
    return d, p


def main(in_dir, out_prefix):
    stats = read_tsv(f"{in_dir}/per_plant_seq_stats.tsv")
    site_of = {r["sample"]: r["site"] for r in read_tsv(f"{in_dir}/Trifolium_location_86.tsv")}
    names = {r["site_id"]: r["short"] for r in read_tsv(f"{in_dir}/site_names.tsv")}
    removed_rel = {r["IID"] for r in read_tsv(f"{in_dir}/relatedness_removed.tsv")}
    retained = {r["sample"]: r["site"] for r in read_tsv(f"{in_dir}/metadata.tsv")}
    with open(f"{in_dir}/samples_missing_stats.tsv", encoding="utf-8") as fh:
        rows = [line.split() for line in fh]
    fmiss = {r[1]: float(r[5]) for r in rows[1:]}

    if len(stats) != len(site_of):
        sys.exit(f"{len(stats)} plants with statistics but {len(site_of)} in the metadata")
    out = []
    for r in stats:
        plant = r["plant"]
        if plant not in site_of:
            sys.exit(f"{plant}: not in the framework's metadata")
        site = site_of[plant]
        if site not in names:
            sys.exit(f"{plant}: site {site} has no short name")
        if plant in retained and retained[plant] != site:
            sys.exit(f"{plant}: site {site} in the input metadata, {retained[plant]} after filtering")
        if plant not in fmiss:
            sys.exit(f"no missingness row for {plant}")
        if plant in retained:
            outcome = "retained"
        elif fmiss[plant] > 0.8:
            outcome = "removed: missing calls > 80%"
        elif plant in removed_rel:
            outcome = "removed: identity by state > 0.999"
        else:
            sys.exit(f"{plant}: neither retained nor removed by a known filter")
        out.append({
            "site": names[site],
            "site_code": site,
            "plant": plant,
            "read_pairs": r["read_pairs"],
            "mapped_pct": r["mapped_pct"],
            "duplicate_pct": r["duplicate_pct"],
            "mean_depth": r["mean_depth_chr1_6"],
            "missing_calls_pct": r["missing_calls_pct"],
            "outcome": outcome,
        })
    missing = set(retained) - {o["plant"] for o in out}
    if missing:
        sys.exit("retained plants without sequencing statistics: " + ", ".join(sorted(missing)))

    out.sort(key=lambda o: (int(o["site"]), o["plant"]))
    with open(f"{out_prefix}.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    kept = [o for o in out if o["outcome"] == "retained"]
    low = [o for o in out if o["outcome"] == "removed: missing calls > 80%"]
    d_all, _ = summary(out)
    d_kept, _ = summary(kept)
    d_low, _ = summary(low)
    lines = [
        "# t_sequencing (provisional Table S12): sequencing of the 86 called plants of *Trifolium israeliticum*",
        "",
        "Draft of 2026-10-25 (user, 2026-10-25: a separate supplementary table, read pairs, depth and missing calls per "
        "plant). Values: `t_sequencing.tsv`. Read pairs, mapped and duplicate percentages from `samtools flagstat` (1.17) "
        "on each plant's BAM; mean depth from `samtools coverage` over the six pseudochromosomes Chr_1–Chr_6, "
        "length-weighted, duplicates, secondary, QC-failed and unmapped reads excluded, no quality threshold; missing calls "
        "from `plink --missing` over the 7,163,793 SNPs the framework read. Built by "
        "`../../scripts/trifolium/per_plant_seq_stats.sh` (on the lab server) and `build_t_sequencing.py`. Sites 1–11 run "
        "from south to north, each plant's site taken from the framework's metadata; outcome = the framework's filters in "
        "their order (missing calls > 80%, then identity by state > 0.999, the member with more missing calls removed).",
        "",
        f"Mean depth: median {statistics.median(d_all):.1f}× (range {min(d_all):.1f}–{max(d_all):.1f}×) over the "
        f"{len(out)} called plants; {statistics.median(d_kept):.1f}× ({min(d_kept):.1f}–{max(d_kept):.1f}×) over the "
        f"{len(kept)} retained; {min(d_low):.1f}–{max(d_low):.1f}× over the {len(low)} removed for missing calls. The 26 "
        "sequenced plants that were not called (all 10 of site 29014 among them) have no BAM and are not listed.",
        "",
        "| Site | Site code | Plant | Read pairs | Mapped, % | Duplicates, % | Mean depth (×) | Missing calls, % "
        "| Outcome in the framework |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for o in out:
        lines.append(f"| {o['site']} | {o['site_code']} | {o['plant']} | {int(o['read_pairs']):,} | {o['mapped_pct']} "
                     f"| {o['duplicate_pct']} | {o['mean_depth']} | {o['missing_calls_pct']} | {o['outcome']} |")
    with open(f"{out_prefix}.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    for label, sel in (("all called", out), ("retained", kept)):
        d, p = summary(sel)
        print(f"{label}: n = {len(sel)}; mean depth median {statistics.median(d):.1f} "
              f"(range {min(d):.1f}-{max(d):.1f}); read pairs median {statistics.median(p):,} "
              f"(range {min(p):,}-{max(p):,})")
    counts = {}
    for o in out:
        counts[o["outcome"]] = counts.get(o["outcome"], 0) + 1
    print("outcomes:", counts)
    per_site = {}
    for o in kept:
        per_site[int(o["site"])] = per_site.get(int(o["site"]), 0) + 1
    print("retained plants per site:", dict(sorted(per_site.items())))
    print("name prefix differs from the site:",
          [(o["plant"], o["site_code"]) for o in out if o["plant"].split("-")[0] != o["site_code"]])


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
