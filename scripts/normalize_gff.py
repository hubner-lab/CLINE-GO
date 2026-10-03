#!/usr/bin/env python3
"""Normalize GFF3 seqids to the pipeline's chromosome-name contract and verify them
against the filtered VCF.

Usage: normalize_gff.py <raw.gff3 | NULL> <filtered.vcf> <normalized.gff3>

The contract (shared with rule filter_vcf in workflow/rules/processing.smk):
  * a leading ``chr`` prefix is stripped in ANY case (Chr1 -> 1, CHR2 -> 2, chr2H -> 2H)
  * sex/organelle contigs stay letters (X, Y, XY, MT); ``M`` is written as ``MT`` because
    plink (which re-emits the VCF with --output-chr MT) spells the mitochondrion that way.
Header/comment lines (``#...``) are copied verbatim.

After writing, the seqid set of the GFF body is compared with the contig set of the
filtered VCF body. No shared contig is an error (exit 1) — every downstream gene/region
join keys on ``chr``, so the annotation would be silently empty otherwise (audit
2026-09-13 B1). A partial overlap is logged as a WARNING with both one-sided lists.

Matching NAMES do not imply matching COORDINATES: a same-names annotation from another
assembly (a Morex v2 GFF against a v3-called VCF, both ``chr1H``..``chr7H``) used to pass
clean and every gene, promoter window, exon count and GO term was then computed at the
wrong positions (finding fd0070). The coordinate extents are therefore compared too, but
only where a disagreement is PHYSICALLY IMPOSSIBLE rather than merely unequal:

  * a GFF feature ending past the VCF's declared ``##contig=<...,length=>``
  * a VCF variant positioned past the GFF's declared ``##sequence-region`` end

Both are reported as WARNINGs. A GFF whose features stop well short of the chromosome
length is NOT flagged — that is an ordinary sparse annotation (SIMDATA's own synthetic
GFF reaches 25 Mb of a 41.7 Mb contig), so a tolerance on the ratio would be noise. The
per-contig extents are printed either way, so the number is there to check by hand.

``NULL`` as the GFF path (no annotation configured) writes an empty output and exits 0.
"""
import re
import sys

_CHR_PREFIX = re.compile(r"^chr", re.IGNORECASE)
# The five codes plink recognises, keyed by UPPER case: plink parses them
# case-insensitively (Mt, mt, x all count) and --output-chr MT writes them back
# as X / Y / XY / MT. A plant annotation spelling the mitochondrion `Mt` (TAIR
# style) would otherwise sit on `Mt` while its own VCF came out as `MT`, and the
# join would silently lose that one contig. Everything else — Pt, 2H, Un,
# scaffold_1 — plink keeps verbatim, and so does this.
_PLINK_CODES = {"X": "X", "Y": "Y", "XY": "XY", "M": "MT", "MT": "MT"}


def canonical(name: str) -> str:
    """The pipeline's canonical chromosome name for a raw contig name."""
    name = _CHR_PREFIX.sub("", name)
    return _PLINK_CODES.get(name.upper(), name)


_CONTIG_HDR = re.compile(r"^##contig=<(?=.*\bID=(?P<id>[^,>]+))(?=.*\blength=(?P<len>[0-9]+))")


def vcf_contigs(path: str):
    """Contig names in the VCF body, their declared header lengths and their max POS.

    Returns (names, declared, max_pos). `declared` holds only the contigs whose
    ``##contig`` header carries a length; `max_pos` is collected from the body, so a
    contig declared in the header but carrying no variant is absent from it.
    """
    seen = set()
    declared: dict = {}
    max_pos: dict = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##contig="):
                    m = _CONTIG_HDR.match(line)
                    if m:
                        declared[m.group("id")] = int(m.group("len"))
                continue
            fields = line.split("\t", 2)
            name = fields[0]
            seen.add(name)
            if len(fields) > 1:
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                if pos > max_pos.get(name, 0):
                    max_pos[name] = pos
    return seen, declared, max_pos


def normalize(gff_in: str, gff_out: str):
    """Write the normalized GFF and return (seqids, max_end, region_end).

    `max_end` is the largest field-5 coordinate seen per canonical seqid; `region_end`
    the end declared by a ``##sequence-region`` directive, when the GFF carries one.
    Both keys are canonicalised, so they compare directly against the VCF's contigs.
    """
    seqids = set()
    max_end: dict = {}
    region_end: dict = {}
    with open(gff_in) as fin, open(gff_out, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            region_end[canonical(parts[1])] = int(parts[3])
                        except ValueError:
                            pass
                fout.write(line)
                continue
            seqid, sep, rest = line.partition("\t")
            if not sep:  # not a 9-column record; keep verbatim
                fout.write(line)
                continue
            seqid = canonical(seqid)
            seqids.add(seqid)
            cols = rest.split("\t")
            if len(cols) >= 4:
                try:
                    end = int(cols[3])
                except ValueError:
                    end = None
                if end is not None and end > max_end.get(seqid, 0):
                    max_end[seqid] = end
            fout.write(seqid + sep + rest)
    return seqids, max_end, region_end


def check_extents(shared, gff_max_end, gff_region_end, vcf_declared, vcf_max_pos):
    """Report coordinate disagreements that no matching assembly can produce.

    Only two comparisons are made, both one-sided and both impossible on a matching
    build: a GFF feature past the VCF's declared contig length, and a VCF variant past
    the GFF's declared sequence-region end. Returns the number of contigs flagged.
    """
    print("INFO: contig extents (GFF max feature end / VCF declared length / VCF max POS):")
    flagged = []
    for name in shared:
        g_end = gff_max_end.get(name)
        v_len = vcf_declared.get(name)
        v_pos = vcf_max_pos.get(name)

        def fmt(v):
            return "-" if v is None else f"{v:,}"

        print(f"INFO:   {name}: {fmt(g_end)} / {fmt(v_len)} / {fmt(v_pos)}")

        if g_end is not None and v_len is not None and g_end > v_len:
            flagged.append(name)
            print(f"WARNING: {name}: a GFF feature ends at {g_end:,} but the VCF declares "
                  f"the contig as {v_len:,} bp long — the annotation cannot belong to this "
                  f"assembly. Every gene, promoter window, exon count and GO term on "
                  f"{name} is computed at the wrong coordinates.")
        r_end = gff_region_end.get(name)
        if r_end is not None and v_pos is not None and v_pos > r_end:
            if name not in flagged:
                flagged.append(name)
            print(f"WARNING: {name}: a VCF variant sits at {v_pos:,} but the GFF declares "
                  f"##sequence-region ending at {r_end:,} — the two inputs are from "
                  f"different assemblies.")

    if flagged:
        print(f"WARNING: {len(flagged)} contig(s) have incompatible GFF/VCF coordinates: "
              f"{', '.join(flagged)}. Names match, coordinates do not — check the "
              f"assembly version of both inputs.")
    return len(flagged)


def main(argv):
    if len(argv) != 4:
        sys.exit(f"usage: {argv[0]} <raw.gff3|NULL> <filtered.vcf> <normalized.gff3>")
    gff_in, vcf, gff_out = argv[1:4]

    if gff_in == "NULL":
        open(gff_out, "w").close()
        print("INFO: No GFF provided, created empty normalized GFF")
        return 0

    seqids, gff_max_end, gff_region_end = normalize(gff_in, gff_out)
    print(f"INFO: Normalized GFF chromosome names (case-insensitive 'chr' strip): "
          f"{len(seqids)} seqids")

    contigs, vcf_declared, vcf_max_pos = vcf_contigs(vcf)
    shared = seqids & contigs
    gff_only = sorted(seqids - contigs)
    vcf_only = sorted(contigs - seqids)

    def show(names, n=10):
        return ", ".join(names[:n]) + (f", ... (+{len(names) - n})" if len(names) > n else "")

    if not shared:
        print("ERROR: the GFF and the filtered VCF share NO chromosome name after "
              "normalisation — every gene/region join would be empty.", file=sys.stderr)
        print(f"ERROR:   GFF seqids ({len(seqids)}): {show(gff_only)}", file=sys.stderr)
        print(f"ERROR:   VCF contigs ({len(contigs)}): {show(vcf_only)}", file=sys.stderr)
        print("ERROR: rename the contigs in one of the inputs so they match "
              "(chr prefix, X/Y/MT spelling).", file=sys.stderr)
        return 1

    print(f"INFO: {len(shared)} chromosome name(s) shared by GFF and VCF")
    if vcf_only:
        print(f"WARNING: {len(vcf_only)} VCF contig(s) have no GFF annotation: {show(vcf_only)}")
    if gff_only:
        print(f"WARNING: {len(gff_only)} GFF seqid(s) absent from the VCF (no SNPs there): "
              f"{show(gff_only)}")

    check_extents(sorted(shared), gff_max_end, gff_region_end, vcf_declared, vcf_max_pos)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
