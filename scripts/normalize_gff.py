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


def vcf_contigs(path: str) -> set:
    seen = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            seen.add(line.split("\t", 1)[0])
    return seen


def normalize(gff_in: str, gff_out: str) -> set:
    seqids = set()
    with open(gff_in) as fin, open(gff_out, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            seqid, sep, rest = line.partition("\t")
            if not sep:  # not a 9-column record; keep verbatim
                fout.write(line)
                continue
            seqid = canonical(seqid)
            seqids.add(seqid)
            fout.write(seqid + sep + rest)
    return seqids


def main(argv):
    if len(argv) != 4:
        sys.exit(f"usage: {argv[0]} <raw.gff3|NULL> <filtered.vcf> <normalized.gff3>")
    gff_in, vcf, gff_out = argv[1:4]

    if gff_in == "NULL":
        open(gff_out, "w").close()
        print("INFO: No GFF provided, created empty normalized GFF")
        return 0

    seqids = normalize(gff_in, gff_out)
    print(f"INFO: Normalized GFF chromosome names (case-insensitive 'chr' strip): "
          f"{len(seqids)} seqids")

    contigs = vcf_contigs(vcf)
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
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
