"""scripts/normalize_gff.py — the GFF half of the chromosome-name contract.

Rule filter_vcf strips `chr` case-insensitively from the VCF body and plink writes
X/Y/MT as letters (--output-chr MT); this script applies the same rule to the GFF
and then refuses a GFF whose seqids share nothing with the filtered VCF's contigs
(audit 2026-09-13 B1: a TAIR10-style `Chr1` GFF used to pass through untouched
while the VCF came out as `1`, so every annotation table was empty with exit 0).
"""
import os
import subprocess
import sys
import tempfile
import unittest

import _support

SCRIPT = _support.script_path("normalize_gff")
mod = _support.load_script("normalize_gff")


def gff_line(chrom, start=1, end=10, attrs="ID=g1"):
    return "\t".join([chrom, "src", "gene", str(start), str(end), ".", "+", ".", attrs])


class CanonicalCase(unittest.TestCase):

    def test_chr_prefix_stripped_in_any_case(self):
        for raw, want in [("chr1", "1"), ("Chr1", "1"), ("CHR2", "2"), ("chr2H", "2H"),
                          ("Chr5H", "5H"), ("chrX", "X"), ("Un", "Un"), ("1", "1")]:
            self.assertEqual(mod.canonical(raw), want, raw)

    def test_mitochondrion_spelled_like_plink(self):
        self.assertEqual(mod.canonical("M"), "MT")
        self.assertEqual(mod.canonical("chrM"), "MT")
        self.assertEqual(mod.canonical("MT"), "MT")

    def test_x_y_stay_letters(self):
        self.assertEqual(mod.canonical("X"), "X")
        self.assertEqual(mod.canonical("Y"), "Y")


class NormalizeGffCase(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.d = self.dir.name
        self.out = os.path.join(self.d, "normalized.gff3")

    def write(self, name, lines):
        p = os.path.join(self.d, name)
        with open(p, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return p

    def vcf(self, contigs):
        header = ["##fileformat=VCFv4.2",
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
        body = [f"{c}\t{100 + i}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1" for i, c in enumerate(contigs)]
        return self.write("filtered.vcf", header + body)

    def invoke(self, gff, vcf):
        return subprocess.run([sys.executable, SCRIPT, gff, vcf, self.out],
                              capture_output=True, text=True)

    def body_seqids(self):
        with open(self.out) as fh:
            return [ln.split("\t")[0] for ln in fh if ln.strip() and not ln.startswith("#")]

    def test_mixed_case_gff_matches_plink_side_names(self):
        gff = self.write("in.gff3", ["##gff-version 3",
                                     gff_line("Chr1"), gff_line("CHR2"), gff_line("chr3"),
                                     gff_line("Chr5H"), gff_line("M"), gff_line("Un")])
        # what filter_vcf + plink --output-chr MT actually emit for the same raw names
        res = self.invoke(gff, self.vcf(["1", "2", "3", "5H", "MT", "Un"]))
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(self.body_seqids(), ["1", "2", "3", "5H", "MT", "Un"])
        self.assertIn("6 chromosome name(s) shared", res.stdout)
        self.assertNotIn("WARNING", res.stdout)

    def test_header_and_comment_lines_are_kept_verbatim(self):
        gff = self.write("in.gff3", ["##gff-version 3", "# a comment", gff_line("chr1")])
        self.invoke(gff, self.vcf(["1"]))
        with open(self.out) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[:2], ["##gff-version 3", "# a comment"])

    def test_disjoint_contig_sets_exit_1(self):
        gff = self.write("in.gff3", [gff_line("Chr1"), gff_line("Chr2")])
        res = self.invoke(gff, self.vcf(["scaffold_1", "scaffold_2"]))
        self.assertEqual(res.returncode, 1)
        self.assertIn("share NO chromosome name", res.stderr)
        self.assertIn("scaffold_1", res.stderr)
        self.assertIn("1, 2", res.stderr)       # the GFF side, post-normalisation

    def test_partial_overlap_warns_with_both_one_sided_lists(self):
        gff = self.write("in.gff3", [gff_line("chr1"), gff_line("chr2"), gff_line("chr9")])
        res = self.invoke(gff, self.vcf(["1", "2", "X"]))
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertIn("2 chromosome name(s) shared", res.stdout)
        self.assertIn("VCF contig(s) have no GFF annotation: X", res.stdout)
        self.assertIn("GFF seqid(s) absent from the VCF (no SNPs there): 9", res.stdout)

    def test_no_gff_writes_an_empty_file_and_exits_0(self):
        res = self.invoke("NULL", self.vcf(["1"]))
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(os.path.getsize(self.out), 0)


if __name__ == "__main__":
    unittest.main()
