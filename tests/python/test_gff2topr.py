"""scripts/gff2topr.py — subprocess only, because it cannot be imported.

The file has no `if __name__` guard and defines zero functions: its whole body runs
at import and calls sys.exit(1) on a bad argument count. So the only testable
surface is the process boundary, which is exactly the Tier 6 smoke-test shape.

Its sole caller is the Shiny app (scripts/clinego.app/R/fct_regions.R:588-598, via
processx::run with a 30 s timeout). No Snakemake rule invokes it.
"""
import os
import subprocess
import sys
import tempfile
import unittest

import _support

SCRIPT = _support.script_path("gff2topr")

HEADER = ["chrom", "gene_start", "gene_end", "gene_symbol", "biotype",
          "exon_chromstart", "exon_chromend"]


def gff_line(chrom, feature, start, end, attrs):
    return "\t".join([chrom, "src", feature, str(start), str(end), ".", "+", ".", attrs])


class Gff2ToprCase(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.d = self.dir.name
        self.out = os.path.join(self.d, "out.tsv")

    def write_gff(self, lines, name="in.gff3"):
        p = os.path.join(self.d, name)
        with open(p, "w") as fh:
            fh.write("##gff-version 3\n")
            fh.write("\n".join(lines) + "\n")
        return p

    def run_script(self, *args):
        return subprocess.run([sys.executable, SCRIPT, *args],
                              capture_output=True, text=True)

    def read_out(self):
        with open(self.out) as fh:
            return [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]


class TestArgumentHandling(Gff2ToprCase):

    def test_wrong_argument_count_exits_1_with_usage(self):
        r = self.run_script("only", "three", "args")
        self.assertEqual(r.returncode, 1)
        self.assertIn("Usage:", r.stdout)

    def test_no_arguments_exits_1(self):
        self.assertEqual(self.run_script().returncode, 1)

    def test_missing_input_file_exits_1(self):
        r = self.run_script(os.path.join(self.d, "nope.gff3"),
                            "gene", "Name", "biotype", self.out)
        self.assertEqual(r.returncode, 1)
        self.assertIn("File not found", r.stdout)


class TestConversion(Gff2ToprCase):

    def test_three_genes_with_exons(self):
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
            gff_line("1", "exon", 120, 300, "ID=e1;Parent=g1"),
            gff_line("1", "exon", 500, 700, "ID=e2;Parent=g1"),
            gff_line("2", "gene", 50, 400, "ID=g2;Name=BETA"),
            gff_line("2", "exon", 60, 120, "ID=e3;Parent=g2"),
            gff_line("3H", "gene", 10, 20, "ID=g3;Name=GAMMA"),
        ])
        r = self.run_script(gff, "gene", "Name", "biotype", self.out)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

        rows = self.read_out()
        self.assertEqual(rows[0], HEADER)
        self.assertEqual(len(rows), 4)                      # header + 3 genes

        by_symbol = {r[3]: r for r in rows[1:]}
        self.assertEqual(by_symbol["ALPHA"][:3], ["1", "100", "900"])
        self.assertEqual(by_symbol["ALPHA"][5], "120,500")  # exon starts, in file order
        self.assertEqual(by_symbol["ALPHA"][6], "300,700")
        self.assertEqual(by_symbol["BETA"][5], "60")
        self.assertEqual(by_symbol["GAMMA"][5], "")         # no exons -> empty field
        # Chromosome names pass through verbatim; normalisation is upstream.
        self.assertEqual(by_symbol["GAMMA"][0], "3H")

    def test_malformed_and_comment_lines_are_skipped(self):
        gff = self.write_gff([
            "# a comment in the middle",
            "1\tsrc\tgene\t100",                             # < 9 columns
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
        ])
        r = self.run_script(gff, "gene", "Name", "biotype", self.out)
        self.assertEqual(r.returncode, 0)
        self.assertEqual(len(self.read_out()), 2)

    def test_feature_without_the_name_field_is_dropped(self):
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1"),        # no Name=
            gff_line("1", "gene", 200, 800, "ID=g2;Name=BETA"),
        ])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        rows = self.read_out()
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1][3], "BETA")

    def test_a_non_gene_feature_type_is_selectable(self):
        gff = self.write_gff([
            gff_line("1", "mRNA", 100, 900, "ID=m1;Name=TRANSCRIPT"),
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
        ])
        self.assertEqual(self.run_script(gff, "mRNA", "Name", "biotype", self.out).returncode, 0)
        rows = self.read_out()
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1][3], "TRANSCRIPT")


class TestPinnedDefects(Gff2ToprCase):
    """Current behaviour, filed in docs/pipeline_improvement_requests.md.

    Each of these asserts what the script DOES, not what it should do. Fixing one
    means changing the assertion in the same commit as the fix.
    """

    def test_exon_before_its_parent_gene_is_silently_dropped(self):
        # gff2topr.py:34 `if parent_id in genes` — an exon is only collected when
        # its gene line was seen FIRST. GFF3 mandates no such ordering.
        gff = self.write_gff([
            gff_line("1", "exon", 120, 300, "ID=e1;Parent=g1"),
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
        ])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        rows = self.read_out()
        self.assertEqual(rows[1][3], "ALPHA")
        self.assertEqual(rows[1][5], "")      # the exon is gone

    def test_biotype_argument_is_accepted_and_ignored(self):
        # gff2topr.py:47,52 — biotype extraction is commented out ("REMOVE for
        # spontaneum") and gene_type is hardcoded. argv[4] has no effect at all.
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA;biotype=lncRNA"),
        ])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        self.assertEqual(self.read_out()[1][4], "protein_coding")

        self.assertEqual(self.run_script(gff, "gene", "Name", "TOTAL_NONSENSE",
                                         self.out).returncode, 0)
        self.assertEqual(self.read_out()[1][4], "protein_coding")

    def test_arguments_are_echoed_to_stdout(self):
        # gff2topr.py:16 prints all five arguments unconditionally, into the log
        # of whatever called it.
        gff = self.write_gff([gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA")])
        r = self.run_script(gff, "gene", "Name", "biotype", self.out)
        self.assertIn("\t".join([gff, "gene", "Name", "biotype", self.out]), r.stdout)


if __name__ == "__main__":
    unittest.main()
