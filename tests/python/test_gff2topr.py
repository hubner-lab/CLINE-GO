"""scripts/gff2topr.py — subprocess only, because it cannot be imported.

The file has no `if __name__` guard: its whole body runs at import and calls
sys.exit(1) on a bad argument count, so the only testable surface is the process
boundary, which is exactly the Tier 6 smoke-test shape.

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

    def test_feature_without_the_name_field_keeps_its_ID_as_the_symbol(self):
        # Was: the gene was dropped from the annotation entirely. A missing
        # display name is not a reason to omit a gene from the regional plot.
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1"),        # no Name=
            gff_line("1", "gene", 200, 800, "ID=g2;Name=BETA"),
        ])
        r = self.run_script(gff, "gene", "Name", "biotype", self.out)
        self.assertEqual(r.returncode, 0)
        rows = self.read_out()
        self.assertEqual(len(rows), 3)
        self.assertEqual([row[3] for row in rows[1:]], ["g1", "BETA"])
        self.assertIn("carry no 'Name' attribute", r.stderr)

    def test_a_non_gene_feature_type_is_selectable(self):
        gff = self.write_gff([
            gff_line("1", "mRNA", 100, 900, "ID=m1;Name=TRANSCRIPT"),
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
        ])
        self.assertEqual(self.run_script(gff, "mRNA", "Name", "biotype", self.out).returncode, 0)
        rows = self.read_out()
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1][3], "TRANSCRIPT")


class TestChildResolution(Gff2ToprCase):
    """Exon/CDS attachment must not depend on line order or on feature depth."""

    def test_exon_before_its_parent_gene_is_still_attached(self):
        # Was pinned the other way round: `if parent_id in genes` collected an
        # exon only when its gene line had been read FIRST, and GFF3 mandates no
        # such ordering. The conversion is two-pass now.
        gff = self.write_gff([
            gff_line("1", "exon", 120, 300, "ID=e1;Parent=g1"),
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
        ])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        rows = self.read_out()
        self.assertEqual(rows[1][3], "ALPHA")
        self.assertEqual(rows[1][5], "120")
        self.assertEqual(rows[1][6], "300")

    def test_exons_reach_the_gene_through_the_mRNA(self):
        # The canonical three-level GFF3 model: exon.Parent is the mRNA, not the
        # gene, so with GFF.feature=gene every exon used to be discarded.
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
            gff_line("1", "mRNA", 100, 900, "ID=g1.1;Parent=g1"),
            gff_line("1", "exon", 120, 300, "ID=e1;Parent=g1.1"),
            gff_line("1", "exon", 500, 700, "ID=e2;Parent=g1.1"),
        ])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        rows = self.read_out()
        self.assertEqual(rows[1][5], "120,500")
        self.assertEqual(rows[1][6], "300,700")

    def test_CDS_is_used_when_the_GFF_has_no_exon_features(self):
        # Same fallback read_gff_exons() applies (scripts/R/lib/gff_parsing.R).
        # SIMDATA's own GFF is CDS-only, so this is the test dataset's shape.
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
            gff_line("1", "mRNA", 100, 900, "ID=g1.1;Parent=g1"),
            gff_line("1", "CDS", 150, 800, "ID=c1;Parent=g1.1"),
        ])
        r = self.run_script(gff, "gene", "Name", "biotype", self.out)
        self.assertEqual(r.returncode, 0)
        self.assertEqual(self.read_out()[1][5], "150")
        self.assertIn("using CDS", r.stderr)

    def test_exon_features_win_over_CDS_when_both_are_present(self):
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
            gff_line("1", "exon", 120, 300, "ID=e1;Parent=g1"),
            gff_line("1", "CDS", 150, 280, "ID=c1;Parent=g1"),
        ])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        self.assertEqual(self.read_out()[1][5], "120")

    def test_an_exon_with_no_resolvable_parent_is_reported(self):
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA"),
            gff_line("1", "exon", 120, 300, "ID=e1;Parent=ghost"),
        ])
        r = self.run_script(gff, "gene", "Name", "biotype", self.out)
        self.assertEqual(r.returncode, 0)
        self.assertEqual(self.read_out()[1][5], "")
        self.assertIn("no 'gene' ancestor", r.stderr)


class TestBiotypeAndLogging(Gff2ToprCase):

    def test_biotype_attribute_is_honoured_when_present(self):
        # Was hardcoded to protein_coding, which made argv[4] — the app's
        # GFF.biotype config value — completely inert.
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA;biotype=lncRNA"),
        ])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        self.assertEqual(self.read_out()[1][4], "lncRNA")

    def test_biotype_falls_back_to_protein_coding_when_absent(self):
        # The H. spontaneum GFF has no biotype attribute; requiring one used to
        # drop every gene, which is why extraction was commented out.
        gff = self.write_gff([gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA")])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        self.assertEqual(self.read_out()[1][4], "protein_coding")

    def test_an_attribute_key_is_matched_whole_not_as_a_suffix(self):
        # 'Name' must not be read out of 'gene_Name='.
        gff = self.write_gff([
            gff_line("1", "gene", 100, 900, "ID=g1;gene_Name=WRONG;Name=RIGHT"),
        ])
        self.assertEqual(self.run_script(gff, "gene", "Name", "biotype", self.out).returncode, 0)
        self.assertEqual(self.read_out()[1][3], "RIGHT")

    def test_arguments_are_echoed_to_stderr_not_stdout(self):
        gff = self.write_gff([gff_line("1", "gene", 100, 900, "ID=g1;Name=ALPHA")])
        r = self.run_script(gff, "gene", "Name", "biotype", self.out)
        self.assertEqual(r.stdout, "")
        self.assertIn("\t".join([gff, "gene", "Name", "biotype", self.out]), r.stderr)


if __name__ == "__main__":
    unittest.main()
