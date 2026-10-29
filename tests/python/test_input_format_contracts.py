"""Input-format contracts every rule and script must keep (release prep, Phase 2b).

Each one is a property of the code text, one call site among many, so it is checked as
text over workflow/rules/*.smk and the scripts (comment lines skipped):

  * 12f364  every plink call reads chromosomes with --chr-set 95. Under plink's default
            human set a contig named 23/24/25/26 IS X/Y/XY/MT (measured, plink 1.90b7.2:
            a 26 even merges with a real MT), so one call without the flag re-emits a
            contig under another name than every other step.
  * 5f1e90  the VCF's sample header is found by an anchored ^#CHROM: an unanchored CHROM
            first matches any earlier ##contig=<ID=CHROMOSOME_I> line.
  * ab75b3  no header read is capped at a line count: `head -1000` never reached a
            #CHROM pushed past line 1000 by one ##contig line per scaffold.
"""
import glob
import os
import re
import unittest

import _support

ROOT = _support.REPO_ROOT


def _files(*patterns):
    found = set()
    for pattern in patterns:
        found.update(glob.glob(os.path.join(ROOT, pattern), recursive=True))
    return sorted(found)


SMK = _files("workflow/rules/*.smk")
CODE = SMK + _files("scripts/*.R", "scripts/R/**/*.R", "scripts/clinego.app/R/*.R",
                    "scripts/*.py", "scripts/*.sh")


def _read(path):
    with open(path) as handle:
        return handle.read()


def _code_lines(path):
    """(line number, line) for every line that is not a comment."""
    return [(i, line) for i, line in enumerate(_read(path).splitlines(), start=1)
            if not line.lstrip().startswith("#")]


def _plink_commands(text):
    """Each shell line of a .smk that starts `plink --`, joined over `\\` continuations
    (a flag must follow: docstrings say "plink canonicalised ...")."""
    lines = text.splitlines()
    commands = []
    for i, line in enumerate(lines):
        if re.match(r"\s*plink\s+--", line):
            command, j = [line], i
            while command[-1].rstrip().endswith("\\") and j + 1 < len(lines):
                j += 1
                command.append(lines[j])
            commands.append(" ".join(command))
    return commands


class TestPlinkChromosomeSet(unittest.TestCase):    # 12f364
    def test_every_rule_plink_call_passes_chr_set_95(self):
        n = 0
        for path in SMK:
            for command in _plink_commands(_read(path)):
                n += 1
                with self.subTest(file=os.path.basename(path), command=command[:90]):
                    self.assertIn("--chr-set 95", command)
        self.assertGreaterEqual(n, 18, "the rule parser found fewer plink calls than exist")

    def test_every_script_that_runs_plink_passes_chr_set_95(self):
        callers = [p for p in CODE if not p.endswith(".smk")
                   and re.search(r"""system2?\(\s*["']plink["']|["']plink --""", _read(p))]
        names = {os.path.basename(p) for p in callers}
        self.assertTrue({"ld_decay_pairs.R", "ld_decay_background.R",
                         "run_haplotype_scan.R"} <= names, names)
        for path in callers:
            with self.subTest(file=os.path.basename(path)):
                self.assertRegex(_read(path), r"""--chr-set["']?,?\s*["']?95""")


class TestVcfHeaderReads(unittest.TestCase):        # 5f1e90, ab75b3
    def test_a_grep_for_the_sample_header_is_anchored(self):
        hits = 0
        for path in CODE:
            for i, line in _code_lines(path):
                if re.search(r"grep\b.*CHROM", line):
                    hits += 1
                    with self.subTest(file=f"{os.path.basename(path)}:{i}"):
                        self.assertIn("^#CHROM", line)
        self.assertGreater(hits, 0)

    def test_no_header_read_is_capped_at_a_line_count(self):
        for path in CODE:
            for i, line in _code_lines(path):
                with self.subTest(file=f"{os.path.basename(path)}:{i}"):
                    self.assertNotRegex(line, r"head\s+-(n\s*)?\d{3,}")


if __name__ == "__main__":
    unittest.main()
