"""Every workflow file parses with Snakemake's own parser (9f4c4a, 2026-10-30).

_smk_harness.py exec()s common.smk as plain Python, so a line Snakemake reads as a keyword
(a variable named `rule`, a directive indented one step off) passed every test here and failed
only at `snakemake -n`. `snakemake -s <file> --print-compilation` runs Snakemake's parser on one
file and prints the Python it compiles to, without executing anything; `include:` is a runtime
call, so the Snakefile and each rule file are parsed one by one. A file Snakemake cannot parse
exits non-zero and prints nothing, so the exit status and a non-empty output are both asserted.

The snakemake executable is part of the image; its absence fails the test rather than skipping
it. The parser writes .snakemake/ into its working directory, so it runs in a temp dir.
"""
import glob
import os
import shutil
import subprocess
import tempfile
import unittest

from _support import REPO_ROOT

WORKFLOW_FILES = [os.path.join(REPO_ROOT, "Snakefile")] + sorted(
    glob.glob(os.path.join(REPO_ROOT, "workflow", "rules", "*.smk")))


def _compile(snakemake, path, cwd):
    return subprocess.run([snakemake, "-s", path, "--print-compilation"],
                          cwd=cwd, capture_output=True, text=True)


class SnakemakeParserTest(unittest.TestCase):
    def setUp(self):
        self.snakemake = shutil.which("snakemake")
        self.assertIsNotNone(self.snakemake, "snakemake is not on PATH: the image is wrong")

    def test_every_workflow_file_parses(self):
        self.assertGreater(len(WORKFLOW_FILES), 10)
        with tempfile.TemporaryDirectory() as cwd:
            for path in WORKFLOW_FILES:
                name = os.path.relpath(path, REPO_ROOT)
                with self.subTest(file=name):
                    res = _compile(self.snakemake, path, cwd)
                    self.assertEqual(res.returncode, 0, res.stderr[-3000:])
                    self.assertTrue(res.stdout.strip(), f"{name}: no compilation printed")
                    compile(res.stdout, name, "exec")

    def test_a_file_snakemake_cannot_parse_is_refused(self):
        # Negative controls for the two shapes above: the first is valid Python that only
        # Snakemake's parser refuses, the second a directive indented one step off.
        broken = {
            "rule_as_variable.smk": 'X = 1\nrule = "a"\n',
            "directive_misindented.smk": ('rule a:\n    output: "x"\n'
                                          '    shell: "touch {output}"\n  params: y = 1\n'),
        }
        with tempfile.TemporaryDirectory() as cwd:
            for name, text in broken.items():
                with self.subTest(file=name):
                    path = os.path.join(cwd, name)
                    with open(path, "w") as handle:
                        handle.write(text)
                    res = _compile(self.snakemake, path, cwd)
                    self.assertNotEqual(res.returncode, 0, res.stdout[-1000:])
                    self.assertEqual(res.stdout.strip(), "")


if __name__ == "__main__":
    unittest.main()
