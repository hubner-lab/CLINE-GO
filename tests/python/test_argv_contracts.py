"""Tests for scripts/lint_argv_contracts.py — the Snakemake<->Rscript argv contract.

Two layers, and they answer different questions:

  * UNIT tests over hand-built fixtures: does each check FIRE when the defect is
    present, and stay quiet when it is not? A linter nobody has seen fail is
    indistinguishable from one that matches nothing.
  * A REPO test: the real tree must be clean. This is the regression gate — it is
    the only thing in the suite that would notice an argument appended to a
    script but not to one of its call sites.

The unit layer includes a MUTATION test: it copies the real `workflow/rules/` and
`scripts/` trees, injects the three historically-observed regressions, and asserts
the linter turns red on each. That is what distinguishes "clean tree" from "check
never fires".
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

from _support import REPO_ROOT, load_script

lint = load_script("lint_argv_contracts")


def write(tmp, rel, text):
    path = os.path.join(tmp, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def mini_tree(tmp, script_body, rule_body):
    """A minimal repo: one argv-reading R script, one .smk rule invoking it."""
    write(tmp, "scripts/toy.R", script_body)
    write(tmp, "workflow/rules/toy.smk", rule_body)
    return lint.lint(tmp)


SCRIPT_3 = """#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
IN_TABLE = args[1]
THRESH   = as.numeric(args[2])
OUT_FILE = args[3]
"""

RULE_OK = '''
rule toy:
    input:
        table = "in.tsv"
    output:
        "out.tsv"
    params:
        thresh = 0.05
    log: "toy.log"
    shell:
        """
        Rscript /pipeline/scripts/toy.R \\
            {input.table} {params.thresh} {output} > {log} 2>&1
        """
'''


class TestScriptContract(unittest.TestCase):
    def test_reads_every_index(self):
        c = lint.ScriptContract("scripts/toy.R", SCRIPT_3)
        self.assertEqual(c.max_index, 3)
        self.assertEqual(c.required, 3)
        self.assertEqual(c.missing, [])

    def test_guarded_position_is_optional(self):
        body = SCRIPT_3 + 'REGIME = if (length(args) >= 4) args[4] else "snp"\n'
        c = lint.ScriptContract("scripts/toy.R", body)
        self.assertEqual(c.max_index, 4)
        self.assertEqual(c.required, 3, "a length(args)-guarded slot is optional")
        self.assertEqual(c.optional, [4])

    def test_unconditional_read_above_a_guard_is_still_required(self):
        # The regression that matters: a new mandatory argument appended past an
        # existing optional one. A guard-minimum heuristic would miss this.
        body = (SCRIPT_3
                + 'REGIME = if (length(args) >= 4) args[4] else "snp"\n'
                + "WHITELIST = args[5]\n")
        c = lint.ScriptContract("scripts/toy.R", body)
        self.assertEqual(c.required, 5)

    def test_sparse_argv_is_detected(self):
        body = SCRIPT_3 + "EXTRA = args[5]\n"
        c = lint.ScriptContract("scripts/toy.R", body)
        self.assertEqual(c.missing, [4])

    def test_output_slot_from_variable_name(self):
        c = lint.ScriptContract("scripts/toy.R", SCRIPT_3)
        self.assertEqual(c.out_slots, {3})
        self.assertEqual(c.in_slots, {1, 2})

    def test_slot_bound_to_both_roles_is_ambiguous_not_guessed(self):
        # stage_custom_climate.R binds args[7] to OUT_RASTER in one mode branch
        # and PRESENT_ALL in the other. Neither classification is safe.
        body = (
            "args = commandArgs(trailingOnly = TRUE)\n"
            "MODE = args[1]\n"
            "if (MODE == 'a') {\n"
            "  OUT_RASTER <- args[2]\n"
            "} else {\n"
            "  PRESENT_ALL <- args[2]\n"
            "}\n"
        )
        c = lint.ScriptContract("scripts/toy.R", body)
        self.assertNotIn(2, c.out_slots)
        self.assertNotIn(2, c.in_slots)

    def test_output_slot_behind_a_guard_still_counts(self):
        # geometric_offset.R:61 — `OUT_DIAGNOSTICS <- if (length(args) >= 19 &&
        # nzchar(args[19])) args[19] else NA`. The name is on the line, not
        # adjacent to args[19].
        body = (SCRIPT_3 + "OUT_DIAG <- if (length(args) >= 4 && "
                           "nzchar(args[4])) args[4] else NA\n")
        c = lint.ScriptContract("scripts/toy.R", body)
        self.assertIn(4, c.out_slots)


class TestShellExtraction(unittest.TestCase):
    """Shell-level quotes must survive; Python-level delimiters must not."""

    def _tokens(self, shell_section):
        rule = lint.Rule("r", "f.smk", 1)
        rule.sections["shell"] = shell_section.splitlines()
        cmd = lint.shell_text(rule)
        i = cmd.index("toy.R") + len("toy.R")
        return lint.tokenize_invocation(cmd[i:])

    def test_triple_quoted_block(self):
        toks = self._tokens(
            '            """\n'
            "            Rscript /pipeline/scripts/toy.R \\\n"
            "                {input.a} {input.b} {output} > {log} 2>&1\n"
            '            """\n')
        self.assertEqual(toks, ["{input.a}", "{input.b}", "{output}"])

    def test_shell_double_quotes_are_preserved(self):
        # This is the load-bearing case: the ONLY thing protecting a space-joined
        # file list from shifting every later argument is this pair of quotes.
        toks = self._tokens(
            "        'Rscript /pipeline/scripts/toy.R '\n"
            "        '\"{params.files_str}\" {output.pvals} {output.qvals} "
            "> {log} 2>&1'\n")
        self.assertEqual(toks[0], '"{params.files_str}"')
        self.assertEqual(len(toks), 3)

    def test_quoted_literal_with_a_space_is_one_token(self):
        toks = self._tokens(
            '            """\n'
            "            Rscript /pipeline/scripts/toy.R \\\n"
            '                {input.a} 1 "Genetic Offset" {output} > {log} 2>&1\n'
            '            """\n')
        self.assertEqual(toks, ["{input.a}", "1", '"Genetic Offset"', "{output}"])

    def test_redirect_ends_the_argument_list(self):
        toks = self._tokens(
            '            """\n'
            "            Rscript /pipeline/scripts/toy.R {input.a} > {log} 2>&1\n"
            "            touch {output}\n"
            '            """\n')
        self.assertEqual(toks, ["{input.a}"])


class TestCallSiteChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp)

    def _checks(self, script_body, rule_body):
        violations, _ = mini_tree(self.tmp, script_body, rule_body)
        return sorted({v.check for v in violations if v.severity == "error"})

    def test_matching_call_site_is_clean(self):
        self.assertEqual(self._checks(SCRIPT_3, RULE_OK), [])

    def test_missing_argument_is_an_arity_mismatch(self):
        # Dropping a middle argument trips BOTH checks, and that is the point:
        # {output} slides down into an input slot at the same time.
        self.assertEqual(
            self._checks(SCRIPT_3, RULE_OK.replace("{params.thresh} ", "")),
            ["arity_mismatch", "slot_role_mismatch"])

    def test_extra_argument_is_an_arity_mismatch(self):
        self.assertEqual(
            self._checks(SCRIPT_3,
                         RULE_OK.replace("{output}", "{output} extra")),
            ["arity_mismatch"])

    def test_unquoted_space_joined_params_is_an_error(self):
        rule = '''
rule toy:
    input:
        tables = ["a.tsv", "b.tsv"]
    output:
        "out.tsv"
    params:
        files_str = lambda wc, input: " ".join(input.tables),
        thresh = 0.05
    log: "toy.log"
    shell:
        """
        Rscript /pipeline/scripts/toy.R \\
            {params.files_str} {params.thresh} {output} > {log} 2>&1
        """
'''
        self.assertEqual(self._checks(SCRIPT_3, rule), ["unquoted_multi_token"])

    def test_quoting_the_same_param_clears_it(self):
        rule = '''
rule toy:
    input:
        tables = ["a.tsv", "b.tsv"]
    output:
        "out.tsv"
    params:
        files_str = lambda wc, input: " ".join(input.tables),
        thresh = 0.05
    log: "toy.log"
    shell:
        """
        Rscript /pipeline/scripts/toy.R \\
            "{params.files_str}" {params.thresh} {output} > {log} 2>&1
        """
'''
        self.assertEqual(self._checks(SCRIPT_3, rule), [])

    def test_bare_output_in_a_multi_output_rule_is_an_error(self):
        rule = '''
rule toy:
    input:
        table = "in.tsv"
    output:
        a = "out_a.tsv",
        b = "out_b.tsv"
    params:
        thresh = 0.05
    log: "toy.log"
    shell:
        """
        Rscript /pipeline/scripts/toy.R \\
            {input.table} {params.thresh} {output} > {log} 2>&1
        """
'''
        # {output} is 2 tokens, so the script would get 4 arguments for 3 slots.
        self.assertIn("arity_mismatch", self._checks(SCRIPT_3, rule))

    def test_output_landing_in_an_input_slot_is_an_error(self):
        rule = RULE_OK.replace("{input.table} {params.thresh} {output}",
                               "{output} {params.thresh} {input.table}")
        self.assertEqual(self._checks(SCRIPT_3, rule), ["slot_role_mismatch"])


class TestRealRepository(unittest.TestCase):
    """The regression gate, plus proof it can go red."""

    def test_repository_has_no_argv_contract_errors(self):
        violations, stats = lint.lint(REPO_ROOT)
        errors = [str(v) for v in violations if v.severity == "error"]
        self.assertEqual(errors, [], "\n".join(errors))
        self.assertGreater(stats["scripts"], 60)

    # Call sites the linter does not parse by design, with their multiplicity:
    # the one write_summary.R call reached through shell() inside a run: block
    # (write_summary.R computes its indices per mode anyway). Its eight shell:
    # call sites ARE parsed and must stay so.
    UNPARSED_BY_DESIGN = {("summary.smk", "write_summary.R"): 1}

    def test_every_rscript_call_site_is_parsed(self):
        """Every textual `Rscript .../scripts/X.R` in the rules must be linted.

        A floor ("more than 80 sites") cannot fail when the parser silently drops
        a couple of rules — which is exactly what an exact +4 section indent did to
        both rules in gea_x_gwas.smk. Equality can.
        """
        import re
        pat = re.compile(r"Rscript\s+(?:\S*/)?scripts/([A-Za-z0-9_]+\.R)")
        rdir = os.path.join(REPO_ROOT, "workflow", "rules")
        textual, parsed = [], []
        for fn in sorted(os.listdir(rdir)):
            if not fn.endswith(".smk"):
                continue
            with open(os.path.join(rdir, fn)) as fh:
                for line in fh:
                    if line.lstrip().startswith("#"):
                        continue
                    textual += [(fn, m) for m in pat.findall(line)]
            for rule in lint.parse_rules(os.path.join(rdir, fn), fn):
                parsed += [(fn, m.group(1))
                           for m in lint._RSCRIPT.finditer(lint.shell_text(rule))]
        # Multisets, not sets: a script called from two rules in one file is
        # "seen" by a set as soon as either rule parses.
        from collections import Counter
        missed = Counter(textual) - Counter(parsed)
        self.assertEqual(dict(missed), self.UNPARSED_BY_DESIGN,
                         "Rscript call sites the parser never saw")
        self.assertEqual(Counter(parsed) - Counter(textual), Counter())

    def test_cli_exits_zero_on_the_real_tree(self):
        p = subprocess.run(
            [sys.executable, os.path.join(REPO_ROOT, "scripts",
                                          "lint_argv_contracts.py"),
             REPO_ROOT, "--quiet"],
            capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_injected_regressions_are_caught(self):
        """Each mutation is one real historical failure mode, re-injected."""
        mutations = [
            # 1. A new mandatory argument appended to a script; call sites stale.
            ("scripts/find_genes_around_regions.R",
             "args = commandArgs(trailingOnly = TRUE)",
             "args = commandArgs(trailingOnly = TRUE)\nBIOTYPE_WL = args[9]",
             "arity_mismatch"),
            # 2. The shell quotes lost from a space-joined file list — the shape
            #    that once made combine_pheno_pvalues.R overwrite an input
            #    (it now hard-stops on length(args); most scripts do not).
            ("workflow/rules/gwas.smk",
             '"{params.files_str}"', "{params.files_str}",
             "unquoted_multi_token"),
            # 3. Two positions transposed, so the output path lands in an input
            #    slot. Arity cannot see this; the slot role can.
            ("workflow/rules/gea.smk",
             "{params.inter_dir} {input.metadata} \\",
             "{params.inter_dir} {output} \\",
             "slot_role_mismatch"),
            # 4. A token dropped in a rule whose sections are indented by 2, not
            #    4 — invisible until the parser stopped assuming +4.
            ("workflow/rules/gea_x_gwas.smk",
             "{params.kbest} {params.plot_dir}", "{params.plot_dir}",
             "arity_mismatch"),
        ]
        for rel, old, new, expected in mutations:
            with self.subTest(mutation=expected, file=rel):
                tmp = tempfile.mkdtemp()
                self.addCleanup(shutil.rmtree, tmp)
                for sub in ("scripts", "workflow"):
                    shutil.copytree(os.path.join(REPO_ROOT, sub),
                                    os.path.join(tmp, sub))
                path = os.path.join(tmp, rel)
                with open(path) as fh:
                    text = fh.read()
                self.assertIn(old, text, f"mutation anchor gone from {rel}")
                with open(path, "w") as fh:
                    fh.write(text.replace(old, new, 1))

                violations, _ = lint.lint(tmp)
                checks = {v.check for v in violations if v.severity == "error"}
                self.assertIn(expected, checks,
                              f"{rel}: linter stayed green on an injected "
                              f"{expected}")


if __name__ == "__main__":
    unittest.main()
