"""EMMAX must read the sNMF-imputed genotype matrix, like every other GEA method.

WHY THIS IS A SEPARATE TEST
---------------------------
Until 2026-10-04 every EMMAX TPED was written from a NON-imputed VCF
(``plink --recode12 transpose --output-missing-genotype 0`` on ``W['vcf_filt']`` /
``W['vcf_ld']``), so each missing call reached ``emmax-intel64`` and
``emmax-kin-intel64`` as ``0 0`` and the binaries filled it with the per-SNP mean,
while LFMM, RDA and GAPIT read the ancestry-informed sNMF imputation
(``impute.R`` -> ``lfmm2vcf.R``). The rewire points every TPED chain at
``W['vcf_imp_full']`` / ``W['vcf_imp']`` instead.

Nothing at run time would notice a regression: a TPED built from the filtered VCF
is a perfectly valid EMMAX input, the scan exits 0, and the p-values merely come
from a different genotype matrix than the rest of the GEA. The contract is a
property of the DAG, so it is checked on the DAG: every rule that writes a TPED
(``--recode12``) or runs an EMMAX scan script must trace its ``vcf`` back, through
any sample-subsetting rules and the ``sNMF.imputed_maf`` filter
(``filter_imputed_maf_*``, which every chain must pass), to an ``lfmm2vcf_*`` rule —
the only producers of an imputed VCF.

The rule graph comes from parsing ``workflow/rules/*.smk`` as text and evaluating
each ``input:`` / ``params:`` expression in the namespace ``_smk_harness.load()``
builds — the same approach ``rule_outputs()`` takes for ``output:`` blocks. Rules
are read whatever their ``if`` gates say, so one config sees every chain.
"""
import ast
import glob
import os
import re
import unittest

import _smk_harness as H

# A rule whose shell writes a TPED, or runs an EMMAX scan off a VCF + TPED.
_EMMAX_SHELL = re.compile(r"--recode12|scripts/emmax\.R|scripts/emmax_phenotypes\.R")
# Keyword names that carry the genotype VCF (input `vcf`, transfer-guard params).
_VCF_KEY = re.compile(r"^(vcf|emmax_vcf_\w+)$")
IMPUTED_PRODUCERS = ("lfmm2vcf_full", "lfmm2vcf_ld")


def _keywords(text):
    """`name = expr` pairs of an input/params/output block, unevaluated."""
    try:
        call = ast.parse("_f(" + text + "\n)", mode="eval").body
    except SyntaxError:
        return {}, []
    return ({k.arg: ast.unparse(k.value) for k in call.keywords if k.arg},
            [ast.unparse(a) for a in call.args])


def _block(lines, i, rule_indent):
    """Text of the section starting at lines[i], up to the next section/rule."""
    section = H._SECTION.match(lines[i])
    indent = len(section.group(1))
    body = [lines[i][section.end():]]
    j = i + 1
    while j < len(lines):
        line = lines[j]
        if line.strip():
            here = H._indent(line)
            if here <= indent and (H._SECTION.match(line) or H._RULE.match(line)
                                   or line.lstrip().startswith('"""')):
                break
            if here <= rule_indent:
                break
        body.append(line)
        j += 1
    return "\n".join(body).strip().rstrip(","), j


def rule_graph(namespace):
    """[{name, smk, input, params, output, shell}] with expressions evaluated."""
    H._install_stubs(namespace)
    rules = []
    for path in sorted(glob.glob(os.path.join(H.REPO_ROOT, "workflow", "rules", "*.smk"))):
        if os.path.basename(path) == "common.smk":
            continue
        with open(path) as handle:
            lines = handle.read().splitlines()
        current = None
        i = 0
        while i < len(lines):
            match = H._RULE.match(lines[i])
            if match:
                envs = [H._local_assigns(lines, i, namespace, e)
                        for e in H._loop_env(lines, i, namespace)]
                current = {"name": match.group(1) or "<anonymous>",
                           "smk": os.path.basename(path), "indent": H._indent(lines[i]),
                           "env": envs[0] if envs else None, "raw": {}, "shell": ""}
                rules.append(current)
                i += 1
                continue
            section = H._SECTION.match(lines[i])
            if section and current is not None and H._indent(lines[i]) > current["indent"]:
                text, j = _block(lines, i, current["indent"])
                current["raw"][section.group(2)] = text
                if section.group(2) == "name":
                    current["name"] = text.strip().strip("f").strip("\"'")
                i = j
                continue
            if current is not None and lines[i].strip() and H._indent(lines[i]) <= current["indent"]:
                current = None
            i += 1

    for rule in rules:
        env = dict(rule["env"] or {})
        if "name" in rule["raw"]:
            try:
                rule["name"] = eval(rule["raw"]["name"], namespace, env)
            except Exception:
                pass
        rule["shell"] = rule["raw"].get("shell", "")
        for section in ("input", "params", "output"):
            named, positional = _keywords(rule["raw"].get(section, ""))
            rule[section] = {}
            for key, expr in list(named.items()) + [(f"_{n}", e) for n, e in enumerate(positional)]:
                try:
                    rule[section][key] = eval(expr, namespace, env)
                except Exception:
                    rule[section][key] = None
    return rules


def _producers(rules):
    by_path = {}
    for rule in rules:
        for value in rule["output"].values():
            for item in (value if isinstance(value, (list, tuple)) else [value]):
                if isinstance(item, str):
                    by_path[item] = rule
    return by_path


def trace(path, by_path, limit=10):
    """Rule names from `path` up its `vcf` input chain, nearest first."""
    chain = []
    while limit and path in by_path:
        rule = by_path[path]
        chain.append(rule["name"])
        if rule["name"] in IMPUTED_PRODUCERS:
            break
        path = rule["input"].get("vcf")
        limit -= 1
    return chain


class TestEmmaxReadsImputedGenotypes(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # GWAS Path A rules sit behind `if PHENO_MISSING != 'DROP'`; the text parse
        # sees them regardless, but their W keys only exist in a MEAN config.
        cls.configs = {
            "DROP": H.make_config(),
            "MEAN": H.make_config(GWAS__missing_strategy="MEAN"),
        }
        cls.graphs = {label: rule_graph(H.load(cfg)) for label, cfg in cls.configs.items()}

    def _emmax_vcfs(self, rules):
        found = []
        for rule in rules:
            # pregea_transfer_guard runs its own script and hands emmax.R the VCFs
            # as `emmax_vcf_*` params, so it is matched on those, not on its shell.
            if not (_EMMAX_SHELL.search(rule["shell"])
                    or any(k.startswith("emmax_vcf_") for k in rule["params"])):
                continue
            for section in ("input", "params"):
                for key, value in rule[section].items():
                    if _VCF_KEY.match(key) and isinstance(value, str) and value != "NULL":
                        found.append((rule["name"], key, value))
        return found

    def test_every_emmax_vcf_traces_to_an_lfmm2vcf_rule(self):
        for label, rules in self.graphs.items():
            by_path = _producers(rules)
            vcfs = self._emmax_vcfs(rules)
            names = {name for name, _k, _v in vcfs}
            # The five TPED rules, GEA/GWAS/preGEA scans and the transfer guard.
            for expected in ("tped_gea", "tped_gea_climate", "tped_gwas_trait",
                             "pregea_tped_pruned", "assoc_emmax_gea_trait",
                             "pregea_emmax_rung_trait", "pregea_transfer_guard"):
                self.assertIn(expected, names, f"[{label}] parser lost rule {expected}")
            if label == "MEAN":
                self.assertIn("tped_gwas", names)
                self.assertIn("gwas_a_emmax_trait", names)
            for name, key, path in vcfs:
                with self.subTest(config=label, rule=name, key=key):
                    self.assertNotIn("__placeholder_", path)
                    chain = trace(path, by_path)
                    self.assertTrue(
                        chain and chain[-1] in IMPUTED_PRODUCERS,
                        f"{name}.{key} = {path} does not come from the sNMF-imputed "
                        f"matrix (producer chain: {chain or 'none'})")
                    self.assertTrue(
                        any(r.startswith("filter_imputed_maf_") for r in chain),
                        f"{name}.{key} = {path} bypasses the sNMF.imputed_maf filter "
                        f"(producer chain: {chain})")

    def test_no_emmax_rule_reads_the_unimputed_vcfs(self):
        for label, cfg in self.configs.items():
            namespace = H.load(cfg)
            raw = {namespace["W"]["vcf_filt"], namespace["W"]["vcf_ld"]}
            for name, key, path in self._emmax_vcfs(self.graphs[label]):
                with self.subTest(config=label, rule=name, key=key):
                    self.assertNotIn(path, raw)

    def test_every_tped_rule_keeps_full_sample_ids_and_letter_contigs(self):
        """--double-id keeps an underscore-bearing sample name whole as FID and IID (plink
        splits WBDC_001 and refuses PG_B_001 without it, d8f924), and --output-chr MT keeps
        X/Y/MT as letters. A run on SIMDATA renamed PG_B_<name> with contig X gave the
        reference's GWAS tables exactly (2026-10-30, afb3cb, 0bf1fb); these flags are what
        made it so, and no SIMDATA run would notice one going missing."""
        rules = self.graphs["MEAN"]
        tped = [r for r in rules if "--recode12" in r["shell"]]
        self.assertEqual({r["name"] for r in tped},
                         {"tped_gea", "tped_gea_climate", "tped_gwas", "tped_gwas_trait",
                          "pregea_tped_pruned"})
        for rule in tped:
            with self.subTest(rule=rule["name"]):
                self.assertIn("--double-id", rule["shell"])
                self.assertIn("--output-chr MT", rule["shell"])


class TestGeaLessGwasHasTheImputedMatrix(unittest.TestCase):
    """GWAS EMMAX now needs the full-matrix imputation chain even when no GEA
    method is configured. add_pheno_association_paths() used to guard those keys
    with `'vcfsnp_full' not in W`, which never fires because common.smk pre-seeds
    every K-dependent key with a placeholder."""

    def test_full_matrix_paths_are_real_without_gea(self):
        # The sNMF.imputed_maf keys ride in the same block as lfmm_imp_full; the gwas
        # summary rule reads imputed_maf_report_full, and the structure one the LD
        # report, so all of them must resolve in a gwas_only project with no GEA.
        for label, cfg in (("standard", H.make_config(GEA__configs=[])),
                           ("gwas_only", H.make_config(GEA__configs=[],
                                                       Regime={"mode": "gwas_only"}))):
            namespace = H.load(cfg)
            for key in ("vcf_imp_full", "lfmm_imp_full", "snmf_full", "lfmm_full",
                        "vcfsnp_full", "gapit_gd", "assoc_kinship", "vcfsnp_imp_full",
                        "imputed_maf_report_full", "lfmm_imp_full_raw", "vcf_imp_full_raw",
                        "imputed_maf_report", "vcfsnp_imp", "lfmm_imp_raw", "vcf_imp_raw"):
                with self.subTest(regime=label, key=key):
                    self.assertNotIn("__placeholder_", namespace["W"][key])

    def test_gwas_only_emmax_chain_is_fully_resolved(self):
        """Walk every input upstream of each GWAS EMMAX rule in a gwas_only project
        with no GEA method: no `__placeholder_` path anywhere, and the genotypes
        still come out of lfmm2vcf_full. Both GWAS paths (MEAN = A, DROP = B)."""
        for strategy, rules_wanted in (("MEAN", ("tped_gwas", "gwas_a_emmax_trait")),
                                       ("DROP", ("tped_gwas_trait", "gwas_b_emmax_trait"))):
            namespace = H.load(H.make_config(
                GEA__configs=[], GWAS__missing_strategy=strategy,
                Regime={"mode": "gwas_only"}))
            rules = rule_graph(namespace)
            by_path = _producers(rules)
            for name in rules_wanted:
                rule = next(r for r in rules if r["name"] == name)
                with self.subTest(strategy=strategy, rule=name):
                    self.assertEqual(trace(rule["input"]["vcf"], by_path)[-1:],
                                     ["lfmm2vcf_full"])
                    seen, stack = set(), [rule]
                    while stack:
                        current = stack.pop()
                        if current["name"] in seen:
                            continue
                        seen.add(current["name"])
                        for key, value in current["input"].items():
                            for item in (value if isinstance(value, (list, tuple)) else [value]):
                                if not isinstance(item, str):
                                    continue
                                self.assertNotIn("__placeholder_", item,
                                                 f"{current['name']}.{key} upstream of {name}")
                                if item in by_path:
                                    stack.append(by_path[item])
                    self.assertIn("impute_full", seen)

    def test_gapit_only_gea_has_its_kinship(self):
        """gapit_gea_trait reads W['assoc_kinship'] (kinship_gea), which used to be
        set only when EMMAX was configured."""
        namespace = H.load(H.make_config(
            GEA__configs=[{"method": "BLINK", "adjust": "bonf", "threshold": "0.05"}]))
        for key in ("assoc_tped", "assoc_kinship", "vcf_imp_full"):
            with self.subTest(key=key):
                self.assertNotIn("__placeholder_", namespace["W"][key])


if __name__ == "__main__":
    unittest.main()
