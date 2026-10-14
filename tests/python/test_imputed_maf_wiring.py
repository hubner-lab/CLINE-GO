"""Every imputed genotype matrix must be labelled by the vcfsnp of the SAME SNP set.

WHY THIS IS A STATIC TEST
-------------------------
``sNMF.imputed_maf`` (filter_imputed_maf_ld / _full) drops SNPs from the imputed
``.lfmm`` / ``.vcf`` and writes a matching ``.vcfsnp``. LFMM, RDA, Gradient Forest,
the offsets and the preGEA ladders read the matrix and attach SNP ids to its
columns BY POSITION from a ``.vcfsnp``. A rule that pairs the filtered matrix with
the unfiltered ``W['vcfsnp']`` / ``W['vcfsnp_full']`` mislabels every column after
the first removed SNP — or, now that lfmm.R / gradient_forest_model.R check the
counts, stops mid-run.

A run cannot be relied on to catch it: at the default threshold SIMDATA removes
no SNP, so the two vcfsnp files are identical and a missed rewire passes. The
contract is a property of the DAG, so it is checked on the DAG, with the same
text-parse rule graph test_emmax_imputed_input.py builds.

Pinned:
  * a rule that reads an imputed lfmm (direct or climate-subset) reads only the
    filtered vcfsnp of that set — never the observed-call vcfsnp;
  * the raw imputed files (``*_raw``) are read only by lfmm2vcf_* and
    filter_imputed_maf_*, so no method can bypass the filter;
  * the climate-subset sources point at the filtered matrices;
  * gradient_forest_frequency is the one deliberate exception: it reads the
    OBSERVED matrix, so it keeps the observed vcfsnp.
"""
import ast
import os
import re
import unittest

import _smk_harness as H
from test_emmax_imputed_input import rule_graph

# Rules that label an imputed matrix's columns with a vcfsnp. Listed so a parser
# regression that loses one fails loudly instead of shrinking the check.
LABELLING_RULES = ("assoc_lfmm_gea_trait", "assoc_rda_gea", "pregea_lfmm_rung_trait",
                   "pregea_rda_setup", "pregea_transfer_guard", "gradient_forest_adaptive",
                   "gradient_forest_random", "geometric_offset", "rda_offset")
RAW_READERS = ("lfmm2vcf_ld", "lfmm2vcf_full", "filter_imputed_maf_ld", "filter_imputed_maf_full")


def _flat(values):
    for value in values:
        for item in (value if isinstance(value, (list, tuple)) else [value]):
            if isinstance(item, str):
                yield item


def _climate_sources(namespace):
    """Evaluate processing.smk's _LFMM_CLIMATE_SOURCES dict in the common.smk namespace."""
    path = os.path.join(H.REPO_ROOT, "workflow", "rules", "processing.smk")
    with open(path) as handle:
        text = handle.read()
    start = text.index("_LFMM_CLIMATE_SOURCES = {")
    tree = None
    end = text.index("}", start)
    while tree is None:
        try:
            tree = ast.parse(text[start:end + 1])
        except SyntaxError:
            end = text.index("}", end + 1)
    return eval(compile(ast.Expression(tree.body[0].value), path, "eval"), namespace)


class TestImputedMatrixLabels(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ns = H.load(H.make_config(sNMF__imputed_maf=0.1))
        cls.W = cls.ns["W"]
        cls.rules = rule_graph(cls.ns)
        W = cls.W
        # imputed lfmm path -> the vcfsnp that labels its columns
        cls.label_of = {
            W["lfmm_imp"]: W["vcfsnp_imp"],
            W["lfmm_imp_climate"]: W["vcfsnp_imp"],
            W["lfmm_imp_full"]: W["vcfsnp_imp_full"],
            W["lfmm_imp_full_climate"]: W["vcfsnp_imp_full"],
        }
        cls.vcfsnps = {W["vcfsnp"], W["vcfsnp_full"], W["vcfsnp_imp"], W["vcfsnp_imp_full"]}

    def _inputs(self, rule):
        return set(_flat(rule["input"].values()))

    def test_the_filtered_keys_are_distinct_files(self):
        # With the filter on, a key mix-up in common.smk would make the checks below vacuous.
        self.assertEqual(len(self.vcfsnps), 4)
        for key in self.label_of:
            self.assertNotIn("__placeholder_", key)

    def test_every_labelling_rule_was_parsed(self):
        names = {r["name"] for r in self.rules}
        for name in LABELLING_RULES + RAW_READERS + ("gradient_forest_frequency",):
            self.assertIn(name, names, f"parser lost rule {name}")

    def test_imputed_matrix_is_labelled_by_its_filtered_vcfsnp(self):
        checked = set()
        for rule in self.rules:
            inputs = self._inputs(rule)
            matrices = inputs & set(self.label_of)
            labels = inputs & self.vcfsnps
            if not (matrices and labels):
                continue
            checked.add(rule["name"])
            allowed = {self.label_of[m] for m in matrices}
            with self.subTest(rule=rule["name"]):
                self.assertTrue(labels <= allowed,
                                f"{rule['name']} reads {sorted(labels - allowed)} next to the "
                                f"filtered imputed matrix; expected one of {sorted(allowed)}")
        missing = set(LABELLING_RULES) - checked
        self.assertEqual(missing, set(), "labelling rules no longer pair a matrix with a vcfsnp")

    def test_raw_imputed_files_feed_only_the_filter(self):
        raw = {self.W[k] for k in ("lfmm_imp_raw", "vcf_imp_raw",
                                   "lfmm_imp_full_raw", "vcf_imp_full_raw")}
        readers = {r["name"] for r in self.rules if self._inputs(r) & raw}
        self.assertEqual(readers, set(RAW_READERS))

    def test_climate_subsets_are_cut_from_the_filtered_matrices(self):
        sources = _climate_sources(self.ns)
        self.assertEqual(sources["lfmm_imp"], self.W["lfmm_imp"])
        self.assertEqual(sources["lfmm_imp_full"], self.W["lfmm_imp_full"])

    def test_gf_frequency_keeps_the_observed_pairing(self):
        rule = next(r for r in self.rules if r["name"] == "gradient_forest_frequency")
        self.assertEqual(rule["input"]["lfmm"], self.W["lfmm_full_climate"])
        self.assertEqual(rule["input"]["vcfsnp"], self.W["vcfsnp_full"])

    def test_filter_rules_read_raw_and_write_the_filtered_set(self):
        W = self.W
        for level, raw_snp in (("", "vcfsnp"), ("_full", "vcfsnp_full")):
            name = "filter_imputed_maf_" + ("ld" if not level else "full")
            rule = next(r for r in self.rules if r["name"] == name)
            with self.subTest(rule=name):
                self.assertEqual(rule["input"]["lfmm_imp"], W[f"lfmm_imp{level}_raw"])
                self.assertEqual(rule["input"]["vcf"], W[f"vcf_imp{level}_raw"])
                self.assertEqual(rule["input"]["vcfsnp"], W[raw_snp])
                self.assertEqual(rule["output"]["lfmm"], W[f"lfmm_imp{level}"])
                self.assertEqual(rule["output"]["vcf"], W[f"vcf_imp{level}"])
                self.assertEqual(rule["output"]["vcfsnp"], W[f"vcfsnp_imp{level}"])
                self.assertTrue(re.search(r"filter_imputed_maf\.R", rule["shell"]))


if __name__ == "__main__":
    unittest.main()
