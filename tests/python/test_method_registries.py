"""workflow/methods/ — only the two places where CODE, not data, can be wrong.

gwas.py is NOT a dict literal: it is a comprehension filtering GEA_METHODS by a
flag (:9-13), so flipping `supports_phenotypes` on one entry silently changes
which methods mode=gwas runs. _gapit() (gea.py:51-68) is a factory that produces
8 of the 11 GEA entries.

Deliberately NOT tested: the GEA_METHODS schema sweep and maladaptation.py's
contents. Those are plain dict literals — a test restating them catches nothing
(edit the dict, edit the test) and CLAUDE.md's "avoid redundancy aggressively"
cuts against it.

workflow/rules/*.smk stays out of scope entirely: common.smk:9 needs the
Snakemake-injected `workflow` global and :66-78 reads `config` and touches the
filesystem, both at import time.
"""
import unittest

import _support

gea = _support.import_methods("gea")
gwas = _support.import_methods("gwas")

GAPIT_MODELS = {"BLINK", "FarmCPU", "MLM", "MLMM", "GLM", "CMLM", "ECMLM", "SUPER"}


class TestGwasDerivation(unittest.TestCase):
    """GWAS_METHODS = {m for m in GEA_METHODS if m.supports_phenotypes}."""

    def test_is_derived_from_gea_methods_not_written_out(self):
        for name, cfg in gwas.GWAS_METHODS.items():
            self.assertIs(cfg, gea.GEA_METHODS[name],
                          f"{name} is a copy, not the GEA_METHODS entry")

    def test_lfmm_and_rda_are_excluded(self):
        # gea.py:95 and :120 set supports_phenotypes False; :118-119 says in so
        # many words that this is what keeps RDA out of GWAS.
        self.assertNotIn("LFMM", gwas.GWAS_METHODS)
        self.assertNotIn("RDA", gwas.GWAS_METHODS)

    def test_emmax_and_every_gapit_model_are_included(self):
        self.assertIn("EMMAX", gwas.GWAS_METHODS)
        self.assertTrue(GAPIT_MODELS.issubset(set(gwas.GWAS_METHODS)))

    def test_exactly_the_flagged_methods_and_no_others(self):
        expected = {n for n, c in gea.GEA_METHODS.items() if c.get("supports_phenotypes")}
        self.assertEqual(set(gwas.GWAS_METHODS), expected)
        self.assertEqual(len(gwas.GWAS_METHODS), 9)

    def test_a_method_without_the_flag_would_be_dropped(self):
        # Proves the filter is a filter: nothing lacking the key survives.
        unflagged = [n for n, c in gea.GEA_METHODS.items()
                     if not c.get("supports_phenotypes", False)]
        self.assertTrue(unflagged)
        for name in unflagged:
            self.assertNotIn(name, gwas.GWAS_METHODS)


class TestGapitFactory(unittest.TestCase):

    def test_produces_a_complete_entry(self):
        entry = gea._gapit("BLINK")
        self.assertEqual(entry["engine"], "gapit")
        self.assertEqual(entry["script"], "scripts/gapit.R")
        self.assertEqual(entry["gapit_model"], "BLINK")
        self.assertTrue(entry["supports_phenotypes"])
        self.assertTrue(entry["supports_drop"])
        self.assertFalse(entry["multivariate"])
        self.assertIsNone(entry["pseudo_trait"])
        self.assertEqual(entry["significance_family"], "univariate_pvalue")

    def test_the_model_name_is_the_only_thing_that_varies(self):
        a, b = gea._gapit("GLM"), gea._gapit("MLM")
        self.assertNotEqual(a.pop("gapit_model"), b.pop("gapit_model"))
        self.assertEqual(a, b)

    def test_n_pcs_defaults_to_the_k_best_sentinel(self):
        # Resolved to sNMF k_best at config-parse time, not here.
        self.assertEqual(gea._gapit("GLM")["params"]["n_pcs"]["default"],
                         gea.K_BEST_SENTINEL)

    def test_every_gapit_registry_entry_came_from_the_factory(self):
        for model in GAPIT_MODELS:
            self.assertEqual(gea.GEA_METHODS[model], gea._gapit(model))

    def test_returns_a_fresh_dict_each_call(self):
        # Shared mutable state here would let one method's params leak into
        # another's at config-parse time.
        a, b = gea._gapit("GLM"), gea._gapit("GLM")
        self.assertIsNot(a, b)
        self.assertIsNot(a["params"], b["params"])


class TestParamSpec(unittest.TestCase):

    def test_p_merges_keywords_into_the_spec(self):
        self.assertEqual(gea.P("int", 3, min=0, max=20, help="h"),
                         {"type": "int", "default": 3, "min": 0, "max": 20, "help": "h"})

    def test_p_without_keywords_is_type_and_default_only(self):
        self.assertEqual(gea.P("bool", False), {"type": "bool", "default": False})


if __name__ == "__main__":
    unittest.main()
