"""workflow/methods/ — the places where CODE, not data, can be wrong, plus the two
assertion families that are not value restatements.

gwas.py is NOT a dict literal: it is a comprehension filtering GEA_METHODS by a
flag (:9-13), so flipping `supports_phenotypes` on one entry silently changes
which methods mode=gwas runs. _gapit() (gea.py:51-68) is a factory that produces
8 of the 11 GEA entries.

Deliberately NOT tested: a GEA_METHODS / MALADAPTATION_METHODS schema sweep —
required-key presence, per-key types, allowed `engine` values. Those are plain dict
literals; a test restating them catches nothing (edit the dict, edit the test) and
CLAUDE.md's "avoid redundancy aggressively" cuts against it.

ADDED 2026-09-12, and the line between the two is worth stating because it is easy
to slide back across. Two assertion families about the same dicts DO have an
independent failure mode, so they are not restatements:

  * TestRegistryScriptsExist — every declared script path is checked against the
    FILESYSTEM. A renamed or deleted R script currently surfaces only as a
    Snakemake error mid-run, arbitrarily far downstream.
  * TestRegistryCrossFieldRules — a capability flag and the script path it implies
    must agree (builds_model without a model_script, and so on). These encode a
    RULE about the registry, not the registry's contents: a new entry that gets it
    wrong goes red without the test being edited.

workflow/rules/*.smk stays out of scope entirely: common.smk:9 needs the
Snakemake-injected `workflow` global and :66-78 reads `config` and touches the
filesystem, both at import time.
"""
import os
import unittest

import _support

gea = _support.import_methods("gea")
gwas = _support.import_methods("gwas")
malad = _support.import_methods("maladaptation")

SCRIPT_KEYS = ("script", "model_script", "offset_script",
               "cumimp_script", "importance_script")

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


class TestRegistryScriptsExist(unittest.TestCase):
    """Every declared script path must resolve to a real file.

    Not a restatement: the assertion is about the filesystem, not about the dict.
    Paths are repo-relative (e.g. "scripts/gapit.R").
    """

    def _check(self, registry, label):
        seen = 0
        for name, cfg in registry.items():
            for key in SCRIPT_KEYS:
                rel = cfg.get(key)
                if rel is None:
                    continue
                path = os.path.join(_support.REPO_ROOT, rel)
                self.assertTrue(os.path.isfile(path),
                                f"{label}[{name}][{key}] -> missing file {rel}")
                seen += 1
        # Guard against a vacuous pass if the registry is ever emptied or a key
        # is renamed out from under SCRIPT_KEYS.
        self.assertGreater(seen, 0, f"{label}: no script paths checked at all")

    def test_gea_method_scripts_exist(self):
        self._check(gea.GEA_METHODS, "GEA_METHODS")

    def test_maladaptation_method_scripts_exist(self):
        self._check(malad.MALADAPTATION_METHODS, "MALADAPTATION_METHODS")


class TestRegistryCrossFieldRules(unittest.TestCase):
    """Rules relating a capability flag to the script path it implies.

    Each encodes an invariant the rule factories depend on, so a new entry that
    violates one goes red without this file being touched.
    """

    def test_multivariate_implies_a_pseudo_trait_and_vice_versa(self):
        # gea.py's header: multivariate=True means ONE p-value column for the whole
        # predictor set, named by pseudo_trait. One without the other leaves the
        # threshold layer with no column name to dispatch on.
        for name, cfg in gea.GEA_METHODS.items():
            self.assertEqual(bool(cfg["multivariate"]),
                             cfg["pseudo_trait"] is not None,
                             f"GEA_METHODS[{name}]: multivariate and pseudo_trait disagree")

    def test_builds_model_false_implies_no_model_script(self):
        # builds_model=False is documented as "one-call: model + offset in a single
        # script"; a model_script alongside it would never be invoked.
        for name, cfg in malad.MALADAPTATION_METHODS.items():
            if not cfg["builds_model"]:
                self.assertIsNone(cfg["model_script"],
                                  f"MALADAPTATION_METHODS[{name}]: "
                                  "builds_model is False but a model_script is declared")

    def test_builds_model_true_requires_a_model_script(self):
        for name, cfg in malad.MALADAPTATION_METHODS.items():
            if cfg["builds_model"]:
                self.assertIsNotNone(cfg["model_script"],
                                     f"MALADAPTATION_METHODS[{name}]: "
                                     "builds_model is True but no model_script")

    def test_cumulative_importance_requires_a_cumimp_script(self):
        for name, cfg in malad.MALADAPTATION_METHODS.items():
            if cfg["supports_cumulative_importance"]:
                self.assertIsNotNone(cfg["cumimp_script"],
                                     f"MALADAPTATION_METHODS[{name}]: claims cumulative "
                                     "importance support with no cumimp_script")

    def test_every_maladaptation_method_computes_an_offset(self):
        # The offset raster is the module's entire output; a method without one
        # would register a target nothing can build.
        for name, cfg in malad.MALADAPTATION_METHODS.items():
            self.assertIsNotNone(cfg["offset_script"],
                                 f"MALADAPTATION_METHODS[{name}]: no offset_script")


if __name__ == "__main__":
    unittest.main()
