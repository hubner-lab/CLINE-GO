"""GEA.snp_sets / GWAS.snp_sets — SNP sets declared in config.

A declared set is the config-side twin of the Shiny GEA tab's "Save SNP set for
maladaptation". common.smk resolves each entry at parse time into exactly the
inputs the app's selection takes (trait::method pairs + one rule per cell);
scripts/build_snp_set.R then runs the app's own selection code on them. So the
correctness of a declared set is decided HERE, in the resolution, and pinned by:

  * TestBaseRule / TestOverrides — which rule each cell ends up with: the three
    base-rule cases, the RDA registry pin (only when a threshold is declared),
    and override precedence cell > trait > method;
  * TestPairs — traits x methods vs an explicit pair list, and RDA's single
    climate_multivariate column;
  * TestRejects — every configuration the parser must refuse at parse time,
    rather than let mode=gea burn compute and fail at the set build;
  * TestTargets — sets are built by mode=gea / mode=gwas only, by rules that can
    produce them, under a wildcard that can never match a Shiny-saved name;
  * TestMaladaptationResolution — `snp_sets: all` = store contents + declared
    names, an unbuilt declared set is an error, an empty set is skipped.
"""
import contextlib
import io
import os
import re
import unittest

import _smk_harness as H

GEA_PAIRS = ["bio_1::EMMAX", "bio_2::EMMAX", "bio_3::EMMAX",
             "bio_1::LFMM", "bio_2::LFMM", "bio_3::LFMM", "climate_multivariate::RDA"]


def gea_cfg(sets, **overrides):
    """BASE_CONFIG with GEA.snp_sets, and RDA's own rule moved off its registry
    pin (bonf 0.01) so the two are distinguishable."""
    cfg = H.make_config(**overrides)
    for c in cfg["GEA"]["configs"]:
        if c["method"] == "RDA":
            c["threshold"] = "0.05"
    cfg["GEA"]["snp_sets"] = sets
    return cfg


def decl(sets, name, **overrides):
    return H.load(gea_cfg(sets, **overrides))["SNP_SET_DECL"][name]


def rules(d):
    return {k: (v["type"], v["value"]) for k, v in d["cells"].items()}


class TestBaseRule(unittest.TestCase):

    def test_no_threshold_uses_each_methods_own_config_rule_rda_included(self):
        d = decl([{"name": "default"}], "gea.default")
        self.assertEqual(d["pairs"], GEA_PAIRS)
        r = rules(d)
        self.assertEqual(r["bio_1::EMMAX"], ("bonf", 0.05))
        # RDA's own GEA.configs entry (0.05 here), NOT the registry pin (0.01):
        # gea.default must match the static selected_snps.tsv rule set.
        self.assertEqual(r["climate_multivariate::RDA"], ("bonf", 0.05))

    def test_defaults_carry_strategy_regime_and_fixed_clumping(self):
        d = decl([{"name": "default"}], "gea.default")
        self.assertEqual((d["strategy"], d["regime"], d["origin"], d["source_module"]),
                         ("Union", "snp", "config", "GEA"))
        self.assertEqual(d["clumping_distance"], 1000000)   # GEA.snp_clumping_distance

    def test_auto_clumping_falls_back_to_the_apps_100kb(self):
        d = decl([{"name": "default"}], "gea.default", GEA__snp_clumping_distance="auto_genome_wide")
        self.assertEqual(d["clumping_distance"], 100000)

    def test_type_only_takes_the_types_default_value(self):
        r = rules(decl([{"name": "s", "threshold_type": "top"}], "gea.s"))
        self.assertEqual(r["bio_1::LFMM"], ("top", 100.0))

    def test_value_only_takes_the_first_configs_type(self):
        cfg_sets = [{"name": "s", "threshold_value": 0.2}]
        r = rules(decl(cfg_sets, "gea.s"))
        self.assertEqual(r["bio_2::EMMAX"], ("bonf", 0.2))
        cfg = gea_cfg(cfg_sets)
        cfg["GEA"]["configs"][0]["adjust"] = "qval"
        cfg["GEA"]["configs"][0]["threshold"] = "0.1"
        r = rules(H.load(cfg)["SNP_SET_DECL"]["gea.s"])
        self.assertEqual(r["bio_2::EMMAX"], ("qval", 0.2))

    def test_a_declared_threshold_pins_rda_to_its_registry_rule(self):
        r = rules(decl([{"name": "s", "threshold_type": "custom",
                         "threshold_value": 1e-6}], "gea.s"))
        self.assertEqual(r["bio_1::EMMAX"], ("custom", 1e-6))
        self.assertEqual(r["climate_multivariate::RDA"], ("bonf", 0.01))

    def test_yaml_string_numbers_are_accepted(self):
        # PyYAML reads 1e-05 (no dot) as a STRING; the parse must still take it.
        r = rules(decl([{"name": "s", "threshold_type": "custom",
                         "threshold_value": "1e-05"}], "gea.s"))
        self.assertEqual(r["bio_3::LFMM"], ("custom", 1e-5))


class TestOverrides(unittest.TestCase):

    def test_precedence_is_cell_over_trait_over_method(self):
        d = decl([{"name": "t", "threshold_type": "bonf", "threshold_value": 0.05,
                   "overrides": {"LFMM": {"value": 0.01},
                                 "bio_2": {"type": "custom", "value": 1e-5},
                                 "bio_2::EMMAX": {"type": "qval", "value": 0.1}}}],
                 "gea.t")
        r = rules(d)
        self.assertEqual(r["bio_1::LFMM"], ("bonf", 0.01))      # method
        self.assertEqual(r["bio_2::LFMM"], ("custom", 1e-5))    # trait beats method
        self.assertEqual(r["bio_2::EMMAX"], ("qval", 0.1))      # cell beats trait
        self.assertEqual(r["bio_1::EMMAX"], ("bonf", 0.05))     # master

    def test_value_only_override_keeps_the_type_it_overrides(self):
        r = rules(decl([{"name": "t", "overrides": {"EMMAX": {"value": 0.01}}}], "gea.t"))
        self.assertEqual(r["bio_1::EMMAX"], ("bonf", 0.01))

    def test_type_only_override_takes_the_types_default(self):
        r = rules(decl([{"name": "t", "overrides": {"bio_3": {"type": "qval"}}}], "gea.t"))
        self.assertEqual(r["bio_3::EMMAX"], ("qval", 0.05))

    def test_an_override_beats_the_rda_pin(self):
        r = rules(decl([{"name": "t", "threshold_value": 0.05,
                         "overrides": {"RDA": {"value": 0.001}}}], "gea.t"))
        self.assertEqual(r["climate_multivariate::RDA"], ("bonf", 0.001))


class TestPairs(unittest.TestCase):

    def test_traits_and_methods_select_their_rectangle(self):
        d = decl([{"name": "p", "methods": ["LFMM", "EMMAX"], "traits": ["bio_2"]}], "gea.p")
        self.assertEqual(d["pairs"], ["bio_2::LFMM", "bio_2::EMMAX"])

    def test_rda_is_selected_through_its_pseudo_trait(self):
        d = decl([{"name": "p", "traits": ["bio_1", "climate_multivariate"],
                   "methods": ["EMMAX", "RDA"]}], "gea.p")
        self.assertEqual(d["pairs"], ["bio_1::EMMAX", "climate_multivariate::RDA"])

    def test_explicit_pairs_pass_through_in_order(self):
        d = decl([{"name": "p", "pairs": ["bio_3::LFMM", "bio_1::EMMAX"]}], "gea.p")
        self.assertEqual(d["pairs"], ["bio_3::LFMM", "bio_1::EMMAX"])
        self.assertEqual(set(d["cells"]), {"bio_3::LFMM", "bio_1::EMMAX"})

    def test_gwas_sets_use_the_phenotype_columns(self):
        cfg = H.make_config()
        cfg["GWAS"]["snp_sets"] = [{"name": "all_traits"}]
        ns = H.load(cfg)
        d = ns["SNP_SET_DECL"]["gwas.all_traits"]
        self.assertEqual(d["source_module"], "GWAS")
        expected = [f"{t}::{m}" for m in ns["GWAS_CONFIGS"] for t in ns["PHENO_TRAITS"]]
        self.assertTrue(ns["PHENO_TRAITS"])
        self.assertEqual(d["pairs"], expected)

    def test_single_method_and_per_trait_strategies_are_accepted(self):
        sets = [{"name": "a", "strategy": "LFMM"},
                {"name": "b", "strategy": "Cross-method per-trait"},
                {"name": "c", "strategy": "Overlap"}]
        ns = H.load(gea_cfg(sets))
        self.assertEqual(ns["SNP_SET_DECL"]["gea.b"]["strategy"], "Cross-method per-trait")


class TestRejects(unittest.TestCase):

    def assertRefused(self, sets, fragment, **overrides):
        with self.assertRaises(ValueError) as cm:
            H.load(gea_cfg(sets, **overrides))
        self.assertIn(fragment, str(cm.exception))

    def test_names(self):
        self.assertRefused([{"name": "gea.x"}], "name must match")
        self.assertRefused([{"name": "a b"}], "name must match")
        self.assertRefused([{"name": "x"}, {"name": "x"}], "duplicate set name")
        self.assertRefused([{"threshold_value": 0.1}], "needs a 'name'")

    def test_shape(self):
        self.assertRefused({"name": "x"}, "must be a list")
        self.assertRefused([{"name": "x", "threshold": 0.1}], "unknown key(s) ['threshold']")

    def test_unknown_method_trait_or_override_key(self):
        self.assertRefused([{"name": "x", "methods": ["BLINK"]}], "not in GEA.configs")
        self.assertRefused([{"name": "x", "traits": ["bio_9"]}], "not emitted by")
        self.assertRefused([{"name": "x", "overrides": {"bio_9": {"value": 0.1}}}],
                           "override key 'bio_9'")
        self.assertRefused([{"name": "x", "traits": ["bio_1"], "methods": ["EMMAX"],
                             "overrides": {"bio_2::EMMAX": {"value": 0.1}}}],
                           "override key 'bio_2::EMMAX'")
        self.assertRefused([{"name": "x", "overrides": {"EMMAX": 0.1}}],
                           "an override is a map")
        self.assertRefused([{"name": "x", "overrides": {"EMMAX": {"vale": 0.1}}}],
                           "an override is a map")

    def test_pairs_are_exclusive_with_traits_and_methods(self):
        self.assertRefused([{"name": "x", "pairs": ["bio_1::EMMAX"], "traits": ["bio_1"]}],
                           "either 'pairs' or 'traits'/'methods'")
        self.assertRefused([{"name": "x", "pairs": ["bio_1::EMMAX"], "methods": ["EMMAX"]}],
                           "either 'pairs' or 'traits'/'methods'")

    def test_malformed_or_unknown_pairs(self):
        self.assertRefused([{"name": "x", "pairs": ["bio_1"]}], "is not trait::method")
        self.assertRefused([{"name": "x", "pairs": ["bio_1::RDA"]}], "emits no column")
        self.assertRefused([{"name": "x", "pairs": ["bio_1::LFMM", "bio_1::LFMM"]}],
                           "duplicate pairs")
        self.assertRefused([{"name": "x", "pairs": []}], "no trait::method pairs")

    def test_a_method_left_without_columns(self):
        self.assertRefused([{"name": "x", "methods": ["EMMAX", "RDA"], "traits": ["bio_1"]}],
                           "method RDA has no column")

    def test_thresholds(self):
        self.assertRefused([{"name": "x", "threshold_type": "fdr"}], "threshold_type must be")
        self.assertRefused([{"name": "x", "threshold_type": "bonf", "threshold_value": 2}],
                           "takes a probability")
        self.assertRefused([{"name": "x", "overrides": {"LFMM": {"type": "top", "value": 0.5}}}],
                           "takes a SNP count")

    def test_strategy(self):
        self.assertRefused([{"name": "x", "strategy": "Intersection"}], "strategy must be")
        # A method name is a strategy only if the set contains that method.
        self.assertRefused([{"name": "x", "methods": ["EMMAX"], "strategy": "LFMM"}],
                           "strategy must be")

    def test_clumping_distance(self):
        self.assertRefused([{"name": "x", "clumping_distance": "auto_genome_wide"}],
                           "clumping_distance must be a positive integer")
        self.assertRefused([{"name": "x", "clumping_distance": 0}],
                           "clumping_distance must be a positive integer")
        self.assertRefused([{"name": "x", "clumping_distance": 1.5}],
                           "clumping_distance must be a positive integer")

    def test_regime(self):
        self.assertRefused([{"name": "x", "regime": "wza"}], "not declarable yet")
        self.assertRefused([{"name": "x", "regime": "window"}], "regime must be 'snp'")

    def test_sets_without_configs(self):
        cfg = H.make_config()
        cfg["GWAS"]["configs"] = []
        cfg["GWAS"]["snp_sets"] = [{"name": "x"}]
        with self.assertRaises(ValueError) as cm:
            H.load(cfg)
        self.assertIn("GWAS.snp_sets is declared but GWAS.configs is empty", str(cm.exception))


class TestTargets(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cfg = gea_cfg([{"name": "default"}, {"name": "strict", "threshold_value": 0.01}])
        cfg["GWAS"]["snp_sets"] = [{"name": "pheno"}]
        cls.ns = H.load(cfg)
        cls.patterns = H.producible_patterns(cls.ns)

    def set_targets(self, mode):
        return [t for t in H.targets(self.ns, mode) if "/snp_sets/" in t]

    def test_gea_builds_the_gea_sets_and_their_records(self):
        got = self.set_targets("gea")
        self.assertEqual(sorted(os.path.relpath(t, self.ns["INTER"]) for t in got), [
            "snp_sets/gea.default/selected_snps.tsv", "snp_sets/gea.default/snp_set.yaml",
            "snp_sets/gea.strict/selected_snps.tsv", "snp_sets/gea.strict/snp_set.yaml"])

    def test_gwas_builds_only_the_gwas_sets(self):
        got = self.set_targets("gwas")
        self.assertEqual(len(got), 2)
        self.assertTrue(all("/gwas.pheno/" in t for t in got))

    def test_no_other_mode_builds_a_set(self):
        for mode in ("processing", "prestructure", "structure", "pregea", "gea_x_gwas"):
            self.assertEqual(self.set_targets(mode), [], mode)

    def test_every_set_target_is_producible_by_the_set_rules(self):
        for t in self.set_targets("gea") + self.set_targets("gwas"):
            rules_ = {rule for regex, _smk, rule in self.patterns if regex.match(t)}
            self.assertEqual(rules_, {"snp_set_build"}, t)
        decl_path = self.ns["snp_set_declaration"]("gea.strict")
        self.assertIn("snp_set_declare",
                      {rule for regex, _smk, rule in self.patterns if regex.match(decl_path)})

    def test_build_inputs_are_the_sets_method_tables(self):
        self.assertEqual(self.ns["snp_set_pvalue_inputs"]("gwas.pheno"),
                         [self.ns["pheno_pvalues"](m) for m in sorted(self.ns["GWAS_CONFIGS"])])

    def test_the_wildcard_never_matches_an_unprefixed_or_foreign_name(self):
        rx = self.ns["SNP_SET_REGEX"]
        for name in ("gea.strict", "gwas.pheno_1"):
            self.assertTrue(re.fullmatch(rx, name), name)
        for name in ("mysel", "GEA.strict", "gea.", "gea.a.b", "gea.a-b", "wza.x", "xgea.a"):
            self.assertIsNone(re.fullmatch(rx, name), name)


class TestMaladaptationResolution(unittest.TestCase):
    """resolve_active_snp_sets() against a per-test store in the harness sandbox."""

    _n = 0

    def ns_with(self, sets, snp_sets="all"):
        TestMaladaptationResolution._n += 1
        cfg = gea_cfg(sets, project_name=f"ssres{TestMaladaptationResolution._n}")
        cfg["Maladaptation"]["snp_sets"] = snp_sets
        return H.load(cfg)

    def store(self, ns, name, rows=1):
        d = H.save_snp_set(ns["PROJECT"], name)
        with open(os.path.join(d, "selected_snps.tsv"), "w") as fh:
            fh.write("SNPID\tchr\tpos\tmin_pvalue\n")
            for i in range(rows):
                fh.write(f"1:{i}\t1\t{i}\t1e-8\n")

    def resolve(self, ns):
        # The resolver reports skipped empty sets on stderr; keep the run quiet.
        with H.patched(), contextlib.redirect_stderr(io.StringIO()):
            return ns["resolve_active_snp_sets"]()

    def test_skipping_an_empty_set_is_reported(self):
        ns = self.ns_with([{"name": "strict"}, {"name": "empty"}])
        self.store(ns, "gea.strict")
        self.store(ns, "gea.empty", rows=0)
        err = io.StringIO()
        with H.patched(), contextlib.redirect_stderr(err):
            ns["resolve_active_snp_sets"]()
        self.assertIn("skips SNP set(s) with zero SNPs: ['gea.empty']", err.getvalue())

    def test_all_is_the_store_plus_the_declared_sets(self):
        ns = self.ns_with([{"name": "strict"}])
        self.store(ns, "mysel")
        self.store(ns, "gea.strict")
        self.assertEqual(self.resolve(ns), ["gea.strict", "mysel"])

    def test_an_unbuilt_declared_set_is_an_error(self):
        ns = self.ns_with([{"name": "strict"}])
        self.store(ns, "mysel")
        with self.assertRaises(ValueError) as cm:
            self.resolve(ns)
        self.assertIn("not built yet", str(cm.exception))
        self.assertIn("gea.strict", str(cm.exception))

    def test_an_empty_set_is_skipped_by_all(self):
        ns = self.ns_with([{"name": "strict"}, {"name": "empty"}])
        self.store(ns, "gea.strict")
        self.store(ns, "gea.empty", rows=0)
        self.assertEqual(self.resolve(ns), ["gea.strict"])

    def test_only_empty_sets_is_the_no_sets_error(self):
        ns = self.ns_with([{"name": "empty"}])
        self.store(ns, "gea.empty", rows=0)
        with self.assertRaises(ValueError) as cm:
            self.resolve(ns)
        self.assertIn("No non-empty curated SNP sets", str(cm.exception))

    def test_an_explicit_list_may_name_declared_sets(self):
        ns = self.ns_with([{"name": "strict"}], snp_sets=["gea.strict"])
        self.store(ns, "gea.strict")
        self.store(ns, "mysel")
        self.assertEqual(self.resolve(ns), ["gea.strict"])

    def test_an_explicit_list_refuses_unbuilt_missing_and_empty(self):
        ns = self.ns_with([{"name": "strict"}], snp_sets=["gea.strict"])
        with self.assertRaises(ValueError) as cm:
            self.resolve(ns)
        self.assertIn("not built yet", str(cm.exception))

        ns = self.ns_with([], snp_sets=["nosuch"])
        with self.assertRaises(ValueError) as cm:
            self.resolve(ns)
        self.assertIn("not found", str(cm.exception))

        ns = self.ns_with([{"name": "empty"}], snp_sets=["gea.empty"])
        self.store(ns, "gea.empty", rows=0)
        with self.assertRaises(ValueError) as cm:
            self.resolve(ns)
        self.assertIn("zero SNPs", str(cm.exception))

    def test_a_stale_directory_stays_in_all(self):
        # Declaration removed from config, directory still on disk: documented.
        ns = self.ns_with([])
        self.store(ns, "gea.removed")
        self.assertEqual(self.resolve(ns), ["gea.removed"])


if __name__ == "__main__":
    unittest.main()
