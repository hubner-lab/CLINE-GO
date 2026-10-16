"""Maladaptation.predictors: the climate predictors the offset models are fitted on.

Pinned:
  * unset / blank / 'auto' resolves to Climate.predictors, so an existing config is
    unchanged;
  * a comma string and a YAML list resolve to the same thing;
  * a predictor outside Climate.predictors, a duplicate, or an empty list fails at
    parse time (fail loudly -- the climate tables only carry Climate.predictors);
  * every maladaptation model / offset / plot rule reads MALA_PREDICTORS, while
    density_plot_future (climate characterization) stays on Climate.predictors.
    Static on the rule graph, because a rule left on PREDICTORS_SELECTED runs fine
    and silently fits the wrong predictor set.
"""
import unittest

import _smk_harness as H
from test_emmax_imputed_input import rule_graph

CLIMATE = "bio_1,bio_3,bio_5,bio_12"
MALA_RULES = ("gradient_forest_adaptive", "gradient_forest_random", "gradient_forest_frequency",
              "gf_sensitivity_check", "gradient_forest_offset", "geometric_offset", "rda_offset",
              "plot_gf_cumimp")


def resolve(value=None, unset=False):
    overrides = {"Climate__predictors": CLIMATE}
    if not unset:
        overrides["Maladaptation__predictors"] = value
    return H.load(H.make_config(**overrides))["MALA_PREDICTORS"]


class TestResolution(unittest.TestCase):
    def test_unset_is_climate_predictors(self):
        self.assertEqual(resolve(unset=True), CLIMATE)

    def test_auto_and_blank_are_climate_predictors(self):
        for value in ("auto", "", " auto ", None):
            with self.subTest(value=value):
                self.assertEqual(resolve(value), CLIMATE)

    def test_subset_string_and_list_agree(self):
        self.assertEqual(resolve("bio_1, bio_5"), "bio_1,bio_5")
        self.assertEqual(resolve(["bio_1", "bio_5"]), "bio_1,bio_5")
        self.assertEqual(resolve("bio_12"), "bio_12")

    def test_unknown_predictor_raises(self):
        with self.assertRaisesRegex(ValueError, "not in Climate.predictors"):
            resolve("bio_1,bio_7")

    def test_duplicate_raises(self):
        with self.assertRaisesRegex(ValueError, "more than once"):
            resolve("bio_1,bio_1")

    def test_empty_list_raises(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            resolve([])


class TestWiring(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ns = H.load(H.make_config(Climate__predictors=CLIMATE,
                                  Maladaptation__predictors="bio_1,bio_5"))
        cls.rules = {r["name"]: r for r in rule_graph(ns)}

    def test_offset_rules_read_mala_predictors(self):
        for name in MALA_RULES:
            with self.subTest(rule=name):
                self.assertIn(name, self.rules, "rule parser lost a maladaptation rule")
                params = self.rules[name]["raw"].get("params", "")
                self.assertRegex(params, r"predictors\s*=\s*MALA_PREDICTORS\b")
                self.assertNotIn("PREDICTORS_SELECTED", params)

    def test_no_other_maladaptation_rule_reads_climate_predictors(self):
        stray = [name for name, r in self.rules.items()
                 if r["smk"] == "maladaptation.smk" and name not in MALA_RULES
                 and name != "density_plot_future"
                 and "PREDICTORS_SELECTED" in r["raw"].get("params", "")]
        self.assertEqual(stray, [])

    def test_future_density_plot_keeps_climate_predictors(self):
        self.assertIn("PREDICTORS_SELECTED", self.rules["density_plot_future"]["raw"]["params"])


if __name__ == "__main__":
    unittest.main()
