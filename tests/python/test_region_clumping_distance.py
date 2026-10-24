"""GEA/GWAS.snp_clumping_distance and the advanced region_clumping_distance.

Pinned:
  * snp_clumping_distance is ONE distance for the agreement between methods, the
    cross-trait overlap annotation and the region clumping; unset
    region_clumping_distance leaves the regions on it (an existing config is unchanged);
  * region_clumping_distance replaces it for the region clumping only, as bp or auto_*;
  * a value that is neither bp nor an auto_* mode, or a negative one, fails at parse
    time (fail loudly);
  * GWAS inherits both keys from GEA unless it sets its own;
  * wiring: the agreement, the significant-SNP and WZA-window rules pass the LD-decay
    arguments that resolve an auto_* spec (finding of 2026-10-13: the agreement step
    used to receive the bare string and turn it into NA), and only the region rule
    reads the region distance.
"""
import os
import re
import unittest

import _smk_harness as H

DOWNSTREAM = os.path.join(H.REPO_ROOT, "workflow", "rules", "_assoc_downstream.smk")


def load(**overrides):
    return H.load(H.make_config(**overrides))


def source(ns, name="GEA"):
    return ns["ASSOC_SOURCES"][name]


class TestResolution(unittest.TestCase):
    def test_unset_region_distance_follows_the_shared_distance(self):
        src = source(load(GEA__snp_clumping_distance=150000))
        self.assertEqual(src["clumping_distance"], 150000)
        self.assertEqual(src["region_distance"], 150000)
        self.assertEqual(src["region_distance_mode"], "fixed")

    def test_unset_region_distance_follows_an_auto_shared_distance(self):
        src = source(load(GEA__snp_clumping_distance="auto_genome_wide"))
        self.assertEqual(src["clumping_distance_mode"], "auto_genome_wide")
        self.assertEqual(src["region_distance"], "auto_genome_wide")
        self.assertEqual(src["region_distance_mode"], "auto_genome_wide")

    def test_blank_region_distance_is_unset(self):
        for value in (None, "", "null"):
            with self.subTest(value=value):
                src = source(load(GEA__snp_clumping_distance=150000,
                                  GEA__region_clumping_distance=value))
                self.assertEqual(src["region_distance"], 150000)

    def test_region_distance_overrides_the_regions_only(self):
        src = source(load(GEA__snp_clumping_distance=5000,
                          GEA__region_clumping_distance=200000))
        self.assertEqual(src["clumping_distance"], 5000)
        self.assertEqual(src["clumping_distance_mode"], "fixed")
        self.assertEqual(src["region_distance"], 200000)
        self.assertEqual(src["region_distance_mode"], "fixed")

    def test_region_distance_may_be_auto_while_agreement_is_fixed(self):
        src = source(load(GEA__snp_clumping_distance=5000,
                          GEA__region_clumping_distance="auto_per_chromosome"))
        self.assertEqual(src["clumping_distance_mode"], "fixed")
        self.assertEqual(src["region_distance_mode"], "auto_per_chromosome")

    def test_invalid_distances_raise(self):
        for key in ("snp_clumping_distance", "region_clumping_distance"):
            for value in ("wide", -5):
                with self.subTest(key=key, value=value):
                    with self.assertRaisesRegex(ValueError, key):
                        load(**{f"GEA__{key}": value})

    def test_summary_reports_the_region_distance(self):
        ns = load(GEA__snp_clumping_distance=5000, GEA__region_clumping_distance=200000)
        self.assertEqual(ns["REGION_DISTANCE"], 200000)
        self.assertEqual(ns["CLUMPING_DISTANCE"], 5000)


class TestWiring(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(DOWNSTREAM) as handle:
            text = handle.read()
        cls.rules = {m.group(1): m.group(2) for m in re.finditer(
            r"^    rule (\w+):\n(.*?)(?=^    rule |\Z)", text, re.S | re.M)}

    def test_agreement_and_annotation_rules_resolve_auto_from_ld_decay(self):
        for name in ("assoc_find_sig_snps", "assoc_combine_selected_snps",
                     "assoc_wza_sig_windows"):
            with self.subTest(rule=name):
                body = self.rules[name]
                self.assertIn('_src(wc.source, "clumping_distance")', body)
                self.assertIn('ld_decay_input(_src(wc.source, "clumping_distance_mode"))', body)
                self.assertIn("{params.ld_decay_path} {params.r2_threshold}", body)
                self.assertNotIn("region_distance", body)

    def test_only_the_region_rule_reads_the_region_distance(self):
        body = self.rules["assoc_create_regions"]
        self.assertIn('_src(wc.source, "region_distance")', body)
        self.assertIn('ld_decay_input(_src(wc.source, "region_distance_mode"))', body)
        self.assertNotIn('_src(wc.source, "clumping_distance")', body)


if __name__ == "__main__":
    unittest.main()
