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
    reads the region distance;
  * clumping_ld_rule (critical | half_decay | r2_threshold, default critical) picks the
    LD-decay column an auto_* distance is read from; GWAS inherits it; a leftover
    clumping_r2_threshold is accepted only at LDdecay.r2_threshold, never picks the
    rule, and says so on stderr;
  * auto_per_chromosome needs LDdecay.scope per_chromosome or both, auto_genome_wide
    needs genome_wide or both (parse-time errors);
  * the LD-decay group reaches every script as one shell word ("Tel Aviv").
"""
import os
import re
import unittest

import _smk_harness as H

DOWNSTREAM = os.path.join(H.REPO_ROOT, "workflow", "rules", "_assoc_downstream.smk")
SUMMARY = os.path.join(H.REPO_ROOT, "workflow", "rules", "summary.smk")


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


class TestLdRule(unittest.TestCase):
    def test_default_rule_is_critical(self):
        ns = load()
        self.assertEqual(source(ns)["clumping_ld_rule"], "critical")
        self.assertEqual(ns["CLUMPING_LD_RULE"], "critical")

    def test_rule_is_read_and_inherited_by_gwas(self):
        ns = load(GEA__clumping_ld_rule="half_decay")
        self.assertEqual(source(ns)["clumping_ld_rule"], "half_decay")
        self.assertEqual(ns["PHENO_CLUMPING_LD_RULE"], "half_decay")

    def test_gwas_can_set_its_own_rule(self):
        ns = load(GEA__clumping_ld_rule="half_decay", GWAS__clumping_ld_rule="r2_threshold")
        self.assertEqual(ns["CLUMPING_LD_RULE"], "half_decay")
        self.assertEqual(ns["PHENO_CLUMPING_LD_RULE"], "r2_threshold")

    def test_an_unknown_rule_fails_at_parse_time(self):
        with self.assertRaises(Exception):
            load(GEA__clumping_ld_rule="extent")

    def test_old_r2_threshold_is_accepted_only_at_the_ld_decay_level(self):
        ns = load(GEA__clumping_r2_threshold=0.2)
        self.assertEqual(ns["CLUMPING_LD_RULE"], "critical")
        with self.assertRaises(ValueError):
            load(GEA__clumping_r2_threshold=0.3)
        ns = load(GEA__clumping_r2_threshold=0.3, LDdecay__r2_threshold=0.3)
        self.assertEqual(ns["LD_DECAY_R2_THRESHOLD"], 0.3)

    def test_old_r2_threshold_never_picks_the_rule_and_says_so(self):
        # the user's rule is critical unless clumping_ld_rule says otherwise; the old key
        # used to warn through a hidden DeprecationWarning only
        import contextlib
        import io
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            ns = load(GEA__clumping_ld_rule="half_decay", GWAS__clumping_r2_threshold=0.2)
        self.assertEqual(ns["PHENO_CLUMPING_LD_RULE"], "half_decay")
        self.assertIn("no longer picks the LD-decay rule", err.getvalue())

    def test_per_chromosome_distance_needs_per_chromosome_ld_rows(self):
        with self.assertRaises(ValueError):
            load(GEA__snp_clumping_distance="auto_per_chromosome", LDdecay__scope="genome_wide")
        ns = load(GEA__snp_clumping_distance="auto_per_chromosome", LDdecay__scope="both")
        self.assertEqual(source(ns)["clumping_distance_mode"], "auto_per_chromosome")

    def test_genome_wide_distance_needs_genome_wide_ld_rows(self):
        with self.assertRaisesRegex(ValueError, "auto_genome_wide needs genome-wide"):
            load(GEA__snp_clumping_distance="auto_genome_wide", LDdecay__scope="per_chromosome")
        ns = load(GEA__snp_clumping_distance="auto_per_chromosome", LDdecay__scope="per_chromosome")
        self.assertEqual(source(ns)["clumping_distance_mode"], "auto_per_chromosome")


class TestLdDecayKeys(unittest.TestCase):
    def test_defaults(self):
        ns = load()
        self.assertEqual(ns["LD_DECAY_MAX_PAIRS"], 2e7)
        self.assertEqual(ns["LD_DECAY_R2_THRESHOLD"], 0.2)
        self.assertEqual(ns["LD_DECAY_SEED"], 1)
        self.assertTrue(ns["LD_DECAY_TAG"].startswith("plink_"))

    def test_invalid_values_fail_at_parse_time(self):
        for key, value in (("LDdecay__max_pairs", 10), ("LDdecay__r2_threshold", 1.5),
                           ("LDdecay__background_markers_per_chr", 1),
                           ("LDdecay__max_distance", 0)):
            with self.subTest(key=key):
                with self.assertRaises(ValueError):
                    load(**{key: value})


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
                self.assertIn("{params.ld_decay_path} {params.ld_rule}", body)
                self.assertIn('_src(wc.source, "clumping_ld_rule")', body)
                self.assertNotIn("region_distance", body)

    def test_only_the_region_rule_reads_the_region_distance(self):
        body = self.rules["assoc_create_regions"]
        self.assertIn('_src(wc.source, "region_distance")', body)
        self.assertIn('ld_decay_input(_src(wc.source, "region_distance_mode"))', body)
        self.assertNotIn('_src(wc.source, "clumping_distance")', body)

    def test_the_ld_decay_group_is_one_shell_word(self):
        # a group named after a site with a space ("Tel Aviv") used to split in two
        # and shift every later argument
        for path in (DOWNSTREAM, SUMMARY):
            with self.subTest(file=os.path.basename(path)):
                with open(path) as handle:
                    text = handle.read()
                self.assertIn('"{params.ld_decay_group}"', text)
                self.assertNotRegex(text, r'(?<!")\{params\.ld_decay_group\}')


if __name__ == "__main__":
    unittest.main()
