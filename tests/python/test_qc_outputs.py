"""Release Phase 4c: Processing's QC tables come from the steps that decided them.

Three wiring facts the R tests cannot see, all properties of the rule graph:

  * 7ca55e  compute_depth_stats reads the RAW input VCF over the metadata keep-list. It read
            W['vcf_filt'], which plink writes with GT only, so every depth was -nan while the
            rule exited 0. vcftools takes the first field of a keep-list line as the ID, so the
            plink list (FID IID) must be cut to its IIDs, or vcftools keeps no sample.
  * 7e11b6  filter_vcf_steps replays filter_vcf's plink call: it must read the same VCF and the
            same keep-list in every configuration (depth filter, het filter, relatedness
            removal), or its step counts describe another filter. Both rules take the
            module-level names _FILTER_VCF_INPUT / _FINAL_SAMPLES of processing.smk, which the
            harness does not execute (the file is not plain Python), so they are evaluated here.
  * 787ffe  mode=processing asks for sample_qc_by_site.tsv, and plot_qc_processing makes it with
            the input metadata (every sample's site, dropped ones included) among its inputs.
"""
import os
import re
import unittest

import _smk_harness as H
from test_emmax_imputed_input import rule_graph

PROCESSING_SMK = os.path.join(H.REPO_ROOT, "workflow", "rules", "processing.smk")
SHARED = ("_FILTER_VCF_INPUT", "_FINAL_SAMPLES")


def load(overrides=None, has_dp=False):
    """common.smk's namespace plus processing.smk's two shared filter inputs, and its rules."""
    ns = H.load(H.make_config(**(overrides or {})))
    # common.smk sets HAS_FORMAT_DP from a bcftools read of the raw VCF's header, which the
    # harness does not have; the depth-filter branch needs it True.
    ns["HAS_FORMAT_DP"] = has_dp
    with open(PROCESSING_SMK) as handle:
        text = handle.read()
    for name in SHARED:
        line = re.search(rf"^{name} = .+$", text, re.M)
        assert line, f"processing.smk no longer defines {name}"
        exec(line.group(0), ns)
    return ns, {r["name"]: r for r in rule_graph(ns)}


# The harness's base config removes related samples (relatedness 0.99, action remove), as
# SIMDATA's does; "keep" isolates each other filter.
KEEP = {"Filter__relatedness_action": "keep"}
CONFIGS = {
    "plain": (KEEP, False),
    "depth filter": ({**KEEP, "Filter__min_depth": 5, "Filter__max_depth": 60}, True),
    "het filter": ({**KEEP, "Filter__het_outlier_sd": 2}, False),
    "relatedness removal": ({"Filter__relatedness": 0.9,
                             "Filter__relatedness_action": "remove"}, False),
}


class TestFilterStepsReplayFilterVcf(unittest.TestCase):

    def test_same_vcf_and_keep_list_in_every_configuration(self):
        for label, (overrides, has_dp) in CONFIGS.items():
            with self.subTest(config=label):
                ns, rules = load(overrides, has_dp)
                replay, real = rules["filter_vcf_steps"], rules["filter_vcf"]
                for key in ("vcf", "samples"):
                    self.assertIsNotNone(real["input"][key])
                    self.assertEqual(replay["input"][key], real["input"][key])
                self.assertEqual(replay["params"]["maf"], ns["MAF"])
                self.assertEqual(replay["params"]["miss"], ns["MISS"])

    def test_the_shared_inputs_follow_the_filters_that_ran(self):
        expected = {
            "plain": ("raw", "samples_filtered"),
            "depth filter": ("vcf_depth_filtered", "samples_filtered"),
            "het filter": ("raw", "samples_het_filtered"),
            "relatedness removal": ("raw", "samples_rel_filtered"),
        }
        for label, (overrides, has_dp) in CONFIGS.items():
            with self.subTest(config=label):
                ns, rules = load(overrides, has_dp)
                vcf_key, samples_key = expected[label]
                vcf = (f"{ns['INDIR']}{ns['VCF_RAW']}" if vcf_key == "raw" else ns["W"][vcf_key])
                self.assertEqual(rules["filter_vcf"]["input"]["vcf"], vcf)
                self.assertEqual(rules["filter_vcf"]["input"]["samples"], ns["W"][samples_key])

    def test_plot_qc_processing_reads_the_replay_and_the_filters_lists(self):
        ns, rules = load({"Filter__het_outlier_sd": 2, "Filter__relatedness": 0.9,
                          "Filter__relatedness_action": "remove"})
        W, O = ns["W"], ns["O"]
        self.assertEqual(rules["filter_vcf_steps"]["output"]["_0"], O["qc_filter_steps"])
        plot = rules["plot_qc_processing"]
        self.assertEqual(plot["input"]["filter_steps"], O["qc_filter_steps"])
        self.assertEqual(plot["input"]["samples_filtered"], W["samples_filtered"])
        self.assertEqual(plot["params"]["samples_het"], W["samples_het_filtered"])
        self.assertEqual(plot["params"]["samples_rel"], W["samples_rel_filtered"])
        ns, rules = load(KEEP)
        plot = rules["plot_qc_processing"]
        self.assertEqual((plot["params"]["samples_het"], plot["params"]["samples_rel"]),
                         ("NULL", "NULL"))


class TestDepthFromTheRawVcf(unittest.TestCase):

    def test_compute_depth_stats_reads_the_input_vcf_over_the_keep_list(self):
        ns, rules = load(has_dp=True)
        depth = rules["compute_depth_stats"]
        self.assertEqual(depth["input"]["vcf"], f"{ns['INDIR']}{ns['VCF_RAW']}")
        self.assertNotEqual(depth["input"]["vcf"], ns["W"]["vcf_filt"])
        self.assertEqual(depth["input"]["samples"], ns["W"]["samples_list"])

    def test_the_keep_list_is_cut_to_iids_and_gz_input_is_read(self):
        shell = load(has_dp=True)[1]["compute_depth_stats"]["shell"]
        self.assertRegex(shell, r"awk '\{\{print \$2\}\}' \{input\.samples\} > \{params\.prefix\}\.keep")
        self.assertEqual(len(re.findall(r"--keep \{params\.prefix\}\.keep", shell)), 2)
        self.assertIn("*.gz) VCF_FLAG=--gzvcf", shell)


class TestSiteTableIsTargeted(unittest.TestCase):

    def test_processing_asks_for_it_and_plot_qc_processing_makes_it(self):
        ns, rules = load()
        path = ns["O"]["qc_site_summary"]
        self.assertTrue(path.endswith("/Processing/tables/sample_qc_by_site.tsv"))
        self.assertIn(path, H.targets(ns, "processing"))
        patterns = H.producible_patterns(ns)
        makers = sorted({rule for regex, _smk, rule in patterns if regex.match(path)})
        self.assertEqual(makers, ["plot_qc_processing"])
        plot = rules["plot_qc_processing"]
        self.assertEqual(plot["output"]["site_qc"], path)
        self.assertEqual(plot["input"]["metadata_raw"], f"{ns['INDIR']}{ns['SAMPLES']}")


if __name__ == "__main__":
    unittest.main()
