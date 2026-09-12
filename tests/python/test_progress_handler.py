"""scripts/snakemake_progress_handler.py — the log handler the Shiny app polls.

Snakemake calls main(msg) per log event; the records land in
{PROJECT}_logs/.progress.jsonl. Two things make this file need care:

  * MODULE-LEVEL MUTABLE STATE at :56-58 (_progress_file, _total_jobs,
    _truncated). Every test resets it in setUp, or the suite becomes
    order-dependent and the second test writes into the first one's tempdir.
  * get_progress_file() reads os.getcwd() AND sys.argv, and calls os.makedirs().
    It is driven here with a real tempdir and a patched argv rather than mocked,
    because the directory creation is part of the contract.
"""
import json
import os
import sys
import tempfile
import unittest

import _support

sph = _support.load_script("snakemake_progress_handler")


class ProgressCase(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.d = self.dir.name

        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        os.chdir(self.d)

        argv = sys.argv
        self.addCleanup(setattr, sys, "argv", argv)

        # The module globals survive an import, so reset them per test.
        sph._progress_file = None
        sph._total_jobs = None
        sph._truncated = False
        self.addCleanup(self.reset_globals)

    def reset_globals(self):
        sph._progress_file = None
        sph._total_jobs = None
        sph._truncated = False

    def records(self):
        with open(sph._progress_file) as fh:
            return [json.loads(ln) for ln in fh if ln.strip()]


class TestGetProgressFile(ProgressCase):

    def test_project_is_taken_from_the_configfile_argument(self):
        sys.argv = ["snakemake", "--configfile", "config_SIMDATA.yaml"]
        got = sph.get_progress_file()
        self.assertEqual(got, os.path.join(self.d, "SIMDATA_logs", ".progress.jsonl"))
        self.assertTrue(os.path.isdir(os.path.join(self.d, "SIMDATA_logs")))

    def test_a_configfile_path_is_reduced_to_its_basename(self):
        sys.argv = ["snakemake", "--configfile", "/somewhere/else/config_WBDC.yaml"]
        self.assertEqual(sph.get_progress_file(),
                         os.path.join(self.d, "WBDC_logs", ".progress.jsonl"))

    def test_bare_config_yaml_falls_through_to_the_scan(self):
        os.makedirs(os.path.join(self.d, "ONLYONE_logs"))
        sys.argv = ["snakemake", "--configfile", "config.yaml"]
        self.assertEqual(sph.get_progress_file(),
                         os.path.join(self.d, "ONLYONE_logs", ".progress.jsonl"))

    def test_no_configfile_and_no_logs_dir_falls_back_to_tmp(self):
        sys.argv = ["snakemake", "-c4"]
        self.assertEqual(sph.get_progress_file(), os.path.join("/tmp", ".progress.jsonl"))

    def test_a_plain_file_named_like_a_logs_dir_is_not_taken(self):
        open(os.path.join(self.d, "NOTADIR_logs"), "w").close()
        sys.argv = ["snakemake"]
        self.assertEqual(sph.get_progress_file(), os.path.join("/tmp", ".progress.jsonl"))

    def test_multiple_logs_dirs_resolve_unsorted_and_arbitrarily(self):
        """PINNED DEFECT, filed in docs/pipeline_improvement_requests.md.

        snakemake_progress_handler.py:48-51 takes candidates[0] from an unsorted
        os.listdir(). With several projects in one working directory the progress
        file lands in an arbitrary project's log dir — and which one can change
        between runs on the same tree, so the Shiny app can poll a file the
        running job never writes to.
        """
        for name in ("ALPHA_logs", "BETA_logs", "GAMMA_logs"):
            os.makedirs(os.path.join(self.d, name))
        sys.argv = ["snakemake"]
        got = sph.get_progress_file()
        listed = [x for x in os.listdir(self.d) if x.endswith("_logs")]
        self.assertEqual(got, os.path.join(self.d, listed[0], ".progress.jsonl"))
        # The contract this test pins is "listdir order", NOT "sorted order".
        # Assert only that the resolution follows listdir, which is the defect.
        self.assertIn(os.path.basename(os.path.dirname(got)),
                      {"ALPHA_logs", "BETA_logs", "GAMMA_logs"})


class TestWrite(ProgressCase):

    def test_truncates_once_then_appends(self):
        sys.argv = ["snakemake", "--configfile", "config_P.yaml"]
        path = sph.get_progress_file()
        with open(path, "w") as fh:
            fh.write('{"stale": true}\n')

        sph._write({"level": "a"})
        sph._write({"level": "b"})
        recs = self.records()
        self.assertEqual([r["level"] for r in recs], ["a", "b"])   # stale line gone

    def test_every_record_is_timestamped(self):
        sys.argv = ["snakemake", "--configfile", "config_P.yaml"]
        sph._write({"level": "a"})
        rec = self.records()[0]
        self.assertTrue(rec["timestamp"].endswith("Z"))


class TestMain(ProgressCase):
    """The four record shapes documented in the module docstring at :13-17."""

    def setUp(self):
        super().setUp()
        sys.argv = ["snakemake", "--configfile", "config_P.yaml"]

    def test_job_info(self):
        sph.main({"level": "job_info", "name": "filter_vcf", "jobid": 7,
                  "log": ["/logs/filter_vcf.log"]})
        rec = self.records()[0]
        self.assertEqual(rec["level"], "job_info")
        self.assertEqual(rec["rule"], "filter_vcf")
        self.assertEqual(rec["jobid"], 7)
        self.assertEqual(rec["log_file"], "/logs/filter_vcf.log")
        self.assertNotIn("total_jobs", rec)          # nothing has supplied it yet

    def test_job_info_falls_back_from_name_to_rule(self):
        sph.main({"level": "job_info", "rule": "ld_prune", "jobid": 1})
        self.assertEqual(self.records()[0]["rule"], "ld_prune")

    def test_job_info_with_no_log_records_none(self):
        sph.main({"level": "job_info", "name": "r", "jobid": 1, "log": []})
        self.assertIsNone(self.records()[0]["log_file"])

    def test_job_info_scrapes_the_total_out_of_the_message_text(self):
        sph.main({"level": "job_info", "name": "r", "jobid": 3,
                  "msg": "rule r:\n    3 of 42 steps (7%) done"})
        self.assertEqual(self.records()[0]["total_jobs"], 42)

    def test_the_scraped_total_persists_to_later_records(self):
        sph.main({"level": "job_info", "name": "a", "jobid": 1, "msg": "1 of 9 steps"})
        sph.main({"level": "job_info", "name": "b", "jobid": 2})
        self.assertEqual(self.records()[1]["total_jobs"], 9)

    def test_job_finished_and_job_error(self):
        sph.main({"level": "job_finished", "name": "snmf", "jobid": 2, "log": ["/l.log"]})
        sph.main({"level": "job_error", "name": "emmax", "jobid": 3})
        recs = self.records()
        self.assertEqual(recs[0], {"level": "job_finished", "rule": "snmf", "jobid": 2,
                                   "log_file": "/l.log",
                                   "timestamp": recs[0]["timestamp"]})
        self.assertEqual(recs[1]["level"], "job_error")
        self.assertEqual(recs[1]["rule"], "emmax")

    def test_dag_info(self):
        sph.main({"level": "dag_info", "total_jobs": 31})
        self.assertEqual(self.records()[0], {"level": "dag_info", "total_jobs": 31,
                                             "timestamp": self.records()[0]["timestamp"]})

    def test_dag_info_without_a_total_writes_nothing(self):
        sph.main({"level": "dag_info"})
        self.assertIsNone(sph._progress_file)

    def test_progress_and_run_info_are_rewritten_as_dag_info(self):
        sph.main({"level": "progress", "total": 12})
        sph.main({"level": "run_info", "total": 12})
        self.assertEqual([r["level"] for r in self.records()], ["dag_info", "dag_info"])

    def test_an_unknown_level_writes_nothing(self):
        sph.main({"level": "d3bug", "msg": "hello"})
        self.assertIsNone(sph._progress_file)

    def test_a_message_with_no_level_writes_nothing(self):
        sph.main({"msg": "hello"})
        self.assertIsNone(sph._progress_file)


if __name__ == "__main__":
    unittest.main()
