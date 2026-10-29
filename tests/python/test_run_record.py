"""Run records (665953, 96ccfb) and the method-parameter tables (ce8cd4, release half).

``workflow/run_record.py`` builds and writes ``_intermediate/run_records/<mode>.json`` and
``history.jsonl``. ``common.smk``'s write_run_record() feeds it from Snakemake's state and is
called by the Snakefile's onstart/onsuccess/onerror handlers. method_param_rows() feeds both
the record's methods block and ``<method>_params_K<k>.tsv`` beside each p-value table.

What is pinned here, beyond "it writes a file":
  * a record never carries a git author/committer, an e-mail or any environment value but
    CLINEGO_IMAGE_ID: results trees get shared;
  * a record is never written at parse time, and an error inside the hook never stops a run;
  * the params table is deterministic text and its name escapes every glob that discovers
    p-value, sig-SNP or WZA tables (the app's fct_paths.R / fct_data_loading.R and
    check_invariants.R), so no reader takes it for one of those tables.
"""
import fnmatch
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from unittest import mock

import _smk_harness as smk
from _support import REPO_ROOT

sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))
import run_record  # noqa: E402

AUTHOR = "SENTINEL_AUTHOR_NAME"
EMAIL = "sentinel.author@example.invalid"


def _git(repo, *args):
    return subprocess.run(["git", "-C", repo, "-c", f"user.name={AUTHOR}", "-c",
                           f"user.email={EMAIL}", *args],
                          check=True, capture_output=True).stdout


def _repo(test):
    """A one-commit checkout in a fresh temp dir."""
    d = tempfile.mkdtemp(prefix="clinego_rr_repo_")
    test.addCleanup(shutil.rmtree, d, ignore_errors=True)
    _git(d, "init", "-q")
    with open(os.path.join(d, "tracked.txt"), "w") as handle:
        handle.write("one\n")
    _git(d, "add", "tracked.txt")
    _git(d, "commit", "-q", "-m", "first")
    return d


def _tmpdir(test):
    d = tempfile.mkdtemp(prefix="clinego_rr_")
    test.addCleanup(shutil.rmtree, d, ignore_errors=True)
    return d


class GitStateTest(unittest.TestCase):
    def test_clean_checkout_gives_its_commit_and_not_dirty(self):
        repo = _repo(self)
        state = run_record.git_state(repo)
        self.assertEqual(state["commit"], _git(repo, "rev-parse", "HEAD").decode().strip())
        self.assertIs(state["dirty"], False)
        self.assertNotIn("changed_paths", state)
        self.assertNotIn("diff_sha256", state)

    def test_an_edited_tracked_file_is_dirty_with_the_hash_of_its_diff(self):
        repo = _repo(self)
        with open(os.path.join(repo, "tracked.txt"), "a") as handle:
            handle.write("two\n")
        state = run_record.git_state(repo)
        self.assertIs(state["dirty"], True)
        self.assertEqual(state["changed_paths"], ["tracked.txt"])
        self.assertEqual(state["n_changed_paths"], 1)
        diff = _git(repo, "diff", "HEAD", "--no-color", "--no-ext-diff")
        self.assertEqual(state["diff_sha256"], hashlib.sha256(diff).hexdigest())

    def test_an_untracked_file_makes_the_checkout_dirty(self):
        repo = _repo(self)
        with open(os.path.join(repo, "new_script.R"), "w") as handle:
            handle.write("x <- 1\n")
        state = run_record.git_state(repo)
        self.assertIs(state["dirty"], True)
        self.assertIn("new_script.R", state["changed_paths"])

    def test_reading_the_state_never_rewrites_the_checkouts_index(self):
        # A tracked file with a newer mtime and the same content: plain `git status` refreshes
        # that entry and rewrites .git/index under index.lock. The record runs at every mode
        # start against the user's live checkout, so it must not.
        repo = _repo(self)
        tracked = os.path.join(repo, "tracked.txt")
        st = os.stat(tracked)
        os.utime(tracked, (st.st_atime + 10, st.st_mtime + 10))
        index = os.path.join(repo, ".git", "index")
        before = os.stat(index).st_mtime_ns
        self.assertIs(run_record.git_state(repo)["dirty"], False)
        self.assertEqual(os.stat(index).st_mtime_ns, before)
        self.assertFalse(os.path.exists(index + ".lock"))
        # Control: the same state does make plain `git status` write the index.
        _git(repo, "status", "--porcelain")
        self.assertNotEqual(os.stat(index).st_mtime_ns, before)

    def test_a_directory_that_is_no_checkout_gives_no_commit_and_a_reason(self):
        state = run_record.git_state(_tmpdir(self))
        self.assertIsNone(state["commit"])
        self.assertIsNone(state["dirty"])
        self.assertTrue(state["note"].startswith("git: "), state["note"])

    def test_a_missing_git_gives_no_commit_and_a_reason(self):
        d = _tmpdir(self)
        with mock.patch.object(run_record.subprocess, "run", side_effect=FileNotFoundError):
            state = run_record.git_state(d)
        self.assertEqual(state, {"commit": None, "dirty": None, "note": "git not found"})


class ImageStateTest(unittest.TestCase):
    def test_an_image_id_passed_in_is_recorded_as_visible(self):
        state = run_record.image_state("sha256:abc", markers=())
        self.assertEqual(state, {"in_container": False, "id": "sha256:abc", "visible": True})

    def test_without_the_variable_the_image_is_not_visible_and_says_how_to_pass_it(self):
        marker = os.path.join(_tmpdir(self), ".dockerenv")
        open(marker, "w").close()
        state = run_record.image_state("  ", markers=(marker,))
        self.assertIs(state["in_container"], True)
        self.assertIsNone(state["id"])
        self.assertIs(state["visible"], False)
        self.assertIn(run_record.IMAGE_ENV, state["note"])


def _record(test, status="success", repo=None, configfiles=None, **kw):
    repo = repo or _repo(test)
    if configfiles is None:
        cfg = os.path.join(_tmpdir(test), "config_X.yaml")
        with open(cfg, "w") as handle:
            handle.write("project_name: X\ncpu: 4\n")
        configfiles = [cfg]
    args = dict(mode="gea", status=status, started="2026-10-29T10:00:00+00:00",
                finished=None if status == "running" else "2026-10-29T10:05:00+00:00",
                argv=["/usr/local/bin/snakemake", "-c4", "--config", "mode=gea",
                      "--configfile", "config_X.yaml"],
                launch_dir="/pipeline", configfiles=configfiles,
                config={"project_name": "X", "cpu": 4, "mode": "gea"}, repo=repo,
                image_id=None, snakemake_version="9.27.1", cpu=4,
                methods={"GEA": {"LFMM": {"K": {"value": "3", "set_by": "sNMF.k_best"}}}},
                snakemake_log=["/pipeline/.snakemake/log/x.snakemake.log"])
    args.update(kw)
    return run_record.build_record(**args)


class RecordTest(unittest.TestCase):
    def test_fields_come_in_the_written_order_with_the_config_last(self):
        record = _record(self)
        self.assertEqual(list(record), [
            "format", "mode", "status", "started", "finished", "command", "launch_dir", "configfiles",
            "git", "image", "software", "cpu", "snakemake_log", "methods", "config"])
        self.assertEqual(record["command"],
                         "snakemake -c4 --config mode=gea --configfile config_X.yaml")
        self.assertEqual(record["software"]["snakemake"], "9.27.1")

    def test_each_configfile_carries_the_sha256_of_its_bytes(self):
        cfg = os.path.join(_tmpdir(self), "config_Y.yaml")
        with open(cfg, "wb") as handle:
            handle.write(b"project_name: Y\n")
        record = _record(self, configfiles=[cfg])
        self.assertEqual(record["configfiles"],
                         [{"path": cfg, "sha256": hashlib.sha256(b"project_name: Y\n").hexdigest()}])

    def test_an_unreadable_configfile_is_noted_not_raised(self):
        missing = os.path.join(_tmpdir(self), "gone.yaml")
        entry = _record(self, configfiles=[missing])["configfiles"][0]
        self.assertIsNone(entry["sha256"])
        self.assertIn("FileNotFoundError", entry["note"])

    def test_no_author_email_or_environment_value_reaches_the_record(self):
        repo = _repo(self)
        with open(os.path.join(repo, "tracked.txt"), "a") as handle:
            handle.write("dirty\n")
        with mock.patch.dict(os.environ, {"CLINEGO_TEST_SECRET": "sentinel-secret-value"}):
            text = json.dumps(_record(self, repo=repo))
        for leak in (AUTHOR, EMAIL, "sentinel-secret-value", "CLINEGO_TEST_SECRET"):
            self.assertNotIn(leak, text)


class WriteTest(unittest.TestCase):
    def test_running_writes_the_record_and_no_history(self):
        inter = _tmpdir(self)
        path = run_record.write(inter, _record(self, status="running"))
        self.assertEqual(path, os.path.join(inter, "run_records", "gea.json"))
        self.assertEqual(os.listdir(os.path.join(inter, "run_records")), ["gea.json"])
        with open(path) as handle:
            loaded = json.load(handle)
        self.assertEqual(loaded["format"], run_record.FORMAT)
        self.assertEqual(loaded["status"], "running")
        self.assertIsNone(loaded["finished"])
        self.assertEqual(loaded["methods"]["GEA"]["LFMM"]["K"]["value"], "3")

    def test_each_finished_run_appends_one_history_line(self):
        inter = _tmpdir(self)
        run_record.write(inter, _record(self, status="running"))
        run_record.write(inter, _record(self, status="success"))
        run_record.write(inter, _record(self, status="error"))
        directory = os.path.join(inter, "run_records")
        self.assertEqual(sorted(os.listdir(directory)), ["gea.json", "history.jsonl"])
        with open(os.path.join(directory, "history.jsonl")) as handle:
            lines = handle.read().splitlines()
        self.assertEqual(len(lines), 2)
        history = [json.loads(line) for line in lines]
        self.assertEqual([h["status"] for h in history], ["success", "error"])
        self.assertEqual(history[0]["mode"], "gea")
        self.assertEqual(len(history[0]["commit"]), 40)
        with open(os.path.join(directory, "gea.json")) as handle:
            self.assertEqual(json.load(handle)["status"], "error")


def _load(**overrides):
    cfg = smk.make_config(**overrides)
    return smk.load(cfg)


LFMM_K2 = [{"method": "EMMAX", "adjust": "bonf", "threshold": "0.05"},
           {"method": "LFMM", "adjust": "bonf", "threshold": "0.05", "params": {"K": 2}},
           {"method": "RDA", "adjust": "bonf", "threshold": "0.01",
            "params": {"condition_pcs": 3, "axes": "auto"}}]


class MethodParamsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ns = _load(GEA__configs=LFMM_K2)

    def rows(self, source, method):
        return {p: (v, s) for p, v, s in self.ns["method_param_rows"](source, method)}

    def test_a_configured_k_is_config_and_the_file_name_k_is_k_best(self):
        rows = self.rows("GEA", "LFMM")
        self.assertEqual(rows["engine"], ("lfmm", "registry"))
        self.assertEqual(rows["adjust"], ("bonf", "config"))
        self.assertEqual(rows["threshold"], ("0.05", "config"))
        self.assertEqual(rows["K"], ("2", "config"))
        self.assertEqual(rows["genomic_control"], ("true", "default"))
        self.assertEqual(rows["file_name_k"], ("3", "sNMF.k_best"))

    def test_an_unset_structure_parameter_is_sNMF_k_best(self):
        rows = self.rows("GEA", "EMMAX")
        self.assertEqual(rows["n_pcs"], ("3", "sNMF.k_best"))
        self.assertEqual(rows["kinship"], ("BN", "default"))

    def test_rda_auto_predictor_set_is_resolved_to_the_climate_predictors(self):
        rows = self.rows("GEA", "RDA")
        self.assertEqual(rows["predictor_set"], ("bio_1,bio_2,bio_3", "Climate.predictors"))
        self.assertEqual(rows["condition_pcs"], ("3", "config"))
        self.assertEqual(rows["axes"], ("auto", "config"))
        self.assertEqual(rows["seed"], ("42", "default"))

    def test_gwas_rows_carry_the_gapit_model_and_the_missing_strategy(self):
        rows = self.rows("GWAS", "BLINK")
        self.assertEqual(rows["engine"], ("gapit", "registry"))
        self.assertEqual(rows["gapit_model"], ("BLINK", "registry"))
        self.assertEqual(rows["n_pcs"], ("3", "sNMF.k_best"))
        self.assertEqual(rows["missing_strategy"], ("DROP", "config"))

    def test_the_table_text_is_deterministic_and_four_columns_wide(self):
        text = self.ns["method_params_tsv"]("GEA", "LFMM")
        self.assertEqual(text, self.ns["method_params_tsv"]("GEA", "LFMM"))
        lines = text.splitlines()
        self.assertEqual(lines[0], "method\tparameter\tvalue\tset_by")
        self.assertTrue(all(len(line.split("\t")) == 4 for line in lines))
        self.assertTrue(all(line.startswith("LFMM\t") for line in lines[1:]))
        self.assertIsNone(re.search(r"20\d\d-\d\d-\d\d", text))

    def test_every_configured_method_gets_its_params_table_as_a_target(self):
        for mode, source in (("gea", "GEA"), ("gwas", "GWAS")):
            targets = smk.targets(self.ns, mode)
            for method in self.ns["ASSOC_SOURCES"][source]["configs"]:
                path = self.ns["method_params_table"](source, method)
                self.assertIn(path, targets)
                self.assertTrue(path.endswith(f"/{source}/tables/methods/{method}/"
                                              f"{method}_params_K3.tsv"), path)

    def test_the_params_table_escapes_every_table_discovery_glob(self):
        # fct_paths.R:401/411 and fct_data_loading.R:834 glob *_pvalues_K*.tsv for the
        # p-value tables, fct_paths.R:191 and check_invariants.R:164 *_sig_snps_*.tsv,
        # fct_data_loading.R:985 *_wza_K*.tsv.
        for source in ("GEA", "GWAS"):
            for method in self.ns["ASSOC_SOURCES"][source]["configs"]:
                name = os.path.basename(self.ns["method_params_table"](source, method))
                for glob in ("*_pvalues_K*.tsv", "*_sig_snps_*.tsv", "*_wza_K*.tsv"):
                    self.assertFalse(fnmatch.fnmatch(name, glob), (name, glob))


class WriteRunRecordTest(unittest.TestCase):
    def test_loading_common_smk_writes_no_run_record(self):
        _load(GEA__configs=LFMM_K2)
        for _root, dirs, _files in os.walk(smk._SANDBOX):
            self.assertNotIn("run_records", dirs)

    def test_a_gea_record_carries_the_gea_method_parameters(self):
        ns = _load(GEA__configs=LFMM_K2)
        inter = _tmpdir(self) + "/"
        ns["INTER"], ns["MODE"] = inter, "gea"
        ns["write_run_record"]("running", ["/pipeline/.snakemake/log/x.snakemake.log"])
        with open(os.path.join(inter, "run_records", "gea.json")) as handle:
            record = json.load(handle)
        self.assertEqual(record["status"], "running")
        self.assertEqual(record["git"]["commit"],
                         _git(REPO_ROOT, "-c", "safe.directory=*", "rev-parse", "HEAD").decode().strip())
        self.assertEqual(record["methods"]["GEA"]["LFMM"]["K"], {"value": "2", "set_by": "config"})
        self.assertEqual(sorted(record["methods"]), ["GEA"])
        self.assertEqual(record["snakemake_log"], ["/pipeline/.snakemake/log/x.snakemake.log"])
        self.assertEqual(record["config"]["GEA"]["configs"][1]["params"], {"K": 2})

    def test_a_mode_without_association_methods_has_no_methods_block(self):
        ns = _load()
        inter = _tmpdir(self) + "/"
        ns["INTER"], ns["MODE"] = inter, "structure"
        ns["write_run_record"]("success", "/pipeline/.snakemake/log/one.snakemake.log")
        with open(os.path.join(inter, "run_records", "structure.json")) as handle:
            record = json.load(handle)
        self.assertNotIn("methods", record)
        self.assertEqual(record["snakemake_log"], ["/pipeline/.snakemake/log/one.snakemake.log"])

    def test_a_failing_write_warns_and_never_raises(self):
        ns = _load()
        blocker = os.path.join(_tmpdir(self), "a_file")
        open(blocker, "w").close()
        ns["INTER"], ns["MODE"] = blocker + "/", "gea"
        err = io.StringIO()
        with redirect_stderr(err):
            ns["write_run_record"]("success", [])
        self.assertIn("WARNING: run record not written", err.getvalue())


class SnakefileHooksTest(unittest.TestCase):
    def test_the_three_handlers_write_the_record_with_their_status(self):
        with open(os.path.join(REPO_ROOT, "Snakefile")) as handle:
            text = handle.read()
        for hook, status in (("onstart", "running"), ("onsuccess", "success"),
                             ("onerror", "error")):
            self.assertRegex(text, rf"(?m)^{hook}:\n    write_run_record\(\"{status}\", log\)$")


if __name__ == "__main__":
    unittest.main()
