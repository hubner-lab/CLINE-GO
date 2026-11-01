"""The clinego launcher (repo root), driven as a subprocess with no docker.

The gate runs this suite inside an image that has no `docker`, and from the release on that
image is the baked one, which sets CLINEGO_IN_IMAGE=1 for every process. So every call here
builds its environment from scratch: PATH starts with stub `docker`, `git` and `snakemake`
scripts, CLINEGO_IN_IMAGE / CLINEGO_IMAGE / CLINEGO_IMAGE_ID are set or removed on purpose,
and the stub git pins the version. The tests pass the same way on the host, in
cline-go:latest and in the baked image.

The stub docker answers `image inspect` with a fixed ID (or fails when STUB_NO_IMAGE=1) and
`--version`; --show output is compared as an argv (shlex.split), never as text, because
bash's printf %q spells some arguments differently from a Python literal.
"""
import os
import shlex
import shutil
import stat
import subprocess
import tempfile
import unittest

from _support import REPO_ROOT

LAUNCHER = os.path.join(REPO_ROOT, "clinego")
MODES = ["processing", "prestructure", "structure", "climate", "traits", "pregea", "gea",
         "gwas", "gea_x_gwas", "maladaptation"]
IMAGE_ID = "sha256:0123456789abcdef"
VERSION = "1.0-3-gabc1234"
COMMIT = "abc1234def5678abc1234def5678abc1234def56"

STUB_DOCKER = """#!/bin/sh
if [ "$1" = image ] && [ "$2" = inspect ]; then
    [ "${STUB_NO_IMAGE:-}" = 1 ] && { echo "Error: No such image" >&2; exit 1; }
    echo "%s"
    exit 0
fi
printf '%%s\\n' "docker $*"
""" % IMAGE_ID

# describe / rev-parse from STUB_DESCRIBE (default a tagged, clean checkout); exits 128
# like git outside a checkout when STUB_NO_GIT=1.
STUB_GIT = """#!/bin/sh
[ "${STUB_NO_GIT:-}" = 1 ] && { echo "fatal: not a git repository" >&2; exit 128; }
for a in "$@"; do
    case $a in
        describe) echo "${STUB_DESCRIBE:-v%s}"; exit 0 ;;
        rev-parse) echo "%s"; exit 0 ;;
    esac
done
exit 1
""" % (VERSION, COMMIT)

STUB_SNAKEMAKE = """#!/bin/sh
for a in "$@"; do printf '%s\\n' "$a"; done
"""


def _write_exe(path, text):
    with open(path, "w") as handle:
        handle.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


class LauncherCase(unittest.TestCase):
    """A temp project folder (cwd), a stub bin dir, and the launcher reached by symlink."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="clinego_cli_")
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.bin = os.path.join(self.tmp, "bin")
        self.project = os.path.join(os.path.realpath(self.tmp), "project")
        os.makedirs(self.bin)
        os.makedirs(self.project)
        _write_exe(os.path.join(self.bin, "docker"), STUB_DOCKER)
        _write_exe(os.path.join(self.bin, "git"), STUB_GIT)
        _write_exe(os.path.join(self.bin, "snakemake"), STUB_SNAKEMAKE)
        # Through a symlink, as ~/.local/bin/clinego and /usr/local/bin/clinego reach it:
        # the launcher must still find its repo.
        self.link = os.path.join(self.bin, "clinego")
        os.symlink(LAUNCHER, self.link)
        with open(os.path.join(self.project, "config_demo.yaml"), "w") as handle:
            handle.write("project_name: demo\n")
        self.user = "%d:%d" % (os.getuid(), os.getgid())

    def env(self, in_image=False, **extra):
        env = {"PATH": self.bin + ":/usr/bin:/bin", "HOME": self.tmp, "LC_ALL": "C"}
        if in_image:
            env["CLINEGO_IN_IMAGE"] = "1"
        env.update(extra)
        return env

    def call(self, *args, in_image=False, launcher=None, cwd=None, **extra):
        return subprocess.run(["bash", launcher or self.link] + list(args),
                              cwd=cwd or self.project, env=self.env(in_image, **extra),
                              capture_output=True, text=True, timeout=30)

    def show(self, *args, **kw):
        result = self.call(*args, **kw)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        return shlex.split(result.stdout)

    def refuse(self, *args, contains=None, **kw):
        result = self.call(*args, **kw)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        lines = result.stderr.splitlines()
        self.assertEqual(len(lines), 1, result.stderr)
        self.assertTrue(lines[0].startswith("clinego: "), lines[0])
        if contains:
            self.assertIn(contains, lines[0])
        return lines[0]


class HelpAndVersionTest(LauncherCase):
    def test_help_names_every_subcommand_and_mode(self):
        for args in ([], ["help"]):
            result = self.call(*args)
            self.assertEqual(result.returncode, 0, result.stderr)
            for sub in ["help", "version", "build", "example", "init", "run", "unlock", "app"]:
                self.assertRegex(result.stdout, r"(?m)^  %s\b" % sub)
            for mode in MODES:
                self.assertIn(mode, result.stdout)

    def test_version_strips_the_v_and_reports_the_image(self):
        result = self.call("version")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("clinego " + VERSION + " ", result.stdout)
        self.assertIn("cline-go:%s: built (%s)" % (VERSION, IMAGE_ID), result.stdout)
        result = self.call("version", STUB_NO_IMAGE="1")
        self.assertIn("not built (run: clinego build)", result.stdout)

    def test_version_works_without_docker(self):
        os.remove(os.path.join(self.bin, "docker"))
        result = self.call("version", PATH=self._path_without_docker())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("clinego " + VERSION, result.stdout)
        self.assertIn("unknown (docker not found)", result.stdout)

    def test_version_falls_back_to_the_version_file(self):
        repo = os.path.join(self.tmp, "clone")
        os.makedirs(repo)
        shutil.copy(LAUNCHER, os.path.join(repo, "clinego"))
        with open(os.path.join(repo, "VERSION"), "w") as handle:
            handle.write("1.0\n")
        result = self.call("version", launcher=os.path.join(repo, "clinego"), STUB_NO_GIT="1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("clinego 1.0 ", result.stdout)

    def test_version_in_image_reads_build_info(self):
        repo = os.path.join(self.tmp, "pipeline")
        os.makedirs(repo)
        shutil.copy(LAUNCHER, os.path.join(repo, "clinego"))
        with open(os.path.join(repo, "BUILD_INFO.json"), "w") as handle:
            handle.write('{"version": "1.0", "commit": "%s", "dirty": false, '
                         '"built": "2026-11-01T00:00:00Z"}' % COMMIT)
        if shutil.which("python3", path="/usr/bin:/bin") is None:
            self.skipTest("no python3 to read BUILD_INFO.json")
        result = self.call("version", launcher=os.path.join(repo, "clinego"),
                           in_image=True, STUB_NO_GIT="1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("clinego 1.0\n", result.stdout)
        self.assertIn(COMMIT + " (dirty: false)", result.stdout)
        self.assertIn("2026-11-01T00:00:00Z", result.stdout)

    def _path_without_docker(self):
        """A PATH holding the tools the launcher needs, and no docker anywhere."""
        tools = os.path.join(self.tmp, "tools")
        os.makedirs(tools)
        for name in ["bash", "sh", "dirname", "basename", "readlink", "head", "tr", "id",
                     "cat", "grep", "sed"]:
            found = shutil.which(name, path="/usr/bin:/bin")
            if found:
                os.symlink(found, os.path.join(tools, name))
        return self.bin + ":" + tools


class ShowTest(LauncherCase):
    def base_run(self, name):
        return ["docker", "run", "--rm", "--init", "--platform", "linux/amd64", "--name", name,
                "--user", self.user, "-v", self.project + ":/data", "-w", "/data",
                "-e", "CLINEGO_IMAGE_ID=" + IMAGE_ID]

    def test_run_defaults(self):
        argv = self.show("run", "processing", "--config", "config_demo.yaml", "--show")
        self.assertNotIn("--cpus", argv)  # docker refuses --cpus above a small VM's count
        self.assertEqual(argv, self.base_run("clinego-config_demo-processing") + [
            "cline-go:" + VERSION, "clinego", "run", "processing",
            "--config", "config_demo.yaml", "--cores", "4"])

    def test_run_with_every_option(self):
        argv = self.show("run", "gea", "--config=config_demo.yaml", "--cores", "16",
                         "--memory", "64g", "--image", "my/img:2", "--show",
                         "--", "-n", "--forcerun", "a rule")
        self.assertEqual(argv, self.base_run("clinego-config_demo-gea") + [
            "--memory", "64g", "--cpus", "16", "my/img:2", "clinego", "run", "gea",
            "--config", "config_demo.yaml", "--cores", "16", "--", "-n", "--forcerun", "a rule"])

    def test_image_order_is_flag_then_env_then_version(self):
        argv = self.show("run", "gwas", "--config", "config_demo.yaml", "--show",
                         CLINEGO_IMAGE="env/img:1")
        self.assertIn("env/img:1", argv)
        argv = self.show("run", "gwas", "--config", "config_demo.yaml", "--show",
                         "--image", "flag/img:1", CLINEGO_IMAGE="env/img:1")
        self.assertIn("flag/img:1", argv)
        self.assertNotIn("env/img:1", argv)

    def test_a_config_in_a_subfolder_is_accepted(self):
        os.makedirs(os.path.join(self.project, "configs"))
        with open(os.path.join(self.project, "configs", "config_b.yml"), "w") as handle:
            handle.write("project_name: b\n")
        argv = self.show("run", "traits", "--config", "configs/config_b.yml", "--show")
        self.assertIn("clinego-config_b-traits", argv)
        self.assertEqual(argv[-3:], ["configs/config_b.yml", "--cores", "4"])

    def test_unlock(self):
        argv = self.show("unlock", "--config", "config_demo.yaml", "--show")
        self.assertEqual(argv, self.base_run("clinego-config_demo-unlock") + [
            "cline-go:" + VERSION, "clinego", "unlock",
            "--config", "config_demo.yaml"])
        argv = self.show("unlock", "--config", "config_demo.yaml", "--image", "x/y:1", "--show")
        self.assertIn("x/y:1", argv)

    def test_app_defaults(self):
        argv = self.show("app", "--show")
        self.assertEqual(argv, [
            "docker", "run", "--rm", "--init", "--platform", "linux/amd64",
            "--name", "clinego-app-3838", "--user", self.user, "-p", "3838:3838",
            "-v", self.project + ":/data", "-w", "/data", "-e", "CLINEGO_PROJECT_ROOT=/data",
            "cline-go:" + VERSION, "clinego", "app"])

    def test_app_with_every_option(self):
        argv = self.show("app", "--port", "3871", "--cores", "16", "--memory", "64g",
                         "--image", "my/img:2", "--show")
        self.assertEqual(argv, [
            "docker", "run", "--rm", "--init", "--platform", "linux/amd64",
            "--name", "clinego-app-3871", "--user", self.user, "-p", "3871:3838",
            "-v", self.project + ":/data", "-w", "/data", "-e", "CLINEGO_PROJECT_ROOT=/data",
            "--memory", "64g", "--cpus", "16", "my/img:2", "clinego", "app"])

    def test_build_defaults(self):
        argv = self.show("build", "--show")
        self.assertEqual(argv, [
            "docker", "build", "--platform", "linux/amd64", "-t", "cline-go:" + VERSION,
            "--build-arg", "CLINEGO_VERSION=" + VERSION,
            "--build-arg", "CLINEGO_COMMIT=" + COMMIT,
            "--build-arg", "CLINEGO_DIRTY=false", REPO_ROOT])

    def test_build_with_every_option_and_a_dirty_tree(self):
        argv = self.show("build", "--tag", "cline-go:lane-packaging", "--no-cache", "--show",
                         STUB_DESCRIBE="a711825-dirty")
        self.assertEqual(argv, [
            "docker", "build", "--platform", "linux/amd64", "-t", "cline-go:lane-packaging",
            "--build-arg", "CLINEGO_VERSION=a711825-dirty",
            "--build-arg", "CLINEGO_COMMIT=" + COMMIT,
            "--build-arg", "CLINEGO_DIRTY=true", "--no-cache", REPO_ROOT])

    def test_build_without_git_passes_empty_commit_and_dirty(self):
        repo = os.path.join(self.tmp, "zip")
        os.makedirs(repo)
        shutil.copy(LAUNCHER, os.path.join(repo, "clinego"))
        with open(os.path.join(repo, "VERSION"), "w") as handle:
            handle.write("1.0\n")
        argv = self.show("build", "--show", launcher=os.path.join(repo, "clinego"),
                         STUB_NO_GIT="1")
        self.assertEqual(argv, [
            "docker", "build", "--platform", "linux/amd64", "-t", "cline-go:1.0",
            "--build-arg", "CLINEGO_VERSION=1.0", "--build-arg", "CLINEGO_COMMIT=",
            "--build-arg", "CLINEGO_DIRTY=", repo])

    def test_example_and_init_run_inside_the_image_on_this_folder(self):
        for args, name in ((["example"], "clinego-example"),
                           (["init", "demo2"], "clinego-init-demo2")):
            argv = self.show(*(args + ["--show"]))
            self.assertEqual(argv, [
                "docker", "run", "--rm", "--init", "--platform", "linux/amd64", "--name", name,
                "--user", self.user, "-v", self.project + ":/data", "-w", "/data",
                "cline-go:" + VERSION, "clinego"] + args)


class RefusalTest(LauncherCase):
    def test_unknown_mode_lists_the_ten(self):
        line = self.refuse("run", "gea_gwas", "--config", "config_demo.yaml")
        for mode in MODES:
            self.assertIn(mode, line)

    def test_config_outside_this_folder(self):
        self.refuse("run", "gea", "--config", "/etc/hostname", contains="relative path")
        self.refuse("run", "gea", "--config", "../project/config_demo.yaml",
                    contains="in this folder")
        self.refuse("run", "gea", "--config", "a/../config_demo.yaml", contains="in this folder")

    def test_missing_config(self):
        self.refuse("run", "gea", "--config", "config_nope.yaml", contains="not found")
        self.refuse("run", "gea", contains="--config")
        self.refuse("unlock", "--config", "config_nope.yaml", contains="not found")

    def test_missing_image(self):
        line = self.refuse("run", "gea", "--config", "config_demo.yaml", STUB_NO_IMAGE="1")
        self.assertEqual(line, "clinego: No image cline-go:%s for this code. Run: clinego build"
                         % VERSION)
        self.refuse("app", STUB_NO_IMAGE="1", contains="clinego build")
        self.refuse("run", "gea", "--config", "config_demo.yaml", "--image", "x:1",
                    STUB_NO_IMAGE="1", contains="No image x:1")

    def test_missing_docker(self):
        os.remove(os.path.join(self.bin, "docker"))
        self.refuse("run", "gea", "--config", "config_demo.yaml",
                    PATH=HelpAndVersionTest._path_without_docker(self), contains="docker not found")

    def test_bad_options(self):
        self.refuse("run", "gea", "--config", "config_demo.yaml", "--bogus", contains="--bogus")
        self.refuse("run", "gea", "--config", "config_demo.yaml", "--port", "1",
                    contains="--port")
        self.refuse("run", "gea", "--config", "config_demo.yaml", "--cores", "0",
                    contains="--cores")
        self.refuse("run", "gea", "--config", "config_demo.yaml", "--memory", "lots",
                    contains="--memory")
        self.refuse("run", "gea", "--config", "config_demo.yaml", "-n", contains="-n")
        self.refuse("run", "gea", "extra", "--config", "config_demo.yaml", contains="after --")
        self.refuse("run", "--config", "config_demo.yaml", contains="needs a mode")
        self.refuse("app", "--port", "70000", contains="--port")
        self.refuse("frobnicate", contains="unknown command")

    def test_a_config_name_that_cannot_name_a_container(self):
        with open(os.path.join(self.project, "my config.yaml"), "w") as handle:
            handle.write("project_name: x\n")
        self.refuse("run", "gea", "--config", "my config.yaml", contains="container")

    def test_build_refuses_in_the_image(self):
        self.refuse("build", in_image=True, contains="host")

    def test_example_refuses_to_overwrite(self):
        os.makedirs(os.path.join(self.project, "test_data"))
        self.refuse("example", contains="never overwrites")
        self.refuse("example", in_image=True, contains="never overwrites")

    def test_init_refuses_to_overwrite_and_bad_names(self):
        self.refuse("init", "demo", contains="never overwrites")
        self.refuse("init", "demo", in_image=True, contains="never overwrites")
        self.refuse("init", "../x", contains="letter or digit")
        self.refuse("init", "a/b", contains="letters, digits")
        self.refuse("init", contains="needs a project name")


class InImageTest(LauncherCase):
    """CLINEGO_IN_IMAGE=1: the stub snakemake prints its argv, one argument per line."""

    def argv(self, *args, **extra):
        result = self.call(*args, in_image=True, **extra)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.splitlines()

    def test_run_builds_the_interface_snakemake_call(self):
        self.assertEqual(self.argv("run", "processing", "--config", "config_demo.yaml"), [
            "-s", os.path.join(REPO_ROOT, "Snakefile"), "-c", "4",
            "--configfile", "config_demo.yaml", "--config", "mode=processing",
            "project_root=" + os.path.realpath(self.project), "--scheduler", "greedy",
            "--rerun-incomplete"])

    def test_run_forwards_cores_and_the_snakemake_args_last(self):
        argv = self.argv("run", "gwas", "--config", "config_demo.yaml", "--cores", "8",
                         "--", "-F", "--forcerun", "a rule")
        self.assertEqual(argv[2:4], ["-c", "8"])
        self.assertEqual(argv[argv.index("--rerun-incomplete") + 1:], ["-F", "--forcerun", "a rule"])

    def test_unlock_carries_unlock_and_a_mode(self):
        argv = self.argv("unlock", "--config", "config_demo.yaml")
        self.assertIn("--unlock", argv)
        self.assertNotIn("--rerun-incomplete", argv)
        self.assertIn("mode=processing", argv)
        self.assertIn("project_root=" + os.path.realpath(self.project), argv)

    def test_host_options_are_refused_in_the_image(self):
        self.refuse("run", "gea", "--config", "config_demo.yaml", "--image", "x:1",
                    in_image=True, contains="host option")
        self.refuse("app", "--port", "3871", in_image=True, contains="host options")

    def test_app_runs_the_package_never_dev_r(self):
        argv = self.show("app", "--show", in_image=True)
        self.assertEqual(argv, ["Rscript", "-e",
                                'clinego.app::run_app(host = "0.0.0.0", port = 3838)'])

    def test_init_writes_the_named_config_from_the_default(self):
        result = self.call("init", "demo2", in_image=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(os.path.join(self.project, "config_demo2.yaml")) as handle:
            text = handle.read()
        self.assertRegex(text, r"(?m)^project_name: demo2$")
        with open(os.path.join(REPO_ROOT, "scripts", "clinego.app", "inst",
                               "config_default.yaml")) as handle:
            default = handle.read()
        self.assertEqual(len(text.splitlines()), len(default.splitlines()))
        self.refuse("init", "demo2", in_image=True, contains="never overwrites")

    def test_example_copies_the_dataset_and_its_config(self):
        example = os.path.join(self.tmp, "example")
        os.makedirs(example)
        for name in ("testdata.vcf", "testdata_metadata.tsv", "config_testdata.yaml"):
            with open(os.path.join(example, name), "w") as handle:
                handle.write(name + "\n")
        result = self.call("example", in_image=True, CLINEGO_EXAMPLE_DIR=example)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sorted(os.listdir(os.path.join(self.project, "test_data"))),
                         ["testdata.vcf", "testdata_metadata.tsv"])
        self.assertTrue(os.path.isfile(os.path.join(self.project, "config_testdata.yaml")))
        self.refuse("example", in_image=True, CLINEGO_EXAMPLE_DIR=example,
                    contains="never overwrites")


if __name__ == "__main__":
    unittest.main()
