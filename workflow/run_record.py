"""Run records: what a mode ran with, kept in the results tree it ran on.

One JSON file per mode, ``{PROJECT}_results/_intermediate/run_records/<mode>.json``, rewritten
by every Snakemake invocation of that mode that has jobs to run, and one line per finished
invocation appended to ``history.jsonl`` beside it (findings 665953, 96ccfb). A results tree
copied with ``cp -a`` keeps both, so the config its upstream modes read, the command, the
pipeline commit and the image stay readable from the tree alone.

The Snakefile's onstart / onsuccess / onerror handlers call common.smk's write_run_record(),
which calls write() here. Snakemake runs those handlers only for a real run with jobs: never
for ``-n``, ``--no-hooks`` or "Nothing to be done". A dry run therefore writes nothing, and a
no-op invocation keeps the record of the run that built the outputs.

A record holds: mode; status (running, then success or error); started and finished, UTC;
the command line; each --configfile with the sha256 of its bytes; the config Snakemake ran
with (the --configfile values merged with --config); the git commit of the checkout holding
the Snakefile, whether it was dirty, the changed paths and a sha256 of ``git diff HEAD``; the
image ID when it is visible; Snakemake and Python versions; the resolved cpu; the mode's
GEA/GWAS method parameters; the path of Snakemake's own log.

It never holds an author, committer or e-mail address, and reads exactly one environment
variable, CLINEGO_IMAGE_ID: a container cannot see its own image, so the ID is visible only
when the caller passes it (``docker run -e CLINEGO_IMAGE_ID=...``).

JSON from the standard library: the image's system python3, which runs the Python tests and
loads common.smk through tests/python/_smk_harness.py, has no PyYAML. Non-table text on
purpose: two otherwise identical trees differ in dates and commits, and
scripts/compare_results.R reports a differing non-table file as changed_bytes, never as a
table difference.
"""
import datetime
import hashlib
import json
import os
import platform
import shlex
import subprocess

FORMAT = "CLINE-GO run record 1"
IMAGE_ENV = "CLINEGO_IMAGE_ID"
RECORD_DIR = "run_records"
HISTORY = "history.jsonl"
MAX_CHANGED_PATHS = 200
_CONTAINER_MARKERS = ("/.dockerenv", "/run/.containerenv")


def utc_now():
    """Now, UTC, to the second, ISO 8601."""
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()


def snakemake_version():
    """The running Snakemake's version; None where snakemake is not importable."""
    try:
        from snakemake import __version__
    except Exception:
        return None
    return __version__


def file_sha256(path):
    """sha256 of a file's bytes."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git(repo, *args, timeout=60):
    """stdout (bytes) of one git command run in `repo`, never writing its index: the record
    runs at every mode start against the user's live checkout, as whatever user the container
    runs as. --no-optional-locks stops `git status` refreshing .git/index under index.lock;
    porcelain `git diff` refreshes it even so (git 2.43), hence diff-index in git_state().
    safe.directory=* because the container user need not own the mounted checkout."""
    return subprocess.run(["git", "--no-optional-locks", "-c", "safe.directory=*", "-C",
                           str(repo), *args],
                          capture_output=True, timeout=timeout, check=True).stdout


def git_state(repo):
    """The commit of the checkout at `repo`, whether it is dirty, the changed paths
    (``git status --porcelain``, untracked included) and a sha256 of the uncommitted diff
    (``git diff-index -p HEAD``: the bytes ``git diff HEAD`` prints, without its index
    refresh). Never author or committer data. No git, no checkout or any git error gives
    commit None and the reason under `note`."""
    try:
        commit = _git(repo, "rev-parse", "HEAD").decode().strip()
        status = _git(repo, "status", "--porcelain", "--untracked-files=normal").decode()
        diff = _git(repo, "diff-index", "-p", "--no-color", "--no-ext-diff", "HEAD")
    except FileNotFoundError:
        return {"commit": None, "dirty": None, "note": "git not found"}
    except subprocess.CalledProcessError as exc:
        lines = exc.stderr.decode(errors="replace").strip().splitlines()
        return {"commit": None, "dirty": None,
                "note": "git: " + (lines[-1] if lines else f"exit status {exc.returncode}")}
    except (subprocess.TimeoutExpired, OSError) as exc:
        return {"commit": None, "dirty": None, "note": f"git: {type(exc).__name__}"}
    changed = [line[3:] for line in status.splitlines() if line.strip()]
    state = {"commit": commit, "dirty": bool(changed)}
    if changed:
        state["n_changed_paths"] = len(changed)
        state["changed_paths"] = changed[:MAX_CHANGED_PATHS]
        state["diff_sha256"] = hashlib.sha256(diff).hexdigest()
    return state


def image_state(image_id, markers=_CONTAINER_MARKERS):
    """The image the run used, when it is visible (`image_id`, the value of
    CLINEGO_IMAGE_ID), and whether the run is inside a container."""
    state = {"in_container": any(os.path.exists(m) for m in markers)}
    image_id = (image_id or "").strip()
    if image_id:
        state.update(id=image_id, visible=True)
    else:
        state.update(id=None, visible=False,
                     note=f"not visible from inside a container; pass "
                          f"-e {IMAGE_ENV}=$(docker images --no-trunc -q <image>)")
    return state


def _configfile_entry(path):
    try:
        return {"path": str(path), "sha256": file_sha256(path)}
    except OSError as exc:
        return {"path": str(path), "sha256": None, "note": f"{type(exc).__name__}: {exc}"}


def _plain(value):
    """Plain JSON types: Snakemake's config may hold Paths or immutables maps."""
    return json.loads(json.dumps(value, default=str))


def build_record(*, mode, status, started, finished, argv, launch_dir, configfiles, config,
                 repo, image_id, snakemake_version, cpu, methods=None, snakemake_log=()):
    """One run record as a dict, in the order it is written."""
    record = {
        "format": FORMAT,
        "mode": mode,
        "status": status,
        "started": started,
        "finished": finished,
        "command": shlex.join([os.path.basename(argv[0]), *argv[1:]]) if argv else None,
        "launch_dir": str(launch_dir) if launch_dir else None,
        "configfiles": [_configfile_entry(p) for p in configfiles],
        "git": git_state(repo),
        "image": image_state(image_id),
        "software": {"snakemake": snakemake_version, "python": platform.python_version()},
        "cpu": cpu,
        "snakemake_log": [str(p) for p in snakemake_log],
    }
    if methods:
        record["methods"] = _plain(methods)
    record["config"] = _plain(config)
    return record


def history_entry(record):
    """The one-line summary history.jsonl keeps of a finished run."""
    git = record.get("git") or {}
    entry = {
        "started": record.get("started"),
        "finished": record.get("finished"),
        "mode": record.get("mode"),
        "status": record.get("status"),
        "commit": git.get("commit"),
        "dirty": git.get("dirty"),
    }
    if git.get("diff_sha256"):
        entry["diff_sha256"] = git["diff_sha256"]
    entry["configfiles"] = {c["path"]: c["sha256"] for c in record.get("configfiles", [])}
    entry["image"] = (record.get("image") or {}).get("id")
    entry["command"] = record.get("command")
    return entry


def _atomic_write(path, text):
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as handle:
        handle.write(text)
    os.replace(tmp, path)


def write(intermediate_dir, record):
    """Write <mode>.json under <intermediate_dir>/run_records/, atomically, and append one
    history.jsonl line when the run has finished. Returns the record's path."""
    directory = os.path.join(intermediate_dir, RECORD_DIR)
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, f"{record['mode']}.json")
    _atomic_write(path, json.dumps(record, indent=2, ensure_ascii=False) + "\n")
    if record.get("status") != "running":
        with open(os.path.join(directory, HISTORY), "a") as handle:
            handle.write(json.dumps(history_entry(record), ensure_ascii=False) + "\n")
    return path
