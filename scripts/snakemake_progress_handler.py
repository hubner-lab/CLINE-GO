#!/usr/bin/env python3
"""
Snakemake log handler script for CLINE-GO pipeline progress tracking.

Used via: snakemake --log-handler-script scripts/snakemake_progress_handler.py

Receives structured log messages from Snakemake and writes JSON-lines to
{PROJECT}_logs/.progress.jsonl so the Shiny app can poll for progress.

Message format received from Snakemake:
    {"level": "job_info", "msg": "...", "name": "rule_name", "jobid": N, ...}

Output JSON-lines written:
    {"timestamp": "...", "level": "job_info",     "rule": "name", "jobid": N, "total_jobs": M}
    {"timestamp": "...", "level": "job_finished", "rule": "name", "jobid": N}
    {"timestamp": "...", "level": "job_error",    "rule": "name", "jobid": N}
    {"timestamp": "...", "level": "dag_info",     "total_jobs": M}
"""

import json
import os
import sys
from datetime import datetime

def get_progress_file():
    """Determine the .progress.jsonl path from the working directory."""
    # Snakemake sets the working directory to the Snakefile directory.
    # Results are in {PROJECT}_results/, logs in {PROJECT}_logs/.
    # We write to the pipeline root's first project logs dir we can find.
    # Fallback: write to /tmp/.progress.jsonl
    wd = os.getcwd()  # /pipeline when running inside Docker
    # Try to infer project from configfile argument in sys.argv
    project = None
    for i, arg in enumerate(sys.argv):
        if arg == "--configfile" and i + 1 < len(sys.argv):
            cf = sys.argv[i + 1]
            basename = os.path.basename(cf)
            if basename.startswith("config_") and basename.endswith(".yaml"):
                project = basename[len("config_"):-len(".yaml")]
            elif basename == "config.yaml":
                project = None
            break

    if project:
        log_dir = os.path.join(wd, f"{project}_logs")
    else:
        # Scan for any *_logs directory
        candidates = [d for d in os.listdir(wd)
                      if d.endswith("_logs") and os.path.isdir(os.path.join(wd, d))]
        log_dir = os.path.join(wd, candidates[0]) if candidates else "/tmp"

    os.makedirs(log_dir, exist_ok=True)
    return os.path.join(log_dir, ".progress.jsonl")


_progress_file = None
_total_jobs = None
_truncated = False

def _get_file():
    global _progress_file, _truncated
    if _progress_file is None:
        _progress_file = get_progress_file()
    if not _truncated:
        # Truncate at start of run so stale data doesn't confuse the app
        with open(_progress_file, "w") as f:
            pass
        _truncated = True
    return _progress_file


def _write(record):
    record["timestamp"] = datetime.utcnow().isoformat() + "Z"
    path = _get_file()
    with open(path, "a") as f:
        f.write(json.dumps(record) + "\n")
        f.flush()


def main(msg):
    """Entry point called by Snakemake for each log event."""
    global _total_jobs

    level = msg.get("level", "")

    if level == "job_info":
        # Snakemake fires job_info before each job starts
        rule = msg.get("name", msg.get("rule", ""))
        jobid = msg.get("jobid")
        log_files = msg.get("log", [])
        log_file = log_files[0] if log_files else None

        # Extract total job count from the message text if present
        # Snakemake includes "N of M" in the msg string
        import re
        text = msg.get("msg", "")
        m = re.search(r"(\d+)\s+of\s+(\d+)", text)
        if m:
            _total_jobs = int(m.group(2))

        record = {"level": "job_info", "rule": rule, "jobid": jobid, "log_file": log_file}
        if _total_jobs is not None:
            record["total_jobs"] = _total_jobs
        _write(record)

    elif level == "job_finished":
        rule = msg.get("name", msg.get("rule", ""))
        jobid = msg.get("jobid")
        log_files = msg.get("log", [])
        log_file = log_files[0] if log_files else None
        _write({"level": "job_finished", "rule": rule, "jobid": jobid, "log_file": log_file})

    elif level == "job_error":
        rule = msg.get("name", msg.get("rule", ""))
        jobid = msg.get("jobid")
        log_files = msg.get("log", [])
        log_file = log_files[0] if log_files else None
        _write({"level": "job_error", "rule": rule, "jobid": jobid, "log_file": log_file})

    elif level == "dag_info":
        # dag_info fires after DAG construction with total job count
        total = msg.get("total_jobs")
        if total is not None:
            _total_jobs = int(total)
            _write({"level": "dag_info", "total_jobs": _total_jobs})

    elif level in ("progress", "run_info"):
        # progress: Snakemake sometimes emits these with counts
        total = msg.get("total")
        if total is not None:
            _total_jobs = int(total)
            _write({"level": "dag_info", "total_jobs": _total_jobs})
