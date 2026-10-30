# CLINE-GO Pipeline
# vim: filetype=python
#
# Modular Snakemake workflow for population genomics.
# Config, paths, and helpers are in workflow/rules/common.smk.
# Rules are organized by pipeline mode in separate .smk files.

include: "workflow/rules/common.smk"
include: "workflow/rules/processing.smk"
include: "workflow/rules/prestructure.smk"
include: "workflow/rules/climate.smk"
include: "workflow/rules/traits.smk"
include: "workflow/rules/pregea.smk"
include: "workflow/rules/structure.smk"
include: "workflow/rules/gea.smk"
if GWAS_CONFIGS:
    include: "workflow/rules/gwas.smk"
if ASSOC_SOURCES:
    include: "workflow/rules/_assoc_downstream.smk"
if GWAS_CONFIGS and GEA_CONFIGS:
    include: "workflow/rules/gea_x_gwas.smk"
include: "workflow/rules/maladaptation.smk"
include: "workflow/rules/summary.smk"

# Job temp files (R's tempdir(), fread(cmd = ) copies, mktemp) go to the results tree, not to
# the container's /tmp, which is the host's root disk (441eee). Snakemake exports its tmpdir
# resource to every job as TMPDIR, and that resource defaults to the tempfile module's
# directory, set here. A TMPDIR already in the environment (an HPC job's scratch) or
# --default-resources "tmpdir='/path'" still wins. Set here, not in common.smk: the Python
# tests exec that file in their own process.
import tempfile
JOB_TMPDIR = f"{INTER}tmp"
if not any(os.environ.get(v) for v in ("TMPDIR", "TEMP", "TMP")):
    tempfile.tempdir = JOB_TMPDIR

# Run records: _intermediate/run_records/<mode>.json + history.jsonl (write_run_record in
# common.smk). Snakemake runs these only for a real run with jobs, never for -n, which is
# also why the job temp folder is made here: a dry run writes nothing.
onstart:
    write_run_record("running", log)
    if tempfile.tempdir == JOB_TMPDIR:
        os.makedirs(JOB_TMPDIR, exist_ok=True)

onsuccess:
    write_run_record("success", log)

onerror:
    write_run_record("error", log)

rule all:
    input: get_targets(MODE)
