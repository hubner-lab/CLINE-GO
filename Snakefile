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

rule all:
    input: get_targets(MODE)
