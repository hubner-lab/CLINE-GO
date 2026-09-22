#!/usr/bin/env bash
# Journal-15 figure set: the pooled SS-Clines corpus (5 blocks, 600 replicates),
# rendered into its OWN directory so the standing figures_ssclines_pooled/ set
# (241 files, subtitled, v1 composite as MAIN_FIGURE_simulation) is untouched.
#
#   MVP_SUBTITLE=0       short titles only; provenance goes into the journal captions
#   MVP_ROW_ORDER=strict panel-B row order must equal the pooled data-derived order
#   SCRIPTS              oracle_stats first (v2 reads its below_by_method.tsv),
#                        block_report before pleiotropy_report (runtime cross-checks),
#                        and NO mvp_main_figure.R, so MAIN_FIGURE_simulation is v2
#
#   bash work/j15_figures.sh          logs: work/blockrep_ssclines_pooled_j15/<script>.log
set -u
cd /mnt/data/eugene/ADAPTOGENE || exit 1
export RUN_TAG=ssclines_pooled_j15
export SCRIPTS="mvp_oracle_stats.R mvp_main_figure_v2.R mvp_block_report.R mvp_dist_panel.R mvp_pleiotropy_report.R"
export MVP_SUBTITLE=0 MVP_ROW_ORDER=strict
exec work/block_report.sh ssclines_pooled \
  ssclines_nvar_mvar,ssclines_ncline_ns,ssclines_nequal_mconst,ssclines_ncline_ctredge,ssclines_nequal_mbreaks \
  600 200 offset12_ssclines_pooled
