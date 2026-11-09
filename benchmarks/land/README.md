# benchmarks/land: scripts behind the three-landscape simulation results

Copies of files that lived only in gitignored working folders, so that the paper's claim of scripts in the repository holds. They are records of what ran; paths inside them still name `work/...`.

- `scratch/`: the per-landscape copies of the scoring and figure scripts, from `work/mvp_land/scratch/`, made by `make_scratch_land.py` and `make_figs_land.py` from the committed scripts, each patch asserted.
- `closeout/`: the builders of the published SS-Clines tables and figures, from `work/closeout_lab/scratch/new/`. `benchmarks/mvp_ms_*.R` predate the LD-window rerun and do not rebuild those items.
- `dk.sh`, `post.sh`, `figs.sh`: the pinned launcher and the per-landscape scoring and figure drivers, from `work/mvp_land/`.

The current builders of the manuscript's tables and figures for all three landscapes are `benchmarks/mvp_land_tables.R` and `benchmarks/mvp_land_figures.R`.

## Data S1 landscapes job (2026-11-08, `f827ce4`)

The seven scripts of the lab's Data S1 job (`docs/manuscript/lab/data_s1_landscapes_LAB_PROMPT.md`; report `benchmarks/mvp_eval/landscape_ms/data_s1_2026-11/REPORT.md`, gitignored): `benchmarks/mvp_land_data_s1.py`, `mvp_land_causal_tested.py`, `mvp_land_data_s1_check_tables.py`, `mvp_land_data_s1_check_strata.py`, `mvp_land_data_s1_trifolium.R`, `mvp_land_figures.R` (`FIG_MODE=merge`) and `mvp_land_tracking_fig.R` (`FP_PCA_OUT`).

Adopted copies (sha256, first 16 hex digits, equal to the file that ran):

- `closeout/build_paper_pack_d150.R` from `work/closeout_lab/results/scripts/patched/`: `afe86083b89bc5c9` (Job F4, panel a of `s_random_sets_climate`).
- `closeout/plot_random_climate_ms_d150.R` from the same folder: `e00ff7378d149189` (Job F4, panel b).

Both are the same bytes as the `work/closeout_lab/scratch/new/` copies. `mvp_land_data_s1_trifolium.R` ran them with the `MINOU` environment of `work/closeout_lab/results/scripts/tools/run_closeout.sh` (steps s17 and s18; that launcher is not in the repository): `AN=ANALYSIS/Trifolium86MainK5d150/fig3/offset_minou PROJ=Trifolium86MainK5d150 PALETTE=minou SITE_NAMES=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/site_names.tsv KCHOICE_PANELS=a DELTA_CORR=1 NO_LINEAR=1 MAIN_MS=1 SUPP_MS=1 MAIN6_MS=1 SUPP2_MS=1`.

Tracked files found byte-identical to what ran (each last changed before the job):

- `benchmarks/ms_style.R` (sourced by the two R scripts, the fig3 scripts and both closeout copies) and `benchmarks/mvp_rank_ties.R` (sourced by `mvp_land_figures.R`).
- `manuscript/trifolium/Trifolium86MainK5ERpc5/fig3/plot_pca_structure_ms.R` and `plot_structure_ms.R`, run from that path (Job F3).
- `manuscript/trifolium/Trifolium86MainK5ERpc5/fig3/fig3_style.R`: the fig3 scripts and `build_paper_pack_d150.R` source `ANALYSIS/Trifolium86MainK5ERpc5/fig3/fig3_style.R`, which has the same bytes (`9711d65aec39aecf`), so no copy was made.
- `closeout/sim_figures_ld.R` and `closeout/rda_correction_ld.R` (above), named in `mvp_land_figures.R` as the source of its copied functions, equal the `patched/` copies of the same name.

Data the scripts read, which stays out of the repository (only `benchmarks/mvp_seeds.tsv` is tracked):

- `mvp_land_data_s1.py`: `benchmarks/mvp_seeds.tsv`; `work/mvp_land/seeds_{ssmtn,estclines}.txt`, `ld_window_by_replicate.tsv` and `config_backup_ld619/`; `.snakemake/log/`; the root `config_MVP{seed}{hd,hdrc,hd_sweep}.yaml`; the `MVP{seed}hd` and `hdrc` trees (`snp_sets/`, `ld_decay_half_distances.tsv`); `work/ld_window_lab/results/phase8/data_s1` and its `data_s1_sha256.txt`; `_archive/mvp_ssclines/configs/`.
- `mvp_land_causal_tested.py`: `benchmarks/mvp_seeds.tsv`, `work/mvp_land/seeds_*.txt`, the hd trees' p-value tables (SS-Clines from `_archive/mvp_ssclines/results/`), `data/mvp/MVP{seed}/{loci_freq,truth_any}.tsv`.
- `mvp_land_data_s1_check_tables.py`, `mvp_land_data_s1_check_strata.py`: `benchmarks/mvp_eval/landscape_ms/tables/` and the job's own outputs.
- `mvp_land_figures.R`, `mvp_land_tracking_fig.R`: the scored tables under `benchmarks/mvp_eval/` (galleries, rank and remeasure folders, `landscape_ms/tables/`, `landscape_ms/explain/`).
- `mvp_land_data_s1_trifolium.R`: `Trifolium86MainK5ERpc5_results` (PCA projections and eigenvalues), `ANALYSIS/Trifolium86MainK5ERpc5/fig3/inputs/clusters_K5.short.tsv` and `site_names.tsv`, the 150 kb offset run `ANALYSIS/Trifolium86MainK5d150/fig3/offset_minou/` and `Trifolium86MainK5d150_results`.
