# benchmarks/land: scripts behind the three-landscape simulation results

Copies of files that lived only in gitignored working folders, so that the paper's claim of scripts in the repository holds. They are records of what ran; paths inside them still name `work/...`.

- `scratch/`: the per-landscape copies of the scoring and figure scripts, from `work/mvp_land/scratch/`, made by `make_scratch_land.py` and `make_figs_land.py` from the committed scripts, each patch asserted.
- `closeout/`: the builders of the published SS-Clines tables and figures, from `work/closeout_lab/scratch/new/`. `benchmarks/mvp_ms_*.R` predate the LD-window rerun and do not rebuild those items.
- `dk.sh`, `post.sh`, `figs.sh`: the pinned launcher and the per-landscape scoring and figure drivers, from `work/mvp_land/`.

The current builders of the manuscript's tables and figures for all three landscapes are `benchmarks/mvp_land_tables.R` and `benchmarks/mvp_land_figures.R`.
