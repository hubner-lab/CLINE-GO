# SS-Clines complete-block runbook

Run **one demography block at a time**. Each block is 120 replicates and a complete,
self-contained sub-design (3 genic levels x 4 architecture sub-levels x 10 replicates), so a
block that finishes is analysable on its own and the run can stop between blocks.

Sized for **160 CPUs / 800 GB**. Estimated **~13 h per block**, ~65 h for all five.

Scope, fixed and not re-decided per block: landscape `SS-Clines`, **2-trait only** (1-trait
excluded on identifiability grounds), **no threshold on any realized statistic** — the block is
complete, so there is no selection rule inside it. Confounding is reported, never selected on.

---

## Block table

Run in this order. Block 1 is the most informative (messiest demography, and the legacy SS-Mtn
corpus has 45 replicates of the same demography for a direct landscape contrast).

| # | `MVP_BLOCK` | `MVP_COHORT` | `MVP_COHORT_TAG` | cell manifest |
|---|---|---|---|---|
| 1 | `SS-Clines:N-variable_m-variable` | `ssclines_b1` | `ssclines_nvar_mvar` | `mvp_sweep_cells_ssclines_nvar_mvar.tsv` |
| 2 | `SS-Clines:N-cline-N-to-S_m-constant` | `ssclines_b2` | `ssclines_ncline_ns` | `mvp_sweep_cells_ssclines_ncline_ns.tsv` |
| 3 | `SS-Clines:N-equal_m-constant` | `ssclines_b3` | `ssclines_nequal_mconst` | `mvp_sweep_cells_ssclines_nequal_mconst.tsv` |
| 4 | `SS-Clines:N-cline-center-to-edge_m-constant` | `ssclines_b4` | `ssclines_ncline_ctredge` | `mvp_sweep_cells_ssclines_ncline_ctredge.tsv` |
| 5 | `SS-Clines:N-equal_m_breaks` | `ssclines_b5` | `ssclines_nequal_mbreaks` | `mvp_sweep_cells_ssclines_nequal_mbreaks.tsv` |

All five cell-manifest filenames are **already registered** in `mvp_build_snp_sets.R`
(the list is filtered by `file.exists`, so naming them before they exist is a no-op).

---

## Set this once per block

```bash
cd /mnt/data/eugene/ADAPTOGENE

# ---- edit these four lines per block, everything below is generic ----
export MVP_BLOCK="SS-Clines:N-variable_m-variable"
export MVP_COHORT="ssclines_b1"
export MVP_COHORT_TAG="ssclines_nvar_mvar"
export CELLS_TSV="benchmarks/mvp_sweep_cells_ssclines_nvar_mvar.tsv"
# ----------------------------------------------------------------------

export MVP_BLOCK_ARM="primary_ssclines"      # keep constant across all five blocks
export OFFSET_DIR="benchmarks/mvp_eval/offset12"
export DK="nix shell nixpkgs#docker-client -c docker"
export UIDGID="$(id -u):$(id -g)"
```

`primary_ssclines` is deliberate: every existing figure script filters `arm == "primary"`, so
the new landscape stays invisible to them until a script opts in. Do not change it to
`primary` — that silently mixes two landscapes into every marginal.

---

## Preflight

```bash
df -h /mnt/data | tail -1                       # need >= 0.5 TB free per block
uptime                                          # load should leave ~160 cores
free -g | head -2                               # need >= 800 GB available
$DK ps --format '{{.Names}}\t{{.Image}}'        # note whose containers are already running
cp benchmarks/mvp_seeds.tsv benchmarks/mvp_seeds.tsv.pre_${MVP_COHORT}
```

**Never** `docker stop $(docker ps -q --filter ancestor=cline-go:latest)` — that matches every
container sharing the image, including other users'. Always `--name` your containers and stop
them individually.

---

## Step 0 — manifest

```bash
$DK run --rm --name mvp-selseeds-${MVP_COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g \
  -e MVP_BLOCK -e MVP_COHORT -e MVP_COHORT_TAG -e MVP_BLOCK_ARM \
  -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_select_seeds.R
```

**Gate.** The log must contain all three:
```
block holds 120 seeds, 0 already in the manifest, adding 120
== 92 existing rows reproduced exactly (22 columns checked) ==     # 212, 332, ... on later blocks
== append-only check passed: N kept + 120 new = N+120 ==
```
Anything else — stop and restore `benchmarks/mvp_seeds.tsv.pre_${MVP_COHORT}`.

Capture this block's seed list (column 22 is `added`):

```bash
export SEEDS=$(awk -F'\t' -v t="$MVP_COHORT_TAG" 'NR>1 && $22==t {print $1}' \
               benchmarks/mvp_seeds.tsv | paste -sd,)
echo "$SEEDS" | tr ',' '\n' | wc -l      # must print 120
```

## Step 1 — fetch

Streams the deposit's three tars (~23.5 GB over the network) and extracts only this block's
members (~500 MB lands). Idempotent; a tar whose members are all present is skipped.

```bash
benchmarks/fetch_mvp.sh 2>&1 | tee "logs_fetch_${MVP_COHORT}.log"
```

**Gate.** `ls data/mvp/raw/genotypes | wc -l` grew by 120; same for `mutations`, `individuals`.
The script verifies by file presence, not tar exit code, so check the count yourself.

## Step 2 — convert

**Pass the seeds explicitly. Never `all`.** The `all` default rewrites finished replicates'
`data/mvp/MVP{seed}/` files, which are Snakemake inputs — a bumped mtime re-runs `processing`
and everything downstream on replicates that were already finished and scored.

```bash
benchmarks/convert_mvp_all.sh 20 "$SEEDS" 2>&1 | tee "logs_convert_${MVP_COHORT}.log"
```

**Gate.** `data/mvp/conversion_summary.tsv` gained 120 rows, zero failures.

## Step 3 — base configs

```bash
$DK run --rm --name mvp-cfg-${MVP_COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_write_configs.R --seeds="$SEEDS"
```

**Gate.** 120 new `config_MVP{seed}.yaml` at the repo root.

## Step 4 — sweep config, one cell

`--cells=1` now means "the default rung only" (`npc = lfmm_K = k_best`), which is the one cell
marker-panel construction reads. Extra cells are LFMM-K / #PC hyperparameter ladders and feed
only the detection arm, which this corpus does not extend.

```bash
$DK run --rm --name mvp-swcfg-${MVP_COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_write_sweep_configs.R \
    --seeds="$SEEDS" --cells=1 --manifest-out="$CELLS_TSV"
```

**Gate — do this before spending any compute:**

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} $h["is_default"]=="TRUE" {n++} END{print n" TRUE rows"}' "$CELLS_TSV"
```
Must equal `120 x (methods per seed)`, i.e. exactly one `TRUE` per (seed, method). A count of
**0** means the ladder patch is missing — stop, do not run the sweep.

## Step 5 — GEA sweep  (~7 h)

Modes per seed, strictly serial: `processing -> prestructure -> structure -> pregea -> gea`.

```bash
CPUS_PER_SEED=8 SNAKE_CORES=4 BLAS_THREADS=2 MEM_PER_SEED=40g \
  benchmarks/mvp_run_sweep.sh "$SEEDS" 20 1 2>&1 | tee "logs_gea_${MVP_COHORT}.log"
```

Positional args are `SEEDS_CSV SEED_JOBS CELLS` — concurrency is the `20`, not an env var.
`CPUS_PER_SEED`, `SNAKE_CORES`, `MEM_PER_SEED`, `BLAS_THREADS` are env. 20 seeds x 8 cores = 160.

`BLAS_THREADS=2` is **required, not tuning** — containers see the host's 192 cores and OpenBLAS
threads to that unless pinned, which makes the work ~40x slower and eats ~80 GB.

**Gate.** 120 `DONE` lines; `ls benchmarks/mvp_eval/params/MVP*/c1 -d | wc -l` grew by 120.

## Step 6 — SNP panels

Writes every set including `all`. That one is only a SNP *list* — free — and five figure
scripts use it as the precision/recall denominator via `panel_pr_recomputed.tsv`. The expensive
part is the Gradient Forest fit, which is excluded at step 8 instead.

```bash
$DK run --rm --name mvp-snpsets-${MVP_COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=64g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_build_snp_sets.R --seeds="$SEEDS"
```

**Gate — pre-screen panels at n >= 3 before going further.** `scripts/rda_offset.R:174` FATALs
below 3 SNPs, and that rule failing aborts the *whole seed's* Snakemake run, killing every panel
for that replicate. Find the offenders now, not after the garden hours:

```bash
for s in ${SEEDS//,/ }; do
  for p in truth union best intersect3 rand_best1 solo_lfmm solo_rda solo_emmax; do
    f="MVP${s}_results/_intermediate/snp_sets/${p}/selected_snps.tsv"
    if [[ -f "$f" ]]; then n=$(( $(wc -l < "$f") - 1 )); else n=0; fi
    (( n < 3 )) && echo -e "MVP${s}\t${p}\t${n}"
  done
done | tee "panel_underfilled_${MVP_COHORT}.tsv"
```
Any row listed here: drop that panel from `--sets` for that seed at step 8, and record it in
the block's exclusion note. An empty panel is a RESULT (`n = 0`, `recall = 0`), not missing data.

## Step 7 — garden fitness truth, then garden env tables

Order matters: `mvp_garden_fitness.R` writes `gardens_{seed}.tsv`, which
`mvp_write_garden_env.R` needs.

```bash
$DK run --rm --name mvp-fitness-${MVP_COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS=4 --cpus=16 --memory=96g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_garden_fitness.R --seeds="$SEEDS" --outdir=/pipeline/$OFFSET_DIR

for s in ${SEEDS//,/ }; do
  $DK run --rm --name mvp-genv-$s --user "$UIDGID" -e USER=adaptogene \
    -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g -v "$PWD":/pipeline cline-go:latest \
    Rscript /pipeline/benchmarks/mvp_write_garden_env.R --seed=$s \
      --gardens=/pipeline/$OFFSET_DIR/gardens_${s}.tsv
done
```

Do **not** pass `--outdir` to `mvp_write_garden_env.R`. Its default is
`data/mvp/MVP{seed}/gardens/`, which is exactly where the sweep configs look; redirecting it
to `$OFFSET_DIR` writes the env tables somewhere nothing reads. `--gardens` **must** be passed,
because it defaults to `offset09/gardens_{seed}.tsv` and this block's fitness tables are in
`offset12`.

**Gate.** 120 each of `gardens_{seed}.tsv` and `garden_fitness_{seed}.tsv` in `$OFFSET_DIR`;
`data/mvp/MVP{seed}/gardens/` holds 112 tables per seed.

## Step 8 — garden sweep configs, panels dropped here

`all` and `neutral_all` are ~3.5 of the 7.2 GB per seed and 92-113 of the ~113 GB peak RSS.
Dropping them is what makes 120 replicates affordable. `rand_best1` is **kept** — it is
size-matched to `best`, so it costs ~24 MB, and it is the floor that makes the headline
interpretable ("beats a size-matched random draw").

```bash
$DK run --rm --name mvp-gcfg-${MVP_COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS=2 --cpus=2 --memory=8g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_write_sweep_config.R --seeds="$SEEDS" \
    --sets=truth,union,best,intersect3,rand_best1,solo_lfmm,solo_rda,solo_emmax
```

**Gate.** `grep -c 'neutral_all\|^\s*- all$' config_MVP*_sweep.yaml` returns nothing for the new
seeds; `grep -A2 emit_plots config_MVP<one new seed>_sweep.yaml` shows `false`.

## Step 9 — garden sweep  (~5 h)

`mode=maladaptation`, 112 gardens, one Snakemake invocation per seed.

```bash
OUT="$PWD/$OFFSET_DIR" LOGDIR="$PWD/$OFFSET_DIR/sweeplogs_${MVP_COHORT}" \
CONFIG_SUFFIX=_sweep BLAS_THREADS=2 HOST_RAM_CEILING_GB=800 \
  benchmarks/mvp_garden_sweep.sh "$SEEDS" 20 8 40g 2>&1 | tee "logs_garden_${MVP_COHORT}.log"
```

Positional args are `SEEDS_CSV SEED_JOBS SNAKE_CORES MEM_PER_SEED`. **Never pass `all`** — it
takes every seed in the manifest, including the finished legacy 92.

**Gate.**
```bash
# 8 panels x 4 methods x 112 gardens = 3584 files per seed
for s in ${SEEDS//,/ }; do
  n=$(find "$OFFSET_DIR/gardens/$s" -name '*.tsv' 2>/dev/null | wc -l); echo "$s $n"
done | awk '$2!=3584 {print "SHORT: "$0}'
```
Empty output = complete. Seeds with dropped panels from step 6 will be short by
`4 x 112 = 448` per dropped panel — expected, cross-check against
`panel_underfilled_${MVP_COHORT}.tsv`.

## Step 10 — score

```bash
$DK run --rm --name mvp-score-${MVP_COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS=8 --cpus=16 --memory=96g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/eval_offset_lind.R --seeds="$SEEDS" --outdir=/pipeline/$OFFSET_DIR
```

**Gate.** `garden_performance.tsv` gained `120 x 112 x panels x 4` rows, zero NA. Check
`scoring_skipped.tsv` — anything there that is not in `panel_underfilled_${MVP_COHORT}.tsv` is a
real problem.

Sign convention is Lind & Lotterhos': a good model gives a **negative** tau. Flip once with
`-tau`, never `abs(tau)`, and never pool the four offset methods in a paired test.

## Step 11 — tables, with the regression gate first

```bash
$DK run --rm --name mvp-tables-${MVP_COHORT} --user "$UIDGID" -e USER=adaptogene \
  -e OPENBLAS_NUM_THREADS=4 --cpus=8 --memory=64g -v "$PWD":/pipeline cline-go:latest \
  Rscript /pipeline/benchmarks/mvp_panel_tables.R \
    --outdir=/pipeline/$OFFSET_DIR --check=/pipeline/benchmarks/mvp_eval/offset11
```

**Gate.** `--check` must report **max abs diff 0** against `offset11`. Non-zero means a scoring
change leaked in and the legacy numbers moved — stop and diagnose before trusting anything new.

## Step 12 — block report

**None of the existing figure scripts will report this block.** They all filter
`arm == "primary"`, and `mvp_oracle_stats.R:84-86` additionally asserts `nrow(PRIM) == 90L`.
That is the intended consequence of the `primary_ssclines` arm — the legacy analysis keeps
working untouched and stays at exactly 90 replicates, and nothing silently absorbs a second
landscape. It also means a **block-report script does not exist yet** and has to be written
before step 12 can run. Until it does, step 12 is:

```bash
# placeholder: block report script pending (see "Still to write" below)
```

Record the block's realized covariates — these are what the paper reports **instead of**
selecting on them: `|tau|(PC1,temp)`, `|tau|(PC1,Env2)`, `meanFst`, `K`, `final_LA`, `n_snps`.

Then append the block outcome to `~/Orthidian/projects/ADAPTOGENE/ADAPTOGENE-simulation.md`
and start the next block.

---

## Still to write (not blocking steps 0-11)

1. **`benchmarks/mvp_block_report.R`** — the block report. Must read
   `offset12/garden_performance.tsv` and the manifest, filter `arm == "primary_ssclines"` and
   the block's `added` tag, and emit: per-cell (genic x sub-level) medians, the realized
   covariate table, and the accuracy-vs-`final_LA` panel. Needed for step 12.
2. **`mvp_lind_compare.R:74-98`** — errors rather than degrades when `all` / `neutral` are
   absent from `garden_performance.tsv` (bare `all` resolves to `base::all`). Guard it before
   running that script against `offset12`, or skip it explicitly.
3. **Panel-denominator authority** — the step-6 `n >= 3` pre-screen makes missing panels a
   *config-time* omission for this arm, while the legacy arm records them *post-hoc* in
   `panel_offset_exclusions.tsv`. Decide which file is authoritative before any cross-arm
   median is computed, or it mixes 90-seed and 120-seed denominators.
4. **Docs** — `docs/methods_simulations.md` and `docs/gea_simulation_benchmarks.md` §8.2 still
   state `R² ∈ [0.30, 0.60]` as the selection rule, which is false for 46 of the legacy 90.
   Restate: legacy arm = band-selected (record the true band, `[0.20,0.75]`), SS-Clines arm =
   complete blocks, no within-scope selection. Also note the deposit ships **Kendall τ**, so
   the "R²" columns are τ², not variance explained.

## Known traps

| trap | consequence | guard |
|---|---|---|
| `--cells=1` before the ladder patch | `is_default` FALSE on every row, `default_cell()` resolves nothing, panels empty after ~7 h | step 4 gate |
| cell manifest not registered in `mvp_build_snp_sets.R` | same empty-panel failure | all five already registered |
| `convert_mvp_all.sh all` | rewrites finished replicates' Snakemake inputs, re-runs everything | pass `$SEEDS` |
| `mvp_garden_sweep.sh all` | re-runs the legacy 92 | pass `$SEEDS` |
| panel with < 3 SNPs | `rda_offset.R:174` FATAL aborts the **whole seed**, not just that panel | step 6 pre-screen |
| `BLAS_THREADS` unset | containers thread to 192 cores; ~40x slower, ~80 GB per job | export it |
| `-e USER` unset | snakemake's `getpass.getuser()` -> `KeyError: getpwuid()` on the first Run | every command above sets it |
| `docker stop --filter ancestor=` | stops other users' containers on this shared host | `--name`, stop individually |
| forgetting `--name` | Ctrl-C kills the client, container keeps the results lock | every command above names it |

## Aborting / resuming

Snakemake's own DAG is the resume mechanism — re-running step 5 or step 9 with the same
`$SEEDS` picks up where it stopped. If a container was killed, clear the lock first:

```bash
$DK run --rm --user "$UIDGID" -e USER=adaptogene -v "$PWD":/pipeline cline-go:latest \
  snakemake -s Snakefile --unlock --config mode=gea --configfile config_MVP<seed>_c1.yaml
```

To abandon a block entirely before step 5: restore
`benchmarks/mvp_seeds.tsv.pre_${MVP_COHORT}` and delete this block's `config_MVP*` files. After
step 5 the results trees exist and abandoning means deleting `MVP{seed}_results/` for the
block's seeds.
