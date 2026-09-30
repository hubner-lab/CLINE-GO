#!/usr/bin/env bash
# =============================================================================
# mvp_run_sweep.sh -- run the GEA parameter-cell sweep across MVP replicates.
#
#   benchmarks/mvp_run_sweep.sh [SEEDS_CSV|all] [SEED_JOBS] [CELLS]
#
# Defaults: all seeds in benchmarks/mvp_seeds.tsv, 7 seeds at a time, 5 cells.
#
# WHAT IT DOES, per seed, strictly serially:
#   upstream:  processing -> prestructure -> structure      (once; Snakemake skips if current)
#   pregea:                                                 (once; supplies the recommender arm)
#   cells:     mode=gea with config_MVP{seed}_c{i}.yaml, i = 1..CELLS
#
# After each cell the four p-value tables and the Manhattan plots are copied out to
#   benchmarks/mvp_eval/params/MVP{seed}/c{i}/
# because the pipeline writes every cell to the SAME path -- {method}_pvalues_K{k_best}.tsv
# encodes k_best, never the swept parameter (common.smk:1308). Harvesting is what turns a
# sequence of overwrites into a ladder.
#
# Cells are NOT forced. Snakemake's params rerun-trigger detects the changed rung on its
# own (verified: "Params have changed since last execution: before: 5 now: 0"), so only the
# method rules whose parameter actually moved are re-fitted. The HARVEST GATE below is what
# guards against that trigger silently failing: a cell whose harvested table is byte-identical
# to the previous cell's is a hard error, not a warning -- an identical table means the rung
# never took effect and the whole ladder would be a fiction.
#
# Seeds run in parallel; cells within a seed never do. Each {PROJECT}_results/ carries its own
# .snakemake/ (locks, metadata) -- verified: the repo-root .snakemake holds only log/ -- so
# cross-seed parallelism is safe, while two modes on ONE project would collide on its lock.
#
# Docker: /snap/bin/docker cannot start on this host (setuid blocked under nosuid+NoNewPrivs),
# so every invocation goes through a nix-provided client. See CLAUDE.local.md.
# -e USER=pipeline is required or snakemake's getpwuid() lookup aborts under --user.
# =============================================================================
set -uo pipefail

PIPELINE_ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
IMAGE="${IMAGE:-cline-go:latest}"
SEEDS_ARG="${1:-all}"
SEED_JOBS="${2:-6}"
CELLS="${3:-5}"

# BLAS_THREADS is not a tuning nicety, it is the difference between this sweep finishing
# overnight and not finishing at all. The container sees the HOST's nproc (192) regardless of
# --cpus, so OpenBLAS-pthread spawns 192 threads that thrash against the container's CPU
# share. Measured in-image on a 1500x1500 crossprod+svd, --cpus=32:
#     unset (192 threads)  52.65 s      1 thread   3.16 s
#     8 threads             1.42 s     16 threads  1.30 s   <- 40x faster than unset
# Unset, one RDA fit on seed 1232548 ran >27 min and climbed to 81 GB resident without
# finishing; the thread-local buffers are also where the memory went. Keep this pinned and
# IDENTICAL across every run in the sweep -- BLAS thread count can reorder floating-point
# summation, so mixing settings within one benchmark would be a reproducibility hole.
BLAS_THREADS="${BLAS_THREADS:-16}"
CPUS_PER_SEED="${CPUS_PER_SEED:-20}"
MEM_PER_SEED="${MEM_PER_SEED:-60g}"
SNAKE_CORES="${SNAKE_CORES:-4}"

MANIFEST="$PIPELINE_ROOT/benchmarks/mvp_seeds.tsv"
EXCLUSIONS="$PIPELINE_ROOT/benchmarks/mvp_method_exclusions.tsv"
PARAMS_DIR="${PARAMS_DIR:-$PIPELINE_ROOT/benchmarks/mvp_eval/params}"
RUNLOG_DIR="${RUNLOG_DIR:-$PIPELINE_ROOT/benchmarks/mvp_eval/runlogs}"

# SWEEP_METHODS / CFG_SUFFIX / PARAMS_DIR are what let journal 06 run an 11-method ladder
# on one replicate WITHOUT overwriting journal 05's harvest. Point CFG_SUFFIX at the
# generated config variant (mvp_write_sweep_configs.R --suffix) and PARAMS_DIR at a fresh
# tree; the defaults reproduce journal 05 exactly.
#
# HARVEST_WZA additionally copies out {method}_wza_K{k}.tsv. mode=gea emits those by
# default (common.smk:1725), but journal 05 never harvested them, so they were overwritten
# by the next cell like everything else.
read -r -a METHODS <<< "${SWEEP_METHODS:-EMMAX LFMM RDA BLINK}"
CFG_SUFFIX="${CFG_SUFFIX:-}"
HARVEST_WZA="${HARVEST_WZA:-0}"
# PROJ_SUFFIX makes this a PARALLEL ARM: project_name becomes MVP{seed}{PROJ_SUFFIX}, so the
# arm gets its own {PROJECT}_results/ tree and its own Snakemake lock, and can therefore run
# CONCURRENTLY with the base arm. That is safe only because the project names differ -- two
# arms sharing a project_name would collide on the lock and overwrite each other's tables.
PROJ_SUFFIX="${PROJ_SUFFIX:-}"
SKIP_PREGEA="${SKIP_PREGEA:-0}"
# HARVEST_DIAG=1 (added 2026-09-29, SS-Clines re-analysis Phase 2a) additionally copies out,
# per cell, what Phase 3 needs to audit one-fit RDA and to read structure beside the p-values.
# See harvest_diag() below.
HARVEST_DIAG="${HARVEST_DIAG:-0}"

mkdir -p "$PARAMS_DIR" "$RUNLOG_DIR"

DOCKER=(nix shell nixpkgs#docker-client -c docker)

# ---- gea gate (added 2026-09-29, SS-Clines re-analysis Phase 2b). Off unless GEA_SLOTS > 0.
# mode=gea holds ~25-29 GB for ~18 of its ~48 min (the rda() fit, 32k-SNP seeds) while every
# upstream mode stays under ~3 GB, so memory is set by how many seeds are in gea AT ONCE, not by
# how many are in flight. The gate lets upstream run wide and admits a seed to gea only when
#   (a) fewer than GEA_SLOTS gea tickets are held, and
#   (b) at least GEA_STAGGER_S seconds have passed since the last admission -- seeds launched
#       together otherwise reach the rda() fit together and their peaks coincide.
# A ticket is a file GATE_DIR/tickets/<proj>. Any driver pointing at the same GATE_DIR shares the
# count, and a ticket written by hand for a seed another process is running counts the same way.
# A ticket is dropped when no container of its project has been seen running for 180 s, so a
# killed driver or a hand-written ticket never leaks a slot. Live overrides, re-read on every
# check: GATE_DIR/{gea_slots,gea_stagger_s,seed_jobs} -- one integer each.
GEA_SLOTS="${GEA_SLOTS:-0}"
GEA_STAGGER_S="${GEA_STAGGER_S:-60}"
GATE_DIR="${GATE_DIR:-$PARAMS_DIR/.gea_gate}"
DOCKER_BIN=""

gate_val() {   # gate_val <name> <default>: GATE_DIR/<name> overrides the env default, live
    local f="$GATE_DIR/$1"
    if [[ -s "$f" ]]; then tr -d '[:space:]' < "$f"; else echo "$2"; fi
}

gate_on() { (( $(gate_val gea_slots "$GEA_SLOTS") > 0 )); }

gea_acquire() {   # gea_acquire <proj> <log>: block until a gea slot is free and the stagger has passed
    gate_on || return 0
    local proj="$1" log="$2" waited=0
    [[ -n "$DOCKER_BIN" ]] || DOCKER_BIN=$(nix shell nixpkgs#docker-client -c sh -c 'command -v docker')
    mkdir -p "$GATE_DIR/tickets"
    while :; do
        if (
            flock 9
            live=$("$DOCKER_BIN" ps --format '{{.Names}}' 2>/dev/null) || exit 1
            now=$(date +%s); held=0
            for t in "$GATE_DIR"/tickets/*; do
                [[ -e "$t" ]] || continue
                p=$(basename "$t")
                if grep -q "^mvp-sweep-${p}-" <<< "$live"; then
                    touch "$t"
                elif (( now - $(date -r "$t" +%s) > 180 )); then
                    rm -f "$t"; continue
                fi
                held=$((held + 1))
            done
            slots=$(gate_val gea_slots "$GEA_SLOTS")
            stagger=$(gate_val gea_stagger_s "$GEA_STAGGER_S")
            last=0; [[ -e "$GATE_DIR/.last_admit" ]] && last=$(date -r "$GATE_DIR/.last_admit" +%s)
            if (( held < slots && now - last >= stagger )); then
                : > "$GATE_DIR/tickets/$proj"; touch "$GATE_DIR/.last_admit"; exit 0
            fi
            exit 1
        ) 9>"$GATE_DIR/.lock"; then break; fi
        (( waited )) || echo "[$proj] gea gate: waiting for a slot  $(date -Is)" | tee -a "$log"
        waited=1; sleep 20
    done
    echo "[$proj] gea gate: admitted  $(date -Is)" | tee -a "$log"
}

gea_release() { gate_on || return 0; rm -f "$GATE_DIR/tickets/$1"; }

if [[ "$SEEDS_ARG" == "all" ]]; then
    mapfile -t SEEDS < <(awk 'NR>1 && NF {print $1}' "$MANIFEST")
else
    IFS=',' read -r -a SEEDS <<< "$SEEDS_ARG"
fi

kbest_of() { awk -v s="$1" 'NR>1 && $1==s {print $11; exit}' "$MANIFEST"; }

snake() {   # snake <seed> <mode> <configfile> <logfile>
    # --name: killing the docker CLIENT leaves the container running and holding the results
    # lock; a named orphan is one `docker stop` away (CLAUDE.md, "Always pass --name").
    "${DOCKER[@]}" run --name "mvp-sweep-MVP${1}${PROJ_SUFFIX}-$2" \
        --user "$(id -u):$(id -g)" --rm -e USER=pipeline \
        -e OPENBLAS_NUM_THREADS="$BLAS_THREADS" -e OMP_NUM_THREADS="$BLAS_THREADS" \
        --cpus="$CPUS_PER_SEED" --memory="$MEM_PER_SEED" \
        -v "$PIPELINE_ROOT:/pipeline" "$IMAGE" \
        snakemake -c"$SNAKE_CORES" -s Snakefile \
            --config mode="$2" --configfile "$3" --scheduler greedy \
            --rerun-incomplete \
        >> "$4" 2>&1
}

harvest_diag() {   # harvest_diag <proj> <res> <k> <dest> <log>
    # Only RDA writes side tables (rda.R: candidates / diagnostics / anova). ld_decay comes
    # from Structure/tables but is pulled into mode=gea by assoc_wza (auto_genome_wide), so
    # it is current for the cell. eigenvalues / tracywidom are the LD-pruned LEA PCA's; the
    # _work/ path encodes the filter and LD parameters, so it is globbed and must resolve to
    # exactly ONE .pca directory. Every file is FATAL when missing -- the completion gate
    # counts them, and a skipped copy would look like a harvested cell.
    local proj="$1" res="$2" k="$3" dest="$4" log="$5" f src
    if [[ " ${METHODS[*]} " == *" RDA "* ]]; then
        for f in candidates diagnostics anova; do
            src="$res/GEA/tables/methods/RDA/RDA_${f}_K${k}.tsv"
            [[ -s "$src" ]] || { echo "[$proj] FATAL: HARVEST_DIAG missing/empty $src" | tee -a "$log"; return 1; }
            cp "$src" "$dest/RDA_${f}.tsv"
        done
    fi
    src="$res/Structure/tables/ld_decay_half_distances.tsv"
    [[ -s "$src" ]] || { echo "[$proj] FATAL: HARVEST_DIAG missing/empty $src" | tee -a "$log"; return 1; }
    cp "$src" "$dest/ld_decay_half_distances.tsv"
    local pca=( "$res"/_work/*/*/*.pca )
    if (( ${#pca[@]} != 1 )) || [[ ! -d "${pca[0]}" ]]; then
        echo "[$proj] FATAL: HARVEST_DIAG expected exactly one _work/*/*/*.pca dir, found: ${pca[*]}" | tee -a "$log"
        return 1
    fi
    for f in eigenvalues tracywidom; do
        src=( "${pca[0]}"/*."$f" )
        [[ ${#src[@]} -eq 1 && -s "${src[0]}" ]] \
            || { echo "[$proj] FATAL: HARVEST_DIAG missing/empty ${pca[0]}/*.$f" | tee -a "$log"; return 1; }
        cp "${src[0]}" "$dest/pca.$f"
    done
}

run_seed() {
    local seed="$1"
    local proj="MVP${seed}${PROJ_SUFFIX}"
    local res="$PIPELINE_ROOT/${proj}_results"
    local k; k="$(kbest_of "$seed")"
    local log="$RUNLOG_DIR/${proj}.log"
    : > "$log"

    if [[ -z "$k" ]]; then echo "[$proj] FATAL: seed not in manifest" | tee -a "$log"; return 1; fi
    echo "[$proj] k_best=$k  cells=$CELLS  $(date -Is)" | tee -a "$log"

    # Clear a stale lock left by a killed run. Safe here specifically because this driver is
    # the only thing that ever runs snakemake against MVP* projects, and it serialises every
    # mode within a seed -- so a lock present at seed start is always a corpse, never a live
    # writer. Do NOT copy this into a context where a project can have a concurrent runner.
    if [[ -d "$res/.snakemake/locks" ]] && [[ -n "$(ls -A "$res/.snakemake/locks" 2>/dev/null)" ]]; then
        echo "[$proj] clearing stale snakemake lock" | tee -a "$log"
        "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -e USER=pipeline \
            -v "$PIPELINE_ROOT:/pipeline" "$IMAGE" \
            snakemake -s Snakefile --unlock --config mode=gea \
                --configfile "config_${proj}_c1${CFG_SUFFIX}.yaml" >> "$log" 2>&1
    fi

    # ---- upstream. Snakemake decides what is stale; MVP1231288 rebuilds by itself here
    # because the linkage-group-boundary fix changed its input VCF's mtime (MVP_README.md:183).
    # Upstream runs against the c1 CELL config, not the base config. Every cell config is
    # the base config with only the GEA: block and LDdecay.scope changed, and upstream reads
    # neither GEA params nor the per-chromosome LD-decay scope -- but the scope patch is what
    # keeps `ld_decay_run_chr` (which segfaults inside PopLDdecay on monomorphic deme x
    # chromosome subsets) out of the DAG. Using the base config here would reintroduce it.
    for mode in processing prestructure structure; do
        echo "[$proj] mode=$mode  $(date -Is)" | tee -a "$log"
        snake "$seed" "$mode" "config_${proj}_c1${CFG_SUFFIX}.yaml" "$log" \
            || { echo "[$proj] FAILED at mode=$mode" | tee -a "$log"; return 1; }
    done

    # ---- PreGEA: the hyperparameter recommender. Its ladders are fitted on the LD-pruned
    # marker set, where roughly half the causal loci are gone (seed 1232548: 30 -> 14 on the
    # temperature axis), so PreGEA NOMINATES a rung and the production cells below are what
    # actually get scored. Never mix a PreGEA recall number into the production surface.
    #
    # SKIP_PREGEA=1 (added 2026-09-28, SS-Clines re-analysis): PreGEA is 50-70 % of per-seed
    # wall time and nothing downstream of the sweep reads it -- the cells take their rungs from
    # the generated config, never from pregea_recommendations.tsv.
    if [[ "$SKIP_PREGEA" == "1" ]]; then
        echo "[$proj] mode=pregea SKIPPED (SKIP_PREGEA=1)" | tee -a "$log"
    else
        echo "[$proj] mode=pregea  $(date -Is)" | tee -a "$log"
        snake "$seed" pregea "config_${proj}_c1${CFG_SUFFIX}.yaml" "$log" \
            || echo "[$proj] WARNING: pregea failed, recommender arm unavailable" | tee -a "$log"
    fi

    # ---- cells
    local prev_dir=""
    for i in $(seq 1 "$CELLS"); do
        local cfg="config_${proj}_c${i}${CFG_SUFFIX}.yaml"
        local dest="$PARAMS_DIR/${proj}/c${i}"
        [[ -f "$PIPELINE_ROOT/$cfg" ]] || { echo "[$proj] FATAL: missing $cfg" | tee -a "$log"; return 1; }

        if [[ -f "$dest/${METHODS[-1]}_pvalues.tsv" ]]; then
            echo "[$proj] c$i already harvested, skipping" | tee -a "$log"
            prev_dir="$dest"; continue
        fi

        gea_acquire "$proj" "$log"
        echo "[$proj] c$i  $(date -Is)" | tee -a "$log"
        local stamp="$res/.cell_start"
        touch "$stamp"
        snake "$seed" gea "$cfg" "$log"; local rc=$?
        gea_release "$proj"
        (( rc == 0 )) || { echo "[$proj] FAILED at c$i" | tee -a "$log"; return 1; }
        echo "[$proj] c$i gea finished  $(date -Is)" | tee -a "$log"

        mkdir -p "$dest/manhattan"
        # Before the p-value tables, so the skip sentinel above ({last method}_pvalues.tsv)
        # still means "cell complete" when a diag file is missing.
        if [[ "$HARVEST_DIAG" == "1" ]]; then
            harvest_diag "$proj" "$res" "$k" "$dest" "$log" || return 1
        fi
        for m in "${METHODS[@]}"; do
            # [changed 2026-08-04] The exclusion test used to live INSIDE the `! -s "$src"`
            # branch, i.e. it was only consulted when the table was missing. But a method
            # declared in mvp_method_exclusions.tsv is dropped from GEA.configs, so mode=gea
            # never regenerates it -- while the table from an EARLIER run is still sitting in
            # {PROJECT}_results/GEA/tables/methods/. The old order therefore found that stale
            # file, copied it into every cell of the ladder, and the byte-identical harvest
            # gate below then (correctly) aborted the seed: c2's table was c1's table.
            # Checking the exclusion FIRST means an excluded method is never harvested at all.
            # Seen on MVP1232568 / MVP1231578 (SUPER) after those seeds were re-run with the
            # method excluded.
            if awk -F'\t' -v s="$seed" -v m="$m" 'NR>1 && $1==s && $2==m {found=1}
                                                  END {exit !found}' "$EXCLUSIONS" 2>/dev/null; then
                echo "[$proj] c$i $m EXCLUDED (declared in mvp_method_exclusions.tsv)" | tee -a "$log"
                continue
            fi
            src="$res/GEA/tables/methods/$m/${m}_pvalues_K${k}.tsv"
            if [[ ! -s "$src" ]]; then
                # A missing table for a NON-excluded method is FATAL: a genuine silent failure
                # must not look like a known gap, which is the class of bug this harness exists
                # to rule out.
                echo "[$proj] FATAL: c$i missing/empty $src" | tee -a "$log"; return 1
            fi
            cp "$src" "$dest/${m}_pvalues.tsv"
            if [[ "$HARVEST_WZA" == "1" ]]; then
                wsrc="$res/GEA/tables/methods/$m/${m}_wza_K${k}.tsv"
                # Non-fatal: supports_wza is a per-method registry flag, so a missing WZA
                # table is a legitimate opt-out, not the silent-failure case the p-value
                # gate above exists to catch.
                [[ -s "$wsrc" ]] && cp "$wsrc" "$dest/${m}_wza.tsv" \
                    || echo "[$proj] c$i $m: no WZA table (supports_wza off?)" >> "$log"
            fi
        done
        cp -r "$res/GEA/plots/manhattan/." "$dest/manhattan/" 2>/dev/null
        cp "$res/GEA/tables/selected_snps.tsv" "$dest/" 2>/dev/null

        # The Manhattans are what a user reads an operating point off, so a stale figure
        # paired with a fresh p-value table would be a silently wrong plot. Any PNG older
        # than this cell's start means Snakemake decided it did not need regenerating.
        local stale
        stale=$(find "$dest/manhattan" -name '*.png' ! -newer "$stamp" | wc -l)
        if (( stale > 0 )); then
            echo "[$proj] FATAL: c$i harvested $stale Manhattan PNG(s) older than the cell" \
                 "start — figures do not match this cell's p-values." | tee -a "$log"
            return 1
        fi

        # ---- HARVEST GATE (see header). An unchanged table means the rung did not apply.
        if [[ -n "$prev_dir" ]]; then
            for m in "${METHODS[@]}"; do
                [[ -s "$dest/${m}_pvalues.tsv" && -s "$prev_dir/${m}_pvalues.tsv" ]] || continue
                if cmp -s "$dest/${m}_pvalues.tsv" "$prev_dir/${m}_pvalues.tsv"; then
                    echo "[$proj] FATAL: c$i $m is byte-identical to the previous cell — the" \
                         "parameter did not take effect. Ladder would be a fiction; aborting." | tee -a "$log"
                    return 1
                fi
            done
        fi
        prev_dir="$dest"
    done

    echo "[$proj] DONE $(date -Is)" | tee -a "$log"
}

# ------------------------------------------------------------- bounded worker pool
echo "INFO: ${#SEEDS[@]} seeds, $SEED_JOBS at a time, $CELLS cells each"
echo "INFO: per seed --cpus=$CPUS_PER_SEED --memory=$MEM_PER_SEED, snakemake -c$SNAKE_CORES"
fail=0
for seed in "${SEEDS[@]}"; do
    # seed_jobs is re-read every 15 s (gate_val), so concurrency can be raised or lowered live.
    while (( $(jobs -rp | wc -l) >= $(gate_val seed_jobs "$SEED_JOBS") )); do sleep 15; done
    run_seed "$seed" &
done
wait
for seed in "${SEEDS[@]}"; do
    # [changed 2026-08-03] was: MVP${seed}.log, with no PROJ_SUFFIX. run_seed writes to
    # MVP{seed}{PROJ_SUFFIX}.log, so on a parallel arm this checked the BASE arm's log --
    # which, once the base arm had run, always contained DONE. The arm could fail on every
    # seed and this loop would still report "all seeds complete".
    proj="MVP${seed}${PROJ_SUFFIX}"
    grep -q "^\[${proj}\] DONE" "$RUNLOG_DIR/${proj}.log" 2>/dev/null \
        || { echo "INCOMPLETE: ${proj} (see $RUNLOG_DIR/${proj}.log)"; fail=1; }
done
echo "INFO: sweep finished, $( [[ $fail -eq 0 ]] && echo 'all seeds complete' || echo 'SOME SEEDS INCOMPLETE' )"
exit $fail
