#!/usr/bin/env bash
# mvp_overnight.sh -- run the journal-09 garden sweep to completion unattended, then score it.
#
# Sequence:
#   1. wait until every base project has its `all`-panel Gradient Forest model (clone lanes
#      hardlink it; without it each lane refits a 1.6 GB forest for ~50 min)
#   2. run clone lanes until all 14 x 112 gardens are harvested, restarting the driver if it
#      exits early (a lane that dies takes its partition with it)
#   3. gap-fill with the serial driver, which walks every garden and skips the complete ones
#   4. score with the Lind & Lotterhos metric and build the comparison tables
#
# Host-RAM guard: this machine is shared. Above HIGH_GB we stop our own p4 then p3 lanes;
# below LOW_GB with headroom we let them back in on the next driver restart. The driver's own
# watchdog still stops a container above 900 GB as a last resort.
set -uo pipefail

ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/CLINE-GO}"
OUT="$ROOT/benchmarks/mvp_eval/offset09"
LOG="$OUT/runlogs/overnight.log"
MANIFEST="$ROOT/benchmarks/mvp_seeds.tsv"
TOTAL=588   # 14 seeds x 42 gardens (30 landscape + all 12 novelty)
HIGH_GB="${HIGH_GB:-820}"
DOCKER=(nix shell nixpkgs#docker-client -c docker)

say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
gardens() { ls -d "$OUT"/gardens/*/*/.complete 2>/dev/null | wc -l; }
models()  { local n=0; for s in $(awk 'NR>1{print $1}' "$MANIFEST"); do
              [[ -s "$ROOT/MVP${s}_results/_intermediate/gradient_forest/all_nospatial/adaptive_model.qs" ]] && n=$((n+1))
            done; echo "$n"; }
host_gb() { free -g | awk '/^Mem:/{print $3}'; }

mkdir -p "$OUT/runlogs"
say "overnight run starting; gardens=$(gardens)/$TOTAL models=$(models)/14"

# ---- 1. wait for the base models -------------------------------------------
while [[ "$(models)" -lt 14 ]]; do sleep 120; done
say "all 14 base all-panel models present"

# ---- 2. clone lanes until every garden is harvested -------------------------
attempt=0
while [[ "$(gardens)" -lt "$TOTAL" ]]; do
    attempt=$((attempt + 1))
    say "lane pass #$attempt; gardens=$(gardens)/$TOTAL host=$(host_gb)G"
    MEM_PER_LANE=24g BLAS_THREADS=2 "$ROOT/benchmarks/mvp_garden_run_par.sh" all 3 4 4 \
        >> "$OUT/runlogs/driver_overnight_${attempt}.log" 2>&1 &
    driver=$!

    # supervise: trim our own lanes if the host climbs
    while kill -0 "$driver" 2>/dev/null; do
        u=$(host_gb)
        if [[ -n "$u" ]] && (( u > HIGH_GB )); then
            for tag in p4 p3; do
                names=$("${DOCKER[@]}" ps --filter "name=gp09" --format '{{.Names}}' | grep "$tag" | tr '\n' ' ')
                if [[ -n "$names" ]]; then
                    say "host ${u}G > ${HIGH_GB}G -- stopping $tag lanes"
                    "${DOCKER[@]}" kill $names >/dev/null 2>&1
                    break
                fi
            done
        fi
        sleep 120
    done
    wait "$driver" 2>/dev/null
    say "lane pass #$attempt ended; gardens=$(gardens)/$TOTAL"
    [[ "$(gardens)" -ge "$TOTAL" ]] && break
    sleep 30
done
say "all gardens harvested: $(gardens)/$TOTAL"

# ---- 3. gap-fill (serial driver walks every garden, skipping complete ones) --
if [[ "$(gardens)" -lt "$TOTAL" ]]; then
    say "gap-fill pass"
    SNAKE_CORES=6 BLAS_THREADS=2 MEM_PER_SEED=40g "$ROOT/benchmarks/mvp_garden_run.sh" all 14 8 all \
        >> "$OUT/runlogs/driver_gapfill.log" 2>&1
    say "gap-fill done: $(gardens)/$TOTAL"
fi

# ---- 4. score ---------------------------------------------------------------
say "scoring with the Lind & Lotterhos metric"
"${DOCKER[@]}" run --rm --user "$(id -u):$(id -g)" -e USER=cline-go \
    -e OPENBLAS_NUM_THREADS=8 --cpus=16 --memory=96g \
    -v "$ROOT:/pipeline" cline-go:latest \
    Rscript /pipeline/benchmarks/eval_offset_lind.R >> "$OUT/runlogs/score.log" 2>&1
say "eval_offset_lind.R exit=$?"

"${DOCKER[@]}" run --rm --user "$(id -u):$(id -g)" -e USER=cline-go \
    -e OPENBLAS_NUM_THREADS=8 --cpus=16 --memory=96g \
    -v "$ROOT:/pipeline" cline-go:latest \
    Rscript /pipeline/benchmarks/mvp_lind_compare.R >> "$OUT/runlogs/score.log" 2>&1
say "mvp_lind_compare.R exit=$?"

say "OVERNIGHT RUN COMPLETE  gardens=$(gardens)/$TOTAL"
