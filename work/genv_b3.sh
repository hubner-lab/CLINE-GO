#!/usr/bin/env bash
set -u
cd /mnt/data/eugene/ADAPTOGENE
SEEDS=$(cat /tmp/seeds_b3.txt)
run_one() {
  local s=$1
  nix shell nixpkgs#docker-client -c docker run --rm --name mvp-genv-$s \
    --user "$(id -u):$(id -g)" \
    -e USER=adaptogene -e PIPELINE_ROOT=/pipeline -e OPENBLAS_NUM_THREADS=2 \
    --cpus=2 --memory=8g -v /mnt/data/eugene/ADAPTOGENE:/pipeline cline-go:latest \
    Rscript /pipeline/benchmarks/mvp_write_garden_env.R --seed=$s \
      --gardens=/pipeline/benchmarks/mvp_eval/offset12/gardens_${s}.tsv \
      > work/genv_b3_${s}.log 2>&1 || echo "FAILED $s"
}
n=0
for s in ${SEEDS//,/ }; do
  run_one "$s" &
  n=$((n+1)); (( n % 12 == 0 )) && wait
done
wait
echo "env tables done $(date -Is)"
