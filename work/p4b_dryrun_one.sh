#!/usr/bin/env bash
# p4b_dryrun_one.sh SEED -- dry-run one seed's sweep config; print "seed rc unexpected_rules"
s="$1"; cd /mnt/data/eugene/ADAPTOGENE
L=benchmarks/mvp_eval/offset13/dryrun/MVP$s.log
nix shell nixpkgs#docker-client -c docker run --rm --name mvp-dry-$s --user "$(id -u):$(id -g)" -e USER=adaptogene \
  -w /pipeline -v "$PWD":/pipeline cline-go:latest snakemake -n -s Snakefile --config mode=maladaptation \
  --configfile config_MVP${s}_sweep.yaml --scheduler greedy > "$L" 2>&1
rc=$?
bad=$(awk '/^Job stats/{f=1;next} f&&/^total/{exit} f&&NF==2&&$1!~/^-/&&$1!="job"{print $1}' "$L" | grep -vxE 'gradient_forest_adaptive|stage_custom_climate_future|geometric_offset|gradient_forest_offset|rda_offset|write_summary|all' | paste -sd,)
echo -e "$s\t$rc\t${bad:-none}"
