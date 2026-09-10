#!/bin/sh
# emmax_run.sh — preflight-guarded launcher for the checked-in EMMAX binaries.
#
# Usage:  scripts/emmax_run.sh /pipeline/scripts/emmax-kin-intel64 -v -d 10 -x <prefix>
#
# Checks that the binary can actually run natively on this host, then exec()s
# it with the arguments given. Every EMMAX call site in the pipeline goes
# through this wrapper so the check cannot be forgotten at a new one.
#
# WHY THIS EXISTS
# ---------------
# scripts/emmax-intel64 and scripts/emmax-kin-intel64 are pre-built ELF
# x86-64 binaries statically linked against Intel MKL and the Intel OpenMP
# runtime (`strings` shows the KMP topology/affinity strings). Under an
# x86-64-on-arm64 emulation layer (Docker Desktop / colima on Apple silicon,
# where the amd64-only image is run through QEMU or Rosetta) the KMP CPU
# topology probe never settles: emmax-kin prints "Identified N individuals",
# then "nex = 0", and then never exits and never crashes. Snakemake waits
# forever, the user kills the docker CLIENT, and the CONTAINER keeps running
# and keeps holding the results-directory lock.
#
# Architecture is a property of the IMAGE, not of the binary, so shipping a
# second, arm64-built EMMAX cannot help: inside an amd64 container `uname -m`
# is always x86_64. The only real fix is a multi-arch image. Until then, the
# correct behaviour is to fail immediately and say why (pipeline philosophy:
# "Fail loudly, not silently") instead of hanging.
#
# Escape hatch: ADAPTOGENE_ALLOW_EMULATED_EMMAX=1 downgrades the emulation
# error to a warning and pins the Intel OpenMP/MKL runtime to a single thread,
# which is the usual cure for KMP-under-emulation spins. Untested here; it is
# an experiment switch, not a supported configuration. It is deliberately NOT
# the default, because forcing a different thread count changes the order of
# the floating-point reductions inside MKL and could perturb the kinship
# matrix on hosts where EMMAX works today.

set -eu

BIN="${1:?usage: emmax_run.sh <path-to-emmax-binary> [args...]}"
shift

[ -x "$BIN" ] || {
    echo "ERROR: EMMAX binary not found or not executable: $BIN" >&2
    exit 1
}

# ELF e_machine, bytes 18-19, little-endian. 62 = EM_X86_64, 183 = EM_AARCH64.
elf_machine=$(od -An -tu1 -j18 -N2 "$BIN" | awk '{print $1 + 256 * $2}')
host_arch=$(uname -m)

case "$elf_machine" in
    62)  bin_arch=x86_64 ;;
    183) bin_arch=aarch64 ;;
    *)   bin_arch="ELF-machine-$elf_machine" ;;
esac

# 1. Hard architecture mismatch: the kernel is not the binary's architecture at
#    all (e.g. an arm64 image that somehow got these binaries). Nothing to try.
if [ "$bin_arch" != "$host_arch" ]; then
    echo "ERROR: $BIN is a $bin_arch binary but this host reports '$host_arch'." >&2
    echo "       EMMAX cannot run here. The image is amd64-only (rocker/shiny:4.5" >&2
    echo "       publishes no arm64 manifest); run EMMAX-based modes on a Linux" >&2
    echo "       x86-64 host, or drop 'emmax' from GEA.configs / GWAS.configs." >&2
    exit 1
fi

# 2. Emulated x86-64: `uname -m` says x86_64, but the kernel underneath is not
#    an x86 kernel, so /proc/cpuinfo carries no x86 'flags'/'vendor_id' lines.
#    This is exactly what an amd64 container sees on an arm64 colima VM — the
#    binary would start, print "nex = 0", and hang forever.
emulated=no
if [ "$host_arch" = "x86_64" ] && [ -r /proc/cpuinfo ]; then
    if ! grep -qE '^(flags|vendor_id)[[:space:]]*:' /proc/cpuinfo; then
        emulated=yes
    fi
fi

if [ "$emulated" = yes ]; then
    if [ "${ADAPTOGENE_ALLOW_EMULATED_EMMAX:-0}" = "1" ]; then
        echo "WARNING: emulated x86-64 detected; ADAPTOGENE_ALLOW_EMULATED_EMMAX=1" >&2
        echo "         is set, so EMMAX will be attempted with the Intel OpenMP/MKL" >&2
        echo "         runtime pinned to one thread. It may still hang." >&2
        OMP_NUM_THREADS=1; MKL_NUM_THREADS=1; KMP_AFFINITY=disabled
        MKL_THREADING_LAYER=SEQUENTIAL
        export OMP_NUM_THREADS MKL_NUM_THREADS KMP_AFFINITY MKL_THREADING_LAYER
    else
        echo "ERROR: this container reports x86_64 but is running on a non-x86 kernel" >&2
        echo "       (/proc/cpuinfo has no 'flags'/'vendor_id' line) — i.e. the amd64" >&2
        echo "       image is being emulated on an arm64 host." >&2
        echo "       EMMAX's Intel OpenMP/MKL runtime HANGS under emulation: it prints" >&2
        echo "       'Identified N individuals' and 'nex = 0' and then never exits, so" >&2
        echo "       Snakemake would wait forever and the results-directory lock would" >&2
        echo "       be held by an orphaned container." >&2
        echo "       Run EMMAX-based modes (gea/gwas/pregea) on a Linux x86-64 host," >&2
        echo "       or drop 'emmax' from GEA.configs / GWAS.configs." >&2
        echo "       To experiment anyway: ADAPTOGENE_ALLOW_EMULATED_EMMAX=1" >&2
        exit 1
    fi
fi

exec "$BIN" "$@"
