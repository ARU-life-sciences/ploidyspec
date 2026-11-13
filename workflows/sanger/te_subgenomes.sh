#!/usr/bin/env bash
# Prototype run of scripts/te_subgenomes.py on a finished species run.
# usage: te_subgenomes.sh <species_id> [extra args for te_subgenomes.py]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$HERE"
sp="$1"; shift
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python3 scripts/te_subgenomes.py "results/$sp" --procs "${LSB_DJOB_NUMPROC:-4}" "$@"
