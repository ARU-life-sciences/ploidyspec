#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

species_id="$1"
manifest="$2"
regex_field="${3:-}"
outdir="results/${species_id}"

extra_args=()
if [ -n "$regex_field" ]; then
    IFS='|' read -ra patterns <<< "$regex_field"
    for pat in "${patterns[@]}"; do
        extra_args+=(--chrom-regex "$pat")
    done
fi

# Multi-k Mash-correction verification sweep (see mash.py / plan doc): only
# kmers/matrix/homeologs are affected by this change, so windowed is skipped
# here on purpose -- it's unmodified and re-running it would just burn cluster
# time re-doing expensive per-window FastK calls for no new information.
K_SWEEP="11,13,15,17,19,23"
THREADS="${LSB_DJOB_NUMPROC:-8}"

# sequences.tsv on disk for these 3 species predates this environment (its
# `source` column has stale machine-local paths from wherever `prepare` was
# last run) -- regenerate it from the current manifest first. Cheap, no FastK
# involved.
echo "=== ${species_id}: prepare (regenerate sequences.tsv from current manifest) ==="
./ploidyspec.sh prepare --manifest "$manifest" --outdir "$outdir" --min-len 1000000 "${extra_args[@]}"

echo "=== ${species_id}: kmers (k sweep ${K_SWEEP}) ==="
./ploidyspec.sh kmers --manifest "$manifest" --outdir "$outdir" --k "$K_SWEEP" --threads "$THREADS" "${extra_args[@]}"

echo "=== ${species_id}: matrix (Mash-corrected, k sweep ${K_SWEEP}) ==="
./ploidyspec.sh matrix --manifest "$manifest" --outdir "$outdir" --k "$K_SWEEP" --threads "$THREADS" "${extra_args[@]}"

echo "=== ${species_id}: homeologs (resolution_limited flag) ==="
./ploidyspec.sh homeologs --manifest "$manifest" --outdir "$outdir" --k "$K_SWEEP" "${extra_args[@]}"

echo "=== ${species_id}: verification done ==="
