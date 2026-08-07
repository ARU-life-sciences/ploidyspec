#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

species_id="$1"
manifest="$2"
regex_field="${3:-}"
hap_regex="${4:-}"
outdir="results/${species_id}"

extra_args=()
if [ -n "$regex_field" ]; then
    IFS='|' read -ra patterns <<< "$regex_field"
    for pat in "${patterns[@]}"; do
        extra_args+=(--chrom-regex "$pat")
    done
fi
if [ -n "$hap_regex" ]; then
    extra_args+=(--hap-regex "$hap_regex")
fi

THREADS="${LSB_DJOB_NUMPROC:-8}"

echo "=== ${species_id}: prepare ==="
./ploidyspec.sh prepare --manifest "$manifest" --outdir "$outdir" --min-len 1000000 "${extra_args[@]}"

echo "=== ${species_id}: kmers (marker-k=13 only) ==="
./ploidyspec.sh kmers --manifest "$manifest" --outdir "$outdir" --k 13 --threads "$THREADS" "${extra_args[@]}"

echo "=== ${species_id}: te-markers ==="
./ploidyspec.sh te-markers --manifest "$manifest" --outdir "$outdir" --threads "$THREADS" "${extra_args[@]}"

echo "=== ${species_id}: te-markers-windowed ==="
./ploidyspec.sh te-markers-windowed --manifest "$manifest" --outdir "$outdir" --threads "$THREADS" "${extra_args[@]}"

echo "=== ${species_id}: verify_te_markers done ==="
