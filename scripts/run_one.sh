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
mkdir -p "$outdir"

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

echo "=== ${species_id}: all ==="
./ploidyspec.sh all \
    --manifest "$manifest" \
    --outdir "$outdir" \
    --k 11,13,15,17,19,23 --window-k 15 --min-len 1000000 \
    --window 250000 \
    --threads "${LSB_DJOB_NUMPROC:-8}" \
    "${extra_args[@]}"

echo "=== ${species_id}: homeologs (bonus, non-fatal) ==="
./ploidyspec.sh homeologs --manifest "$manifest" --outdir "$outdir" --k 11,13,15,17,19,23 "${extra_args[@]}" || true
./ploidyspec.sh windowed-homeologs --manifest "$manifest" --outdir "$outdir" \
    --window-k 15 --window 250000 --threads "${LSB_DJOB_NUMPROC:-8}" "${extra_args[@]}" || true

echo "=== ${species_id}: done ==="
