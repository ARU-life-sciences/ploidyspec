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

echo "=== ${species_id}: te-markers ==="
./ploidyspec.sh te-markers --manifest "$manifest" --outdir "$outdir" \
    --marker-k 13 --min-count 100 --min-ratio 2.0 --min-len 1000000 \
    --threads "${LSB_DJOB_NUMPROC:-8}" "${extra_args[@]}"

echo "=== ${species_id}: subgenome-report ==="
./ploidyspec.sh subgenome-report --outdir "$outdir"

echo "=== ${species_id}: report (refresh) ==="
./ploidyspec.sh report --outdir "$outdir"

echo "=== ${species_id}: done ==="
