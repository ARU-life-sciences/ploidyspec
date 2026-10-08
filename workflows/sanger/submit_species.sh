#!/usr/bin/env bash
# Submit `ploidyspec all` (run_one.sh) for one species from jobs.tsv, with every
# argument shell-quoted (a chromosome pattern containing `|` must survive bsub).
# Usage: submit_species.sh SPECIES ["extra flags for all"]
# Env: QUEUE (default normal), MEM (MB, default 32000), THREADS (default 8),
#      WAIT_FOR (an LSF job id to wait for).
# Always adds --with-te-markers --with-windowed-homeologs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$HERE"
sp="$1"; extra="${2:-}"
line=$(awk -F'\t' -v s="$sp" '$1==s' workflows/sanger/jobs.tsv)
[ -n "$line" ] || { echo "no jobs.tsv row for $sp" >&2; exit 1; }
IFS=$'\t' read -r _ manifest regex hap <<< "$line"
mkdir -p "results/$sp" logs
mem="${MEM:-32000}"; threads="${THREADS:-8}"
cmd="PLOIDYSPEC_ALL_EXTRA=$(printf '%q' "--with-te-markers --with-windowed-homeologs $extra") \
workflows/sanger/run_one.sh $(printf '%q ' "$sp" "$manifest" "$regex" "$hap")"
dep=()
[ -n "${WAIT_FOR:-}" ] && dep=(-w "done(${WAIT_FOR})")
bsub -q "${QUEUE:-normal}" -n "$threads" -M "$mem" -R "span[hosts=1] select[mem>=$mem] rusage[mem=$mem]" \
    "${dep[@]}" -o "results/$sp/lsf.%J.out" -e "results/$sp/lsf.%J.err" -J "ploidyspec.$sp" \
    "bash -c $(printf '%q' "$cmd")" | tee -a logs/new_species_jobs.txt
