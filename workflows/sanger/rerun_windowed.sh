#!/usr/bin/env bash
# Re-run the position-free windowed stage (2026-10-07) and everything that reads
# it, for every species in jobs.tsv that already has a windowed/ directory.
# One LSF job per species; threads and memory sized from the longest
# chromosome (~14 bytes/bp per concurrent FastK profile, see windowed.py).
# Usage: rerun_windowed.sh [species ...]   (default: all eligible)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$HERE"
mkdir -p logs
: > logs/rerun_windowed_jobs.txt

want=" $* "
while IFS=$'\t' read -r sp manifest regex hap; do
    [ $# -eq 0 ] || [[ "$want" == *" $sp "* ]] || continue
    d="results/$sp"
    [ -d "$d/windowed" ] && [ -f "$d/sequences.tsv" ] || continue
    maxlen=$(awk -F'\t' 'NR>1 && $5>m {m=$5} END {print m+0}' "$d/sequences.tsv")
    threads=8
    while [ $threads -gt 2 ] && [ $(( threads * 14 * maxlen / 1000000000 )) -gt 48 ]; do
        threads=$(( threads / 2 ))
    done
    mem=$(( threads * 14 * maxlen / 1000000 + 4000 ))
    [ $mem -lt 8000 ] && mem=8000
    args=(--manifest "$manifest" --outdir "$d")
    if [ -n "$regex" ]; then
        IFS='|' read -ra pats <<< "$regex"
        for p in "${pats[@]}"; do args+=(--chrom-regex "$p"); done
    fi
    [ -n "$hap" ] && args+=(--hap-regex "$hap")
    q=$(printf '%q ' "${args[@]}")
    job="source /etc/profile.d/modules.sh; module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0; \
export MPLCONFIGDIR=\${TMPDIR:-/tmp}; set -e; \
./ploidyspec.sh windowed $q --window-k 15 --window 250000 --threads $threads; \
./ploidyspec.sh windowed-homeologs $q --window-k 15 --window 250000 --threads $threads || true; \
./ploidyspec.sh structure --outdir $d; \
./ploidyspec.sh rediploidization $q --threads $threads; \
./ploidyspec.sh report --outdir $d"
    dep=()
    [ -n "${WAIT_FOR:-}" ] && [[ " ${WAIT_FOR_SPECIES:-} " == *" $sp "* ]] && dep=(-w "done($WAIT_FOR)")
    bsub -q long -n "$threads" -M "$mem" -R "span[hosts=1] select[mem>=$mem] rusage[mem=$mem]" \
        "${dep[@]}" -o "$d/lsf.windowed.%J.out" -e "$d/lsf.windowed.%J.err" \
        -J "windowed.$sp" "bash -c $(printf '%q' "$job")" | tee -a logs/rerun_windowed_jobs.txt
    echo "  $sp: max chrom $((maxlen / 1000000)) Mb, $threads threads, ${mem} MB"
done < workflows/sanger/jobs.tsv
