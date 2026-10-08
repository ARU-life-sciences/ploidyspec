#!/usr/bin/env bash
# Add the unrelated-chromosome control tracks for the residual-tetrasomy test
# (2026-10-08) to every species in jobs.tsv that already has windowed-homeolog
# tracks, then re-run the stages that read them. One LSF job per species;
# threads and memory sized as in rerun_windowed.sh.
# Usage: rerun_residual_controls.sh [species ...]   (default: all eligible)
# Queue: $QUEUE (default normal).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$HERE"
mkdir -p logs
: > logs/rerun_residual_controls_jobs.txt

want=" $* "
while IFS=$'\t' read -r sp manifest regex hap; do
    [ $# -eq 0 ] || [[ "$want" == *" $sp "* ]] || continue
    d="results/$sp"
    [ -s "$d/homeologs/windowed_homeologs_all.tsv" ] && [ -f "$d/sequences.tsv" ] || continue
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
./ploidyspec.sh windowed-homeologs $q --controls-only --window-k 15 --window 250000 --threads $threads; \
./ploidyspec.sh rediploidization $q --threads $threads; \
./ploidyspec.sh summary --outdir $d; \
./ploidyspec.sh report --outdir $d"
    bsub -q "${QUEUE:-normal}" -n "$threads" -M "$mem" -R "span[hosts=1] select[mem>=$mem] rusage[mem=$mem]" \
        -o "$d/lsf.controls.%J.out" -e "$d/lsf.controls.%J.err" \
        -J "controls.$sp" "bash -c $(printf '%q' "$job")" | tee -a logs/rerun_residual_controls_jobs.txt
    echo "  $sp: max chrom $((maxlen / 1000000)) Mb, $threads threads, ${mem} MB"
done < workflows/sanger/jobs.tsv
