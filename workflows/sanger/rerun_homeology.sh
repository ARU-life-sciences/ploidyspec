#!/usr/bin/env bash
# Re-run windowed-homeologs with the genome-wide homeology map (2026-10-09),
# then the stages that read it (rediploidization, summary, report). One LSF
# job per jobs.tsv species with k-mer tables; threads and memory sized from
# the longest chromosome as in rerun_residual_controls.sh.
# Usage: rerun_homeology.sh [species ...]   (default: all eligible)
# Env: QUEUE (default normal); MAP_ONLY=1 computes only the map.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$HERE"
mkdir -p logs
log=logs/rerun_homeology_jobs.txt
: > "$log"

active=" $(bjobs -noheader -o job_name 2>/dev/null | sed 's/^[^.]*\.//' | tr '\n' ' ') "
want=" $* "
while IFS=$'\t' read -r sp manifest regex hap; do
    [ $# -eq 0 ] || [[ "$want" == *" $sp "* ]] || continue
    d="results/$sp"
    [ -f "$d/sequences.tsv" ] && [ -f "$d/matrix/whole_chrom_distance_matrix.csv" ] || continue
    if [[ "$active" == *" $sp "* ]]; then echo "  $sp: job active, skipped"; continue; fi
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
    if [ "${MAP_ONLY:-0}" = 1 ]; then
        stages="./ploidyspec.sh windowed-homeologs $q --map-only --window 250000 --threads $threads"
    else
        stages="./ploidyspec.sh windowed-homeologs $q --window-k 15 --window 250000 --threads $threads; \
./ploidyspec.sh rediploidization $q --threads $threads; \
./ploidyspec.sh summary --outdir $d; \
./ploidyspec.sh report --outdir $d"
    fi
    job="source /etc/profile.d/modules.sh; module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0; \
export MPLCONFIGDIR=\${TMPDIR:-/tmp}; set -e; $stages"
    bsub -q "${QUEUE:-normal}" -n "$threads" -M "$mem" -R "span[hosts=1] select[mem>=$mem] rusage[mem=$mem]" \
        -o "$d/lsf.homeology.%J.out" -e "$d/lsf.homeology.%J.err" \
        -J "homeology.$sp" "bash -c $(printf '%q' "$job")" | tee -a "$log"
    echo "  $sp: max chrom $((maxlen / 1000000)) Mb, $threads threads, ${mem} MB"
done < workflows/sanger/jobs.tsv
