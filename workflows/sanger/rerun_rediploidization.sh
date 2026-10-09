#!/usr/bin/env bash
# Re-run the cheap stages that read existing windowed tracks (rediploidization,
# summary, report) after a change to how they are read -- e.g. the residual-
# tetrasomy recalibration (2026-10-08). One small LSF job per jobs.tsv species
# that has windowed-homeolog tracks; species with a job of that name already
# queued or running are skipped (they pick up the current code themselves).
# Usage: rerun_rediploidization.sh [species ...]   (default: all eligible)
# Queue: $QUEUE (default normal).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$HERE"
mkdir -p logs
: > logs/rerun_rediploidization_jobs.txt
active=" $(bjobs -noheader -o job_name 2>/dev/null | sed 's/^[^.]*\.//' | tr '\n' ' ') "

want=" $* "
while IFS=$'\t' read -r sp manifest regex hap; do
    [ $# -eq 0 ] || [[ "$want" == *" $sp "* ]] || continue
    d="results/$sp"
    [ -s "$d/homeologs/windowed_homeologs_all.tsv" ] || continue
    if [[ "$active" == *" $sp "* ]]; then echo "  $sp: job active, skipped"; continue; fi
    args=(--manifest "$manifest" --outdir "$d")
    if [ -n "$regex" ]; then
        IFS='|' read -ra pats <<< "$regex"
        for p in "${pats[@]}"; do args+=(--chrom-regex "$p"); done
    fi
    [ -n "$hap" ] && args+=(--hap-regex "$hap")
    q=$(printf '%q ' "${args[@]}")
    job="source /etc/profile.d/modules.sh; module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0; \
export MPLCONFIGDIR=\${TMPDIR:-/tmp}; set -e; \
./ploidyspec.sh rediploidization $q --threads 2; \
./ploidyspec.sh summary --outdir $d; \
./ploidyspec.sh report --outdir $d"
    bsub -q "${QUEUE:-normal}" -n 2 -M 8000 -R "span[hosts=1] select[mem>=8000] rusage[mem=8000]" \
        -o "$d/lsf.rediploidization.%J.out" -e "$d/lsf.rediploidization.%J.err" \
        -J "rediploidization.$sp" "bash -c $(printf '%q' "$job")" | tee -a logs/rerun_rediploidization_jobs.txt
done < workflows/sanger/jobs.tsv
