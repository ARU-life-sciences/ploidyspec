#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

STAGE_JOB_ID="${1:-}"  # optional LSF job ID to wait on (the curated-fasta re-bgzip job)

mkdir -p logs
while IFS=$'\t' read -r species_id manifest regex_field hap_regex; do
    outdir="results/${species_id}"
    mkdir -p "$outdir"
    dep_args=()
    if [ -n "$STAGE_JOB_ID" ]; then
        dep_args=(-w "done(${STAGE_JOB_ID})")
    fi
    # bsub flattens multiple job-command arguments into one space-joined string
    # that gets re-parsed by a shell on the exec host -- any shell metacharacter
    # inside an argument (e.g. the literal `|` separating chrom-regex patterns)
    # becomes live syntax again unless pre-escaped with %q here.
    job_cmd=$(printf '%q ' scripts/run_one.sh "$species_id" "$manifest" "$regex_field" "$hap_regex")
    bsub -q long -n 8 -M 16000 -R "span[hosts=1] select[mem>=16000] rusage[mem=16000]" \
        "${dep_args[@]}" \
        -o "${outdir}/lsf.%J.out" -e "${outdir}/lsf.%J.err" \
        -J "ploidyspec.${species_id}" \
        "$job_cmd"
done < scripts/jobs.tsv
