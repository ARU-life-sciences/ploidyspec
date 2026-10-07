#!/usr/bin/env bash
# Progress of the panel-wide rediploidization run (ROADMAP.md, Progress log).
# Usage: scripts/rediploidization_status.sh   -- one line per species + totals
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
printf "%-14s %-6s %s\n" species state summary
while read -r sp job; do
    if bjobs "$job" 2>/dev/null | grep -qE 'PEND|RUN'; then
        state=$(bjobs -noheader -o stat "$job")
    elif grep -qs "Successfully completed" "logs/redip_$sp.$job.out"; then
        state=DONE
    elif [ -s "logs/redip_$sp.$job.out" ]; then
        state=FAIL   # see logs/redip_$sp.$job.err
    else
        state=UNKN   # no LSF output yet
    fi
    summary=""
    if [ "$state" = DONE ]; then
        summary=$(awk -F'\t' '$1 ~ /copy_state:|n_distinct_fusions/ && $2 != 0 {sub("copy_state:", "", $1); printf "%s=%s ", $1, $2}' \
            "results/$sp/rediploidization/rediploidization_summary.tsv")
    fi
    printf "%-14s %-6s %s\n" "$sp" "$state" "$summary"
done < logs/redip_panel_jobs.txt | tee /dev/stderr | awk '{n[$2]++} END {for (s in n) printf "%s=%d ", s, n[s]; print ""}' >&2
