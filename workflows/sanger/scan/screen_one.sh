#!/usr/bin/env bash
# k-mer screen of one ToL release: prepare, kmers and matrix at k = 23 only,
# homeologs and structure; then the k-mer tables and per-chromosome FASTAs are
# deleted (cleanup). Writes <outdir>/screen.done or <outdir>/screen.failed,
# and <outdir>/screen.mode: how chromosomes were read (names, auto, hap1, none).
# Usage: screen_one.sh OUTDIR THREADS  (OUTDIR holds manifest.tsv)
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$HERE"
d=$1; threads=$2
source /etc/profile.d/modules.sh
module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0
export MPLCONFIGDIR=${TMPDIR:-/tmp}
rm -f "$d/screen.done" "$d/screen.failed" "$d/screen.mode"
n_files=$(grep -vc '^#' "$d/manifest.tsv")
haps_read() { awk -F'\t' 'NR>1{print $2}' "$d/sequences.tsv" 2>/dev/null | sort -u | wc -l; }
# prepare: header names; if a haplotype comes out empty (its scaffolds carry no
# chromosome names, xgLitLitt3 hap2), number by matching (--chrom-naming auto);
# if that fails too (a contig-level haplotype, daBudDavi1 hap2), hap1 alone
mode=""
for try in names auto hap1; do
    rm -f "$d/sequences.tsv" "$d/unplaced.tsv"
    m="$d/manifest.tsv"; naming=(--chrom-naming detect)
    [ $try = auto ] && naming=(--chrom-naming auto)
    if [ $try = hap1 ]; then
        { head -1 "$d/manifest.tsv"; grep -v '^#' "$d/manifest.tsv" | head -1; } > "$d/manifest_hap1.tsv"
        m="$d/manifest_hap1.tsv"; n_files=1
    fi
    if ./ploidyspec.sh prepare --manifest "$m" --outdir "$d" "${naming[@]}" \
        && [ "$(haps_read)" -eq "$n_files" ]; then
        mode=$try; break
    fi
done
echo "${mode:-none}" > "$d/screen.mode"
args=(--manifest "$m" --outdir "$d" --k 23 --threads "$threads")
if [ -n "$mode" ] \
    && ./ploidyspec.sh kmers "${args[@]}" \
    && ./ploidyspec.sh matrix "${args[@]}" \
    && ./ploidyspec.sh homeologs "${args[@]}" \
    && ./ploidyspec.sh structure --outdir "$d"; then
    touch "$d/screen.done"
else
    touch "$d/screen.failed"
fi
./ploidyspec.sh cleanup --outdir "$d" || true
