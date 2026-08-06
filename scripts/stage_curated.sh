#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module load samtools/1.20--h50ea8bc_0

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

# curated (plain-gzip, read-only) assembly paths, derived straight from meta/meta.tsv
# rather than a /tmp scratch file -- compute nodes don't share /tmp with the login node.
python3 - <<'PYEOF' > /tmp/curated_files.$$.tsv
import os
with open("meta/meta.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        species = parts[1]
        for p in parts[3:]:
            if not p.strip():
                continue
            path = p if os.path.exists(p) else (
                p[:-len(".fasta.gz")] + ".fa.gz" if p.endswith(".fasta.gz") else
                p[:-len(".fa.gz")] + ".fasta.gz" if p.endswith(".fa.gz") else p
            )
            if "/assembly/curated/" in path:
                print(f"{species}\t{path}")
PYEOF

while IFS=$'\t' read -r species path; do
    base=$(basename "$path")
    dest="data/${species}/${base}"
    mkdir -p "data/${species}"
    if [ -s "$dest" ]; then
        echo "skip (exists): $dest"
        continue
    fi
    echo "staging: $path -> $dest"
    tmp="${dest}.tmp"
    zcat "$path" | bgzip -@4 -c > "$tmp"
    mv "$tmp" "$dest"
    samtools faidx "$dest"
done < "/tmp/curated_files.$$.tsv"
rm -f "/tmp/curated_files.$$.tsv"
echo "STAGING DONE"
