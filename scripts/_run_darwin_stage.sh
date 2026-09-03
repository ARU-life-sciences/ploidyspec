#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module load samtools/1.20--h50ea8bc_0
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
python3 scripts/_setup_darwin_batch.py
