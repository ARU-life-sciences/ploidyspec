#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0
cd /lustre/scratch122/tol/teams/blaxter/users/mb39/ploidyspec
./ploidyspec.sh te-markers --manifest manifests/lpTriTurg1_AB.tsv --outdir results/lpTriTurg1_AB --marker-k 13 --min-count 100 --min-ratio 2.0 --threads 8
