#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0
cd /lustre/scratch122/tol/teams/blaxter/users/mb39/ploidyspec
K_SWEEP="11,13,15,17,19,23"
./ploidyspec.sh kmers --manifest manifests/lpElePalu1.tsv --outdir results/lpElePalu1 --k "$K_SWEEP" --threads 8
./ploidyspec.sh matrix --manifest manifests/lpElePalu1.tsv --outdir results/lpElePalu1 --k "$K_SWEEP" --threads 8
./ploidyspec.sh homeologs --manifest manifests/lpElePalu1.tsv --outdir results/lpElePalu1
