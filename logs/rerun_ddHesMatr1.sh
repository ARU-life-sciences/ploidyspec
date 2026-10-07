#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module load fastk/1.2-c1 samtools/1.20--h50ea8bc_0
cd /lustre/scratch122/tol/teams/blaxter/users/mb39/ploidyspec
export MPLCONFIGDIR=${TMPDIR:-/tmp}
m=manifests/ddHesMatr1.tsv; d=results/ddHesMatr1; T=${LSB_DJOB_NUMPROC:-8}
R=(--chrom-regex 'chromosome:?\s*(\d+)\b' --chrom-regex 'SUPER[_-](\d+)_HAP\d+\b' --chrom-regex 'SUPER[_-](\d+)\b(?!_)' --hap-regex '(?:SUPER_\d+_)?HAP(\d+)(?:_SUPER)?')
K=(--k 11,13,15,17,19,23)
./ploidyspec.sh matrix --manifest $m --outdir $d "${K[@]}" --threads $T "${R[@]}"
./ploidyspec.sh homeologs --manifest $m --outdir $d "${K[@]}" "${R[@]}"
./ploidyspec.sh te-markers --manifest $m --outdir $d --threads $T "${R[@]}"
./ploidyspec.sh subgenome-report --manifest $m --outdir $d "${R[@]}"
./ploidyspec.sh structure --outdir $d
./ploidyspec.sh rediploidization --manifest $m --outdir $d "${K[@]}" --threads $T "${R[@]}"
./ploidyspec.sh report --outdir $d
