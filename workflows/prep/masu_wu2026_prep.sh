#!/usr/bin/env bash
# Wu et al. 2026 (Sci Data) phased Hi-C masu assembly, replacing our
# self-scaffolded pair (masu_salmon_prep_and_run.sh). HAP1 = GCA_052849675.1
# (200 scaffolds), HAP2 = GCA_052849655.1 (50 scaffolds), as in the paper.
# Download: download_masu_wu2026.sh. NCBI headers are kept; chromosomes are
# read with --chrom-regex 'chromosome (\d+),' as for fSalAlp3.
set -euo pipefail
source /etc/profile.d/modules.sh
module load samtools/1.20--h50ea8bc_0
cd /lustre/scratch122/tol/teams/blaxter/users/mb39/ploidyspec/data/salmonid/wu2026
for pair in HAP1:GCA_052849675.1_ASM5284967v1 HAP2:GCA_052849655.1_ASM5284965v1; do
  h=${pair%%:*}; acc=${pair#*:}; out=masu_wu2026_${h,,}.fa
  [ -s $out.fai ] && continue
  zcat ${acc}_genomic.fna.gz > $out.part && mv $out.part $out
  samtools faidx $out
  echo "$out: $(grep -c 'chromosome [0-9]*,' $out) chromosomes"
done
