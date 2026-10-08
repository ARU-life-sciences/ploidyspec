#!/usr/bin/env bash
set -euo pipefail
source /etc/profile.d/modules.sh
module load minimap2/2.28--he4a0461_1 samtools/1.20--h50ea8bc_0 fastk/1.2-c1
cd /lustre/scratch122/tol/teams/blaxter/users/mb39/ploidyspec
d=data/salmonid
T=${LSB_DJOB_NUMPROC:-16}
export MPLCONFIGDIR=${TMPDIR:-/tmp}
# 1. place HAP1 contigs on HAP2's chromosomes
if [ ! -s $d/hap1_vs_hap2.paf ]; then
  minimap2 -x asm5 -c -t $T $d/masu_hap2.fa $d/masu_hap1.fa > $d/hap1_vs_hap2.paf.part
  mv $d/hap1_vs_hap2.paf.part $d/hap1_vs_hap2.paf
fi
python3 workflows/prep/scaffold_by_reference.py --contigs $d/masu_hap1.fa \
  --reference-fai $d/masu_hap2.fa.fai --paf $d/hap1_vs_hap2.paf \
  --out $d/masu_hap1.scaffolded.fa --placements $d/hap1_placements.tsv
samtools faidx $d/masu_hap1.scaffolded.fa
# 2. HAP2: add 'chromosome: N' so both files use the same naming
awk '/^>/{ if (match($1, /^>group([0-9]+)$/, m)) print $0" chromosome: "m[1]; else print; next } {print}' \
  $d/masu_hap2.fa > $d/masu_hap2.named.fa
samtools faidx $d/masu_hap2.named.fa
# 3. ploidyspec
mkdir -p results/OncMaso1
./ploidyspec.sh all --manifest manifests/OncMaso1.tsv --outdir results/OncMaso1 \
  --threads $T --with-te-markers --with-windowed-homeologs
