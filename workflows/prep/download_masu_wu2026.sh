#!/usr/bin/env bash
set -euo pipefail
cd /lustre/scratch122/tol/teams/blaxter/users/mb39/ploidyspec/data/salmonid/wu2026
for p in GCA_052849675.1_ASM5284967v1 GCA_052849655.1_ASM5284965v1; do
  a=${p:4:3}/${p:7:3}/${p:10:3}
  for f in ${p}_genomic.fna.gz ${p}_assembly_report.txt; do
    [ -s $f ] || wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/$a/$p/$f
  done
done
md5sum *.gz
