# ploidyspec roadmap

Written 2026-10-07 after a full review of the code, docs and the 47-run panel.
Three goals, in priority order:

1. **Software anyone can clone/install and run** on their own haplotype assemblies.
2. **A locked-down method** with a small, defensible set of core outputs, calibrated
   against anchors.
3. **A write-up** (paper + methods docs) that makes the whole thing easy to hold in
   your head.

Two scientific axes run through all of it, and every core output should say which
one it measures:

- **Origin** — auto ↔ allo (the polyploid continuum, Twyford et al. 2025).
- **Rediploidization** — how far, and *where*, the genome has returned to disomic
  (diploid-like) inheritance: ancient pairs resolved, lineages split by fusion,
  arm-level/local resolution, still-tetrasomic chromosomes.

Rediploidization is a core feature, not a side reading. Right now it's "a pattern
read across three lines of evidence" (`INTERPRETATION.md`) and fusion detection is
a manual cross-reference. Making it a first-class, automated output is part of
Phase 1, not a later extra.

---

## Current state (summary of the review)

**Software**
- Core package (`ploidyspec/`, ~3.6k lines, stdlib + numpy + matplotlib) is
  clean, resumable, 91 unit tests passing (`python3 -m unittest discover -s tests`).
- Not installable: no `pyproject.toml`, no entry point, undeclared Python deps,
  FastK must be found by hand.
- About half the analysis lives in `scripts/` and is not run by `ploidyspec all`:
  `auto_allo_spectrum.py` (partition_consistency, distance_ratio_cv, pair_depth_cv,
  flip_rate), `genome_partition.py`, `reorder_heatmap.py`, `poly_space_pca.py`.
- Sanger-specific glue (`module load`, `bsub`, lustre paths, one-off `_rerun_*`)
  is mixed with the tool; manifests use absolute lustre paths.
- Defaults differ from what the panel used (README `--k 15` vs panel `11..23` sweep).
- No end-to-end test (real FastK on a tiny FASTA), no CI.
- Chromosome/haplotype naming still regex-driven.
- `results/` is 1.6 TB (mostly k-mer tables); no cleanup option; memory needs
  undocumented (`ddHesMatr1` needed 40 GB).

**Science / results**
- 44 species (47 result dirs, wheat split A/B/AB). `te-markers` run everywhere it
  can be — the 7 missing species all have only one usable haplotype.
- Stale `te_frac` (prior assembly): `dcCerAlpi1`, `ddHypMacu1`, `drLytSali1`,
  `lpElePalu1`.
- **No usable diploid anchor.** `ddAraThal4`'s alternate haplotype is fragmentary,
  so only HAP1 (5 units) is analysed — no within-pair baseline possible.
  Diploid-looking 2-copy species span `te_frac` 0.03–0.135, so the TE floor is
  species-specific; `ddMalSylv1` (lit. diploid, 2 copies) is the best candidate.
- `daEupConf1` has one usable haplotype, so the planned Tetmer head-to-head needs
  a different species (or Tetmer on reads).
- ~10 metrics, some failed (run-length/`flip_rate`), some confounded
  (`distance_ratio_cv` denominator). Needs pruning to a core set.
- Open threads: ~~inverted direction of the fusion-wave distance ranking~~ —
  resolved by the rediploidization stage: fused-lineage divergence
  (`fusions.tsv` `lineage_divergence`) ranks `SchYoun1`'s five fusions in the
  paper's wave order; `SchYoun1` lineage split on fused scaffolds (fused copies
  are detected but not yet fed into te-markers/windowed); `drAriEdul1` chr01
  (~13 Mb breakpoint); `ddSalTria1` and `drRosSpin1` ploidy conflicts; PCA
  imputation.
- `INTERPRETATION.md` (82 KB) is a lab notebook — methods, results and
  corrections interleaved.

---

## Phase 0 — housekeeping

- [x] Commit the outstanding work (poly-space PCA, field-positioning section,
      per-species READMEs).
- [x] Gitignore `logs/`.
- [x] Move one-off scripts (`_rerun_*`, `_repair_*`, `_run_*`, `_setup_darwin_batch.py`)
      to `scripts/archive/`.

## Phase 1 — software runnable by others

1. [x] **Packaging**: `pyproject.toml` with a `ploidyspec` console script;
       `environment.yml` (python, numpy, matplotlib, samtools, FastK); README
       install section rewritten around `pip install` / conda. Container later.
2. [x] **Rediploidization as a core stage** (see design below).
3. [ ] **Fold per-species metrics into the package**: partition_consistency,
       distance_ratio_cv / pair_depth_cv, genome partition, reordered heatmap become
       stages run by `all`; cross-species tables + PCA become `ploidyspec panel`.
4. [ ] **Separate cluster glue**: `workflows/sanger/` for LSF/module scripts;
       manifests accept paths relative to the manifest file.
5. [ ] **Defaults = what the panel used** (k sweep 11–23 etc.); `--cleanup` to
       drop k-mer tables after a run; document runtime/memory.
6. [ ] **Simulated test data + end-to-end test + CI**: small synthetic genomes —
       diploid, autotetraploid, allotetraploid, and a *partially rediploidized*
       autotetraploid (one fused/disomic chromosome, one arm-level disomic region).
       Doubles as the paper's validation figure.
7. [ ] **Auto-detect naming** so manifests/regexes become optional (long-term
       goal: point it at raw FASTA(s)).

### Rediploidization stage — design sketch

Everything below already exists as data; this makes it automatic and puts it in
one table.

- **Fusion/fission detection** (currently manual, `INTERPRETATION.md` §"Detecting
  chromosome fusion"): per haplotype set, copies at ~2× sibling length + a
  chromosome number missing from exactly those haplotypes → candidate fusion;
  confirm by windowed match confined to one block. Mirror logic for fission.
  Validate: recover all 5 `SchYoun1` fusions and `SchCurv1` chr19+22 without
  header hints.
- **Per-chromosome lineage structure**: the distance-based 2-way split + TE
  `split_ratio` already in `te_marker_fraction_by_lineage.tsv`, extended to the
  windowed track so it reports *where* along the chromosome the split holds
  (whole chromosome vs one arm vs a block — chr19 vs chr17 vs `drAriEdul1` chr01).
- **Ancient-pair resolution**: `distance_ratio` per homeolog pair and its spread
  (one synchronous event vs asynchronous resolution).
- **Output**: `rediploidization/rediploidization_by_chrom.tsv` — one row per
  chromosome number: copies, fusion status, lineage split (distance + TE), extent
  (whole/arm/local/none), ancient partner + distance_ratio, and a state label:
  *tetrasomic-like* / *partially resolved* / *resolved into disomic lineages* /
  *anciently resolved pair*. Plus a genome-level summary (fraction of genome per
  state, spread of resolution).
- Labels are descriptive readings of the data, not proof of inheritance mode
  (single individual, no segregation data) — said explicitly in the output docs.

## Phase 2 — lock down the method

- [ ] Choose core outputs, one question each:
      copy number · ancient pairing · copy divergence · lineage/inheritance
      structure · repeat-content divergence · **rediploidization state**.
      Everything else → supplementary.
- [ ] Diploid anchor with two chromosome-scale haplotypes (`ddMalSylv1` after
      literature check, or a new outbred diploid).
- [ ] Rediploidization anchors: snow carps (fusions + chr17 arm-level);
      `drLytSali1` (tetrasomic by classical genetics); ideally a salmonid, where
      residual-tetrasomy regions are already mapped from segregation data —
      `windowed-homeologs` should recover them even from one haplotype.
- [ ] Re-run the 4 stale species on their Darwin reassemblies.
- [ ] Close or explicitly park the open threads listed above.

## Phase 3 — write-up

Split `INTERPRETATION.md` into: methods doc (from the code), results, and an
archived lab notebook. Paper skeleton:

- **Intro**: polyploid continuum; rediploidization as the time axis of that
  continuum; need for reference-free, single-individual profiling at ToL scale.
- **Methods**: one subsection per core output.
- **Validation**: simulations; allo anchors (wheat A/B, `daGleHede1`,
  `drTriRepe1`, `drSorDevo1`, `drMyrSpic1`); auto anchors (snow carps — blind
  recovery of fusions, chr19 lineage split, chr17 — and `drLytSali1`); diploid
  anchor.
- **Results**: panel survey on both axes; case studies — snow carp
  fusion-triggered rediploidization, cryptic young allopolyploids
  (`dmRanRepe1`, `dcCerAlpi1`, `daInuConz1`), segmental `ddHesMatr1`,
  demi-duplication `llColAutu1`, holocentric `lpElePalu1`.
- **Discussion**: divergence magnitude vs inheritance mode (diploidized autos look
  allo — Twyford Box 2, Xie 2026); absolute vs within-genome readings; assembly
  curation artefacts; single-individual limits; comparison with Tetmer.

## Phase 4 — more clades, targeted

Hold broad expansion until Phase 2. Prioritise species that add *anchors*:
confirmed diploids with two chromosome-scale haplotypes, confirmed autopolyploids,
and systems with mapped rediploidization (salmonids, other fusion-rich
autopolyploid lineages).
