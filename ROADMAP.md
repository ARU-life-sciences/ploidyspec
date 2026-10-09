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
3. [x] **Fold per-species metrics into the package**: partition_consistency,
       distance_ratio_cv / pair_depth_cv, genome partition, reordered heatmap become
       stages run by `all`; cross-species tables + PCA become `ploidyspec panel`.
4. [x] **Separate cluster glue**: `workflows/sanger/` for LSF/module scripts;
       manifests accept paths relative to the manifest file.
5. [x] **Defaults = what the panel used** (k sweep 11–23 etc.); `--cleanup` to
       drop k-mer tables after a run; document runtime/memory.
       *Done with one deliberate deviation:* `--k` defaults to `15,23`, since
       k=23 was chosen for all 53,309 panel pairs (identical distances, ~1/3 the
       k-mer work). README documents memory (~40 MB per Mb of longest chromosome),
       runtime and disk.
6. [x] **Simulated test data + end-to-end test + CI**: small synthetic genomes —
       diploid, autotetraploid, allotetraploid, and a *partially rediploidized*
       autotetraploid (one fused/disomic chromosome, one arm-level disomic region).
       Doubles as the paper's validation figure.
7. [~] **Auto-detect naming** so manifests/regexes become optional (long-term
       goal: point it at raw FASTA(s)). Done (2026-10-08): chromosome-scale
       sequences found by length and numbered without names
       (`--chrom-naming auto`, default fallback), numbers matched across
       haplotypes by k-mer one-to-one matching; end-to-end test on unnamed,
       shuffled simulated assemblies. Contig-level haplotypes: placement on a
       scaffolded haplotype with `workflows/prep/scaffold_by_reference.py`.
       Still open: separating haplotypes inside one unlabelled file (manifest
       labels still needed), and folding the scaffolding into `prepare`.

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

- [x] Choose core outputs, one question each (option B, below):
      copy number · ancient pairing · copy divergence · lineage/inheritance
      structure · repeat-content divergence · **rediploidization state**.
      Everything else → supplementary.
- [x] Diploid anchor: `ddMalSylv1` (confirmed diploid by the user). Allele
      distance 0.0038, cross-chromosome 31x that, no pairs, no partition;
      per-chromosome te_frac up to 0.25 (see INTERPRETATION.md corrections).
- [ ] Rediploidization anchors: snow carps (fusions + chr17 arm-level);
      `drLytSali1` (tetrasomic by classical genetics); ideally a salmonid, where
      residual-tetrasomy regions are already mapped from segregation data —
      `windowed-homeologs` should recover them even from one haplotype.
- [x] Re-run the 4 stale species on their Darwin reassemblies (submitted
      2026-10-07 from scratch, plus `ddLepDrab1`; jobs in
      `logs/rerun4_jobs.txt`, old dirs in `superseded/`).
- [x] Fix `windowed` (2026-10-07): position-free. Each copy's windows are
      looked up in each other copy's whole-chromosome k-mer table (FastK
      `-p:table` + Profex), distance -ln(c)/k; rows directional. Lineage
      splits are read along every copy and the median copy reported.
- [x] Option B in code (2026-10-07): `structure` writes `genome_partition.tsv`
      and `pair_synchrony.tsv` only (other inheritance metrics retired, last
      values in `meta/archive/`); `panel` writes `panel_summary.tsv` plus
      `supplementary/` (TE-marker table, PCA over core numbers); report puts
      core sections first and TE markers under Supplementary.
- [x] Simulations now include what broke on real data: a deletion and an
      inversion in every haplotype but HAP1 (non-collinear copies), a
      `mislabelled` scenario (shifted chromosome numbering in one file, a
      chr05/chr06 swap in another), and a re-run test on a changed assembly
      in the same output directory. Cache signatures now include the source
      file's size and mtime.
- [x] Panel-wide windowed re-run (`workflows/sanger/rerun_windowed.sh`), then
      re-read `split_extent`-dependent states and regional calls.
- [x] `ddHesMatr1`'s chr01/chr02 swap: mutual swaps now need 2.0x per unit
      instead of 3.0x (user decision 2026-10-07). Matrix onward re-running
      (LSF 332545).
- [x] Core outputs: option B (user decision 2026-10-07) — copy number, copy
      divergence, ancient pairing, lineage split with where-along-the-
      chromosome, genome partition, fusions, rediploidization state; TE
      markers as a within-genome contrast in supplementary; drop the
      inheritance metrics except one ancient-pair synchrony measure; PCA
      supplementary. Needs the windowed fix.
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

---

## Progress log

Newest first. Check a running batch with the command given for it.

- **2026-10-08 — residual-tetrasomy test recalibrated.** The panel rerun with
  controls gave five new hits; three were calibration failures. Four-copy
  assemblies with a divergent copy (`ddLepDrab1`, `ddHypMacu1`) had an
  "allelic" level of 0.017-0.021, so 3x that sat at the homeolog distance;
  near-homozygous `laPotCris1` had a shuffled null of 3 windows. Fixes: the
  allelic level is the closest other copy, median over the chromosome's
  copies; a near-allelic window must also be < 0.25x the species' median
  homeolog distance (masu's near-allelic windows: median 0.0045 vs 0.048
  homeolog); the null is floored at 8 windows; terminal-only (>= 90%) hits are
  low confidence in the summary. Masu unchanged (8/10, 54 Mb), Salix 0;
  Lepidium, Hypericum, Potamogeton now 0; Lathraea chr03/chr04 kept (interior,
  26 Mb); Glechoma down to 2 small interior segments.

- **2026-10-08 — residual-tetrasomy test for two-haplotype assemblies.**
  `ploidyspec/residual.py`, run by the rediploidization stage when
  windowed-homeologs has run. Per window, the closest homeolog copy against
  the copy's typical allelic distance (< 3x = near-allelic), unless the window
  is equally close to an unrelated control chromosome
  (`windowed_homeolog_controls.tsv`, new; `windowed-homeologs --controls-only`
  adds it to old runs). Runs tested against a genome-wide shuffle null and
  decided per homeolog pair by more than half its copies. Pooled pairs keeping
  residual tetrasomy become `partially_resolved`; the summary reads them.
  - Masu salmon (`OncMaso1`): 8/10 paired chromosomes, 54 Mb, 62% at
    chromosome ends; chr05-chr28 diverged throughout. Lien et al. 2016
    (Atlantic salmon) is the reference for what to expect.
  - The control was needed: without it `ddSalPent1` (Salix, salicoid WGD)
    showed 10 chromosomes with end-to-end near-identity (0.001, below allelic);
    411 of its 521 near-allelic windows are shared with the control, leaving 0.
  - Simulation `residual_tetrasomy` (new): chr1/chr8 terminal stretch found
    from all 4 copies, a shared satellite at every chromosome start rejected
    by the control, the allotetraploid's pairs clean.
  - Still to do: re-run the panel with controls
    (`workflows/sanger/rerun_residual_controls.sh`); an arm-level homeolog map
    (salmonid homeology is per arm, so whole-chromosome pairing finds 5 masu
    pairs).
- **2026-10-08 — non-plant Darwin Tree of Life batch submitted.** 12 species
  in `jobs.tsv` (job IDs in `logs/new_species_jobs.txt`): Arctic charr
  (`fSalAlp3`) and whitefish (`fCorLav1`), salmonids; four Stylommatophora
  land snails (`xgCepNemo3`, `xgHygCinc1`, `xgMonCant1`, `xgDauRufa1`;
  McHale et al. report a WGD in the clade) with three non-stylommatophoran
  gastropod controls (`xgPomEleg1`, `xgLitLitt3`, `xgStaPalu1`); gar and
  bowfin (`fLepOcu1`, `fAmiCal2`), before the teleost WGD; and the weevil
  `icStrMela3`, triploid AAB by ToL's k-mer ploidy plot, assembled as two sets
  in hap1 (chr1-10, chr11-20) and one in hap2. Weevils `icOtiRugo1` and
  `icPolImpr1` were checked and are diploid. `xgCatTriz1` dropped
  (contig-level). `xgMonCant1` and `xgLitLitt3` run with
  `--chrom-naming auto` (unnamed hap2), which `jobs.tsv` does not record.

- **2026-10-08 — automatic chromosome numbering (1.7, part).** `prepare`
  numbers chromosome-scale sequences without names when headers carry none
  (or with `--chrom-naming auto`); the matrix stage matches the other
  haplotypes to the reference by k-mer distance. Unnamed, shuffled simulated
  allo- and autotetraploids recover the true grouping and states. Masu
  salmon (`OncMaso1`) added: HAP1 contigs placed on HAP2's 33 chromosomes
  with `workflows/prep/scaffold_by_reference.py`; first job died at 48 GB in
  minimap2 (`-c` base-level alignment), re-running without it.

- **2026-10-08 — windowed splits tested against a null.** After the
  position-free fix, single short segments (2-3 Mb, the minimum) on
  chromosomes with whole-chromosome split ~1.0 were inflating
  `one_divergent_copy`/`candidate`. The windowed scan is now repeated for
  groupings of the copies that cut across the observed one; pooled over the
  genome they give a per-species noise level, and a split extent counts only
  above that null's 95th percentile and above the chromosome's own crossing
  groupings. Removed 26 calls in six species (`daBudDavi1` 10, `daGalBore1`
  4, `ddSalTria1` 4, `drAriEdul1` 4, `SchCurv1` 3, `SchYoun1` 1), all towards
  `tetrasomic_like`/`candidate`; snow carp chr19 and chr17 and the simulated
  regional split pass.
  - `one_divergent_copy` now needs chromosome-wide evidence (distance or TE
    split); a regional segment isolating one copy is ordinary haplotype
    structure in a polysomic genome. `drLytSali1` (re-run on its new assembly)
    goes from 12/15 one_divergent_copy to 13/15 tetrasomic_like, matching its
    known tetrasomic inheritance; 32 calls change in nine species.
  - `all --with-windowed-homeologs` no longer stops when a species has no
    homeolog pairs (it aborted `dcCerAlpi1` before structure/te-markers/
    rediploidization); finishing job submitted.
  - Re-runs landed: `ddLepDrab1` (41 units relabelled; 8 pairs at uniform
    depth, CV 0.06; mostly tetrasomic-like within each number: AAAABBBB-like
    on x=8), `ddHesMatr1` (swap fixed; the chr01-chr02 "pair" is gone),
    `drLytSali1`, `lpElePalu1` (38.5 GB peak), `ddHypMacu1`.

- **2026-10-07 — Phase 2 start: hand-checks and re-runs.** Details in
  INTERPRETATION.md "Phase 2 hand-check: corrections".
  - `ddLepDrab1`'s allo signal was HAP4 chromosome-numbering left
    uncorrected: the relabelling check now matches units to chromosome
    numbers one-to-one. Clean re-run submitted.
  - `ddHesMatr1` HAP1 chr01/chr02 swap found by hand; below the auto
    threshold. Its only homeolog pair is an artefact. Decision pending.
  - Cached per-unit FASTAs/k-mer tables were reused across assemblies:
    `drLytSali1`'s September run used the old assembly. `kmers` now
    rebuilds a unit whose source, name or length changed. `dcCerAlpi1`,
    `ddHypMacu1`, `drLytSali1`, `lpElePalu1` re-running from scratch (LSF
    330623–330626, `ddLepDrab1` 330864). Check with
    `bjobs -w | grep rerun`.
  - Windowed tracks compare equal coordinates and are out of register for
    most of most chromosomes, panel-wide. Fix proposed above.
  - `ddMalSylv1` adopted as diploid anchor.

- **2026-10-07 — rediploidization: sibling-enrichment fusion check and
  partition pooling.**
  - Long-copy fusion candidates now need the partner at >= 2.5x
    (`--sibling-excess`) its per-k-mer rate in the copy's least-enriched
    sibling. Closes follow-up 3: `drMyrSpic1` HAP1_chr13 (homeologs at 0.9x
    next to a 15 Mb fragment) and both `daPilAura1` candidates (0.6–1.9x)
    drop out; `SchCurv1`'s fusion is 3.4x. No panel species has a
    `candidate_partner_present` row any more.
  - Species with no accepted pairs but a significant, balanced k=2 genome
    partition (z >= 10, `--partition-z`) pool each chromosome with its
    reciprocal best match across the groups (`partition_pool`):
    `daInuConz1` 16/16, `dmRanRepe1` 16/16, `dcCerAlpi1` 32/36, all
    `resolved_lineages`. No other species changed.
  - **Caveat found on the way:** pooling any two chromosome numbers in the
    diploid `ddMalSylv1` also reads `resolved_lineages` (22–47x). A pooled
    `resolved_lineages` therefore means only "not interchangeable"; the
    partition z (or the homeolog FDR) carries the evidence. Documented in
    `OUTPUTS.md`. For the paper, report pooled states as a test for
    tetrasomy, not as evidence of allopolyploidy.

- **2026-10-07 — panel refresh on the merged code.** Homeolog pooling
  (2-copy chromosomes pooled with their accepted partner), the
  `one_divergent_copy` state, ROADMAP 1.3–1.5 and the homeolog min-effect
  floor (`--min-effect 0.02`; stops chance pairs in a simulated diploid,
  changes no accepted pair on the panel) are all on `main`. All 8 end-to-end
  tests pass. `structure` and `rediploidization` re-run for all 45
  species with no failures; `meta/` regenerated with `ploidyspec panel`.
  - `auto_allo_spectrum.tsv` reproduces exactly. Partition z-values moved
    by up to ~±2 now that the null is seeded per species. `ddHesMatr1`'s
    weak k=3 partition and `daLatClan1`'s k=2 drop out; `dcCerAlpi1`'s best k
    is now 3 (z=20.8) on the same 16/20 core. Citations in
    `INTERPRETATION.md` updated.
  - With pooling, every allo anchor reads `resolved_lineages`, and so do the
    *Potamogeton* species and `llColAutu1`. `daPilAura1` reads one
    divergent copy in four on all 9 base chromosomes (AAAB-like).
    `lpElePalu1` is mixed.
  - Still `not_assessable`: two-haplotype species with no pairs, including
    the cryptic allo candidates. Next step for the stage: pool by partition
    group, not only by pair.
  - Write-up: panel results are in `INTERPRETATION.md` under "Detecting
    chromosome fusion".

- **2026-10-07 — first review of the panel rediploidization run** (40/43
  species done; `ddSalCine1`, `drIngLaur1` running, wheat resubmitted at 32 GB
  after `TERM_MEMLIMIT` at 4 GB).
  - **Fusions:** confirmed only in the snow carps. Three
    `candidate_partner_present` hits, both explained: `daPilAura1` HAP1_chr01
    (273 Mb vs 185 Mb sibling) contains its very young homeolog chr14 and
    chr06/chr07 material (HAP2_chr07 is short), likely a scaffolding difference
    between haplotypes; `drMyrSpic1` HAP1_chr13 is only "long" because HAP2_chr13
    is a 15 Mb fragment, and its chr03/chr09 hits look like the third member of
    the chr03–chr09 homeolog set, consistent with the AABBCC allohexaploid
    origin. The partner-missing rule kept all three from being called fusions.
  - **`drLytSali1`: 11/15 chromosomes `tetrasomic_like`**, matching the
    tetrasomic inheritance known from classical genetics — the magnitude
    metrics read this species as allo-like, this one doesn't.
  - **Recurring outlier haplotypes:** `ddHypMacu1` HAP1 on 8/8 (the known
    HAP1 assembly-quality problem — validates the summary metric);
    `ddLepDrab1` HAP4 on 15/16 — this is what drove its 0.94
    `partition_consistency`, so its "strongest allo candidate" status needs
    re-reading as one consistently divergent haplotype (assembly or a
    different-origin copy) rather than a two-subgenome split; `ddSalTria1`
    HAP3 on 8/19.
  - **Partial structure on many chromosomes:** `drAriEdul1` (5 partial, 3
    candidate of 17) and `ddEmpNigr1` (6 partial, 4 candidate of 13) — a
    partly diploidized autopolyploid reading for `ddEmpNigr1` would explain
    its literature-auto vs high-`te_frac` tension. Not yet checked by hand.
  - **Most of the panel is `not_assessable`** (two-haplotype assemblies).
    But `ploidy_ancestry_summary.tsv` already shows a pattern there:
    `daPilAura1`, `lpElePalu1` (some pairs) have `distance_ratio` ~1 — the
    homeolog partner is as close as the other copy of the same number, i.e.
    the four copies across the two numbers look interchangeable, while the
    allo anchors sit at 3.4–4.6. See follow-ups.
  - **Follow-ups for the stage** (do in the `simulate` worktree, test on
    simulations):
    1. Assess two-haplotype assemblies by pooling each homeolog pair's copies
       into one group (copies of chrA + chrB) and running the same lineage
       reading on it. This would turn most `not_assessable` rows into a real
       state.
    2. Triploids (`ddSalTria1`): with 3 copies every split is 1-vs-2, so
       `single_copy_outlier` is the wrong label when the same haplotype
       recurs (possible AAB structure). Give 3-copy splits their own state.
    3. Long-copy fusion test: don't call a copy "long" when its only sibling
       is a fragment (`drMyrSpic1` HAP2_chr13 at 15 Mb).
    Status: 1 done (homeolog pooling, plus partition pooling), 2 done
    (`one_divergent_copy`), 3 done (sibling-enrichment check) — see the
    entries above.
- **2026-10-07 — panel-wide rediploidization run submitted** (Phase 1.2
  follow-up). 43 LSF jobs (`redip_<species>`, 2 cores, 4 GB), one per species
  with a matrix; snow carps were already done. Job IDs in
  `logs/redip_panel_jobs.txt`. Check progress with
  `workflows/sanger/rediploidization_status.sh` (per-species state + non-zero
  `copy_state` counts + distinct fusions; totals on the last line). Results land
  in `results/<species>/rediploidization/` and each `report.html`.
  Species with many chromosome-length unplaced scaffolds (`ddSalCine1` 38,
  `ddPopNigr1`, `drIngLaur1`, the wheat runs) take longest: ~2.7 min per scaffold.
  **To do when finished:** review fusions and `candidate_partner_present` hits
  and false-positive states across the panel, then write the results into
  `INTERPRETATION.md`.
- **2026-10-07 — rediploidization stage added** (`60c333b`). Validated blind
  on the snow carps: `SchCurv1` chr19+22 (HAP3/HAP4) and the regional chr17
  split; all five `SchYoun1` fusions, ranked by `lineage_divergence` in the
  published wave order.
- **2026-10-07 — Phase 0 done, packaging done** (`d127213`, `4f756ba`,
  `4d0ca74`).
