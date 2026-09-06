# ploidyspec — what it measures, and where rediploidization fits

A working reference for collaborators: what each pipeline stage computes, how
auto/allopolyploidy and rediploidization progress show up in the data, where
the method's edges are, and a species-by-species summary of the panel so far.
Column-by-column file reference lives in [`OUTPUTS.md`](OUTPUTS.md); this
document is the conceptual layer above it.

## In one paragraph

Given two or more haplotype-phased genome assemblies from *one individual*,
ploidyspec compares canonical k-mer sets between every chromosome copy —
genome-wide, in sliding windows along each chromosome, and in their
repeat/TE content specifically. No alignment or reference genome is needed.
Auto/allopolyploidy, rediploidization progress, and chromosome fusion are
different *readings* of that same underlying distance data, not separate
measurements — which is the main source of confusion, so each section below
says explicitly which reading it is.

## The pipeline

Six stages, run roughly in this order: `prepare` → `kmers` → `matrix` →
`homeologs` → `windowed` → `report`, with `windowed-homeologs`,
`te-markers`, and `te-markers-windowed` as opt-in branches feeding back into
the report.

**`matrix` — contemporary structure.**
*Answers: which chromosomes have how many copies, and how close are those
copies to each other?* Mash-corrected k-mer distance, all-vs-all, between
every chromosome-scale sequence. `ploidy_summary.tsv` gives copy count per
chromosome number (`n_haplotype_copies`) — raw copy number, not inheritance
mode. The heatmap's *shape* is the auto/allo tell: a flat, uniform block
(all copies equally close) reads auto-like or homogenized; a split into
sub-blocks (some copies much closer to each other than to the rest) reads
allo-like — two or more independently-diverging lineages.

**`homeologs` — ancient retained-duplicate search.**
*Answers: does this chromosome still look like the sister it was duplicated
from in the original WGD, despite everything renumbering since?* Every
cross-chromosome-number pair tested against the genome-wide background
distance (z-test, Benjamini–Hochberg FDR). The key rediploidization number
is `distance_ratio` — the pair's distance divided by each chromosome's own
contemporary-copy distance. Near 1 means the "ancient" partner is barely
more diverged than a chromosome is from its own copies: not yet resolved
into two independent lineages. In the hundreds or thousands means long
since resolved.

**`windowed` / `windowed-homeologs` — fine-scale divergence.**
*Answers: is the whole-chromosome average hiding something local?* Same
metric, computed in sliding windows instead of one number per pair. A pair
that's mostly diverged but carries one block of near-identity is a
candidate recent or ongoing homeologous exchange (HE) — or, as we learned
below, a signature of a physically fused chromosome. The heatmap view
(rows = haplotype pairs, columns = position, color = distance) makes mosaic
patterns visible across many chromosomes in one image.

**`te-markers` / `te-markers-windowed` — repeat-content divergence.**
*Answers: have these two copies built independent repeat/TE histories, or
are they still sharing one?* Finds k-mers that are high-copy (repeat-derived)
and private to one copy versus the other. `te_marker_fraction` is the share
of a chromosome's repeat content that's subgenome-private rather than
shared — addressed on its own below. The windowed version paints each
window by which side's markers it matches, giving a direct picture of HE
tracks.

## Where rediploidization fits

Not one output column — a pattern read across three lines of evidence, all
already sitting in the files above:

1. **Contemporary shape** (`matrix`) — flat/uniform copies = not yet
   separated, or already fully homogenized back together; split blocks =
   separated.
2. **Ancient-pair resolution** (`homeologs`, `distance_ratio`) — per
   chromosome, and its *spread* across the genome. Tight and high = one
   synchronized event, long finished. Wide spread = pairs resolving at
   different times.
3. **Local and structural signal** (`windowed` heatmaps, `te-markers-windowed`,
   and literal karyotype change) — a block of near-identity inside an
   otherwise-resolved pair, or a chromosome fusion between haplotype copies
   within one individual — the most direct evidence available.

**Reading "has undergone" vs. "is undergoing":** complete rediploidization
looks like uniform, high `distance_ratio` across most or all chromosomes, no
local low-distance blocks in the windowed data, stable `te_marker_fraction`
genome-wide (`daGleHede1` is the clean case in this panel: 18/18 chromosomes
paired, all at ratio 3.2–4.6×). In-progress rediploidization looks like
heterogeneous `distance_ratio` within one genome, local exchange-looking
blocks in windowed heatmaps, and — the cleanest signal available — visible
chromosome fusion between haplotype copies *within the same individual*
(see the snow carp section below).

**Counterintuitive but important:** near-total *loss* of ancient-pair signal
doesn't mean nothing's happening. It can mean sequence-level divergence has
already outrun what k-mer comparison can detect, even while the karyotype is
still actively restructuring via fusion. Sequence-level and
structural-level rediploidization can be at visibly different stages in the
same genome — this is exactly what we see in the snow carp pair below.

## Are TE markers only useful for allopolyploids?

Mechanically, `te_marker_fraction` asks one question: have these two copies
built *independent* repeat histories? That's inherently shaped like an
allopolyploid question — two once-separate genomes, two separate TE
lineages since they split. But it's a real, informative test for
autopolyploids too, just read in the other direction:

- **High value** → evidence *for* independent divergence — allopolyploid, or
  an autopolyploid old enough to have stopped fully homogenizing.
- **Low value** → *consistent with* autopolyploid (one ancestral genome,
  copies still homogenizing via exchange), but not proof — a young
  allopolyploid, or one whose subgenomes have been heavily homogenized, can
  look identical. This was the panel's open calibration gap (see below —
  the snow carps now close it).

So: useful on both sides, as a continuous index rather than an allo-only
switch. For salmonid work specifically, where regions of residual tetrasomy
are already mapped from segregation data, `te_marker_fraction` computed per
chromosome arm gives an independent, sequence-only cross-check, and the
windowed painting gives a physical view of exactly where homeologous
exchange is (or isn't) still active.

**Now run on both snow carp species, and it immediately exposed a real
limitation in how the index was being read.** Genome-wide means came back
at 0.209 (`SchCurv1`) and 0.223 (`SchYoun1`) — indistinguishable from
`daGleHede1`'s confirmed-allo 0.223. Taken at face value, that says
`te_marker_fraction` can't tell a confirmed autopolyploid from a confirmed
allopolyploid, which would undercut the whole index. It can't, **at that
resolution** — a flat genome-wide mean is the wrong level to read a genome
that's only partly rediploidized. Broken down per chromosome and, where a
chromosome has ≥3 haplotype copies, by an automatic 2-way split of those
copies by whole-chromosome distance (new: `te_marker_fraction_by_lineage.tsv`,
described in `OUTPUTS.md`), the real signal appears immediately:
`SchCurv1`'s `chr19` (the confirmed chromosome fusion) shows within-lineage
`te_marker_fraction` of 0.11 versus 0.59 across lineages — a 5.6×
split invisible in either the flat genome-wide mean or that chromosome's
own flat pairwise average (~0.35). Every other chromosome sits at
split_ratio 0.8–2.2 (no comparable structure), except `chr17` (split_ratio
3.2 in `SchCurv1`, 1.9 independently in `SchYoun1`) — which the paper
separately documents as a *second* rediploidization mechanism at that exact
chromosome (a centromeric inversion producing a partial, short-arm-only
disomic pattern), found by this tool with no prior knowledge of that
result. See `OUTPUTS.md` for the caveats on treating `split_ratio` as a
diagnostic rather than a classifier.

**Applying the same lineage-split analysis to the rest of the panel** (any
species with `te-markers` already run, no prior knowledge of fusion
structure needed) surfaced comparably strong splits at specific chromosomes
in 4 other species: `daBudDavi1` chr05 (6.13), `ddHypMacu1` chr06 (4.73),
`drAriEdul1` chr01 (4.66), `drLytSali1` chr07 (3.67). Digging into
`daBudDavi1`'s chr05 (the strongest of the four) found it was **not** real
biology — `HAP2_chr05` is 32% short and has a third of its three siblings'
high-copy k-mer count, tracking them normally on every other chromosome, and
every pair involving it is lopsided in the same way `ddHypMacu1`'s `HAP1`
was. Same confound, different species, caught by inspecting length and
per-pair marker counts by hand. That's now automated: `te_marker_fraction_by_lineage.tsv`
carries a `flagged_units` column marking exactly this pattern (see
`OUTPUTS.md`) — deliberately built to *not* flag a real 2-vs-2 split like
`SchCurv1`'s `chr19` fusion (both lineages have >1 member, so no majority-vs-
minority signature exists to catch), only a clear majority-vs-minority split
like `daBudDavi1`'s.

Digging into the remaining three found a second failure mode. `ddHypMacu1`'s
chr06 has *normal* length but ~2× the repeat density of its three
mutually-agreeing siblings — invisible to a length/raw-count check, since a
longer sequence isn't the issue here. Corroborating evidence: the same
haplotype (`HAP1`) shows 1.9–5.1× elevated density on 3 other chromosomes
and 0.43× (deflated) on a fourth — scattered in both directions, consistent
with those being assembly noise rather than four real biological signals.
`HAP1` does independently show 570 excluded scaffold fragments vs. 2–3 for
`HAP3`/`HAP4` in `unplaced.tsv` — but that specific asymmetry turns out to
be *expected* for any >2-haplotype curated assembly, not evidence `HAP1`
was sequenced or assembled to a lower standard: ToL's curation pipeline
only runs unlocalized-sequence placement on the primary `HAP1`/`HAP2` pair,
so `HAP3`/`HAP4`+ files are built only from already-placed scaffolds and
will always show near-zero unplaced counts by construction. The scattered
density anomalies are still the more direct evidence for treating `HAP1`'s
signal cautiously here (see `OUTPUTS.md`'s assembly-quality-confound entry
for the full explanation and why the correlation, not just the asymmetry,
is what matters). Also now automated: a separate
`high_density_units` column, using repeat *density* (high-copy k-mers per
bp) rather than raw count — raw count alone would false-flag a real fusion
too, since a longer fused sequence simply contains more total repeat
content (confirmed against `SchCurv1`'s `chr19`: ~2.3–2.6× the raw count in
the fused lineage but only ~1.35–1.4× the density, correctly unflagged).

Final status of all four: `daBudDavi1` chr05, `drLytSali1` chr07, and
`ddHypMacu1` chr06 are all flagged as assembly artifacts, none of them real
signal. `drAriEdul1`'s chr01 is the one that survived — see below.

### `drAriEdul1` chr01 — the one that looks real

Passes both automated checks (`flagged_units` and `high_density_units` both
empty), and manual digging didn't find an artifact explanation either.
Whole-chromosome Mash-corrected distances tell a clean, quantitative story:

```
HAP1–HAP2: 0.0118    HAP2–HAP3: 0.0057 (closest pair)
HAP1–HAP3: 0.0118    HAP2–HAP4: 0.0268
HAP1–HAP4: 0.0260    HAP3–HAP4: 0.0266
```

`HAP1`/`HAP2`/`HAP3` form a genuinely tight trio (mutual distances
0.006–0.012); `HAP4` sits 2.2–4.7× further from all three than they are
from each other. The windowed heatmap looks at first glance like `HAP1` is
equally diverged from everyone — that's a scale-compression artifact of the
*raw*, uncorrected windowed track: its color scale gets stretched by
`HAP2`/`HAP3`'s unusually tight relationship (down to 0.14), making
everything else look uniformly dark by comparison even though `HAP1` is
genuinely closer to `HAP2`/`HAP3` than to `HAP4` on the calibrated metric.
Read correctly, the windowed track shows real spatial structure: `HAP2`,
`HAP3`, and `HAP4` sit relatively close for the first ~13Mb, then `HAP4`
breaks away sharply from both for the rest of the chromosome — a specific
transition point, not diffuse noise.

Why this one reads as real rather than artifactual, unlike the other three:
- `HAP4` is *not* the haplotype `unplaced.tsv` flags as fragmented — that's
  `HAP1`/`HAP2` (89 excluded fragments each, vs. 2 each for `HAP3`/`HAP4`).
  This is a 4-haplotype species, so per the curation-pipeline explanation in
  `OUTPUTS.md` (only the primary `HAP1`/`HAP2` pair goes through unlocalized-
  sequence placement), that asymmetry is expected by construction and isn't
  itself evidence `HAP1`/`HAP2` are lower quality — but it does mean `HAP4`
  is one of the two haplotypes *not* carrying that structural asymmetry
  either way, which is one less thing to explain away here.
- Length (+15%) and density (+66%, just under the 1.75× threshold) are
  both modestly and consistently elevated together — a different pattern
  from `ddHypMacu1`'s pure-density spike or `daBudDavi1`'s pure length
  deficit.
- There's a specific breakpoint (~13Mb), not uniform noise across the whole
  chromosome.

This species' own literature entry already flags taxonomic uncertainty —
the sequenced sample may actually be the tetraploid *A. wyensis* rather
than diploid *A. edulis* — and whitebeams (*Sorbus*/*Aria*) are a textbook
genus for hybrid/polyploid complexity, so a real, partial differentiation
signal at one chromosome is a plausible outcome, not a surprising one. Not
proven the way the snow carp findings were (no independent literature
confirmation here) — the strongest surviving candidate of the four, not a
confirmed finding. Resolving it further would need alignment/synteny data,
beyond what k-mer comparison alone can do.

## Detecting chromosome fusion (and fission)

Currently a manual cross-reference of existing files, not its own pipeline
stage:

1. Compare a chromosome's length across its own haplotype copies
   (`sequences.tsv`). One copy at roughly double (or some other clean
   multiple) of its siblings is a candidate fusion.
2. Cross-check `ploidy_summary.tsv`: is another chromosome number missing
   entirely from exactly the same haplotype(s)? That's the fusion-partner
   candidate.
3. Confirm with the windowed heatmap between the candidate scaffold and both
   candidate partners. A real fusion shows the low-distance match *confined
   to one block* — the part of the scaffold that really is the other
   chromosome — absent from haplotype copies that don't carry the fusion.

**Fission** would be the mirror image: a copy unusually short relative to
its siblings, paired with an extra chromosome number present only in that
haplotype set.

## The snow carps: what we found, what the paper found, and a correction

Two Schizothoracine species were run through the pipeline: *Schizothorax
curvilabiatus* (`SchCurv1`, 4 haplotypes, 25 chromosomes) and
*Schizopygopsis younghusbandi* (`SchYoun1`, 4 haplotypes split from parent-
of-origin-labelled `M1`/`M2`/`P1`/`P2`, 25 chromosomes). This was done
*before* reading Xie et al. 2026, *Nature*, "Chromosomal fusions trigger
rediploidization of autopolyploid genomes"
(doi:10.1038/s41586-026-10439-1) — which turns out to describe exactly this
system, with the same two species.

**What tallies:**
- **Chromosome counts.** We independently derived 98 chromosome-scale units
  for `SchCurv1` (4×25 minus 2, since chr22 is missing from 2 of 4
  haplotypes) and effectively 90 distinguishable units for `SchYoun1` — the
  paper states the cytogenetically-observed and assembled counts as exactly
  98 and 90.
- **The fusion set.** We found 5 chromosome fusions in `SchYoun1` purely
  from header-naming analysis of the raw FASTA (`Chr04+15`, `Chr08+16`,
  `Chr11+14`, `Chr19+22`, `Chr20+23`), and one in `SchCurv1` (`chr19+22`,
  found independently via length/absence detection, not header naming).
  The paper's own synteny/Hi-C/FISH-confirmed fusion set is exactly these
  same 5 pairs, and it independently confirms `chr19+22` as the one shared,
  ancestral fusion present in both species — matching our finding that it's
  the *only* fusion signature in `SchCurv1`, while `SchYoun1` carries all
  five.
- **Auto, not allo.** The paper establishes unambiguous autotetraploid
  origin for the whole subfamily (cytogenetics, outgroup synteny, meiotic
  FISH showing quadrivalent pairing). Our own read of the whole-chromosome
  heatmaps — flat, uniform block structure, no split subgenome pattern — is
  exactly the auto-like signature, arrived at independently. This is the
  panel's first confirmed-autopolyploid ground truth (previously we only
  had confirmed allopolyploids plus one true diploid).
- **The lineage-split TE-marker signal at `chr19`.** Once `te-markers` was
  run, splitting `SchCurv1`'s `chr19` copies into the unfused (`HAP1`/`HAP2`)
  and fused (`HAP3`/`HAP4`) lineages showed within-lineage
  `te_marker_fraction` of 0.11 versus 0.59 across lineages — the two
  post-fusion lineages have built independent repeat histories from each
  other, not just diverged in bulk sequence. This is the same
  fusion-splits-a-tetrasomic-quartet-into-two-lineages model the paper
  describes, recovered with a completely different method (differential
  repeat-content markers, not Ks/gene trees) at the same locus.
- **`chr17`, unprompted.** The same lineage-split analysis, run automatically
  across every chromosome (not just the known fusion), flagged `chr17` as
  the next-most-structured chromosome in both species (`split_ratio` 3.2 in
  `SchCurv1`, 1.9 in `SchYoun1`) — well below `chr19`'s fusion-strength
  signal, but clearly above the rest of the genome (0.8–2.2). `chr17` is
  *not* one of the five known fusions. The paper separately documents a
  distinct rediploidization mechanism there: a centromeric inversion
  producing a partial, short-arm-only shift to disomic inheritance, with the
  long arm staying tetrasomic — exactly consistent with a real-but-partial
  signal weaker than a full chromosome fusion. Not independently confirmed
  here (would need the same due-diligence check as any homeolog hit), but a
  strong candidate that fell out of the analysis without being told to look
  for it.

**A correction we made along the way:** we initially misread `SchCurv1`'s
`chr19↔chr22` ancient-homeolog hit as a pipeline artifact — a chimeric,
mis-scaffolded sequence to be discarded — because the matching windowed
heatmap block was localized rather than genome-wide. That reasoning was
backwards. Under the paper's own model, an **unbalanced chromosome fusion**
(only 2 of 4 ancestral copies fuse; the fused pair then pairs preferentially
with itself at meiosis instead of the original four-way pairing) is exactly
what triggers disomic inheritance and ohnologue divergence at that locus —
so a match confined to one block, present only in the haplotypes that carry
the fusion, is the *correct* signature of a real fused ohnologue pair, not
evidence against one. `results/SchCurv1/homeologs/ploidy_ancestry_summary.tsv`
and `OUTPUTS.md` have both been corrected to reflect this.

**What doesn't cleanly tally yet, worth a real conversation:** ranking our
raw distance/z-score for all 5 `SchYoun1` fusion pairs groups them into
exactly the paper's three waves (wave 1 alone, wave 2's three pairs
clustered, wave 3 alone) — a correspondence unlikely to be chance. But the
*direction* is inverted from naive expectation: `chr19↔chr22` (the paper's
most-diverged, most-resolved pair) shows the *lowest* raw k-mer distance of
the five; `chr11↔chr14` (least resolved) shows the *highest*. Since k-mer
Jaccard distance should increase monotonically with real divergence, this
needs unpacking — plausibly an artefact of averaging across a whole
chromosome when the paper's own results show divergence concentrated near
the fusion site and expanding outward incompletely, meaning most of a
partially-rediploidized chromosome's length may still show close to zero
signal and dilute the average. The windowed data (which we already have)
is the natural next place to check this properly, rather than the
whole-chromosome average.

**Also worth carrying forward:** the paper explicitly warns that some
downstream consequences of fusion-triggered rediploidization — biased gene
retention, biased expression — "mimic subgenome dominance in
allopolyploids," and that care is needed before reading asymmetric
divergence patterns as proof of allopolyploid origin on that basis alone.
That caveat applies beyond the fish — worth revisiting any panel species
called "allo" primarily from asymmetry rather than independent (e.g.
phylogenomic) evidence.

## Caveats to carry into any discussion

- No confirmed diploid-vs-autopolyploid calibration anchor existed until the
  snow carps — now the panel has confirmed allopolyploids, one confirmed
  autopolyploid pair, and one true diploid.
- Everything comes from a single individual's assembly. No pedigree or
  segregation data, so any inheritance-mode claim (multisomic vs. disomic)
  is indirect evidence, not direct proof.
- FDR correction gets more conservative as chromosome count grows. "No
  significant pair" on a many-chromosome genome doesn't mean no signal —
  check `homeolog_candidates_ranked.tsv` for near-misses (this is exactly
  what happened with `SchYoun1`'s `chr19↔chr22`: real signal, real q=0.26).
- Fused/chimeric scaffolds can produce genuine ancient-homeolog hits, not
  just artifacts — always check chromosome length + `n_haplotype_copies` +
  the windowed heatmap to tell which it is before concluding either way.
- k-mer similarity can't fully separate true shared ancestry from
  convergent repeat content (e.g. shared satellite DNA). A real,
  non-artifact, localized block of similarity could still be either.

## 2026-09 Darwin assembly batch

A batch of 15 species with new/updated Darwin Tree of Life `hap1`/`hap2`(+)
assemblies was run through `all` (+ `homeologs`/`windowed-homeologs`) in
early September 2026: 4 already in the panel with newer assembly versions
(`dcCerAlpi1`, `ddHypMacu1`, `drLytSali1`, `lpElePalu1`) and 11 new species
(`daSenVulg1`, `drSorDevo1`, `drMyrSpic1`, `drMyrVert1`, `daPilAura1`,
`ddLepDrab1`, `ddHesMatr1`, `ddSalTria1`, `drSorDevo1`, `dmRanRepe1`,
`daLatClan1`). `drSorAngl1` was excluded — only `hap1` exists for it
anywhere on the farm, so it can't be run. All 14 runnable species are now
complete and verified (`matrix/ploidy_summary.tsv` +
`homeologs/homeolog_pairs.tsv` populated, LSF "Successfully completed").
`ddHesMatr1` (tetraploid, 6 chromosomes) needed three resubmissions at
increasing memory (16GB → 28GB → 40GB) to clear `TERM_MEMLIMIT` during
k-mer table building — its actual peak usage tracked each new ceiling
closely (18.6GB, then 28.16GB) before finally completing under 40GB,
which may reflect this species' genuine k-mer/repeat-content memory
footprint rather than a fluke. **`te-markers` + `subgenome-report` have
since been run on all 10 new-species assemblies from this batch**
(2026-09-06) — the `te_frac` values below are real for those 10. The `te_frac`
values already in the table for the 4 re-run species are still carried
over from their *previous* (pre-Darwin-batch) assembly and are marked
stale below; treat those as historical context only, not as measurements
of the new assembly.

Two standout results from this te-markers run: **`ddLepDrab1`'s `te_frac`
is 0.637 — the highest of the 10, and higher than any confirmed allo
anchor in the panel** (`daGleHede1` 0.223, `drTriRepe1` 0.169). That's now
a *third* independent metric (alongside `partition_consistency` 0.94 and
`distance_ratio_cv` 0.132) converging on allo-like for this species —
the strongest multi-metric case in the panel outside the confirmed
anchors. **`drMyrVert1`'s `te_frac` is 0.521**, the second-highest of the
10 — a genuinely new finding, since no literature source addresses this
species' own origin (it's only ever used as an outgroup in *M. spicatum*
studies). `daSenVulg1`'s `te_frac` is very low (0.041) — its apparent
`distance_ratio_cv` conflict (0.903) turned out to be a metric artifact,
not a real disagreement (traced below in its species-table row); the low
`te_frac` on its own remains the one thing genuinely at odds with the
literature's allo lean.

The batch splits into three readable groups:

- **`dcCerAlpi1` and `dmRanRepe1` reconfirm/add clean, complete
  no-homeolog results** — `dcCerAlpi1` 0/36 chromosomes paired (same as the
  prior assembly), `dmRanRepe1` 0/16 — the strongest candidates in the
  panel for a genuine diploid/no-WGD data point, which is exactly what
  [[project_autoallo_calibration_gap]] flagged as missing. Neither's
  literature ploidy status has been checked yet, so these are data-only
  observations, not confirmed calibration anchors.
- **Near-fully-resolved ancient pairing** (most or all chromosomes find a
  partner): `drSorDevo1` (34/34, full), `drMyrVert1` (14/14, full),
  `daPilAura1` (18/18, full, tight/consistent divergence 0.031–0.044 —
  single-age WGD signature), `ddLepDrab1` (16/16, full, tetraploid —
  cleanly resolves into 8 chromosome pairs), `lpElePalu1` (18/19, unchanged
  from prior assembly), `drLytSali1` (10/15, unchanged from prior
  assembly), `daLatClan1` (18/21), `drMyrSpic1` (16/21), `daSenVulg1`
  (16/20).
- **Sparse pairing** (most chromosome copies look interchangeable rather
  than distinctly diverged): `ddHypMacu1` (2/8, unchanged in count from the
  prior assembly), `ddSalTria1` (4/19, triploid), `ddHesMatr1` (2/6,
  tetraploid — only `chr01↔chr02` accepted out of 6 chromosomes, the
  species that needed 40GB to complete), `drLytSali1` also qualifies here
  relative to its ploidy level (only 10/15 vs. `ddLepDrab1`'s full 16/16 at
  the same 4-copy ploidy).

All four re-run species reproduced the same accepted-pair *count* as their
prior-assembly runs (`dcCerAlpi1` 0/36, `ddHypMacu1` 2/8, `lpElePalu1`
18/19, `drLytSali1` 10/15) — a useful stability check, though pair
*identity* wasn't cross-checked chromosome-by-chromosome against the old
runs (the old `results/` output was overwritten in place before this
comparison was made).

## A continuous auto/allo spectrum from inheritance-mode metrics

Prompted by `drLytSali1`: `te_marker_fraction` and bulk k-mer divergence
both measure *magnitude* of divergence, which turns out to be a poor
auto/allo discriminator on its own — the reason `drLytSali1` looks
allo-like on both is that they can be confounded by the *same* mechanism
(tetrasomic recombination homogenizes bulk sequence identity but doesn't
homogenize TE insertion polymorphisms, so a genuine autopolyploid can show
near-zero bulk distance *and* elevated `te_marker_fraction` at once — see
its row below). What actually distinguishes allopolyploidy (disomic
inheritance: two subgenomes recombine only within themselves, staying
stable) from autopolyploidy (tetrasomic inheritance: all copies exchange
freely) is *consistency of lineage identity*, not divergence magnitude.

Four new statistics, computed by `scripts/auto_allo_spectrum.py` and
written to `meta/auto_allo_spectrum.tsv` — all reanalysis of data the
pipeline already produces, no new FastK/k-mer work (the fourth,
`mean_run_length_windows`/`flip_rate`, is described further down —
it's the direct version of what `mean_windowed_cv` only proxies, and it
did not hold up empirically):

- **`partition_consistency`**: for every chromosome with ≥3 haplotype
  copies, splits them into two lineages by whole-chromosome distance
  (reusing `subgenome_report.bipartition_by_distance`), then checks whether
  the *same* partition recurs across every chromosome. A fixed partition
  every time (e.g. always `{HAP1,HAP2}` vs `{HAP3,HAP4}`) is a disomic/allo
  signature; a partition whose "odd one out" rotates between different
  haplotype labels chromosome-to-chromosome is a tetrasomic/auto signature.
  Only defined for species with ≥3 copies and ≥2 splittable chromosomes —
  most of the panel's 2-copy species can't be evaluated on this axis at
  all.
- **`mean_windowed_cv`**: mean coefficient of variation of the raw windowed
  divergence track across every accepted homeolog pair. Low = uniform
  divergence along the whole pair (one clean historical split); high =
  patchy/mosaic (frequent local exchange or introgression).
- **`pair_depth_cv`**: coefficient of variation of `homeolog_pairs.tsv`'s
  `mean_distance` across a species' accepted pairs. Low = every pair
  diverged to about the same depth (one historical event); high = wide
  spread (complex/asynchronous history).
- **`distance_ratio_cv`**: the *actual* cross-chromosome divergence-depth
  consistency test — coefficient of variation of `homeologs/
  ploidy_ancestry_summary.tsv`'s `distance_ratio` (pair distance normalized
  to each pair's own within-chromosome baseline) across a species' accepted
  pairs, distinct from `pair_depth_cv`'s raw-`mean_distance` version. Low =
  every pair diverged to the same relative depth (one clean historical
  event, allo-like); high = wide, inconsistent spread (messier/multi-event
  history, auto-like). This is what actually formalizes the `distance_ratio`
  spread already flagged ad hoc elsewhere in this document (`daLatSqua1`
  20×, `llColAutu1` 37×) into a proper statistic — see results below.

**`singleton_artifact_suspected`**: cross-references the modal partition's
recurring "odd one out" haplotype against `unplaced.tsv` fragment counts.
**This turned out to be a much blunter instrument than intended** — every
>2-haplotype species in the panel structurally shows near-zero unplaced
fragments on `HAP3`/`HAP4`+ by default (the curation-pipeline artifact
documented above), so this check fires on that same expected background
pattern almost everywhere (`daGalBore1`, `ddHypPerf1`, `ddSalTria1`,
`ddHesMatr1`, `drLytSali1` all flag "yes" for exactly this reason) — it is
**not**, on its own, evidence of anything unusual for those species. It
only carries real weight where the asymmetry is extreme *even relative to
that expected baseline*: `ddHypMacu1`'s `HAP1` (570 unplaced fragments vs.
2-3 for siblings, previously confirmed to correlate with inflated
`te_marker_fraction`) is the one case here with independent corroboration.
Treat every other "yes" as "the routine default," and treat a "no" as the
more informative result — `ddLepDrab1` is the only species whose singleton
(`HAP4`) is *not* explained by this pattern, making its near-total
partition consistency (see below) the most credible real-structure finding
in this analysis. The real gold-standard check remains what was done for
`ddHypMacu1` and `drAriEdul1`: verify whether the asymmetry actually
correlates with the divergence signal, not just whether it exists.

**`mean_run_length_windows` / `flip_rate` — the direct run-length test.**
`mean_windowed_cv` (above) is a magnitude-patchiness *proxy* for "does the
divergence stay uniform along the chromosome"; it isn't literally a
measurement of how long a stable local block lasts before switching. This
is the direct version: for each chromosome's whole-chromosome partition,
walk the within-chromosome windowed track and ask, window by window,
whether the local self-vs-cross ordering still agrees with the global
grouping. A real disomic/allo split should hold up almost everywhere (long
runs, low flip rate — subgenomes recombine internally, rarely exchange
with each other); a tetrasomic/auto split should flip locally far more
often (short runs, high flip rate — routine multivalent recombination
continually reshuffling which copies resemble which).

**This did not survive contact with the data as a discriminator.** Across
the 12 species it's computable for (all need ≥3 haplotype copies),
`SchCurv1` and `SchYoun1` — both independently confirmed autopolyploids
(Xie et al. 2026) — land at opposite ends of the observed range (4.12 vs.
5.31 mean run length), so whatever this signal is dominated by at 250kb
window resolution, it isn't primarily inheritance mode; more likely it
tracks per-species noise characteristics (genome size, repeat content,
raw window count — `ddHesMatr1` alone contributes 12,763 of the ~34,000
windows tested panel-wide). Deliberately **not** folded into
`combined_allo_score` for this reason — reported as informational only.
The one genuine standout is `ddHesMatr1` (2.25 windows, flip_rate 0.44 —
roughly double every other species), which is a good match for its
literature description as a "segmental allotetraploid" (mixed
bivalent/quadrivalent meiotic pairing — a real patchwork of disomic- and
tetrasomic-like regions is exactly what unusually frequent local flipping
would look like), but one matching outlier out of twelve is not evidence
the metric works in general — treat it as a single corroborating data
point for `ddHesMatr1` specifically, not validation of the method.

**`distance_ratio_cv` results, across all 24 species with ≥2 accepted
pairs — this one behaved well.** Confirmed allo anchors sit low, as
expected: `drTriRepe1` 0.037 (tightest in the panel — fits its recent,
15–28kya, single-event origin), `daGleHede1` 0.110 (the primary anchor).
`ddLepDrab1` also sits low (0.132) — a *third* independent metric now
pointing the same direction as `partition_consistency` (0.94) and the
routine-artifact check (cleared), strengthening it as the panel's best
allo candidate. `drMyrSpic1` sits moderately high (0.315) in a way that
makes direct mechanistic sense: its confirmed origin is a documented
**stepwise** allohexaploidy ((AB)C, two sequential events at different
times), and pairs from different historical events should diverge to
different relative depths — exactly what elevated CV here means.
`daSenVulg1`'s apparent 0.903 outlier was traced to source (see its
species-table row below): the ancient-divergence signal itself is
uniform, the ratio just has a noisy denominator (116× within-chromosome
heterozygosity variation unrelated to the WGD) — a real, generalizable
caveat about this metric, not evidence about this species specifically.
`lpElePalu1` is the single highest genuine value in the entire panel
(1.302), which fits its
already-documented genus-level agmatoploidy/symploidy (chromosome
fission-fusion) — a genus where chromosome number changes via fission and
fusion rather than clean WGD would produce exactly this kind of wildly
inconsistent pair-to-pair divergence depth, reinforcing the prior read
that `lpElePalu1`'s pairing signal is likely real chromosome-number
biology rather than unresolved WGD. As with the other magnitude-based
metrics, `drLytSali1` lands low (0.083) via the same tetrasomic-
homogenization confound as `mean_windowed_cv`/`pair_depth_cv` — consistent
with, not contradicting, everything already established about that
species, and the reason this metric is also excluded from
`combined_allo_score`.

**Reading agreement and conflict across `te_marker_fraction`,
`partition_consistency`, `distance_ratio_cv`/`pair_depth_cv`, and
run-length/`flip_rate`.** These four measure genuinely different things —
repeat content, lineage identity, event-count consistency, and local
switching rate, respectively — so "conflict" between them is informative,
not just noise, *if* you know which of five patterns you're looking at.
Every one of these has now shown up at least once in this panel:

1. **True agreement (all available metrics point the same way) — highest
   confidence, either direction.** `drTriRepe1` (allo: low `distance_ratio_
   cv` 0.037, `te_frac` 0.169) and `drLytSali1`/`ddEmpNigr1` (auto: low
   `partition_consistency` 0.33–0.40) are the clean cases.
2. **Tetrasomic homogenization: `partition_consistency` reads auto (low),
   every magnitude-based metric reads allo (low CV/tight).** Ongoing
   multivalent recombination equalizes bulk divergence and depth ratios
   genome-wide, but doesn't touch lineage *identity* the way
   `partition_consistency` tests it. `drLytSali1` is the clean example —
   trust `partition_consistency` here, not the CV metrics.
3. **Assembly-quality artifact: `partition_consistency` reads allo (high),
   but the singleton is a known fragmented/incomplete haplotype, not a
   real lineage.** `ddHypMacu1`'s HAP1 (570 unplaced fragments vs. 2–3 for
   siblings) drives its 1.00 consistency; `singleton_artifact_suspected`
   plus a direct `unplaced.tsv` check catches this. Discard the
   `partition_consistency` reading entirely here, don't average it in.
4. **Genuine mosaic/segmental biology: metrics disagree because the
   organism itself is heterogeneous, not because one metric is wrong.**
   `ddHesMatr1` — literature calls it a "segmental allotetraploid" (mixed
   bivalent/quadrivalent meiosis); its run-length flip_rate (0.44, by far
   the panel's highest) and low `te_frac` (0.089) both fit a genome that's
   part-disomic, part-tetrasomic depending on region. Don't force a single
   auto/allo label onto a species where the real answer is "both, in
   different places."
5. **Multi-event history inflates the "messy" metrics without implying
   auto.** `drMyrSpic1`'s confirmed stepwise (AB)C allohexaploid origin
   produces an elevated `distance_ratio_cv` (0.315) and a huge te_frac
   range (0.018–0.971) — high spread here means "more than one historical
   event," not "not allo." `llColAutu1`'s elevated `distance_ratio_cv`
   (0.586) likely reflects the same principle via a different named
   mechanism (Chacon & Renner 2014's "demi-duplication" — repeated fusion
   of mismatched-ploidy gametes is inherently multi-event). Never read
   "high spread" as auto by default — check whether a documented
   multi-event or non-WGD mechanism explains it first.
6. **Metric-construction artifact: the *ratio* is unreliable even though
   the underlying signal is clean.** `daSenVulg1`'s `distance_ratio_cv`
   (0.903) looked like a real conflict with its very uniform raw pair
   distances (`pair_depth_cv` 0.031) — traced to the ratio's denominator
   (ordinary within-chromosome heterozygosity) varying 116× for reasons
   unrelated to the ancient WGD. When `distance_ratio_cv` and
   `pair_depth_cv` disagree sharply, check the raw `own_mean_distance`
   spread in `ploidy_ancestry_summary.tsv` before trusting either number.

**Practical trust ordering, in light of all of this**: `partition_consistency`
first (most mechanistically direct, but always cross-check
`singleton_artifact_suspected`) → `te_marker_fraction` second (direct
repeat-content signal, but check the HAP1/2-vs-HAP3/4+ curation asymmetry
and, for >2-copy species, prefer the lineage-split value over the flat
genome-wide mean) → `distance_ratio_cv`/`pair_depth_cv` third (useful for
spotting multi-event or non-WGD histories, not a direct auto/allo signal
on its own, and check the two against each other before trusting either)
→ run-length/`flip_rate` last (mostly noise, only informative as a
segmental/mosaic-architecture detector). A single remaining disagreement
after this process — like `daSenVulg1`'s low `te_frac` against its
literature lean — is a genuine open question, not something to force an
answer onto.

**Selected results** (full panel in `meta/auto_allo_spectrum.tsv`):

| species | partition_consistency | singleton (artifact?) | run length / flip rate | mean_windowed_cv | pair_depth_cv | read |
|---|---|---|---|---|---|---|
| `drLytSali1` | 0.40 | HAP4 (routine default) | 3.72 / 0.263 | 0.012 | 0.046 | **auto-consistent** — rotating partition matches the strong classical-genetics literature (tetrasomic inheritance, double reduction); the low CVs are the tetrasomic-homogenization confound, not allo evidence |
| `ddEmpNigr1` | 0.38 | n/a | 4.61 / 0.210 | – | – | **auto-consistent** — matches lit call via an independent mechanism from `te_marker_fraction` |
| `SchCurv1` | 0.42 | cleared (median-robust check) | 4.12 / 0.236 | 0.035 | 0.150 | **auto-consistent** — matches confirmed auto call (Xie et al. 2026) |
| `SchYoun1` | 0.33 | n/a | 5.31 / 0.181 | – | – | **auto-consistent** — matches confirmed auto call, despite the highest run-length of the confirmed-auto pair (see caveat above — this metric doesn't discriminate reliably) |
| `daBudDavi1` | 0.26 | n/a | 5.35 / 0.179 | – | – | reads auto-like on partition_consistency, but longest run-length in the panel (allo-like by that axis) — internally inconsistent, corroborating the existing lit-vs-data tension rather than resolving it |
| `ddHypPerf1` | 0.38 | routine default | 4.62 / 0.210 | – | – | reads auto-like — new data point for a species the literature calls genuinely contested |
| `ddHypMacu1` | 1.00 | **HAP1, confirmed real artifact** | 4.12 / 0.235 | 0.021 | – (1 pair) | **not trustworthy** — full consistency is the same known assembly-fragmentation confound that inflated its raw `te_marker_fraction`, not real disomic structure; run-length sits with the confirmed autos, not with a stable-partition read |
| `ddLepDrab1` | 0.94 | HAP4, **not explained by the default pattern** | 4.49 / 0.213 | 0.035 | 0.102 | **best candidate for genuine allo-like structural stability** in the panel on partition_consistency, but its run-length sits *below* `SchYoun1` (confirmed auto) — the two metrics disagree here, treat the allo read as tentative pending literature |
| `ddSalTria1` | 0.63 | routine default (triploid) | 5.48 / 0.165 | 0.033 | 0.019 | ambiguous — triploid, partial consistency, singleton explained by the routine pattern; longest run-length in the panel |
| `ddHesMatr1` | 0.67 | routine default | **2.25 / 0.444 (outlier)** | 0.021 | – (1 pair) | segmental allotetraploid per lit — the run-length outlier is the one place this metric adds real, corroborating signal |
| `drAriEdul1` | 0.24 | HAP4, real signal (not artifact, per prior manual check) | 5.34 / 0.180 | – | 0.028 | ambiguous — long run-length here doesn't cleanly resolve the identity/cytotype question |
| `daGleHede1`/wheat `AB` | n/a (2-copy) | n/a | n/a (needs ≥3 copies) | 0.027 / – | 0.042 / – | metric 1 and 4 don't apply; both already anchors via the other metrics |

**On a single combined number** — the script computes one
(`combined_allo_score`, an average of the three metrics mapped to a common
[0,1] direction), but **it is not reliable as a headline figure** and
shouldn't be quoted in isolation. Worked counter-example: `drLytSali1`
scores 0.40 on `partition_consistency` (correctly reads auto) but 0.78 on
the naive combined average, because the other two metrics get dragged
allo-ward by the exact same tetrasomic-homogenization effect that produces
its low bulk divergence in the first place — averaging dilutes the one
metric that actually saw through the confound. **Use
`partition_consistency` as the primary signal whenever it's defined** (≥3
copies), and read the other two as corroborating/context, not as
equal-weighted inputs to one score.

## Diffuse, genome-wide partitions (species with no accepted homeolog pairs)

`homeolog_pairs.tsv`'s FDR test looks for discrete 1:1 ancestral-duplicate
*pairs* — it can miss a real, weaker, genome-wide structure where no single
pair is individually significant but the chromosome set as a whole still
factors into subgenome-sized groups. `scripts/genome_partition.py` tests
for exactly this: average-linkage agglomerative clustering over the
*complete* all-pairs distance matrix (`homeolog_candidates_ranked.tsv`,
every cross-chromosome pair already tested, no new k-mer work), recording
the partition at every possible k in one merge pass, then testing each k's
separation ratio (mean between-group distance ÷ mean within-group
distance) against a permutation null of ~500 random partitions with the
same group-size profile. Output in `meta/genome_partition.tsv`.

**Validated flagship case: `daInuConz1`** (*Inula conyza*, now
*Pentanema conyzae*). Filed as a clean diploid floor reference — lowest
`te_marker_fraction` in the panel (0.017), 0/16 accepted pairs. Wrong
resolution, not wrong instinct: `genome_partition.py` finds a highly
significant (z=10.5) 8-vs-8 chromosome-number bipartition, independently
confirmed by a clean, **non-overlapping** repeat-density split between the
two groups (1354 vs. 1074 high-copy k-mers/Mb — length alone doesn't
separate them, but density does). Cytology confirms this species is
tetraploid, 2n=4x=32, base number **x=8** — replicated across 7+ independent
karyological studies since 1977 — meaning the discovered 8/8 split matches
the confirmed base chromosome number exactly, found *before* that
literature check and with no way to have been influenced by it. The clade
(*Pentanema conyzae* group, Inuleae) is independently documented as showing
"reticulate patterns... best explained by hybridization and introgression"
with allopolyploidization explicitly proposed for close congeners — so this
looks like a real, allo-leaning founding tetraploidy event, structurally
invisible to every pairwise-FDR test but recoverable at the genome-wide
level. `te_marker_fraction` measured only ordinary same-chromosome-number
heterozygosity here (genuinely low) — it was never testing the actual
subgenome axis for this species at all, exactly as with `daGleHede1`'s
partition-consistency discussion above.

**Full panel, sorted by strength** (see `meta/genome_partition.tsv` for
every significant k per species, not just the best one — several species
show nested structure at multiple k simultaneously):

- **Previously "zero signal," now real structure**: `dcCerAlpi1` (k=2,
  uneven 16/20 split, z=21.7 — the highest-confidence new finding in the
  panel, and a plausible explanation for its previously-unexplained
  moderate `te_frac` of 0.198 despite 0 accepted pairs), `dmRanRepe1` (k=2,
  8/8, z=11.1), `daInuConz1` (above).
- **Extends/refines already-known pairing into a coarser layer**:
  `drLytSali1` (k=5, five groups of 3 — each one is a known accepted pair
  plus one extra chromosome), `ddLepDrab1` (k=8 exactly reproduces the 8
  known pairs at z=11.1, but **k=7 is more significant at z=10.8**, merging
  two known pairs into a quartet — direct support for a nested,
  allo-octaploid-like block structure, persisting down through k=6→2),
  `daLatSqua1` (k=5, groups of 4 rather than known pairs of 2 — same
  block pattern as `ddLepDrab1`), `daSenVulg1` (k=10, all-pairs — extends
  the previously-known 8 accepted pairs to 10).
- **Uneven group sizes, mechanistically consistent with "demi-duplication"
  rather than clean WGD**: `llColAutu1` (k=10, sizes 3,4,4,5,5,6,6,6,6,6 —
  irregular rather than clean pairs, matching the repeated
  mismatched-ploidy-gamete-fusion mechanism already proposed in the
  literature for this genus).
- **Weak, statistically fragile — flag but don't trust yet**: `daSonOler1`
  (z=5.4, notable because this was previously the panel's *cleanest*
  diploid reference), `ddSalTria1` (z=3.2), `ddHesMatr1` (z=2.1). **Real
  caveat**: every k from 2 to N−2 is scanned per species (10–30 values for
  most of the panel), so z>2 alone will turn up somewhere by chance for
  some species even under a true null — the z>10 findings above are robust
  to this, these three are not.
- **Everything else** in the panel reproduces its already-known accepted
  pairs almost exactly (e.g. `drSorDevo1`, `ddSalPent1`, `drTriRepe1`, the
  *Potamogeton* species) — a useful consistency check between the two
  independent methods, not new information on its own.

**Balanced vs. degenerate k=2 splits — a further diagnostic, and a
correction to two species previously filed as diploid.** Not every
significant k=2 partition means the same thing. Checking whether each
chromosome has a *unique reciprocal best match* on the other side (the same
check that validated `daInuConz1`'s 8/8 split) separates real whole-genome
bipartitions from degenerate clustering:

- **Balanced, with clean reciprocal 1:1 matching — real candidate block
  structure**: `daInuConz1` (8/8, z=11.4, 7/8 pairs cleanly reciprocal),
  `dcCerAlpi1` (16/20, z=21.5 — the highest in the panel — 16/20
  chromosomes reciprocally clean, 4 unresolved), and, most strikingly,
  **`dmRanRepe1` (8/8, z=10.7, all 8 pairs perfectly reciprocal, zero
  exceptions — cleaner than `daInuConz1`'s own matching)**. `lpElePalu1`
  (9/10, z=5.8) also qualifies but more weakly, and this genus has a
  documented alternative explanation (fission-fusion/agmatoploidy) that
  could produce block-like clustering without two distinct parental
  origins — treat it with more caution than the other three.
- **Degenerate — one tight leftover pair swallowed into one large blob,
  not a genuine parental split**: every other species with a significant
  k=2 row, including `daGleHede1` (see the worked comparison above) and
  the confirmed-allo `drTriRepe1`, `drSorDevo1`, and all five *Potamogeton*
  species. A degenerate k=2 doesn't mean "not allo" — several of these are
  confirmed or strongly-evidenced allopolyploids by other means — it means
  the genome-wide parental-origin fingerprint (as opposed to the
  chromosome-by-chromosome pairwise relationship) is no longer detectable,
  whether from age, from faster-homogenizing processes like biased gene
  conversion/TE turnover between subgenomes, or from the two parental
  lineages having been less differentiated from each other to begin with.

**`dcCerAlpi1` and `dmRanRepe1` no longer belong in the diploid bucket.**
Both were previously filed as the panel's best candidates for a genuine
diploid/no-WGD calibration anchor (0 accepted pairs, "clean diploid
signature"). The balanced, near-complete reciprocal bipartite structure
found here is the same *quality* of evidence that supports calling
`daInuConz1` a likely allotetraploid — `dmRanRepe1`'s matching is if
anything cleaner (perfect 8/8, no exceptions at all). Both also show
notably *lower* cross-group distances (`dcCerAlpi1` 0.080–0.088,
`dmRanRepe1` 0.070–0.073) than `daInuConz1`'s (0.093–0.099), suggesting —
if real — these may be among the **youngest, least-diverged candidate
allopolyploid signals in the entire panel**: parental genome-of-origin
block structure fully intact, not yet worn down the way `daGleHede1`'s
appears to be. This is a genuine revision to [[project_autoallo_calibration_gap]]
territory: two species we were treating as diploid anchors now look like
they may instead be *cryptic, undescribed* allopolyploids — worsening,
not filling, the panel's diploid-anchor gap until literature can confirm
or rule this out for *Cerastium alpinum* and *Ranunculus repens*
specifically. Neither has a species-level ploidy/origin literature check
yet (both genera have documented polyploid complexes generally, which
doesn't confirm anything for these species specifically) — this is the
single highest-value literature gap this analysis has surfaced.

## The "checkerboard" pattern, and reordering by known pair membership

Several species with complete or near-complete accepted homeolog pairing
(`drTriRepe1`, `drRosSpin1`, `daSenVulg1`, `drSorDevo1`, `llColAutu1`, and
others) show a scattered, checkerboard-like pattern in the default
chromosome-number-ordered heatmap, distinct from the smooth two-block
pattern of the diffuse-partition species above. `scripts/
reorder_heatmap.py` tests whether this is a numbering-order artifact of
already-strong, individually significant pairs (it should resolve into a
clean pattern once reordered by the actual `homeolog_pairs.tsv`
assignment) or something unresolved. Reuses `homeolog_candidates_ranked.tsv`,
no new k-mer work. Output: `results/<species>/homeologs/
homeolog_pairs_reordered_heatmap.png`, run panel-wide.

**For most species, it's exactly a numbering artifact.** `drTriRepe1`,
`drRosSpin1`, `drSorDevo1`, and the majority of the panel collapse cleanly:
every accepted pair becomes one sharp, tight cell immediately off the
(masked) diagonal, with an otherwise flat, uniform background everywhere
else. Chromosome numbers here don't track a clean size rank or any
subgenome-aware scheme (checked directly on `drTriRepe1`: `chr01` is
55.9Mb, smaller than 14 of its 15 siblings, breaking any simple
size-order story) — pair partners just end up scattered at arbitrary
numeric distance from each other, which is what produces the checkerboard
look in raw chromosome-number order. Once reordered, there's nothing left
to explain.

**Two species don't fully collapse, and that's real signal, not a
failure of the reordering.**

- **`llColAutu1`** keeps visible extra structure even after reordering by
  its 24 pairs — several groups of 4–8 chromosomes (spanning multiple
  pairs) stay noticeably tighter with each other than with the rest of the
  genome (e.g. `chr18/28/25/26/10/17/30/35`), and `chr01`/`chr02` show
  unusually elevated distance reaching into otherwise-unrelated regions.
  This directly corroborates the uneven, non-pairwise group sizes
  `genome_partition.py` already found for this species (3,4,4,5,5,6,6,6,6,6)
  and the literature's proposed "demi-duplication" mechanism (repeated
  fusion of mismatched-ploidy gametes) — a single clean WGD should produce
  the flat-background-plus-pairs pattern every other species shows;
  `llColAutu1`'s persistent extra structure is consistent with a messier,
  multi-event history layered on top of the primary pairing.
- **`ddLepDrab1`** shows a clear, visible tighter sub-block spanning
  `chr05,chr16,chr03,chr06,chr10,chr12` — three of its eight pairs
  clustering with each other more than with the rest. This is a direct
  visual confirmation of the nested quartet structure `genome_partition.py`
  found underneath the primary 8 pairs (k=7 merges exactly `chr05,chr16`
  with `chr10,chr12` into one group, more significant than the primary
  k=8 pairing itself) — corroborating evidence for the allo-octaploid-block
  hypothesis from two independent methods now.
- **`daSenVulg1`**: the 4 chromosomes left unpaired by the strict FDR test
  (`chr03,chr04,chr11,chr19`) are visibly tighter with each other than with
  the paired set — consistent with the `genome_partition.py` k=10 result
  that extended this species' 8 known pairs to 10 by picking up 2 more
  that didn't individually clear significance.

**On making this the default view**: not a full replacement, for two
concrete reasons. First, it's undefined for any species with zero accepted
pairs — exactly the diffuse-partition species in the section above, which
still need the general-purpose clustering-based heatmap since there's no
pair table to sort by. Second, it deliberately drops to chromosome-number
resolution (aggregating away haplotype copies), which is coarser than the
existing unit-level `whole_chrom_distance_heatmap_contrast.png` for ≥3-copy
species — that's exactly the resolution that caught `drLytSali1`'s
rotating-singleton pattern, which a chromosome-number-only view can't see.
Recommended instead: generate both. The pair-reordered view is the clearer
default *read* whenever accepted pairs exist (which is most of the panel),
but the general clustering-based heatmap remains the one actually capable
of *discovering* structure the FDR test missed, as demonstrated repeatedly
in this document.

## Species summary

Every species run through the pipeline so far. `copies` = distinct
haplotype-copy counts seen across chromosomes; `accepted pairs` = ancient
homeolog pairs surviving FDR; `te_frac` = mean `te_marker_fraction` where
`te-markers` has been run. Calls marked *(lit)* come from literature
cross-referencing (`meta/literature_ploidy.tsv`); the rest are data-only
reads from this panel and haven't been checked against literature —
treat those as provisional.

| species | copies | accepted pairs | te_frac | read |
|---|---|---|---|---|
| `daGleHede1` (*Glechoma hederacea*) | 2 | 18/18 | 0.223 | **allo** *(lit)* — primary calibration anchor; clean, fully-resolved, uniform subgenome split |
| `daBudDavi1` (*Buddleja davidii*) | 4 | 0/19 | 0.167 | **allo** *(lit, phylogenomic)*, but heatmap alone reads auto — see caveat above; one localized chr08 block found in windowed heatmap, candidate incomplete homogenization. Yang et al. 2023 documents a specific cytonuclear conflict for this species (nests with diploids in the plastid tree, polyploids in ASTRAL/nrDNA) — plausibly *why* our raw-sequence read is auto-like despite the allo call; reframed from unexplained tension to corroborated complexity |
| `daGalBore1` (*Galium boreale*) | 4 | 4/11 | 0.129 | unresolved *(lit: ploidy confirmed — part of a documented 2n=22/44/66 polyploid complex — origin not addressed)* |
| `daEupConf1` (*Euphrasia confusa*) | 1 | 22/22 | – | **allo** *(lit, pre-existing annotation)*; only 1 usable haplotype but strong ancient-pairing signal within it |
| `daInuConz1` (*Inula conyza*, now *Pentanema conyzae*) | 2 | 0/16 | 0.017 | **not diploid** — confirmed tetraploid, 2n=4x=32, x=8 (7+ karyological studies since 1977); `te_frac`/0 accepted pairs only reflect ordinary heterozygosity, the real subgenome signal is a diffuse 8-vs-8 chromosome-number split (`genome_partition.py`, z=10.5) matching x=8 exactly, confirmed by non-overlapping repeat density between the two groups — allo-leaning per clade-level hybridization/reticulation evidence (Gutiérrez-Larruscain et al. 2018/2019), see dedicated section above |
| `daLatClan1` | 2 | 18/21 | 0.254 | near-fully-resolved ancient pairing; hexaploid genus (2n=42, PHYA gene duplication confirmed in this species too), but auto/allo origin not addressed by any source found — moderate te_frac, similar magnitude to several confirmed allos, but on its own not decisive given the literature gap. 2026-09 Darwin batch. |
| `daLatSqua1` (*Lathraea squamaria*) | 2 | 16/18 | 0.320 | strong ancient-pairing + high te_frac — allo-leaning; wide `distance_ratio` spread (20×; `distance_ratio_cv` = **1.01, the widest in the panel** except `lpElePalu1`) is now explicit auto-vs-allo evidence, not just "candidate asynchronous resolution" — this level of pair-to-pair inconsistency is more consistent with a messier/multi-event history than one clean allopolyploidization. Hexaploid (2n=42) and PHYA gene-duplication confirmed independently across nearly a century of sources, but origin (auto vs allo) genuinely unaddressed in the literature found |
| `daPilAura1` | 2 | 18/18 (full) | 0.207 | fully-resolved ancient pairing, all 18 chromosomes paired into 9 pairs, tight/consistent divergence (0.031–0.044) — consistent single-age WGD signature, and te_frac (0.207) sits close to `daGleHede1`'s confirmed-allo value (0.223). Tetraploid (2n=36), facultatively apomictic, extensively studied as a hybridizing *parent* of further crosses (e.g. *P. rubra*), but its own tetraploidy's origin is never directly addressed. 2026-09 Darwin batch. |
| `daSenVulg1` | 2 | 16/20 | 0.041 | near-fully-resolved ancient pairing, uniform divergence (~0.06–0.07 across all 8 pairs, raw `pair_depth_cv` = 0.031, tight) — consistent single-age WGD signature by the raw distance. **Traced the apparent `distance_ratio_cv` = 0.90 "conflict" to its source: it's not messy ancient divergence, it's the *normalization denominator*** — each chromosome's own within-chromosome (HAP1-vs-HAP2) baseline distance varies 116× across the genome (0.000016–0.001859), for reasons unrelated to the ancient WGD (ordinary heterozygosity/recombination-rate variation). The ancient signal itself stays uniform; only the ratio to a noisy denominator looks inconsistent. So this isn't really a metric conflict — `pair_depth_cv` and the ancient-pairing uniformity agree with each other and with an allo-consistent read; `distance_ratio_cv` is simply the wrong tool here, and `te_frac` (0.041, low) is the one genuine outstanding disagreement, on a completely different axis (repeat content, not divergence depth). Genuinely contested historically (auto-from-*S.-vernalis*, Kadereit 1984, later refuted by isozyme evidence), but more recent molecular work (Chapman & Abbott 2010's RAY2b homeolog finding; Kim et al. 2008, *Science*) leans **allo** — the low te_frac is the one thing still at odds with that lean. 2026-09 Darwin batch. |
| `daSonOler1` (*Sonchus oleraceus*) | 2 | 0/16 | 0.039 | diploid-looking, low te signal |
| `dcCerAlpi1` (*Cerastium alpinum*) | 2 | 0/36 | 0.198 *(stale, prior assembly)* | **no longer treated as diploid** — moderate te signal despite no accepted ancient pairs is now explained: `genome_partition.py` finds the panel's strongest (z=21.5) 16-vs-20 diffuse bipartition, 16/20 chromosomes reciprocally clean, unusually low cross-group distance (0.080–0.088) suggesting a young/undegraded candidate allopolyploid — see dedicated section above. No species-level literature check exists yet; reconfirmed 0/36 on 2026-09 Darwin reassembly, see batch note above |
| `dcHonPepl1` (*Honckenya peploides*) | 2 | 0/34 | 0.028 | mostly diploid-looking; some tetraploid-like clusters noted visually earlier, not yet reconciled with this low te value |
| `ddAraThal4` (*Arabidopsis thaliana*) | 1 | 0/5 | – | **true diploid** *(lit, confirmed)* — panel's diploid anchor |
| `ddEmpNigr1` (*Empetrum nigrum*) | 4 | 2/13 | 0.239 | **auto** *(lit)*, but te_frac sits at/above the confirmed-allo range — unresolved tension, flagged in `meta/literature_ploidy.tsv` |
| `ddHesMatr1` (*Hesperis matronalis*) | 4 (tetraploid) | 2/6 | 0.089 | sparse pairing — only `chr01↔chr02` accepted out of 6 chromosomes; needed three resubmissions (16GB→28GB→40GB) to clear `TERM_MEMLIMIT` during k-mer table building. Low te_frac (0.089) is consistent with mostly-undiverged copies. **Literature directly describes this species as a "segmental allotetraploid"** (Francis et al. 2009: both bivalent and quadrivalent meiotic pairing observed — partial subgenome homology, a genuine third category between clean auto and clean allo) — maps remarkably well onto our own sparse-pairing + low-te_frac result (mostly no divergence signal, one real pair), the best-supported nuanced case found in this literature pass. 2026-09 Darwin batch. |
| `ddHypMacu1` (*Hypericum maculatum*) | 4 | 2/8 | 0.291 (corrected; raw 0.445 was assembly-quality-inflated) *(stale, prior assembly)* | **auto** *(lit)*; corrected te_frac still elevated vs. confirmed autos — same tension as above; reconfirmed 2/8 (chr01↔chr03) on 2026-09 Darwin reassembly, see batch note above |
| `ddHypPerf1` (*H. perforatum*) | 4 | 0/8 | 0.408 | contested in the literature (both auto and allo hypotheses actively cited); high te_frac leans allo but unconfirmed |
| `ddLepDrab1` (*Lepidium draba*) | 4 | 16/16 (full) | **0.637 — highest in the panel, above any confirmed allo anchor** | cleanest allotetraploid-like structure seen in the panel — all 16 chromosomes resolve into 8 pairs; corroborated by **three** independent metrics now (`partition_consistency` 0.94, not explained by the routine HAP1/2-vs-HAP3/4 artifact; `distance_ratio_cv` 0.13, one of the tightest in the panel; and now `te_frac` 0.637, the highest of any species in this analysis) — the panel's strongest multi-metric case outside the confirmed anchors, though the run-length metric still disagrees (see spectrum section above). Literature corroborates at the tribe level: Lepidieae (containing *L. draba*) is specifically flagged as falling within the "hybrid or highly polyploid classes," attributed to rampant hybridization that "often precedes allopolyploidization" — medium confidence, not yet a species-specific confirmed call, but this is now the single best-evidenced *candidate* allo species in the entire panel without a formal literature confirmation. 2026-09 Darwin batch. |
| `ddMalSylv1` (*Malva sylvestris*) | 2 | 0/21 | 0.135 | diploid *(lit)* — "a dead diploid" per visual inspection |
| `ddPopNigr1` (*Populus nigra*) | 1 | 0/19 | – | diploid-looking (alt haplotype too fragmentary) |
| `ddSalCine1` (*Salix cinerea*) | 1 | 0/18 | – | diploid-looking (alt haplotype too fragmentary) |
| `ddSalPent1` (*Salix pentandra*) | 2 | 38/38 | 0.200 | strong, complete ancient-pairing signal — allo-leaning. Gulyaev et al. 2022 gives real (if hedged, "might have arisen") phylogenomic evidence: the polyploid group containing *S. pentandra* shows nuclear/chloroplast tree conflict plus admixture-analysis support for a hybrid origin between the *Salix* and *Vetrix* clades — medium confidence, consistent with our complete pairing signal |
| `ddSalTria1` (*Salix triandra*) | 3 (triploid) | 4/19 | 0.247 | sparse pairing — most chromosome copies look interchangeable rather than distinctly diverged, consistent with a fairly homogeneous triploid rather than three diverged subgenomes; te_frac (0.247) sits moderate, close to several confirmed allos, but a wide range (0.055–0.816) suggests very uneven divergence across pairs — hard to interpret cleanly given the ploidy-identity conflict below. **Only literature found (Blackburn & Harrison 1924) reports this species as diploid** (haploid n=19 or 22) — a direct conflict with our assembly's triploid structure, with no modern source addressing the discrepancy; worth treating the triploid call itself with some caution pending a second look at header/manifest parsing for this species. 2026-09 Darwin batch. |
| `dmRanRepe1` (*Ranunculus repens*) | 2 | 0/16 | 0.066 | **no longer treated as diploid** — no chromosome pair shows significant divergence signal at strict FDR, but `genome_partition.py` finds a *perfect* 8-vs-8 reciprocal bipartite match (z=10.7, zero exceptions, cleaner than `daInuConz1`'s own matching) with the lowest cross-group distance of any candidate block structure in the panel (0.070–0.073) — the low te_frac now reads as consistent with a genuinely *young* allopolyploid rather than no polyploidy at all. Currently the strongest "young cryptic allopolyploid" candidate in the whole panel; no species-level literature check exists yet (genus *Ranunculus* has documented polyploid complexes generally). 2026-09 Darwin batch. |
| `drAriEdul1` (*Aria edulis*) | 4 | 16/17 | 0.219/0.209 (no inflation found) | ambiguous per *A. edulis* diploid literature; allo-consistent if sample is actually the tetraploid relative *A. wyensis* — Green 2024's own hypothesis is that *A. wyensis* is an allotetraploid derived from *A. edulis* × a member of the *A. porrigentiformis* group. Misidentification under "*S. aria*"-type names is well documented (Pellicer et al. 2012: 14 specimens collected as *S. aria* were actually triploid), and the *Aria* complex contains *both* confirmed allo- and autotetraploid species (Dickinson 2018) — origin still hinges on resolving exact species ID first. `chr01` lineage split (4.66×) passes both artifact checks and manual digging — `HAP4` (a clean assembly) sits 2.2–4.7× further from `HAP1`/`HAP2`/`HAP3` than they do from each other, with a specific ~13Mb breakpoint; strongest surviving real-signal candidate in the panel, unconfirmed |
| `drIngLaur1` (*Inga laurina*) | 1 | 0/13 | – | diploid-looking |
| `drLytSali1` (*Lythrum salicaria*) | 4 | 10/15 | 0.307 *(stale, prior assembly)* | **auto** *(lit)*, but high te_frac and substantial ancient-pairing signal both lean allo — worth revisiting given the fusion-mimics-dominance caveat above; reconfirmed 10/15 on 2026-09 Darwin reassembly, see batch note above |
| `drMyrSpic1` (*Myriophyllum spicatum*) | 2 | 16/21 | 0.359 (range 0.018–0.971, huge spread) | near-fully-resolved ancient pairing (7 tight pairs, distance 0.057–0.061, + 1 much weaker pair at 0.128) — bimodal divergence hints at two duplication events or a recent rearrangement on an older WGD. Elevated te_frac and its huge pair-to-pair range fit the confirmed origin well. **Confirmed allohexaploid (AABBCC) via a 2026 chromosome-scale genome** (Wang et al. 2026, *Plant Journal*), with a documented stepwise (AB)C origin — maps remarkably well onto our own bimodal pairing pattern (7 tight pairs ≈ the older A-B split, 1 weaker pair ≈ the later, more divergent C lineage) and the wide te_frac spread (pairs from different historical events diverging to different repeat-content depths). One of the best-evidenced allo calls found in this whole literature pass, on par with `daGleHede1`/wheat. 2026-09 Darwin batch. |
| `drMyrVert1` (*Myriophyllum verticillatum*) | 2 | 14/14 (full) | **0.521 — second-highest in this batch** | fully-resolved ancient pairing, all 14 chromosomes paired, and now a strikingly high te_frac. Unlike its close relative *M. spicatum* (confirmed allo, see above), this species is used almost exclusively as an outgroup in *M. spicatum* studies and its own origin is never directly addressed in any source found — even its tetraploid classification (2n=28) has some historical dispute. This high te_frac is a genuinely new finding suggesting *M. verticillatum* may itself be allopolyploid like its better-studied relative, worth flagging for a dedicated literature check. 2026-09 Darwin batch. |
| `drRosSpin1` (*Rosa spinosissima*) | 2 | 14/14 | 0.122 | complete ancient-pairing signal, moderate te_frac — allo-leaning per most (older) sources, **but a 2026 chromosome-level genome paper reports this species as diploid (2n=14)**, directly contradicting the majority tetraploid consensus (2n=4x=28) that our own structure (14 loci × 2 copies, all fully paired into 7 ancient-duplicate groups) is actually more consistent with. A genuine, recent, unresolved contradiction similar in kind to `drAriEdul1`'s identity issue — worth checking whether the Darwin ToL individual matches the 2026 genome paper's sample before trusting either the ploidy level or the allo call here |
| `drSorDevo1` (*Sorbus devoniensis*) | 2 | 34/34 (full) | 0.153 | fully-resolved ancient pairing, all 34 chromosomes paired into 17 pairs — strongest complete-pairing signal in the panel. te_frac is moderate rather than dramatic, but the structural evidence alone is already very strong. **High-confidence allotetraploid**, extensively confirmed: hybridogenic origin from *S. torminalis* (sexual diploid, confirmed maternal parent via chloroplast/microsatellite data across many studies) × subgenus *Aria*/*S. aria* s.l. (paternal), via a triploid-bridge mechanism; peroxidase isozyme banding directly combines patterns from both parents (Proctor et al. 1989). One of the best-evidenced allo calls in the entire panel, matching our own strongest complete-pairing result extremely well. 2026-09 Darwin batch. |
| `drTriDubi3` (*Trifolium dubium*) | 1 | 14/15 | – | strong ancient-pairing despite single usable haplotype |
| `drTriRepe1` (*T. repens*, white clover) | 2 | 16/16 | 0.169 | **allo** *(lit, high confidence — upgraded)* — extensively confirmed allotetraploid (2n=4x=32) across ≥9 independent studies (GISH, ITS/cpDNA, chromosome-level genome assembly), diploid progenitors *T. occidentale* + *T. pallescens*; now on par with `daGleHede1`/wheat as a calibration-quality anchor. Tight, uniform `distance_ratio` (1.1× spread) is consistent with one clean, synchronized event |
| `drUrtDioi1` (*Urtica dioica*) | 1 | 0/13 | – | diploid-looking |
| `laPotLuce1`/`laPotNata1` (*P. lucens*, *P. natans*) | 2 | ~24–26 (near-complete) | 0.25–0.36 | **"probably allotetraploid"** — the only two of the five *Potamogeton* species with a direct, species-specific origin statement (Wang et al. 2007), medium confidence; both are also well-documented parents of named natural hybrids |
| `laPotPerf1` (*P. perfoliatus*) | 2 | ~24–26 | 0.25–0.36 | "presumed allopolyploid... definitive evidence needed" (Ganie et al. 2020) — hedged by the source itself, low-medium confidence; the most-hybridized of the five (parent to ≥3 named hybrids) |
| `laPotNodo1` (*P. nodosus*) | 2 | ~24–26 | 0.25–0.36 | tetraploid (2n=52) confirmed, but no species-specific origin statement found — only indirect support via its role as a hybridization partner (P.×schreberi) |
| `laPotCris1` (*P. crispus*) | 2 | ~24–26 | 0.25–0.36 | the most cytologically unsettled of the five — reported chromosome number varies wildly across sources (2n=50/52/78/84); no direct origin statement, treat its pairing result with extra caution given the unstable base ploidy count |
| `llColAutu1` (*Colchicum autumnale*) | 2 | 48/51 | 0.050 | near-complete ancient pairing but low te_frac — inconsistent pair, worth a second look; wide `distance_ratio` spread in the panel (37× max, `distance_ratio_cv` = 0.59) — now explicit evidence of an inconsistent, non-single-event history, not just "worth a second look." The most authoritative source found (Chacon & Renner 2014) states directly that whether past allo- or autopolyploidy explains *Colchicum*'s chromosome-number lability "remains an open question... no experimental crosses" — but does infer frequent **"demi-duplication"** (fusion of gametes of different ploidy) as the dominant mechanism, a genuine third pathway distinct from classic auto/allo. Demi-duplication plausibly explains our own odd signal directly: strong copy-number/structural pairing (48/51) without necessarily merging two long-diverged genomes (low te_frac), and the high `distance_ratio_cv` fits a mechanism that repeatedly fuses mismatched-ploidy gametes rather than one clean WGD |
| `lpElePalu1` (*Eleocharis palustris*) | 2 | 18/19 | 0.083 *(stale, prior assembly)* | holocentric chromosomes, documented agmatoploidy/symploidy (fission/fusion) in the genus — many pairs sit at `distance_ratio` ≈1, plausibly real biology rather than unresolved WGD signal. **`distance_ratio_cv` = 1.30, the highest in the entire panel** — a strong, mechanistically apt confirmation: a genus where chromosome number changes via fission/fusion rather than clean WGD should produce exactly this kind of wildly inconsistent pair-to-pair divergence depth, reinforcing that this species' pairing signal is real chromosome-number biology, not an unresolved/messy WGD; reconfirmed 18/19 on 2026-09 Darwin reassembly, see batch note above |
| `lpTriTurg1_A` / `_B` (wheat, A/B subgenomes) | 2 | 0 | 0.012–0.015 (within-subgenome) | **allo** *(confirmed)* — structurally-guaranteed low within-subgenome baseline |
| `lpTriTurg1_AB` (wheat, A×B cross-subgenome) | – | – | ~0.42 | **allo** *(confirmed)* — cross-subgenome anchor, the panel's other primary calibration point alongside `daGleHede1` |
| `SchCurv1` (*Schizothorax curvilabiatus*) | 2, 4 | 1/25 (real) | 0.209 genome-wide, but 0.11 within-lineage vs 0.59 cross-lineage at `chr19` | **auto** *(confirmed, Xie et al. 2026)* — 1 ancestral fusion (`chr19+22`), sequence-level ancient signal mostly lost; genome-wide te_frac mean looks allo-like but that's the wrong resolution — see lineage-split analysis above |
| `SchYoun1` (*Schizopygopsis younghusbandi*) | 2, 4 | 0/25 (1 near-miss) | 0.223 genome-wide | **auto** *(confirmed, Xie et al. 2026)* — 5 fusions, clean `M1`+`P1` (unfused) vs `M2`+`P2` (fused) subgenome-like split; fused scaffolds not yet re-processed for their own lineage-split te_frac (see next steps) |

---

*ploidyspec: k-mer-based subgenome/ploidy divergence profiling, no reference
or alignment required. Full column-by-column reference in
[`OUTPUTS.md`](OUTPUTS.md).*
