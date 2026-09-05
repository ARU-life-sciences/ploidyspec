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
footprint rather than a fluke. **`te-markers` has not been (re-)run on
any of these 15's new assemblies** — it's a separate, expensive stage not
included in `all`. The `te_frac` values already in the table for the 4
re-run species are carried over from their *previous* (pre-Darwin-batch)
assembly and are marked stale below; treat them as historical context only,
not as measurements of the new assembly.

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
| `daBudDavi1` (*Buddleja davidii*) | 4 | 0/19 | 0.167 | **allo** *(lit, phylogenomic)*, but heatmap alone reads auto — see caveat above; one localized chr08 block found in windowed heatmap, candidate incomplete homogenization |
| `daGalBore1` (*Galium boreale*) | 4 | 4/11 | 0.129 | unresolved *(lit: ploidy confirmed, origin not addressed)* |
| `daEupConf1` (*Euphrasia confusa*) | 1 | 22/22 | – | **allo** *(lit, pre-existing annotation)*; only 1 usable haplotype but strong ancient-pairing signal within it |
| `daInuConz1` (*Inula conycafé*) | 2 | 0/16 | 0.017 | diploid-looking, low te signal; not independently literature-checked |
| `daLatClan1` | 2 | 18/21 | – *(not yet run)* | near-fully-resolved ancient pairing; data-only, not literature-checked. 2026-09 Darwin batch. |
| `daLatSqua1` (*Lathraea squamaria*) | 2 | 16/18 | 0.320 | strong ancient-pairing + high te_frac — allo-leaning; wide `distance_ratio` spread (20×), candidate asynchronous resolution, not yet artifact-checked pair-by-pair |
| `daPilAura1` | 2 | 18/18 (full) | – *(not yet run)* | fully-resolved ancient pairing, all 18 chromosomes paired into 9 pairs, tight/consistent divergence (0.031–0.044) — consistent single-age WGD signature; data-only, not literature-checked. 2026-09 Darwin batch. |
| `daSenVulg1` | 2 | 16/20 | – *(not yet run)* | near-fully-resolved ancient pairing, uniform divergence (~0.06–0.07 across all 8 pairs) — consistent single-age WGD signature; data-only, not literature-checked. 2026-09 Darwin batch. |
| `daSonOler1` (*Sonchus oleraceus*) | 2 | 0/16 | 0.039 | diploid-looking, low te signal |
| `dcCerAlpi1` (*Cerastium alpinum*) | 2 | 0/36 | 0.198 *(stale, prior assembly)* | moderate te signal despite no accepted ancient pairs — worth a closer look; reconfirmed 0/36 on 2026-09 Darwin reassembly, see batch note above |
| `dcHonPepl1` (*Honckenya peploides*) | 2 | 0/34 | 0.028 | mostly diploid-looking; some tetraploid-like clusters noted visually earlier, not yet reconciled with this low te value |
| `ddAraThal4` (*Arabidopsis thaliana*) | 1 | 0/5 | – | **true diploid** *(lit, confirmed)* — panel's diploid anchor |
| `ddEmpNigr1` (*Empetrum nigrum*) | 4 | 2/13 | 0.239 | **auto** *(lit)*, but te_frac sits at/above the confirmed-allo range — unresolved tension, flagged in `meta/literature_ploidy.tsv` |
| `ddHesMatr1` | 4 (tetraploid) | 2/6 | – *(not yet run)* | sparse pairing — only `chr01↔chr02` accepted out of 6 chromosomes; needed three resubmissions (16GB→28GB→40GB) to clear `TERM_MEMLIMIT` during k-mer table building, peak usage tracked each ceiling closely, possibly a genuine high-memory-footprint species; data-only, not literature-checked. 2026-09 Darwin batch. |
| `ddHypMacu1` (*Hypericum maculatum*) | 4 | 2/8 | 0.291 (corrected; raw 0.445 was assembly-quality-inflated) *(stale, prior assembly)* | **auto** *(lit)*; corrected te_frac still elevated vs. confirmed autos — same tension as above; reconfirmed 2/8 (chr01↔chr03) on 2026-09 Darwin reassembly, see batch note above |
| `ddHypPerf1` (*H. perforatum*) | 4 | 0/8 | 0.408 | contested in the literature (both auto and allo hypotheses actively cited); high te_frac leans allo but unconfirmed |
| `ddLepDrab1` | 4 | 16/16 (full) | – *(not yet run)* | cleanest allotetraploid-like structure seen in the panel — all 16 chromosomes resolve into 8 pairs; data-only, not literature-checked. 2026-09 Darwin batch. |
| `ddMalSylv1` (*Malus sylvestris*) | 2 | 0/21 | 0.135 | diploid *(lit)* — "a dead diploid" per visual inspection |
| `ddPopNigr1` (*Populus nigra*) | 1 | 0/19 | – | diploid-looking (alt haplotype too fragmentary) |
| `ddSalCine1` (*Salix cinerea*) | 1 | 0/18 | – | diploid-looking (alt haplotype too fragmentary) |
| `ddSalPent1` (*Salix pentandra*) | 2 | 38/38 | 0.200 | strong, complete ancient-pairing signal — allo-leaning, not literature-checked |
| `ddSalTria1` | 3 (triploid) | 4/19 | – *(not yet run)* | sparse pairing — most chromosome copies look interchangeable rather than distinctly diverged, consistent with a fairly homogeneous triploid rather than three diverged subgenomes; data-only, not literature-checked. 2026-09 Darwin batch. |
| `dmRanRepe1` | 2 | 0/16 | – *(not yet run)* | clean, complete no-homeolog result — no chromosome pair shows significant divergence signal, same pattern as `dcCerAlpi1`; candidate diploid/no-WGD data point, see batch note above. 2026-09 Darwin batch. |
| `drAriEdul1` (*Aria edulis*) | 4 | 16/17 | 0.219/0.209 (no inflation found) | ambiguous per *A. edulis* diploid literature; allo-consistent if sample is actually the tetraploid relative *A. wyensis*. `chr01` lineage split (4.66×) passes both artifact checks and manual digging — `HAP4` (a clean assembly) sits 2.2–4.7× further from `HAP1`/`HAP2`/`HAP3` than they do from each other, with a specific ~13Mb breakpoint; strongest surviving real-signal candidate in the panel, unconfirmed |
| `drIngLaur1` (*Inga laurina*) | 1 | 0/13 | – | diploid-looking |
| `drLytSali1` (*Lythrum salicaria*) | 4 | 10/15 | 0.307 *(stale, prior assembly)* | **auto** *(lit)*, but high te_frac and substantial ancient-pairing signal both lean allo — worth revisiting given the fusion-mimics-dominance caveat above; reconfirmed 10/15 on 2026-09 Darwin reassembly, see batch note above |
| `drMyrSpic1` | 2 | 16/21 | – *(not yet run)* | near-fully-resolved ancient pairing (7 tight pairs, distance 0.057–0.061, + 1 much weaker pair at 0.128) — bimodal divergence hints at two duplication events or a recent rearrangement on an older WGD; data-only, not literature-checked. 2026-09 Darwin batch. |
| `drMyrVert1` | 2 | 14/14 (full) | – *(not yet run)* | fully-resolved ancient pairing, all 14 chromosomes paired; data-only, not literature-checked. 2026-09 Darwin batch. |
| `drRosSpin1` (*Rosa spinosissima*) | 2 | 14/14 | 0.122 | complete ancient-pairing signal, moderate te_frac — allo-leaning, not literature-checked |
| `drSorDevo1` | 2 | 34/34 (full) | – *(not yet run)* | fully-resolved ancient pairing, all 34 chromosomes paired into 17 pairs — strongest complete-pairing signal in the panel; data-only, not literature-checked. 2026-09 Darwin batch. |
| `drTriDubi3` (*Trifolium dubium*) | 1 | 14/15 | – | strong ancient-pairing despite single usable haplotype |
| `drTriRepe1` (*T. repens*, white clover) | 2 | 16/16 | 0.169 | **suspected allo** *(lit)* — tight, uniform `distance_ratio` (1.1× spread), consistent with one clean, synchronized event like `daGleHede1` |
| `drUrtDioi1` (*Urtica dioica*) | 1 | 0/13 | – | diploid-looking |
| `laPotCris1/Luce1/Nata1/Nodo1/Perf1` (*Potamogeton*, 5 spp.) | 2 | 24–26 (near-complete) | 0.25–0.36 | allo-leaning across the genus (high te_frac, near-complete ancient pairing); heterogeneous `distance_ratio` within each species (1.6–5×) — candidates for mosaic resolution, matches earlier visual read of mixed diploid/tetraploid-like blocks in the heatmaps |
| `llColAutu1` (*Colchicum autumnale*) | 2 | 48/51 | 0.050 | near-complete ancient pairing but low te_frac — inconsistent pair, worth a second look; widest `distance_ratio` spread in the panel (37×), strong asynchronous-resolution candidate |
| `lpElePalu1` (*Eleocharis palustris*) | 2 | 18/19 | 0.083 *(stale, prior assembly)* | holocentric chromosomes, documented agmatoploidy/symploidy (fission/fusion) in the genus — many pairs sit at `distance_ratio` ≈1, plausibly real biology rather than unresolved WGD signal; see prior deep-dive; reconfirmed 18/19 on 2026-09 Darwin reassembly, see batch note above |
| `lpTriTurg1_A` / `_B` (wheat, A/B subgenomes) | 2 | 0 | 0.012–0.015 (within-subgenome) | **allo** *(confirmed)* — structurally-guaranteed low within-subgenome baseline |
| `lpTriTurg1_AB` (wheat, A×B cross-subgenome) | – | – | ~0.42 | **allo** *(confirmed)* — cross-subgenome anchor, the panel's other primary calibration point alongside `daGleHede1` |
| `SchCurv1` (*Schizothorax curvilabiatus*) | 2, 4 | 1/25 (real) | 0.209 genome-wide, but 0.11 within-lineage vs 0.59 cross-lineage at `chr19` | **auto** *(confirmed, Xie et al. 2026)* — 1 ancestral fusion (`chr19+22`), sequence-level ancient signal mostly lost; genome-wide te_frac mean looks allo-like but that's the wrong resolution — see lineage-split analysis above |
| `SchYoun1` (*Schizopygopsis younghusbandi*) | 2, 4 | 0/25 (1 near-miss) | 0.223 genome-wide | **auto** *(confirmed, Xie et al. 2026)* — 5 fusions, clean `M1`+`P1` (unfused) vs `M2`+`P2` (fused) subgenome-like split; fused scaffolds not yet re-processed for their own lineage-split te_frac (see next steps) |

---

*ploidyspec: k-mer-based subgenome/ploidy divergence profiling, no reference
or alignment required. Full column-by-column reference in
[`OUTPUTS.md`](OUTPUTS.md).*
