# drLytSali1 — *Lythrum salicaria* (purple loosestrife)

<!--
SKELETON — same five sections for every species in the panel:

1. Species overview: ToL ID, scientific/common name, assembly source and
   stage, basic structural stats (chromosome count, haplotype-copy count).
   Just the facts needed to orient a reader who's never seen this species
   before -- no interpretation yet.

2. Literature context: what's published about this species' ploidy level
   and origin, independent of anything this pipeline found. Cite sources.
   State confidence honestly, including "nothing found" or "contested."
   Pull from meta/literature_ploidy.tsv -- this section should never
   contain a claim that isn't traceable to a citation there.

3. Direct data outputs: what the pipeline actually measured, reported
   plainly, metric by metric, before any synthesis. If a metric doesn't
   apply (e.g. partition_consistency needs >=3 copies) say so and why,
   don't omit it silently. This section is a report of numbers, not an
   argument.

4. Process inference: synthesize (2) and (3) into a mechanistic story --
   what evolutionary process(es) actually produced the numbers in (3).
   Name which of INTERPRETATION.md's "reading agreement and conflict"
   patterns this species falls into (true agreement / tetrasomic
   homogenization / assembly artifact / mosaic biology / multi-event
   history / metric-construction noise), and say why. This is where
   apparent contradictions between metrics get resolved or flagged as
   still-open.

5. Auto/allo/in-between statement: one paragraph, one clear call --
   allo / auto / genuine third category (name it) / unresolved -- with a
   confidence level and the single strongest piece of evidence behind it.
   End with what would change this call (a literature check, a specific
   further analysis) if it's not fully settled.

Delete this comment block once a species' sections are actually filled in;
keep the section headers and ordering identical across every species so
they stay easy to compare side by side.
-->

## 1. Species overview

- **ToL ID**: `drLytSali1`
- **Species**: *Lythrum salicaria* L. (purple loosestrife), family Lythraceae
- **Assembly**: curated stage, 2 source FASTA files (`drLytSali1.hap1.1`,
  `drLytSali1.hap2.1`), each a grab-bag mixing multiple haplotype tags —
  manifest uses `AUTO` haplotype detection with a unified `HAP(\d+)` regex,
  correctly recovering **4 haplotype copies** per chromosome from the
  header tags rather than the 2 the raw file count would naively suggest
- **Structure**: 15 chromosome numbers, 4 copies each (60 chromosome-scale
  units total)

## 2. Literature context

**Ploidy level**: tetraploid, native Eurasian range has both diploid and
autotetraploid cytotypes; invasive North American populations are
uniformly autotetraploid (Kubátová et al. 2008, *Journal of Biogeography*;
Balogh 2018).

**Origin — autopolyploid, high confidence.** This is one of the
best-evidenced calls in the whole panel, resting on classical genetics
evidence rather than modern genomics, replicated independently since the
1940s:

- **Tetrasomic inheritance** directly demonstrated via segregation ratios
  (Little 1958, *The Botanical Review*: "chromosome pairing and disjunction
  are completely at random in tetraploid *Lythrum salicaria*"; Lewis &
  Jones 1992, citing Fisher & Mather 1943's original linkage analysis of
  the style-length locus).
- **Double reduction** — a segregation phenomenon specific to polysomic
  (auto-)inheritance, essentially absent under disomic (allo-)
  inheritance — documented for this species since the 1940s and still
  cited as a reference case (Lamon et al. 2026, bioRxiv, reviewing double
  reduction across polyploids).
- **Breeding evidence**: full fertility in experimental hybrids, and
  successful backcrossing to *L. salicaria* from cultivar hybrids with the
  related *L. virgatum*, consistent with autopolyploidy rather than the
  reduced/asymmetric fertility expected from a fixed allopolyploid hybrid
  (Houghton-Thompson 2001, PhD thesis).

Key citations: Houghton-Thompson 2001; Lewis & Jones 1992; Balogh 2018,
2024; Chun, Nason & Moloney 2009 (*Molecular Ecology*); Little 1958;
Kubátová et al. 2008; Ford 1945; Lamon et al. 2026. Full detail in
`meta/literature_ploidy.tsv`.

## 3. Direct data outputs

| metric | value | note |
|---|---|---|
| Structural pairing | 10/15 chromosomes (5 pairs) | `chr02↔chr03`, `chr13↔chr14`, `chr11↔chr15`, `chr04↔chr07`, `chr05↔chr10` |
| Whole-genome bulk distance | ~0.02–0.03 across all 15 chromosomes | uniform, close to a single-genome haplotype-pair baseline |
| `te_marker_fraction` | mean 0.307, range 0.117–0.670 (n=90 pairs) | *elevated* — above both confirmed allo anchors (`daGleHede1` 0.223, `drTriRepe1` 0.169) |
| `partition_consistency` | 0.40 (rotating) | modal singleton `HAP4`, but rotates between `HAP2`/`HAP3`/`HAP4` depending on chromosome — not artifact-explained |
| `distance_ratio_cv` | 0.083 (tight) | — |
| `mean_run_length_windows` / `flip_rate` | 3.72 / 0.263 | shortest run-length among the species this is computable for |

`ddHesMatr1`-style singleton-artifact check: not applicable (rotates
rather than fixing on one haplotype, so this isn't the `ddHypMacu1`-style
fragmented-assembly confound).

## 4. Process inference

This is the clearest case in the panel of **Pattern 2: tetrasomic
homogenization** (see `INTERPRETATION.md`'s "reading agreement and
conflict" section) — magnitude-based metrics read allo-like, the
identity-based metric reads auto-like, and the mechanism explains why both
are correct readings of different things:

- **Why the magnitude metrics (`distance_ratio_cv`, bulk distance,
  run-length) all read low/tight/allo-like**: tetrasomic inheritance means
  all 4 copies pair and recombine in random multivalent configurations
  every generation, not strict copy-vs-copy. That's a continuous sequence
  *homogenizer* — any point mutation arising on one copy has a real chance
  of being shuffled onto or off the others via crossover and gene
  conversion, which is exactly why bulk distance sits at ~0.02–0.03,
  barely above a single-genome's own haplotype-pair baseline.
- **Why `te_marker_fraction` reads high/allo-like anyway**: TE insertions
  are discrete presence/absence events added by transposition, not
  substitution — removing or spreading one requires an exact excision or a
  crossover breakpoint landing precisely at that locus, a far less
  efficient homogenizing mechanism than the gene conversion that erases
  point mutations. Independent, stochastic post-polyploidization TE
  activity per copy is well documented in autopolyploids, not just
  allopolyploids, and needs no separate parental origin to produce a real,
  elevated signal.
- **Why `partition_consistency` is the one metric that sees through this**:
  it tests lineage *identity* stability across chromosomes, not
  divergence *magnitude* — and it's the one axis tetrasomic homogenization
  doesn't touch. The rotating singleton (not a fixed haplotype every time)
  is itself additional evidence *against* a fixed two-parent structure:
  a real allopolyploid should show the *same* two-way partition on every
  chromosome, not one whose "odd one out" changes identity.

## 5. Auto/allo statement

**Autopolyploid, high confidence.** The classical genetics literature
(tetrasomic inheritance, double reduction, breeding/backcross fertility,
replicated independently since the 1940s) is definitive on its own and
doesn't need this pipeline's confirmation — but `partition_consistency`'s
rotating-singleton result independently corroborates it via a completely
different, sequence-based mechanism, and the apparent contradiction from
`te_marker_fraction`/bulk-distance is fully explained (§4), not
unresolved. Nothing here would change this call; the te_marker_fraction
tension that originally motivated deeper investigation of this species is
now a resolved, mechanistically-understood feature of autopolyploidy under
tetrasomic inheritance, not evidence against it.
