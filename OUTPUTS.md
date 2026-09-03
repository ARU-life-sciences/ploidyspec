# ploidyspec output reference

What every file means, how to read the numbers, and worked examples using
real output already on disk (`daGleHede1` — *Glechoma hederacea*, a known
allopolyploid — and the `lpTriTurg1_AB` durum wheat validation run, whose
A/B subgenomes are independently known ground truth). **Both of these are
confirmed allopolyploids**, not an auto/allo pair — `daGleHede1`'s
haplotype-copy distance (0.0099) is used below purely as a within-individual
heterozygosity reference point, not as an autopolyploid-like baseline. The
panel currently has no confirmed autopolyploid or diploid calibration point
for the low end of the te-marker auto/allo scale (see `subgenomes/`
section) — that's a real gap, not an oversight.

## Directory layout

```
results/<species>/
  sequences.tsv, unplaced.tsv     prepare's output -- every stage below reads sequences.tsv
  matrix/                         whole-chromosome distance, ploidy, homologous-chromosome reports
  windowed/                       sliding-window divergence between haplotype copies of the same chromosome
  homeologs/                      ancient (paleopolyploid) homeolog pairs and their windowed tracks
  subgenomes/                     fossil-TE marker output, windowed subgenome painting, auto/allo index
```

Everything under `matrix/`/`windowed/`/`homeologs/`/`subgenomes/` is derived,
in that dependency order — `windowed`, `homeologs` and `subgenomes` all need
`matrix` to have run first (and `subgenomes` needs `windowed` too, for the
`windowed_distance_cv` column). `sequences.tsv`/`unplaced.tsv` sit at the top
level since everything else depends on them.

## `sequences.tsv` / `unplaced.tsv` (from `prepare`)

`sequences.tsv`: one row per chromosome-scale unit that made it into the
analysis — `unit_id` (e.g. `HAP1_chr01`), `hap`, `chrom`, `seq_id`, `length`,
`source`, `desc`. This is the master list every other stage groups by
`chrom` (to find haplotype copies of the same chromosome) or reads directly.

`unplaced.tsv`: every sequence that got *excluded*, with why (`reason`:
`below-min-len`, `no-hap-match`, `no-chrom-match`). Read this first whenever
a species produces fewer units than you expected — it's almost always a
header-format/regex mismatch, not a bug.

## `matrix/` — whole-chromosome divergence, ploidy, homologous chromosomes

**`whole_chrom_distance_matrix.csv`**: the square all-vs-all distance
matrix, values are Mash-corrected divergence (see below), not raw k-mer
Jaccard. Feeds the heatmaps and `homeologs/`.

**`whole_chrom_pairs.tsv`**: long-form version, one row per pair, with the
full detail: `k_sweep`, `chosen_k` (the largest swept k that cleared the
noise floor for that pair — see `resolution_limited`), `jaccard_at_chosen_k`,
`containment_at_chosen_k`, `distance` (the Mash-corrected value, what's in
the matrix), `containment_p` (a secondary containment-based divergence
estimate, cross-check), `resolution_limited`, `k_consistency_spread`.

**How to read `distance`**: it's an estimated per-site substitution rate
(Mash correction, Ondov et al. 2016), *not* a bounded fraction — a bulk
`1-Jaccard` distance conflates k-mer-level and base-level divergence (one
substitution knocks out up to *k* neighboring k-mers), so raw Jaccard
overstates true divergence non-linearly. Real anchors from `daGleHede1`:

| Comparison | `distance` |
|---|---|
| True homolog pair (`HAP1_chr01` vs `HAP2_chr01`, same individual) | **0.0099** (~1%) |
| Ancient homeolog pair (`chr02` vs `chr03`, retained WGD duplicates) | **0.047** (~5%) |
| Unrelated chromosome pair (`HAP1_chr01` vs `HAP1_chr02`) | **0.11–0.13** (~11–13%) |

**`resolution_limited`**: `True` means no k in the sweep cleared the
chance-collision noise floor for that pair (see `mash.py`) — the distance
value exists but is closer to noise than signal, treat it as a lower bound /
"beyond what this k-mer sweep can resolve," not a trustworthy number.

**`whole_chrom_multi_k.tsv`**: the full per-(pair, k) sweep behind the
consensus `distance` — every k's raw jaccard/containment/distance and
whether *that specific k* was resolution-limited. This is what
`k_resolution_diagnostic.png` plots.

**`whole_chrom_distance_heatmap.png`**: clustered heatmap, color scale
calibrated to the matrix's own max value (not a fixed 0–1 — most real
comparisons land under 0.2, so a fixed scale looks washed out).
**`whole_chrom_distance_heatmap_contrast.png`**: the diagonal-masked,
2nd–98th-percentile-stretched variant (the Hi-C-style trick) — use this one
specifically to look for subgenome block structure in the *background*,
since the standard heatmap's dynamic range is dominated by the near-zero
true-homolog cells.

**`ploidy_summary.tsv`**: one row per chromosome number — `n_haplotype_copies`
(the ploidy level *for that chromosome*, e.g. 4 for a tetraploid) and
`haplotype_labels`. Worth checking even when you already "know" the
species' ploidy: `daGalBore1`'s manifest nominally listed 2 files, but
AUTO-detection recovered 4 real haplotype copies per chromosome from header
tags, matching *Galium boreale*'s known tetraploid cytology (2n=4x=44) —
this file is where that shows up per-chromosome.

**`homologous_chromosomes.tsv`**: the same-chromosome-number subset of
`whole_chrom_pairs.tsv` — i.e. only the *true* haplotype/homologous pairs
(as opposed to `homeologs/homeolog_pairs.tsv`'s cross-chromosome-number
*ancestral* pairs), pulled into its own view so you don't have to filter the
full all-vs-all table to answer "what are chr01's copies and how different
are they."

**`k_resolution_diagnostic.png`**: per-pair distance-vs-k, red where a pair
falls below the noise floor at that k — the empirical answer to "how far
back can this k-mer sweep actually see."

## `windowed/` — divergence along the chromosome

**`windowed_chrNN.tsv`** (one per chromosome number) and **`windowed_all.tsv`**:
`jaccard_distance` per window per haplotype-copy pair. **This track is raw
Jaccard, not Mash-corrected** — `windowed` stayed uncorrected deliberately
(re-deriving Mash correction per window would multiply the already-dominant
cost of this stage). Use it for *relative* shape along the chromosome
(where does divergence spike/dip), not as a calibrated divergence estimate —
for that, use `matrix/`'s `distance`.

**`windowed_genome_overview.png`**: small-multiples of every chromosome's
track on one page, x-axis normalized to % of chromosome length so different
lengths are comparable.

**`windowed_chrNN_heatmap.png`** (one per chromosome) and
**`windowed_genome_overview_heatmap.png`**: the same per-window
`jaccard_distance` data as a heatmap instead of a line plot — one row per
haplotype-copy pair, one column per window, color = distance (viridis,
2nd–98th percentile contrast-stretched, matching `matrix/`'s heatmaps). More
readable than the line-plot version once a chromosome has more than two or
three haplotype copies (line plots overlap into noise; a mosaic of
lighter/darker patches along one row is still legible) — the overview
version uses one shared color scale across every chromosome so patches are
comparable genome-wide, not just within one chromosome. A patchy row (short
stretches that break from an otherwise-uniform distance) is a candidate
partial-rediploidization or homeologous-exchange region for that pair.

## `homeologs/` — ancient (paleopolyploid) chromosome pairs

**`homeolog_pairs.tsv`**: candidate retained-duplicate chromosome pairs from
an ancestral whole-genome duplication, found by testing every
cross-chromosome-number pair's distance against the empirical background of
*all* cross-chromosome-number distances (`z_score`, one-sided lower-tail
test against a sigma-clipped normal fit to that background), then
Benjamini-Hochberg FDR correction (`q_value`, default threshold 0.05, `--fdr-alpha` to
change it) across every candidate tested simultaneously, then a greedy 1:1
match (each chromosome gets at most one ancestral partner).

**Read `z_score` alongside `p_value`/`q_value`, not instead of them**: once
`|z|` gets much past ~9, `p`/`q` underflow to exactly `0.0` in float64 and
stop conveying *how much* evidence there is — `z` stays finite. Real
example, `daGleHede1` (all 9 pairs land far past that point):

```
chrom_a  chrom_b  mean_distance  z_score   p_value  q_value
chr02    chr03    0.047344       -20.686   0        0
chr14    chr17    0.054208       -18.610   0        0
```

Both are "p≈0, definitely significant," but `-20.7` vs `-18.6` still tells
you `chr02`↔`chr03` sits further from the background than `chr14`↔`chr17` —
information the p/q columns alone can no longer show at this depth.

`n_haplotype_copy_pairs`, `resolution_limited`, `n_resolution_limited_of_total`:
same meaning as in `matrix/`, aggregated across every haplotype-copy
combination between the two chromosome numbers. If `resolution_limited` is
`True` here, the pairing's acceptance is resting on a distance estimate that
was itself below the noise floor — treat it as lower-confidence even if `q`
looks small.

**`homeolog_candidates_ranked.tsv`**: every cross-chromosome-number pair
tested, not just the FDR-accepted subset in `homeolog_pairs.tsv` — ranked by
ascending distance, each row carrying its own `z_score`/`p_value`/`q_value`
and an `accepted` flag. Exists because FDR correction gets more conservative
as chromosome count grows (more simultaneous tests), so a real, moderate
effect size can fail significance on a many-chromosome species for the same
reason a weak one would on a small species — check this file, not just
`homeolog_pairs.tsv`, before concluding "no ancient signal" for a species
with many chromosomes.

**`ploidy_ancestry_summary.tsv`**: one row per chromosome number, combining
the contemporary haplotype-copy count (from `matrix/ploidy_summary.tsv`)
with its accepted ancient partner (if any) — the two signals that otherwise
require cross-referencing `matrix/ploidy_summary.tsv` and
`homeolog_pairs.tsv` by hand. `distance_ratio` (homeolog-pair distance
divided by the mean of both chromosomes' own within-chromosome distance) is
a diagnostic, **not a classifier**: a ratio near 1 means the "ancient"
partner is actually about as close as a chromosome's own contemporary copy
— worth checking by hand whether that reflects a real biological signal
(e.g. a genuinely young/tetrasomic-like duplication, or genus-specific
chromosome fission/fusion biology) rather than assuming it means the same
thing a ratio in the hundreds does (deep, clearly-ancient divergence, as
seen in some species in this project's own panel).

**Fused-chromosome homeolog pairs — a real category, not automatically an
artifact**: a homeolog pair can come out significant because one
haplotype's copy of a chromosome isn't a clean single unit — it's a scaffold
containing another chromosome's sequence too, usually visible as roughly
double the length of that chromosome's other haplotype copies, paired with
that other chromosome being entirely *absent* from the same haplotype(s).
**This can be a genuine, biologically real chromosome fusion, not a
mislabeling bug — check before concluding either way.** Case in point,
`SchCurv1` (`chr19`↔`chr22`): `HAP3_chr19`/`HAP4_chr19` are ~69.5Mb versus
`HAP1_chr19`/`HAP2_chr19`'s ~38–40Mb (close to `chr19`+`chr22` combined),
`chr22` is missing entirely from `HAP3`/`HAP4`, and the windowed heatmap
shows the match to `chr22` confined to specific blocks within `HAP3`/`HAP4`
only. An earlier version of this note called this an assembly artifact to
be discarded — that was wrong. `Xie et al. 2026, Nature`
(doi:10.1038/s41586-026-10439-1) independently confirms, via synteny to an
outgroup, Hi-C, and meiotic FISH, that `chr19`+`chr22` is a real ancestral
**unbalanced chromosome fusion** shared across all snow carp genera:
only 2 of the 4 ancestral chromosome copies fuse, the fused pair then
pairs preferentially with itself at meiosis (bivalent) instead of the
original 4-way (tetravalent) pairing, and *that* is what triggers disomic
inheritance and ohnologue divergence at that locus — this is literally the
mechanism the paper's title describes. Under that model, a match **confined
to a block** rather than spread genome-wide is the *correct* signature for a
physically fused ohnologue pair, not evidence against it — only the
fused-in segment should match, and it should only match in the haplotype
copies that carry the fusion. Don't read "localized, not genome-wide" as an
automatic artifact flag; read the length/absence pattern first, and treat
"confined to a block" as consistent with a real fusion rather than
disqualifying one.

Before trusting *any* `homeolog_pairs.tsv` hit, still check whether either
chromosome's `n_haplotype_copies` is lower than its partner's (visible in
this same file) and whether the short haplotype(s) are anomalously long
relative to their siblings — that tells you whether you're looking at a
straightforward retained duplicate or a fusion-derived one, which changes
how you read `distance_ratio` (it now describes resolution *at the fused
locus specifically*, not the whole chromosome) but doesn't by itself mean
discard the pair. `ploidy_ancestry_summary.tsv` may carry a manually-added
`notes` column (present for `SchCurv1`) recording this kind of
investigated, per-pair verdict — it's not part of the automated schema, so
its absence elsewhere in the panel means "not yet checked," not "confirmed
clean."

**`windowed_chrAAxBB.tsv/.png`, `windowed_chrAAxBB_heatmap.png`,
`windowed_homeologs_all.tsv`, `windowed_homeologs_overview.png`,
`windowed_homeologs_overview_heatmap.png`**: same shape as `windowed/`'s
files (including the heatmap variants), but tracking divergence along the
*ancestral* pairing instead of true haplotype copies — same raw-Jaccard
caveat applies.

## `subgenomes/` — differential fossil-TE markers (the resolver/phaser)

Only present for species where `te-markers`/`te-markers-windowed` have
actually been run (opt-in via `all --with-te-markers`/`--with-te-markers-windowed`,
or as standalone stages — not run by plain `all`, since `te-markers-windowed`
in particular is expensive, ~8h for wheat's large chromosomes). Implements
the Jaron/Cerca method
(github.com/KamilSJaron/k-mer-approaches-for-biodiversity-genomics, linked
in `README.md`): isolates k-mers that are both **high-copy**
(repetitive/TE-like, count ≥ `--min-count`, default 100) and
**differentially represented** (≥ `--min-ratio`, default 2x) between two
haplotype copies of the same chromosome — independent repeat-family
expansion in each parental lineage is a much more targeted allopolyploidy
signal than bulk divergence.

**`te_markers_<a>x<b>.tsv`**: every differential marker k-mer for that pair
— `kmer, count_a, count_b, assigned` (`a` or `b`, whichever side it's
enriched in).

**`te_markers_summary.tsv`**: per-pair counts — `n_markers_a`/`n_markers_b`
(differential markers favoring each side) and `n_highcopy_a`/`n_highcopy_b`
(*total* high-copy k-mers per side, markers included). This is where the
scale of the signal first becomes visible:

| | `n_markers_a` + `n_markers_b` (per chromosome) |
|---|---|
| `daGleHede1` haplotype copies (confirmed allopolyploid) | ~600–5,700 |
| Wheat's known A vs B subgenomes | ~330,000–605,000 |

That 100–1000x gap in raw counts is driven mostly by genome/repeat-content
size (wheat's genome and total repeat load dwarf `daGleHede1`'s), not
directly comparable across species — see `te_marker_fraction` below for the
size-normalized version, which tells a different story.

**`te_markers_windowed_<a>x<b>.tsv/.png`** (the phaser's output): each
haplotype copy's own windows, tiled and checked against *both* marker sets.
`assigned` is `a`/`b`/`ambiguous`/`none`. A haplotype's windows should mostly
match its *own* side (`HAP1`'s windows matching the `a`-side markers that
were themselves derived by comparing `HAP1` against `HAP2`) — a window that
instead matches the *other* side is a candidate homeologous-exchange or
introgression site, not noise. On wheat's known subgenomes this recovers
88–95% self-matching with near-zero cross-matching, visually a near-solid
two-color split across the whole chromosome — the strongest validation this
method has. On `daGleHede1` most windows are correctly `ambiguous` (most of
a genome isn't repeat-dense enough to carry a strong signal either way),
with real signal concentrated in specific repeat-rich regions.

**`auto_allo_index.tsv`**: `te_marker_fraction` =
`(n_markers_a + n_markers_b) / (n_highcopy_a + n_highcopy_b)`, bounded [0,1]
by construction (markers are a subset of high-copy k-mers) — "what fraction
of each haplotype's high-copy repeat content is subgenome-differential."
Near 0: nearly identical repeat content (autopolyploid-like). Near 1:
almost entirely non-overlapping (allopolyploid-like). Report and compare the
continuous value; don't treat any specific number as a verdict.

`daGleHede1`'s range (0.11–0.42 across its 18 chromosomes, mean ~0.22) and
wheat's (0.40–0.44) overlap on their high ends — consistent evidence the
metric picks up real allopolyploid signal in both, at different intensities
(plausibly different subgenome divergence times/degrees). The panel's first
confirmed autopolyploids, `SchCurv1`/`SchYoun1` (snow carp, Xie et al. 2026),
sit at genome-wide means of 0.209/0.223 — indistinguishable from
`daGleHede1` at that resolution. **This isn't a failure of the metric; it's
the wrong resolution to read it at.** See
`te_marker_fraction_by_lineage.tsv` below.

**`te_marker_fraction_by_lineage.tsv`**: for every chromosome with ≥3
haplotype copies, an automatic 2-way split of the copies by whole-chromosome
distance (`matrix/whole_chrom_distance_matrix.csv`, reused — no new k-mer
work), then `te_marker_fraction` averaged separately **within** each group
versus **across** them (`split_ratio` = cross ÷ within). Exists because a
flat, unweighted genome-wide mean can hide real structure in a genome that's
only *partly* resolved into two lineages — confirmed case: `SchCurv1`'s
`chr19` (the one confirmed ancestral chromosome fusion in that species) has
within-lineage fraction 0.11 but cross-lineage 0.59 (`split_ratio` 5.6),
buried inside a flat per-chromosome mean around 0.35 and invisible in the
genome-wide mean entirely. Every other chromosome in `SchCurv1` sits at
`split_ratio` 0.8–2.2, i.e. no comparable structure — except `chr17`
(`split_ratio` 3.2, and 1.9 independently in `SchYoun1`), which is *not*
one of the five known chromosome fusions. The paper documents a **second,
separate rediploidization mechanism** at exactly this chromosome — a
centromeric inversion producing a partial, short-arm-only disomic pattern
— which plausibly explains a real but weaker split_ratio than a full
fusion like `chr19`'s. Not independently confirmed here (would need the
same length/windowed-heatmap check as any homeolog hit — see the fused-
scaffold caveat above), but a strong, unprompted candidate worth checking
before assuming it's noise. The grouping (`group_a`/`group_b` columns) is a
simple farthest-pair-seeded 2-way split, not a rigorous clustering method —
treat `split_ratio` as a diagnostic, same as `distance_ratio` in
`homeologs/`, not a classifier: most chromosomes in most species will show
`split_ratio` near 1 (no real structure), and that's the expected, correct
result, not a bug.

**Assembly-quality confound, checked and found in one species so far**:
a haplotype assembly that's more fragmented than its siblings can inflate
`te_marker_fraction` for every pair it's in, independent of any real
biology — fragmented/incomplete regions can look like spurious
"differential" content. Confirmed for `ddHypMacu1`: `unplaced.tsv` shows its
`HAP1` file has 570 excluded scaffold fragments vs. 2-3 for `HAP3`/`HAP4`
(essentially perfectly chromosome-scale), and `HAP1`-involving pairs average
`te_marker_fraction` 0.598 vs. 0.291 for `HAP2`/`HAP3`/`HAP4`-only pairs —
almost double. **Use 0.291, not the raw panel mean (0.445, which is
`HAP1`-inflated), as `ddHypMacu1`'s trustworthy value.** Checked
panel-wide (comparing `unplaced.tsv` counts across each species' source
files) for the same pattern: `drAriEdul1` also has a lopsided assembly
(`HAP1`/`HAP2` far more fragmented than `HAP3`/`HAP4`) but its
`HAP1`/`HAP2`-involving pairs average *lower* (0.209) than
`HAP3`/`HAP4`-only pairs (0.270) — no inflation there, its reported value
stands as-is. No other species in the panel shows a comparably lopsided
`unplaced.tsv` distribution among source files that had `te-markers` run.
Before trusting any single-species `te_marker_fraction` headline number,
check `unplaced.tsv` for this asymmetry and, if present, verify (as done
here) whether it actually correlates with elevated fraction values before
assuming inflation.

`windowed_distance_cv` (coefficient of variation of the
raw-Jaccard `windowed/` track for that pair) is a secondary, corroborating
column — higher means more patchy/heterogeneous divergence along the
chromosome, consistent with (but not proof of) mosaic subgenome structure.

**`subgenome_windows_summary.tsv`**: per-unit window-assignment counts —
`n_self`/`n_other`/`n_ambiguous`/`n_none` and their percentages. The
directly readable version of "how much of this haplotype copy phases
cleanly as itself."

**`subgenome_anomalous_windows.tsv`**: every window where a haplotype's own
sequence matched the *other* side's markers — candidate
homeologous-exchange/introgression coordinates, ready to cross-reference
against `windowed/`'s divergence track or an assembly viewer.
