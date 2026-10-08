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
  windowed/                       divergence along each copy against the other copies of the same chromosome
  homeologs/                      ancient (paleopolyploid) homeolog pairs and their windowed tracks
  subgenomes/                     fossil-TE marker output, windowed subgenome painting, auto/allo index
  structure/                      inheritance-mode metrics, diffuse genome-wide chromosome partitions
  rediploidization/               chromosome fusions, per-chromosome lineage structure and rediploidization state
```

Everything under `matrix/`/`windowed/`/`homeologs/`/`subgenomes/` is derived,
in that dependency order — `windowed`, `homeologs` and `subgenomes` all need
`matrix` to have run first (and `subgenomes` needs `windowed` too, for the
`windowed_distance_cv` column). `sequences.tsv`/`unplaced.tsv` sit at the top
level since everything else depends on them.

## `sequences.tsv` / `unplaced.tsv` (from `prepare`)

`sequences.tsv`: one row per chromosome-scale unit that made it into the
analysis — `unit_id` (e.g. `HAP1_chr01`), `hap`, `chrom`, `seq_id`, `length`,
`source`, `desc`, `numbering`. This is the master list every other stage
groups by `chrom` (to find haplotype copies of the same chromosome) or reads
directly. `numbering` is `names` when the chromosome number came from the
header, `auto` when ploidyspec assigned it (see README, "Naming your
chromosome-scale sequences"). With automatic numbering, one haplotype is
numbered by length and the others are matched to it in the matrix stage; a
sequence that matched nothing keeps a provisional number from 1001
(`HAP2_chr1001`).

`unplaced.tsv`: every sequence that got *excluded*, with why (`reason`:
`below-min-len`, `no-hap-match`, `no-chrom-match`, or with automatic numbering
`auto-not-chromosome-scale`: shorter than a tenth of the median length of the
haplotype's larger sequences). Read this first whenever
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
one row per window of one copy (`unit_a`, positions along *its own*
coordinates) against another copy (`unit_b`). Each window's k-mers are looked
up anywhere in `unit_b`'s whole chromosome (FastK profile against `unit_b`'s
k-mer table, `--window-k`), so no alignment or shared coordinates are needed:
- `kmers_a`: valid k-mer positions in the window (assembly gaps excluded).
- `shared`: how many of them occur anywhere in `unit_b`.
- `containment` = `shared`/`kmers_a`; `distance` = -ln(`containment`)/k
  (Mash-style; 1.0 when nothing is shared).

Rows are directional: A→B and B→A are separate tracks along each copy.
Window distances run lower than `matrix/`'s whole-chromosome `distance`
(about 0.6x on `ddEmpNigr1`, because a window is matched against the whole
other copy), with the same ordering of pairs. Use them for shape along the
chromosome, and `matrix/` for calibrated divergence.

Before 2026-10-07 this stage compared windows at *equal coordinates* in every
copy (column `jaccard_distance`). That assumed collinear assemblies and fell
out of register after the first indel or gap larger than a window: a pair
0.004 apart over the whole chromosome read ~0.99 Jaccard distance in most
windows. Outputs with a `jaccard_distance` column are from that version and
should be re-run.

Memory: FastK's profile holds the other copy's table, ~14 bytes per bp of
the longest chromosome per thread (the stage logs its estimated peak).

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
*ancestral* pairing instead of true haplotype copies (each copy's windows
against the other chromosome number's copies).

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

**`flagged_units` column**: automatically marks units that look like an
incomplete/fragmented assembly rather than a real second lineage — a unit
whose sequence length *or* total high-copy k-mer count sits below 75% of a
*majority-agreed* baseline among its same-chromosome siblings (median of
whichever units mutually agree with each other; requires that agreeing set
to be a strict majority, not just the largest cluster — see `flag_low_content_units`
in `subgenome_report.py`). The majority requirement matters: a real 2-vs-2
chromosome fusion also produces a length/high-copy outlier by a naive
median comparison (the group median gets pulled toward the longer, fused
side), but splits the group evenly — neither side is a majority — so it's
never mistaken for this. A real fusion or subgenome split doesn't reduce a
unit's own length or repeat content, only how much of it is shared with
specific other copies, so this is a genuinely different check from the
bipartition itself, computed independently and reported alongside it.
Confirmed case: `daBudDavi1`'s `chr05` split at 6.13× —
comparable in size to `SchCurv1`'s real `chr19` fusion signal — but as a
lopsided **1-vs-3** split (`HAP2` alone), not a clean 2-vs-2. `HAP2_chr05`
is 32% short (25.5Mb vs. siblings' 37–38Mb) and has a third of their
high-copy k-mer count, tracking its siblings normally (within 3–5%) on
every *other* chromosome — chromosome-specific, not a systemic assembly
problem with that haplotype. Every pair involving it is lopsided (~10,000+
markers favoring the other side, `HAP2` contributing 8–18 of its own),
exactly the `ddHypMacu1` pattern below. `flagged_units` catches this
automatically now: `chr05`'s row carries `HAP2_chr05` in that column. A
clean 2-vs-2 split with no flagged units (`SchCurv1`'s `chr19`, `chr17` in
both snow carps) is the shape a real signal takes; a lopsided split with a
flagged unit is the shape this specific artifact takes — but `flagged_units`
being empty doesn't guarantee the split is real, only that this one known
failure mode has been ruled out.

**`high_density_units` column**: a separate, more cautious check for a
different failure mode — a unit with *normal* length but anomalously high
repeat **density** (high-copy k-mers per bp, not raw count) relative to a
majority-agreed baseline among its siblings (>1.75× the baseline; see
`flag_high_density_units`). Density, not raw high-copy count, matters here:
a real chromosome fusion also elevates a unit's raw high-copy count (a
longer, fused sequence simply contains more total repeat content), which a
raw-count version of this check would wrongly flag — confirmed against
`SchCurv1`'s real `chr19` fusion, whose fused lineage has ~2.3–2.6× the raw
high-copy count of the unfused lineage but only ~1.35–1.4× the *density*,
roughly proportional to the length increase rather than a density anomaly,
and correctly unflagged. `chr17` (candidate, associated with a documented
centromeric-inversion rediploidization mechanism per Xie et al. 2026) sits
at ~1.4–1.6× density, also unflagged — deliberately: unlike missing
content, which is never real, elevated density *can* be genuine biology
(independent TE activity in a real second lineage is exactly what
`te_marker_fraction` exists to detect), so this column is a caution to go
check further, not a verdict the way `flagged_units` is closer to being.
Confirmed case: `ddHypMacu1`'s `chr06` (split_ratio 4.73, initially read as
a clean candidate since `flagged_units` — a length/raw-count check — didn't
catch it) has normal length but ~2× the density of its three
mutually-agreeing siblings, clearing the threshold. Corroborating evidence:
the same haplotype (`HAP1`) shows 1.9–5.1× elevated density on 3 other
chromosomes and 0.43× (deflated) on a fourth — scattered in both
directions, consistent with `HAP1` being independently confirmed as the
lower-quality assembly for this species (570 excluded scaffold fragments
below) rather than four coincidental real biological signals.

**Assembly-quality confound, now checked automatically (`flagged_units`
above), confirmed in two species**: a haplotype file that's more
fragmented/incomplete than its siblings can inflate `te_marker_fraction`
for every pair it's in, independent of any real biology — missing content
trivially looks like spurious "differential" content by construction.

**Why this specifically hits `HAP1`/`HAP2` and not `HAP3`/`HAP4`+ — a
curation-pipeline artifact, not necessarily a sequencing/assembly-quality
difference.** For any species with more than 2 haplotype copies (an
auto-tetraploid-style curated assembly), ToL's curation process only runs
the unlocalized-sequence placement step on the primary `HAP1`/`HAP2` pair.
Any additional haplotype file (`HAP3`, `HAP4`, …) is built only from
already-placed `SUPER` (chromosome-scale) scaffolds plus whatever
unlocalized sequences trivially came with them — it never goes through that
placement step. So `HAP3`/`HAP4`+ will *structurally* show near-zero
`unplaced.tsv` entries for any such species, by construction, not because
their underlying assembly is more complete or higher quality. Expect the
`HAP1`/`HAP2` vs. `HAP3`/`HAP4`+ asymmetry in `unplaced.tsv` fragment counts
as the *default* for every >2-haplotype species in the panel — it is not on
its own evidence of a problem, and shouldn't be read as one without checking
whether it actually correlates with an inflated `te_marker_fraction` (below).

Confirmed for `ddHypMacu1`: `unplaced.tsv` shows its `HAP1` file has 570
excluded scaffold fragments vs. 2-3 for `HAP3`/`HAP4`, and `HAP1`-involving
pairs average `te_marker_fraction` 0.598 vs. 0.291 for `HAP2`/`HAP3`/`HAP4`-only
pairs — almost double, so in this case the asymmetry *does* correlate with
real inflation (plausibly because content that never got placed into
`HAP1`'s chromosome-scale scaffolds is genuinely absent from what
`te-markers` compares, not just an unplaced-file bookkeeping difference).
**Use 0.291, not the raw panel mean (0.445, which is `HAP1`-inflated), as
`ddHypMacu1`'s trustworthy value.** Confirmed a second time for
`daBudDavi1`'s `chr05`/`HAP2` (above, a 2-haplotype species — same failure
mode, different cause, since the HAP1/HAP2-curation explanation above
doesn't apply there). Checked panel-wide (comparing `unplaced.tsv` counts
across each species' source files) for the same correlation: `drAriEdul1`
also has a lopsided assembly (`HAP1`/`HAP2` far more fragmented than
`HAP3`/`HAP4`) but its `HAP1`/`HAP2`-involving pairs average *lower* (0.209)
than `HAP3`/`HAP4`-only pairs (0.270) — the structural asymmetry is present
(as expected for any >2-haplotype species) but doesn't correlate with
inflation here, so its reported value stands as-is. **The lesson: don't
skip the check just because the asymmetry itself is expected — always
verify whether it correlates with an actual `te_marker_fraction` difference
before deciding whether to correct for it.** Before trusting any
single-species `te_marker_fraction` headline number, check `unplaced.tsv`
for this asymmetry and, if present, verify (as done here) whether it
actually correlates with elevated fraction values before assuming
inflation.

`windowed_distance_cv` (coefficient of variation of the
`windowed/` track of `unit_a`'s windows against `unit_b`) is a secondary, corroborating
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

## `rediploidization/` — fusions, lineage structure, rediploidization state

Pulls together the rediploidization evidence the other stages compute, plus one
new test (fusion detection), into one per-chromosome table. Run after `matrix`,
`homeologs` and `windowed`; uses `te-markers` output when present. Only the
unplaced-scaffold fusion test does new k-mer work, and only for the handful of
unplaced sequences long enough to be whole chromosomes.

**`fusions.tsv`**: one row per scaffold with a fusion signal. Two kinds of
scaffold are tested:
- `placed_long_copy`: a placed copy at >= `--long-ratio` (1.4) x the median
  length of its same-chromosome siblings (e.g. `SchCurv1`'s 69.5 Mb `HAP3_chr19`
  vs 38-40 Mb siblings). Uses the matrix stage's shared-k-mer counts.
- `unplaced_scaffold`: a sequence in `unplaced.tsv` that matched no chromosome
  number but is >= `--orphan-min-frac` (0.5) x the median chromosome length
  (e.g. `SchYoun1`'s `Sy_Chr04_15_M2`). New k-mer table, intersected with
  every placed copy. Results are cached in `orphan_containment.tsv`.

`containment` = fraction of a chromosome copy's k-mers found in the scaffold
(best copy per chromosome number). A chromosome counts as a component when its
containment is >= `--containment-z` (10) robust z-scores above the species'
background and >= 2x the background median. `SchCurv1`: chr22 in the fused chr19
copies = 0.44-0.45, background median 0.063, z ~ 48.

`sibling_excess` (placed long copies only): for each component, the share of the
long copy's own k-mers found in that chromosome, divided by the same share in its
least-enriched sibling copy. A component is kept only at >= `--sibling-excess`
(2.5). Normalising by each copy's own k-mer count means a fragmentary sibling
reads the same per-k-mer rate as a complete one. A fused copy carries the
partner and its unfused sibling doesn't: `SchCurv1` chr22 in `HAP3_chr19` =
3.4x. Ancient homeology is shared by every copy: `drMyrSpic1` `HAP1_chr13`,
"long" only because its sibling is a 15 Mb fragment, has its homeologs chr03 and
chr09 at 0.9x. `daPilAura1`'s long `HAP1_chr01` has chr06/chr07/chr14 at
0.6-1.9x. All of these were `candidate_partner_present` before this check.
Hits that fail it still appear in `sibling_excess` but not in `components`.

`lineage_divergence` = -ln(containment)/k: a Mash-style divergence between the
fused lineage and the unfused copies of each component. Use it to rank fusions
by age. `SchYoun1` gives chr19+22 ~0.053, chr04+15 / chr08+16 / chr20+23
~0.026-0.028, and chr11+14 ~0.015. That is the three-wave grouping and order of
Xie et al. 2026, oldest first, recovered without prior knowledge.

`status`:
- `fusion`: >= 2 component chromosomes, and every partner is missing from that
  haplotype (it was absorbed into this scaffold, not duplicated).
- `candidate_partner_present`: >= 2 components but the partner chromosome is
  still present in the same haplotype. This can be a translocation, a
  duplication, or (often, in two-haplotype assemblies where one haplotype is
  fragmented) a long copy that shares repeats or ancient homeology with several
  chromosomes. Not used as fusion evidence downstream.

**`rediploidization_by_chrom.tsv`**: one row per chromosome number.
- `n_copies`, `haps_missing`: copies placed under this number, and which
  haplotypes lack one. In a fusion these are the haplotypes carrying it.
- `fusion`: confirmed fusions involving this chromosome (`with chr22 in HAP3,HAP4`).
- `group_a`/`group_b`: the 2-way split of copies by whole-chromosome distance
  (same method as `te_marker_fraction_by_lineage.tsv`; >= 3 copies needed).
- `dist_split`: mean cross-group / mean within-group whole-chromosome distance.
  ~1 = no lineage structure. `SchCurv1` chr19 = 8.1, chr17 = 1.2, most
  chromosomes 1.0-1.1.
- `te_split`: the TE-marker `split_ratio` for the same chromosome, if te-markers ran.
- `window_split_frac`, `split_extent`, `split_segments`: where along the
  chromosome the split holds. Read along each copy that has another copy in its
  own group: per window, mean distance to the other group's copies over mean
  distance to its own group's (from `windowed/`). Windows >= `--window-split`
  (1.25) count as split. Runs of split windows (one-window gaps bridged, >=
  `--min-segment-bp` 2 Mb) are segments; `split_extent` is `whole` (segments
  cover >= half the chromosome), `regional`, or `none`. A real lineage split is
  seen from every copy, so the reported reading is the median copy's
  (`window_split_frac` is the median over copies), and `split_segments` are on
  that copy's coordinates, prefixed with its haplotype (`HAP2:0.0-17.8Mb`).
  `ddEmpNigr1` chr03: `whole`, 92% of windows split.
- `split_extent_raw`, `window_covered`, `window_null_q95`: the windowed split is
  tested against a null before it counts. Window distances carry real noise, so
  a run of split windows at the 2 Mb minimum can appear by chance
  (`daBudDavi1`: single 2-3 Mb segments on chromosomes with a whole-chromosome
  split of ~1.0). The same scan is run for groupings of the copies that cut
  across the observed one (a balanced split against the other balanced splits;
  a one-copy split against the splits isolating each other copy). Their covered
  fractions, pooled over the genome, give the species' noise level.
  `split_extent` is `whole`/`regional` only when `window_covered` (fraction of
  the chromosome in split segments, median copy) exceeds both
  `window_null_q95` (that null's 95th percentile) and every crossing grouping of
  the same chromosome. Otherwise it is `none`; `split_extent_raw` keeps the
  untested reading. On the panel this removed 26 split calls in six species,
  all on short segments, and kept `SchCurv1` chr19 (0.89 vs null 0.07) and
  chr17 (0.24).
- `copy_state`: the reading for this chromosome number:
  - `fusion_lineages`: fused in some haplotypes, not others, which splits the
    copies into a fused and an unfused lineage (the snow carp mechanism, Xie et al. 2026).
  - `resolved_lineages`: balanced split, `dist_split` >= 2x threshold, and a
    whole-chromosome windowed split.
  - `partially_resolved`: >= 2 of {`dist_split` >= 1.25, `te_split` >= 2,
    windowed segment}.
  - `candidate`: exactly one of those.
  - `tetrasomic_like`: none. The copies are interchangeable, as expected under
    polysomic inheritance (or complete homogenization).
  - `one_divergent_copy`: the split isolates one copy from the rest
    (`outlier_hap` names it), on chromosome-wide evidence (distance split >=
    1.25 or TE split >= 2). A windowed segment alone does not count here: one
    of four copies carrying a divergent block in some region is ordinary
    haplotype structure in a polysomic genome (`drLytSali1`, tetrasomic by
    classical genetics, had 12/15 chromosomes in this state on regional
    segments alone; now 1/15). It has two readings, so check which:
    - One haplotype is odd on most chromosomes (`most_frequent_outlier_hap`):
      likely assembly quality, see the assembly-quality confound above
      (`ddHypMacu1` HAP1, 8/8).
    - The isolated copy changes haplotype between chromosomes, at a consistent
      magnitude: a genuinely divergent genome copy (AAAB-like; `daPilAura1`, one
      copy at split 1.3-1.7 on all 9 base chromosomes).

    With 3 copies every split is one-vs-two, so triploids land here by
    construction.

    Presence/absence differences count too. In simulations at 0.2% allelic
    divergence, giving every copy but one an 8-14% deletion made the complete
    copy the divergent one (split 1.1-1.4) and put a `regional` split over the
    missing region. A copy assembled less completely than its siblings can
    therefore read as divergent; compare lengths and `unplaced.tsv` before
    reading it as biology.
  - `not_assessable`: < 3 copies and no homeolog or partition partner to pool with.
- `state_basis`, `pooled_with`: what the state was read from:
  - `copies`: the chromosome's own >= 3 copies (all three lines of evidence).
  - `fusion`: a confirmed fusion.
  - `partition_pool`: as `homeolog_pool`, for species with no accepted homeolog
    pairs but a significant genome-wide bipartition from the `structure` stage
    (`structure/genome_partition.tsv`, k=2 row, z >= `--partition-z` (10), smaller
    group >= a quarter of the chromosomes). Each chromosome is pooled with its
    reciprocal best match in the other group (`pooled_with`); chromosomes without
    a reciprocal match stay `not_assessable`. `daInuConz1` (16/16, z=10.6),
    `dmRanRepe1` (16/16, z=10.9), `dcCerAlpi1` (32/36, z=18.9).
    **The evidence here is the partition test, not the pooled state.** Any two
    different chromosome numbers are far more distant than two alleles: pooling
    each chromosome of the diploid `ddMalSylv1` with its nearest other number
    gives `dist_split` 22-47, i.e. `resolved_lineages`. So `resolved_lineages`
    from either kind of pool means "the copies do not interchange across the two
    numbers"; it does not separate allopolyploidy from a long-diploidized
    autopolyploid, and with a partition pool it adds nothing beyond the
    partition's z. What pooling *can* show is the opposite case, a partner as
    close as the allele (`tetrasomic_like`, most of `lpElePalu1`).
  - `homeolog_pool`: the chromosome has < 3 copies (two-haplotype assemblies),
    so its copies were pooled with those of its FDR-accepted homeolog partner
    (`pooled_with`) and split as one group. A two-haplotype tetraploid has 2
    copies per number but 4 per homeolog pair. `group_a`/`group_b`/`dist_split`
    then describe the pooled group. Only whole-chromosome distance exists across
    different numbers, so the state is `resolved_lineages` (`dist_split` >= 2x
    threshold: the two numbers are separate lineages, e.g. every allo anchor,
    `daGleHede1` 3.2-4.6), `candidate` (1-2x), `tetrasomic_like` (< 1x: the
    homeolog is as close as the same-number copy, so the copies are
    interchangeable across numbers, e.g. most of `lpElePalu1`), or
    `one_divergent_copy`.
  - `none`: nothing to read it from.
- `ancient_partner`, `distance_ratio`, `ancient_state`: the FDR-accepted
  homeolog partner from `homeologs/` and its `distance_ratio`
  (`ploidy_ancestry_summary.tsv`). `ancient_state` is `paired`, `unpaired`, or
  `fusion_partner` when a detected fusion explains the pair (`SchCurv1`
  chr19<->chr22). Those pairs are fused chromosomes, not retained WGD duplicates.

**`rediploidization_summary.tsv`**: genome-level counts.
- Chromosome numbers per `copy_state`; `n_distinct_fusions` (distinct
  chromosome combinations) and `n_fused_scaffolds` (one per haplotype carrying
  one: `SchCurv1` = 1 and 2).
- `partition_pool_z`: present when chromosomes were pooled across the genome
  partition; the z of the k=2 partition used.
- `most_frequent_outlier_hap`: one haplotype isolated on many chromosomes is a
  haplotype-level issue, not per-chromosome biology.
- Ancient-paired chromosome count, with the median and coefficient of variation
  of `distance_ratio`: tight and high = one synchronized, long-finished event;
  wide spread = pairs resolving at different times.

**Caveats.** Every state is a descriptive reading of one individual's
assemblies, not a measurement of inheritance mode (that needs segregation
data). The thresholds are provisional, tuned on the snow carp anchors
(`SchCurv1`: chr19+22 fusion, chr17 regional split, the rest tetrasomic-like),
and are all CLI flags. Simulated test genomes (ROADMAP Phase 1.6) are the
planned calibration.

## `structure/` — genome-wide partitions and pair synchrony

Read from `homeologs/` (no new k-mer work).

**`pair_synchrony.tsv`**: one row. `n_accepted_pairs`, `pair_depth_median` and
`pair_depth_cv` of the accepted homeolog pairs' `mean_distance`. One
whole-genome duplication diverges every pair to about the same depth (low CV,
`daGleHede1` 0.04); several events, or pairs resolving at different times,
spread them out.

**`genome_partition.tsv`**: every k at which the chromosome numbers factor into
k groups (average-linkage clustering of the chromosome-number distance matrix)
more than random partitions with the same group sizes, at |z| >= 2. The null is
500 permutations with a fixed per-species seed, so z-scores are reproducible
for a species. They differ by permutation noise (about ±1) from the numbers
quoted from the old panel-wide script, which shared one random stream across
species.

Retired 2026-10-07 (`inheritance_metrics.tsv`): `partition_consistency`, run
length / `flip_rate`, `mean_windowed_cv`, `distance_ratio_cv` and
`combined_allo_score`. `partition_consistency` and the windowed metrics were
driven by haplotype-labelling and window-registration artefacts
(INTERPRETATION.md, "Phase 2 hand-check: corrections"), and `distance_ratio_cv`
is confounded by tetrasomic homogenization. Their last panel values are in
`meta/archive/auto_allo_spectrum.tsv`.

## `ploidyspec panel` — cross-species tables

`ploidyspec panel --results results --outdir meta --categories meta/species_categories.tsv`
collects every species' outputs. Core tables:
- `panel_summary.tsv`: one row per species: chromosome numbers and modal copies
  per chromosome; median allele (same-number) and cross-number distance; accepted
  homeolog pairs and their synchrony; the k=2 and strongest genome partition z;
  distinct fusions; and the number of chromosomes in each rediploidization
  state.
- `genome_partition.tsv`: the `structure/` partitions stacked, computed in memory
  for species that haven't run the stage.
- `rediploidization_panel.tsv`: one row per species from
  `rediploidization_summary.tsv`.

Supplementary (`supplementary/`):
- `te_markers_panel.tsv`: `te_split_median` (the within-genome contrast:
  median TE-marker `split_ratio` over chromosomes with >= 3 copies) and the
  absolute `te_marker_fraction_mean`, for reference only. A confirmed diploid
  (`ddMalSylv1`) reaches 0.25 per chromosome, above the allo anchor
  `daGleHede1`'s 0.22.
- `poly_space_features.tsv` / `poly_space_loadings.tsv` / `poly_space_pca.png`:
  an exploratory PCA over seven core numbers (allele distance, paired fraction,
  pair-depth CV, best partition z, resolved and tetrasomic fractions, TE split).
  Missing values are mean-imputed, so species with few real features
  (`n_features_present`) sit near the origin by construction.

`--categories` only colours the plot.
