// ploidyspec handbook: algorithms, outputs, interpretation and every species
// Build (from the repository root):
//   python3 docs/guide/make_handbook_species.py
//   typst compile --root . docs/guide/ploidyspec_handbook.typ

#import "handbook_lib.typ": *

#set document(title: "ploidyspec handbook", author: "ploidyspec project")
#set page(
  paper: "a4",
  margin: (x: 2.2cm, top: 2.4cm, bottom: 2.4cm),
  header: context {
    if counter(page).get().first() > 1 [
      #set text(8pt, fill: grey)
      ploidyspec handbook #h(1fr) October 2026
    ]
  },
  footer: context {
    if counter(page).get().first() > 1 [
      #set text(8pt, fill: grey)
      #h(1fr) #counter(page).display() #h(1fr)
    ]
  },
)
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 1.0em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): it => {
  pagebreak(weak: true)
  v(0.6em)
  set text(17pt, fill: accent, weight: "bold")
  it
  v(0.4em)
}
#show heading.where(level: 2): it => {
  v(0.5em)
  set text(12.5pt, fill: accent)
  it
  v(0.15em)
}
#show heading.where(level: 3): it => {
  v(0.3em)
  set text(11pt, style: "italic")
  it
}
#show raw.where(block: false): it => box(fill: rgb("#f1f3f5"), inset: (x: 2pt), outset: (y: 2pt), radius: 2pt, text(size: 9pt, it))
#show raw.where(block: true): it => block(fill: rgb("#f7f8fa"), inset: 8pt, radius: 4pt, width: 100%, text(size: 8.5pt, it))
#show figure.caption: it => text(size: 9pt, it)
#set table(stroke: 0.4pt + rgb("#cbd5e0"), inset: 5pt, align: left)
#show table: set text(size: 8.8pt)
#show table: set par(justify: false)

#let note(title: none, body) = block(
  fill: soft, inset: 10pt, radius: 4pt, width: 100%, stroke: (left: 3pt + accent),
  [#if title != none [#text(weight: "bold", fill: accent, title) \ ]#body],
)
#let caution(title: none, body) = block(
  fill: softwarm, inset: 10pt, radius: 4pt, width: 100%, stroke: (left: 3pt + warm),
  [#if title != none [#text(weight: "bold", fill: warm, title) \ ]#body],
)
// an algorithm box: numbered steps in a framed block
#let algo(title, body) = block(
  stroke: 0.6pt + accent, radius: 4pt, inset: 10pt, width: 100%, breakable: true,
  [#text(weight: "bold", fill: accent, size: 9.5pt)[Algorithm: #title] #v(-0.2em) #set text(size: 9.6pt); #body],
)
#let params(..rows) = table(
  columns: (1.6fr, 0.9fr, 3fr),
  table.header([*Parameter (flag)*], [*Default*], [*Why this value*]),
  ..rows,
)
#let stage(name, body) = box(
  fill: soft, stroke: 0.6pt + accent, radius: 4pt, inset: 6pt, width: 100%,
  [#text(weight: "bold", size: 9pt, fill: accent, name) \ #text(size: 8pt, body)],
)
#let arrow = align(center + horizon, text(14pt, fill: accent, sym.arrow.r))

// ---------------------------------------------------------------- title page
#v(3cm)
#text(30pt, weight: "bold", fill: accent)[ploidyspec]
#v(0.1cm)
#text(16pt)[Handbook: algorithms, outputs, interpretation, \ and every species so far]
#v(0.8cm)
#line(length: 100%, stroke: 1pt + accent)
#v(0.5cm)
#grid(columns: (1fr, 1fr), gutter: 0.8cm,
  note(title: "What this is")[
    The technical companion to the _ploidyspec guide_. The guide tells the
    story: why the project exists, how to run it, what we found first. This
    handbook gives each algorithm step by step, every output and how to read
    it, four worked walkthroughs on real species, how the readings relate to
    auto- and allopolyploidy and to rediploidization, and a summary of every
    species run so far. It reflects the code and results as of 10 October 2026.
  ],
  caution(title: "What changed since the guide")[
    - *Odd ploidy*: a haplotype file holding two chromosome sets is split
      into an extra copy (the AAB weevil).
    - *Homeology map*: homeology found per block, not only per chromosome.
    - *Residual tetrasomy* and *homeologous exchange* along homeolog pairs.
    - *Renumbering cycles* fixed (whitefish).
    - Summary rules made stricter (pairs, not chromosomes; minimum evidence).
    - Masu's residual tetrasomy shrank from 54 to 17.7 Mb once scaffolding
      artefacts were excluded.
    The guide's chapters 5–7 and its species appendix are superseded by
    chapters 3–9 here.
  ],
)

#pagebreak()
#outline(title: [Contents], indent: auto, depth: 2)

// ================================================================ 1
= Orientation

== The pipeline in one picture

#figure(
  grid(
    columns: (1fr, auto, 1fr, auto, 1fr, auto, 1fr),
    gutter: 4pt,
    stage("1 prepare")[manifest → chromosome copies (units); naming or automatic numbering],
    arrow,
    stage("2 kmers")[FastK table per unit at each k],
    arrow,
    stage("3 matrix")[all-vs-all distances; renumbering and extra-set check],
    arrow,
    stage("4 homeologs")[whole-chromosome homeolog pairs (FDR)],
  ) + v(6pt) + grid(
    columns: (1fr, auto, 1fr, auto, 1fr, auto, 1fr),
    gutter: 4pt,
    stage("5 windowed")[distance along each copy vs its same-number copies],
    arrow,
    stage("6 windowed-homeologs")[homeology map (blocks) → tracks between homeolog copies + controls],
    arrow,
    stage("7 structure")[genome partitions; pair synchrony],
    arrow,
    stage("8 rediploidization")[fusions, lineage states, residual tetrasomy, exchange],
  ) + v(6pt) + grid(
    columns: (1fr, auto, 1fr, 1fr),
    gutter: 4pt,
    stage("9 summary, 10 report")[four answers with confidence and evidence; HTML page],
    arrow,
    stage("optional: te-markers")[repeat k-mers private to one copy],
    stage("many species: panel")[cross-species tables],
  ),
  caption: [`ploidyspec all` runs these in order. Stage 6 runs with
  `--with-windowed-homeologs` (recommended), te-markers with
  `--with-te-markers`. Each stage reads only files written by earlier stages,
  so any stage can be re-run on its own.],
)

== Vocabulary used throughout

/ Unit: one assembled copy of one chromosome, named `HAP2_chr05`: haplotype `HAP2`, chromosome number 5.
/ Copies: the units sharing a chromosome number. Their count is set by the assembly, not inferred.
/ Same-number (allelic) distance: distance between copies of one chromosome number. In a diploid this is heterozygosity.
/ Homeologs: different chromosome numbers that descend from one ancestral chromosome through a duplication.
/ Duplicated sets: chromosome numbers that fall into groups of homeologs (pairs, triples) across most of the genome.
/ Pooling: reading a chromosome together with its homeolog partner when each number has only two copies.
/ State: the per-chromosome reading (`tetrasomic_like`, `resolved_lineages`, ...), Section 4.9.
/ Block: a stretch of a chromosome whose closest other chromosome is consistently the same one (Section 4.6).

== How to read this handbook

If you want the *ideas* first, read Chapter 2 and then the walkthroughs in
Chapter 6. If you want to *check how a number was made*, go to its stage in
Chapter 4: each section gives the algorithm as numbered steps, the parameters
with the reason for each default, and the failure it guards against. Chapter 5
says what is in every output file. Chapter 9 summarizes each species.

// ================================================================ 2
= Auto, allo and rediploidization: what the data can show

== Two axes, not one label

A polyploid genome has two copies or more of a whole chromosome set. The
classical labels sort polyploids by *origin*: an *autopolyploid* duplicates one
genome, an *allopolyploid* joins two diverged genomes. They differ in what
happens at meiosis. Autopolyploid copies are interchangeable and pair at
random (*polysomic*, in a tetraploid *tetrasomic*, inheritance); allopolyploid
copies pair with their own parent's partner (*disomic* inheritance), so the
parental sets persist as *subgenomes*.

Twyford et al. (2025) argue this is a continuum, and their key point for us is
that genomes move along it. After a duplication, an autopolyploid gradually
*rediploidizes*: chromosomes stop pairing at random, region by region, and the
copies diverge into two lineages that eventually look exactly like the
subgenomes of an allopolyploid. So ploidyspec reports two separate axes:

+ *Origin-like structure*: today, do the copies of each chromosome form one
  interchangeable pool, or separate lineages at a consistent depth?
+ *Rediploidization*: is there evidence of the transition from the first state
  to the second: where, by what mechanism, and how far it has gone?

== What each process leaves in the sequence

The table below is the conceptual core of the tool: every reading in
ploidyspec is one of these signatures.

#figure(
  table(
    columns: (1.25fr, 2.2fr, 1.5fr),
    table.header([*Process*], [*What it leaves in one individual's assembly*], [*Where ploidyspec sees it*]),
    [Recent autopolyploidy, still tetrasomic],
    [All copies of a chromosome about equally close; no lineage structure; repeat content shared.],
    [`copy_state = tetrasomic_like`; pooled pairs tetrasomic-like],

    [Allopolyploidy (two diverged parents)],
    [Copies fall into two groups at a consistent depth on every chromosome; each group carries its parent's repeat bursts.],
    [`resolved_lineages` throughout; homeolog pairs with low `pair_depth_cv`; TE split],

    [Rediploidization by fusion (snow carps)],
    [Some copies of a chromosome fused to another chromosome; fused and unfused copies form lineages; divergence highest at the fusion.],
    [`fusions.tsv` (`status = fusion`); `fusion_lineages`; `lineage_divergence`],

    [Rediploidization by inversion or gradual divergence],
    [A lineage split along part of a chromosome only.],
    [`partially_resolved` with `split_extent = regional`],

    [Residual tetrasomy (salmonids)],
    [In an old, mostly rediploidized duplication, a few regions where homeologs still pair and stay as close as alleles.],
    [`residual_*` columns; `residual_tetrasomy.tsv`],

    [Homeologous exchange (allopolyploids)],
    [One copy carries its homeolog's sequence over a stretch: closer there to the homeolog than to its own homolog.],
    [`exchange_*` columns],

    [Karyotype change since the duplication],
    [A chromosome whose arms are duplicated on two different chromosomes.],
    [homeology map: `multi_partner_list`],

    [Asynchronous resolution],
    [Homeolog pairs or blocks diverged to different depths.],
    [`pair_depth_cv`; `block_dist_cv`],

    [An extra divergent set (AAB, AAAB)],
    [One copy of every chromosome stands apart from the others.],
    [`one_divergent_copy`; extra-set split (Section 4.3)],
  ),
  caption: [Processes, their signatures, and the outputs that read them.],
)

== What one individual cannot tell you

Four limits shape every answer the tool gives.

+ *A rediploidized autopolyploid looks allopolyploid.* Once the lineages have
  separated everywhere, nothing in today's sequence records whether they
  started as one genome or two. ploidyspec says "allo-like" and adds that a
  long-rediploidized autopolyploid looks the same. Only leftovers of the
  transition (residual tetrasomy, fusions, a mix of states) point to an
  autopolyploid origin.
+ *Divergence records isolation, not duplication* (Gaynor et al., 2026):
  copies that still recombine stay similar. A tetrasomic autopolyploid leaves
  little to see beyond equal distances.
+ *Exchange and residual tetrasomy are mirror images.* Both make homeologs
  locally near-identical. The difference is which copies agree: in residual
  tetrasomy all copies are close to each other; in an exchange one copy carries
  its homeolog's sequence and is far from its own homolog. ploidyspec now
  tests both (Section 4.10) because the first version confused them.
+ *k-mers see homeology only up to about 8–10% divergence.* Glechoma's
  homeologs (0.042 apart) are mapped over 78% of the genome; the ~90-My-old
  salmonid duplication only over 6–10% (the least-diverged regions).
  "Not detected" never means "not there".

#caution(title: "Assembly artefacts look like biology")[
  Missing sequence reads as divergence; a contig placed on the wrong homeolog
  reads as exchange or residual tetrasomy; one haplotype file numbering its
  chromosomes differently reads as a broken genome. Much of this handbook's
  algorithmic detail is about guarding against these. The masu walkthrough
  (Section 6.4) shows what happens when a guard is missing.
]

== How the four answers are built

#figure(
  table(
    columns: (1fr, 2.6fr),
    table.header([*Question*], [*Built from*]),
    [Ploidy], [copies per number (assembly) × duplicated sets (homeolog pairs, genome partition, or homeology map), with odd copy numbers flagged],
    [Origin-like structure], [per-chromosome states; for two-copy genomes the pooled homeolog pairs, residual tetrasomy (≥ 2 pairs) and the partition; one-copy-apart patterns],
    [TE markers], [the within-genome contrast `te_split`, never the absolute fraction],
    [Rediploidization], [fusions; mixtures of states; regional splits; residual tetrasomy; exchange; rearrangements and block-divergence spread from the map; pair synchrony],
  ),
  caption: [Section 4.11 gives the full decision logic.],
)

// ================================================================ 3
= Measurement primitives

Every number in ploidyspec comes from comparing sets of k-mers (words of k
bases). There is no alignment, no annotation and no reference genome.

== k-mer tables

FastK (Myers) counts every canonical k-mer of a sequence (a k-mer and its
reverse complement count as one) into a sorted table. `kmers` builds one table
per unit and per k; `Logex` intersects two tables; `Histex` counts distinct
k-mers; `Profex` reads a per-position profile (below).

== Whole-chromosome distance: Mash with a noise floor

#algo("pairwise distance between two units A and B")[
  + For each k in the sweep (default 15 and 23): $n_A$, $n_B$ = distinct k-mers;
    $s$ = shared k-mers; Jaccard $J = s \/ (n_A + n_B - s)$.
  + Mash distance $d = -1/k ln(2J \/ (1 + J))$, the per-base substitution rate
    under which a k-mer survives with probability $(1-d)^k$ (Ondov et al., 2016).
  + Chance floor: two unrelated sequences share $E = n_A n_B \/ (4^k \/ 2)$
    k-mers by chance. Resolution $z = (s - E) \/ sqrt(E)$; k is usable if $z ≥ 5$.
  + Report $d$ at the *largest usable k* (most specific). If no k is usable,
    the pair is `resolution_limited` and reported without a distance.
]
On every one of the 53,309 chromosome pairs of the first panel, k = 23 was
usable, so k = 23 is the working scale. Typical values: alleles 0.0001–0.04,
homeologs 0.04–0.10, unrelated chromosomes 0.07–0.16.

== Window distance: containment

For a window W of one copy and a whole other copy B, containment
$c$ = fraction of W's k-mer *positions* whose k-mer occurs anywhere in B, and
$d = -ln(c) \/ k$ (1.0 if nothing is shared). Containment is asymmetric and
position-free: it asks whether W's sequence exists anywhere in B, not at the
same coordinate.

== Position-free windows

#algo("windowed track of copy A against copy B")[
  + Profile A's whole sequence against B's k-mer table (`FastK -p:B`): for every
    k-mer position of A, 1 if its k-mer is in B.
  + Profile A against its *own* table: the valid positions (Ns and gaps drop out
    of the denominator instead of reading as divergence).
  + `Profex -z` streams runs of present positions; an awk program sums them into
    fixed bins of gcd(window, step) positions, so the per-position profile never
    reaches Python.
  + Sum bins into windows (default 250 kb, step = window) along A's own
    coordinates; $c$ = present / valid; $d = -ln(c) \/ k$.
]
This replaced comparing windows at equal coordinates, which fell out of
register after the first indel larger than a window (the guide's Section 8.1).
Memory is about 14 bytes per base of B per concurrent profile.

== Statistical tools used repeatedly

- *Robust z*: (value − median) / (1.4826 × MAD), with the MAD floored. Used for
  fusion containment and homeology-map windows; resistant to the very outliers
  being looked for.
- *Shuffle nulls*: labels (near-allelic, exchange, same-partner) are shuffled
  across the whole genome and the longest run recorded, 100–500 times; a real
  run must beat the 95th percentile. This adapts automatically to each
  species' noise level and window density.
- *Floors on nulls*: when a genome is nearly homozygous, near-allelic windows
  are rare and the shuffled maximum falls to 3–4 windows, so short clusters
  would pass. Run nulls are floored (8 windows by default).
- *Benjamini–Hochberg FDR* for homeolog pairs (Section 4.4).
- *Votes across copies*: real lineage structure or residual tetrasomy is seen
  from every copy; an assembly error is usually seen from the one or two
  copies it touches.

// ================================================================ 4
= The stages, algorithm by algorithm

== prepare: from a manifest to chromosome copies

#algo("prepare")[
  + Read the manifest: FASTA path and haplotype label (`AUTO` = read the
    haplotype from headers with `--hap-regex`).
  + Index each FASTA (`samtools faidx`; plain gzip must be recompressed with
    bgzip). Sequences shorter than `--min-len` (1 Mb) go to `unplaced.tsv`.
  + Unlocalised sequences (`unloc`, `_random`) are excluded: NCBI names them
    "chromosome N" too, which would give two units the same number.
  + Chromosome number from the header (`chromosome: N`, `SUPER_N`, or
    `--chrom-regex`). `--chrom-naming detect` falls back to automatic numbering
    when no header has a number.
  + *Automatic numbering*: per haplotype, keep sequences ≥ 10% of the median
    length of the longer half (chromosome-scale); the reference haplotype (most
    such sequences, excluding contig-level ones) is numbered by length; the
    others get provisional numbers (≥ 1001) that the matrix stage replaces by
    their closest reference chromosome.
  + Write `sequences.tsv` (one row per unit) and `unplaced.tsv`.
]

== kmers

Extracts each unit to `chroms/` and builds its FastK table per k
(`ktabs_kNN/`). Each table records its source file, size and modification
time, and the sequence name and length; any change rebuilds it (a re-run on a
renumbered assembly once silently reused old sequences).

== matrix: distances, renumbering, extra sets

The all-vs-all distance matrix (Section 3.2), then three checks on chromosome
numbering, all using that matrix and all logged in
`matrix/chrom_label_corrections.tsv`. Every check renames the unit's files on
disk, so all later stages see the corrected numbers.

#algo("renumbering between haplotype files")[
  + Reference = the source file with most units. Only other files are ever
    changed (correcting both sides of a swap independently can merge two
    chromosomes).
  + For each other unit: distance to its declared number's reference copies vs
    to the best other number. Ratio ≥ 1.5 makes it a candidate.
  + Match the file's units one-to-one to reference numbers, closest first
    (two units cannot claim one number; `ddLepDrab1` HAP4 needed this).
  + Apply a move when the matched number is ≥ 3× closer than the declared one;
    a mutual swap needs only ≥ 2× on both sides (`ddHesMatr1`, 3.1× and 2.3×).
  + A move into a slot whose occupant stays would duplicate a unit name. If the
    occupant is no closer to its own label than to the slot the matching gave
    it, it is *displaced* there (the label had no support). Otherwise the whole
    file's batch is rejected as `ambiguous`, never partly applied.
]
The displacement rule was added for the European whitefish (`fCorLav1`): its
HAP2 numbering differs from HAP1 on 30 of 40 chromosomes in long cycles (7 → 9
→ 10 → 11 → 13 → 12 → 7, at 10–32× evidence), but HAP2 chr33 has no allele in
HAP1 at all (0.11–0.12 to everything), stayed put, blocked the cycle through
chr33 and the whole correction was rejected. With displacement, 30 corrections
plus 1 displacement apply cleanly.

#algo("extra chromosome set in one haplotype (odd ploidy)")[
  + Private numbers: numbers present in one haplotype only.
  + For each private unit: mean distance to every *shared* number's copies;
    background = the median of those.
  + Match private units one-to-one to shared numbers, closest first; keep a
    match if background / distance ≥ 1.3.
  + If at least 3 units and at least half the private units match, they are an
    extra set: they become a new haplotype (`HAP1` → `HAP1B`) carrying the
    matched numbers (`status = extra_set`).
]
The weevil `icStrMela3` (Section 6.3) is the case: HAP1 holds A as chr1–10 and
B as chr11–20; each B chromosome is 0.047–0.066 from one A number against a
median 0.093–0.099 to the rest (1.5–2.4×). A single private chromosome (the
Monacha fusion leftover) is never a set.

== homeologs: whole-chromosome pairs

#algo("homeolog pairs")[
  + For every pair of chromosome numbers (i, j): the mean distance over all copy
    combinations, so each number is one node whatever its copy count.
  + Fit the background: mean and SD of all pair distances, σ-clipping low values
    (2σ, 3 rounds) so a genome full of homeologs does not hide its own signal.
  + z = (d − mean) / SD; one-sided normal p; Benjamini–Hochberg q.
  + Accept pairs with q ≤ 0.05 *and* d ≤ 0.98 × median background (`--min-effect`
    2%), greedily one partner per number, closest first.
]
The 2% effect floor exists because a tight background makes tiny differences
significant: on simulated diploids, a pair 1.25% below the median reached
z = −4. Outputs: `homeolog_pairs.tsv`, every candidate in
`homeolog_candidates_ranked.tsv`, and per chromosome in
`ploidy_ancestry_summary.tsv` (`distance_ratio` = homeolog distance / own
allelic distance).

== windowed: along each chromosome, same number

For each chromosome number, every copy is windowed against every other copy of
that number (Section 3.4, k = 15). These tracks carry the lineage splits along
chromosomes (Section 4.9) and each copy's own-homolog distance per window,
which the residual and exchange tests need (Section 4.10).

== windowed-homeologs: the homeology map, then pair tracks

This stage first maps homeology genome-wide, then computes windowed tracks
between the copies of every homeolog pair (whole-chromosome pairs plus block
pairs from the map), plus control tracks.

#algo("homeology map")[
  + Take the haplotype with the most chromosomes. Window every chromosome
    (250 kb) against every *other* chromosome of that haplotype, at k = 23
    (Section 3.4; n² profiles).
  + Per window: the closest chromosome, the median distance to all others, and
    robust z = (median − closest) / (1.4826 MAD). Homeolog-like if z ≥ 3.
  + Label each window with its closest chromosome if homeolog-like, else none.
    A *block* is a run of one label, bridging up to 2 windows of anything else.
  + Null: shuffle the labels over the whole genome (shuffling within a
    chromosome would keep a chromosome homeologous end to end as long as
    itself); a block must be longer than the 95th percentile of the longest
    shuffled run, at least 6 windows and `--min-segment-bp` (2 Mb).
  + *Reciprocity*: a block A → B counts only if B has a block back on A.
  + Pairs of numbers with ≥ 2 Mb of reciprocal blocks join the pair list;
    chromosomes with reciprocal blocks on two or more partners are reported.
]
Why k = 23: at k = 15 about 40% of a window's k-mers occur by chance in any
100 Mb chromosome, and salmonid homeologs barely stand out (masu 0.049 vs 0.057
for unrelated chromosomes). Why reciprocity: homeology runs both ways, repeats
need not. Without it, Glechoma gained two short blocks onto chr16 at a third of
the homeolog distance (shared repeats) and the diploid _Malva sylvestris_ had
seven blocks; with it, Glechoma has exactly its nine known pairs (78% of the
genome) and _Malva_ none.

#algo("pair tracks and controls")[
  + For each pair (i, j): window every copy of i against every copy of j and the
    reverse (k = 15), giving per window the distance to each partner copy.
  + Control: each paired chromosome against one unrelated chromosome (a paired
    number that is none of its partners), same windows.
]

== structure: partitions and synchrony

#algo("genome partition")[
  + Cluster chromosome numbers by average linkage on the number-to-number
    distances, recording every k from n groups down to 2.
  + Separation = mean between-group / mean within-group distance.
  + Null: 500 random groupings with the same group sizes (seeded per species).
    z = (observed − null mean) / null SD; report k with |z| ≥ 2.
]
A strong balanced two-group partition (`daInuConz1` 8 vs 8, z = 10.6) is a
candidate subgenome structure that no single pair test finds; the
rediploidization stage pools chromosomes across it (Section 4.9).
`pair_depth_cv` = SD / mean of the accepted pairs' distances: one duplication
event diverges every pair to about the same depth (`daGleHede1` 0.042).

== te-markers (optional)

For each pair of copies of a chromosome, at k = 13: high-copy k-mers (count
≥ 100, mostly repeats) and markers (high-copy k-mers ≥ 2× more abundant in one
copy). `te_marker_fraction` = markers / high-copy k-mers of both copies.
For chromosomes with three or more copies, the copies are split in two by
distance and `split_ratio` = cross-group fraction / within-group fraction. The
absolute fraction does not separate diploid, auto and allo (the diploid _Malva_
spans 0.05–0.25); the within-genome split does.

== rediploidization: fusions and per-chromosome states

=== Fusions

#algo("fusion detection")[
  + Candidates: placed copies ≥ 1.4× the median length of their sibling copies,
    and unplaced scaffolds ≥ half the median chromosome length.
  + For a candidate, containment of each *other* chromosome's best copy in it
    (shared k-mers / that copy's k-mers).
  + Background: containment between every pair of ordinary units of different
    numbers; robust z.
  + A component needs z ≥ 10, ≥ 2× the background median, and ≥ ¼ of the
    strongest hit.
  + Placed copies only: the partner must be ≥ 2.5× more enriched in the
    candidate than in its least-enriched unfused sibling (`sibling_excess`).
    Ancient homeology is shared by every copy; a fusion is not.
  + `status = fusion` when every partner number is *missing* from that
    haplotype (absorbed, not duplicated); otherwise
    `candidate_partner_present`.
  + `lineage_divergence` = $-ln(c) \/ k$ of the partner containment: how long the
    fused lineage has been diverging.
]
A fusion on one haplotype in a genome with no other duplication (the land snail
_Monacha_, chr23+chr24 joined in HAP1 only) is a fusion polymorphism or a
scaffolding join, not rediploidization; the summary says so.

=== Lineage split within a chromosome number

#algo("split, extent and null (three or more copies)")[
  + Bipartition: seed with the two most distant copies, assign the rest to the
    closer seed. `dist_split` = mean cross-group / mean within-group distance;
    *balanced* if both groups have ≥ 2 copies.
  + Along every copy that has a partner in its own group: per window, mean
    distance to the other group / mean distance to its own group. Split windows
    have ratio ≥ 1.25; segments are runs (1-window gaps bridged) ≥ 2 Mb.
  + The *median copy*'s covered fraction is the reading (a real split shows from
    every copy). `whole` if ≥ 50%, `regional` if any segment.
  + Null: repeat for groupings that *cut across* the observed one (other
    balanced splits; or isolating each other copy), pooled genome-wide. The
    extent counts only if it beats that null's 95th percentile and every
    crossing grouping of its own chromosome.
]

=== Pooling for two-copy genomes

With two copies per number nothing can be split. The chromosome is pooled with
its whole-chromosome homeolog partner (four copies), or else with its
reciprocal best match across a significant (z ≥ 10, balanced) two-group
partition. Pooled groups are read on whole-chromosome distance only.

=== States

#figure(
  table(
    columns: (1.15fr, 2.7fr),
    table.header([*State*], [*Rule*]),
    [`fusion_lineages`], [the chromosome is a component of a detected fusion],
    [`resolved_lineages`], [balanced split, `dist_split` ≥ 2.5 and extent `whole`; pooled: `dist_split` ≥ 2.5 (any two different numbers are far apart, so for pooled pairs this mostly says "these are homeologs, not alleles")],
    [`partially_resolved`], [two or more of: `dist_split` ≥ 1.25, `te_split` ≥ 2, a windowed segment; pooled: resolved but carrying residual tetrasomy],
    [`candidate`], [exactly one of those],
    [`tetrasomic_like`], [none: the copies are interchangeable],
    [`one_divergent_copy`], [one copy apart, with chromosome-wide (distance or TE) evidence; a windowed segment alone does not count],
    [`not_assessable`], [fewer than three copies and nothing to pool with],
  ),
  caption: [Per-chromosome states (`rediploidization_by_chrom.tsv`).],
)

== Residual tetrasomy and homeologous exchange

Both read the pair tracks of Section 4.6 against each copy's own-homolog track
from Section 4.5, for every homeolog pair (whole-chromosome or block).

#algo("residual tetrasomy")[
  + *Allelic level* per copy: median over windows of the distance to its closest
    own-number copy; per chromosome the median over its copies; capped at 1.5×
    the genome median (a chromosome carrying a large exchange must not raise its
    own yardstick).
  + Per window of a copy, *near-allelic* when all of:
    - closest partner copy < 3 × max(allelic level, 0.001);
    - closest partner copy < 0.25 × the species' median homeolog distance;
    - the copy's *own* homolog is also within 3 × allelic level (in residual
      tetrasomy every copy is close to every other);
    - the unrelated control is not as close (shared repeats are close to every
      chromosome carrying them).
  + Runs of near-allelic windows (1-window gaps bridged) count when longer than
    the 95th percentile of the genome-wide longest run under shuffling (200
    permutations), floored at 8 windows, and ≥ 2 Mb.
  + *Pair vote*: the pair has residual tetrasomy when more than half the copies
    of its two numbers show a counted run (3 of 4 with two haplotypes).
  + Per chromosome: the best-covered copy's segments, summed over partners;
    segments within 10% of an end are `start`/`end`.
]
#algo("homeologous exchange")[
  + Per window of a copy: *exchange-like* when the closest partner copy is less
    than half the distance to the copy's own homolog, and the own homolog is
    beyond 3 × allelic level.
  + Runs against their own shuffle null (same floor and length).
  + No pair vote: an exchange is usually carried by one copy;
    `exchange_support` says how many.
]
#caution(title: "Why the own-homolog condition matters")[
  The first version counted any window close to a homeolog copy. An exchange in
  one copy then passed the pair vote from the three copies that see it (the
  simulated `segmental_homeology` scenario failed this way), and so did contigs
  swapped between homeologs in masu's reference-scaffolded HAP2. 6 of masu's 8
  earlier residual chromosomes rested on windows whose own homolog was
  0.03–0.085 away, homeolog distance, not missing sequence.
]

== summary: the decision logic

The summary reads the files above and answers four questions, each with a
confidence and the evidence it used. In order of evaluation:

*Ploidy.* Copies = the modal copy number per chromosome (high confidence: it is
what the assembly contains). Duplicated sets, first match wins: a genome
partition into groups of ≥ 2 numbers (k ≥ n/4, z ≥ 10); homeolog pairs covering
≥ 50% of the numbers; a balanced k = 2 partition; reciprocal homeology-map
blocks covering ≥ 30% of the genome with ≥ 2 pairs (then x is not read from the
numbers, since the sets are reshuffled). With ≥ 3 pairs but no sets, an older
duplication "shows on part of the genome". Odd copy numbers and split-out extra
sets are stated.

*Origin-like structure.*
+ A fusion on one haplotype with no other duplication: not assessable.
+ Two copies and no duplication evidence (no sets and fewer than 3 pairs by
  either method): not assessable (gar and bowfin pass one pair each).
+ One copy apart on most chromosomes: the same haplotype every time is read as
  an assembly or phasing problem, unless that haplotype is a split-out extra
  set, which reads as AAB-like (high); changing haplotype reads as A…AB-like.
+ Two-copy (pooled) genomes: tetrasomic-like on most pairs → auto-like;
  residual tetrasomy on ≥ 2 independent pairs → "auto-like origin, mostly
  rediploidized"; ≥ 60% split → allo-like (high if pair-depth CV ≤ 0.06).
+ Three or more copies: mostly tetrasomic → auto-like; mostly split →
  allo-like; a mixture → mixed.

*TE markers.* Only `te_split` is read: ≥ 2 on chromosomes that also split by
distance means separate repeat histories.

*Rediploidization.* Fusions (between copies of a duplicated genome) → yes
(high); mixtures of split and tetrasomic-like chromosomes → likely under way;
two-copy genomes with residual tetrasomy → partly; no residual tetrasomy but
chromosomes with two or more homeolog partners → advanced; otherwise diverged
throughout. Exchange candidates, rearrangements and block-divergence spread
(CV ≥ 0.15) are appended as "Also:". Residual tetrasomy ≥ 90% at chromosome
ends is low confidence (homeolog-specific subtelomeric repeats pass the
control).

== report, panel, simulate, cleanup

`report` writes one self-contained HTML page per species. `panel` stacks all
species into `panel_summary.tsv`, `rediploidization_panel.tsv`,
`species_summaries.tsv`. `simulate` writes synthetic assemblies with a known
answer (Chapter 7). `cleanup` deletes k-mer tables and keeps every result.

// ================================================================ 5
= What every output means

`OUTPUTS.md` in the repository documents every column. This chapter is the
map: what each file is for and the columns you will actually read.

#table(
  columns: (1.6fr, 2.8fr),
  table.header([*File*], [*What to read in it*]),
  [`sequences.tsv`], [one row per unit: `unit_id`, `hap`, `chrom`, `seq_id`, `length`, `source`; `numbering` (`names`/`auto`)],
  [`unplaced.tsv`], [every excluded sequence and why (`below-min-len`, `unlocalised`, `no-chrom-match`, ...)],
  [`matrix/whole_chrom_pairs.tsv`], [per unit pair: `chosen_k`, k-mer counts, `jaccard_at_chosen_k`, `distance`, `resolution_limited`],
  [`matrix/ploidy_summary.tsv`], [copies per chromosome number and which haplotypes],
  [`matrix/chrom_label_corrections.tsv`], [`status` (`corrected`, `displaced`, `extra_set`, `ambiguous`), old/new number, the distances and ratio, `new_hap`],
  [`matrix/*heatmap*.png`], [the full matrix; the contrast version shows block structure],
  [`homeologs/homeolog_pairs.tsv`], [accepted pairs: `mean_distance`, `z_score`, `q_value`],
  [`homeologs/ploidy_ancestry_summary.tsv`], [per chromosome: partner, `distance_ratio`],
  [`homeologs/homeology_map.tsv`], [per window: `best`, `best_dist`, `median_dist`, `z`, `homeolog_like`],
  [`homeologs/homeolog_blocks.tsv`], [per block: `chrom`, `partner`, `start`–`end`, `median_dist`, `position`, `reciprocal`],
  [`homeologs/homeology_map_summary.tsv`], [`duplicated_frac`, `n_partner_pairs`, `multi_partner_list`, `block_dist_cv`],
  [`homeologs/homeology_map.png`], [each chromosome painted by its partner; hatched = unmatched],
  [`homeologs/windowed_chrAAxBB.png`], [window distances between the copies of a homeolog pair],
  [`windowed/windowed_all.tsv`], [per window of `unit_a` vs `unit_b` (same number): `kmers_a`, `shared`, `distance`],
  [`structure/genome_partition.tsv`], [partitions with `group_sizes`, `separation_ratio`, `z_score`, `groups`],
  [`structure/pair_synchrony.tsv`], [`pair_depth_median`, `pair_depth_cv`],
  [`subgenomes/te_marker_fraction_by_lineage.tsv`], [`split_ratio` per chromosome],
  [`rediploidization/fusions.tsv`], [per tested scaffold: `components`, `containment`, `sibling_excess`, `lineage_divergence`, `status`],
  [`rediploidization/rediploidization_by_chrom.tsv`], [per chromosome: groups, `dist_split`, `te_split`, `split_extent`, `split_segments`, `copy_state`, `state_basis`, `pooled_with`, residual and exchange columns, `block_partners`],
  [`rediploidization/residual_tetrasomy.tsv`], [per copy and partner: windows, near-allelic and shared-with-control counts, `allele_level`, runs vs null, segments, exchange segments],
  [`rediploidization/rediploidization_summary.tsv`], [genome-level counts, residual and exchange totals, `map_*` metrics],
  [`summary.tsv`, `summary.md`], [the four answers, confidence, evidence],
  [`report.html`], [everything on one page],
)

#note(title: "Reading a per-chromosome row in practice")[
  Start with `copy_state` and `state_basis` (`copies`, `homeolog_pool`,
  `partition_pool`, `fusion`). For `copies`, look at `dist_split` and
  `split_extent`; for a pool, at `pooled_with` and `dist_split`. Then
  `residual_support` (e.g. `chr17:3/4` = 3 of the pair's 4 copies show a
  counted run), `residual_segments` (`HAP1:16.8-19.5Mb(interior)~chr17` = on
  HAP1, 16.8–19.5 Mb, mid-chromosome, against chr17), and the exchange columns
  the same way.
]

// ================================================================ 6
= Walkthroughs

== Glechoma hederacea: a two-haplotype allotetraploid

*Input.* Two haplotypes, 18 chromosomes each, numbered by name.

*Matrix.* HAP1 chr02 vs HAP2 chr02: 0.0128 (alleles). HAP1 chr02 vs HAP1 chr03:
0.0474. HAP1 chr01 vs HAP1 chr02: 0.113 (unrelated). So chr02 and chr03 are
four times further apart than alleles but less than half as far as unrelated
chromosomes: homeologs.

*Homeologs.* Nine pairs, all at q ≈ 0 and z between −18.6 and −20.7:
chr02–03 (0.047), chr01–05, chr04–10, chr06–09, chr07–08, chr12–13, chr11–15,
chr16–18, chr14–17 (0.054). `distance_ratio` 3.5 (median). Pair-depth CV 0.042:
one synchronous event. The genome partition into nine pairs has z = 23.8.

*Homeology map* (Figure below). 78% of the genome in 39 reciprocal blocks, the
same nine pairs, no chromosome with two partners, block CV 0.047.

#figure(image("../../results/daGleHede1/homeologs/homeology_map.png", width: 100%),
  caption: [Glechoma: every chromosome is duplicated on exactly one partner
  along most of its length. Grey gaps are regions where the homeolog is not
  detected (diverged beyond k-mer reach, or lost).])

*Rediploidization stage.* Each chromosome is pooled with its partner (four
copies); every pool splits by `dist_split` 3.5–4.6, so 16 chromosomes are
`resolved_lineages`. The residual test finds one pair, chr14/chr17: 2.75 Mb
on HAP1 chr14 (16.8–19.5 Mb) and 6.75 Mb on HAP2 chr17 (16.0–22.8 Mb),
mid-chromosome, support 3 of 4. Both chromosomes become `partially_resolved`.

#figure(image("../../results/daGleHede1/homeologs/windowed_chr14x17.png", width: 100%),
  caption: [chr14 × chr17, all eight copy directions. Outside ~14–23 Mb the
  homeologs are 0.03–0.05 apart; inside they drop to ~0.003, as close as
  alleles. A mid-chromosome block like this could be retained pairing, a
  recent homeologous exchange shared by both haplotypes, or a centromeric
  repeat specific to this pair: the data cannot tell which.])

*Summary.* Ploidy: 2 copies × 18 numbers in sets of 2, x ≈ 9, up to 4x.
Origin: allo-like (disomic-like), high: one residual pair is not enough to call
an auto-like origin (the rule needs two independent pairs). Rediploidization:
partly, 9.5 Mb of residual tetrasomy on one pair, plus one exchange candidate.

*What to take away.* Glechoma is the textbook allo-like genome: a single
synchronous duplication, every chromosome paired, and almost no shared
stretches. The one near-identical block is worth a look but does not change the
reading.

== Schizothorax curvilabiatus: rediploidization by fusion

*Input.* Four haplotypes × 25 chromosome numbers (an autotetraploid snow carp).
HAP3 and HAP4 chr19 are 69.5 Mb, against 37.8 and 40.2 Mb for HAP1 and HAP2
chr19; chr22 is missing from HAP3 and HAP4.

*Fusions.* HAP3 chr19 and HAP4 chr19 are long-copy candidates. Each contains
chr22's k-mers far above background (containment 0.452 and 0.437, z ≈ 48), 3.4×
more than their unfused siblings do; chr22 is absent from both haplotypes:
`status = fusion`, components chr19+chr22, `lineage_divergence` 0.035–0.036.

*Lineages.* chr19's copies split HAP1+HAP2 vs HAP3+HAP4 at `dist_split` 8.1
(0.007 within, 0.055 across), `te_split` 5.6, and along the whole chromosome
(89% of windows, against a genome null of 7%): `fusion_lineages`.
chr17 splits only over its last 8 Mb (24% covered, beats the null),
`dist_split` 1.21, `te_split` 3.2: `partially_resolved`, regional. Most
chromosomes (17 of 25) are `tetrasomic_like`.

#figure(image("fig/snowcarp_windows.png", width: 100%),
  caption: [`SchCurv1` windows along HAP1. chr19: the fused copies stay far all
  along. chr17: the copies separate only towards the end.])

*Summary.* 4 copies × 25, no older sets. Origin: mixed (tetrasomic-like and
split chromosomes). Rediploidization: yes, one fusion between copies (the snow
carp mechanism), high. This is the published result of Xie et al. (2026),
recovered without being told where to look.

== Strophosoma melanogrammum: an AAB triploid in two files

*Input.* HAP1 has 20 chromosomes, HAP2 has 10. The ToL k-mer ploidy plot says
triploid, AAB.

*Matrix.* chr1–10 are shared; chr11–20 are private to HAP1. Each private
chromosome is closest to one shared number: HAP1 chr11 to chr01 at 0.049 against
a median 0.098 to the other numbers (1.97×), chr12 to chr03 (1.99×), chr13 to
chr04, ... chr19 to chr09 (2.37×). Ten matches out of ten: an extra set.
HAP1 chr11–20 become `HAP1B_chr01`–`HAP1B_chr10`.

*Copies.* Every chromosome now has three copies. On chr01, HAP1 and HAP2 are
0.011 apart (the two A copies); HAP1B is 0.049 from both (B). The bipartition
isolates HAP1B with `dist_split` 4.6 and `te_split` 3.1: `one_divergent_copy`,
odd haplotype HAP1B, on 9 of 10 chromosomes.

*Summary.* Ploidy: 3 copies × 10 numbers, HAP1B split out, odd copy number.
Origin: AAB-like, one divergent chromosome set, high. Rediploidization: not
assessable (a divergent set is not a lineage that split off).

*The odd one out.* On chr09 the outlier is HAP2, not HAP1B: HAP1 chr09 sits
between A and B (0.030 to HAP2 chr09, 0.027 to B). A homeologous exchange, or a
phasing mix-up in HAP1 chr09; worth checking in the assembly.

== Oncorhynchus masou: when a guard is missing

*Input.* Two haplotypes × 33. HAP2 was contig-level and was scaffolded on HAP1
for this project (`workflows/prep/scaffold_by_reference.py`): its contigs were
placed where they align best on HAP1.

*The first result.* Residual tetrasomy on 8 of 10 paired chromosomes, 54 Mb,
62% terminal: an apparently textbook salmonid. Lien et al. (2016) describe such
regions in Atlantic salmon.

*What was wrong.* In a region where two homeologs are similar, a HAP2 contig
from chromosome B aligns almost as well to A and can be placed there. HAP1 A is
then near-identical to "partner" HAP2 B (which carries A's sequence), and far
from its own HAP2 A (which carries B's). Three of the four copies see a
near-allelic stretch and the pair vote passed. The tell: those windows' own
homolog distance was 0.03–0.085, homeolog distance.

*The current result.* With the own-homolog condition: residual tetrasomy on 4
chromosomes in 2 pairs (chr17/chr29 and chr28/chr33, both seen from all 4
copies), 17.7 Mb, 70% at chromosome ends; exchange candidates on 5 chromosomes, each in
1 of 2 copies (the swapped contigs, now read as what they are). The homeology
map finds 10% of the genome in blocks with two chromosomes carrying two
partners (chr10: chr12+chr13; chr28: chr05+chr33).

*What to take away.* Our own preprocessing produced the strongest
"biological" signal in the panel. Masu remains a useful salmonid, but its
residual numbers carry that history; Arctic charr and whitefish, assembled
with both haplotypes scaffolded, are the cleaner salmonid tests, and the
k-mer method only reaches their least-diverged 6–17%.

#figure(image("../../results/fSalAlp3/homeologs/homeology_map.png", width: 82%),
  caption: [Arctic charr: the salmonid duplication is mostly beyond k-mer reach.
  Reciprocal blocks (solid) cover 6% of the genome. Hatched blocks on chr12 and
  chr24 point at chr01 (the largest chromosome, 117 Mb) without a block back:
  consistent with chr01 being a fused chromosome whose arms are homeologous to
  chr12 and chr24, but too weak from chr01's side to count.])

// ================================================================ 7
= Validation

== Simulated genomes

`ploidyspec simulate` writes random-sequence genomes with transposable-element
families at high copy number; every haplotype except HAP1 carries a
transposition and an inversion, so copies are not collinear. The end-to-end
test (21 tests) runs the whole pipeline on each scenario.

#table(
  columns: (1.3fr, 2.7fr, 1.4fr),
  table.header([*Scenario*], [*What it contains*], [*Checked*]),
  [diploid], [two haplotypes], [no pairs, no blocks],
  [autotetraploid], [four equivalent haplotypes], [all `tetrasomic_like`, no blocks],
  [autotetraploid, 2 hap], [the same, each base chromosome numbered twice], [pairs (i, i+n), pooled `tetrasomic_like`],
  [allotetraploid], [two progenitors, own TE bursts, 2 haplotypes], [`resolved_lineages`, no residual tetrasomy],
  [rediploidized], [whole and 40% lineage splits, a placed and an unplaced fusion], [states and both fusions],
  [mislabelled], [one file shifted by one, a swap in another], [all relabelled],
  [residual tetrasomy], [diverged homeologs except the tails of two pairs; a satellite at every chromosome start], [residual on exactly those four chromosomes; satellite excluded by the control; "auto-like origin"],
  [allotriploid, one file], [AAB with A + B in HAP1's file], [HAP1B split out; AAB-like],
  [segmental homeology], [a fused chromosome whose arms are homeologous to two others; an exchange in one haplotype], [two partners on chr01; exchange on HAP1 chr04 only],
)

== Real anchors

#table(
  columns: (1.1fr, 1.6fr, 2.3fr),
  table.header([*Anchor*], [*Known*], [*ploidyspec*]),
  [`daGleHede1`], [allopolyploid], [9 pairs, CV 0.042, allo-like (high); map: exactly 9 pairs, 78%],
  [`ddMalSylv1`], [diploid], [no pairs, 0% duplicated; not assessable],
  [`fLepOcu1`, `fAmiCal2`], [no duplication since 2R], [one pair each, 0% in blocks; not assessable],
  [`SchCurv1`, `SchYoun1`], [autotetraploid snow carps, fusions (Xie et al. 2026)], [all published fusions, in wave order; mixed with fusions],
  [`drLytSali1`], [tetrasomic by classical genetics], [13/15 `tetrasomic_like`],
  [`icStrMela3`], [triploid, AAB (k-mer ploidy plot)], [3 copies, AAB-like (high)],
)

// ================================================================ 8
= Interpreting results, and the limits

== A reading order for any species

+ *Is the assembly what you think?* `sequences.tsv`, `unplaced.tsv`,
  `chrom_label_corrections.tsv`, `ploidy_summary.tsv`. Many corrections,
  `ambiguous` rows, or uneven copy counts come first.
+ *Is there a duplication?* Homeolog pairs, the partition, the homeology map.
  If none: two-copy genomes are not assessable for origin.
+ *How do the copies relate?* States per chromosome; for two-copy genomes the
  pooled pairs and residual tetrasomy.
+ *What traces of transition?* Fusions, regional splits, residual tetrasomy,
  exchange, chromosomes with two partners, divergence spread.
+ *Then* the summary sentences, with their evidence.

== Patterns and what they mean

#table(
  columns: (1.6fr, 2.5fr),
  table.header([*Pattern*], [*Most likely reading*]),
  [all copies equally close; no lineages], [tetrasomic: recent autopolyploid, not yet rediploidized],
  [two lineages at one depth everywhere, synchronous pairs], [allo-like: allopolyploid, or a long-rediploidized autopolyploid],
  [mixture of tetrasomic and split chromosomes, fusions], [rediploidization in progress (snow carps)],
  [diverged pairs with near-allelic stretches on ≥ 2 pairs], [autopolyploid origin with residual tetrasomy (salmonid-like); check terminal fraction and assembly history],
  [chromosomes with blocks on two partners], [karyotype rearranged since the duplication],
  [one copy apart, same haplotype always], [assembly/phasing problem, or a split-out extra set (AAB)],
  [one copy apart, haplotype changing], [AAAB-like, a divergent extra genome],
  [exchange in one copy], [homeologous exchange, or a contig on the wrong homeolog],
)

== Limits

- One individual: structure, not inheritance; not origin.
- The assembly is part of the measurement (Section 2.3).
- Old duplications fade below k-mer detection (~8–10% divergence).
- Thresholds are provisional (tuned on simulations, snow carps, Glechoma,
  _Malva_) and every one is a command-line flag (Appendix A).
- Holocentric genomes (_Eleocharis_) break and fuse chromosomes readily, so
  chromosome numbers need not track ancestral units.

== Open questions

- Confirm the Glechoma chr14/chr17 block, the Lathraea chr03/chr04 and the
  Pilosella chr04/chr05 residual segments with dotplots.
- Check weevil HAP1 chr09 (between A and B).
- Gene-based methods for old duplications (salmonids, snails).
- The `SchCurv1` chr17 block against Xie et al.'s coordinates.

// ================================================================ 9
= Every species so far

Each card below is generated from the result files
(`docs/guide/make_handbook_species.py`). The table gives the overview; the
cards give, per species: what is known, the assembly and any numbering fixes,
distances, duplication by each method, states, fusions, residual tetrasomy and
exchange with their support, the four answers with confidence and evidence,
and caveats. "Dup. (map)" is the share of the genome in reciprocal homeolog
blocks; "–" means the map has not run (one haplotype) or the run predates it.

#include "handbook_species.typ"

// ================================================================ A
= Appendix A: parameters

#params(
  [`--k`], [15,23], [largest k above the chance floor is used; k = 23 on every panel pair],
  [`--min-len`], [1 Mb], [chromosome-scale sequences only],
  [rename ratio / swap], [3× / 2×], [real mismatches were 4.6–96×; `ddHesMatr1`'s swap 3.1× and 2.3×],
  [extra-set ratio, min units], [1.3×, 3], [weevil B set 1.5–2.4×; a fusion leftover is one unit],
  [`--fdr-alpha`, `--min-effect`], [0.05, 2%], [a 1.25% difference reached z = −4 on simulated diploids],
  [`--window`, `--window-k`], [250 kb, 15], [allelic resolution along chromosomes],
  [map k, z, min block], [23, 3, 6 windows + null], [chance hits negligible at k = 23],
  [map pair minimum], [2 Mb reciprocal], [one short block is not a pair],
  [`--dist-split`, `--te-split`, `--window-split`], [1.25, 2.0, 1.25], [tuned on snow carps],
  [`--min-segment-bp`], [2 Mb], [shortest segment or block],
  [`--residual-ratio`], [3], [near-allelic = within 3× the allelic level],
  [`--residual-homeolog-frac`], [0.25], [masu residual windows 0.1× homeolog; false ones 0.7–0.9×],
  [`--residual-min-run`], [8 windows], [near-homozygous genomes had nulls of 3],
  [allelic cap], [1.5× genome median], [a large exchange hid itself at 5×],
  [exchange ratio], [2×], [homeolog clearly closer than own homolog],
  [fusion: long ratio, z, sibling excess], [1.4×, 10, 2.5×], [SchCurv1 fusion 3.4×, homeology 0.9×],
  [partition z for pooling], [10], [daInuConz1 10.6],
  [summary: residual pairs, duplication pairs, map fraction], [2, 3, 30%], [one pair is not an origin; gar and bowfin pass one pair],
)

= Appendix B: glossary

/ Allelic distance: distance between copies of the same chromosome number.
/ Block: a run of windows whose closest other chromosome is consistently the same.
/ Containment: fraction of one sequence's k-mers present in another.
/ Disomic / tetrasomic: pairing with a fixed partner / at random among four.
/ Extra set: a second chromosome set held in one haplotype file under other numbers.
/ Homeolog: a chromosome derived from the same ancestral chromosome through polyploidy.
/ Homeologous exchange: a stretch of one homeolog replaced by the other's sequence.
/ Pooling: reading two homeologous numbers' copies together.
/ Reciprocal block: a block whose partner has a block back on it.
/ Residual tetrasomy: regions of an otherwise rediploidized genome where homeologs still behave as one pool.
/ Rediploidization: the return of a polyploid to disomic behaviour, region by region.

= Appendix C: references

#set par(justify: false, hanging-indent: 1.2em)
#set text(size: 9.3pt)

Deb, S. K., Edger, P. P., Pires, J. C. & McKain, M. R. (2023). Patterns, mechanisms, and consequences of homoeologous exchange in allopolyploid angiosperms. _New Phytologist_ 238: 2284–2304.

Gaynor, M. L., Feng, K., Soltis, D. E., Soltis, P. S. & Smith, S. A. (2026). All detectable ancient whole-genome duplications involve hybridization. _bioRxiv_ preprint, not peer reviewed.

Lien, S. et al. (2016). The Atlantic salmon genome provides insights into rediploidization. _Nature_ 533: 200–205.

Myers, G. FastK: a k-mer counter for high-fidelity shotgun datasets. github.com/thegenemyers/FASTK

Ondov, B. D. et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. _Genome Biology_ 17: 132.

Twyford, A. D. et al. (2025). The polyploid continuum and the landscape of polyploid genomic variation. _American Journal of Botany_ 112: e70121.

Xie, C. et al. (2026). Chromosomal fusions trigger rediploidization of autopolyploid genomes. _Nature_ 654: 706–713.
