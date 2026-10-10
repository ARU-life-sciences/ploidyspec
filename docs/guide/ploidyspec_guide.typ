// ploidyspec: a guide to the tool, its outputs and what we have found so far
// Build: typst compile ploidyspec_guide.typ

#let accent = rgb("#2b6cb0")
#let warm = rgb("#c05621")
#let soft = rgb("#ebf4ff")
#let softwarm = rgb("#fffaf0")
#let grey = rgb("#4a5568")

#set document(title: "ploidyspec: a guide", author: "ploidyspec project")
#set page(
  paper: "a4",
  margin: (x: 2.3cm, top: 2.4cm, bottom: 2.4cm),
  header: context {
    if counter(page).get().first() > 1 [
      #set text(8pt, fill: grey)
      ploidyspec: a guide #h(1fr) October 2026
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

#let note(title: none, body) = block(
  fill: soft, inset: 10pt, radius: 4pt, width: 100%, stroke: (left: 3pt + accent),
  [#if title != none [#text(weight: "bold", fill: accent, title) \ ]#body],
)
#let caution(title: none, body) = block(
  fill: softwarm, inset: 10pt, radius: 4pt, width: 100%, stroke: (left: 3pt + warm),
  [#if title != none [#text(weight: "bold", fill: warm, title) \ ]#body],
)
#let stage(name, body) = box(
  fill: soft, stroke: 0.6pt + accent, radius: 4pt, inset: 6pt, width: 100%,
  [#text(weight: "bold", size: 9pt, fill: accent, name) \ #text(size: 8pt, body)],
)
#let arrow = align(center + horizon, text(14pt, fill: accent, sym.arrow.r))

// ---------------------------------------------------------------- title page
  #v(3.5cm)
  #text(30pt, weight: "bold", fill: accent)[ploidyspec]
  #v(0.2cm)
  #text(15pt)[Reading polyploidy, its origin and its return to diploidy \ from haplotype-resolved genome assemblies]
  #v(1cm)
  #line(length: 100%, stroke: 1pt + accent)
  #v(0.5cm)
  #text(11pt)[A guide to the method, the tool, its outputs, \ and what the first 45 assemblies show]
  #v(3cm)
  #grid(columns: (1fr, 1fr), gutter: 1cm,
    note(title: "What this document is")[
      A single place to understand the whole project: why it exists, how the
      software works, what every output means, how to use the outputs to answer
      concrete questions about a genome, and what we have found so far,
      including what we got wrong and corrected.
    ],
    note(title: "Status")[
      *Superseded in part (10 October 2026):* the companion
      `ploidyspec_handbook.pdf` gives every algorithm step by step, all outputs
      with walkthroughs, the new residual-tetrasomy, exchange, homeology-map
      and odd-ploidy readings, and current summaries of all 59 species. Its
      chapters 3–9 replace chapters 5–7 and Appendix A here.
      Written 7-8 October 2026. The software is on `main`. The panel was
      re-run after the fixes described in Chapter 8, so the results chapter
      reflects the current method.
    ],
  )

#pagebreak()

// ---------------------------------------------------------------- contents
#outline(title: [Contents], indent: auto, depth: 2)

// ================================================================ 1
= At a glance

Give ploidyspec the haplotype-resolved assembly of one individual (two or
more haplotype FASTA files, or one file with haplotypes named in the
headers) and it answers four questions, each with the evidence behind it:

#table(
  columns: (1.15fr, 2.2fr, 1.4fr),
  align: (left, left, left),
  table.header([*Question*], [*How ploidyspec answers it*], [*Where to look*]),
  [*1. What is the ploidy?*],
  [Counts the haplotype copies of each chromosome, then looks for chromosome
   numbers that are themselves duplicated (homeolog pairs). Copies × duplicated
   sets gives the ploidy relative to the base number.],
  [`matrix/ploidy_summary.tsv`, `homeologs/homeolog_pairs.tsv`, `panel_summary.tsv`],

  [*2. Does it behave like an auto- or an allopolyploid, and how sure are we?*],
  [Measures whether the copies of each chromosome are interchangeable
   (tetrasomic-like) or fall into two lineages that do not mix
   (disomic-like), and whether that pattern is consistent across
   chromosomes. Several independent lines of evidence give the confidence.],
  [`rediploidization/`, `structure/`, Section #link(<sec-q2>)[6.2]],

  [*3. What do the TE markers say?*],
  [Counts repeat (transposable element) k-mers that are abundant in one
   copy but not another: a fossil record of repeat activity in separate
   lineages. Read as a contrast within a genome, not as an absolute score.],
  [`subgenomes/`, Section #link(<sec-te>)[6.3]],

  [*4. Is there evidence of rediploidization?*],
  [Looks for chromosome fusions between copies, chromosomes whose copies have
   split into lineages along all or part of their length, and ancient pairs
   that resolved at different times.],
  [`rediploidization/`, Section #link(<sec-q4>)[6.4]],
)

For each species the `summary` stage writes these four answers, with a
confidence and the evidence, to `summary.md`. Appendix #link(<app-species>)[A]
has them for every species in the panel.

#v(0.4em)
It does all of this from DNA sequence alone: no gene annotation, no
reference genome, no read data, no alignment. The core measurement is the
k-mer distance between chromosome copies, whole and in windows.

#caution(title: "The one thing to remember")[
  Sequence from a single individual shows how its chromosome copies relate
  to each other today. It shows inheritance-like structure, not inheritance
  itself (that needs crosses or populations), and not origin. A long
  rediploidized autopolyploid can look exactly like an allopolyploid. The
  outputs are worded accordingly: "tetrasomic-like", "resolved lineages",
  never "is an autopolyploid".
]

// ================================================================ 2
= Background: why this question, why now

== Polyploidy is common, and its categories are blurry

Whole-genome duplication (WGD) has happened again and again in the history
of plants and animals, and nearly all flowering plants descend from
polyploid ancestors. Polyploids are traditionally sorted into two kinds.
*Autopolyploids* double a genome within one species; their chromosome
copies are near-identical and pair freely at meiosis (*polysomic*, or
*tetrasomic* in a tetraploid, inheritance). *Allopolyploids* combine the
genomes of two diverged species; each chromosome has a preferred partner
from its own parent (*disomic* inheritance), and the two parental sets
remain as *subgenomes*.

Twyford et al. (2025) argue that this binary is better seen as a
*continuum*, or a multidimensional "poly-space", because the definitions
(by origin, by meiotic pairing, by inheritance) often disagree, and because
genomes move through that space over time. Their Box 2 makes the point
that matters most here: after duplication a genome undergoes
*diploidization* (or *rediploidization*), and a diploidized autopolyploid
can have diverged duplicate chromosomes "indistinguishable from those in a
classic allopolyploid". They suggest measurable axes to place genomes in
poly-space: heterozygosity, sequence divergence between duplicated sets,
gene-level divergence, and structural variation, ideally from
haplotype-resolved genomes. ploidyspec measures several of these directly
(Section #link(<sec-polyspace>)[8.3]).

== Rediploidization: the time axis of the continuum

Xie et al. (2026) caught rediploidization near its start. Snow carps
(Schizothoracinae) share one autotetraploid origin. In _Schizopygopsis
younghusbandi_ most of the 25 chromosome quartets are still tetrasomic:
their four copies are about equally similar (Ks ≈ 0.01). But on ten
chromosomes, two of the four copies have fused with another chromosome
(five fusions). On those chromosomes, fused copies pair with fused copies
and unfused with unfused, so the quartet has split into two lineages that
diverge. Divergence is highest at the fusion site and spreads outwards
along the arms, and the fusions happened in three "waves" of different
ages. A separate mechanism, a centromeric inversion, made the short arm of
chromosome 17 disomic while the long arm stayed tetrasomic. The authors
also note that downstream consequences of rediploidization, such as biased
gene loss and expression between the two lineages, mimic subgenome
dominance in allopolyploids.

This is the pattern ploidyspec was built to read: per chromosome, and
along each chromosome, whether the copies still behave as one pool or have
separated. The snow carp genomes are part of our panel, and we recover
their fusions and their wave order without being told where to look
(Section #link(<sec-snowcarp>)[7.2]).

== Homoeologous exchange blurs things from the other side

Allopolyploids do not stay neatly separated either. Deb et al. (2023)
review homoeologous exchange (HE): recombination between chromosomes of
different subgenomes. HE can be reciprocal, or "with replacement", where one
subgenome's segment is overwritten by a copy of the other's, homogenizing
that region. HE is common in young allopolyploids, concentrated towards
chromosome ends, often biased towards one subgenome, and leaves a mosaic.
For ploidyspec this means an allopolyploid can show local blocks where
homeologs are near-identical. Those blocks look like the regional splits
of a rediploidizing autopolyploid seen in a mirror.

== Divergence records isolation, not duplication

A preprint by Gaynor et al. (2026, bioRxiv; not yet peer reviewed) makes
an argument that applies to any divergence-based method, including ours.
Using haplotype-phased polyploid assemblies, they find that synonymous
divergence (Ks) builds up only between gene copies that do not recombine.
Copies on chromosomes that pair at meiosis stay similar, whatever the
ploidy. Duplicate-gene divergence therefore detects meiotic isolation, or
divergence inherited from hybridizing parents, rather than genome
doubling as such. That is why tetrasomic autopolyploid copies stay close.

ploidyspec's distances behave the same way, and we use this deliberately.
Low distance among all copies of a chromosome means the copies still mix.
Two clusters at a consistent depth mean two lineages have been isolated
from each other. Divergence alone cannot say whether the isolation came
from hybridization (allo) or evolved later (rediploidized auto). It also
explains a blind spot: an old autopolyploidy that stayed tetrasomic leaves
little divergence to find.

== Other ways to tell auto from allo, and how ploidyspec differs

- *Gene-loss balance (P-index; Wang et al., 2019).* After duplication, an
  allopolyploid tends to lose genes unequally between subgenomes (biased
  fractionation), and an autopolyploid equally. Wang et al. turn this into
  a 0-1 index from gene collinearity against an outgroup genome, with
  about 0.3 separating auto from allo. It needs gene annotation, synteny,
  and a suitable outgroup. Xie et al. found biased gene loss towards the
  fused copies in an autopolyploid, so biased fractionation can also arise
  after rediploidization.
- *Duplicated BUSCOs and gene clusters.* McHale et al. (2025) detect the
  ancient land snail and slug WGD from the proportion of duplicated BUSCO
  genes (16-21% in Stylommatophora versus 0.5-1.5% in other molluscs). They
  follow it in the patchwork retention of duplicated Hox clusters,
  mirroring vertebrates. These approaches suit old duplications whose copies
  have diverged beyond the reach of k-mer comparison.
- *k-mer spectra from reads* (for example GenomeScope 2.0 and Smudgeplot,
  Ranallo-Benavidez et al., 2020; Tetmer, Becher et al., 2020) estimate
  ploidy and subgenome divergence without an assembly, but summarize the
  whole genome in one number.
- *Allele dosage and segregation.* Gaynor et al. and Xie et al. use
  genotype classes (AAAB/ABBB versus AABB) or resequencing of many
  individuals to infer disomic versus polysomic inheritance directly.

ploidyspec sits between these. It uses chromosome-scale, haplotype-resolved
assemblies of one individual, as produced at scale by the Darwin Tree of
Life project, and reads structure chromosome by chromosome and along
chromosomes, without annotation. It complements gene-based methods rather
than replacing them.

// ================================================================ 3
= What we are trying to achieve

== The aim

Given a haplotype-resolved assembly, describe where the genome sits on two
axes:

+ *Origin-like structure:* do the chromosome copies form one interchangeable
  pool (auto-like), or separate lineages at a consistent depth across the
  genome (allo-like), or something in between?
+ *Rediploidization:* has a once-tetrasomic genome started to split into
  lineages, where, by what mechanism (fusion, inversion, gradual
  divergence), and how far has it gone?

The tool should run on any species without tuning, report evidence rather
than labels, and say when a question cannot be answered from the data.

== Design principles

- *Reference-free and annotation-free.* Everything comes from k-mers of the
  assembled sequence, so a newly assembled species can be analysed the day
  its assembly is released.
- *Per chromosome, then per window.* Genome-wide summaries hide exactly the
  heterogeneity that rediploidization creates.
- *Descriptive states, not verdicts.* Each chromosome gets a state with the
  numbers that produced it.
- *Several independent lines of evidence.* Whole-chromosome distance,
  windowed distance, repeat markers, chromosome fusions and genome-wide
  partitions each fail differently, so agreement between them carries
  weight.
- *Validated on genomes with a known answer.* A simulator builds synthetic
  assemblies (diploid, auto, allo, rediploidizing, mislabelled) and an
  end-to-end test checks the whole pipeline against the truth.

// ================================================================ 4
= The tool

== Input

ploidyspec needs, for one individual, every haplotype at chromosome scale,
described by a *manifest*: a two-column table of FASTA path and haplotype
label.

```
# fasta_path                         hap_label
data/sp/sp.hap1.1.primary.fa.gz      HAP1
data/sp/sp.hap2.1.primary.fa.gz      HAP2
```

Use `AUTO` instead of a label when one file holds several haplotypes named
in the headers. Each chromosome-scale sequence (default ≥ 1 Mb) needs a
chromosome number, read from the header: `chromosome: 4` in the
description, or `SUPER_4` in the name. Other naming schemes are handled with
`--chrom-regex`/`--hap-regex`. Shorter or unnumbered sequences are listed in
`unplaced.tsv`, and the large ones are still tested as possible fused
chromosomes.

#note(title: "What counts as 'a copy'")[
  ploidyspec counts what the assembly contains. A tetraploid can be
  assembled as four haplotypes (four copies of each chromosome number) or as
  two haplotypes in which each base chromosome appears under two numbers
  (two copies each, arranged in homeolog pairs). Both are handled, and the
  answers to the four questions take the layout into account.
]

== Installing and running

```
conda env create -f environment.yml && conda activate ploidyspec
pip install .
ploidyspec all --manifest manifests/sp.tsv --outdir results/sp --threads 8
ploidyspec all ... --with-te-markers --with-windowed-homeologs   # optional stages
ploidyspec panel --results results --outdir meta                 # many species
```

The only external tools are `samtools` and FastK (`FastK`, `Logex`,
`Histex`, `Tabex`, `Profex`). Every stage can be run on its own, and each
skips work already done, so an interrupted run resumes. `report.html` in
the output directory collects tables and plots from whatever has run.

*Resources.* Memory scales with the longest chromosome, not with genome
size: about 40 MB per Mb for the matrix stage, and about 14 bytes per base
of the longest chromosome per thread for the windowed stage. Runtime on our
panel was 0.2-13 hours per species. k-mer tables dominate disk use;
`--cleanup` removes them and keeps every result.

== The measurement underneath everything

A *k-mer* is a word of k bases. ploidyspec builds a FastK table of every
k-mer in each chromosome copy and compares tables. For two sequences A and
B with Jaccard index $J = |A inter B| \/ |A union B|$ of their k-mer sets,
the Mash distance (Ondov et al., 2016) estimates the per-base divergence:

$ d = -1/k ln((2J)/(1+J)) $

Short k-mers collide by chance and saturate, long ones are too sensitive.
So the matrix stage computes $d$ at several k (default 15 and 23) and keeps
the largest k that clears the chance-collision floor. On all 53,309
chromosome pairs of our panel that was k = 23. For windows, ploidyspec
uses *containment* $c$, the fraction of one sequence's k-mers present in
the other, and the matching distance $d = -ln(c) \/ k$.

#note(title: "Scale to keep in mind")[
  Typical values in our panel: two alleles of the same chromosome, 0.0001
  to 0.04; ancient homeologs, 0.05 to 0.10; unrelated chromosomes, 0.07 to
  0.16. The confirmed diploid _Malva sylvestris_ (`ddMalSylv1`) has alleles
  0.0038 apart and unrelated chromosomes 0.116 apart, 31 times further.
]

// ================================================================ 5
= The modules, one by one

#figure(
  grid(
    columns: (1fr, auto, 1fr, auto, 1fr, auto, 1fr),
    gutter: 4pt,
    stage("1 prepare")[manifest → one FASTA per chromosome copy],
    arrow,
    stage("2 kmers")[FastK table per copy, at each k],
    arrow,
    stage("3 matrix")[all-vs-all copy distances; relabelling check],
    arrow,
    stage("4 homeologs")[duplicated chromosome numbers (FDR)],
  ) + v(6pt) + grid(
    columns: (1fr, auto, 1fr, auto, 1fr, auto, 1fr),
    gutter: 4pt,
    stage("5 windowed")[distance along each copy, position-free],
    arrow,
    stage("6 structure")[genome partitions, pair synchrony],
    arrow,
    stage("7 rediploidization")[fusions, lineage states per chromosome],
    arrow,
    stage("8 summary, 9 report")[four answers with confidence; one HTML page],
  ) + v(6pt) + grid(
    columns: (1fr, 1fr, 1fr),
    gutter: 6pt,
    stage("optional: te-markers")[repeat k-mers private to one copy],
    stage("optional: windowed-homeologs")[windows between homeolog pairs],
    stage("many species: panel")[cross-species tables, exploratory PCA],
  ),
  caption: [The pipeline. `ploidyspec all` runs stages 1-9 in order; each stage
  reads the outputs of the stages before it.],
)

== prepare

Reads the manifest, finds each chromosome-scale sequence, gives it a
haplotype and a chromosome number, and writes `sequences.tsv` (one row per
*unit*, a single chromosome copy such as `HAP2_chr05`) and `unplaced.tsv`
(everything else, with the reason). If many sequences land in
`unplaced.tsv`, the naming regex needs adjusting.

== kmers

Extracts each unit to its own FASTA and builds a FastK table for each k.
These are cached and reused by name. Each cached unit also records where it
came from (assembly file, its size and modification time, sequence name
and length), and is rebuilt if any of these change. This was added after a
re-run on a new assembly silently reused the old sequences
(Section #link(<sec-corrections>)[8.1]).

== matrix: copy distances and the relabelling check

Computes the distance between every pair of units, writes the full
matrix, heatmaps, `ploidy_summary.tsv` (copies per chromosome number) and
`homologous_chromosomes.tsv`.

It also checks *labelling*. Different haplotype assemblies are sometimes
numbered differently: one file's "chr2" can be another file's "chr4". The
check takes the file with the most units as reference. Each other file's
units are matched one-to-one to the reference chromosome numbers, closest
first. A unit is relabelled when its matched number is at least 3× closer
than its declared one. A pair of units that each match the other's number
(a swap) needs 2× on both sides. Ambiguous cases are listed in
`matrix/chrom_label_corrections.tsv`, and relabelled units are renamed on
disk so every later stage sees the corrected numbers.

== homeologs: are some chromosome numbers duplicates of each other?

After a WGD, some chromosome numbers are each other's ancient copies
(*homeologs*). For every pair of chromosome numbers this stage takes the
mean distance over all copy combinations. It compares that mean to the
background of all pairs using an empirical p-value, and controls the false
discovery rate at 5%. A pair is accepted when it is significant and at
least 2% closer than the background median. Outputs:
`homeolog_pairs.tsv` (accepted pairs), `homeolog_candidates_ranked.tsv`
(every pair, ranked) and `ploidy_ancestry_summary.tsv`. The last gives
each chromosome's partner and `distance_ratio`, the homeolog distance
divided by the allele distance.

== windowed: what happens along each chromosome

Each copy is cut into 250 kb windows along its own coordinates. Each window
is looked up, with FastK's profile mode, in each other copy's whole
chromosome. That gives the fraction of the window's k-mer positions found
anywhere in the other copy, and the distance $-ln(c) \/ k$ (k = 15). No
alignment or shared coordinates are needed, so indels, inversions and gaps
between assemblies do not matter. The result is one track per ordered pair
of copies, along the first copy.

#figure(
  image("fig/windowed_fix.png", width: 100%),
  caption: [Why the windowed stage was rewritten. Left: the original method
  compared windows at the same coordinates in each copy. Once two assemblies
  differ by an indel larger than a window, every later window compares
  different sequence. Two copies of `ddEmpNigr1` chr03 that are 0.004 apart
  over the whole chromosome read as almost completely different (Jaccard
  ≈ 0.99) in almost every window. Right: the position-free method on the same
  chromosome. HAP3, the close copy, stays near zero all along, and HAP2 and
  HAP4 form the other lineage. The spike near 29 Mb is sequence present in
  HAP1 and absent from HAP2.],
)

== structure: genome-wide partitions and pair synchrony

`genome_partition.tsv` asks whether the chromosome numbers fall into groups
that are closer within than between, more than chance, even when no single
homeolog pair is significant. Chromosome numbers are clustered (average
linkage) and, at every number of groups k, the between/within separation is
compared with 500 random groupings of the same sizes. The randomness is
seeded per species, so results reproduce. A strong, balanced two-group
split (`daInuConz1`: 8 vs 8, z = 10.6) is a candidate subgenome structure
that pairwise tests miss.

`pair_synchrony.tsv` gives the spread of divergence depth across the accepted
homeolog pairs (`pair_depth_cv`). One duplication event diverges every pair
to about the same depth (low CV, `daGleHede1` 0.04). Several events, or
pairs that stopped recombining at different times, spread them out.

== te-markers (optional): fossil repeat activity

See Section #link(<sec-te>)[6.3] for the idea. For each pair of copies of a
chromosome, at k = 13, the stage finds *high-copy* k-mers (count ≥ 100,
mostly transposable elements and other repeats). It then finds *markers*:
high-copy k-mers at least 2× more abundant in one copy than the other. The
summary is

$ "te_marker_fraction" = (m_A + m_B) / (h_A + h_B) $

with $m$ the markers and $h$ the high-copy k-mers of copies A and B.
`te_marker_fraction_by_lineage.tsv` splits chromosomes with three or more
copies into two groups by distance and reports `split_ratio`: the
cross-group fraction over the within-group fraction. Copies with too little
or too much repeat content for their length are flagged as assembly
problems. The method follows K. Jaron's fossil-TE approach for separating
subgenomes.

== rediploidization: fusions and per-chromosome states

The stage that combines everything. It has three parts.

*Fusion detection.* It tests two kinds of sequence:
- placed copies at least 1.4× the median length of their sibling copies;
- unplaced scaffolds at least half the median chromosome length.

For each, it measures how much of every other chromosome's k-mer content
it contains. A chromosome counts as a fused component when its containment
is far above background (robust z ≥ 10, ≥ 2× the background median, and
≥ a quarter of the top hit). For placed copies, the partner must also be
≥ 2.5× more enriched, per k-mer, than in the copy's unfused siblings
(`sibling_excess`). This rules out ancient homeology, which every copy
shares, and siblings that are merely fragments. The scaffold is called a
`fusion` only if the partner chromosome is missing from that haplotype.
`lineage_divergence` $= -ln(c) \/ k$ dates the fused lineage against the
unfused copies.

*Lineage structure per chromosome.* For each chromosome number with three or
more copies it reads:

- `dist_split`: split the copies into the two most separated groups;
  mean cross-group over mean within-group distance;
- `te_split`: the TE-marker `split_ratio` for the same chromosome;
- `split_extent`: along each copy, which windows show the split. Reported
  as `whole`, `regional` (with the segments) or `none`, from the median
  copy. It counts only if it beats a null: the same scan for groupings of
  the copies that cut across the observed one, pooled over the genome
  (Section #link(<sec-q4>)[6.4]).

*Pooling.* A two-haplotype assembly has only two copies of each chromosome
number, too few to split. ploidyspec pools a chromosome with its accepted
homeolog partner (four copies). Failing that, it uses its reciprocal best
match across a significant two-group genome partition, and reads the
pooled group.

Each chromosome then gets one `copy_state`:

#table(
  columns: (1.15fr, 2.6fr),
  table.header([*State*], [*Meaning*]),
  [`fusion_lineages`], [Fused in some haplotypes and not others: the fused and unfused copies form separate lineages (the snow carp mechanism).],
  [`resolved_lineages`], [Balanced split, whole-chromosome split ≥ 2.5×, holding along most of the chromosome; or, for a pooled group, the two numbers are separate lineages.],
  [`partially_resolved`], [Two or more of: distance split ≥ 1.25, TE split ≥ 2, a windowed split segment.],
  [`candidate`], [Exactly one of those.],
  [`tetrasomic_like`], [None: the copies are interchangeable. For a pooled group, the homeolog is as close as the allele.],
  [`one_divergent_copy`], [One copy stands apart from the rest along the chromosome (distance or TE split; a windowed segment alone does not count). This is an odd haplotype (assembly quality, or a divergent extra genome, AAAB-like), not two lineages. Every split of a triploid looks like this.],
  [`not_assessable`], [Fewer than three copies and nothing to pool with.],
)

== summary, report, panel, simulate, cleanup

`summary` answers the four questions of Chapter 6 for the species, each with
a confidence (high, medium, low or not assessable) and the evidence behind
it, following the decision tables of that chapter. It writes `summary.tsv`
and `summary.md`, and the report opens with it. Appendix A lists the
summaries for every species in the panel.

`report` writes one self-contained HTML page per species, with core
sections first and TE markers under "Supplementary". `panel` stacks every
species into `panel_summary.tsv`, `genome_partition.tsv` and
`rediploidization_panel.tsv`. It also writes a TE table and an exploratory
PCA over the core numbers in `supplementary/`. `simulate` writes synthetic
assemblies with a known answer, and `cleanup` deletes the large
intermediate files.

// ================================================================ 6
= Answering the four questions

This chapter is the practical guide: given a finished run, how to read
the outputs to answer each question, and how much to trust the answer.

== Question 1: what is the ploidy?

#figure(
  table(
    columns: (1.3fr, 1fr, 1fr, 2fr),
    table.header([*Assembly layout*], [*Copies per number*], [*Homeolog pairs*], [*Reading*]),
    [Two haplotypes, no pairs], [2], [none], [Diploid-like. (Or a polyploid whose duplicate sets have diverged beyond detection; check `genome_partition.tsv`.)],
    [Two haplotypes, numbers in pairs], [2], [most numbers], [Tetraploid assembled as two haplotypes: 2 copies × 2 sets. `daGleHede1`: 18 numbers in 9 pairs.],
    [Four haplotypes, no pairs], [4], [none], [Tetraploid, all four copies in view.],
    [Four haplotypes, numbers in pairs], [4], [most numbers], [Octoploid on half the numbers: `ddLepDrab1`, 4 × 16 = 64 chromosomes, 8 pairs, 2n = 8x on x = 8.],
    [Three haplotypes], [3], [-], [Triploid (`ddSalTria1`).],
  ),
  caption: [Reading ploidy from copies and homeolog pairs.],
)

Base number x is the number of chromosome numbers divided by the sets per
haplotype, and ploidy is copies × sets. *Confidence is high for the copy count*, which is set by
the assembly. That is also its limit: a collapsed or missing haplotype
lowers the count. `ddAraThal4` has only one usable haplotype, so its copy
count says nothing about ploidy. *Confidence in the number of sets* depends
on how many chromosome numbers fall into accepted pairs, and on whether
`ploidy_summary.tsv` is uniform. Check counts against the literature
(`meta/literature_ploidy.tsv`) when possible.

#caution(title: "Mislabelled chromosomes make ploidy look strange")[
  If one haplotype file numbers chromosomes differently, every comparison
  of "chr2" mixes two chromosomes. Always check
  `matrix/chrom_label_corrections.tsv`. Corrected rows are fixed
  automatically; an `ambiguous` row is a warning to look by hand.
]

== Question 2: auto-like or allo-like, and with what confidence? <sec-q2>

#figure(
  image("fig/three_shapes.png", width: 100%),
  caption: [The three shapes that matter, from real chromosomes (whole-chromosome
  distance between each pair of copies, × 1000). Left: all six distances
  similar, so the copies are one interchangeable pool. Middle: two tight pairs
  at a uniform distance from each other, so two lineages. Right: two lineages
  whose split coincides with a fusion carried by two of the four copies.],
)

Read the evidence in this order.

+ *Within each chromosome* (three or more copies): are the copies
  interchangeable or split into two lineages? This is `copy_state`, with
  `dist_split` behind it.
+ *Across chromosomes:* is the split at the same depth on every chromosome?
  An allo-like genome splits everywhere, at a consistent depth (one origin).
  A rediploidizing genome mixes tetrasomic-like and split chromosomes, at
  different depths.
+ *Between chromosome numbers* (two-haplotype assemblies): accepted homeolog
  pairs and their `pair_depth_cv`; strong balanced genome partitions.
+ *Repeat history:* a TE-marker split that agrees with the distance split
  (`te_split` ≥ 2 on the same chromosomes).

#figure(
  table(
    columns: (1.1fr, 2.2fr, 1.1fr),
    table.header([*Reading*], [*Evidence pattern*], [*Confidence*]),
    [*Allo-like (disomic-like) throughout*],
    [Most chromosomes `resolved_lineages`, at similar depth; or (two haplotypes) most numbers in homeolog pairs with low `pair_depth_cv` (≈ 0.03-0.05), resolved when pooled; TE split agrees.],
    [*High* when two or more lines agree across most chromosomes.],

    [*Auto-like (tetrasomic-like) throughout*],
    [Most chromosomes `tetrasomic_like`: all copies about equally close; no lineage split by distance, TE or windows.],
    [*High* with four or more copies. Not available from two haplotypes alone.],

    [*Mixed / transitional*],
    [Some chromosomes tetrasomic-like, others split or partially split; fusions; splits at different depths.],
    [*Medium*: the per-chromosome pattern is solid, the mechanism needs checking.],

    [*Cryptic subgenome structure*],
    [No accepted pairs, but a strong balanced genome partition (z ≥ 10).],
    [*Medium*: one line of evidence; confirm with karyology.],

    [*Not assessable*],
    [Two haplotypes, no pairs, no partition; or one usable haplotype.],
    [-],
  ),
  caption: [How to read origin-like structure, with confidence. "Lines of evidence"
  are whole-chromosome distance, windowed extent, TE markers, homeolog pairs and
  genome partitions.],
)

#caution(title: "Four things that look allo-like but are not")[
  - *A rediploidized autopolyploid.* It splits into lineages exactly like an
    allopolyploid: the snow carp fused chromosomes; Twyford et al. Box 2.
  - *One odd haplotype.* A copy assembled less completely than its siblings,
    or a mislabelled haplotype, isolates one copy (`one_divergent_copy`). If
    the same haplotype is odd on most chromosomes, suspect the assembly
    (`most_frequent_outlier_hap` in the summary).
  - *Pooled `resolved_lineages`.* Any two different chromosome numbers are
    far more distant than two alleles. So pooling a chromosome with its
    partner reads `resolved_lineages` even in a diploid (we checked on
    `ddMalSylv1`: splits of 22-47×). For pooled groups the evidence is the
    homeolog test or the partition z, and the informative outcome is the
    opposite one, `tetrasomic_like`.
  - *Absolute `te_marker_fraction`.* See Section #link(<sec-te>)[6.3].
]

== Question 3: TE markers, explained <sec-te>

*The idea.* Transposable elements (TEs) copy themselves around the genome in
bursts. Two lineages that evolve apart, such as the parents of an
allopolyploid, have separate bursts. Each lineage ends up with its own
recently amplified repeats, and those stay behind as a "fossil" signature.
In an allopolyploid, a subgenome carries its parent's repeats and the other
subgenome largely does not. So many abundant repeat k-mers are much more
common in one copy than in the other: high `te_marker_fraction`. In an
autopolyploid, all copies come from one genome with one repeat history, so
few should be.

*What we found.* The confirmed diploid `ddMalSylv1` has a per-chromosome
`te_marker_fraction` of 0.05-0.25 (median 0.14). Ordinary allelic
differences produce markers too: TE insertions present in one allele but
not the other. This range overlaps the allopolyploid anchors (`daGleHede1`
median 0.22, `drTriRepe1` 0.15). The autotetraploid snow carp `SchCurv1`
spans 0.07 to 0.61.

#figure(
  image("fig/te_fraction.png", width: 100%),
  caption: [Absolute `te_marker_fraction` does not separate diploid, auto- and
  allopolyploid. Each point is one chromosome copy pair. The shaded band is the
  range in the confirmed diploid. `daInuConz1` is low because it compares
  allelic copies within each subgenome; its subgenome split lies between
  chromosome numbers, which this index never compares.],
)

*How to use TE markers.* Use them as a *contrast within one genome*. For a
chromosome with three or more copies, split the copies by distance and ask
whether markers are concentrated between the two groups (`split_ratio` ≥ 2)
rather than within them. When the TE split and the distance split agree,
the lineages have separate repeat histories, an independent line of
evidence. `SchCurv1` chr19, the fused chromosome, has a TE split of 5.6
against 0.8-2.2 for most of its genome.

#note(title: "Reading te_split in practice")[
  - `te_split` ≈ 1: repeats do not track the lineages.
  - `te_split` ≥ 2 on the same chromosomes as `dist_split` ≥ 1.25: two
    lineages with separate repeat histories.
  - High absolute fractions with `te_split` ≈ 1: allelic TE polymorphism,
    or an assembly artefact. Check `flagged_units`.
]

== Question 4: is there evidence of rediploidization? <sec-q4>

Rediploidization leaves four kinds of trace that ploidyspec can see:

+ *Chromosome fusions between copies*: `fusions.tsv` rows with `status =
  fusion`. These are the strongest single evidence: a copy that has absorbed
  another chromosome cannot pair freely with its unfused siblings. Their
  `lineage_divergence` orders the fusions in time.
+ *Chromosomes whose copies have split into lineages* in a genome where
  others have not: a mixture of `resolved_lineages`/`fusion_lineages` and
  `tetrasomic_like`.
+ *Splits along part of a chromosome*: `partially_resolved` with a
  `regional` extent. This is the pattern of a centromeric inversion (snow
  carp chr17), or of divergence spreading out from a fusion site.
+ *Ancient pairs at different depths*: a high `pair_depth_cv`. Duplicate
  pairs that stopped recombining at different times, such as `drMyrSpic1`
  at 0.33 against 0.03-0.05 for single-event anchors.

#figure(
  table(
    columns: (1.2fr, 2.3fr),
    table.header([*Reading*], [*Pattern*]),
    [Rediploidization under way], [Fusions and/or a mix of split and tetrasomic-like chromosomes in a genome with three or more copies.],
    [Completed (or allo-like)], [Every chromosome split at a consistent depth; pairs synchronous.],
    [Not started], [All chromosomes tetrasomic-like, no fusions.],
    [Partial or ongoing on a few chromosomes], [Regional splits on a minority of chromosomes. Check that segments are not noise (see the caution below).],
  ),
  caption: [Reading rediploidization.],
)

#note(title: "How regional splits are kept honest")[
  Window distances are real measurements with real noise, so a run of
  split windows can occur by chance. Two rules guard against it.

  - *A null for each species.* The window scan is repeated for groupings of
    the copies that cut across the observed one: for a balanced split, the
    other balanced splits; for a one-copy split, the splits isolating each
    other copy. Their coverage, pooled over the genome, gives the species'
    noise level. A split's extent counts only if it covers more of the
    chromosome than that null's 95th percentile and more than every crossing
    grouping of the same chromosome (`window_covered`, `window_null_q95`,
    `split_extent_raw` show the working). In `daBudDavi1` this removed ten
    calls resting on single 2-3 Mb segments. Snow carp chr19 (0.89 against a
    null of 0.07) and chr17 (0.24) pass, as does the simulated regional split.
  - *One copy apart needs chromosome-wide evidence.* In a polysomic genome,
    one of four copies often carries a divergent haplotype block somewhere.
    That is ordinary polymorphism, not a lineage. So `one_divergent_copy`
    requires a distance or TE split across the chromosome. On regional
    segments alone, `drLytSali1`, tetrasomic by classical genetics, read
    12/15 chromosomes as `one_divergent_copy`; it now reads 13/15
    `tetrasomic_like`.
]

// ================================================================ 7
= What we have found so far

== Validation on simulated genomes

`ploidyspec simulate` builds synthetic assemblies with a known answer:
random sequence with TE families inserted at high copy number. Every
haplotype except HAP1 carries a transposition and an inversion, so the
copies are not collinear, as with real independently assembled haplotypes.
The end-to-end test runs the whole pipeline on each and checks the result.

#figure(
  table(
    columns: (1.4fr, 2.6fr, 1fr),
    table.header([*Scenario*], [*What it contains*], [*Chromosomes correct*]),
    [diploid], [two haplotypes of one genome], [7 / 7],
    [autotetraploid], [four equivalent haplotypes], [7 / 7],
    [autotetraploid, 2 haplotypes], [the same, assembled as two haplotypes with each chromosome numbered twice], [14 / 14],
    [allotetraploid], [two diverged progenitors, each with its own TE burst], [14 / 14],
    [rediploidized], [autotetraploid with a whole-chromosome lineage split, a split over 40% of one chromosome, one placed and one unplaced fusion], [7 / 7, both fusions found],
    [mislabelled], [autotetraploid with one haplotype numbered shifted by one, and a swap in another], [7 / 7, all 9 relabelled],
  ),
  caption: [End-to-end validation (11 tests, all passing). On the rediploidized
  scenario, the pre-fix windowed code calls chr07 `one_divergent_copy` and
  loses chr02's regional split, so the test now guards against that bug.],
)

== The snow carps: blind recovery of a published result <sec-snowcarp>

- *Fusions.* In `SchYoun1` all five published fusions are recovered from
  unplaced scaffolds. In `SchCurv1`, chr19+chr22 is fused in HAP3 and HAP4
  (containment 0.44-0.45, z ≈ 48, 3.4× enriched over unfused siblings).
- *Their order in time.* `lineage_divergence` ranks the `SchYoun1` fusions
  as chr19+22 (≈ 0.053), then chr04+15, chr08+16 and chr20+23
  (≈ 0.026-0.028), then chr11+14 (≈ 0.015). This is the three-wave order
  Xie et al. inferred, found without knowing it.
- *Lineage split.* The fused chromosome splits along its whole length,
  in 92% of windows.
- *chr17.* The analysis flagged chr17 unprompted. With the corrected windows,
  every copy places the split in the last ≈ 8 Mb of the chromosome. Xie et al.
  report a short-arm-only disomic shift caused by a centromeric inversion;
  whether that 8 Mb block is the short arm still has to be checked against
  their coordinates.

#figure(
  image("fig/snowcarp_windows.png", width: 100%),
  caption: [Window distances along HAP1 of `SchCurv1`. chr19: the fused copies
  (HAP3, HAP4) are further from HAP1 than HAP2 is, along the whole chromosome.
  chr17: the copies separate only towards the end.],
)

== The panel: 45 assemblies <sec-panel>

#figure(
  image("fig/panel_states.png", width: 100%),
  caption: [Rediploidization state of each chromosome number, per species (label:
  copies × chromosome numbers), after the corrections in Section
  #link(<sec-corrections>)[8.1]. `daInuConz1`, `dmRanRepe1` and `dcCerAlpi1`
  read `resolved_lineages` through partition pooling: the evidence there is
  the partition z, not the bar.],
)

*Two-haplotype allopolyploid anchors* (`daGleHede1`, `drTriRepe1`,
`drSorDevo1`, `ddSalPent1`, `drMyrVert1`, `drRosSpin1`, the _Potamogeton_
species) have most chromosome numbers in homeolog pairs, synchronous pair
depths (`pair_depth_cv` 0.03-0.05 for the clearest), strong genome
partitions, and resolved pooled lineages. That is what a single, old,
disomic duplication looks like.

*Cryptic subgenome structure.* `daInuConz1` (_Pentanema conyzae_) has no
accepted homeolog pair but an 8 vs 8 partition at z = 10.6. That matches
its documented tetraploid count 2n = 4x = 32 on x = 8. `dmRanRepe1` (8 vs 8,
z = 10.9) and `dcCerAlpi1` (16 vs 20, z = 18.9) show the same kind of
signal.

*Autopolyploid anchors* read as they should. `drLytSali1`, tetrasomic by
classical genetics, is 13/15 `tetrasomic_like` on its current assembly. The
*snow carps* combine fusions and whole-chromosome lineage splits with a
tetrasomic-like majority (`SchCurv1` 17/25, `SchYoun1` 14/25 plus 10
fusion-lineage chromosomes): rediploidization in progress.

*Two lineages per chromosome.* `ddEmpNigr1` shows, on most chromosomes, two
lineages at a uniform cross distance (≈ 0.020). That suggests a single
origin with disomic-like structure; it is now 2 `resolved` and 8
`partially_resolved` of 13.

*One divergent copy.* `daPilAura1` shows one copy in four apart on all 9
base chromosomes, with the odd copy changing haplotype. That is
AAAB-like. `ddHypMacu1` has the same haplotype (HAP1) apart on all 8
chromosomes, matching a known assembly problem in that haplotype.

*An octoploid, read correctly at last.* `ddLepDrab1` has four haplotypes
of 16 chromosome numbers in 8 homeolog pairs at a uniform depth (CV 0.06).
Within each number the four copies are mostly tetrasomic-like. That fits two
x = 8 sets each present four times (AAAABBBB-like), matching 2n = 64 on
x = 8. Its earlier "strongest allopolyploid candidate" signals came from one
haplotype file numbering 14 chromosomes differently (Section
#link(<sec-corrections>)[8.1]).

#figure(
  table(
    columns: (1.25fr, 0.6fr, 0.95fr, 0.75fr, 0.8fr, 0.8fr, 2.1fr),
    table.header([*Species*], [*Copies × nos.*], [*Allele dist. (median)*], [*Pairs*], [*Pair depth CV*], [*Partition z (k)*], [*Reading*]),
    [`daGleHede1`], [2 × 18], [0.014], [9], [0.042], [23.8 (9)], [allo anchor; one synchronous event],
    [`drSorDevo1`], [2 × 34], [0.011], [17], [0.030], [42.5 (17)], [allo anchor],
    [`ddSalPent1`], [2 × 38], [0.009], [19], [0.038], [47.2 (19)], [allo anchor],
    [`drMyrSpic1`], [2 × 21], [0.009], [8], [0.332], [24.7 (7)], [hexaploid; asynchronous pairs],
    [`daInuConz1`], [2 × 16], [0.0007], [0], [-], [10.6 (2)], [cryptic 8 vs 8 split, x = 8],
    [`dmRanRepe1`], [2 × 16], [0.023], [0], [-], [10.9 (2)], [cryptic 8 vs 8 split],
    [`dcCerAlpi1`], [2 × 36], [0.008], [0], [-], [20.8 (3)], [cryptic 16 vs 20 split],
    [`ddMalSylv1`], [2 × 21], [0.0038], [0], [-], [none], [diploid anchor],
    [`drLytSali1`], [4 × 15], [0.025], [5], [0.047], [10.5 (5)], [auto anchor; 13/15 tetrasomic-like],
    [`SchCurv1`], [4 × 25], [0.013], [2], [-], [-], [auto; 1 fusion; 17/25 tetrasomic-like],
    [`SchYoun1`], [4 × 25], [0.017], [0], [-], [-], [auto; 5 fusions; 14/25 tetrasomic-like],
    [`ddLepDrab1`], [4 × 16], [0.021], [8], [0.062], [18.8 (8)], [octoploid, AAAABBBB-like on x = 8],
    [`ddEmpNigr1`], [4 × 13], [0.020], [1], [-], [-], [two lineages per chromosome],
    [`daPilAura1`], [2 × 18], [0.039], [9], [0.105], [22.7 (9)], [one divergent copy in four],
  ),
  caption: [Selected species from `panel_summary.tsv`. Pairs:
  accepted homeolog pairs. Partition: strongest genome partition. Allele
  distance: median distance between copies of the same chromosome number.],
)

// ================================================================ 8
= Caveats, corrections and open questions

== What we got wrong, and fixed <sec-corrections>

Checking results by hand during Phase 2 found four problems. All four are
now fixed and covered by tests. They are listed here because
earlier write-ups quoted numbers that depended on them.

#table(
  columns: (1.2fr, 2.4fr, 1.6fr),
  table.header([*Problem*], [*What happened*], [*Fix*]),
  [Window registration],
  [Windows were compared at equal coordinates. After the first indel larger
   than a window they compared different sequence, so across the panel only
   0-34% of windows were usable. Every windowed-derived number was affected:
   regional splits, `split_extent`, run-length and windowed-CV metrics.],
  [Position-free windows (FastK profiles). Panel re-run.],

  [Chromosome relabelling],
  [`ddLepDrab1` HAP4 numbers 14 of 16 chromosomes differently. The check
   rejected the whole batch because two units tied, so HAP4 stayed
   mislabelled. Its "strongest allo candidate" signals were all artefacts.
   `ddHesMatr1` has a HAP1 chr01/chr02 swap below the old threshold; its
   only homeolog pair was the swap.],
  [One-to-one matching; mutual swaps at 2×. Both re-run.],

  [Stale cache],
  [Per-copy FASTAs and k-mer tables were reused by name. A re-run of
   `drLytSali1` on a new, renumbered assembly kept the old sequences.],
  [Cache records its source; rebuilt on any change.],

  [Window noise],
  [After the windowed fix, chance runs of split windows and ordinary
   haplotype blocks were read as lineage splits, inflating
   `one_divergent_copy` and `candidate`.],
  [Null from crossing groupings; one copy apart needs chromosome-wide
   evidence (Section #link(<sec-q4>)[6.4]).],
)

Several inheritance-mode metrics from earlier work (`partition_consistency`,
run length, `flip_rate`, the windowed CVs and a combined score) were
retired, because labelling and window registration drove them. Their last
values are archived in `meta/archive/` with a note.

== Limits of the approach

- *One individual.* Inheritance and origin are not observed, only structure.
- *Assembly quality is part of the measurement.* Missing or extra sequence
  in one haplotype reads as divergence: a copy missing 10% of a chromosome
  reads as the divergent copy, with a regional split. Compare lengths and
  `unplaced.tsv` before reading biology into one odd haplotype.
- *Old duplications fade.* Homeolog distances approach the background of
  unrelated chromosomes, and k-mer methods lose them. Ancient WGDs such as
  the land snail event (McHale et al., 2025) are better seen with gene-based
  approaches (duplicated BUSCOs, synteny, Ks).
- *Stay-tetrasomic autopolyploids leave little divergence* (Gaynor et al.,
  2026). They are visible only when the assembly shows all copies.
- *Homoeologous exchange* can make allopolyploid chromosomes look locally
  auto-like (Deb et al., 2023); regional splits in an allopolyploid deserve
  that reading.
- *Thresholds are provisional*, tuned on the snow carps and on simulations.
  They are all exposed as command-line options.

== Where ploidyspec sits in poly-space <sec-polyspace>

Twyford et al. (2025) propose placing polyploids in a multidimensional space
of measurable genomic features. ploidyspec provides several axes from one
assembly:

#table(
  columns: (1.4fr, 2fr),
  table.header([*Twyford et al. axis*], [*ploidyspec measurement*]),
  [Heterozygosity / segregating variation], [`allelic_distance_median`: distance between copies of the same chromosome number],
  [Sequence divergence between duplicated sets], [homeolog pair distance, `distance_ratio`, `dist_split`],
  [Structural divergence], [fusions between copies; regional lineage splits],
  [Temporal heterogeneity (diploidization)], [mix of states across chromosomes; `pair_depth_cv`; fusion `lineage_divergence`],
  [TE distribution], [`te_split` (within-genome contrast)],
)

`ploidyspec panel` builds an exploratory PCA over these numbers
(`meta/supplementary/poly_space_pca.png`). It is a picture, not a
classifier.

== Open questions

- Confirm the `SchCurv1` chr17 block against the published coordinates.
- Literature checks for the cryptic candidates (`dmRanRepe1`, `dcCerAlpi1`).
- Assess two-haplotype genomes for tetrasomy more directly: the pooled test
  can show tetrasomy, but cannot show allopolyploidy.
- Add anchors with mapped rediploidization, ideally a salmonid.
- Auto-detect chromosome and haplotype naming, so manifests need no regexes.

// ================================================================ 9
= Output reference

#table(
  columns: (2.3fr, 2.2fr),
  table.header([*File*], [*What it holds*]),
  [`sequences.tsv`, `unplaced.tsv`], [one row per chromosome copy (unit); everything else, with the reason],
  [`matrix/whole_chrom_pairs.tsv`], [distance (and k-mer counts) for every pair of units],
  [`matrix/ploidy_summary.tsv`], [copies per chromosome number and their haplotypes],
  [`matrix/chrom_label_corrections.tsv`], [relabelled and ambiguous units],
  [`matrix/*heatmap*.png`], [the full distance matrix],
  [`homeologs/homeolog_pairs.tsv`], [accepted homeolog pairs (FDR 5%, ≥ 2% below background)],
  [`homeologs/homeolog_candidates_ranked.tsv`], [every pair of chromosome numbers, ranked, with z, p, q],
  [`homeologs/ploidy_ancestry_summary.tsv`], [each chromosome's partner and `distance_ratio`],
  [`homeologs/homeolog_pairs_reordered_ heatmap.png`], [chromosome numbers clustered; pairs outlined],
  [`windowed/windowed_all.tsv`], [per window of `unit_a` against `unit_b`: `kmers_a`, `shared`, `containment`, `distance`],
  [`structure/genome_partition.tsv`], [significant k-group partitions of chromosome numbers, with z],
  [`structure/pair_synchrony.tsv`], [`n_accepted_pairs`, `pair_depth_median`, `pair_depth_cv`],
  [`subgenomes/auto_allo_index.tsv`], [`te_marker_fraction` per pair of copies],
  [`subgenomes/te_marker_fraction_ by_lineage.tsv`], [within/cross-lineage fractions and `split_ratio`],
  [`rediploidization/fusions.tsv`], [tested scaffolds: components, containment, `sibling_excess`, `lineage_divergence`, `status`],
  [`rediploidization/rediploidization_by_chrom.tsv`], [per chromosome: splits, extent, segments, `copy_state`, `state_basis`, partner],
  [`rediploidization/rediploidization_summary.tsv`], [counts per state, fusions, outlier haplotype, pair depths],
  [`summary.tsv`, `summary.md`], [the four answers with confidence and evidence],
  [`report.html`], [everything above, one page],
  [`meta/panel_summary.tsv` (panel)], [one row per species: the core numbers],
  [`meta/species_summaries.tsv` (panel)], [every species' four answers],
)

Full column-by-column documentation with worked examples is in `OUTPUTS.md`
in the repository.

// ================================================================ 10
= Glossary

/ Allele distance: k-mer distance between copies of the same chromosome number; a measure of heterozygosity.
/ Containment: fraction of one sequence's k-mers found in another.
/ Copy / unit: one assembled copy of one chromosome (`HAP2_chr05`).
/ Disomic inheritance: each chromosome pairs with a preferred partner; copies form fixed pairs.
/ Fossil TE marker: a repeat k-mer much more abundant in one copy than another, left by separate bursts of TE activity.
/ Fusion: two chromosomes joined into one in some haplotypes.
/ Homeolog: a chromosome derived from the same ancestral chromosome through polyploidy, in a different set.
/ Homoeologous exchange: recombination between homeologs, swapping or replacing segments.
/ Lineage split: copies of one chromosome falling into two groups that are closer within than between.
/ Mash distance: per-base divergence estimated from shared k-mers.
/ Ohnolog: a gene duplicated by whole-genome duplication.
/ Pooling: reading a chromosome together with its homeolog partner when each has too few copies on its own.
/ Rediploidization: return of a polyploid genome to diploid-like (disomic) behaviour, chromosome by chromosome.
/ Tetrasomic inheritance: all four copies pair at random; copies stay interchangeable.

// ================================================================ appendix
= Appendix A: species summaries <app-species>

What `ploidyspec summary` says for each of the 45 assemblies in the panel,
generated from `meta/species_summaries.tsv`. Each card gives the four answers
of Chapter 6 with their confidence and the evidence behind them. Where the
literature gives a chromosome count it is shown top right, for comparison
with the ploidy answer. Remember what these are: readings of structure in one
individual's assembly, not inheritance and not origin. A two-haplotype
assembly in particular can show that duplicated sets do not interchange, but
not why.

#include "species_summaries.typ"

// ================================================================ refs
= References

#set par(justify: false, hanging-indent: 1.2em)
#set text(size: 9.3pt)

Becher, H. et al. (2020). Tetmer: k-mer spectrum model for tetraploid genomes. As cited in Twyford et al. (2025).

Deb, S. K., Edger, P. P., Pires, J. C. & McKain, M. R. (2023). Patterns, mechanisms, and consequences of homoeologous exchange in allopolyploid angiosperms: a genomic and epigenomic perspective. _New Phytologist_ 238(6): 2284-2304. doi:10.1111/nph.18927

Gaynor, M. L., Feng, K., Soltis, D. E., Soltis, P. S. & Smith, S. A. (2026). All detectable ancient whole-genome duplications involve hybridization. _bioRxiv_ preprint, posted 28 August 2026, not peer reviewed. doi:10.64898/2026.08.25.747000

Jaron, K. S. Separate sub-genomes of an allopolyploid. k-mer approaches for biodiversity genomics (wiki). github.com/KamilSJaron/k-mer-approaches-for-biodiversity-genomics

Li, Z. et al. (2021). Patterns and processes of diploidization in land plants. _Annual Review of Plant Biology_ 72: 387-410.

Lv, Z., Addo Nyarko, C., Ramtekey, V., Behn, H. & Mason, A. S. (2024). Defining autopolyploidy: cytology, genetics, and taxonomy. _American Journal of Botany_ 111(8): e16292.

Mason, A. S. & Wendel, J. F. (2020). Homoeologous exchanges, segmental allopolyploidy, and polyploid genome evolution. _Frontiers in Genetics_ 11: 1014.

McHale, F., Mulhair, P. O. & Holland, P. W. H. (2025). Evolution of duplicated Hox gene clusters in land snails and slugs. _Journal of Experimental Zoology Part B_ 344(6): 363-368. doi:10.1002/jez.b.23322

Myers, G. FastK: a k-mer counter for high-fidelity shotgun datasets. github.com/thegenemyers/FASTK

Ondov, B. D. et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. _Genome Biology_ 17: 132.

Ranallo-Benavidez, T. R., Jaron, K. S. & Schatz, M. C. (2020). GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. _Nature Communications_ 11: 1432.

Twyford, A. D., Conover, J. L., Doyle, J. J., Mason, A. S., Soltis, D. E., Soltis, P. S. & Wendel, J. F. (2025). The polyploid continuum and the landscape of polyploid genomic variation. _American Journal of Botany_ 112(11): e70121. doi:10.1002/ajb2.70121

Wang, J., Qin, J., Sun, P., et al. (2019). Polyploidy index and its implications for the evolution of polyploids. _Frontiers in Genetics_ 10: 807. doi:10.3389/fgene.2019.00807

Xie, C., Ma, Z., Zhou, C., Ma, K., et al. (2026). Chromosomal fusions trigger rediploidization of autopolyploid genomes. _Nature_ 654: 706-713. doi:10.1038/s41586-026-10439-1
