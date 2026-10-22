// ploidyspec guide: the method, the tool, every output and plot, and every species
// Build (from the repository root):
//   python3 docs/guide/make_guide_species.py
//   typst compile --root . docs/guide/ploidyspec_guide.typ

#import "guide_lib.typ": *

#set document(title: "ploidyspec guide", author: "ploidyspec project")
#set page(
  paper: "a4",
  margin: (x: 2.2cm, top: 2.4cm, bottom: 2.4cm),
  header: context {
    if counter(page).get().first() > 1 [
      #set text(8pt, fill: grey)
      ploidyspec guide #h(1fr) October 2026
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
#text(15pt)[Reading polyploidy, its origin and its return to diploidy \ from haplotype-resolved genome assemblies]
#v(0.8cm)
#line(length: 100%, stroke: 1pt + accent)
#v(0.5cm)
#text(11pt)[The method, the tool, every output and diagnostic plot, \ worked examples, and a summary of every species run so far]
#v(2cm)
#grid(columns: (1fr, 1fr), gutter: 0.8cm,
  note(title: "What this is")[
    One place to understand the whole project: why it exists (Chapters 1–3),
    how to run it (4), how each number is computed (5–6), what every output
    file and plot means (7–8), four walkthroughs on real species (9), how it
    was validated and corrected (10–11), how to interpret a result (12),
    what it says about each of the 59 assemblies (13), and how it stands
    against the literature (14).
  ],
  note(title: "Status")[
    Reflects the code and results of 12 October 2026, after the panel was
    re-run with every correction in Chapter 11. Generated parts (Chapter 13)
    are rebuilt from the result files by `docs/guide/make_guide_species.py`.
  ],
)

#pagebreak()
#outline(title: [Contents], indent: auto, depth: 2)


// ================================================================ 1
= At a glance

Give ploidyspec the haplotype-resolved assembly of one individual (two or
more haplotype FASTA files, or one file with haplotypes named in the
headers) and it answers four questions, each with the evidence behind it and
a confidence:

#table(
  columns: (1.15fr, 2.3fr, 1.3fr),
  table.header([*Question*], [*How ploidyspec answers it*], [*Where*]),
  [*1. What is the ploidy?*],
  [Counts the haplotype copies of each chromosome (an extra chromosome set held
   in one file is split out), then looks for chromosome numbers that are
   themselves duplicated: homeolog pairs, genome partitions, and blocks of
   segmental homeology.],
  [Sections 6.3–6.7, 12],
  [*2. Auto-like or allo-like?*],
  [Whether the copies of each chromosome are interchangeable or fall into
   lineages that do not mix, consistently across chromosomes; for two-haplotype
   assemblies, through homeolog pairs and residual tetrasomy.],
  [Sections 6.9–6.11],
  [*3. What do the TE markers say?*],
  [Whether repeat (transposable element) k-mers private to one copy track the
   lineages: a within-genome contrast, not an absolute score.],
  [Sections 6.8, 8.11],
  [*4. Is there evidence of rediploidization?*],
  [Fusions between copies and the divergence spreading from them, chromosomes
   split into lineages along all or part of their length, residual tetrasomy,
   homeologous exchange and lineage switches, chromosomes rearranged since the
   duplication, and pairs or blocks diverged to different depths.],
  [Sections 6.9–6.11],
)

It does this from DNA sequence alone: no gene annotation, no reference
genome, no read data, no alignment. The core measurement is k-mer distance
between chromosome copies, whole and in windows. Chapter 14 compares this with
published methods and checks the panel's conclusions against the literature.

#caution(title: "The one thing to remember")[
  Sequence from a single individual shows how its chromosome copies relate
  to each other today. It shows inheritance-like structure, not inheritance
  itself (that needs crosses or populations), and not origin. A long
  rediploidized autopolyploid can look exactly like an allopolyploid. The
  outputs are worded accordingly: "tetrasomic-like", "resolved lineages",
  "allo-like", never "is an autopolyploid".
]


= Background: why this question, why now

This chapter summarizes the literature the method rests on. Chapter 14 goes
further: it sets ploidyspec beside every comparable method, checks our
conclusions against the published record, and lists what the literature changed
in the tool. All works are listed in Appendix C.

== Polyploidy is common, and its categories are blurry

Whole-genome duplication (WGD) has happened again and again in the history
of plants and animals. Most lineages have been through several cycles of
duplication and diploidization, and all seed plants descend from polyploid
ancestors @Doyle2017-ii @Van-de-Peer2009-hf. Phylogenomic surveys keep finding
new WGDs; duplicate genes from them are lost roughly exponentially, with a
half-life of about 21.6 million years in angiosperms @Ren2018-pj @Qiao2019-wp.
Established WGDs may be rare successes among many failed polyploids, but they
coincide with radiations and novelty @Van-de-Peer2009-hf, and much basic
biology of polyploids is still open @Soltis2010-jk.

Polyploids are traditionally sorted into two kinds. *Autopolyploids* double a
genome within one species; their chromosome copies are near-identical and
pair freely at meiosis (*polysomic*, or *tetrasomic* in a tetraploid,
inheritance). *Allopolyploids* combine the genomes of two diverged species;
each chromosome has a preferred partner from its own parent (*disomic*
inheritance), and the two parental sets remain as *subgenomes*.
#tc(<Stebbins1947-ei>) already added a third class, the *segmental
allopolyploid*, whose parents were only partly diverged, so that some
chromosomes pair disomically and some polysomically.

The two kinds are defined in two incompatible ways. *Taxonomic* (mode-of-origin)
definitions ask whether one species or two contributed the genomes; *genetic*
or *cytological* definitions ask whether inheritance is polysomic or disomic
@Doyle2017-ii. The two often disagree. #tc(<Lv2024-yo>) found pure or
predominantly polysomic inheritance in 39 of 43 putative autopolyploids
(91%), but also in seven species that are taxonomically allopolyploid, with a
further five partly polysomic; inheritance can relabel species. Complex
polyploids with more than two contributing lineages, and introgression between
cytotypes, blur the picture further @Leal2024-va.

#tc(<Twyford2025-xy>) argue that the binary is better seen as a
*continuum*, or a multidimensional "poly-space", because the definitions
disagree and because genomes move through that space over time. Their Box 2
makes the point that matters most here: after duplication a genome undergoes
*diploidization* (or *rediploidization*), and a diploidized autopolyploid can
have diverged duplicate chromosomes "indistinguishable from those in a classic
allopolyploid". They suggest measurable axes to place genomes in poly-space:
heterozygosity, sequence divergence between duplicated sets, gene-level
divergence, and structural variation, ideally from haplotype-resolved genomes.
Phased genomes are what now make these questions answerable: which part of the
genome came from which lineage, and whether the lineages have mixed
@Thomas2026-kv. ploidyspec measures several of these axes directly
(Section 12.4).

== Meiosis decides what the sequence records

What a polyploid genome looks like after some time depends on how its
chromosomes pair. In a new polyploid, related chromosomes can form multivalents
or pair with any partner; stable polyploids usually evolve diploid-like,
bivalent pairing, under genetic control that is now partly known
@Grandont2013-or @Soares2021-rh. The transition from polysomic to disomic
inheritance need not wait for dedicated pairing genes: in simulations, drift
together with a threshold for pairing fidelity is enough, inversions and
reduced chiasma numbers help, and sub- or neofunctionalization speeds it up
@Le-Comber2010-au.

The transition is rarely genome-wide at once. In tetraploid rose, pairing is
*partially preferential*, differs between chromosomes, varies along chromosome
arms, and depends on the genotype @Bourke2017-hd. In apomictic _Brachiaria
decumbens_, homologs pair preferentially across the genome but one chromosome
segregates multisomically @Worthington2016-jx. The Adriatic sturgeon shows
tetrasomic segregation at 11 of 14 loci and disomic at three
@Dalle-Palle2022-vr, and an AABB banana shows both disomic and polysomic
inheritance cytologically @Jeridi2012-vy. Mixed inheritance within one genome
is therefore expected, not exceptional. For ploidyspec this is why every
reading is made per chromosome and per window.

Inheritance itself is measured from crosses, populations or cytology:
segregation in mapping populations @Bourke2017-hd @Bicknell2023-qk, allele
dosage in families @Delaval2024-dr, approximate Bayesian computation on
population samples @Roux2015-wy, and the population-genetic tools reviewed by
#tc(<Dufresne2014-rt>). One assembled individual cannot measure inheritance;
it records the sequence consequences of past pairing.

== Rediploidization: the time axis of the continuum

#tc(<Xie2026-aq>) caught rediploidization near its start (commentary in
#tc(<Xie2026-zb>)). Snow carps (Schizothoracinae) share one autotetraploid
origin. In _Schizopygopsis younghusbandi_ most of the 25 chromosome quartets
are still tetrasomic: their four copies are about equally similar
(Ks ≈ 0.01). But on ten chromosomes, two of the four copies have fused with
another chromosome (five fusions). On those chromosomes, fused copies pair
with fused copies and unfused with unfused, so the quartet has split into two
lineages that diverge. Divergence is highest at the fusion site and spreads
outwards along the arms, and the fusions happened in three "waves" of
different ages. A separate mechanism, a centromeric inversion, made the short
arm of chromosome 17 disomic while the long arm stayed tetrasomic. The authors
also note that downstream consequences of rediploidization, such as biased
gene loss and expression between the two lineages, mimic subgenome dominance
in allopolyploids.

At the other end of the time axis, salmonids show what a duplication about 80–95
million years old looks like. Their ancestor was autotetraploid, and
rediploidization was slow and regional @Lien2016-lf. About a quarter of each
salmonid genome rediploidized only after the lineages had split
(*lineage-specific ohnologue resolution*, LORe), more than 50 million years
after the WGD, so these regions have diverged less than the rest
@Robertson2017-ya. In Atlantic salmon the most similar duplicated regions
(above 90% sequence similarity, 25.6% of the anchored sequence) cluster in
seven pairs of chromosome arms, each a large collinear block that includes a
telomere @Lien2016-lf; in
resequencing data from charr, whitefish and salmon they are where deviant SNPs
concentrate @Dallaire2023-bn. The same pattern holds on deeper timescales:
across 74 teleost genomes some duplicated regions stayed tetraploid for more
than 60 million years after the teleost WGD, which points to an autopolyploid
origin @Parey2022-ic; and sturgeon and paddlefish share one WGD that was masked
for decades because the two lineages separated before rediploidization was half
done @Redmond2023-fm. Rediploidization also reshapes karyotypes: the
"pseudotriploid" _Cardamine cordifolia_ (2n = 24) is a diploidized tetraploid
whose chromosome number fell from 32 by terminal translocations
@Mandakova2016-ul, and the mangrove _Sonneratia_ rearranged its chromosomes
extensively after a whole-genome triplication @Feng2024-ew. Gene loss during
rediploidization proceeds mainly by DNA deletion rather than pseudogenization
@Crombez2025-iz, and allele loss is already measurable in a young
autotetraploid @Yu2023-rd; #tc(<Li2021-dp>) review the patterns in plants.

This is the pattern ploidyspec was built to read: per chromosome, and
along each chromosome, whether the copies still behave as one pool or have
separated. The snow carp genomes are part of our panel, and we recover
their fusions and their wave order without being told where to look
(Section 9.2); three salmonids are in the panel too (Sections 9.4, 14.4).

== Homoeologous exchange blurs things from the other side

Allopolyploids do not stay neatly separated either. Homoeologous exchange (HE)
is recombination between chromosomes of different subgenomes
@Gaeta2010-va @Mason2020-mz @Deb2023-bo. HE can be reciprocal, or "with
replacement", where one subgenome's segment is overwritten by a copy of the
other's, homogenizing that region. It is common in young allopolyploids,
concentrated towards chromosome ends, often biased towards one subgenome, and
leaves a mosaic @Deb2023-bo. In a nascent wheat allotetraploid, HEs are biased
to low-copy, subtelomeric regions, behave like homologous recombination
hotspots, and often fall within genes @Zhang2020-ms. Repeated rounds of HE
with replacement act as a "polyploid ratchet" @Gaeta2010-va and produce
presence/absence variation @Hurgobin2018-nw; fixed exchanges can leave a
genome that looks partly auto- and partly allopolyploid, a modern route to
segmental allopolyploidy @Mason2020-mz. Segmental allopolyploid genomes have
been described in peanut @Leal-Bertioli2018-gp, sour cherry
@Wohner2023-pu and the decaploid _Thinopyrum ponticum_ @Deng2026-gt;
allotriploid AAB bananas carry frequent exchanges between A and B
@Xie2024-lk; and allopolyploid dodders show homoeologous recombination between
subgenomes @Segawa2025-tu. HE is not universal: three independent origins of
_Brachypodium hybridum_ show none @Mu2023-qq, and the young allopolyploid
_Arabidopsis suecica_ shows no "genome shock" but gradual adaptation of meiosis
@Burns2021-gz @Burns2024-md.

For ploidyspec this means an allopolyploid can show local blocks where one
homeolog carries the other's sequence. ploidyspec reads these as exchange
candidates and as lineage switches, and keeps them apart from residual
tetrasomy, where homeologs still pair in an old autopolyploid.

== Repeats and centromeres turn over fastest

The fastest-changing parts of a genome are repeats. Centromeric satellites have
diverged between the subgenomes of wheat @Su2019-kj; in _Brachypodium_,
centromeres turned over between homologous pairs from diploid to tetraploid
@Chen2024-qk; in cotton, centromeric repeats of the D progenitor invaded the A
subgenome after polyploidization @Yan2025-bz, and TE landscapes evolved
asymmetrically between subgenomes @Shi2026-ar. Repeats active in one
progenitor before the merger are therefore a subgenome label that survives
after unique sequence has diverged beyond alignment. This is the basis of
TE-signature subgenome phasing @Gordon2019-in @Session2023-vh @Jia2022-bo, and
subgenome-specific TE content is used to support allopolyploid origins in new
assemblies @Song2024-am @Robbins2026-qo @Nieuwenhuis2024-hf.

Two consequences for ploidyspec. In durum wheat our subgenome partition
recovers the A and B subgenomes exactly from the same signal those methods
use, and every chromosome is closer to the others of its own subgenome than
to its homeolog. Homeolog pairs therefore only stand out against the
cross-subgenome background, where all seven A–B pairs are found
(Section 6.4). And shared satellites
masquerade as near-identical sequence: Glechoma's apparent residual
tetrasomy was a centromeric satellite (Section 6.10).

== Divergence records isolation, not duplication

A preprint by #tc(<Gaynor2026-rc>) (bioRxiv; not yet peer reviewed) makes
an argument that applies to any divergence-based method, including ours.
Using haplotype-phased polyploid assemblies, they find that synonymous
divergence (Ks) builds up only between gene copies that do not recombine.
Copies on chromosomes that pair at meiosis stay similar, whatever the
ploidy. Duplicate-gene divergence therefore detects meiotic isolation, or
divergence inherited from hybridizing parents, rather than genome
doubling as such. That is why tetrasomic autopolyploid copies stay close.
Simulations make the same point for Ks histograms: the peak can mark parental
divergence, the last migration between parents, or the ancestral coalescence,
and can be hidden at zero @Dunn2025-qy. Gene conversion between WGD duplicates
is highest just after polyploidization and declines with time @Qiao2019-wp.

ploidyspec's distances behave the same way, and we use this deliberately.
Low distance among all copies of a chromosome means the copies still mix.
Two clusters at a consistent depth mean two lineages have been isolated
from each other. Divergence alone cannot say whether the isolation came
from hybridization (allo) or evolved later (rediploidized auto). It also
explains a blind spot: an old autopolyploidy that stayed tetrasomic leaves
little divergence to find.

== Other ways to tell auto from allo, and how ploidyspec differs

- *Gene-loss balance (P-index; #tc(<Wang2019-zs>)).* After duplication, an
  allopolyploid tends to lose genes unequally between subgenomes (biased
  fractionation), and an autopolyploid equally. Wang et al. turn this into
  a 0–1 index from gene collinearity against an outgroup genome, with
  about 0.3 separating auto from allo. It needs gene annotation, synteny,
  and a suitable outgroup. Xie et al. found biased gene loss towards the
  fused copies in an autopolyploid, and #tc(<Parey2022-ic>) biased retention
  in the autopolyploid teleost WGD, so biased fractionation can also arise
  after rediploidization.
- *Duplicated BUSCOs and gene clusters.* #tc(<McHale2025-wt>) detect the
  ancient land snail and slug WGD from the proportion of duplicated BUSCO
  genes (16–21% in Stylommatophora versus 0.5–1.5% in other molluscs) and
  follow it in the patchwork retention of duplicated Hox clusters, mirroring
  vertebrates; #tc(<Farhat2023-fp>) infer a WGD in neogastropods from synteny,
  Hox duplication and Ks. These approaches suit old duplications whose copies
  have diverged beyond the reach of k-mer comparison.
- *Gene trees and ancestral karyotypes* (#tc(<Robertson2017-ya>);
  #tc(<Parey2022-ic>); #tc(<Redmond2023-fm>)) date when each duplicated region
  stopped recombining, across many species. They need annotated genomes of
  many related species.
- *k-mer spectra from reads* (GenomeScope 2.0 and Smudgeplot,
  #tc(<Ranallo-Benavidez2020-uv>); Tetmer, #tc(<Becher2020-go>)) estimate
  ploidy, heterozygosity and subgenome divergence without an assembly, but
  summarize the whole genome in one number. Tetmer showed tetraploid
  _Euphrasia_ to be allopolyploid with about 5% subgenome divergence
  @Becher2020-go; the same read-based profiling detects polyploidy across
  Microsporidia @Khalaf2024-xo. #tc(<Jenike2025-uu>) review k-mer methods for
  biodiversity genomics.
- *Subgenome phasing* assigns chromosomes of a known allopolyploid to
  subgenomes: by repeat signatures @Gordon2019-in @Session2023-vh, by
  subgenome-specific k-mers (SubPhaser, @Jia2022-bo, which fails on
  autopolyploids by design), or by differential divergence against related
  genomes (M3, @Lv2025-te). #tc(<Zhang2023-va>) benchmark these;
  #tc(<Session2024-sl>) reviews the field. They assume allopolyploidy and
  answer "which subgenome", not "auto or allo".
- *Allele dosage and segregation.* Gaynor et al. and Xie et al. use
  genotype classes (AAAB/ABBB versus AABB) or resequencing of many
  individuals to infer disomic versus polysomic inheritance directly.

ploidyspec sits between these. It uses chromosome-scale, haplotype-resolved
assemblies of one individual, as produced at scale by the Darwin Tree of
Life project, and reads structure chromosome by chromosome and along
chromosomes, without annotation. It complements gene-based methods rather
than replacing them. Section 14.1 sets them side by side.

== Haplotype-resolved assemblies make this possible

Polyploid genomes were long assembled as collapsed mosaics
@Kyriakidou2018-xv. Long reads, Hi-C and new phasing methods now separate
haplotypes in polyploids @Zhao2025-yc: a tetraploid potato with four
haplotypes, where about 50% of the genome turned out identical by descent and
about 20% structurally rearranged @Sun2022-tc; hexaploid sweet potato
@Yang2017-zi; the amphitriploid (AAABBB) gibel carp @Wang2022-am; AAB bananas
@Xie2024-lk; octoploid strawberry @Song2024-am; autopolyploid snow carps
@Xie2026-aq; a phased masu salmon @Wu2026-zm. The Earth BioGenome Project
aims to sequence all eukaryotes @Exposito-Alonso2020-qh and has set standards
for phased reference assemblies @Howe2026-nw. ploidyspec is designed for
exactly this stream of assemblies, and treats the assembly as part of the
measurement: reference-guided scaffolding @Alonge2019-fr can place a contig on
the wrong homeolog, which is how masu's first residual-tetrasomy signal arose
(Section 9.4).

== What we are trying to achieve

=== The aim

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

=== Design principles

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

= Auto, allo and rediploidization: what the data can show

== Two axes, not one label

A polyploid genome has two copies or more of a whole chromosome set. The
classical labels sort polyploids by *origin*: an *autopolyploid* duplicates one
genome, an *allopolyploid* joins two diverged genomes. They differ in what
happens at meiosis. Autopolyploid copies are interchangeable and pair at
random (*polysomic*, in a tetraploid *tetrasomic*, inheritance); allopolyploid
copies pair with their own parent's partner (*disomic* inheritance), so the
parental sets persist as *subgenomes*.

#tc(<Twyford2025-xy>) argue this is a continuum, and their key point for us is
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

    [Delayed (LORe-like) rediploidization],
    [Homeolog blocks far less diverged than the rest but not near-allelic: arms that stopped recombining late (salmonid LORe).],
    [homeology map: `n_lore_blocks`, `lore_bp`, `lore_pairs`; block `low_divergence`],

    [An extra divergent set (AAB, AAAB)],
    [One copy of every chromosome stands apart from the others.],
    [`one_divergent_copy`; extra-set split (Section 6.3)],
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
+ *Divergence records isolation, not duplication* @Gaynor2026-rc:
  copies that still recombine stay similar. A tetrasomic autopolyploid leaves
  little to see beyond equal distances. A pair of number sets that is *not*
  detected is therefore one of two opposite things: copies so similar they
  still recombine (then they show as near-identical copies, not as missing
  pairs), or copies diverged beyond k-mer reach.
+ *Exchange and residual tetrasomy are mirror images.* Both make homeologs
  locally near-identical. The difference is which copies agree: in residual
  tetrasomy all copies are close to each other; in an exchange one copy carries
  its homeolog's sequence and is far from its own homolog. ploidyspec now
  tests both (Section 6.10) because the first version confused them.
+ *k-mers see homeology only up to about 8–10% divergence.* Glechoma's
  homeologs (0.042 apart) are mapped over 79% of the genome; the ~90-My-old
  salmonid duplication over 16–30%, and Cerastium's homeologs (~0.08 apart)
  are found by the whole-chromosome test but not by the window map. In durum
  wheat, whose genome is mostly fast-turning-over transposons, the seven A–B
  pairs (0.080–0.082 apart) are found only at whole-chromosome scale, and
  only against the cross-subgenome background; the window map finds no
  blocks.
  "Not detected" never means "not there".

#caution(title: "Assembly artefacts look like biology")[
  Missing sequence reads as divergence; a contig placed on the wrong homeolog
  reads as exchange or residual tetrasomy; one haplotype file numbering its
  chromosomes differently reads as a broken genome. Much of this guide's
  algorithmic detail is about guarding against these. The masu walkthrough
  (Section 9.4) shows what happens when a guard is missing.
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
  caption: [Section 6.11 gives the full decision logic.],
)

// ================================================================ 3

= Running ploidyspec

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

== The pipeline

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
/ State: the per-chromosome reading (`tetrasomic_like`, `resolved_lineages`, ...), Section 6.9.
/ Block: a stretch of a chromosome whose closest other chromosome is consistently the same one (Section 6.6).


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
    under which a k-mer survives with probability $(1-d)^k$ @Ondov2016-jo.
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
register after the first indel larger than a window (the guide's Section 12.1).
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
- *Benjamini–Hochberg FDR* for homeolog pairs (Section 6.4).
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

The all-vs-all distance matrix (Section 5.2), then three checks on chromosome
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
    it, it is *displaced* there (the label had no support). If the move vacates
    the occupant's matched slot, costs it ≤ 5% and improves the two units'
    summed distance ≥ 1.5×, the pair is *swapped* (evidence on one side only:
    `lpElePalu1` HAP2 chr12 is 3.2× closer to chr11 while HAP2 chr11 is tied).
    Otherwise the whole file's batch is rejected as `ambiguous`, never partly
    applied.
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
The weevil `icStrMela3` (Section 9.3) is the case: HAP1 holds A as chr1–10 and
B as chr11–20; each B chromosome is 0.047–0.066 from one A number against a
median 0.093–0.099 to the rest (1.5–2.4×). A single private chromosome (the
Monacha fusion leftover) is never a set.

== homeologs: whole-chromosome pairs

#algo("homeolog pairs")[
  + For every pair of chromosome numbers (i, j): the mean distance over all copy
    combinations, so each number is one node whatever its copy count.
  + *Size adjustment*: regress these distances on log(L_i) + log(L_j) and
    |log L_i − log L_j| (three rounds, dropping pairs beyond 3 robust SD so
    homeologs do not shape the fit) and keep each pair's residual around the
    mean. Pairs are tested, and the partition and pooling are run, on these
    `adjusted_distance` values.
  + Fit the background: mean and SD of all pair distances, σ-clipping low values
    (2σ, 3 rounds) so a genome full of homeologs does not hide its own signal.
  + z = (d − mean) / SD; one-sided normal p; Benjamini–Hochberg q.
  + Accept pairs with q ≤ 0.05 *and* d ≤ 0.98 × median background (`--min-effect`
    2%), greedily one partner per number, closest first.
  + *Subgenome partition*: cut the average-linkage tree of the raw
    number-to-number distances at two groups, and fit the size model with a
    same-set term added. The partition is accepted when the same-set
    coefficient (in robust SD) is at least 10 SD above 200 random splits of
    the same group sizes (z ≥ 10) and the smaller set holds ≥ 25% of the
    numbers.
  + *Cross-subgenome test*: with an accepted partition, remove only the size
    terms (keeping the set effect) and test the cross-set pairs against the
    cross-set background, as above. The result replaces the genome-wide test
    only if it finds more pairs (`pair_test` in `subgenome_partition.tsv`).
]
The size adjustment exists because whole-chromosome distance depends on the
chromosomes' sizes: it correlates 0.8 with the size difference in Cerastium
and in the diploid _Malva_. Without it, clustering sorted chromosomes into
large and small; adjusted, Cerastium resolves into 18 pairs.

But the size adjustment alone also removed real subgenomes. In an
allopolyploid each subgenome carries its own repeats, so non-homologous
chromosomes of one subgenome can be closer than homeologs across them
(@fig-ranunculus). Tested against every other pair in the genome, the
homeologs then do not stand out. When one subgenome is also larger, the size
regression absorbs the split. _Ranunculus repens_ (2 × 16, 4x on x = 8
@Karbstein2025-ou) is the case: A (chr01–08, 4.7 Gb) and B (chr09–16, 3.8 Gb),
within-set distances 0.064–0.068, homeologs 0.070–0.073, other A–B pairs
0.078–0.081. The genome-wide test found no pair. The subgenome partition
(z = 109) and the cross-subgenome test find all eight (1–9, 2–15, 3–10,
4–12, 5–13, 6–11, 7–16, 8–14; z −13 to −17), and every A chromosome is larger
than its homeolog. The same change recovers durum wheat's seven A–B pairs
(sets exactly A and B, z = 179) and eight pairs in _Inula conyzae_ (z = 136).
No species without subgenomes reaches z = 5 (_Malva_ 2.6), and Cerastium
(z = 35, sets misplacing two chromosomes) keeps the genome-wide test's 18 pairs.

#figure(
  image("fig/ranunculus_subgenomes.png", width: 100%),
  caption: [_Ranunculus repens_: number-to-number distances ordered by
  subgenome (white boxes) with homeolog pairs (red boxes), and the 8 × 8
  A-versus-B block. Within a subgenome chromosomes are closer than homeologs
  across subgenomes; each A chromosome is closest to one B partner.],
) <fig-ranunculus>

The 2% effect floor exists because a tight background makes tiny
differences significant: on simulated diploids, a pair 1.25% below the median
reached z = −4. Outputs: `homeolog_pairs.tsv`, every candidate in
`homeolog_candidates_ranked.tsv`, and per chromosome in
`ploidy_ancestry_summary.tsv` (`distance_ratio` = homeolog distance / own
allelic distance).

== windowed: along each chromosome, same number

For each chromosome number, every copy is windowed against every other copy of
that number (Section 5.4, k = 15). These tracks carry the lineage splits along
chromosomes (Section 6.9) and each copy's own-homolog distance per window,
which the residual and exchange tests need (Section 6.10).

== windowed-homeologs: the homeology map, then pair tracks

This stage first maps homeology genome-wide, then computes windowed tracks
between the copies of every homeolog pair (whole-chromosome pairs plus block
pairs from the map), plus control tracks.

#algo("homeology map")[
  + Take the haplotype with the most chromosomes. Window every chromosome
    (250 kb) against every *other* chromosome of that haplotype, at k = 23
    (Section 5.4; n² profiles).
  + *Size correction*: a larger chromosome contains more of any window's k-mers.
    Fit one slope of window distance on log(target length) over all (window,
    target) pairs, dropping outliers, and subtract it.
  + *Cross-subgenome mode*: when the homeologs stage used the cross-subgenome
    test, each window is scored only against the other subgenome's
    chromosomes, so its own subgenome's shared repeats cannot win.
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
  + *Low-divergence (LORe-like) blocks*: with at least 6 reciprocal blocks,
    a block is LORe-like when its median distance lies between 0.25× and
    0.6× the median of all reciprocal blocks, and the pair is low in both
    directions. Below 0.25× it is `near_allelic`, left to the
    residual-tetrasomy test and its repeat guards.
]
Why the size correction: in the repeat-rich _Hesperis matronalis_ (six
chromosomes of up to 875 Mb, every window 0.022–0.026 from every chromosome)
the closest chromosome was simply the largest, producing chr01/chr02 "blocks"
over 21% of the genome; in Arctic charr, ten unmatched blocks pointed at its
largest chromosome, which alignment showed shares only scattered short hits
with them. A per-target baseline would also remove this, but absorbs real
homeology where many chromosomes are related (_Lathraea_), so a single slope is
used. Why k = 23: at k = 15 about 40% of a window's k-mers occur by chance in any
100 Mb chromosome, and salmonid homeologs barely stand out (masu 0.049 vs 0.057
for unrelated chromosomes). Why reciprocity: homeology runs both ways, repeats
need not. Without it, Glechoma gained two short blocks onto chr16 at a third of
the homeolog distance (shared repeats) and the diploid _Malva sylvestris_ had
seven blocks; with it, Glechoma has exactly its nine known pairs (78% of the
genome) and _Malva_ none.

Why a LORe-like class: in salmonids about a quarter of the genome
rediploidized tens of millions of years after the rest @Robertson2017-ya
@Lien2016-lf. Those arms stopped recombining late, so they are less diverged
than other homeologs, but they are no longer near-identical, so the
residual-tetrasomy test (Section 6.10) is blind to them by design. The map
sees them: in whitefish 13 reciprocal blocks (243 Mb, 40% of the block span)
sit at 0.011–0.026 against a median block distance of 0.045. The band's lower bound keeps
near-allelic blocks out (they are residual tetrasomy or repeats:
_Lathraea clandestina_'s chr03/chr04 blocks at 0.002 against a median of 0.016
are repeats, Section 6.10), and the two-way rule drops one-way low blocks. A
second, younger duplication looks the same, so the reading names both.

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
candidate subgenome structure; the rediploidization stage pools chromosomes
across it (Section 6.9). For _Inula_ the genome-wide pair test found nothing,
but the size-independent subgenome partition of the homeologs stage
(Section 6.4) gives the same sets and eight cross-set pairs.
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

=== Lineage switches

#algo("lineage switches (three or more copies)")[
  + Along each copy, label every window with its closest other copy of the same
    number; runs of one label (2-window gaps bridged) covering ≥ 10% of the
    chromosome and ≥ 2 Mb are kept.
  + Two consecutive runs with neighbours X then Y are a *switch* when over the
    first the copy is ≥ 2.5× closer to X than to Y (median over windows), over
    the second the reverse, and the far copy sits ≥ 3× the genome's allelic
    level away (per chromosome the most allelic copy's median closest
    distance; median over chromosomes).
]
The exchange readout compares different chromosome numbers; when subgenome
copies share a number, an exchange shows as a switch. The weevil's HAP1 chr09
is A to 40 Mb and B after it. In polysomic genomes the nearest copy changes in
short blocks between equally close copies, which the contrast and lineage-level
rules exclude.

=== Divergence along fused chromosomes

Each confirmed fusion is windowed along the fused sequence against the unfused
copies of its components (`d_unfused`, labelled by component) and the other
fused copy (`d_fused`). The junction is the boundary that best separates the
two components' windows; `junction_gradient` = median `d_unfused` within the
inner 20% of each arm over the outer 20%. Above 1 means divergence peaks at the
fusion site @Xie2026-aq.

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

Both read the pair tracks of Section 6.6 against each copy's own-homolog track
from Section 6.5, for every homeolog pair (whole-chromosome or block).

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
      chromosome carrying them);
    - the window is not low-complexity: its sequence compresses (zlib) to at
      least 0.75 of the genome's median ratio.
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
#caution(title: "Why the complexity guard matters")[
  Position-free containment cannot tell shared repeats from shared sequence.
  Glechoma's chr14/chr17 block (16–22 Mb) read as near-allelic in all eight
  copy directions, but it is a centromeric satellite: it compresses to
  0.06–0.13 against 0.32 for the rest of the chromosome and aligns at only
  94–96%. Homeologous centromeres inherit one satellite family, so the
  unrelated control does not remove it. _Lathraea clandestina_'s chr03/chr04
  segment is repeat-dense (0.19–0.27 against a median 0.29) and aligns at
  98–99% against 99.9–100% for the alleles. Pilosella's segments, by contrast,
  align at the allelic level (99.6–99.8%) and are kept.
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
+ Two copies and no duplication evidence (no sets, fewer than 3 whole-
  chromosome pairs at |z| ≥ 5, and fewer than 3 map pairs): not assessable
  (gar and bowfin pass one pair each; Stagnicola one strong and two marginal).
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
(CV ≥ 0.15) and LORe-like blocks (at least one reciprocal pair, ≥ 4 Mb) are
appended as "Also:". Residual tetrasomy ≥ 90% at chromosome ends is low
confidence (homeolog-specific subtelomeric repeats pass the control), although
real delayed rediploidization and homoeologous exchange are also biased
towards chromosome ends @Lien2016-lf @Zhang2020-ms, so terminal signal is
flagged, not discarded.

== report, panel, simulate, cleanup

`report` writes one self-contained HTML page per species. `panel` stacks all
species into `panel_summary.tsv`, `rediploidization_panel.tsv`,
`species_summaries.tsv`. `simulate` writes synthetic assemblies with a known
answer (Chapter 10). `cleanup` deletes k-mer tables and keeps every result.

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
  [`homeologs/subgenome_partition.tsv`], [the two subgenomes, their z against random splits, and which pair test was used],
  [`homeologs/homeology_map.tsv`], [per window: `best`, `best_dist`, `median_dist`, `z`, `homeolog_like`],
  [`homeologs/homeolog_blocks.tsv`], [per block: `chrom`, `partner`, `start`–`end`, `median_dist`, `position`, `reciprocal`, `low_divergence`],
  [`homeologs/homeology_map_summary.tsv`], [`duplicated_frac`, `n_partner_pairs`, `multi_partner_list`, `block_dist_cv`, `n_lore_blocks`, `lore_bp`],
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

= Reading the diagnostic plots

#include "plots_chapter.typ"


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
copies); every pool splits by `dist_split` 3.5–4.6, so all 18 are
`resolved_lineages`. The residual test finds nothing. An earlier version called
chr14/chr17 residual-tetrasomic over a mid-chromosome block (figure below);
that block is a centromeric satellite and is now excluded (Section 6.10).

#figure(image("../../results/daGleHede1/homeologs/windowed_chr14x17.png", width: 100%),
  caption: [chr14 × chr17, all eight copy directions. Outside ~14–23 Mb the
  homeologs are 0.03–0.05 apart; inside, k-mer containment says ~0.003, as close
  as alleles. Alignment disagrees (94–96% identity there, for the allele as well
  as the homeolog), and the block compresses 3–5× better than the rest of the
  chromosome: a satellite array shared by the two homeologous centromeres, not
  retained pairing.])

*Summary.* Ploidy: 2 copies × 18 numbers in sets of 2, x ≈ 9, up to 4x.
Origin: allo-like (disomic-like), high. Rediploidization: diverged throughout.

*What to take away.* Glechoma is the textbook allo-like genome: a single
synchronous duplication, every chromosome paired, no stretch where homeologs
are as close as alleles. Its one apparent exception taught the tool to ignore
repeat-dense windows.

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
carp mechanism), high. This is the published result of #tc(<Xie2026-aq>),
recovered without being told where to look. Aligning HAP3 chr17 to HAP1 shows
why chr17 splits only at its end: HAP3 is inverted from ~27 Mb to the end, the
inversion mechanism Xie et al. describe for this chromosome.

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

*The odd one out.* On chr09 the outlier is HAP2, not HAP1B. Along HAP1 chr09
the windows show why: it matches the other A copy over 0–40 Mb (0.004 against
0.015 to B) and the B copy from 40 Mb to its end (0.004 against 0.018 to A), so
chr09 reads AAB on one half and ABB on the other. The lineage-switch readout
(Section 6.9) reports it as `HAP1:0.0-40.0Mb~HAP2;40.0-81.4Mb~HAP1B`: a
homeologous exchange with replacement, or a phasing switch between the two sets
held in HAP1's file (Hi-C would tell).

== Oncorhynchus masou: when a guard is missing

*Input.* Two haplotypes × 33. HAP2 was contig-level and was scaffolded on HAP1
for this project (`workflows/prep/scaffold_by_reference.py`): its contigs were
placed where they align best on HAP1.

*The first result.* Residual tetrasomy on 8 of 10 paired chromosomes, 54 Mb,
62% terminal: an apparently textbook salmonid. #tc(<Lien2016-lf>) describe such
regions in Atlantic salmon.

*What was wrong.* In a region where two homeologs are similar, a HAP2 contig
from chromosome B aligns almost as well to A and can be placed there. HAP1 A is
then near-identical to "partner" HAP2 B (which carries A's sequence), and far
from its own HAP2 A (which carries B's). Three of the four copies see a
near-allelic stretch and the pair vote passed. The tell: those windows' own
homolog distance was 0.03–0.085, homeolog distance.

*The current result.* With the own-homolog condition and the repeat guard:
residual tetrasomy on one pair, chr17/chr29 (8.6 Mb, seen from all 4 copies);
exchange candidates on 3 chromosomes, each in 1 of 2 copies (the swapped
contigs, now read as what they are). The size-corrected homeology map finds
16.5% of the genome in blocks, 17 chromosome pairs, and 8 chromosomes with
blocks on two or three others (chr32: chr01+chr04+chr33): the reshuffled
salmonid karyotype. One residual pair is not enough to call an auto-like
origin, so the origin reads allo-like (medium).

*What to take away.* Our own preprocessing produced the strongest
"biological" signal in the panel. Masu remains a useful salmonid, but its
residual numbers carry that history; Arctic charr and whitefish, assembled
with both haplotypes scaffolded, are the cleaner salmonid tests. The k-mer map
reaches 18–30% of their genomes; neither shows residual tetrasomy, but both
show LORe-like blocks (Section 6.6): whitefish 243 Mb in five chromosome pairs,
charr 8.5 Mb on one pair (chr09/chr11), at chromosome ends, as #tc(<Lien2016-lf>)
found the delayed regions of salmon to include telomeres.

*Since then.* A phased Hi-C assembly of masu with both haplotypes at
chromosome level has been published @Wu2026-zm (GenBank GCA_052849675.1 and
GCA_052849655.1). It replaces our self-scaffolded pair as `OncMaso1_wu2026`
(`workflows/prep/masu_wu2026_prep.sh`). Its HAP2 numbers its chromosomes
differently from HAP1: renumbering corrected 23, displaced 2 (HAP2 chr17 →
31, chr23 → 25) and left HAP2 chr32 ambiguous (2.9×). The first renumbering
run blocked these corrections because a few chromosomes had no allele in the
other haplotype (Section 6.3).

On the published assembly: residual tetrasomy on none of 28 tested
chromosomes (the self-scaffolded chr17/chr29 signal does not reappear), no
LORe-like blocks, and 6 near-allelic blocks that the residual test's guards
reject. The homeology map finds 22% of the genome in 63 blocks, 17
chromosome pairs, and 7 chromosomes with blocks on two or three partners
(chr08: chr07+chr16+chr26). Exchange candidates sit on 13 chromosomes
(166 Mb), almost all in one of the two copies; on a Hi-C-phased assembly
they are either real exchanges or phasing switches between homeologs, and
need alignment before either reading. The origin reads allo-like (medium) and
rediploidization advanced. Set against Atlantic salmon, where about a quarter
of the genome rediploidized late @Lien2016-lf, the k-mer map does not see
masu's delayed regions as a distinct low-divergence class. Chapter 13 shows
both masu assemblies.

#figure(image("../../results/fSalAlp3/homeologs/homeology_map.png", width: 82%),
  caption: [Arctic charr after the size correction: reciprocal blocks (solid)
  cover 18% of the genome in 18 chromosome pairs, and 6 chromosomes carry
  blocks on two homeologs (salmonid fusions since the duplication). Before the
  correction, ten unmatched blocks pointed at chr01, the largest chromosome;
  alignment showed only scattered short hits, a target-size artefact.])

// ================================================================ 7

= Validation

== Simulated genomes

`ploidyspec simulate` writes random-sequence genomes with transposable-element
families at high copy number; every haplotype except HAP1 carries a
transposition and an inversion, so copies are not collinear. The end-to-end
test (22 tests) runs the whole pipeline on each scenario.

#table(
  columns: (1.3fr, 2.7fr, 1.4fr),
  table.header([*Scenario*], [*What it contains*], [*Checked*]),
  [diploid], [two haplotypes], [no pairs, no blocks],
  [autotetraploid], [four equivalent haplotypes], [all `tetrasomic_like`, no blocks],
  [autotetraploid, 2 hap], [the same, each base chromosome numbered twice], [pairs (i, i+n), pooled `tetrasomic_like`],
  [allotetraploid], [two progenitors, own TE bursts, 2 haplotypes], [`resolved_lineages`, no residual tetrasomy],
  [rediploidized], [whole and 40% lineage splits, a placed and an unplaced fusion], [states, both fusions, and the fusion junction at the middle of each fused copy],
  [mislabelled], [one file shifted by one, a swap in another], [all relabelled],
  [residual tetrasomy], [diverged homeologs except the last 60% of two pairs; a satellite at every chromosome start], [residual on exactly those four chromosomes; satellite excluded; "auto-like origin"],
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
  [`SchCurv1`, `SchYoun1`], [autotetraploid snow carps, fusions @Xie2026-aq], [all published fusions, in wave order; mixed with fusions],
  [`drLytSali1`], [tetrasomic by classical genetics], [13/15 `tetrasomic_like`],
  [`icStrMela3`], [triploid, AAB (k-mer ploidy plot)], [3 copies, AAB-like (high)],
  [`lpTriTurg1`], [durum wheat, AABB], [the subgenome partition splits the 14 automatically numbered chromosomes into exactly the 7 A and 7 B (z = 179), and the cross-subgenome test finds all seven pairs, 1A–1B … 7A–7B (0.080–0.082, z −11 to −15), blind to the names; no window-map blocks; allo-like],
  [`dmRanRepe1`], [4x on x = 8 @Karbstein2025-ou], [two subgenomes of 8 (z = 109, A 1.23× B), eight one-to-one pairs; window map 7 of 8 pairs, 6% of the genome; allo-like],
)

// ================================================================ 8
== Fusion profiles: divergence spreads from the fusion site

#tc(<Xie2026-aq>) report that in snow carps divergence between fused and
unfused copies is highest at the fusion site and spreads along the arms. The
fusion profiles (Section 6.9) measure this directly, without being told where
the junctions are. Every one of the eleven fused copies is placed with its
junction found at 94–100% agreement, and the gradient follows fusion age:

#table(
  columns: (0.9fr, 1.3fr, 1.2fr, 1.2fr),
  table.header([*Wave (Xie et al.)*], [*Fusion*], [*Overall divergence*], [*Junction gradient*]),
  [oldest], [chr19+chr22 (`SchYoun1`, `SchCurv1`)], [0.035–0.054], [2.7–4.1],
  [middle], [chr04+15, chr08+16, chr20+23], [0.026–0.029], [2.1–2.6],
  [youngest], [chr11+chr14], [0.015], [0.6–0.7 (flat)],
  [-], [_Monacha_ chr23+24 (one haplotype)], [0.009], [0.8 (flat)],
)

Older fusions have both diverged further and concentrated that divergence at
the junction; the youngest wave has not yet. The fused copies themselves
stay within 0.001–0.008 of each other along their whole length.


= What we got wrong, and fixed


Checking results by hand found problems in three rounds (7 October, 9–11
October, and 21 October). All are fixed and covered by tests; they are listed because earlier
write-ups quoted numbers that depended on them.

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
   evidence (Section 6.9).],

  [Renumbering cycles],
  [The whitefish's HAP2 numbering differs from HAP1 on 30 of 40 chromosomes in
   long cycles; one unit with no allele blocked a cycle and the whole correction
   was rejected. Eleocharis HAP2 chr11/chr12 were swapped with evidence on one
   side only.],
  [Displacement and one-sided swaps (Section 6.3).],

  [A haplotype holding two sets],
  [The AAB weevil's HAP1 holds A and B under different numbers, so it read as
   2 copies of 20 chromosomes and "a phasing problem".],
  [Extra sets split out as a third copy (Section 6.3).],

  [Residual tetrasomy from artefacts],
  [Masu's 54 Mb rested on homeolog-swapped contigs in its reference-scaffolded
   HAP2; Glechoma's block was a centromeric satellite; Lathraea's was
   repeat-dense; one pair (Glechoma) was enough to call an auto-like origin.],
  [Own homolog must be near; repeat-dense windows excluded; origin needs two
   pairs (Sections 6.10, 6.11).],

  [Chromosome size],
  [Whole-chromosome distance tracks chromosome size (correlation 0.8 with the
   size difference), so partitions sorted chromosomes by size; in the homeology
   map the largest chromosome won every window in repeat-rich genomes.],
  [Size-adjusted pair distances; one size slope removed in the map (Sections
   6.4, 6.6).],

  [Size adjustment hiding subgenomes],
  [The size regression absorbed a split into a larger and a smaller
   subgenome, and the genome-wide pair test compared homeologs against
   within-subgenome pairs that share repeats and are closer. Ranunculus,
   durum wheat and _Inula_ lost every pair (the first version had found
   Ranunculus's eight).],
  [Size-independent subgenome partition and a cross-subgenome pair test;
   map windows scored against the other subgenome only (Section 6.4).],

  [Exchange between copies of one number],
  [The exchange readout only compared different chromosome numbers.],
  [Lineage switches (Section 6.9).],

  [Silent on delayed rediploidization],
  [Salmonid-style delayed (LORe) regions are less diverged than other
   homeologs but not near-allelic, so neither the residual test nor the summary
   mentioned them, although the map showed them (Chapter 14).],
  [LORe-like block class in the homeology map (Section 6.6).],

  [Masu assembly],
  [Our masu HAP2 was contig-level and scaffolded on HAP1 by us; its swapped
   contigs produced most of the first residual signal.],
  [Re-run on the published phased assembly @Wu2026-zm (Section 9.4).],
)

Several inheritance-mode metrics from earlier work (`partition_consistency`,
run length, `flip_rate`, the windowed CVs and a combined score) were
retired, because labelling and window registration drove them. Their last
values are archived in `meta/archive/` with a note.


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
  [homeolog blocks far less diverged than the rest, both directions], [delayed (LORe-like) rediploidization of those regions, or a younger duplication],
  [two subgenomes, chromosomes closer within a subgenome than to their homeologs], [allopolyploid: each subgenome carries its progenitor's repeats (durum wheat, Ranunculus, _Inula_); pairs show only against the cross-subgenome background],
  [a strong partition but no detectable homeolog pairs], [subgenomes separated by repeat history after unique sequence diverged beyond k-mer reach],
  [two copies, no pairs, no partition, in a known polyploid], [the second set is diverged beyond k-mer reach, or the assembly holds one set per haplotype; not evidence of recombining copies],
  [one copy apart, same haplotype always], [assembly/phasing problem, or a split-out extra set (AAB)],
  [one copy apart, haplotype changing], [AAAB-like, a divergent extra genome],
  [exchange in one copy], [homeologous exchange, or a contig on the wrong homeolog],
)

== Limits

- One individual: structure, not inheritance; not origin.
- The assembly is part of the measurement (Section 3.3).
- Old duplications fade below k-mer detection (~8–10% divergence).
- Thresholds are provisional (tuned on simulations, snow carps, Glechoma,
  _Malva_) and every one is a command-line flag (Appendix A).
- Holocentric genomes (_Eleocharis_) break and fuse chromosomes readily
  @Hipp2009-hz @Bures2013-zm @Mandrioli2020-hc, so chromosome numbers need not
  track ancestral units, and arm- and junction-based readings are weaker.
- Inheritance is per chromosome and can be partial or genotype-dependent
  @Bourke2017-hd @Dalle-Palle2022-vr; "allo-like" on a chromosome means its
  copies do not mix today, not that the species is disomic everywhere.
- Repeat signals (TE split, partition in wheat) separate subgenomes by the
  same principle as TE-signature phasing @Session2023-vh; they are not
  independent of repeat turnover at centromeres @Su2019-kj @Chen2024-qk.

== Where ploidyspec sits in poly-space

#tc(<Twyford2025-xy>) propose placing polyploids in a multidimensional space
of measurable genomic features. ploidyspec provides several axes from one
assembly:

#table(
  columns: (1.4fr, 2fr),
  table.header([*Twyford et al. axis*], [*ploidyspec measurement*]),
  [Heterozygosity / segregating variation], [`allelic_distance_median`: distance between copies of the same chromosome number],
  [Sequence divergence between duplicated sets], [homeolog pair distance, `distance_ratio`, `dist_split`],
  [Structural divergence], [fusions between copies; regional lineage splits],
  [Temporal heterogeneity (diploidization)], [mix of states across chromosomes; `pair_depth_cv`; fusion `lineage_divergence`; LORe-like blocks],
  [TE distribution], [`te_split` (within-genome contrast)],
)

`ploidyspec panel` builds an exploratory PCA over these numbers
(Section 8.13). It is a picture, not a classifier.

== Open questions

Settled on 10–11 October (by alignment, read k-mer spectra, or a fix):
Glechoma chr14/chr17 (satellite), Lathraea chr03/chr04 (repeats; 98–99% vs
99.9% alleles), Pilosella chr04/chr05 (real: homeolog = allele identity),
weevil chr09 (A to 40 Mb, B after: a lineage switch), snow carp chr17 (an
inversion in HAP3 from ~27 Mb), Monacha (HAP2 chr24 + chr23 joined at
~37 Mb), drAriEdul1 chr01 (a 19 Mb divergent HAP4 block), _Salix triandra_
(triploid) and _Rosa spinosissima_ (tetraploid, AABB), Hesperis chr01/chr02
and charr chr01 (size artefacts), Cerastium (18 large–small homeolog pairs:
subgenomes of different chromosome size). Settled on 21 October: Ranunculus
(two subgenomes of 8 and eight homeolog pairs, hidden by the size adjustment
and by subgenome-specific repeats until the cross-subgenome test).

Still open:
- Xie et al.'s exact coordinates for the `SchCurv1` chr17 inversion.
- Primary chromosome counts for Cerastium and Ranunculus.
- Hi-C to tell biology from assembly for the weevil chr09 switch and the
  Monacha join.
- Gene-based methods for duplications beyond k-mer reach (salmonids, snails).
- The open questions raised by the literature (Section 14.6).

// ================================================================ 9

= Every species so far

Each card below is generated from the result files
(`docs/guide/make_guide_species.py`). The table gives the overview; the
cards give, per species: what is known, the assembly and any numbering fixes,
distances, duplication by each method, states, fusions, residual tetrasomy and
exchange with their support, the four answers with confidence and evidence,
and caveats. "Dup. (map)" is the share of the genome in reciprocal homeolog
blocks; "–" means the map has not run (one haplotype) or the run predates it.

#include "guide_species.typ"

// ================================================================ 14

= The literature: where ploidyspec stands

#include "literature_chapter.typ"

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
  [subgenome partition z, smaller set, permutations], [10, 25%, 200], [subgenome species 109–179; no other species above 5 (_Malva_ 2.6)],
  [summary: residual pairs, duplication pairs, map fraction], [2, 3, 30%], [one pair is not an origin; gar and bowfin pass one pair],
  [LORe-like band, min blocks], [0.25–0.6× median block distance, 6], [below 0.25× is the residual test's (near-allelic); whitefish LORe-like blocks 0.011–0.026 against a median of 0.045],
  [LORe-like reading], [≥ 1 two-way pair, ≥ 4 Mb], [one short low block is noise],
)


= Appendix B: glossary

/ Allelic distance: distance between copies of the same chromosome number.
/ Block: a run of windows whose closest other chromosome is consistently the same.
/ Containment: fraction of one sequence's k-mers present in another.
/ Disomic / tetrasomic: pairing with a fixed partner / at random among four.
/ Extra set: a second chromosome set held in one haplotype file under other numbers.
/ Homeolog: a chromosome derived from the same ancestral chromosome through polyploidy.
/ Homeologous exchange: a stretch of one homeolog replaced by the other's sequence.
/ LORe: lineage-specific ohnologue resolution; regions that rediploidized late, after speciation (salmonids).
/ Mixed inheritance: some chromosomes or loci segregate disomically, others polysomically, in one species.
/ Pooling: reading two homeologous numbers' copies together.
/ Reciprocal block: a block whose partner has a block back on it.
/ Residual tetrasomy: regions of an otherwise rediploidized genome where homeologs still behave as one pool.
/ Rediploidization: the return of a polyploid to disomic behaviour, region by region.
/ Segmental allopolyploid: a polyploid whose subgenomes are only partly diverged, so pairing is partly preferential (Stebbins 1947).


= Appendix C: references

All works cited in this guide and the further literature reviewed for
Chapter 14 (86 works).

#set text(size: 9pt)
#bibliography("references.bib", title: none, full: true, style: "apa")
