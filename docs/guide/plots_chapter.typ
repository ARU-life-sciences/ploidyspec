// Reading the diagnostic plots -- included by ploidyspec_guide.typ
#import "guide_lib.typ": *

Every stage writes plots next to its tables, and `report.html` collects them.
This chapter takes them one at a time: what is drawn, how to read it, the
patterns that matter, and the traps. Examples are from real runs.

== Conventions shared by every heatmap

- *Colour runs from yellow (close) to dark purple (far)* (matplotlib's
  `viridis_r`). Yellow always means "similar sequence".
- *The scale is stretched* to the 2nd–98th percentile of the values in that
  plot (stated in the title), so colours are comparable within a plot, not
  between plots. Read the colour bar.
- *Light grey means no data*: a copy shorter than the others ends earlier, or
  a cell is masked on purpose (the diagonal of a contrast heatmap).
- *Distances are k-mer distances*, not alignment identity: roughly the
  per-base divergence where sequence is shared, but containment-based window
  distances saturate in repeat-dense sequence (Section 6.10).

== Whole-chromosome distance heatmap (`matrix/whole_chrom_distance_heatmap.png`)

*What it shows.* Every chromosome copy against every other, at the chosen k
(Section 5.2), clustered so that similar copies sit together.

#figure(image("../../results/daGleHede1/matrix/whole_chrom_distance_heatmap.png", width: 78%),
  caption: [Glechoma (two haplotypes × 18 numbers). Yellow 2 × 2 squares on the
  diagonal: the two alleles of each chromosome (~0.01). Teal 2 × 2 squares next
  to them: a homeolog pair, e.g. chr10 with chr04 (~0.05). Everything else is
  purple: unrelated chromosomes (~0.11–0.12).])

*How to read it.* Look for the block sizes along the diagonal:
- *Yellow squares the size of the copy number* (2 × 2 with two haplotypes,
  4 × 4 with four): the copies of each chromosome. In a tetrasomic genome the
  4 × 4 block is uniformly yellow; in a genome with two lineages per chromosome
  it splits into two yellow 2 × 2 squares inside a lighter 4 × 4.
- *Teal blocks pairing two numbers*: homeologs. If most numbers have one, the
  genome has duplicated sets (Glechoma: nine).
- *A copy whose row is lighter or darker than its siblings'*: one divergent
  copy (AAAB-like, or an odd haplotype).

*Traps.* A yellow cell far from the diagonal between two different numbers
means a numbering mismatch: check `chrom_label_corrections.tsv` (the matrix
stage fixes most). Chromosome size also shifts the background colour slightly
(Section 6.4), which is why the homeolog test uses adjusted distances.

== Contrast heatmap (`whole_chrom_distance_heatmap_contrast.png`)

The same matrix with the diagonal masked (white) and the colour stretched over
the remaining values. Masking the near-zero self cells lets the background use
the whole colour range, so faint structure among unrelated chromosomes, for
example two groups of chromosomes that share repeat families (wheat's A and B,
Section 10.2), becomes visible. Read it for *structure in the background*, not
for the copies themselves.

#figure(image("../../results/daGleHede1/matrix/whole_chrom_distance_heatmap_contrast.png", width: 62%),
  caption: [Glechoma, contrast version: the homeolog squares stand out against a
  flat background; there is no further grouping of chromosome numbers.])

== k-resolution diagnostic (`matrix/k_resolution_diagnostic.png`)

*What it shows.* For every pair of copies (one grey line each), the corrected
distance at each k of the sweep; points are red where that k is below the
chance-collision floor (Section 5.2).

#figure(image("../../results/daGleHede1/matrix/k_resolution_diagnostic.png", width: 72%),
  caption: [Glechoma. Three bands: alleles (~0.01–0.02), homeologs (~0.05) and
  unrelated chromosomes (~0.11–0.13). At k = 11 everything collapses to ~0:
  short k-mers are shared by chance, so k = 11 cannot separate anything.])

*How to read it.* The distances should level off as k grows: the plateau
(here k ≥ 17) is the real divergence, and the matrix uses the largest usable
k. Distinct bands are the classes of comparison; a genome with homeologs has a
middle band. Red points at large k mean those pairs are beyond what k-mers can
measure.

== Homeolog candidates (`homeologs/homeolog_pairs.png`)

*What it shows.* Every pair of chromosome numbers, ranked by mean distance;
accepted homeolog pairs in red and labelled.

#figure(image("../../results/daGleHede1/homeologs/homeolog_pairs.png", width: 85%),
  caption: [Glechoma: nine accepted pairs (red, ~0.05) far below a continuous
  background (~0.11–0.12). A gap like this is the clearest picture of a
  duplicated genome.])

*How to read it.* A clean duplication shows a *gap* between a low cluster and
the rest. A gradual slope without a gap, with a few red points at its foot,
means weak evidence: read the z-scores in `homeolog_candidates_ranked.tsv`.
The plot shows raw distances; acceptance is tested on size-adjusted ones
(Section 6.4), so a red point can sit slightly above a blue one.

== Reordered chromosome-number heatmap (`homeolog_pairs_reordered_heatmap.png`)

*What it shows.* Chromosome numbers (each averaged over its copies) clustered
by average linkage, accepted pairs outlined in red.

#figure(image("../../results/daGleHede1/homeologs/homeolog_pairs_reordered_heatmap.png", width: 62%),
  caption: [Glechoma: nine yellow 2 × 2 squares, each outlined, on a flat
  background: two sets of nine.])

*How to read it.* The square sizes give the number of sets: 2 × 2 for sets
of two (tetraploid on a two-haplotype assembly), 3 × 3 for three (hexaploid
history, e.g. `drMyrSpic1`). Larger blocks of moderate colour without
outlines are diffuse group structure (the genome partition, Section 6.7).

== Homeology map (`homeologs/homeology_map.png`)

*What it shows.* Each chromosome of the mapped haplotype as a bar, painted
along its length by the chromosome its homeolog block lies on (labels give the
partner number). Grey: no homeolog detected. Hatched and faded: a block
without a block back from the partner (not counted).

#figure(image("../../results/daGleHede1/homeologs/homeology_map.png", width: 100%),
  caption: [Glechoma: almost every chromosome painted by one partner over most
  of its length (chr01 by chr05, chr02 by chr03, ...): 79% of the genome in
  homeolog blocks, every chromosome with one partner.])

*How to read it.* One colour per chromosome along most of its length: a
duplicated genome with intact chromosomes (Glechoma, Salix, Trifolium). Two or
more colours on one chromosome: its arms are duplicated on different
chromosomes, i.e. fusions, fissions or translocations since the duplication
(salmonids, Potamogeton). Mostly grey: no detectable duplication, or one too
old for k-mers (Section 3.3). The charr map is in Section 9.4.

== Windowed tracks along a chromosome (`windowed/windowed_chrNN.png`, `_heatmap.png`)

*What it shows.* Along each copy of one chromosome number, in 250 kb windows,
the distance to each other copy (position-free, Section 5.4). The line plot
draws one line per directed pair; the heatmap one row per directed pair
("HAP3 vs HAP1" = windows of HAP3 against the whole of HAP1).

#figure(image("../../results/SchCurv1/windowed/windowed_chr19_heatmap.png", width: 100%),
  caption: [Snow carp `SchCurv1` chr19. HAP1 vs HAP2 and HAP3 vs HAP4 are
  yellow all along: two pairs of close copies. Every cross row (HAP1/HAP2
  against HAP3/HAP4) is dark: two lineages over the whole chromosome. HAP3 and
  HAP4 run to 69 Mb because they carry the fused chr22; the HAP1 and HAP2
  rows end at 38–40 Mb (grey).])

#figure(image("../../results/SchCurv1/windowed/windowed_chr19.png", width: 100%),
  caption: [The same chromosome as lines: the cross pairs (0.04–0.06 on the
  chr19 part) sit far above the within-lineage pairs (near 0).])

*How to read it.*
- *All rows similar and low*: interchangeable copies (tetrasomic-like).
- *Rows in two groups*, low within and high across, all along: two lineages
  (`resolved_lineages`, `fusion_lineages`).
- *Groups that appear only over part of the chromosome*: a regional split
  (snow carp chr17, the last 8 Mb: an inversion in HAP3).
- *One copy's rows dark everywhere*: one divergent copy.
- *A row that is low then high, while another goes high then low at the same
  point*: a copy changing lineage (a lineage switch, Section 6.9).

*Traps.* Isolated dark spikes in one row are usually sequence present in one
copy and missing from the other (gaps, unplaced contigs), not divergence.
Short scattered blocks in a polysomic genome are ordinary haplotype structure;
the tool only counts splits that beat a null (Section 6.9).

== Genome overview (`windowed/windowed_genome_overview*.png`)

Every chromosome's heatmap side by side on one shared colour scale, x axis in
Mb. Use it to find the chromosomes worth opening.

#figure(image("../../results/SchCurv1/windowed/windowed_genome_overview_heatmap.png", width: 82%),
  caption: [`SchCurv1`: most chromosomes are speckled yellow-green (four
  interchangeable copies with scattered haplotype blocks). Three stand out:
  chr19 (two lineages along its whole length, the fusion), chr11 (banded rows
  over 5–28 Mb) and chr17 (dark rows at its end, the inverted region).])

*How to read it.* Look for *banding*: rows that stay dark or light over long
stretches. Speckle is noise or ordinary polymorphism; bands are structure.

== Tracks between homeologs (`homeologs/windowed_chrAAxBB*.png`, overview)

The same plots for the copies of a homeolog pair: windows of each copy of one
number against each copy of the other. The level is the homeolog distance; a
stretch where it falls to the allelic level is a residual-tetrasomy or
exchange candidate (Section 6.10).

#figure(image("../../results/daGleHede1/homeologs/windowed_chr14x17.png", width: 100%),
  caption: [Glechoma chr14 × chr17: homeologs ~0.03–0.05 apart, falling to
  ~0.003 over 14–23 Mb in every direction. This is the trap: the block is a
  centromeric satellite shared by the two homeologs (it compresses 3–5×
  better than the rest and aligns at only 94–96%), not retained pairing.])

*How to read it.* A drop seen from *all* copies of both numbers, in sequence
that is not repeat-dense, is residual tetrasomy (Pilosella chr04 × chr05). A
drop in one copy only, while its own homolog stays far, is an exchange. A
drop exactly where the chromosome is repeat-rich, typically mid-chromosome or
at the ends, is most likely shared repeats: the tool now excludes those
windows, but the plot still shows them.

== TE-marker painting (`subgenomes/te_markers_windowed_*.png`)

*What it shows* (with `--with-te-markers-windowed`). For a pair of copies, each
window coloured by whose fossil-TE markers it matches (Section 6.8): blue for
the first copy, red for the second, grey when markers of both or neither
dominate (ambiguous), light grey when the window has no markers.

#figure(image("../../results/daGleHede1/subgenomes/te_markers_windowed_HAP1_chr01xHAP2_chr01.png", width: 100%),
  caption: [Glechoma chr01, its two alleles: almost all grey (ambiguous), with
  a few short blue and red patches. The two alleles share one repeat history,
  so this is what "no subgenome difference" looks like.])

*How to read it.* Between copies from different subgenomes or lineages, each
copy should be painted mostly in its own colour along long stretches: separate
repeat histories. Between alleles of one lineage (as above) the painting is
grey. Mixed painting with long switches between colours along one copy points
to exchange between lineages. Read it together with `te_split`: the painting
shows where, the split ratio how strong.

== Fusion profile (`rediploidization/fusion_profile_<scaffold>.png`)

*What it shows.* Along a fused chromosome, each window's distance to the
unfused copies of the component it belongs to (coloured by component), and to
the other fused copy (grey). The title gives the junction and the gradient
(Section 6.9).

#figure(image("../../results/SchYoun1/rediploidization/fusion_profile_Sy_Chr19_22_P2.png", width: 100%),
  caption: [`SchYoun1`, the chr19+chr22 fused scaffold of HAP4. The chr22 part
  (orange) and the chr19 part (blue) meet at 33 Mb. Divergence from the
  unfused copies is ~0.04–0.05 around the junction and falls to ~0.01 at both
  ends (gradient 3.7); the two fused copies (grey) stay ~0.003 apart
  throughout.])

*How to read it.* A peak at the junction falling towards the ends: recombination
with the unfused copies stopped first next to the fusion and divergence is
spreading outwards (Xie et al. 2026). A flat profile: a young fusion (the
youngest snow carp wave, gradient 0.6–0.7) or an old one whose divergence has
already spread. The grey points show the fused lineage's own diversity.

== Panel PCA (`meta/supplementary/poly_space_pca.png`)

*What it shows.* An exploratory principal component analysis of six numbers
per species (allelic distance, homeology-map duplicated fraction, paired
fraction, best partition z, resolved and tetrasomic fractions), each scaled to
unit variance; species missing any number are left out.

#figure(image("../../meta/supplementary/poly_space_pca.png", width: 88%),
  caption: [52 species. Left: duplicated, paired, allo-like genomes (Salix,
  Sorbus, Glechoma, Potamogeton). Bottom right: assemblies with four copies
  where many chromosomes are tetrasomic-like (Buddleja, Galium, Lythrum). Top
  right: genomes without detectable duplication (snails, gar, bowfin, _Malva_).])

*How to read it.* PC1 (50%) runs from duplicated genomes (negative:
duplicated fraction, paired fraction and partition z all load −0.48 to −0.52)
to undivided or tetrasomic ones (positive: tetrasomic fraction +0.25,
allelic distance +0.21). PC2 (27%) separates genomes whose chromosomes split
into lineages (positive: resolved fraction +0.44) from high-heterozygosity,
tetrasomic ones (negative: allelic distance −0.56, tetrasomic fraction −0.50).
It is a picture of where species sit relative to each other, not a classifier.
