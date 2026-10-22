// The literature: where ploidyspec stands -- included by ploidyspec_guide.typ
#import "guide_lib.typ": *
#show figure.where(kind: table): set block(breakable: true)

This chapter reads ploidyspec against the 86 works in Appendix C (reviewed
12 October 2026). It answers four questions: which methods already address
parts of our problem, and what each needs (14.1); where ploidyspec is new and
where it is not (14.2); whether our conclusions hold up against the published
record (14.3, 14.4); and what the literature changed in the tool and still
leaves open (14.5, 14.6). Section 14.7 indexes every work by theme.

#note(title: "The short version")[
  The literature supports our framing (a continuum with two axes, read per
  chromosome) and several of our results quantitatively: the snow carp fusion
  waves, the _Euphrasia_ subgenome divergence, the wheat subgenome partition.
  Our genuinely new contribution is reading both axes, per chromosome and
  along chromosomes, from one individual's phased assembly with no annotation,
  reference, reads or population. Subgenome partitioning by repeats is *not*
  new; we converge with TE-signature phasing. Three of our readings needed
  correcting: salmonids were missing delayed (LORe) rediploidization, masu
  rested on our own scaffolding, and "allo-like" for rose was worded too
  strongly. _Ranunculus_'s "no structure" was a fault of ours: its two
  subgenomes and eight homeolog pairs were hidden by the size adjustment and
  by subgenome-specific repeats, and are now found.
]

== The methods landscape

#figure(
  table(
    columns: (1.25fr, 1.35fr, 1.6fr, 1.6fr),
    table.header([*Approach*], [*Needs*], [*Reads*], [*Limit for our question*]),
    [Read k-mer spectra: GenomeScope 2.0, Smudgeplot @Ranallo-Benavidez2020-uv; Tetmer @Becher2020-go; review @Jenike2025-uu],
    [short reads],
    [genome-wide ploidy, heterozygosity, subgenome divergence @Khalaf2024-xo],
    [one number per genome; no chromosome or regional resolution; auto-like and young allo can look alike],

    [Ks distributions and WGD dating @Ren2018-pj @Qiao2019-wp @Dunn2025-qy],
    [annotated genes; often transcriptomes],
    [timing of duplications from duplicate-gene divergence],
    [Ks accumulates only between non-recombining copies @Gaynor2026-rc; the peak can reflect parental divergence, migration or coalescence @Dunn2025-qy],

    [Fractionation balance: P-index @Wang2019-zs],
    [annotation, synteny, an outgroup genome],
    [auto vs allo from balanced vs biased gene loss (threshold ~0.3)],
    [biased loss also arises after autopolyploid rediploidization @Parey2022-ic @Xie2026-aq],

    [Gene trees and ancestral karyotypes @Robertson2017-ya @Parey2022-ic @Redmond2023-fm],
    [many annotated genomes of related species],
    [when each duplicated region stopped recombining (AORe vs LORe); auto vs allo origin],
    [needs a clade of genomes; not per individual],

    [Duplicated BUSCOs, Hox clusters @McHale2025-wt @Farhat2023-fp],
    [annotation],
    [old WGDs beyond k-mer reach],
    [presence of a WGD only; nothing on inheritance],

    [Subgenome phasing by repeats or k-mers: PolyCRACKER @Gordon2019-in, TE signatures @Session2023-vh, SubPhaser @Jia2022-bo; benchmark @Zhang2023-va; review @Session2024-sl],
    [an assembly of a known allopolyploid; homeolog groups assigned first],
    [which chromosomes belong to which progenitor],
    [assumes allopolyploidy; fails on autopolyploids by design @Jia2022-bo; says nothing about rediploidization],

    [Differential divergence against relatives: M3 @Lv2025-te],
    [the polyploid plus related diploid genomes],
    [subgenome assignment, fast],
    [needs references from relatives],

    [Inheritance from crosses, families, populations @Bourke2017-hd @Bicknell2023-qk @Delaval2024-dr @Roux2015-wy @Dufresne2014-rt @Dalle-Palle2022-vr],
    [mapping populations, pedigrees or population samples; markers],
    [polysomic vs disomic inheritance, per chromosome or locus; origin by model choice],
    [the direct answer, but expensive and rare for wild species],

    [Resequencing, allele composition, meiotic cytology @Xie2026-aq @Gaynor2026-rc @Grandont2013-or],
    [a phased assembly plus population data or cytology],
    [tetrasomic vs disomic regions, mechanisms],
    [population data needed],

    [*ploidyspec*],
    [*one individual's haplotype-resolved, chromosome-scale assembly*],
    [*copies, duplicated sets, per-chromosome and per-window lineage structure, fusions, residual tetrasomy, exchange, switches, LORe-like blocks, karyotype change*],
    [*structure, not inheritance or origin; blind beyond ~8–10% divergence; the assembly is part of the measurement*],
  ),
  caption: [Methods that address parts of the auto/allo and rediploidization question.],
)

== Where ploidyspec is new, and where it is not

*New, as far as the reviewed literature shows:*

+ *Both axes from one individual's assembly.* Every other route to
  rediploidization (gene trees across a clade, resequenced populations,
  crosses) needs data beyond one genome. ploidyspec reads origin-like structure
  and rediploidization per chromosome and along chromosomes from the phased
  assembly alone, with no annotation, reference, outgroup or reads. That is
  what the stream of phased reference genomes from the Darwin Tree of Life
  and Earth BioGenome projects provides @Howe2026-nw @Exposito-Alonso2020-qh.
+ *Fusion-driven rediploidization recovered blind.* #tc(<Xie2026-aq>) needed
  haplotype-resolved genomes, resequencing, meiotic cytology and gene trees to
  show that fusions trigger rediploidization and that divergence spreads from
  the fusion site. ploidyspec finds every published snow carp fusion, places
  each junction at 94–100% agreement, and orders the waves by junction
  gradient (oldest 2.7–4.1, middle 2.1–2.6, youngest 0.6–0.7), from the
  assemblies alone (Section 10.3).
+ *Regional readouts with artefact guards.* Residual tetrasomy, homoeologous
  exchange, lineage switches and LORe-like blocks are each read per window and
  per pair, with guards against the failure modes the literature documents:
  centromeric satellites @Su2019-kj @Chen2024-qk, shared repeats, swapped
  contigs from reference-guided scaffolding @Alonge2019-fr, chromosome size,
  and numbering errors between haplotype files.
+ *A homeology map without genes.* Genome-wide segmental homeology from k-mer
  windows, which shows karyotype change since the duplication (chromosomes
  with blocks on two partners) and regions of different age.
+ *Robust to real assemblies.* Odd ploidy and a second chromosome set inside
  one haplotype file (AAB), renumbering cycles, and size artefacts are
  handled automatically; published methods mostly assume clean input.
+ *One scale across the tree of life.* Plants, fish, snails and a weevil on
  the same distances and the same rules.

*Not new, and we should say so:*

- Partitioning subgenomes by repeat history. Our durum wheat partition and the
  TE-marker split use the principle of TE-signature phasing @Gordon2019-in
  @Session2023-vh @Jia2022-bo; they are a by-product, not a claim.
- Ploidy from k-mers @Ranallo-Benavidez2020-uv @Becher2020-go; we read copies
  from the assembly and use read spectra only as a check.
- The concepts themselves: the continuum @Twyford2025-xy, divergence as a
  record of isolation @Gaynor2026-rc, LORe @Robertson2017-ya, fusion-triggered
  rediploidization @Xie2026-aq.

== Our conclusions checked against the literature

#figure(
  table(
    columns: (1.5fr, 1.7fr, 0.75fr, 1.6fr),
    table.header([*Our conclusion*], [*Literature*], [*Verdict*], [*Consequence*]),
    [Snow carp fusions in three waves; divergence peaks at the junction],
    [#tc(<Xie2026-aq>): unbalanced fusions trigger rediploidization, spreading from fusion sites in three waves],
    [supported],
    [an independent replication from assemblies alone],

    [`SchCurv1` chr17 splits only towards its end: an inversion in HAP3],
    [Xie et al.: a centromeric inversion made chr17's short arm disomic],
    [supported],
    [coordinates still to compare (Section 12.5)],

    [`daEupConf1`: duplicated sets ~0.057 apart (one haplotype)],
    [Tetmer: tetraploid _Euphrasia_ is allopolyploid, subgenome divergence ~5% @Becher2020-go],
    [supported],
    [a quantitative match between read- and assembly-based estimates],

    [Durum wheat: exact A/B partition and all seven A–B pairs, blind to the names],
    [subgenome-specific repeats and satellites @Su2019-kj @Gordon2019-in @Jia2022-bo @Session2023-vh],
    [supported],
    [present the partition as convergent with TE-signature phasing; the pairs show only against the cross-subgenome background, because each subgenome's own repeats make its chromosomes closer to each other than to their homeologs],

    [Glechoma residual block is a centromeric satellite],
    [centromeric satellites turn over and homogenize between and within subgenomes @Su2019-kj @Chen2024-qk @Yan2025-bz],
    [supported],
    [keep the complexity guard],

    [Allo-like ≠ allopolyploid; a rediploidized autopolyploid looks the same],
    [#tc(<Twyford2025-xy>), #tc(<Lv2024-yo>), #tc(<Doyle2017-ii>), #tc(<Parey2022-ic>)],
    [supported],
    [wording already says so],

    [Divergence is read as isolation of copies],
    [#tc(<Gaynor2026-rc>); #tc(<Dunn2025-qy>)],
    [supported],
    [--],

    [Mixed per-chromosome states are informative],
    [mixed and partial inheritance in rose, sturgeon, banana, _Brachiaria_ @Bourke2017-hd @Dalle-Palle2022-vr @Jeridi2012-vy @Worthington2016-jx; theory @Le-Comber2010-au],
    [supported],
    [present mixtures as expected, not anomalous],

    [Gar and bowfin: no duplication since 2R],
    [the outgroups to the teleost WGD @Parey2022-ic],
    [supported],
    [negative controls],

    [Pilosella: residual tetrasomy on 4 pairs (79 Mb), aligned at allelic identity],
    [_P. piloselloides_ is autopolyploid with tetrasomic inheritance @Bicknell2023-qk],
    [consistent],
    [an auto-like reading with residual tetrasomy fits the genus],

    [Weevil HAP1 chr09: A to 40 Mb, B after (a lineage switch)],
    [frequent exchange between A and B in AAB bananas @Xie2024-lk; gene conversion and sporadic recombination in a clonal amphitriploid fish @Wang2022-am; clonal polyploidy common in weevils @Stenberg2003-rt],
    [plausible],
    [if _Strophosoma_ is clonal, mitotic recombination or gene conversion, not meiotic exchange; or a phasing switch (Hi-C needed)],

    [Salmonids: "no residual tetrasomy" read as fully rediploidized],
    [~25% of each salmonid genome rediploidized late (LORe), >90% similar, telomere-including arms @Robertson2017-ya @Lien2016-lf; deviant SNPs concentrate there @Dallaire2023-bn],
    [incomplete],
    [LORe regions are not near-allelic, so the residual test misses them by design: LORe-like block class added (Section 6.6)],

    [Masu: residual tetrasomy on one pair (chr17/chr29)],
    [reference-guided scaffolding can misplace contigs @Alonge2019-fr; a phased Hi-C assembly exists @Wu2026-zm],
    [assembly-dependent],
    [re-run on Wu et al.'s assembly (Section 9.4)],

    [Rosa spinosissima: allo-like (high), 7 pairs],
    [rose pairing is partially preferential, varies along arms and between genotypes @Bourke2017-hd; segmental allopolyploidy @Stebbins1947-ei @Mason2020-mz],
    [overstated],
    [read "copies do not mix today on any of the 7 pairs"; segmental allopolyploidy is as likely as strict allopolyploidy],

    [Cerastium alpinum: 18 large–small homeolog pairs],
    [2n = 72 @Boscaiu1999-lh, octoploid throughout its range on x = 9 @Brysting2000-cm; an allopolyploid can look diploid-like in k-mers @Nieuwenhuis2024-hf],
    [partial],
    [the pairs are one level of an 8x genome; the deeper level (sets of four) is beyond k-mer reach],

    [Ranunculus repens: two subgenomes of 8 and eight homeolog pairs (first reported as no structure)],
    [diploid relatives have x = 8 @Karbstein2025-ou, so 2n = 32 is 4x; subgenome-specific repeats are the TE-phasing signal @Session2023-vh],
    [supported (after a fix)],
    [the earlier "no structure" was our size adjustment removing a larger and a smaller subgenome; read as allo-like with asymmetric subgenomes (A 1.23× B)],

    [Eleocharis: auto-like origin, residual tetrasomy on 10 of 19 numbers],
    [holocentric chromosomes: frequent fission and fusion, restricted chiasmata, inverted meiosis @Hipp2009-hz @Bures2013-zm @Mandrioli2020-hc],
    [caution],
    [chromosome numbers need not track ancestral units; treat arm- and number-based readings as provisional],

    [Hesperis: AAAB-like on 6 numbers],
    [Brassicaceae chromosome numbers change by descending dysploidy after WGD; number is no proxy for history @Mandakova2016-ul],
    [consistent],
    [do not infer x from the count],

    [Residual tetrasomy ≥ 90% terminal is low confidence],
    [real delayed regions include telomeres @Lien2016-lf; exchanges are subtelomeric @Zhang2020-ms],
    [partly contradicted],
    [keep the flag (repeats are also terminal), but terminal signal is expected biology, not presumed artefact],
  ),
  caption: [Each conclusion of the panel that the literature bears on.],
)

== The panel in the light of the literature

#table(
  columns: (1.1fr, 2.2fr, 2.2fr),
  table.header([*Species*], [*What the literature says*], [*What ploidyspec reads*]),
  [`SchCurv1`, `SchYoun1` snow carps],
  [single ancestral autotetraploidy; rediploidization started by fusions, incomplete at speciation @Xie2026-aq @Xie2026-zb],
  [4 × 25; mixed states; every fusion and the wave order; chr17 inversion],

  [`fSalAlp3`, `fCorLav1`, `OncMaso1` salmonids],
  [Ss4R autotetraploidy ~80–95 Mya; LORe over ~25% @Lien2016-lf @Robertson2017-ya; charr and whitefish deviant SNPs in recently rediploidized regions @Dallaire2023-bn; phased masu @Wu2026-zm],
  [map reaches 16–30% of the genome; karyotype reshuffled; LORe-like blocks: whitefish 243 Mb (40% of the block span, five pairs), charr 8.5 Mb (one pair, at chromosome ends)],

  [`fLepOcu1`, `fAmiCal2` gar, bowfin],
  [non-duplicated outgroups of the teleost WGD @Parey2022-ic],
  [one pair each, 0% in blocks: not assessable, as expected],

  [`lpTriTurg1` durum wheat],
  [AABB; subgenome-specific centromeric satellites @Su2019-kj; TE-based phasing separates wheat subgenomes @Gordon2019-in @Jia2022-bo; HE in nascent AADD @Zhang2020-ms],
  [partition exactly A vs B (z = 179) and all seven A–B pairs; none in the window map],

  [`icStrMela3` weevil],
  [AAB precedent with frequent exchange @Xie2024-lk; clonal polyploids keep gene conversion @Wang2022-am; weevil clonality and polyploidy @Stenberg2003-rt],
  [3 × 10, AAB-like (high); chr09 lineage switch],

  [`daEupConf1` _Euphrasia_],
  [tetraploid _Euphrasia_ allopolyploid, ~5% subgenome divergence @Becher2020-go],
  [one usable haplotype; 11 groups of 2 numbers, ~0.057 apart],

  [`daLatClan1`, `daLatSqua1` _Lathraea_],
  [parasitic Orobanchaceae: genome inflation and gene loss @Lyko2021-se; parasitic dodders are allopolyploid with homoeologous recombination @Segawa2025-tu],
  [allo-like with rearranged karyotypes; _L. clandestina_'s near-allelic block is repeats; one short terminal LORe-like pair (chr05/chr07, 5 Mb), possibly subtelomeric repeats],

  [`daPilAura1` _Pilosella_],
  [tetrasomic autopolyploidy in _P. piloselloides_; apomixis loci @Bicknell2023-qk],
  [residual tetrasomy on 4 pairs; one copy apart on most numbers; LORe-like blocks on two pairs (115 Mb)],

  [`dcCerAlpi1` _Cerastium_],
  [2n = 72, octoploid on x = 9 @Boscaiu1999-lh @Brysting2000-cm],
  [18 large–small pairs (x ≈ 18): one level of the duplication],

  [`dmRanRepe1` _Ranunculus_],
  [_Ranunculus_: polyploidy, apomixis and reticulation; diploid x = 8, 2.69 Gb @Karbstein2025-ou],
  [2 × 16: two subgenomes of 8 (A 1.23× B) and eight homeolog pairs (Section 14.3)],

  [`drRosSpin1` _Rosa_],
  [partial preferential pairing in tetraploid rose @Bourke2017-hd],
  [7 pairs, synchronous (CV 0.03), no residual tetrasomy],

  [`lpElePalu1` _Eleocharis_],
  [holocentric Cyperaceae; fission and fusion dominate number change @Hipp2009-hz @Bures2013-zm],
  [residual tetrasomy on 10 of 19 numbers, 310 Mb; 13 near-allelic map blocks, no LORe-like ones],

  [`ddHesMatr1` _Hesperis_],
  [Brassicaceae dysploidy and pseudo-ploidy @Mandakova2016-ul],
  [4 × 6, AAAB-like; one lineage switch],

  [`xgCepNemo3` and the land snails],
  [stylommatophoran WGD shown by duplicated BUSCOs and Hox clusters @McHale2025-wt; WGD in neogastropods @Farhat2023-fp],
  [little or no k-mer-detectable homeology: the duplication is beyond k-mer reach],
)

== What the literature changed in ploidyspec

#table(
  columns: (1.4fr, 2.4fr),
  table.header([*Change*], [*Why*]),
  [LORe-like block class in the homeology map, with a near-allelic lower bound and a two-way rule (Section 6.6)],
  [salmonid delayed rediploidization @Robertson2017-ya @Lien2016-lf was visible in the map but unreported],

  [Masu re-run on the published phased assembly],
  [@Wu2026-zm; our reference-scaffolded HAP2 produced artefacts (Section 9.4)],

  [Wording: positions on two axes, "copies do not mix today", never a species label],
  [@Twyford2025-xy @Lv2024-yo @Doyle2017-ii @Gaynor2026-rc],

  [Limits: holocentric genomes, partial inheritance, repeat-based signals (Section 12.3)],
  [@Hipp2009-hz @Bourke2017-hd @Session2023-vh],

  [Terminal residual signal flagged, not presumed artefact],
  [@Lien2016-lf @Zhang2020-ms],

  [Reading of "no pairs" in a known polyploid (Section 12.2)],
  [@Gaynor2026-rc],
)

Not yet in the code: a holocentric flag in the summary, the segmental
alternative in the "allo-like" sentence, and a confidence rule for terminal
residual tetrasomy that accounts for real terminal biology.

== Open questions the literature raises

- *Validate the LORe-like class.* Map charr and whitefish LORe-like blocks onto
  Atlantic salmon's seven delayed arm pairs @Lien2016-lf through synteny.
- *Whitefish versus charr.* Whitefish shows far more LORe-like sequence than
  charr; deviant SNPs in both @Dallaire2023-bn could test whether this is
  biology or assembly.
- *Masu on a clean assembly.* Does residual tetrasomy on chr17/chr29 survive
  on @Wu2026-zm?
- *Weevil reproduction.* Whether _Strophosoma melanogrammum_ is clonal decides
  whether its chr09 switch can be meiotic exchange @Stenberg2003-rt
  @Wang2022-am.
- *Rose ploidy and pairing.* Our read spectrum says tetraploid; a 2026 genome
  report says diploid; pairing in tetraploid rose is genotype-dependent
  @Bourke2017-hd.
- *The deeper Cerastium level* (8x on x = 9) and a primary count for this
  _Ranunculus repens_ individual @Brysting2000-cm @Karbstein2025-ou.
- *Head-to-head comparisons.* SubPhaser and M3 on wheat and Glechoma
  @Jia2022-bo @Lv2025-te; GenomeScope, Smudgeplot and Tetmer on the panel's
  reads @Ranallo-Benavidez2020-uv @Becher2020-go; Ks against our distances on
  the snow carps @Gaynor2026-rc @Dunn2025-qy.
- *Beyond k-mer reach.* Combine with duplicated-BUSCO or synteny methods for
  the snails and the deepest salmonid homeology @McHale2025-wt @Farhat2023-fp.

== The literature by theme

#table(
  columns: (1.3fr, 3fr),
  table.header([*Theme*], [*Works*]),
  [Definitions and the continuum], [@Stebbins1947-ei @Doyle2017-ii @Soltis2010-jk @Twyford2025-xy @Lv2024-yo @Leal2024-va @Thomas2026-kv],
  [WGD prevalence and consequences], [@Van-de-Peer2009-hf @Ren2018-pj @Qiao2019-wp @Crombez2025-iz @Yu2023-rd @Li2021-dp],
  [Meiosis and inheritance], [@Grandont2013-or @Soares2021-rh @Le-Comber2010-au @Bourke2017-hd @Dalle-Palle2022-vr @Jeridi2012-vy @Worthington2016-jx @Dufresne2014-rt @Roux2015-wy @Delaval2024-dr @Bicknell2023-qk],
  [Rediploidization and LORe], [@Robertson2017-ya @Lien2016-lf @Dallaire2023-bn @Parey2022-ic @Redmond2023-fm @Xie2026-aq @Xie2026-zb @Mandakova2016-ul @Feng2024-ew],
  [Homoeologous exchange and segmental allopolyploidy], [@Gaeta2010-va @Mason2020-mz @Deb2023-bo @Zhang2020-ms @Hurgobin2018-nw @Mu2023-qq @Leal-Bertioli2018-gp @Wohner2023-pu @Deng2026-gt @Burns2021-gz @Burns2024-md],
  [Repeats, TEs and centromeres], [@Su2019-kj @Chen2024-qk @Yan2025-bz @Shi2026-ar @Gordon2019-in @Session2023-vh],
  [Subgenome phasing], [@Jia2022-bo @Zhang2023-va @Session2024-sl @Lv2025-te @Song2024-am @Robbins2026-qo],
  [Divergence and Ks], [@Gaynor2026-rc @Dunn2025-qy @Wang2019-zs],
  [k-mer methods], [@Ondov2016-jo @Ranallo-Benavidez2020-uv @Becher2020-go @Jenike2025-uu @Khalaf2024-xo],
  [Assembly of polyploids and standards], [@Kyriakidou2018-xv @Zhao2025-yc @Sun2022-tc @Yang2017-zi @Alonge2019-fr @Howe2026-nw @Exposito-Alonso2020-qh],
  [Phased polyploid genomes], [@Wang2022-am @Xie2024-lk @Nieuwenhuis2024-hf @Segawa2025-tu @Wu2026-zm],
  [Holocentric chromosomes], [@Hipp2009-hz @Bures2013-zm @Mandrioli2020-hc],
  [Taxa in the panel], [@Boscaiu1999-lh @Brysting2000-cm @Karbstein2025-ou @Lyko2021-se @Stenberg2003-rt @McHale2025-wt @Farhat2023-fp],
)
