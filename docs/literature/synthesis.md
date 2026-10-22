# ploidyspec and the literature: synthesis (67 papers, Oct 2026)

## 1. How the literature shapes what we are doing

**Frame: the polyploid continuum (Twyford 2025).** Twyford et al. argue that
"auto vs allo" is a continuum. It runs from autopolyploidy, through segmental
allopolyploidy (Stebbins 1947), to strict allopolyploidy, and genomes move along
it over time as they diploidize. ploidyspec measures where one individual sits on
that continuum, along two axes:

1. **Origin-like structure:** how distinct the chromosome sets are (partition z,
   homeolog pairs, duplicated fraction).
2. **Rediploidization:** how far the copies have gone from recombining to
   diverging (residual tetrasomy, exchange, lineage switches, fusions).

The tool should say this explicitly. It reports positions on axes, not a binary
label. Lv 2024 shows why the binary label is unsafe:

- most autopolyploids (39 of 43, 91%) are polysomic;
- several taxonomic *allo*polyploids are polysomic too;
- Bourke 2017 shows pairing behaviour in tetraploid rose varies by chromosome
  and by genotype.

**What k-mer divergence can and cannot see (Gaynor 2026).** Ks, and our Mash
distance likewise, only accumulates between copies that do *not* recombine. A
young autopolyploid with polysomic pairing is therefore invisible to divergence
by construction. "No structure" means "copies still recombine, or have only just
stopped". It does not mean "auto origin". This is the same logic Gaynor uses to
argue that all *detectable* ancient WGDs involve hybridization. Our "auto-like"
readings are statements about current recombination, not about history. The
guide should use that wording.

**Rediploidization is asynchronous and regional (Lien 2016; Robertson 2017;
Parey 2022; Redmond 2023; Dallaire 2023).**

- In salmonids about 25% of the genome (7–8 arm pairs, >90% similarity) is in
  delayed rediploidization, i.e. LORe.
- Teleosts and sturgeon/paddlefish show the same pattern on deeper timescales.
- Dallaire 2023 places ongoing residual tetrasomy in specific arms.

Per-chromosome and per-window readouts, which we already have, are the right
granularity.

**Exchange is distal and recombination-like (Zhang 2020 wheat; Deb 2023
review).** Homoeologous exchanges are biased to subtelomeric, low-copy regions.
Our exchange readout, and the "terminal-only" guard on residual tetrasomy
(RESIDUAL_TERMINAL_ONLY), have to be read with this in mind. Terminal signal
is expected biology, not an artefact.

**Centromeres turn over fastest (Chen 2024 Brachypodium; Yan 2025).**
Centromeric satellites diverge and are replaced quickly between subgenomes.
This supports our satellite and low-complexity guard: the Glechoma false
residual call was a centromeric satellite. It also predicts that
subgenome-specific TE or satellite content can make chromosomes of one
subgenome closer to each other than to their homeologs. That is what we see in
lpTriTurg1, dmRanRepe1 and daInuConz1, where the homeolog pairs show only
against the cross-subgenome background.

**Karyotype change drives rediploidization (Mandáková 2016; Xie 2026; Hipp
2009).** Descending dysploidy after WGD produces "pseudo"-ploidy:
Cardamine is a tetraploid that reads as a pseudotriploid. Fusions lock
recombination and create divergence gradients away from the junction.
Holocentric lineages (Carex, and our Eleocharis) fission and fuse freely,
which weakens arm-level and junction readings.

## 2. Where the novelty is (and isn't)

| Capability | Prior art | ploidyspec |
|---|---|---|
| Ploidy / heterozygosity from reads | GenomeScope2 / Smudgeplot (Ranallo-Benavidez 2020) | Not our claim; we read assemblies, and use reads only to check odd ploidy (Salix 28/60/87) |
| Subgenome phasing by TE signatures | Session & Rokhsar 2023; PolyCRACKER; SubPhaser (Jia 2022) | **Convergent, not novel.** Our wheat A/B partition reproduces this. Those methods need homeolog groups first and assume allopolyploidy; SubPhaser fails on autos |
| Mash-based subgenome partition | M3 (Lv 2025) | M3 needs outgroup/progenitor references; ours is reference-free and single-individual |
| WGD detection and dating | Ks distributions, P-index (Wang 2019), Van de Peer 2009 | Ks needs annotation and sees only non-recombining copies (Gaynor); we work on raw sequence |
| Delayed rediploidization maps | Gene-tree methods (Robertson 2017; Parey 2022; Redmond 2023); resequencing and linkage (Xie 2026; Dallaire 2023; Bourke 2017) | We get regional readings from **one individual's haplotype-resolved assembly**, with no annotation, population or cross |

**Defensible novelty:**

1. **Annotation-free, reference-free, single-individual reading of both axes**,
   per chromosome and per window. It works from the haplotype-resolved
   assemblies that Tree of Life and the Earth BioGenome Project now produce in
   bulk (Howe 2026 standards). No other tool takes "phased FASTA in, continuum
   position out".
2. **Fusion-driven rediploidization from the assembly alone.** Junction
   divergence gradients replicate Xie 2026's snow-carp waves:
   - oldest wave 2.7–4.1, middle 2.1–2.6, youngest 0.6–0.7;
   - the chr17 short-arm inversion is consistent too.

   Xie needed resequencing populations and gene trees to get this.
3. **Residual tetrasomy, exchange and lineage-switch readouts with explicit
   artefact guards** (satellite / low-complexity, control-shared, repeat-level
   identity, size). Each of the guards answers a failure documented in the
   literature.
4. **Practical robustness for real assemblies:** odd ploidy and extra sets
   (AAB), chromosome-numbering errors (swaps, displacement), and size
   artefacts. Methods papers rarely deal with these, and they matter at scale.
5. **Comparative panel** (more than 50 species, plants, fish and invertebrates)
   on one scale.

**Don't claim:** subgenome phasing per se, ploidy estimation from reads, or
dating WGDs.

## 3. Are our inferences right? Conclusion-by-conclusion

| Our conclusion | Literature verdict | Action |
|---|---|---|
| Snow carp fusion waves / gradients | **Supported** (Xie 2026) | Cite as validation |
| SchCurv1 chr17 HAP3 inversion | **Supported** (Xie's short-arm centromeric inversion) | Cite |
| Wheat: exact A/B partition and all 7 A–B pairs, blind (after the 21 Oct subgenome fix) | **Supported** (Session & Rokhsar; Chen 2024 centromere turnover) | Reword the partition as "convergent with TE-signature phasing" |
| Glechoma residual = centromeric satellite | **Supported** (satellite dynamics) | Keep guard |
| Weevil chr09 lineage switch | **Plausible.** Matches frequent homologous exchange in AAB banana (Xie 2024) | Cite as a biological precedent |
| Cerastium: 18 large–small pairs, size-asymmetric subgenomes | **Plausible.** 2n=72 confirmed (Boşcaiu 1999); okra shows an allopolyploid can look diploid-like in k-mers (Nieuwenhuis 2024) | Keep; call it "asymmetric subgenomes" without assigning progenitors |
| Rosa spinosissima "AABB / allo-like" | **Overstated.** Bourke 2017: rose pairing is partially preferential and varies by chromosome and genotype | Reword: "segmental / partially preferential". Partition strength ≠ disomy |
| Ranunculus: no structure once size removed | **Wrong, fixed 21 Oct.** 4x on x=8 (Karbstein 2025): two subgenomes of 8 (A 1.23× B) and 8 one-to-one homeolog pairs. The size adjustment and subgenome-specific repeats had hidden them | Read as allo-like with asymmetric subgenomes |
| "No residual tetrasomy" in salmonids ⇒ rediploidized | **Incomplete.** LORe regions are low-divergence but *not* near-allelic, so our residual test (< 0.25× homeolog distance) misses them by design | **New check (below).** Add a LORe-like class |
| Masu residual on chr17/chr29 only (8.6 Mb) | **Under-powered.** Our masu HAP2 was self-scaffolded (RaGOO-like; Alonge 2019); Wu 2026 published a Hi-C phased assembly (JBNQWP/JBNQWQ000000000) | Rerun masu on Wu 2026 |
| Eleocharis 10/19 chromosomes residual, 310 Mb | **Caution.** Holocentric (Hipp 2009; Mandrioli & Manicardi 2020); fission/fusion is free and arm logic does not apply | Flag holocentric lineages in the summary |
| Pilosella residual tetrasomy real | Consistent with apomictic / polysomic Hieracium-group biology (Lv 2024) | Keep |
| Mixed per-chromosome states (some pairs disomic, some not) | **Supported**: mixed inheritance in sturgeon (Redmond 2023) and rose (Bourke 2017); segmental allopolyploidy (Deng 2026; Leal-Bertioli 2018) | Present it as the expected outcome, not an anomaly |
| Potato-like autos with partial structure | Sun 2022: about 50% IBD between haplotypes, about 20% SV | IBD blocks can mimic "near-allelic"; we guard via own-homolog rule, so mention it |

### New check: LORe is already visible in our homeology maps

Blocks with median homeolog distance < 0.6× the genome median homeolog-block
distance:

| Species | Low blocks | Mb (of mapped) | Where |
|---|---|---|---|
| fCorLav1 (whitefish) | 19 | 266 of 683 (39%) | whole arms (chr04/37, chr03/13, chr40/06, chr35/38 at 0.002) |
| fSalAlp3 (charr) | 6 | 23 of 490 | chromosome starts (chr09/11, chr34/36) |
| OncMaso1 (masu) | 9 | 28 of 466 | mostly chromosome ends (chr30/31, chr32/33, chr03/15) |

These are LORe-like: arm-scale, reciprocal and low-divergence, but not
near-allelic. The whitefish share fits the ~25% salmonid expectation (Lien
2016) and overshoots it. Masu's are terminal and short, which may be the
self-scaffolded assembly losing them. Charr is between the two. The summary
currently says nothing about these blocks.

**Proposed:** a "delayed-rediploidization (LORe-like)" block class in the
homeology map summary. The rule is block distance below a fraction of the
genome median homeolog distance, but above the residual threshold. Validate it
against Lien 2016's arm list via the charr/salmon synteny.

## 4. Recommended changes before the paper

1. Reword "auto/allo" as positions on the continuum. Say "copies recombining"
   rather than "autopolyploid" (Twyford 2025; Gaynor 2026; Lv 2024).
2. Add the LORe-like block class (above), and report its Mb for salmonids.
3. Rerun masu on the Wu 2026 phased assembly.
4. Reword rose (and other "allo-like" calls with partial pairs) as
   segmental / partially preferential.
5. Add a holocentric caveat flag (Eleocharis, and Juncaceae/Cyperaceae in
   general).
6. Cite convergent prior art for the wheat partition (Session & Rokhsar,
   SubPhaser, M3) and position our partition as a by-product, not the claim.
7. Note that a terminal bias in exchange is expected (Zhang 2020).

## 5. Paperpile mismatches

- The file named "Becher 2020" contains **Deng 2026** (segmental allopolyploid
  with HE). The Becher/Tetmer paper itself is missing, and it is a key methods
  comparison.
- The file named "Bureš 2013" is **Mandrioli & Manicardi 2020** (holocentric).
