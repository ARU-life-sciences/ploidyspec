# daBudDavi1 — *Buddleja davidii* (butterfly bush)

## 1. Species overview

- **ToL ID**: `daBudDavi1`
- **Species**: *Buddleja davidii* Franch. (butterfly bush), family
  Scrophulariaceae
- **Assembly**: release stage, 2 source FASTA files (`daBudDavi1.hap1.2`,
  `daBudDavi1.hap2.2`); `AUTO` haplotype detection recovers **4 haplotype
  copies** per chromosome from header tags
- **Structure**: 19 chromosome numbers, 4 copies each (76 chromosome-scale
  units total) — matches the literature chromosome count exactly (see §2)

## 2. Literature context

**Ploidy level**: tetraploid, 2n=76, base number x=19 — matches this
assembly's structure exactly (19 chromosome numbers × 4 copies).

**Origin — allopolyploid, medium-high confidence, with a documented
complication.** Yang et al. 2023 (*Annals of Botany*) used PhyloNet
network inference across the Asian *Buddleja* clade and identified several
allopolyploid speciation events, with *B. davidii* clustering among them.
Critically, the same study found a **cytonuclear conflict specific to this
species**: *B. davidii* clusters with five other polyploid taxa in the
ASTRAL/nrDNA (nuclear) trees, but nests with *diploid* species in the
plastid (maternal) tree — a real, published signature of complex hybrid
ancestry, not a data-quality artifact.

Key citation: Yang, Ge, Guo, Olmstead & Sun 2023, *Annals of Botany*
(reticulate evolution of Asian *Buddleja*). Full detail in
`meta/literature_ploidy.tsv`.

## 3. Direct data outputs

| metric | value | note |
|---|---|---|
| Structural pairing | 0/19 chromosomes | no ancient-homeolog pair clears FDR significance |
| `te_marker_fraction` | mean 0.167, range 0.052–0.573 (n=114 pairs) | moderate-low, closer to the auto end of the panel's range |
| `partition_consistency` | 0.26 (rotating/low) | no fixed HAP-label-to-subgenome grouping recurs across the 19 chromosomes |
| Diffuse chromosome-number partition (`genome_partition.py`) | no significant structure at any k | z < 2 for every tested grouping — unlike `daInuConz1`/`dcCerAlpi1`/`dmRanRepe1`, there's no detectable "some chromosomes resemble each other more" signal here either |
| Hierarchically-reordered heatmap | smooth gradient, no blocks | distance ranges narrowly (0.100–0.106) with a continuous gradient from `chr19` (most divergent) to `chr07`/`chr11` (least), not discrete clusters |
| `distance_ratio_cv` / `pair_depth_cv` | n/a | require ≥2 accepted pairs; none exist |

## 4. Process inference

Unlike most other species in this panel, `daBudDavi1` doesn't fall neatly
into one of `INTERPRETATION.md`'s described conflict patterns — every
metric here reads weak-to-neutral rather than strongly contradictory.
There's no discrete pairwise signal, no diffuse block structure, and no
consistent within-chromosome subgenome grouping; the reordered heatmap
shows a smooth continuum rather than the sharp pairs or clean bipartitions
seen elsewhere in the panel. That absence of *any* strong structural
signal is itself consistent with, rather than contradicting, the
literature's own finding: a genuine nuclear-vs-plastid cytonuclear
conflict is exactly the kind of history (recent introgression, incomplete
lineage sorting, or a hybridization event whose nuclear and organellar
genomes tell different stories) that would leave the raw sequence-content
signal this pipeline measures muted and hard to resolve into a clean
auto/allo structural signature, even when the origin genuinely is
allopolyploid.

## 5. Auto/allo statement

**Allopolyploid, medium-high confidence — from the literature, not from
this pipeline's own data.** The published phylogenomic evidence (network
inference placing *B. davidii* among confirmed allopolyploid speciation
events in its clade) stands on its own. This pipeline's structural and
repeat-content metrics neither confirm nor contradict that call — they
read as uniformly weak/neutral, which is itself plausibly explained by the
same cytonuclear conflict Yang et al. (2023) already documented for this
exact species, rather than being an independent problem with the allo
call. **For a genome note**: *"Chromosome-level k-mer divergence analysis
found no discrete retained-duplicate homeolog pairs and no diffuse
subgenome partition, consistent with the cytonuclear conflict reported by
Yang et al. (2023), who found B. davidii nesting among confirmed
allopolyploids in nuclear gene trees but with diploid taxa in the plastid
tree."*
