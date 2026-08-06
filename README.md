# Ploidyspec

This will be a bash tool which runs a series of FastK commands (searches in path/provide binary) to:

1. Split each chromosome to its own FASTA (one per hap1 chromosome, one per hap2 copy) — samtools faidx or seqkit split on the existing assemblies.
2. FastK -k15 -t1 -Nchr1 chr1.fa per chromosome → produces a .ktab (the distinct k-mer set + counts for that sequence alone).
3. Logex (FastK's set-algebra tool) does boolean ops directly between .ktab files — union/intersect/diff — so you can get "k-mers in copy A ∩ copy B" as a table, then read off its size from the resulting .hist. That's your pairwise intersection count without writing custom code.
  - I'd sanity-check the exact Logex syntax against --help/the GitHub README when you get to it — I'm confident in the building blocks (FastK/Tabex/Profex/Logex/Histex) but not 100% on current flag syntax, it's changed across versions.
4. One caveat: use canonical k-mers (min of k-mer and its reverse complement) — chromosomes/copies could be assembled on opposite strands, and FastK does support canonical counting by default I believe, but worth confirming in --help.

For the windowed/dotplot version, neither tool does this out of the box — you'd want: bin each chromosome into windows (e.g. 50–100kb), and for every k-mer shared between copy A and copy B, look up which window it falls in on each side, and tally into a (window_A × window_B) co-occurrence matrix, then render as a heatmap. That's the same mechanism a classical dotplotter (Gepard/dotter) uses internally — a "dot" is just a shared k-mer, binned instead of point-plotted. FastK's .prof output (per-position profile, via -p) plus Profex gives you the position info needed to build this; alternatively, since these are single chromosomes (tens of Mb), I could just do the k-mer extraction + windowed intersection directly in Python/numpy (2-bit-pack the k-mers, canonicalize, sort+intersect) — reusing the exact np.intersect1d pattern already in analyze_substitution_sharing.py, no new tool dependency, and no cluster install needed.

## Implementation

The plan above is implemented as `ploidyspec.sh`, a thin wrapper around the `ploidyspec` Python
package (stdlib + numpy/matplotlib/scipy for matrix assembly and plotting). It shells out directly
to `samtools`, `FastK`, `Logex` and `Histex` (auto-detected from sibling `samtools*`/`FASTK*`
directories, or `--samtools`/`--fastk-dir`, or `$PATH`) — no k-mer math is reimplemented in Python.
Instead of the window×window dotplot co-occurrence matrix sketched above, the windowed stage
computes a 1-D **windowed distance track**: at each window it re-runs FastK on just that slice and
intersects it against the *same window coordinates* in every other haplotype copy of that
chromosome. That's a much cheaper way to answer "how auto/allopolyploid is this genome" — it
directly shows, along the chromosome, where haplotype copies are homogeneous (autopolyploid-like)
vs. sharply diverged (allopolyploid subgenomes), and where the split changes (candidate
homeologous exchange / introgression breakpoints) — at the cost of assuming rough colinearity
between the haplotype copies of a given chromosome (true here since these are chromosome-scale
phased assemblies of the same karyotypic chromosome; see the ≥20% length-difference warning the
tool prints if that assumption looks shaky for a given chromosome).

### Pipeline stages

1. **prepare** — parse a manifest of haplotype assembly FASTAs into per-chromosome "units". Each
   sequence is assigned a haplotype label (explicit per file, or `AUTO`-detected per sequence from
   its FASTA description) and a chromosome number (via regex against the description), filtered to
   sequences ≥ `--min-len` (default 1 Mb) so unplaced scaffolds/contigs are dropped. Output:
   `sequences.tsv` (+ `unplaced.tsv` listing everything excluded and why).
2. **kmers** — `samtools faidx` extracts each unit into its own FASTA, then `FastK -k<K> -t1`
   builds a canonical k-mer table (`.ktab`/`.hist`) per chromosome, in parallel.
3. **matrix** — all-vs-all whole-chromosome k-mer Jaccard distance
   (`1 - |A∩B|/|A∪B|`, via `Logex -H 'x=A &. B'` + `Histex -A`) across every unit, batched up to 8
   tables per `Logex` call (its letter-argument limit) to cut process overhead. Output: a full
   distance matrix CSV, a long-form pairs TSV (shared/union/Jaccard/containment), and a
   clustered heatmap PNG.
4. **windowed** — units are grouped by chromosome number (i.e. every haplotype's copy of "chr7"
   together); within each group, tumbling windows (default 250 kb) are re-extracted and
   re-`FastK`'d, and Jaccard distance is computed between every pair of haplotype copies at each
   window (again batched into one `Logex` call per window when ≤8 copies). Output: per-chromosome
   TSV + line plot of distance vs. position for every haplotype pair, plus a genome-wide small-
   multiples overview PNG.
5. **all** — runs the four stages above in order; each stage skips work it's already done, so `all`
   is safe to re-run after an interruption.

### Usage

```
./ploidyspec.sh all \
  --manifest manifests/daBudDavi1.tsv \
  --outdir results/daBudDavi1 \
  --k 15 --min-len 1000000 \
  --window 250000 \
  --threads 8
```

Or run stages individually (`prepare`, `kmers`, `matrix`, `windowed`) with the same flags.

The manifest is a TSV of `<fasta_path>\t<haplotype_label|AUTO>`, one row per assembly file — see
`manifests/daBudDavi1.tsv` for a real example: a tetraploid *Buddleja davidii* assembly where one
file holds haplotype 1 (headers say `chromosome: N`) and a second file holds haplotypes 2–4
together (headers say `HAP2_SUPER_N` / `HAP3_SUPER_N` / `HAP4_SUPER_N`, so that file's haplotype
label is `AUTO`-detected per sequence). `--chrom-regex` / `--hap-regex` override the default header
patterns for genomes with different naming conventions.


https://github.com/KamilSJaron/k-mer-approaches-for-biodiversity-genomics/wiki/Separate-sub-genomes-of-an-allopolyploid
