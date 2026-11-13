# ploidyspec

K-mer-based subgenome/ploidy divergence profiling for polyploid genomes: given two or
more haplotype-scale genome assemblies of the same individual, ploidyspec measures how
diverged each chromosome's copies are (Mash-corrected k-mer distance), searches for
retained ancestral (paleopolyploid) chromosome pairs, and — where enough repeat content
exists — estimates a continuous auto&harr;allo index from differential fossil-TE k-mer
markers between subgenomes. No reference genome or alignment required.

See [`OUTPUTS.md`](OUTPUTS.md) for what every output file/column means, with worked
examples from real runs. This README covers installation and how to run it.

## Quick start

```
conda env create -f environment.yml && conda activate ploidyspec && pip install .
ploidyspec check                                   # finds samtools, FastK and Python deps
ploidyspec simulate --scenario allotetraploid --outdir demo
ploidyspec all --fasta demo/allotetraploid/*.fa --outdir demo/out --threads 4
cat demo/out/summary.md                            # compare with demo/allotetraploid/truth.tsv
```

That takes about 7 minutes on 4 cores. For your own genome, pass one
FASTA per haplotype copy (plain, gzip or bgzip), in order:

```
ploidyspec all --fasta hap1.fa.gz hap2.fa.gz --outdir results/mysp --threads 8
```

and open `results/mysp/report.html`. Add `--full` for the further readings (residual
tetrasomy, homeologous exchange, the homeology map, the TE-marker auto/allo index; see
below for the cost). `--full` leaves out the windowed TE painting
(`--with-te-markers-windowed`), which is the most expensive stage and only informative
when the subgenomes carry clearly different repeats.

## Install / setup

With conda (recommended — brings samtools and FastK from bioconda):

```
git clone https://github.com/ARU-life-sciences/ploidyspec.git && cd ploidyspec
conda env create -f environment.yml
conda activate ploidyspec
pip install .
ploidyspec check
```

Or, if you already have `samtools` and FastK, `pip install
git+https://github.com/ARU-life-sciences/ploidyspec.git` (Python ≥3.9; the only Python
dependencies are `numpy` and `matplotlib`). Cloning the repository downloads about 1 GB,
almost all of it the panel's results; the pip line fetches the same, so for the tool
alone `git clone --depth 1` is quicker. Running from a checkout without installing also
works: `./ploidyspec.sh` is a thin wrapper around `python3 -m ploidyspec.cli`.
`ploidyspec check` lists anything missing and how to install it; `ploidyspec --version`
prints the version.

External tools this shells out to: `samtools`, and FastK's `FastK`/`Logex`/`Histex`/
`Tabex` (github.com/thegenemyers/FASTK). Ploidyspec finds them three ways, tried in
order:
1. `--samtools <path>` / `--fastk-dir <path>` (a directory containing the FastK
   binaries) passed explicitly.
2. Directories inside the repo checkout matching `samtools*` / `FASTK*` (e.g.
   `samtools-1.24/samtools`, `FASTK-1.2/FastK` — see `.gitignore`).
3. `$PATH` (what the conda environment provides).

## Naming your chromosome-scale sequences

Ploidyspec needs to know, per sequence, which **haplotype copy** it is and which
**chromosome number** it represents. The haplotype comes from the manifest. The chromosome
number can come from the headers, or ploidyspec can work it out:

- **No names needed** (`--chrom-naming auto`, and the default whenever no header carries a
  chromosome number). In each haplotype, sequences of at least `--min-len` that are also
  at least a tenth of the median length of its larger sequences count as chromosome-scale;
  the rest go to `unplaced.tsv`. The haplotype with the most chromosome-scale sequences is
  numbered by length (1 = longest), and every other haplotype's sequences take the number
  of the reference chromosome they share most k-mers with (one-to-one matching in the
  `matrix` stage). Tested on simulated assemblies with every name removed and the
  sequences shuffled. A sequence that matches nothing keeps a provisional number from
  1001 and is reported in the log.
- **From headers** (`--chrom-naming names`, or the default when headers carry numbers),
  read via regex as below. Get your assembly's naming close to one of these patterns and
  you won't need to configure anything.

Each haplotype must still be chromosome-scale. A haplotype assembled only to contigs
can be placed on another, scaffolded haplotype of the same individual first with
`workflows/prep/scaffold_by_reference.py` (minimap2).

**Chromosome number** — any of:
- `chromosome: <N>` (or `chromosome<N>`, case-insensitive) anywhere in the description.
- `chromosome: <N><L>` or `chromosome: <L><N>` with a subgenome letter A–I (wheat `3B`,
  bird `1A`, *Brassica* `C09`, cotton `D13`). These are numbered 100 × letter + number
  (A = 1), so `3B` becomes chromosome 203 (`chr203` in unit names) and `C09` becomes 309.
- `SUPER_<N>_HAP<M>` (chromosome-then-haplotype, as several ToL assemblies name them).
- `SUPER_<N>` or `SUPER-<N>` in the sequence name, as long as nothing else immediately
  follows the number (so `SUPER_22_unloc_1` does not match).

**Haplotype number** — only used when a manifest row is marked `AUTO` (see below):
`HAP<N>_SUPER` or `HAP<N>-SUPER` in the sequence name (haplotype number comes first).

If your assembly's headers match none of these, override with
`--chrom-regex`/`--hap-regex` (a regex with one capturing group for the number; patterns
combined with `|` may each have their own group, and the one that matched is used). A
real example from this project's `workflows/sanger/jobs.tsv`:

```
# drLytSali1: manifest rows labelled AUTO, scaffolds named "SUPER_<N>_HAP<M>", so the
# haplotype is read from the name, in the opposite order to the default --hap-regex:
--hap-regex '(?:SUPER_\d+_)?HAP(\d+)(?:_SUPER)?'
```

`--chrom-regex` may be passed multiple times (each tried in order, first match wins) or
combined into one regex with `|` alternation. After a `prepare` run, check
`unplaced.tsv` in the output directory — if far more sequences ended up there than
expected, your regex isn't matching and needs adjusting before continuing.

**Input files.** Each haplotype must be one FASTA file, plain, gzip or bgzip. Plain
gzip cannot be indexed by samtools, so it is recompressed with bgzip (or, without bgzip,
decompressed) into `<outdir>/inputs/`; a FASTA in a directory you cannot write to is
linked from there so its index can be written. With `--fasta hap1.fa hap2.fa ...` the
files are labelled HAP1, HAP2, ... in order (written to `<outdir>/manifest.tsv`). If one
haplotype's headers carry chromosome numbers and another's do not, all are numbered
automatically.

A **manifest** (`--manifest`, instead of `--fasta`) is a TSV of
`<fasta_path>\t<haplotype_label|AUTO>`, one row per assembly file — see `manifests/daBudDavi1.tsv` for a real example. Use an explicit label
(`HAP1`, `HAP2`, ...) when one file is exactly one haplotype copy; use `AUTO` when a
single file contains multiple haplotype copies distinguished by header (matched via
`--hap-regex`/the default above).

## Testing on simulated genomes

`ploidyspec simulate --outdir sims` writes five synthetic assemblies with a known
answer — diploid, autotetraploid (as four haplotypes, and as two haplotypes with each
chromosome numbered twice), allotetraploid, and a partially rediploidized
autotetraploid (a whole-chromosome lineage split, a regional split, one placed and
one unplaced fusion) — each with a `manifest.tsv` and the expected results in
`truth.tsv`. Run any stage on them like a real species. The end-to-end test runs the
whole pipeline on all five and checks it recovers the truth (needs samtools + FastK,
about 25 minutes):

```
python3 -m unittest discover -s tests                         # unit tests, seconds
PLOIDYSPEC_E2E=1 python3 -m unittest tests.test_end_to_end    # end-to-end
```

## Pipeline stages

Cheap stages, always run by `all`:
1. **prepare** — parse the manifest into per-chromosome units (`sequences.tsv`) using
   the naming rules above, filtered to sequences &ge; `--min-len` (default 1 Mb).
2. **kmers** — canonical FastK k-mer table per unit.
3. **matrix** — all-vs-all whole-chromosome Mash-corrected k-mer distance, ploidy
   summary, homologous-chromosome-pair report, clustered + contrast heatmaps. Includes an
   automatic correctness check (`chrom_reconcile.py`) for cross-file chromosome-number
   mislabeling — see `OUTPUTS.md` if a run reports corrections. A haplotype file
   holding two chromosome sets under different numbers (an AAB triploid assembled as
   A+B in one file) is split into an extra copy, so odd ploidy reads as three copies.
4. **homeologs** — FDR-controlled search for retained ancestral (paleopolyploid)
   chromosome pairs among *different* chromosome numbers, plus the full ranked-candidate
   evidence (not just the FDR-accepted subset).
5. **windowed** — divergence along each chromosome: every copy's windows are looked up
   in each other copy's whole chromosome, so assemblies needn't be collinear.
6. **structure** — diffuse genome-wide chromosome partitions, and how uniform
   divergence depth is across the accepted homeolog pairs. No new k-mer work.
7. **rediploidization** — chromosome fusions between haplotype copies (k-mer
   containment), lineage structure along each chromosome, and a per-chromosome
   rediploidization state (`fusion_lineages`, `partially_resolved`,
   `tetrasomic_like`, ...) alongside ancient pairing. Uses `te-markers` output
   if present.
8. **summary** — plain-language answers to the four questions (ploidy, auto- or
   allo-like structure, TE markers, rediploidization), each with a confidence and the
   evidence behind it (`summary.tsv`, `summary.md`; also at the top of the report).
9. **report** — self-contained HTML report (`report.html` in the output directory),
   embedding whatever plots/tables the stages that ran produced.

Opt-in stages (moderate-to-expensive; flags on `all`, or run individually):
- `--with-windowed-homeologs` — a genome-wide homeology map first (every window of
  every chromosome of one haplotype against every other chromosome, k=23: which
  chromosome each block is duplicated on, so homeology that runs per arm is found),
  then sliding-window divergence between the pairs `homeologs` and the map found.
  These tracks feed the residual-tetrasomy and homeologous-exchange readings of
  `rediploidization`. The map's cost grows with the square of the chromosome count
  (about 2 h on 8 threads for 40 chromosomes of ~100 Mb).
- `--with-te-markers` — differential fossil-TE k-mer markers between same-chromosome
  haplotype copies (the auto&harr;allo index), plus `subgenome-report`'s summary.
- `--with-te-markers-windowed` — the most expensive stage: paints windows along each
  haplotype copy by which subgenome's markers they match (implies
  `--with-te-markers`; not part of `--full`). It separates subgenomes window by window
  only where their repeats differ strongly (durum wheat A vs B: 87% of windows on their
  own side, 0.3% on the other). Elsewhere the labels mostly follow which copy has more
  repeat k-mers, so it does not confirm or rule out the homeologous-exchange calls.

## Usage

```
ploidyspec all \
  --manifest manifests/daBudDavi1.tsv \
  --outdir results/daBudDavi1 \
  --threads 8
```

(or `--fasta hap1.fa hap2.fa` in place of `--manifest`).

Defaults are the settings the project's panel used: chromosome-scale = >= 1 Mb,
250 kb windows, window k = 15, marker k = 13. The one difference is the k-mer sweep
(`--k`, default `15,23`). The panel ran `11,13,15,17,19,23`, but the distance is taken
from the largest k that clears the chance-collision floor, and that was k = 23 for all
53,309 chromosome pairs in the panel. The default gives identical distances for about a
third of the k-mer counting. Pass the six-value sweep to reproduce the panel's runs
exactly.

Add `--with-te-markers` (or `--with-te-markers-windowed`) to also get the auto&harr;allo
index in the same run; `--full` turns on `--with-windowed-homeologs` and
`--with-te-markers` together. Every stage is also runnable individually with the same flags
(`./ploidyspec.sh <stage> --manifest ... --outdir ...`) — each stage skips work it's
already done, so re-running after an interruption (or after adding `--with-*` flags to
an already-completed run) picks up where it left off rather than redoing everything.

For several species, `ploidyspec panel --results results --outdir meta` builds the
cross-species tables (per-species summary, genome partitions, rediploidization states)
and, under `supplementary/`, the TE-marker table and an exploratory PCA.

See [`OUTPUTS.md`](OUTPUTS.md) for what every output file and column means, with worked
examples from real species in this project's panel.

## Runtime, memory and disk

From the panel's runs (`all` with the six-k sweep and the windowed stage, 8 threads,
Sanger farm):

- **Memory** scales with the longest chromosome, not total genome size: roughly
  40 MB of RAM per Mb of the largest chromosome-scale sequence, with ~1-5 GB for
  typical 30-100 Mb chromosomes. Larger examples:
  - 359 Mb (`daInuConz1`): 13.7 GB
  - 458 Mb (`lpElePalu1`): 18.6 GB
  - 875 Mb (`ddHesMatr1`): 26.5 GB
  - 842 Mb (wheat): 32.6 GB

  Wheat ran out of memory at 4 GB in the `rediploidization` stage alone. Size jobs
  from the longest chromosome.
- **Runtime**: 0.2-13 hours per species, mostly k-mer counting and the windowed
  stage. Scales with genome size and copy number. The 15,23 default roughly halves
  the k-mer part.
- **Disk**: k-mer tables and per-chromosome FASTAs dominate (the 45-species panel's
  `results/` is ~1.6 TB). `ploidyspec all --cleanup`, or `ploidyspec cleanup
  --outdir <dir>` afterwards, deletes them and keeps every table and plot. Stages
  re-run later rebuild what they need (re-run `kmers` before `matrix`).

---

Method reference for the fossil-TE marker stages:
https://github.com/KamilSJaron/k-mer-approaches-for-biodiversity-genomics/wiki/Separate-sub-genomes-of-an-allopolyploid

## Licence

MIT; see [`LICENSE`](LICENSE).
