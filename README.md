# ploidyspec

K-mer-based subgenome/ploidy divergence profiling for polyploid genomes: given two or
more haplotype-scale genome assemblies of the same individual, ploidyspec measures how
diverged each chromosome's copies are (Mash-corrected k-mer distance), searches for
retained ancestral (paleopolyploid) chromosome pairs, and — where enough repeat content
exists — estimates a continuous auto&harr;allo index from differential fossil-TE k-mer
markers between subgenomes. No reference genome or alignment required.

See [`OUTPUTS.md`](OUTPUTS.md) for what every output file/column means, with worked
examples from real runs. This README covers installation and how to run it.

## Install / setup

Python dependencies are stdlib plus `numpy` and `matplotlib` only (no `scipy`).

External tools this shells out to: `samtools`, and FastK's `FastK`/`Logex`/`Histex`/
`Tabex` (github.com/thegenemyers/FASTK). Ploidyspec auto-detects these three ways, tried
in order:
1. `--samtools <path>` / `--fastk-dir <path>` (a directory containing the FastK
   binaries) passed explicitly.
2. Sibling directories next to the repo checkout matching `samtools*` / `FASTK*` (e.g.
   `../samtools-1.24/samtools`, `../FASTK-1.2/FastK` — this is how the binaries are laid
   out on this project's cluster; see `.gitignore`'s `/samtools-1.24` and `/FASTK-1.2`
   entries).
3. `$PATH`.

No install step beyond having those binaries reachable one of those three ways —
`ploidyspec.sh` is a thin wrapper (`exec python3 -m ploidyspec.cli "$@"`), nothing to
build.

## Naming your chromosome-scale sequences

Ploidyspec needs to know, per sequence, which **haplotype copy** it is and which
**chromosome number** it represents. Both are read from the FASTA header/description via
regex — get your assembly's naming close to one of the patterns below and you won't need
to configure anything.

**Chromosome number** — either of:
- `chromosome: <N>` (or `chromosome<N>`, case-insensitive) anywhere in the description.
- `SUPER_<N>` or `SUPER-<N>` in the sequence name, as long as nothing else immediately
  follows the number (so a compound name like `SUPER_4_HAP2` does *not* match this one —
  see below).

**Haplotype number** — only used when a manifest row is marked `AUTO` (see below):
`HAP<N>_SUPER` or `HAP<N>-SUPER` in the sequence name (haplotype number comes first).

If your assembly's headers don't match either pattern — e.g. the opposite naming order
(`SUPER_<N>_HAP<M>`, chromosome-then-haplotype, which several species in this project's
own panel actually use) — override with `--chrom-regex`/`--hap-regex` (each takes a
regex with exactly one capturing group, the number). Two real examples from this
project's `scripts/jobs.tsv`:

```
# daGalBore1: headers mix "chromosome: N" free text with plain "SUPER_N" scaffold
# names -- one combined --chrom-regex covering both, tried against each header in order:
--chrom-regex 'chromosome:?\s*(\d+)\b|SUPER[_-](\d+)_HAP\d+\b|SUPER[_-](\d+)\b(?!_)'

# drLytSali1: scaffold names are "SUPER_<N>_HAP<M>" (chromosome-number-first, the
# opposite order the default --hap-regex expects), so both are overridden:
--chrom-regex 'chromosome:?\s*(\d+)\b|SUPER[_-](\d+)_HAP\d+\b|SUPER[_-](\d+)\b(?!_)'
--hap-regex '(?:SUPER_\d+_)?HAP(\d+)(?:_SUPER)?'
```

`--chrom-regex` may be passed multiple times (each tried in order, first match wins) or
combined into one regex with `|` alternation as above. After a `prepare` run, check
`unplaced.tsv` in the output directory — if far more sequences ended up there than
expected, your regex isn't matching and needs adjusting before continuing.

A **manifest** is a TSV of `<fasta_path>\t<haplotype_label|AUTO>`, one row per assembly
file — see `manifests/daBudDavi1.tsv` for a real example. Use an explicit label
(`HAP1`, `HAP2`, ...) when one file is exactly one haplotype copy; use `AUTO` when a
single file contains multiple haplotype copies distinguished by header (matched via
`--hap-regex`/the default above).

## Pipeline stages

Cheap stages, always run by `all`:
1. **prepare** — parse the manifest into per-chromosome units (`sequences.tsv`) using
   the naming rules above, filtered to sequences &ge; `--min-len` (default 1 Mb).
2. **kmers** — canonical FastK k-mer table per unit.
3. **matrix** — all-vs-all whole-chromosome Mash-corrected k-mer distance, ploidy
   summary, homologous-chromosome-pair report, clustered + contrast heatmaps. Includes an
   automatic correctness check (`chrom_reconcile.py`) for cross-file chromosome-number
   mislabeling — see `OUTPUTS.md` if a run reports corrections.
4. **homeologs** — FDR-controlled search for retained ancestral (paleopolyploid)
   chromosome pairs among *different* chromosome numbers, plus the full ranked-candidate
   evidence (not just the FDR-accepted subset).
5. **windowed** — sliding-window divergence along each chromosome between its haplotype
   copies.
6. **report** — self-contained HTML report (`report.html` in the output directory),
   embedding whatever plots/tables the stages that ran produced.

Opt-in stages (moderate-to-expensive; flags on `all`, or run individually):
- `--with-windowed-homeologs` — sliding-window divergence between the ancestral pairs
  `homeologs` found.
- `--with-te-markers` — differential fossil-TE k-mer markers between same-chromosome
  haplotype copies (the auto&harr;allo index), plus `subgenome-report`'s summary.
- `--with-te-markers-windowed` — the most expensive stage: paints windows along each
  haplotype copy by which subgenome's markers they match (implies
  `--with-te-markers`).

## Usage

```
./ploidyspec.sh all \
  --manifest manifests/daBudDavi1.tsv \
  --outdir results/daBudDavi1 \
  --k 15 --min-len 1000000 \
  --window 250000 \
  --threads 8
```

Add `--with-te-markers` (or `--with-te-markers-windowed`) to also get the auto&harr;allo
index in the same run. Every stage is also runnable individually with the same flags
(`./ploidyspec.sh <stage> --manifest ... --outdir ...`) — each stage skips work it's
already done, so re-running after an interruption (or after adding `--with-*` flags to
an already-completed run) picks up where it left off rather than redoing everything.

See [`OUTPUTS.md`](OUTPUTS.md) for what every output file and column means, with worked
examples from real species in this project's panel.

---

Method reference for the fossil-TE marker stages:
https://github.com/KamilSJaron/k-mer-approaches-for-biodiversity-genomics/wiki/Separate-sub-genomes-of-an-allopolyploid
