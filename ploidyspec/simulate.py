"""
Synthetic haplotype-resolved assemblies with a known answer, for end-to-end
testing and for calibrating thresholds against ground truth.

Each scenario writes one FASTA per haplotype (`chromosome: N` in placed
sequences' headers, the format `prepare` reads by default), a manifest, and
`truth.tsv` listing what the pipeline should recover. Sequence is random
(no real genome needed) with transposable-element-like repeat families
inserted at high copy number so the te-markers stage has repeat content to
work with.

Scenarios:
- diploid: 2 haplotypes of one genome.
- autotetraploid: 4 haplotypes, all equally related -- polysomic,
  every chromosome tetrasomic_like.
- allotetraploid: two progenitor genomes diverged from a common ancestor, each
  with its own TE burst, numbered as separate chromosomes (1..n from A,
  n+1..2n from B) in 2 haplotypes -- the usual layout of a curated
  allopolyploid assembly. Homeolog pairs (i, i+n) are the expected finding.
- rediploidized: an autotetraploid where HAP3/HAP4 have partly split off as a
  second lineage, one mechanism per chromosome:
    chr1  whole-chromosome lineage split           -> resolved_lineages
    chr2  split over the first 40% only            -> partially_resolved
    chr3+chr4 fused in HAP3/HAP4, placed as chr3    -> fusion_lineages (placed long copy)
    chr5+chr6 fused in HAP3/HAP4, unplaced scaffold -> fusion_lineages (unplaced scaffold)
    remaining chromosomes untouched                 -> tetrasomic_like
"""

import csv
import os

import numpy as np

from .common import log

SCENARIOS = ("diploid", "autotetraploid", "allotetraploid", "rediploidized")
LETTERS = np.frombuffer(b"ACGT", dtype=np.uint8)

DEFAULTS = dict(
    n_chrom=8,
    chrom_len=1_500_000,
    het=0.002,  # per-haplotype substitution rate from its lineage
    lineage_div=0.02,  # divergence of a rediploidized lineage
    allo_div=0.03,  # each progenitor's divergence from the common ancestor
    te_families=2,
    te_len=1000,
    te_copies=120,  # per chromosome per family: clears te-markers' --min-count 100
    te_copy_div=0.01,
)


class Simulator:
    def __init__(self, seed, **params):
        self.rng = np.random.default_rng(seed)
        self.p = {**DEFAULTS, **params}

    def random_seq(self, n):
        return self.rng.choice(4, size=n, p=[0.3, 0.2, 0.2, 0.3]).astype(np.uint8)

    def mutate(self, seq, rate, start=0, end=None):
        """Substitutions at `rate` within seq[start:end]; returns a new array."""
        out = seq.copy()
        end = len(seq) if end is None else end
        hit = np.flatnonzero(self.rng.random(end - start) < rate) + start
        out[hit] = (out[hit] + self.rng.integers(1, 4, size=len(hit))) % 4
        return out

    def te_family(self):
        return self.random_seq(self.p["te_len"])

    def insert_copies(self, seq, family, copies):
        """Insert `copies` slightly diverged copies of `family` at random positions."""
        cuts = np.sort(self.rng.integers(0, len(seq), size=copies))
        parts, prev = [], 0
        for c in cuts:
            parts.append(seq[prev:c])
            parts.append(self.mutate(family, self.p["te_copy_div"]))
            prev = c
        parts.append(seq[prev:])
        return np.concatenate(parts)

    def genome(self, n_chrom, families):
        chroms = []
        for _ in range(n_chrom):
            s = self.random_seq(self.p["chrom_len"])
            for fam in families:
                s = self.insert_copies(s, fam, self.p["te_copies"])
            chroms.append(s)
        return chroms

    def burst(self, seq):
        """Lineage-specific TE burst: a new family, private to this lineage."""
        return self.insert_copies(seq, self.te_family(), self.p["te_copies"])

    def haplotypes(self, chroms, n_hap):
        return [[self.mutate(c, self.p["het"]) for c in chroms] for _ in range(n_hap)]


def to_text(seq):
    return LETTERS[seq].tobytes().decode()


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, desc, seq in records:
            f.write(f">{name} {desc}\n" if desc else f">{name}\n")
            text = to_text(seq)
            for i in range(0, len(text), 80):
                f.write(text[i : i + 80] + "\n")


def build_scenario(name, seed, params):
    """Returns ({hap: [(seq_name, description, seq)]}, [truth rows])."""
    sim = Simulator(seed, **params)
    p = sim.p
    n = p["n_chrom"]
    ancestral = [sim.te_family() for _ in range(p["te_families"])]

    def placed(hap, chrom, seq):
        return (f"{hap}_chr{chrom}", f"chromosome: {chrom}", seq)

    haps, truth = {}, []
    if name in ("diploid", "autotetraploid"):
        n_hap = 2 if name == "diploid" else 4
        genome = sim.genome(n, ancestral)
        for i, chroms in enumerate(sim.haplotypes(genome, n_hap), 1):
            haps[f"HAP{i}"] = [placed(f"HAP{i}", c, s) for c, s in enumerate(chroms, 1)]
        state = "not_assessable" if n_hap == 2 else "tetrasomic_like"
        truth = [dict(chrom=f"chr{c:02d}", copy_state=state, fusion="", homeolog="")
                 for c in range(1, n + 1)]

    elif name == "allotetraploid":
        ancestor = sim.genome(n, ancestral)
        progenitors = []
        for _ in range(2):
            chroms = [sim.burst(sim.mutate(c, p["allo_div"])) for c in ancestor]
            progenitors.append(chroms)
        genome = progenitors[0] + progenitors[1]
        for i, chroms in enumerate(sim.haplotypes(genome, 2), 1):
            haps[f"HAP{i}"] = [placed(f"HAP{i}", c, s) for c, s in enumerate(chroms, 1)]
        for c in range(1, 2 * n + 1):
            partner = c + n if c <= n else c - n
            truth.append(dict(chrom=f"chr{c:02d}", copy_state="not_assessable",
                              fusion="", homeolog=f"chr{partner:02d}"))

    elif name == "rediploidized":
        if n < 6:
            raise SystemExit("rediploidized scenario needs --n-chrom >= 6")
        genome = sim.genome(n, ancestral)
        d = p["lineage_div"]
        L = p["chrom_len"]
        # lineage X (HAP1/HAP2) is the ancestral genome; lineage Y (HAP3/HAP4)
        # differs from it on chr1-chr6 only
        y = list(genome)
        y[0] = sim.burst(sim.mutate(genome[0], d))
        y[1] = sim.mutate(genome[1], d, 0, int(0.4 * L))
        fused_34 = sim.burst(sim.mutate(np.concatenate([genome[2], genome[3]]), d))
        fused_56 = sim.burst(sim.mutate(np.concatenate([genome[4], genome[5]]), d))
        for i in (1, 2):
            chroms = [sim.mutate(c, p["het"]) for c in genome]
            haps[f"HAP{i}"] = [placed(f"HAP{i}", c, s) for c, s in enumerate(chroms, 1)]
        for i in (3, 4):
            hap = f"HAP{i}"
            recs = [placed(hap, 1, sim.mutate(y[0], p["het"])),
                    placed(hap, 2, sim.mutate(y[1], p["het"])),
                    placed(hap, 3, sim.mutate(fused_34, p["het"])),
                    (f"{hap}_fused_5_6", "", sim.mutate(fused_56, p["het"]))]
            recs += [placed(hap, c, sim.mutate(genome[c - 1], p["het"])) for c in range(7, n + 1)]
            haps[hap] = recs
        expected = {1: ("resolved_lineages", ""), 2: ("partially_resolved", ""),
                    3: ("fusion_lineages", "chr04"), 4: ("fusion_lineages", "chr03"),
                    5: ("fusion_lineages", "chr06"), 6: ("fusion_lineages", "chr05")}
        for c in range(1, n + 1):
            state, partner = expected.get(c, ("tetrasomic_like", ""))
            truth.append(dict(chrom=f"chr{c:02d}", copy_state=state, fusion=partner, homeolog=""))
    else:
        raise SystemExit(f"unknown scenario {name!r}; choose from {', '.join(SCENARIOS)}")
    return haps, truth


def simulate(scenario, outdir, seed, params):
    names = SCENARIOS if scenario == "all" else (scenario,)
    for i, name in enumerate(names):
        sdir = os.path.abspath(os.path.join(outdir, name))
        os.makedirs(sdir, exist_ok=True)
        log(f"simulating {name} -> {sdir}")
        haps, truth = build_scenario(name, seed + i, params)
        with open(os.path.join(sdir, "manifest.tsv"), "w") as m:
            m.write("# fasta_path\thap_label\n")
            for hap, records in haps.items():
                path = os.path.join(sdir, f"{name}.{hap}.fa")
                write_fasta(path, records)
                m.write(f"{path}\t{hap}\n")
        with open(os.path.join(sdir, "truth.tsv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["chrom", "copy_state", "fusion", "homeolog"],
                               delimiter="\t")
            w.writeheader()
            w.writerows(truth)
    return [os.path.join(outdir, n) for n in names]
