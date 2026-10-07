"""
Per-species genome-structure outputs read from the matrix and homeologs stages
(no new k-mer work). The cross-species tables are built by `ploidyspec panel`.

Writes, in results/<species>/structure/:
- genome_partition.tsv: diffuse genome-wide chromosome partitions -- whether the
  chromosome numbers factor into k groups more than chance, at every k with
  |z| >= 2 (INTERPRETATION.md, "Diffuse, genome-wide partitions").
- pair_synchrony.tsv: how uniform divergence depth is across the accepted
  homeolog pairs (pair_depth_cv). One whole-genome duplication diverges every
  pair to about the same depth; several events, or pairs resolving at
  different times, spread them out.

The earlier inheritance-mode metrics (partition_consistency, run length /
flip rate, windowed CVs, distance_ratio_cv, combined_allo_score) were retired
on 2026-10-07 (ROADMAP Phase 2, option B): partition_consistency and the
windowed metrics were driven by haplotype-labelling and window-registration
artefacts (INTERPRETATION.md, "Phase 2 hand-check: corrections"), and
distance_ratio_cv is confounded by tetrasomic homogenization. Their last values
are archived in meta/archive/auto_allo_spectrum.tsv.
"""
import csv
import os
import random
import statistics

from .common import homeologs_dir, log

N_PERMUTATIONS = 500
MIN_Z = 2.0  # only report k's whose observed separation clears this
PARTITION_SEED = 20260906


def structure_dir(outdir):
    d = os.path.join(outdir, "structure")
    os.makedirs(d, exist_ok=True)
    return d


def metric_pair_depth_cv(species_dir):
    pairs_tsv = os.path.join(homeologs_dir(species_dir), "homeolog_pairs.tsv")
    if not os.path.exists(pairs_tsv):
        return None, 0
    with open(pairs_tsv) as f:
        pairs = list(csv.DictReader(f, delimiter="\t"))
    depths = [float(r["mean_distance"]) for r in pairs if r.get("mean_distance")]
    if len(depths) < 2:
        return None, len(depths)
    mean = statistics.mean(depths)
    if mean == 0:
        return None, len(depths)
    return statistics.pstdev(depths) / mean, len(depths)


def load_chrom_distance_matrix(species_dir):
    path = os.path.join(species_dir, "homeologs", "homeolog_candidates_ranked.tsv")
    if not os.path.exists(path):
        return None, None
    dist = {}
    chroms = set()
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            a, b = r["chrom_a"], r["chrom_b"]
            d = float(r["distance"])
            dist[(a, b)] = d
            dist[(b, a)] = d
            chroms.add(a)
            chroms.add(b)
    return dist, sorted(chroms, key=lambda x: int(x[3:]))


def cluster_mean_dist(dist, c1, c2):
    vals = [dist[(a, b)] for a in c1 for b in c2 if (a, b) in dist]
    if not vals:
        return None
    return sum(vals) / len(vals)


def agglomerative_merge_sequence(dist, chroms):
    """Returns a list of (k, clusters) from k=len(chroms) down to k=1,
    average-linkage, greedy nearest-pair merge at every step."""
    clusters = [[c] for c in chroms]
    sequence = [(len(clusters), [list(c) for c in clusters])]
    while len(clusters) > 1:
        best = None
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d = cluster_mean_dist(dist, clusters[i], clusters[j])
                if d is None:
                    continue
                if best is None or d < best[0]:
                    best = (d, i, j)
        if best is None:
            break
        _, i, j = best
        clusters[i] = clusters[i] + clusters[j]
        del clusters[j]
        sequence.append((len(clusters), [list(c) for c in clusters]))
    return sequence


def separation_ratio(dist, clusters):
    within = []
    between = []
    for idx, c in enumerate(clusters):
        for a_i, a in enumerate(c):
            for b in c[a_i + 1 :]:
                if (a, b) in dist:
                    within.append(dist[(a, b)])
        for other in clusters[idx + 1 :]:
            for a in c:
                for b in other:
                    if (a, b) in dist:
                        between.append(dist[(a, b)])
    if not within or not between:
        return None
    return statistics.mean(between) / statistics.mean(within)


def null_distribution(dist, chroms, sizes, rng, n=N_PERMUTATIONS):
    ratios = []
    pool = list(chroms)
    for _ in range(n):
        rng.shuffle(pool)
        clusters = []
        idx = 0
        for s in sizes:
            clusters.append(pool[idx : idx + s])
            idx += s
        r = separation_ratio(dist, clusters)
        if r is not None:
            ratios.append(r)
    return ratios


def analyze_species(species_dir, species_name, seed=PARTITION_SEED):
    """Every k with a valid partition, strongest |z| first. The null uses its own
    RNG seeded per species, so a species' z-scores don't depend on which other
    species were analysed before it (the old panel script shared one RNG)."""
    rng = random.Random(seed)
    dist, chroms = load_chrom_distance_matrix(species_dir)
    if dist is None or len(chroms) < 4:
        return []

    sequence = agglomerative_merge_sequence(dist, chroms)
    results = []
    for k, clusters in sequence:
        if k < 2 or k > len(chroms) - 2:
            continue
        sizes = sorted(len(c) for c in clusters)
        if min(sizes) < 2:
            continue
        obs = separation_ratio(dist, clusters)
        if obs is None:
            continue
        null = null_distribution(dist, chroms, sizes, rng)
        if len(null) < 10:
            continue
        null_mean = statistics.mean(null)
        null_std = statistics.pstdev(null)
        if null_std == 0:
            continue
        z = (obs - null_mean) / null_std
        results.append(
            {
                "species": species_name,
                "k": k,
                "group_sizes": ",".join(str(s) for s in sizes),
                "separation_ratio": obs,
                "null_mean_ratio": null_mean,
                "z_score": z,
                "groups": clusters,
            }
        )
    results.sort(key=lambda r: -abs(r["z_score"]))
    return results



PARTITION_FIELDS = ["species", "k", "group_sizes", "separation_ratio", "null_mean_ratio",
                    "z_score", "groups"]


SYNCHRONY_FIELDS = ["species", "n_accepted_pairs", "pair_depth_median", "pair_depth_cv"]


def pair_synchrony(species_dir, species):
    cv, n = metric_pair_depth_cv(species_dir)
    median = None
    pairs_tsv = os.path.join(homeologs_dir(species_dir), "homeolog_pairs.tsv")
    if os.path.exists(pairs_tsv):
        with open(pairs_tsv) as f:
            depths = [float(r["mean_distance"]) for r in csv.DictReader(f, delimiter="\t")
                      if r.get("mean_distance")]
        median = statistics.median(depths) if depths else None
    return {
        "species": species,
        "n_accepted_pairs": n,
        "pair_depth_median": f"{median:.4f}" if median is not None else "",
        "pair_depth_cv": f"{cv:.4f}" if cv is not None else "",
    }


def partition_rows(species_dir, species):
    rows = []
    for r in analyze_species(species_dir, species):
        if abs(r["z_score"]) < MIN_Z:
            continue
        rows.append({
            "species": species,
            "k": r["k"],
            "group_sizes": r["group_sizes"],
            "separation_ratio": f"{r['separation_ratio']:.4f}",
            "null_mean_ratio": f"{r['null_mean_ratio']:.4f}",
            "z_score": f"{r['z_score']:.2f}",
            "groups": " | ".join(",".join(g) for g in r["groups"]),
        })
    return rows


def write_rows(path, fields, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def compute_structure(outdir):
    species = os.path.basename(os.path.normpath(outdir))
    sdir = structure_dir(outdir)
    write_rows(os.path.join(sdir, "pair_synchrony.tsv"), SYNCHRONY_FIELDS,
               [pair_synchrony(outdir, species)])
    parts = partition_rows(outdir, species)
    write_rows(os.path.join(sdir, "genome_partition.tsv"), PARTITION_FIELDS, parts)
    stale = os.path.join(sdir, "inheritance_metrics.tsv")
    if os.path.exists(stale):
        os.remove(stale)
    log(f"wrote pair_synchrony.tsv and genome_partition.tsv "
        f"({len(parts)} partition(s) with |z| >= {MIN_Z}) in {sdir}")
