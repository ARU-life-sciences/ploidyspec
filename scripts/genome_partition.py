#!/usr/bin/env python3
"""
Detects diffuse, genome-wide chromosome partitions -- structure where no
single pair clears FDR significance (homeolog_pairs.tsv is empty or sparse),
but the chromosomes as a whole split into >=2 groups that are mutually more
similar within than between. This is a different signal from homeolog_
pairs.tsv's discrete 1:1 ancestral-duplicate test: it's asking whether the
whole chromosome set factors into subgenome-sized blocks at all, at any
resolution k, not whether any specific pair is a retained duplicate.

Reuses homeologs/homeolog_candidates_ranked.tsv's `distance` column, which
already has every cross-chromosome-number pair tested (not just the FDR-
accepted subset) -- no new k-mer work.

Method: average-linkage agglomerative clustering over the complete distance
matrix, recording the partition at every k from N-1 down to 2 in a single
merge pass. At each k, computes the separation ratio (mean between-group
distance / mean within-group distance) for the greedy partition, then tests
it against a null of ~500 random partitions with the SAME group-size
profile, giving a z-score. High |z| at some k means the chromosomes
genuinely factor into k groups of that shape more than chance predicts --
whether that's k=2 (a simple diffuse bipartition, e.g. daInuConz1),
k=N/3 (e.g. drLytSali1's reported 5 groups of 3), or multiple significant k
simultaneously (nested structure, e.g. ddLepDrab1's pairs of already-
accepted pairs, a candidate allo-octaploid block structure).

Usage: python3 scripts/genome_partition.py [species ...] > meta/genome_partition.tsv
       (no species args = whole panel)
"""
import csv
import os
import random
import statistics
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(REPO, "results")

N_PERMUTATIONS = 500
MIN_Z = 2.0  # only report k's whose observed separation clears this


def load_distance_matrix(species_dir):
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


def null_distribution(dist, chroms, sizes, n=N_PERMUTATIONS):
    ratios = []
    pool = list(chroms)
    for _ in range(n):
        random.shuffle(pool)
        clusters = []
        idx = 0
        for s in sizes:
            clusters.append(pool[idx : idx + s])
            idx += s
        r = separation_ratio(dist, clusters)
        if r is not None:
            ratios.append(r)
    return ratios


def analyze_species(species_dir, species_name):
    dist, chroms = load_distance_matrix(species_dir)
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
        null = null_distribution(dist, chroms, sizes)
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
    return results


def main():
    random.seed(20260906)  # fixed seed: z-scores are cited by value in INTERPRETATION.md
    species_args = sys.argv[1:]
    if species_args:
        species_list = species_args
    else:
        species_list = sorted(
            d for d in os.listdir(RESULTS_DIR) if os.path.isdir(os.path.join(RESULTS_DIR, d))
        )

    fieldnames = [
        "species",
        "k",
        "group_sizes",
        "separation_ratio",
        "null_mean_ratio",
        "z_score",
        "groups",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for species in species_list:
        species_dir = os.path.join(RESULTS_DIR, species)
        if not os.path.isdir(species_dir):
            continue
        results = analyze_species(species_dir, species)
        results.sort(key=lambda r: -abs(r["z_score"]))
        for r in results:
            if abs(r["z_score"]) < MIN_Z:
                continue
            writer.writerow(
                {
                    "species": r["species"],
                    "k": r["k"],
                    "group_sizes": r["group_sizes"],
                    "separation_ratio": f"{r['separation_ratio']:.4f}",
                    "null_mean_ratio": f"{r['null_mean_ratio']:.4f}",
                    "z_score": f"{r['z_score']:.2f}",
                    "groups": " | ".join(",".join(g) for g in r["groups"]),
                }
            )


if __name__ == "__main__":
    main()
