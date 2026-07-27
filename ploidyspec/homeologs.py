import csv
import itertools
import os
import statistics
from collections import defaultdict

from .common import log
from .kmer_tables import load_sequences


def load_distance_matrix(outdir):
    path = os.path.join(outdir, "whole_chrom_distance_matrix.csv")
    with open(path) as f:
        r = list(csv.reader(f))
    ids = r[0][1:]
    mat = {ids[i]: {ids[j]: float(v) for j, v in enumerate(row[1:])} for i, row in enumerate(r[1:])}
    return ids, mat


def chrom_groups(units):
    groups = defaultdict(list)
    for u in units:
        groups[int(u["chrom"])].append(u["unit_id"])
    return groups


def cross_chrom_distances(groups, mat):
    """Mean k-mer distance between every pair of *different* chromosome numbers,
    averaged over every haplotype-copy combination between the two groups -- i.e.
    treats each chromosome *number* as a single node, independent of how many
    haplotype copies of it happen to be assembled."""
    chrom_nums = sorted(groups)
    pair_dist = {}
    for i, j in itertools.combinations(chrom_nums, 2):
        vals = [mat[a][b] for a in groups[i] for b in groups[j]]
        pair_dist[(i, j)] = statistics.mean(vals)
    return pair_dist


def detect_homeolog_pairs(pair_dist, chrom_nums, gap_search_frac=0.5):
    """
    Look for a retained ancestral (paleopolyploid/WGD) subgenome pairing among
    *different* chromosome numbers: chromosome pairs whose whole-chromosome k-mer
    distance sits well below the random cross-chromosome background, consistent
    with shared descent from a whole-genome duplication rather than coincidence.

    Approach: take the best (lowest-distance) partner candidates, find the largest
    gap in their sorted distances (real ancestral pairs cluster well below the
    random background, which is unimodal and much higher), keep everything below
    that gap, then greedily resolve it into a 1:1 matching.
    """
    sorted_pairs = sorted(pair_dist.items(), key=lambda kv: kv[1])
    n_candidates = max(len(chrom_nums) // 2 + 1, int(len(sorted_pairs) * gap_search_frac))
    search_space = sorted_pairs[:n_candidates]

    best_gap = -1.0
    cutoff_idx = len(search_space)
    for k in range(1, len(search_space)):
        gap = search_space[k][1] - search_space[k - 1][1]
        if gap > best_gap:
            best_gap = gap
            cutoff_idx = k
    pool = search_space[:cutoff_idx]

    matched = {}
    accepted = []
    for (i, j), d in pool:
        if i in matched or j in matched:
            continue
        matched[i] = j
        matched[j] = i
        accepted.append((i, j, d))

    unmatched = [c for c in chrom_nums if c not in matched]
    background = [d for _, d in sorted_pairs[cutoff_idx:]]
    return accepted, unmatched, background, best_gap


def run(seq_tsv, outdir):
    units = load_sequences(seq_tsv)
    ids, mat = load_distance_matrix(outdir)
    groups = chrom_groups(units)
    chrom_nums = sorted(groups)
    if len(chrom_nums) < 4:
        log(f"only {len(chrom_nums)} chromosome numbers -- too few to search for a paleopolyploid pairing")
        return []

    pair_dist = cross_chrom_distances(groups, mat)
    accepted, unmatched, background, gap = detect_homeolog_pairs(pair_dist, chrom_nums)

    path = os.path.join(outdir, "homeolog_pairs.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom_a", "chrom_b", "mean_distance", "n_haplotype_copy_pairs"])
        for i, j, d in sorted(accepted, key=lambda x: x[2]):
            w.writerow([f"chr{i:02d}", f"chr{j:02d}", f"{d:.6f}", len(groups[i]) * len(groups[j])])

    bg_mean = statistics.mean(background) if background else float("nan")
    bg_min = min(background) if background else float("nan")
    log(
        f"candidate ancient homeolog pairs: {len(accepted)} pairs across {2 * len(accepted)} of "
        f"{len(chrom_nums)} chromosomes (gap={gap:.4f}, background min={bg_min:.4f} mean={bg_mean:.4f})"
    )
    for i, j, d in sorted(accepted, key=lambda x: x[2]):
        log(f"  chr{i:02d} <-> chr{j:02d}: mean distance={d:.4f}")
    if unmatched:
        log(f"  no significant ancestral partner found for: {', '.join(f'chr{c:02d}' for c in unmatched)}")

    plot_ranked_distances(outdir, pair_dist, accepted)
    return [(i, j) for i, j, _ in accepted]


def plot_ranked_distances(outdir, pair_dist, accepted):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sorted_pairs = sorted(pair_dist.items(), key=lambda kv: kv[1])
    accepted_keys = {(min(i, j), max(i, j)) for i, j, _ in accepted}

    xs = list(range(len(sorted_pairs)))
    ys = [d for _, d in sorted_pairs]
    colors = ["crimson" if (min(k[0], k[1]), max(k[0], k[1])) in accepted_keys else "steelblue" for k, _ in sorted_pairs]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(xs, ys, c=colors, s=14)
    for idx, ((i, j), d) in enumerate(sorted_pairs):
        if (min(i, j), max(i, j)) in accepted_keys:
            ax.annotate(f"{i:02d}-{j:02d}", (idx, d), fontsize=6, xytext=(2, 4), textcoords="offset points")
    ax.set_xlabel("cross-chromosome pair rank (ascending mean distance)")
    ax.set_ylabel("mean whole-chromosome k-mer Jaccard distance")
    ax.set_title("Candidate ancestral (paleopolyploid) homeolog pairs -- red = accepted")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "homeolog_pairs.png"), dpi=150)
    plt.close(fig)
