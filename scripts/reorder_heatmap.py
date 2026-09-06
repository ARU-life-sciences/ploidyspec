#!/usr/bin/env python3
"""
Reorders each species' chromosome-number-level distance matrix (from
homeolog_candidates_ranked.tsv -- every cross-chromosome pair already
tested, aggregated across haplotype copies) so that each accepted ancient-
homeolog pair's two members sit adjacent to each other, instead of the
pipeline's raw chromosome-number order. No new k-mer work.

Tests whether an apparent "checkerboard" pattern in the default-ordered
heatmap is a numbering-order artifact of already-known, strong,
individually-significant 1:1 pairs (it should resolve into a clean block-
diagonal pattern once reordered by known pair membership) or a genuine,
unresolved scatter (it would stay checkerboard-like even after reordering).

Pair order (which pair-block comes before which) is by ascending mean_
distance -- tightest/most significant pairs first -- purely for a readable
gradient, not a claim about anything biological. Unpaired chromosomes (for
species with incomplete pairing) are appended at the end, sorted by number.

Usage: python3 scripts/reorder_heatmap.py [species ...]
       (no args = every species with an accepted homeolog_pairs.tsv)
Writes results/<species>/homeologs/homeolog_pairs_reordered_heatmap.png
"""
import csv
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(REPO, "results")


def load_all_chroms(species_dir):
    path = os.path.join(species_dir, "matrix", "ploidy_summary.tsv")
    with open(path) as f:
        return [r["chrom"] for r in csv.DictReader(f, delimiter="\t")]


def load_distance_matrix(species_dir):
    path = os.path.join(species_dir, "homeologs", "homeolog_candidates_ranked.tsv")
    if not os.path.exists(path):
        return None
    dist = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            a, b = r["chrom_a"], r["chrom_b"]
            d = float(r["distance"])
            dist[(a, b)] = d
            dist[(b, a)] = d
    return dist


def load_accepted_pairs(species_dir):
    path = os.path.join(species_dir, "homeologs", "homeolog_pairs.tsv")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [
            (r["chrom_a"], r["chrom_b"], float(r["mean_distance"]))
            for r in csv.DictReader(f, delimiter="\t")
        ]


def build_pair_order(all_chroms, pairs):
    pairs_sorted = sorted(pairs, key=lambda p: p[2])
    order = []
    paired = set()
    for a, b, _ in pairs_sorted:
        order.append(a)
        order.append(b)
        paired.add(a)
        paired.add(b)
    unpaired = sorted((c for c in all_chroms if c not in paired), key=lambda x: int(x[3:]))
    order.extend(unpaired)
    return order, len(pairs_sorted)


def plot_reordered(species, species_dir, order, n_pairs, dist):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    n = len(order)
    mat = np.zeros((n, n))
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            if i == j:
                mat[i, j] = np.nan
            else:
                mat[i, j] = dist.get((a, b), np.nan)

    off_diag = mat[~np.eye(n, dtype=bool)]
    off_diag = off_diag[~np.isnan(off_diag)]
    if off_diag.size == 0:
        return False
    vmin, vmax = np.percentile(off_diag, [2, 98])

    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("white")

    fig_w = max(6, n * 0.28)
    fig, ax = plt.subplots(figsize=(fig_w, fig_w))
    im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(order, rotation=90, fontsize=max(4, 9 - n // 15))
    ax.set_yticklabels(order, fontsize=max(4, 9 - n // 15))

    # mark pair-block boundaries every 2 chromosomes, up to n_pairs*2
    for k in range(1, n_pairs):
        ax.axhline(2 * k - 0.5, color="white", linewidth=0.6, alpha=0.6)
        ax.axvline(2 * k - 0.5, color="white", linewidth=0.6, alpha=0.6)
    if n_pairs * 2 < n:
        ax.axhline(n_pairs * 2 - 0.5, color="red", linewidth=1.0, alpha=0.8)
        ax.axvline(n_pairs * 2 - 0.5, color="red", linewidth=1.0, alpha=0.8)

    ax.set_title(
        f"{species}: chromosome-number distance, reordered by accepted "
        f"homeolog pair membership\n({n_pairs} accepted pairs, sorted "
        f"tightest-first; red line = start of unpaired chromosomes, if any)\n"
        f"2nd-98th pct: {vmin:.3f}-{vmax:.3f}",
        fontsize=8,
        wrap=True,
    )
    fig.colorbar(im, ax=ax, shrink=0.7, label="distance")
    fig.tight_layout()
    out_path = os.path.join(
        species_dir, "homeologs", "homeolog_pairs_reordered_heatmap.png"
    )
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


def main():
    species_args = sys.argv[1:]
    if species_args:
        species_list = species_args
    else:
        species_list = sorted(
            d for d in os.listdir(RESULTS_DIR) if os.path.isdir(os.path.join(RESULTS_DIR, d))
        )

    for species in species_list:
        species_dir = os.path.join(RESULTS_DIR, species)
        pairs = load_accepted_pairs(species_dir)
        if not pairs:
            continue
        dist = load_distance_matrix(species_dir)
        if dist is None:
            continue
        all_chroms = load_all_chroms(species_dir)
        order, n_pairs = build_pair_order(all_chroms, pairs)
        ok = plot_reordered(species, species_dir, order, n_pairs, dist)
        print(f"{species}: {'wrote' if ok else 'skipped'} reordered heatmap "
              f"({n_pairs} pairs, {len(order)} chromosomes)")


if __name__ == "__main__":
    main()
