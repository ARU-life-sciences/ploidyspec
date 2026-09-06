#!/usr/bin/env python3
"""
Reorders each species' chromosome-number-level distance matrix (from
homeolog_candidates_ranked.tsv -- every cross-chromosome pair already
tested, aggregated across haplotype copies) by proper hierarchical
(average-linkage) clustering, reusing the exact same `average_linkage_
order` function the pipeline's own whole_chrom_distance_heatmap already
uses at the unit level. No new k-mer work.

v2: the first version of this script sorted the FDR-accepted pairs
(homeolog_pairs.tsv) as a flat, independent list, ordered only by each
pair's own distance. That's wrong whenever pairs relate to EACH OTHER --
it destroyed real multi-chromosome block structure (llColAutu1's uneven
3-6-chromosome demi-duplication groups got scattered into independent
pairs) and gave no way to see nested structure (ddLepDrab1's pairs-of-
pairs). Hierarchical clustering fixes both: it nests pairs inside quartets
inside larger related groups automatically, at whatever level the data
actually supports, with no k to hand-pick. FDR-accepted pairs are still
marked (thick border) so you can see which adjacencies are individually
significant vs. which come from the coarser, unsupervised clustering.

Usage: python3 scripts/reorder_heatmap.py [species ...]
       (no args = every species with a homeolog_candidates_ranked.tsv)
Writes results/<species>/homeologs/homeolog_pairs_reordered_heatmap.png
"""
import csv
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(REPO, "results")
sys.path.insert(0, REPO)

from ploidyspec.whole_matrix import average_linkage_order  # noqa: E402


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


def load_accepted_pairs(species_dir):
    path = os.path.join(species_dir, "homeologs", "homeolog_pairs.tsv")
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        return {
            frozenset([r["chrom_a"], r["chrom_b"]])
            for r in csv.DictReader(f, delimiter="\t")
        }


def plot_reordered(species, species_dir, chroms, dist, accepted_pairs):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    n = len(chroms)
    mat = np.zeros((n, n))
    for i, a in enumerate(chroms):
        for j, b in enumerate(chroms):
            mat[i, j] = 0.0 if i == j else dist.get((a, b), np.nan)

    order = average_linkage_order(mat)
    ordered_chroms = [chroms[i] for i in order]
    ordered_mat = mat[np.ix_(order, order)].copy()

    off_diag = ordered_mat[~np.eye(n, dtype=bool)]
    off_diag = off_diag[~np.isnan(off_diag)]
    if off_diag.size == 0:
        return False
    vmin, vmax = np.percentile(off_diag, [2, 98])
    np.fill_diagonal(ordered_mat, np.nan)

    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("white")

    fig_w = max(6, n * 0.28)
    fig, ax = plt.subplots(figsize=(fig_w, fig_w))
    ax.imshow(ordered_mat, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(ordered_chroms, rotation=90, fontsize=max(4, 9 - n // 15))
    ax.set_yticklabels(ordered_chroms, fontsize=max(4, 9 - n // 15))

    # outline cells that are individually FDR-accepted pairs, so you can see
    # which adjacencies are statistically confirmed vs. just clustering-implied
    for i, a in enumerate(ordered_chroms):
        for j, b in enumerate(ordered_chroms):
            if i != j and frozenset([a, b]) in accepted_pairs:
                ax.add_patch(
                    plt.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1, fill=False,
                        edgecolor="red", linewidth=1.4,
                    )
                )

    im = ax.images[0]
    fig.colorbar(im, ax=ax, shrink=0.7, label="distance")
    ax.set_title(
        f"{species}: chromosome-number distance, hierarchically reordered "
        f"(average-linkage)\nred outline = individually FDR-accepted pair "
        f"({len(accepted_pairs)} total)\n2nd-98th pct: {vmin:.3f}-{vmax:.3f}",
        fontsize=8,
        wrap=True,
    )
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
        dist, chroms = load_distance_matrix(species_dir)
        if dist is None or not chroms or len(chroms) < 3:
            continue
        accepted_pairs = load_accepted_pairs(species_dir)
        ok = plot_reordered(species, species_dir, chroms, dist, accepted_pairs)
        print(f"{species}: {'wrote' if ok else 'skipped'} reordered heatmap "
              f"({len(accepted_pairs)} accepted pairs, {len(chroms)} chromosomes)")


if __name__ == "__main__":
    main()
