#!/usr/bin/env python3
"""
An empirical answer to Twyford et al. 2025's call for a multidimensional
"poly-space" -- a PCA over multiple genomic-divergence axes, to see where
real polyploids actually cluster rather than assuming a binary. Builds a
species x metric feature matrix from data this project already computed
(meta/auto_allo_spectrum.tsv, meta/genome_partition.tsv, each species'
subgenomes/auto_allo_index.tsv) -- no new k-mer work.

Features (one row per species, all standardized before PCA):
- te_marker_fraction: mean, from auto_allo_index.tsv (repeat-content signal)
- partition_consistency: from auto_allo_spectrum.tsv (lineage-identity signal)
- distance_ratio_cv: from auto_allo_spectrum.tsv (event-count-consistency signal)
- pair_depth_cv: from auto_allo_spectrum.tsv (raw depth-consistency signal)
- mean_windowed_cv: from auto_allo_spectrum.tsv (within-pair patchiness)
- flip_rate: from auto_allo_spectrum.tsv (local switching rate)
- genome_partition_best_z: max |z_score| across all k in genome_partition.tsv
  (diffuse whole-genome block-structure strength; 0 if no significant row)

Missing values are mean-imputed per column (most metrics only apply to a
subset of species -- see INTERPRETATION.md for why) and this is a real
limitation of this exploratory analysis, not a hidden detail: a species
with mostly-imputed features sits near the origin by construction, not
because it's "intermediate" in any biological sense. n_features_present
is reported per species so this can be read correctly.

PCA is plain NumPy SVD on standardized (z-scored) features -- no sklearn
dependency, consistent with the rest of this repo's dependency-light style.

Usage: python3 scripts/poly_space_pca.py
Writes meta/poly_space_features.tsv and meta/poly_space_pca.png
"""
import csv
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(REPO, "results")
META_DIR = os.path.join(REPO, "meta")

FEATURE_COLS = [
    "te_marker_fraction",
    "partition_consistency",
    "distance_ratio_cv",
    "pair_depth_cv",
    "mean_windowed_cv",
    "flip_rate",
    "genome_partition_best_z",
]

# Hand-assigned read, from this project's own literature pass
# (meta/literature_ploidy.tsv) and INTERPRETATION.md -- for colouring the
# plot only, not used anywhere in the PCA computation itself.
CATEGORY = {
    "daGleHede1": "allo (confirmed)",
    "lpTriTurg1_A": "allo (confirmed)",
    "lpTriTurg1_B": "allo (confirmed)",
    "lpTriTurg1_AB": "allo (confirmed)",
    "drTriRepe1": "allo (confirmed)",
    "drSorDevo1": "allo (confirmed)",
    "drMyrSpic1": "allo (confirmed)",
    "SchCurv1": "auto (confirmed)",
    "SchYoun1": "auto (confirmed)",
    "drLytSali1": "auto (confirmed)",
    "ddEmpNigr1": "auto (lit, tension)",
    "ddHypMacu1": "auto (lit, tension)",
    "ddHesMatr1": "third category (segmental)",
    "llColAutu1": "third category (demi-dup)",
    "ddHypPerf1": "contested",
    "daInuConz1": "cryptic allo candidate",
    "dcCerAlpi1": "cryptic allo candidate",
    "dmRanRepe1": "cryptic allo candidate",
    "daBudDavi1": "allo (lit, weak data)",
}


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def get_te_marker_fraction(species):
    path = os.path.join(RESULTS_DIR, species, "subgenomes", "auto_allo_index.tsv")
    if not os.path.exists(path):
        return None
    vals = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            v = r.get("te_marker_fraction")
            if v:
                vals.append(float(v))
    if not vals:
        return None
    return sum(vals) / len(vals)


def get_genome_partition_best_z(rows_by_species, species):
    rows = rows_by_species.get(species, [])
    if not rows:
        return None
    return max(abs(float(r["z_score"])) for r in rows)


def main():
    spectrum_rows = {r["species"]: r for r in read_tsv(os.path.join(META_DIR, "auto_allo_spectrum.tsv"))}

    partition_rows = {}
    for r in read_tsv(os.path.join(META_DIR, "genome_partition.tsv")):
        partition_rows.setdefault(r["species"], []).append(r)

    species_list = sorted(spectrum_rows)
    raw = {}
    for sp in species_list:
        srow = spectrum_rows[sp]
        row = {
            "te_marker_fraction": get_te_marker_fraction(sp),
            "partition_consistency": float(srow["partition_consistency"]) if srow.get("partition_consistency") else None,
            "distance_ratio_cv": float(srow["distance_ratio_cv"]) if srow.get("distance_ratio_cv") else None,
            "pair_depth_cv": float(srow["pair_depth_cv"]) if srow.get("pair_depth_cv") else None,
            "mean_windowed_cv": float(srow["mean_windowed_cv"]) if srow.get("mean_windowed_cv") else None,
            "flip_rate": float(srow["flip_rate"]) if srow.get("flip_rate") else None,
            "genome_partition_best_z": get_genome_partition_best_z(partition_rows, sp),
        }
        raw[sp] = row

    # only keep species with at least 2 non-missing features -- below that,
    # a point is pure column-mean imputation and tells us nothing
    kept = [sp for sp in species_list if sum(v is not None for v in raw[sp].values()) >= 2]

    n_present = {sp: sum(v is not None for v in raw[sp].values()) for sp in kept}

    # mean-impute per column, over `kept` species only
    col_means = {}
    for col in FEATURE_COLS:
        vals = [raw[sp][col] for sp in kept if raw[sp][col] is not None]
        col_means[col] = sum(vals) / len(vals) if vals else 0.0

    X = np.array(
        [[raw[sp][col] if raw[sp][col] is not None else col_means[col] for col in FEATURE_COLS] for sp in kept]
    )

    # write the feature table actually used, imputed or not, for transparency
    with open(os.path.join(META_DIR, "poly_space_features.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["species", "n_features_present"] + FEATURE_COLS)
        for i, sp in enumerate(kept):
            w.writerow([sp, n_present[sp]] + [f"{x:.4f}" for x in X[i]])

    # standardize (z-score) each column, then PCA via SVD
    mu = X.mean(axis=0)
    sigma = X.std(axis=0)
    sigma[sigma == 0] = 1.0
    Xz = (X - mu) / sigma

    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    scores = U * S  # n_species x n_components
    explained = (S**2) / np.sum(S**2)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 9))
    categories = sorted(set(CATEGORY.get(sp, "data-only / unresolved") for sp in kept))
    cmap = plt.get_cmap("tab10")
    color_for = {cat: cmap(i % 10) for i, cat in enumerate(categories)}

    for i, sp in enumerate(kept):
        cat = CATEGORY.get(sp, "data-only / unresolved")
        marker = "o" if n_present[sp] >= 4 else "x"
        ax.scatter(scores[i, 0], scores[i, 1], color=color_for[cat], marker=marker, s=60)
        ax.annotate(sp, (scores[i, 0], scores[i, 1]), fontsize=7, xytext=(3, 3), textcoords="offset points")

    for cat, c in color_for.items():
        ax.scatter([], [], color=c, label=cat)
    ax.scatter([], [], color="gray", marker="x", label="<4 features present (mostly imputed)")
    ax.legend(fontsize=8, loc="best")

    ax.set_xlabel(f"PC1 ({explained[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1]*100:.1f}% variance)")
    ax.set_title(
        "Poly-space: PCA over inheritance-mode/divergence metrics\n"
        "(empirical answer to Twyford et al. 2025's call for a multidimensional\n"
        "polyploid landscape -- colours are this project's own prior reads, not PCA input)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(META_DIR, "poly_space_pca.png"), dpi=150)
    plt.close(fig)

    # loadings: which raw metrics drive PC1/PC2
    print("Explained variance:", [f"{e*100:.1f}%" for e in explained])
    print("\nPC1 loadings:")
    for col, load in sorted(zip(FEATURE_COLS, Vt[0]), key=lambda x: -abs(x[1])):
        print(f"  {col:28s} {load:+.3f}")
    print("\nPC2 loadings:")
    for col, load in sorted(zip(FEATURE_COLS, Vt[1]), key=lambda x: -abs(x[1])):
        print(f"  {col:28s} {load:+.3f}")
    print(f"\n{len(kept)} species included (of {len(species_list)} with any spectrum row)")


if __name__ == "__main__":
    main()
