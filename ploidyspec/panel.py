"""
Cross-species tables and the poly-space PCA, built from each species' own
outputs under a results directory (one subdirectory per species). Replaces the
panel-wide scripts (auto_allo_spectrum.py, genome_partition.py,
poly_space_pca.py). Species that haven't run the `structure` stage get their
structure metrics computed in memory here (no files written under results/).

Writes to --outdir:
- auto_allo_spectrum.tsv: one row of structure/inheritance_metrics.tsv per species
- genome_partition.tsv: every species' significant diffuse partitions
- rediploidization_panel.tsv: one row per species from rediploidization_summary.tsv
- poly_space_features.tsv / poly_space_loadings.tsv / poly_space_pca.png: an
  exploratory PCA over seven metrics (Twyford et al. 2025's "poly-space").
  Missing values are mean-imputed, so species with few features sit near the
  origin by construction; n_features_present says how many were real.
"""

import csv
import os

import numpy as np

from .common import log
from .structure import (
    INHERITANCE_FIELDS,
    PARTITION_FIELDS,
    inheritance_metrics,
    partition_rows,
)

FEATURE_COLS = [
    "te_marker_fraction",
    "partition_consistency",
    "distance_ratio_cv",
    "pair_depth_cv",
    "mean_windowed_cv",
    "flip_rate",
    "genome_partition_best_z",
]


def read_tsv(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return list(csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t"))


def write_tsv(path, fields, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def species_dirs(results_dir):
    out = []
    for d in sorted(os.listdir(results_dir)):
        path = os.path.join(results_dir, d)
        if os.path.exists(os.path.join(path, "matrix", "whole_chrom_distance_matrix.csv")):
            out.append((d, path))
    return out


def species_structure(species, path):
    sdir = os.path.join(path, "structure")
    inh = read_tsv(os.path.join(sdir, "inheritance_metrics.tsv"))
    parts = read_tsv(os.path.join(sdir, "genome_partition.tsv"))
    if inh is None or parts is None:
        return inheritance_metrics(path, species), partition_rows(path, species)
    return inh[0], parts


def mean_te_marker_fraction(path):
    rows = read_tsv(os.path.join(path, "subgenomes", "auto_allo_index.tsv")) or []
    vals = [float(r["te_marker_fraction"]) for r in rows if r.get("te_marker_fraction")]
    return sum(vals) / len(vals) if vals else None


def rediploidization_row(species, path):
    rows = read_tsv(os.path.join(path, "rediploidization", "rediploidization_summary.tsv"))
    if rows is None:
        return None
    return {"species": species, **{r["metric"]: r["value"] for r in rows}}


def poly_space(spectrum, partitions, te_frac, categories, outdir):
    best_z = {}
    for r in partitions:
        best_z[r["species"]] = max(best_z.get(r["species"], 0.0), abs(float(r["z_score"])))

    def num(v):
        return float(v) if v not in (None, "") else None

    raw = {}
    for sp, row in spectrum.items():
        raw[sp] = {
            "te_marker_fraction": te_frac.get(sp),
            "partition_consistency": num(row.get("partition_consistency")),
            "distance_ratio_cv": num(row.get("distance_ratio_cv")),
            "pair_depth_cv": num(row.get("pair_depth_cv")),
            "mean_windowed_cv": num(row.get("mean_windowed_cv")),
            "flip_rate": num(row.get("flip_rate")),
            "genome_partition_best_z": best_z.get(sp),
        }
    # below 2 real features a point is pure column-mean imputation
    kept = [sp for sp in sorted(raw) if sum(v is not None for v in raw[sp].values()) >= 2]
    if len(kept) < 3:
        log(f"poly-space: only {len(kept)} species with >= 2 features -- skipping PCA")
        return
    n_present = {sp: sum(v is not None for v in raw[sp].values()) for sp in kept}
    col_means = {}
    for col in FEATURE_COLS:
        vals = [raw[sp][col] for sp in kept if raw[sp][col] is not None]
        col_means[col] = sum(vals) / len(vals) if vals else 0.0
    X = np.array([[raw[sp][c] if raw[sp][c] is not None else col_means[c] for c in FEATURE_COLS]
                  for sp in kept])

    with open(os.path.join(outdir, "poly_space_features.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["species", "n_features_present"] + FEATURE_COLS)
        for i, sp in enumerate(kept):
            w.writerow([sp, n_present[sp]] + [f"{x:.4f}" for x in X[i]])

    sigma = X.std(axis=0)
    sigma[sigma == 0] = 1.0
    Xz = (X - X.mean(axis=0)) / sigma
    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    scores = U * S
    explained = (S**2) / np.sum(S**2)

    with open(os.path.join(outdir, "poly_space_loadings.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["component", "explained_variance"] + FEATURE_COLS)
        for c in range(min(3, len(S))):
            w.writerow([f"PC{c + 1}", f"{explained[c]:.4f}"] + [f"{v:+.3f}" for v in Vt[c]])

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 9))
    cats = sorted({categories.get(sp, "data-only / unresolved") for sp in kept})
    cmap = plt.get_cmap("tab10")
    colour = {c: cmap(i % 10) for i, c in enumerate(cats)}
    for i, sp in enumerate(kept):
        cat = categories.get(sp, "data-only / unresolved")
        ax.scatter(scores[i, 0], scores[i, 1], color=colour[cat],
                   marker="o" if n_present[sp] >= 4 else "x", s=60)
        ax.annotate(sp, (scores[i, 0], scores[i, 1]), fontsize=7, xytext=(3, 3),
                    textcoords="offset points")
    for c, col in colour.items():
        ax.scatter([], [], color=col, label=c)
    ax.scatter([], [], color="gray", marker="x", label="<4 features present (mostly imputed)")
    ax.legend(fontsize=8, loc="best")
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}% variance)")
    ax.set_title("Poly-space: PCA over inheritance-mode/divergence metrics\n"
                 "(colours are prior reads from --categories, not PCA input)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "poly_space_pca.png"), dpi=150)
    plt.close(fig)
    log(f"poly-space PCA over {len(kept)} species: PC1 {explained[0] * 100:.1f}%, "
        f"PC2 {explained[1] * 100:.1f}%")


def build_panel(results_dir, outdir, categories_path=None):
    os.makedirs(outdir, exist_ok=True)
    spectrum, partitions, redip, te_frac = {}, [], [], {}
    for species, path in species_dirs(results_dir):
        inh, parts = species_structure(species, path)
        spectrum[species] = inh
        partitions.extend(parts)
        te = mean_te_marker_fraction(path)
        if te is not None:
            te_frac[species] = te
        row = rediploidization_row(species, path)
        if row:
            redip.append(row)
    log(f"panel: {len(spectrum)} species from {results_dir}")

    write_tsv(os.path.join(outdir, "auto_allo_spectrum.tsv"), INHERITANCE_FIELDS,
              [spectrum[s] for s in sorted(spectrum)])
    write_tsv(os.path.join(outdir, "genome_partition.tsv"), PARTITION_FIELDS, partitions)
    if redip:
        fields = ["species"] + sorted({k for r in redip for k in r if k != "species"},
                                      key=lambda k: (not k.startswith("n_chrom"), k))
        write_tsv(os.path.join(outdir, "rediploidization_panel.tsv"), fields, redip)

    categories = {}
    if categories_path:
        categories = {r["species"]: r["category"] for r in read_tsv(categories_path) or []}
    poly_space(spectrum, partitions, te_frac, categories, outdir)
    log(f"wrote auto_allo_spectrum.tsv, genome_partition.tsv"
        f"{', rediploidization_panel.tsv' if redip else ''} and poly_space_* in {outdir}")
