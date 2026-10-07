"""
Cross-species tables built from each species' own outputs under a results
directory (one subdirectory per species). Species that haven't run the
`structure` stage get its outputs computed in memory here (no files written
under results/).

Writes to --outdir (core outputs, ROADMAP Phase 2 option B):
- panel_summary.tsv: one row per species -- copy number, copy divergence,
  ancient pairing and its synchrony, the strongest genome partition, and the
  rediploidization state counts.
- genome_partition.tsv: every species' significant diffuse partitions.
- rediploidization_panel.tsv: one row per species from rediploidization_summary.tsv.

and to --outdir/supplementary/:
- te_markers_panel.tsv: TE-marker readings. The within-genome contrast
  (te_split_median) is the usable one; absolute te_marker_fraction overlaps
  between a confirmed diploid (ddMalSylv1, up to 0.25 per chromosome) and the
  allo anchors, so it is reported for reference only.
- poly_space_features.tsv / poly_space_loadings.tsv / poly_space_pca.png: an
  exploratory PCA over the core numbers (Twyford et al. 2025's "poly-space").
  Missing values are mean-imputed, so species with few features sit near the
  origin by construction; n_features_present says how many were real.
"""

import csv
import os
import statistics
from collections import Counter

import numpy as np

from .common import log
from .structure import PARTITION_FIELDS, pair_synchrony, partition_rows

STATES = ("tetrasomic_like", "candidate", "partially_resolved", "resolved_lineages",
          "fusion_lineages", "one_divergent_copy", "not_assessable")

SUMMARY_FIELDS = (
    ["species", "n_chromosome_numbers", "copies_per_chromosome",
     "allelic_distance_median", "cross_chrom_distance_median",
     "n_accepted_pairs", "pair_depth_median", "pair_depth_cv",
     "partition_k2_z", "partition_best_k", "partition_best_z",
     "n_distinct_fusions"]
    + [f"n_{s}" for s in STATES]
)

FEATURE_COLS = [
    "allelic_distance_median",
    "paired_fraction",
    "pair_depth_cv",
    "partition_best_z",
    "resolved_fraction",
    "tetrasomic_fraction",
    "te_split_median",
]

TE_FIELDS = ["species", "n_chromosomes", "te_split_median", "te_marker_fraction_mean"]


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
    sync = read_tsv(os.path.join(sdir, "pair_synchrony.tsv"))
    parts = read_tsv(os.path.join(sdir, "genome_partition.tsv"))
    if sync is None or parts is None:
        return pair_synchrony(path, species), partition_rows(path, species)
    return sync[0], parts


def copy_divergence(path):
    """(median same-chromosome-number distance, median cross-number distance)
    from the matrix stage's pairs."""
    rows = read_tsv(os.path.join(path, "matrix", "whole_chrom_pairs.tsv")) or []
    same, cross = [], []
    for r in rows:
        ca, cb = r["unit_a"].rsplit("_chr", 1)[1], r["unit_b"].rsplit("_chr", 1)[1]
        (same if ca == cb else cross).append(float(r["distance"]))
    med = lambda v: statistics.median(v) if v else None  # noqa: E731
    return med(same), med(cross)


def modal_copies(path):
    rows = read_tsv(os.path.join(path, "matrix", "ploidy_summary.tsv")) or []
    counts = Counter(r["n_haplotype_copies"] for r in rows)
    return (counts.most_common(1)[0][0], len(rows)) if counts else ("", 0)


def te_markers_row(species, path):
    sdir = os.path.join(path, "subgenomes")
    index = read_tsv(os.path.join(sdir, "auto_allo_index.tsv"))
    if index is None:
        return None
    fracs = [float(r["te_marker_fraction"]) for r in index if r.get("te_marker_fraction")]
    lineage = read_tsv(os.path.join(sdir, "te_marker_fraction_by_lineage.tsv")) or []
    splits = [float(r["split_ratio"]) for r in lineage if r.get("split_ratio")]
    return {
        "species": species,
        "n_chromosomes": len({r["chrom"] for r in index}),
        "te_split_median": f"{statistics.median(splits):.3f}" if splits else "",
        "te_marker_fraction_mean": f"{statistics.mean(fracs):.4f}" if fracs else "",
    }


def rediploidization_row(species, path):
    rows = read_tsv(os.path.join(path, "rediploidization", "rediploidization_summary.tsv"))
    if rows is None:
        return None
    return {"species": species, **{r["metric"]: r["value"] for r in rows}}


def fmt(v, nd=4):
    return "" if v is None else f"{v:.{nd}f}"


def summary_row(species, path, sync, parts, redip):
    copies, n_chrom = modal_copies(path)
    same, cross = copy_divergence(path)
    k2 = [float(r["z_score"]) for r in parts if str(r["k"]) == "2"]
    best = max(parts, key=lambda r: float(r["z_score"]), default=None)
    row = {
        "species": species,
        "n_chromosome_numbers": n_chrom,
        "copies_per_chromosome": copies,
        "allelic_distance_median": fmt(same),
        "cross_chrom_distance_median": fmt(cross),
        "n_accepted_pairs": sync.get("n_accepted_pairs", ""),
        "pair_depth_median": sync.get("pair_depth_median", ""),
        "pair_depth_cv": sync.get("pair_depth_cv", ""),
        "partition_k2_z": fmt(k2[0], 2) if k2 else "",
        "partition_best_k": best["k"] if best else "",
        "partition_best_z": best["z_score"] if best else "",
        "n_distinct_fusions": (redip or {}).get("n_distinct_fusions", ""),
    }
    for s in STATES:
        row[f"n_{s}"] = (redip or {}).get(f"copy_state:{s}", "")
    return row


def features(summary, te_rows, redip):
    def num(v):
        return float(v) if v not in (None, "") else None

    out = {}
    for sp, row in summary.items():
        r = redip.get(sp, {})
        n = num(r.get("n_chromosome_numbers"))
        assessable = (n - (num(r.get("copy_state:not_assessable")) or 0)) if n else None

        def frac(*states):
            if not assessable:
                return None
            return sum(num(r.get(f"copy_state:{s}")) or 0 for s in states) / assessable

        paired = num(r.get("ancient_paired_chromosomes"))
        out[sp] = {
            "allelic_distance_median": num(row["allelic_distance_median"]),
            "paired_fraction": paired / n if paired is not None and n else None,
            "pair_depth_cv": num(row["pair_depth_cv"]),
            "partition_best_z": num(row["partition_best_z"]),
            "resolved_fraction": frac("resolved_lineages", "fusion_lineages"),
            "tetrasomic_fraction": frac("tetrasomic_like"),
            "te_split_median": num((te_rows.get(sp) or {}).get("te_split_median")),
        }
    return out


def poly_space(raw, categories, outdir):
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
        w.writerow(["feature", "PC1", "PC2"])
        for j, col in enumerate(FEATURE_COLS):
            w.writerow([col, f"{Vt[0, j]:.4f}", f"{Vt[1, j]:.4f}"])

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
    ax.set_title("Poly-space (exploratory): PCA over the core per-species numbers\n"
                 "(colours are prior reads from --categories, not PCA input)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "poly_space_pca.png"), dpi=150)
    plt.close(fig)
    log(f"poly-space PCA over {len(kept)} species: PC1 {explained[0] * 100:.1f}%, "
        f"PC2 {explained[1] * 100:.1f}%")


def build_panel(results_dir, outdir, categories_path=None):
    supp = os.path.join(outdir, "supplementary")
    os.makedirs(supp, exist_ok=True)
    summary, synchrony, partitions, redip, te_rows = {}, [], [], {}, {}
    for species, path in species_dirs(results_dir):
        sync, parts = species_structure(species, path)
        synchrony.append(sync)
        partitions.extend(parts)
        r = rediploidization_row(species, path)
        if r:
            redip[species] = r
        t = te_markers_row(species, path)
        if t:
            te_rows[species] = t
        summary[species] = summary_row(species, path, sync, parts, r)
    log(f"panel: {len(summary)} species from {results_dir}")

    write_tsv(os.path.join(outdir, "panel_summary.tsv"), SUMMARY_FIELDS,
              [summary[s] for s in sorted(summary)])
    write_tsv(os.path.join(outdir, "genome_partition.tsv"), PARTITION_FIELDS, partitions)
    if redip:
        rows = [redip[s] for s in sorted(redip)]
        fields = ["species"] + sorted({k for r in rows for k in r if k != "species"},
                                      key=lambda k: (not k.startswith("n_chrom"), k))
        write_tsv(os.path.join(outdir, "rediploidization_panel.tsv"), fields, rows)
    write_tsv(os.path.join(supp, "te_markers_panel.tsv"), TE_FIELDS,
              [te_rows[s] for s in sorted(te_rows)])

    categories = {}
    if categories_path:
        categories = {r["species"]: r["category"] for r in read_tsv(categories_path) or []}
    poly_space(features(summary, te_rows, redip), categories, supp)
    log(f"wrote panel_summary.tsv, genome_partition.tsv"
        f"{', rediploidization_panel.tsv' if redip else ''} in {outdir}; "
        f"te_markers_panel.tsv and poly_space_* in {supp}")
