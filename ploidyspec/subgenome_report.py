"""
Consolidates te_markers.py's per-pair output into a continuous auto<->allo
index and a subgenome-assignment summary. Pure aggregation of already-written
files (te_markers_summary.tsv, te_markers_windowed_*.tsv, windowed_*.tsv) --
no new FastK/Tabex calls, so this is fast even though te-markers itself isn't.
"""

import csv
import os
import statistics
from collections import Counter

from .common import log, subgenomes_dir, windowed_dir


def te_marker_fraction(n_markers_a, n_markers_b, n_highcopy_a, n_highcopy_b):
    """
    Fraction of each haplotype's high-copy (repetitive/TE-like) k-mer content
    that's subgenome-differential. Markers are a subset of high-copy k-mers
    by construction, so this is bounded [0,1]. Near 0: repeat content nearly
    identical between haplotype copies (autopolyploid-like, shared origin).
    Near 1: repeat content almost entirely non-overlapping (allopolyploid-
    like, independent parental origins). No calibrated threshold exists yet
    (see OUTPUTS.md) -- report the continuous value.
    """
    denom = n_highcopy_a + n_highcopy_b
    if denom <= 0:
        return 0.0
    return (n_markers_a + n_markers_b) / denom


def windowed_distance_cv(outdir, chrom_str, unit_a, unit_b):
    """
    Coefficient of variation (std/mean) of the raw-Jaccard windowed
    divergence track for this haplotype pair -- a secondary corroborating
    signal alongside the primary te-marker-based index. Uses windowed/'s
    output as-is (raw Jaccard, not Mash-corrected: windowed deliberately
    stayed uncorrected for cost reasons), so this is about relative
    variability/patchiness across the chromosome, not an absolute calibrated
    divergence estimate. None if the windowed track isn't available.
    """
    path = os.path.join(windowed_dir(outdir), f"windowed_{chrom_str}.tsv")
    if not os.path.exists(path):
        return None
    hap_a = unit_a.rsplit("_chr", 1)[0]
    hap_b = unit_b.rsplit("_chr", 1)[0]
    values = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if {row["hap_a"], row["hap_b"]} == {hap_a, hap_b}:
                values.append(float(row["jaccard_distance"]))
    if len(values) < 2:
        return None
    mean = statistics.mean(values)
    if mean == 0:
        return None
    return statistics.pstdev(values) / mean


def summarize_windowed_assignment(windowed_path, unit_a, unit_b):
    """
    Per-unit window-assignment counts/percentages from a
    te_markers_windowed_<a>x<b>.tsv, plus the list of "anomalous" windows
    where a haplotype's own sequence matched the OTHER haplotype's markers --
    candidate homeologous-exchange/introgression sites, not noise (see
    te_markers.py::paint_unit_windows).
    """
    self_side = {unit_a: "a", unit_b: "b"}
    counts = {unit_a: Counter(), unit_b: Counter()}
    anomalous = []
    with open(windowed_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            unit = row["unit"]
            assigned = row["assigned"]
            counts[unit][assigned] += 1
            self_ = self_side[unit]
            other = "b" if self_ == "a" else "a"
            if assigned == other:
                anomalous.append(
                    dict(
                        unit=unit,
                        win_start=row["win_start"],
                        win_end=row["win_end"],
                        matched_side=assigned,
                    )
                )

    rows = []
    for unit in (unit_a, unit_b):
        c = counts[unit]
        total = sum(c.values())
        self_ = self_side[unit]
        other = "b" if self_ == "a" else "a"
        rows.append(
            dict(
                unit=unit,
                total_windows=total,
                n_self=c.get(self_, 0),
                n_other=c.get(other, 0),
                n_ambiguous=c.get("ambiguous", 0),
                n_none=c.get("none", 0),
                pct_self=(c.get(self_, 0) / total * 100) if total else 0.0,
                pct_other=(c.get(other, 0) / total * 100) if total else 0.0,
            )
        )
    return rows, anomalous


def compute_subgenome_report(outdir):
    summary_path = os.path.join(subgenomes_dir(outdir), "te_markers_summary.tsv")
    if not os.path.exists(summary_path):
        log(
            f"no te_markers_summary.tsv in {outdir} -- run `te-markers` first, "
            f"skipping subgenome report"
        )
        return

    with open(summary_path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        log(
            "te_markers_summary.tsv has no pairs (species has <2 haplotype copies "
            "per chromosome number) -- nothing to report"
        )
        return
    if "n_highcopy_a" not in rows[0]:
        log(
            "te_markers_summary.tsv predates the n_highcopy_a/b columns -- "
            "rerun `te-markers` to get the auto/allo index"
        )
        return

    index_rows = []
    windows_summary_rows = []
    anomalous_rows = []
    for r in rows:
        ua, ub = r["unit_a"], r["unit_b"]
        n_a, n_b = int(r["n_markers_a"]), int(r["n_markers_b"])
        nh_a, nh_b = int(r["n_highcopy_a"]), int(r["n_highcopy_b"])
        fraction = te_marker_fraction(n_a, n_b, nh_a, nh_b)
        cv = windowed_distance_cv(outdir, r["chrom"], ua, ub)
        index_rows.append(
            dict(
                unit_a=ua,
                unit_b=ub,
                chrom=r["chrom"],
                n_markers_a=n_a,
                n_markers_b=n_b,
                n_highcopy_a=nh_a,
                n_highcopy_b=nh_b,
                te_marker_fraction=fraction,
                windowed_distance_cv=cv,
            )
        )

        windowed_path = os.path.join(
            subgenomes_dir(outdir), f"te_markers_windowed_{ua}x{ub}.tsv"
        )
        if os.path.exists(windowed_path):
            per_unit, anomalous = summarize_windowed_assignment(windowed_path, ua, ub)
            windows_summary_rows.extend(per_unit)
            anomalous_rows.extend(anomalous)

    write_index_tsv(outdir, index_rows)
    if windows_summary_rows:
        write_windows_summary_tsv(outdir, windows_summary_rows)
        write_anomalous_tsv(outdir, anomalous_rows)
        log(
            f"wrote auto_allo_index.tsv ({len(index_rows)} pairs), "
            f"subgenome_windows_summary.tsv, subgenome_anomalous_windows.tsv "
            f"({len(anomalous_rows)} anomalous windows)"
        )
    else:
        log(
            f"wrote auto_allo_index.tsv ({len(index_rows)} pairs) -- no "
            f"te-markers-windowed output found, skipping window-level reports"
        )


def write_index_tsv(outdir, rows):
    path = os.path.join(subgenomes_dir(outdir), "auto_allo_index.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "unit_a",
                "unit_b",
                "chrom",
                "n_markers_a",
                "n_markers_b",
                "n_highcopy_a",
                "n_highcopy_b",
                "te_marker_fraction",
                "windowed_distance_cv",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    r["unit_a"],
                    r["unit_b"],
                    r["chrom"],
                    r["n_markers_a"],
                    r["n_markers_b"],
                    r["n_highcopy_a"],
                    r["n_highcopy_b"],
                    f"{r['te_marker_fraction']:.6f}",
                    (
                        ""
                        if r["windowed_distance_cv"] is None
                        else f"{r['windowed_distance_cv']:.6f}"
                    ),
                ]
            )


def write_windows_summary_tsv(outdir, rows):
    path = os.path.join(subgenomes_dir(outdir), "subgenome_windows_summary.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "unit",
                "total_windows",
                "n_self",
                "n_other",
                "n_ambiguous",
                "n_none",
                "pct_self",
                "pct_other",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    r["unit"],
                    r["total_windows"],
                    r["n_self"],
                    r["n_other"],
                    r["n_ambiguous"],
                    r["n_none"],
                    f"{r['pct_self']:.2f}",
                    f"{r['pct_other']:.2f}",
                ]
            )


def write_anomalous_tsv(outdir, rows):
    path = os.path.join(subgenomes_dir(outdir), "subgenome_anomalous_windows.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["unit", "win_start", "win_end", "matched_side"])
        for r in rows:
            w.writerow([r["unit"], r["win_start"], r["win_end"], r["matched_side"]])
