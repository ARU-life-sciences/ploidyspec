"""
Consolidates te_markers.py's per-pair output into a continuous auto<->allo
index and a subgenome-assignment summary. Pure aggregation of already-written
files (te_markers_summary.tsv, te_markers_windowed_*.tsv, windowed_*.tsv) --
no new FastK/Tabex calls, so this is fast even though te-markers itself isn't.
"""

import csv
import os
import statistics
from collections import Counter, defaultdict

from .common import log, matrix_dir, subgenomes_dir, windowed_dir


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


def load_unit_lengths(outdir):
    """unit_id -> sequence length, from sequences.tsv (written by `prepare`)."""
    path = os.path.join(outdir, "sequences.tsv")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {r["unit_id"]: int(r["length"]) for r in csv.DictReader(f, delimiter="\t")}


def collect_unit_highcopy(index_rows):
    """unit_id -> n_highcopy, read off whichever auto_allo_index.tsv row it
    first appears in. This is a fixed property of a unit's own sequence at
    the marker k (confirmed empirically: a unit's n_highcopy is identical
    across every pairwise row it appears in), not something that depends on
    which partner it's compared against."""
    hc = {}
    for r in index_rows:
        hc.setdefault(r["unit_a"], r["n_highcopy_a"])
        hc.setdefault(r["unit_b"], r["n_highcopy_b"])
    return hc


def _majority_cluster(values, tolerance):
    """
    values: {unit: positive number}. Finds the largest set of units whose
    values are all mutually within `tolerance` (relative) of each other,
    returning (cluster_set, median_of_cluster) -- or (None, None) if no
    strict majority (more than half the units) agrees. Deliberately strict:
    a real chromosome fusion also produces a length/highcopy outlier by
    this measure (fused copies ~2x an unfused one), but splits a 4-unit
    group evenly (2-vs-2) -- neither side is a majority, so this returns
    (None, None) and flag_low_content_units correctly stays silent. Only a
    clear majority-vs-minority split (e.g. 3-vs-1) produces a baseline.
    """
    items = [(u, v) for u, v in values.items() if v > 0]
    best = []
    for _, v in items:
        cluster = [u2 for u2, v2 in items if abs(v2 - v) / v <= tolerance]
        if len(cluster) > len(best):
            best = cluster
    if not best or len(best) <= len(items) / 2:
        return None, None
    return set(best), statistics.median(values[u] for u in best)


def flag_low_content_units(units, lengths, highcopy, tolerance=0.15, threshold=0.75):
    """
    Flags units whose sequence length or total high-copy k-mer count sits
    well below a *majority-agreed* baseline among their same-chromosome
    siblings -- the signature of an incomplete/fragmented haplotype
    assembly, not real biological divergence. Requires a majority cluster
    (see _majority_cluster) before flagging anyone, specifically so a real
    2-vs-2 chromosome fusion or subgenome split (SchCurv1's chr19: two
    unfused copies at ~38-40Mb, two fused copies at ~69.5Mb, neither side a
    majority) is never mistaken for this. Confirmed against a real
    3-vs-1 case: daBudDavi1's chr05 HAP2 (32% short, a third of its three
    siblings' high-copy count, siblings mutually agree) is flagged; ddHypMacu1's
    HAP1 (570 excluded scaffold fragments vs. 2-3 for siblings) is the other
    confirmed case, both driving a spurious te_marker_fraction split with
    nothing to do with subgenome or fusion structure. A real fusion or
    subgenome split doesn't reduce a unit's own length or repeat content,
    only how much of it is shared with specific other copies -- so this is a
    different, complementary check to the lineage split itself, not a
    duplicate of it. Needs >=3 units with known length/highcopy; returns []
    otherwise.
    """
    lens = {u: lengths[u] for u in units if u in lengths}
    hcs = {u: highcopy[u] for u in units if u in highcopy}
    if len(lens) < 3 or len(hcs) < 3:
        return []

    len_cluster, len_baseline = _majority_cluster(lens, tolerance)
    hc_cluster, hc_baseline = _majority_cluster(hcs, tolerance)

    flagged = []
    for u in units:
        low_len = (
            len_cluster is not None and u not in len_cluster and u in lens
            and lens[u] / len_baseline < threshold
        )
        low_hc = (
            hc_cluster is not None and u not in hc_cluster and u in hcs
            and hcs[u] / hc_baseline < threshold
        )
        if low_len or low_hc:
            flagged.append(u)
    return flagged


def load_whole_chrom_distances(outdir):
    """unit-id-pair -> whole-chromosome Mash-corrected distance, from
    matrix/whole_chrom_distance_matrix.csv. Always numeric (unlike
    whole_chrom_pairs.tsv's distance column, which is blank when a pair is
    resolution_limited across every k)."""
    path = os.path.join(matrix_dir(outdir), "whole_chrom_distance_matrix.csv")
    if not os.path.exists(path):
        return {}
    dist = {}
    with open(path) as f:
        rows = list(csv.reader(f))
    header = rows[0][1:]
    for row in rows[1:]:
        uid = row[0]
        for j, val in enumerate(row[1:]):
            dist[(uid, header[j])] = float(val)
    return dist


def bipartition_by_distance(units, dist):
    """
    Split >=3 haplotype copies of one chromosome into two lineages using
    their whole-chromosome distance to each other: seed with the most-distant
    pair, then assign every other copy to whichever seed it's closer to. A
    dependency-free stand-in for 2-means on a distance matrix -- not meant to
    be a rigorous clustering method, just good enough to recover a real split
    when one exists (allopolyploid subgenomes, fusion-derived lineages). When
    no real split exists the two groups end up no more separated than any
    other bipartition would be, which shows up downstream as split_ratio
    near 1 -- read it as a diagnostic, not a classifier, same as
    distance_ratio in homeologs.py.
    """
    if len(units) < 3:
        return None
    best = None
    for i, u in enumerate(units):
        for v in units[i + 1 :]:
            d = dist.get((u, v))
            if d is None:
                continue
            if best is None or d > best[0]:
                best = (d, u, v)
    if best is None:
        return None
    _, seed_a, seed_b = best
    group_a, group_b = [seed_a], [seed_b]
    for u in units:
        if u in (seed_a, seed_b):
            continue
        da, db = dist.get((u, seed_a)), dist.get((u, seed_b))
        if da is None or db is None:
            continue
        (group_a if da < db else group_b).append(u)
    if not group_a or not group_b:
        return None
    return group_a, group_b


def compute_lineage_te_fractions(outdir, index_rows):
    """
    For each chromosome with >=3 haplotype copies, splits the copies into two
    lineages by whole-chromosome distance and reports mean te_marker_fraction
    within each lineage versus across them. Exists because the flat,
    unweighted genome-wide mean in auto_allo_index.tsv can look uninformative
    for a genome that's only partly resolved into two lineages -- confirmed
    case: SchCurv1's chr19 (one confirmed chromosome fusion, Xie et al. 2026)
    has within-lineage te_marker_fraction of 0.07-0.14 but 0.57-0.61 across
    lineages, a signal invisible in that chromosome's flat pairwise mean
    (~0.35) and easy to miss unless you already know the fusion structure.
    This makes the check automatic instead of requiring that prior knowledge.

    Also flags units that look like an incomplete/fragmented assembly rather
    than a real second lineage (see flag_low_content_units) -- confirmed
    case: daBudDavi1's chr05 split 6.13x, but as a lopsided 1-vs-3 (HAP2
    alone, 32% short and a third of its siblings' high-copy count) rather
    than a real 2-vs-2 fusion-like split. The split_ratio is still reported
    in that case, but flagged_units marks it as likely an assembly artifact
    rather than biology, without requiring the same manual length/marker
    digging every time.
    """
    own_dist = load_whole_chrom_distances(outdir)
    lengths = load_unit_lengths(outdir)
    highcopy = collect_unit_highcopy(index_rows)
    by_chrom = defaultdict(list)
    for r in index_rows:
        by_chrom[r["chrom"]].append(r)

    out_rows = []
    for chrom, pairs in sorted(by_chrom.items()):
        units = sorted(set(p["unit_a"] for p in pairs) | set(p["unit_b"] for p in pairs))
        split = bipartition_by_distance(units, own_dist)
        if split is None:
            continue
        group_a, group_b = split
        set_a = set(group_a)

        within_vals, cross_vals = [], []
        for p in pairs:
            same_group = (p["unit_a"] in set_a) == (p["unit_b"] in set_a)
            (within_vals if same_group else cross_vals).append(p["te_marker_fraction"])
        if not within_vals or not cross_vals:
            continue

        within_mean = sum(within_vals) / len(within_vals)
        cross_mean = sum(cross_vals) / len(cross_vals)
        flagged_units = flag_low_content_units(units, lengths, highcopy)
        out_rows.append(
            dict(
                chrom=chrom,
                group_a=",".join(group_a),
                group_b=",".join(group_b),
                n_within=len(within_vals),
                n_cross=len(cross_vals),
                within_te_marker_fraction=within_mean,
                cross_te_marker_fraction=cross_mean,
                split_ratio=(cross_mean / within_mean) if within_mean > 0 else None,
                flagged_units=",".join(flagged_units),
            )
        )
    return out_rows


def write_lineage_tsv(outdir, rows):
    path = os.path.join(subgenomes_dir(outdir), "te_marker_fraction_by_lineage.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "chrom",
                "group_a",
                "group_b",
                "n_within",
                "n_cross",
                "within_te_marker_fraction",
                "cross_te_marker_fraction",
                "split_ratio",
                "flagged_units",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    r["chrom"],
                    r["group_a"],
                    r["group_b"],
                    r["n_within"],
                    r["n_cross"],
                    f"{r['within_te_marker_fraction']:.6f}",
                    f"{r['cross_te_marker_fraction']:.6f}",
                    "" if r["split_ratio"] is None else f"{r['split_ratio']:.2f}",
                    r["flagged_units"],
                ]
            )


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

    lineage_rows = compute_lineage_te_fractions(outdir, index_rows)
    if lineage_rows:
        write_lineage_tsv(outdir, lineage_rows)
        log(
            f"wrote te_marker_fraction_by_lineage.tsv ({len(lineage_rows)} chromosomes "
            f"with >=3 copies split into two lineages)"
        )

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
