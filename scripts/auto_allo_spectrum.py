#!/usr/bin/env python3
"""
Three inheritance-mode-based metrics for the auto<->allo spectrum, all
derived from data the pipeline already computes (no new FastK/k-mer work):

1. partition_consistency: does the same lineage-partition of >=3 haplotype
   copies recur across every chromosome (disomic/allo signature), or does
   the "odd one out" rotate between different haplotype labels on different
   chromosomes (tetrasomic/auto signature)? Reuses
   subgenome_report.bipartition_by_distance on matrix/whole_chrom_distance_
   matrix.csv, generalizing the by-hand check that caught drLytSali1's
   rotating singleton.

2. mean_windowed_cv: average coefficient-of-variation of the raw windowed
   jaccard-distance track across every accepted homeolog pair's windowed
   file (homeologs/windowed_chrAAxBB.tsv). Low = uniform divergence along
   the whole chromosome pair (one clean historical split, allo-like). High
   = patchy/mosaic (frequent local exchange, auto-like or introgression).

3. pair_depth_cv: coefficient-of-variation of homeolog_pairs.tsv's
   mean_distance across all of a species' accepted pairs. Low = every pair
   diverged to about the same depth (single historical event, allo-like).
   High = wide spread (auto-like/complex history).

4. mean_run_length_windows / flip_rate: the DIRECT run-length/switching-rate
   test (as opposed to mean_windowed_cv's magnitude-patchiness proxy) --
   for each chromosome's whole-chromosome partition (same grouping as
   metric 1), walks the within-chromosome windowed track and asks, window
   by window, whether the local self-vs-cross distance ordering still
   agrees with the global partition. Long runs/low flip rate = the split
   holds up almost everywhere (rare exchange, allo-like); short runs/high
   flip rate = the local pattern keeps flipping (routine recombination,
   auto-like). EMPIRICALLY, across the 12 species this is computable for,
   this does NOT cleanly separate confirmed autos from allo-leaning
   species -- SchCurv1 (confirmed auto, Xie et al. 2026) and SchYoun1
   (also confirmed auto) sit at opposite ends of the observed range
   (4.12 vs 5.31 mean run length), so whatever this signal mostly reflects
   at 250kb window resolution (likely per-species noise level: genome
   size, repeat content, window count), it isn't primarily inheritance
   mode. The one genuine standout is ddHesMatr1 (2.25, flip_rate 0.44,
   roughly double any other species) -- consistent with its literature
   description as a "segmental allotetraploid" (mixed bivalent/quadrivalent
   meiotic pairing, i.e. a real patchwork of disomic- and tetrasomic-like
   regions, which is exactly what very frequent local flipping would look
   like). Deliberately NOT folded into combined_allo_score given this null
   result on the panel as a whole -- report as informational only.

5. distance_ratio_cv: coefficient of variation of homeolog_pairs' `distance_
   ratio` (homeolog-pair distance / mean of both chromosomes' own within-
   chromosome distance, from homeologs/ploidy_ancestry_summary.tsv) across
   a species' accepted pairs -- the actual "cross-chromosome divergence-
   depth consistency" test, distinct from pair_depth_cv's raw mean_distance
   version. Low = every pair diverged to the same depth relative to its own
   baseline (one clean historical event, allo-like); high = wide,
   inconsistent spread (messier/multi-event history, auto-like). Confirms
   the daLatSqua1 (CV 1.01) / llColAutu1 (CV 0.59) outlier status already
   flagged ad hoc as `distance_ratio` spread in INTERPRETATION.md, now a
   formal statistic. Also lands very low (0.08) for drLytSali1 -- the SAME
   tetrasomic-homogenization confound as metrics 2/3 (multivalent
   recombination equalizes divergence across the whole genome, not just
   within one pair), so this is deliberately NOT folded into
   combined_allo_score either, for the same reason.

Metrics 1-3 are each directionally mapped to an "allo-likeness" score in
(0,1] (1/(1+x) for the two CV-based metrics, direct value for partition
consistency), then averaged over whichever are actually available for that
species into a single heuristic combined_score. This combined score is
explicitly NOT a classifier -- report it alongside the raw metrics and
te_marker_fraction, never in place of them.

Usage: python3 scripts/auto_allo_spectrum.py > meta/auto_allo_spectrum.tsv
"""
import csv
import os
import re
import statistics
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from ploidyspec.common import homeologs_dir, matrix_dir, windowed_dir  # noqa: E402
from ploidyspec.kmer_tables import load_sequences  # noqa: E402
from ploidyspec.subgenome_report import (  # noqa: E402
    bipartition_by_distance,
    load_whole_chrom_distances,
)

RESULTS_DIR = os.path.join(REPO, "results")


def hap_of(unit_id):
    return unit_id.rsplit("_chr", 1)[0]


def canonical_partition(group_a, group_b):
    """Canonical, order-independent representation of a 2-way hap-label
    split, so the same partition compares equal regardless of which side
    bipartition_by_distance happened to call 'a' vs 'b'."""
    a = frozenset(hap_of(u) for u in group_a)
    b = frozenset(hap_of(u) for u in group_b)
    return frozenset([a, b])


def metric_partition_consistency(species_dir):
    """Returns (consistency, n_chroms_split, modal_singleton_hap). The third
    value is set only when the modal partition is a strict 1-vs-rest split,
    naming the recurring singleton haplotype -- a high consistency score
    driven by one haplotype label that's a KNOWN fragmented/low-quality
    assembly (check unplaced.tsv) is an assembly-quality artifact, not a
    real disomic/allo signature. Confirmed case: ddHypMacu1's HAP1 (570
    unplaced scaffolds vs 2-3 for siblings) is the singleton on all 8/8
    chromosomes, giving a spurious consistency of 1.0."""
    seq_tsv = os.path.join(species_dir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        return None, 0, None
    units = load_sequences(seq_tsv)
    by_chrom = {}
    for u in units:
        by_chrom.setdefault(int(u["chrom"]), []).append(u["unit_id"])

    dist = load_whole_chrom_distances(species_dir)
    if not dist:
        return None, 0, None

    partitions = []
    for chrom, unit_ids in by_chrom.items():
        if len(unit_ids) < 3:
            continue
        split = bipartition_by_distance(unit_ids, dist)
        if split is None:
            continue
        group_a, group_b = split
        partitions.append(canonical_partition(group_a, group_b))

    n = len(partitions)
    if n < 2:
        return None, n, None

    counts = {}
    for p in partitions:
        counts[p] = counts.get(p, 0) + 1
    modal_partition = max(counts, key=counts.get)
    modal_count = counts[modal_partition]

    singleton_hap = None
    sides = list(modal_partition)
    if len(sides) == 2:
        for side in sides:
            if len(side) == 1:
                singleton_hap = next(iter(side))

    return modal_count / n, n, singleton_hap


def metric_run_length(species_dir):
    """
    The actual run-length/switching-rate test, as opposed to metric_windowed_cv's
    magnitude-patchiness proxy: for each chromosome with a valid >=3-copy
    partition (group_a, group_b from bipartition_by_distance, based on
    WHOLE-chromosome distance), walk the within-chromosome windowed track
    (windowed/windowed_chrNN.tsv) and ask, window by window, whether the
    LOCAL pattern still agrees with the GLOBAL partition -- is the mean
    within-group ("self") distance still lower than the mean between-group
    ("cross") distance in this specific window, or has it locally inverted
    (a candidate local-exchange/recombination event)?

    Returns (mean_run_length_windows, flip_rate, n_windows_total). A real
    disomic/allo split should hold up almost everywhere (long runs, low
    flip rate -- rare exchange between subgenomes that otherwise stay
    separate). A tetrasomic/auto split should show the local pattern
    flipping often (short runs, high flip rate -- routine multivalent
    recombination continually reshuffling which copies look alike).
    """
    seq_tsv = os.path.join(species_dir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        return None, None, 0
    units = load_sequences(seq_tsv)
    by_chrom = {}
    for u in units:
        by_chrom.setdefault(int(u["chrom"]), []).append(u["unit_id"])

    dist = load_whole_chrom_distances(species_dir)
    if not dist:
        return None, None, 0

    all_run_lengths = []
    all_flips = 0
    all_windows = 0

    for chrom, unit_ids in by_chrom.items():
        if len(unit_ids) < 3:
            continue
        split = bipartition_by_distance(unit_ids, dist)
        if split is None:
            continue
        group_a, group_b = split
        set_a, set_b = set(group_a), set(group_b)

        self_pairs = set()
        cross_pairs = set()
        for g in (set_a, set_b):
            for i, u in enumerate(sorted(g)):
                for v in sorted(g)[i + 1 :]:
                    self_pairs.add(frozenset([u, v]))
        for u in set_a:
            for v in set_b:
                cross_pairs.add(frozenset([u, v]))
        if not self_pairs or not cross_pairs:
            continue

        chrom_str = f"chr{chrom:02d}"
        path = os.path.join(windowed_dir(species_dir), f"windowed_{chrom_str}.tsv")
        if not os.path.exists(path):
            continue

        by_window = {}
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                key = frozenset([r["unit_a"], r["unit_b"]])
                win = int(r["win_start"])
                by_window.setdefault(win, {})[key] = float(r["jaccard_distance"])

        agreements = []
        for win in sorted(by_window):
            vals = by_window[win]
            self_vals = [vals[p] for p in self_pairs if p in vals]
            cross_vals = [vals[p] for p in cross_pairs if p in vals]
            if not self_vals or not cross_vals:
                continue
            agreements.append(statistics.mean(self_vals) < statistics.mean(cross_vals))

        if len(agreements) < 4:
            continue

        run_lengths = []
        current = 1
        flips = 0
        for i in range(1, len(agreements)):
            if agreements[i] == agreements[i - 1]:
                current += 1
            else:
                run_lengths.append(current)
                current = 1
                flips += 1
        run_lengths.append(current)

        all_run_lengths.extend(run_lengths)
        all_flips += flips
        all_windows += len(agreements)

    if not all_run_lengths or all_windows == 0:
        return None, None, 0
    return statistics.mean(all_run_lengths), all_flips / all_windows, all_windows


def metric_windowed_cv(species_dir):
    pairs_tsv = os.path.join(homeologs_dir(species_dir), "homeolog_pairs.tsv")
    if not os.path.exists(pairs_tsv):
        return None, 0
    with open(pairs_tsv) as f:
        pairs = list(csv.DictReader(f, delimiter="\t"))
    if not pairs:
        return None, 0

    cvs = []
    for row in pairs:
        a = row["chrom_a"].replace("chr", "").zfill(2)
        b = row["chrom_b"].replace("chr", "").zfill(2)
        label = f"chr{a}x{b}"
        path = os.path.join(homeologs_dir(species_dir), f"windowed_{label}.tsv")
        if not os.path.exists(path):
            continue
        by_pair = {}
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                key = frozenset([r["hap_a"], r["hap_b"]])
                by_pair.setdefault(key, []).append(float(r["jaccard_distance"]))
        for key, values in by_pair.items():
            if len(values) < 3:
                continue
            mean = statistics.mean(values)
            if mean == 0:
                continue
            cvs.append(statistics.pstdev(values) / mean)

    if not cvs:
        return None, 0
    return statistics.mean(cvs), len(cvs)


_HAP_TAG_RE = re.compile(r"(?:SUPER_\d+_)?HAP(\d+)(?:_SUPER|[_x])?")


def check_singleton_artifact(species_dir, singleton_hap):
    """Cross-references the modal partition's recurring singleton haplotype
    against unplaced.tsv fragment counts. Two different naming conventions
    in this panel mean a hap label can't be read off unplaced.tsv one way:
    AUTO-labeled (curated-stage) species carry a parseable HAP tag directly
    in seq_id (e.g. "SUPER_2_HAP1_unloc_1") even when several haps share one
    source fasta (ddHesMatr1's grab-bag files); positionally-labeled
    (release-stage) species like ddHypMacu1 have plain INSDC seq_ids with no
    HAP tag at all, so the source fasta path (one-to-one with hap there) is
    the only way to know. Tries the seq_id regex first, falls back to the
    source-path mapping from sequences.tsv's placed units.

    Compares the singleton's count against the MEDIAN (not mean) of the
    other haps -- robust to a real, unrelated outlier among the *other*
    haps corrupting the comparison for a singleton that isn't actually
    unusual itself. Confirmed case: SchCurv1's modal singleton is HAP1
    (count 0, same as HAP2/HAP3), but HAP4 alone carries 2006 unplaced
    fragments (real fusion-boundary complexity, not an artifact of HAP1) --
    a mean-based comparison falsely flagged HAP1 by getting dragged off by
    HAP4; the median comparison correctly clears it.

    A haplotype that's a genuine outlier (either much MORE fragmented, like
    ddHypMacu1's HAP1 with 570 unplaced scaffolds vs 2-3 for siblings, or
    much LESS fragmented/uniquely complete, like a species where only
    HAP1/HAP2 get unlocalized-sequence placement so HAP3+ has ~0 unplaced by
    construction -- see OUTPUTS.md) is a candidate assembly-composition
    artifact driving partition_consistency, not necessarily a real lineage
    signal -- but this is a first-pass screen, not a verdict: it flags an
    asymmetry exists, not that the asymmetry actually explains the
    partition (drAriEdul1's HAP4 asymmetry is flagged here too, but a prior,
    more careful investigation checking whether it correlates with the
    actual divergence signal ruled it out as an artifact -- see
    INTERPRETATION.md). Returns True/False/None (None = couldn't check)."""
    if not singleton_hap:
        return None
    seq_tsv = os.path.join(species_dir, "sequences.tsv")
    unplaced_tsv = os.path.join(species_dir, "unplaced.tsv")
    if not os.path.exists(seq_tsv) or not os.path.exists(unplaced_tsv):
        return None

    source_to_hap = {}
    all_haps = set()
    with open(seq_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            source_to_hap.setdefault(r["source"], r["hap"])
            all_haps.add(r["hap"])

    if len(all_haps) < 3:
        return None

    counts = {h: 0 for h in all_haps}
    with open(unplaced_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            m = _HAP_TAG_RE.search(r["seq_id"])
            hap = f"HAP{m.group(1)}" if m else source_to_hap.get(r["source"])
            if hap in counts:
                counts[hap] += 1

    singleton_count = counts.get(singleton_hap, 0)
    others = [c for h, c in counts.items() if h != singleton_hap]
    if not others:
        return None
    others_median = statistics.median(others)
    if others_median == 0 and singleton_count == 0:
        return None
    ratio = (singleton_count + 1) / (others_median + 1)
    return ratio >= 3.0 or ratio <= (1.0 / 3.0)


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


def metric_distance_ratio_cv(species_dir):
    """
    The actual "cross-chromosome divergence-depth consistency" test, as
    distinct from metric_pair_depth_cv above: reads homeologs/
    ploidy_ancestry_summary.tsv's `distance_ratio` column (homeolog-pair
    distance divided by the mean of both chromosomes' own within-chromosome
    distance -- already used ad hoc in INTERPRETATION.md to flag daLatSqua1's
    20x and llColAutu1's 37x spreads), dedupes to one value per accepted
    pair (the file lists each pair twice, once per chromosome side), and
    reports its coefficient of variation across a species' accepted pairs.

    A real single-event allopolyploidization should diverge every homeolog
    pair to roughly the same depth relative to each pair's own baseline
    (tight distance_ratio spread, low CV). A messier/auto or multi-event
    history should show wide, inconsistent spread (high CV) -- this reframes
    distance_ratio spread explicitly as auto-vs-allo evidence rather than
    the "candidate asynchronous resolution, uncertain" framing used before
    this metric existed as a formal statistic. Needs >=2 accepted pairs.
    """
    path = os.path.join(homeologs_dir(species_dir), "ploidy_ancestry_summary.tsv")
    if not os.path.exists(path):
        return None, 0
    seen = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            partner = r.get("homeolog_partner")
            ratio = r.get("distance_ratio")
            if not partner or not ratio:
                continue
            key = frozenset([r["chrom"], partner])
            seen[key] = float(ratio)
    ratios = list(seen.values())
    if len(ratios) < 2:
        return None, len(ratios)
    mean = statistics.mean(ratios)
    if mean == 0:
        return None, len(ratios)
    return statistics.pstdev(ratios) / mean, len(ratios)


def main():
    species_list = sorted(
        d for d in os.listdir(RESULTS_DIR) if os.path.isdir(os.path.join(RESULTS_DIR, d))
    )

    fieldnames = [
        "species",
        "partition_consistency",
        "n_chroms_split",
        "modal_singleton_hap",
        "singleton_artifact_suspected",
        "mean_run_length_windows",
        "flip_rate",
        "n_windows_for_run_length",
        "mean_windowed_cv",
        "n_pairs_cv",
        "pair_depth_cv",
        "distance_ratio_cv",
        "n_pairs_for_ratio_cv",
        "n_accepted_pairs",
        "combined_allo_score",
        "n_metrics_available",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for species in species_list:
        species_dir = os.path.join(RESULTS_DIR, species)
        if not os.path.exists(os.path.join(matrix_dir(species_dir), "whole_chrom_distance_matrix.csv")):
            continue

        pc, n_split, singleton_hap = metric_partition_consistency(species_dir)
        artifact = check_singleton_artifact(species_dir, singleton_hap)
        mrl, flip_rate, n_windows = metric_run_length(species_dir)
        wcv, n_wcv = metric_windowed_cv(species_dir)
        dcv, n_pairs = metric_pair_depth_cv(species_dir)
        rcv, n_rcv = metric_distance_ratio_cv(species_dir)

        scores = []
        if pc is not None:
            scores.append(pc)
        if wcv is not None:
            scores.append(1.0 / (1.0 + wcv))
        if dcv is not None:
            scores.append(1.0 / (1.0 + dcv))

        combined = statistics.mean(scores) if scores else None

        writer.writerow(
            {
                "species": species,
                "partition_consistency": f"{pc:.4f}" if pc is not None else "",
                "n_chroms_split": n_split,
                "modal_singleton_hap": singleton_hap or "",
                "singleton_artifact_suspected": (
                    "" if artifact is None else ("yes" if artifact else "no")
                ),
                "mean_run_length_windows": f"{mrl:.3f}" if mrl is not None else "",
                "flip_rate": f"{flip_rate:.4f}" if flip_rate is not None else "",
                "n_windows_for_run_length": n_windows,
                "mean_windowed_cv": f"{wcv:.4f}" if wcv is not None else "",
                "n_pairs_cv": n_wcv,
                "pair_depth_cv": f"{dcv:.4f}" if dcv is not None else "",
                "distance_ratio_cv": f"{rcv:.4f}" if rcv is not None else "",
                "n_pairs_for_ratio_cv": n_rcv,
                "n_accepted_pairs": n_pairs,
                "combined_allo_score": f"{combined:.4f}" if combined is not None else "",
                "n_metrics_available": len(scores),
            }
        )


if __name__ == "__main__":
    main()
