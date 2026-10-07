"""
Per-species genome-structure metrics, read from the matrix, homeologs and
windowed outputs (no new k-mer work). Formerly scripts/auto_allo_spectrum.py and
scripts/genome_partition.py, which ran over the whole panel at once; the
cross-species tables are now built by `ploidyspec panel` from each species'
structure/ outputs.

Writes, in results/<species>/structure/:
- inheritance_metrics.tsv: one row of auto/allo-spectrum metrics --
  partition_consistency, run length / flip rate, windowed CV, pair-depth CV,
  distance-ratio CV and a heuristic combined score. See OUTPUTS.md and
  INTERPRETATION.md ("A continuous auto/allo spectrum") for what each means and
  which ones failed as discriminators on the panel.
- genome_partition.tsv: diffuse genome-wide chromosome partitions -- whether the
  chromosome numbers factor into k groups more than chance, at every k with
  |z| >= 2 (INTERPRETATION.md, "Diffuse, genome-wide partitions").

Metric details (from the original panel script):

1. partition_consistency: does the same lineage-partition of >=3 haplotype
   copies recur across every chromosome (disomic/allo signature), or does
   the "odd one out" rotate between different haplotype labels on different
   chromosomes (tetrasomic/auto signature)? When the recurring (modal)
   partition is one copy vs the rest (modal_singleton_hap is set), it measures
   one consistently divergent haplotype -- an assembly artefact or an AAAB-like
   copy, not two subgenomes (ddLepDrab1: HAP4 on 15/16 chromosomes). Read it
   that way; it still enters combined_allo_score, since dropping it there leaves
   only the near-constant CV terms.
2. mean_windowed_cv: average coefficient of variation of the windowed
   distance track across accepted homeolog pairs. Low = uniform divergence.
3. pair_depth_cv: CV of accepted homeolog pairs' mean_distance.
4. mean_run_length_windows / flip_rate: window-by-window agreement with the
   whole-chromosome partition. Did NOT separate confirmed autos from allos on
   the panel -- informational only, not in the combined score.
5. distance_ratio_cv: CV of accepted pairs' distance_ratio. Confounded by
   tetrasomic homogenization (drLytSali1) -- not in the combined score.

combined_allo_score averages the allo-direction mapping of metrics 1-3 that
are available (1/(1+x) for the CVs). It is a heuristic, NOT a classifier.
"""
import csv
import os
import random
import re
import statistics

from .common import homeologs_dir, log, matrix_dir, window_rows, windowed_dir
from .kmer_tables import load_sequences
from .subgenome_report import bipartition_by_distance, load_whole_chrom_distances

N_PERMUTATIONS = 500
MIN_Z = 2.0  # only report k's whose observed separation clears this
PARTITION_SEED = 20260906


def structure_dir(outdir):
    d = os.path.join(outdir, "structure")
    os.makedirs(d, exist_ok=True)
    return d


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

        if min(len(set_a), len(set_b)) < 2:
            continue

        chrom_str = f"chr{chrom:02d}"
        path = os.path.join(windowed_dir(species_dir), f"windowed_{chrom_str}.tsv")
        if not os.path.exists(path):
            continue

        # each copy's windows along its own coordinates: {anchor: {win: {other: d}}}
        by_anchor = {}
        for r, d in window_rows(path):
            by_anchor.setdefault(r["unit_a"], {}).setdefault(int(r["win_start"]), {})[r["unit_b"]] = d

        for anchor, by_window in by_anchor.items():
            own = set_a if anchor in set_a else set_b
            agreements = []
            for win in sorted(by_window):
                vals = by_window[win]
                self_vals = [d for u, d in vals.items() if u in own]
                cross_vals = [d for u, d in vals.items() if u not in own]
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
        for r, d in window_rows(path):
            by_pair.setdefault((r["unit_a"], r["unit_b"]), []).append(d)
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


def load_chrom_distance_matrix(species_dir):
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


def null_distribution(dist, chroms, sizes, rng, n=N_PERMUTATIONS):
    ratios = []
    pool = list(chroms)
    for _ in range(n):
        rng.shuffle(pool)
        clusters = []
        idx = 0
        for s in sizes:
            clusters.append(pool[idx : idx + s])
            idx += s
        r = separation_ratio(dist, clusters)
        if r is not None:
            ratios.append(r)
    return ratios


def analyze_species(species_dir, species_name, seed=PARTITION_SEED):
    """Every k with a valid partition, strongest |z| first. The null uses its own
    RNG seeded per species, so a species' z-scores don't depend on which other
    species were analysed before it (the old panel script shared one RNG)."""
    rng = random.Random(seed)
    dist, chroms = load_chrom_distance_matrix(species_dir)
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
        null = null_distribution(dist, chroms, sizes, rng)
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
    results.sort(key=lambda r: -abs(r["z_score"]))
    return results



INHERITANCE_FIELDS = [
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
PARTITION_FIELDS = ["species", "k", "group_sizes", "separation_ratio", "null_mean_ratio",
                    "z_score", "groups"]


def inheritance_metrics(species_dir, species):
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

    def f(v, nd):
        return f"{v:.{nd}f}" if v is not None else ""

    return {
        "species": species,
        "partition_consistency": f(pc, 4),
        "n_chroms_split": n_split,
        "modal_singleton_hap": singleton_hap or "",
        "singleton_artifact_suspected": "" if artifact is None else ("yes" if artifact else "no"),
        "mean_run_length_windows": f(mrl, 3),
        "flip_rate": f(flip_rate, 4),
        "n_windows_for_run_length": n_windows,
        "mean_windowed_cv": f(wcv, 4),
        "n_pairs_cv": n_wcv,
        "pair_depth_cv": f(dcv, 4),
        "distance_ratio_cv": f(rcv, 4),
        "n_pairs_for_ratio_cv": n_rcv,
        "n_accepted_pairs": n_pairs,
        "combined_allo_score": f(combined, 4),
        "n_metrics_available": len(scores),
    }


def partition_rows(species_dir, species):
    rows = []
    for r in analyze_species(species_dir, species):
        if abs(r["z_score"]) < MIN_Z:
            continue
        rows.append({
            "species": species,
            "k": r["k"],
            "group_sizes": r["group_sizes"],
            "separation_ratio": f"{r['separation_ratio']:.4f}",
            "null_mean_ratio": f"{r['null_mean_ratio']:.4f}",
            "z_score": f"{r['z_score']:.2f}",
            "groups": " | ".join(",".join(g) for g in r["groups"]),
        })
    return rows


def write_rows(path, fields, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def compute_structure(outdir):
    species = os.path.basename(os.path.normpath(outdir))
    sdir = structure_dir(outdir)
    write_rows(os.path.join(sdir, "inheritance_metrics.tsv"), INHERITANCE_FIELDS,
               [inheritance_metrics(outdir, species)])
    parts = partition_rows(outdir, species)
    write_rows(os.path.join(sdir, "genome_partition.tsv"), PARTITION_FIELDS, parts)
    log(f"wrote inheritance_metrics.tsv and genome_partition.tsv "
        f"({len(parts)} partition(s) with |z| >= {MIN_Z}) in {sdir}")
