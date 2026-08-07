import csv
import itertools
import math
import os
import statistics
from collections import defaultdict

from .common import homeologs_dir, log, matrix_dir
from .kmer_tables import load_sequences


def load_distance_matrix(outdir):
    path = os.path.join(matrix_dir(outdir), "whole_chrom_distance_matrix.csv")
    with open(path) as f:
        r = list(csv.reader(f))
    ids = r[0][1:]
    mat = {
        ids[i]: {ids[j]: float(v) for j, v in enumerate(row[1:])}
        for i, row in enumerate(r[1:])
    }
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


def load_pair_resolution_flags(outdir):
    """{frozenset({unit_a, unit_b}): resolution_limited bool} from whole_chrom_pairs.tsv.
    Empty dict (not an error) if that file predates the Mash-correction columns."""
    path = os.path.join(matrix_dir(outdir), "whole_chrom_pairs.tsv")
    flags = {}
    if not os.path.exists(path):
        return flags
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "resolution_limited" not in (reader.fieldnames or []):
            return flags
        for row in reader:
            flags[frozenset((row["unit_a"], row["unit_b"]))] = row[
                "resolution_limited"
            ] in ("True", "true", "1")
    return flags


def chrom_pair_resolution_status(groups, flags, i, j):
    """Resolution status of a chromosome-number pair, aggregated over every
    haplotype-copy combination between the two groups. Returns None if flags
    aren't available (e.g. pre-Mash-correction pairs.tsv) for any combination,
    else (any_limited, n_limited, n_total)."""
    keys = [frozenset((a, b)) for a in groups[i] for b in groups[j]]
    statuses = [flags.get(k) for k in keys]
    if any(s is None for s in statuses):
        return None
    n_limited = sum(1 for s in statuses if s)
    return n_limited > 0, n_limited, len(statuses)


def _fit_background(distances, n_iter=3, clip_sigma=2.0):
    """Iteratively refit mean/stdev of the background, excluding low outliers
    (candidate homeolog pairs) each round so a genome where most chromosomes
    retain a homeolog partner doesn't contaminate its own null estimate."""
    sample = list(distances)
    mean = statistics.mean(sample)
    stdev = statistics.pstdev(sample) or 1e-9
    for _ in range(n_iter):
        kept = [d for d in distances if d >= mean - clip_sigma * stdev]
        if len(kept) < 2 or len(kept) == len(sample):
            break
        sample = kept
        mean = statistics.mean(sample)
        stdev = statistics.pstdev(sample) or 1e-9
    return mean, stdev


def empirical_pvalues(pair_dist):
    """
    Empirical p-value per cross-chromosome-number pair: a one-sided lower-tail
    z-test against a normal approximation fit to the background of all
    cross-chromosome-number distances (sigma-clipped so candidate homeolog
    pairs don't inflate the estimate -- see _fit_background). Not a literal
    reshuffle-and-recompute permutation test (that would need re-deriving
    k-mer distances under permuted labels, real new compute this codebase has
    no machinery for), but the same underlying question: is this pair's
    distance unusually small relative to what unrelated-chromosome
    comparisons look like.

    A naive "fraction of the same n distances that are <= mine" p-value
    (i.e. rank/n) was tried first and rejected after verification caught a
    real degeneracy: it's mathematically identical to feeding
    Benjamini-Hochberg its own expected-value-under-the-null as the observed
    p-value, so q_(k) = p_(k)*n/k = (k/n)*n/k = 1 for every single pair
    regardless of the data -- confirmed empirically on synthetic data with an
    obvious 9-pair signal (0/9 accepted). The p-value has to be able to
    resolve finer than 1/n for FDR correction to do anything at all; a
    parametric fit to the background achieves that.

    Returns (p_values, z_scores). z_scores is reported alongside p in every
    output because p underflows to exactly 0.0 in float64 once |z| gets much
    past ~9 (confirmed on real data: daGleHede1's homeolog pairs sit ~20
    background-stdevs below the mean, since that background happens to be
    unusually tight) -- at that point p and q genuinely are indistinguishable
    from 0 and stop conveying how much evidence there is, while z stays a
    finite, comparable number.
    """
    distances = list(pair_dist.values())
    mean, stdev = _fit_background(distances)
    z_scores = {pair: (d - mean) / stdev for pair, d in pair_dist.items()}
    p_values = {
        pair: 0.5 * (1 + math.erf(z / math.sqrt(2))) for pair, z in z_scores.items()
    }
    return p_values, z_scores


def bh_qvalues(p_by_key):
    """Benjamini-Hochberg FDR-adjusted p-values (q-values), stdlib only --
    q_(rank) = min(1, min_{rank' >= rank} p_(rank') * n / rank'), computed as a
    running minimum from the largest rank down to the smallest so it comes
    out monotone non-decreasing as p increases, per the standard BH procedure."""
    items = list(p_by_key.items())
    n = len(items)
    order = sorted(range(n), key=lambda idx: items[idx][1])
    q_values = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        p = items[idx][1]
        running_min = min(running_min, p * n / rank)
        q_values[idx] = running_min
    return {items[idx][0]: q_values[idx] for idx in range(n)}


def detect_homeolog_pairs(pair_dist, chrom_nums, fdr_alpha=0.05):
    """
    Look for a retained ancestral (paleopolyploid/WGD) subgenome pairing among
    *different* chromosome numbers: chromosome pairs whose whole-chromosome
    k-mer distance is unusually low relative to the empirical background of
    every other cross-chromosome-number comparison (see empirical_pvalues),
    with Benjamini-Hochberg FDR correction across all candidates tested
    simultaneously. Pairs with q <= fdr_alpha are accepted, then greedily
    resolved into a 1:1 matching (ascending distance order, each chromosome
    number gets at most one partner) since a chromosome can only have one
    true ancestral partner.
    """
    p_values, z_scores = empirical_pvalues(pair_dist)
    q_values = bh_qvalues(p_values)

    sorted_pairs = sorted(pair_dist.items(), key=lambda kv: kv[1])
    matched = {}
    accepted = []
    for (i, j), d in sorted_pairs:
        if q_values[(i, j)] > fdr_alpha:
            continue
        if i in matched or j in matched:
            continue
        matched[i] = j
        matched[j] = i
        accepted.append(
            (i, j, d, p_values[(i, j)], q_values[(i, j)], z_scores[(i, j)])
        )

    unmatched = [c for c in chrom_nums if c not in matched]
    accepted_keys = {(i, j) for i, j, *_ in accepted}
    background = [d for pair, d in pair_dist.items() if pair not in accepted_keys]
    return accepted, unmatched, background


def run(seq_tsv, outdir, fdr_alpha=0.05):
    units = load_sequences(seq_tsv)
    ids, mat = load_distance_matrix(outdir)
    groups = chrom_groups(units)
    chrom_nums = sorted(groups)
    if len(chrom_nums) < 4:
        log(
            f"only {len(chrom_nums)} chromosome numbers -- too few to search for a paleopolyploid pairing"
        )
        return []

    pair_dist = cross_chrom_distances(groups, mat)
    accepted, unmatched, background = detect_homeolog_pairs(
        pair_dist, chrom_nums, fdr_alpha
    )
    flags = load_pair_resolution_flags(outdir)
    status_by_pair = {
        (i, j): chrom_pair_resolution_status(groups, flags, i, j)
        for i, j, *_ in accepted
    }

    path = os.path.join(homeologs_dir(outdir), "homeolog_pairs.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "chrom_a",
                "chrom_b",
                "mean_distance",
                "z_score",
                "p_value",
                "q_value",
                "n_haplotype_copy_pairs",
                "resolution_limited",
                "n_resolution_limited_of_total",
            ]
        )
        for i, j, d, p, q, z in sorted(accepted, key=lambda x: x[2]):
            status = status_by_pair[(i, j)]
            limited_str = "NA" if status is None else status[0]
            n_limited_str = "NA" if status is None else f"{status[1]}/{status[2]}"
            w.writerow(
                [
                    f"chr{i:02d}",
                    f"chr{j:02d}",
                    f"{d:.6f}",
                    f"{z:.3f}",
                    f"{p:.6g}",
                    f"{q:.6g}",
                    len(groups[i]) * len(groups[j]),
                    limited_str,
                    n_limited_str,
                ]
            )

    bg_mean = statistics.mean(background) if background else float("nan")
    bg_min = min(background) if background else float("nan")
    log(
        f"candidate ancient homeolog pairs: {len(accepted)} pairs across {2 * len(accepted)} of "
        f"{len(chrom_nums)} chromosomes (FDR alpha={fdr_alpha}, background min={bg_min:.4f} mean={bg_mean:.4f})"
    )
    for i, j, d, p, q, z in sorted(accepted, key=lambda x: x[2]):
        status = status_by_pair[(i, j)]
        warn = ""
        if status is not None and status[0]:
            warn = (
                f" [WARNING: {status[1]}/{status[2]} underlying haplotype-pair distances "
                f"are resolution-limited -- this pairing rests on noise-floor estimates]"
            )
        log(
            f"  chr{i:02d} <-> chr{j:02d}: mean distance={d:.4f} z={z:.2f} q={q:.4g}{warn}"
        )
    if unmatched:
        log(
            f"  no significant ancestral partner found for: {', '.join(f'chr{c:02d}' for c in unmatched)}"
        )

    plot_ranked_distances(outdir, pair_dist, accepted)
    return [(i, j) for i, j, *_ in accepted]


def plot_ranked_distances(outdir, pair_dist, accepted):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sorted_pairs = sorted(pair_dist.items(), key=lambda kv: kv[1])
    accepted_keys = {(min(i, j), max(i, j)) for i, j, *_ in accepted}

    xs = list(range(len(sorted_pairs)))
    ys = [d for _, d in sorted_pairs]
    colors = [
        (
            "crimson"
            if (min(k[0], k[1]), max(k[0], k[1])) in accepted_keys
            else "steelblue"
        )
        for k, _ in sorted_pairs
    ]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(xs, ys, c=colors, s=14)
    for idx, ((i, j), d) in enumerate(sorted_pairs):
        if (min(i, j), max(i, j)) in accepted_keys:
            ax.annotate(
                f"{i:02d}-{j:02d}",
                (idx, d),
                fontsize=6,
                xytext=(2, 4),
                textcoords="offset points",
            )
    ax.set_xlabel("cross-chromosome pair rank (ascending mean distance)")
    ax.set_ylabel("mean whole-chromosome k-mer Jaccard distance")
    ax.set_title(
        "Candidate ancestral (paleopolyploid) homeolog pairs -- red = accepted"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(homeologs_dir(outdir), "homeolog_pairs.png"), dpi=150)
    plt.close(fig)
