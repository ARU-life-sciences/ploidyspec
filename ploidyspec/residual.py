"""
Residual tetrasomy between homeologs: stretches of a chromosome where its
homeolog copies are as close to it as its own allelic copy.

After an autopolyploidization, rediploidization proceeds chromosome by
chromosome and region by region. Where the duplicated copies have stopped
recombining they diverge, and their k-mer distance grows well beyond the
allelic distance; where they still pair and exchange (residual tetrasomy, as
at the ends of many salmonid chromosomes) they stay as close as alleles. A
two-haplotype assembly of such a genome numbers the duplicated copies as
different chromosomes, and the whole-chromosome distance only shows the
diverged majority; this test reads the windowed-homeologs tracks to find the
regions that have not diverged.

Per window of each anchor copy (positions along that copy, position-free
comparison; see windowed.py): the distance to the closest copy of the
homeolog partner, over the anchor's typical allelic distance (median over its
windows of the distance to its own chromosome's other copies). A window is
*near-allelic* when that ratio is below `ratio` (default 3) and the anchor
is not also near-allelic to an unrelated control chromosome in that window
(windowed.control_pairs). Shared satellites and subtelomeric repeats make a
window close to every chromosome carrying them -- ddSalPent1's chromosome ends
read 0.001 against their salicoid homeologs, below the allelic distance --
while residual tetrasomy is specific to the homeolog. Runs from before the
controls existed are read without them and flagged `controlled = no`.

Null: near-allelic windows also turn up singly (shared repeats, a misplaced
contig). Runs are scanned with gaps of one window bridged; the near-allelic
labels are shuffled across all tracks of the species and the longest run
recorded, repeated `n_perm` times. A run counts only if it has more
near-allelic windows than the 95th percentile of that genome-wide maximum, and
is at least `min_segment_bp` long. A homeolog pair counts as residual-
tetrasomic when more than half the copies of its two numbers show such a run
(3 of 4 with two haplotypes): a real tetrasomic region shows from every copy,
but a contig scaffolded onto the wrong homeolog -- its tetrasomic sequence
aligns equally to both -- only from the two copies it touches. Each
chromosome then reports its best-covered copy.

Needs `windowed` and `windowed-homeologs` outputs; returns nothing without them.
"""

import os
import random
import statistics
from collections import defaultdict

from .common import homeologs_dir, window_rows, windowed_dir
from .windowed import CONTROLS_TSV

RESIDUAL_RATIO = 3.0
RESIDUAL_NULL_QUANTILE = 0.95
RESIDUAL_PERMUTATIONS = 200
ALLELE_FLOOR = 1e-3  # an almost homozygous chromosome must not inflate every ratio
TERMINAL_FRAC = 0.1  # a segment within this fraction of either end is terminal
MAX_GAP = 1

FIELDS = ["chrom", "partner", "anchor", "control", "n_windows", "n_near_allelic", "n_shared_with_control",
          "allele_level",
          "longest_run", "null_max_run", "segments", "covered_bp", "covered_frac"]


def load_tracks(path):
    """{(anchor, partner_chrom): {win_start: (win_end, {other_unit: distance})}}"""
    out = defaultdict(dict)
    if not os.path.exists(path):
        return out
    for r, d in window_rows(path):
        key = (r["unit_a"], int(r["chrom_b"]))
        entry = out[key].setdefault(int(r["win_start"]), (int(r["win_end"]), {}))
        entry[1][r["unit_b"]] = d
    return out


def allele_levels(path):
    """{anchor: median over its windows of the mean distance to its own other copies}"""
    per = defaultdict(list)
    if not os.path.exists(path):
        return {}
    by_win = defaultdict(list)
    for r, d in window_rows(path):
        by_win[(r["unit_a"], r["win_start"])].append(d)
    for (anchor, _), ds in by_win.items():
        per[anchor].append(statistics.mean(ds))
    return {a: statistics.median(v) for a, v in per.items()}


def near_allelic_track(windows, allele, ratio, control=None):
    """[(start, end, near_allelic, shared_with_control)] along one anchor. A
    window close to the control chromosome as well is not near-allelic."""
    cut = ratio * max(allele, ALLELE_FLOOR)
    out = []
    for s in sorted(windows):
        near = min(windows[s][1].values()) < cut
        ctl = control.get(s) if control else None
        shared = bool(ctl and min(ctl[1].values()) < cut)
        out.append((s, windows[s][0], near and not shared, near and shared))
    return out


def runs(flags, max_gap=MAX_GAP):
    """Runs of True in a list of booleans, bridging up to max_gap Falses.
    Returns [(first_index, last_index, n_true)]."""
    out, first, last, n, gap = [], None, None, 0, 0
    for i, f in enumerate(flags):
        if f:
            if first is None:
                first = i
            last, n, gap = i, n + 1, 0
        elif first is not None:
            gap += 1
            if gap > max_gap:
                out.append((first, last, n))
                first, n, gap = None, 0, 0
    if first is not None:
        out.append((first, last, n))
    return out


def null_max_run(lengths, n_true, n_perm, quantile, seed=0):
    """The `quantile` of the genome-wide longest run when n_true near-allelic labels
    are scattered at random over tracks of the given lengths."""
    rng = random.Random(seed)
    total = sum(lengths)
    if not total or not n_true:
        return 0
    maxima = []
    for _ in range(n_perm):
        hits = set(rng.sample(range(total), n_true))
        best, offset = 0, 0
        for L in lengths:
            flags = [offset + i in hits for i in range(L)]
            best = max([best] + [r[2] for r in runs(flags)])
            offset += L
        maxima.append(best)
    maxima.sort()
    return maxima[int(quantile * (len(maxima) - 1))]


def position(start, end, span_start, span_end):
    lo = span_start + TERMINAL_FRAC * (span_end - span_start)
    hi = span_end - TERMINAL_FRAC * (span_end - span_start)
    if start <= lo:
        return "start"
    if end >= hi:
        return "end"
    return "interior"


def residual_tetrasomy(outdir, chrom_of, thresholds):
    """
    chrom_of: {unit_id: chromosome number}. Returns (rows per anchor,
    {chrom: reading of its median anchor}, null threshold) -- empty when the
    windowed-homeologs tracks are missing.
    """
    tracks = load_tracks(os.path.join(homeologs_dir(outdir), "windowed_homeologs_all.tsv"))
    if not tracks:
        return [], {}, None
    alleles = allele_levels(os.path.join(windowed_dir(outdir), "windowed_all.tsv"))
    controls = load_tracks(os.path.join(homeologs_dir(outdir), CONTROLS_TSV))
    control_of = {anchor: (c, w) for (anchor, c), w in controls.items()}
    ratio = thresholds.get("residual_ratio", RESIDUAL_RATIO)
    min_bp = thresholds["min_segment_bp"]

    built = {}
    for (anchor, partner), windows in sorted(tracks.items()):
        if anchor not in alleles:
            continue
        built[(anchor, partner)] = near_allelic_track(windows, alleles[anchor], ratio,
                                                      control_of.get(anchor, (None, None))[1])
    lengths = [len(t) for t in built.values()]
    n_true = sum(t[2] for track in built.values() for t in track)
    null = null_max_run(lengths, n_true, thresholds.get("residual_permutations", RESIDUAL_PERMUTATIONS),
                        RESIDUAL_NULL_QUANTILE)

    rows, by_chrom = [], defaultdict(list)
    for (anchor, partner), track in built.items():
        span_start, span_end = track[0][0], track[-1][1]
        all_runs = runs([t[2] for t in track])
        segs = []
        for i, j, n in all_runs:
            s, e = track[i][0], track[j][1]
            if n > null and e - s + 1 >= min_bp:
                segs.append((s, e, position(s, e, span_start, span_end)))
        covered = sum(e - s + 1 for s, e, _ in segs)
        frac = covered / (span_end - span_start + 1)
        chrom = chrom_of[anchor]
        ctl = control_of.get(anchor, (None, None))[0]
        row = dict(chrom=f"chr{chrom:02d}", partner=f"chr{partner:02d}", anchor=anchor,
                   control=f"chr{ctl:02d}" if ctl is not None else "",
                   n_windows=len(track), n_near_allelic=sum(t[2] for t in track),
                   n_shared_with_control=sum(t[3] for t in track),
                   allele_level=f"{alleles[anchor]:.5f}",
                   longest_run=max((n for _, _, n in all_runs), default=0), null_max_run=null,
                   segments=";".join(f"{s / 1e6:.2f}-{e / 1e6:.2f}Mb({p})" for s, e, p in segs),
                   covered_bp=covered, covered_frac=f"{frac:.3f}")
        rows.append(row)
        by_chrom[chrom].append((covered, frac, anchor, segs, partner))

    # decided per homeolog pair, over the copies of both numbers: real residual
    # tetrasomy shows from every copy, a contig placed on the wrong homeolog (its
    # tetrasomic sequence aligns equally to both) from only the two it touches
    pair_reads = defaultdict(list)
    for chrom, reads in by_chrom.items():
        for covered, frac, anchor, segs, partner in reads:
            pair_reads[tuple(sorted((chrom, partner)))].append(covered > 0)
    readings = {}
    for chrom, reads in by_chrom.items():
        partner = reads[0][4]
        votes = pair_reads[tuple(sorted((chrom, partner)))]
        passes = sum(votes) > len(votes) / 2
        covered, frac, anchor, segs, _ = max(reads) if passes else (0, 0.0, "", [], partner)
        readings[chrom] = dict(
            residual_bp=covered, residual_frac=frac, partner=partner,
            residual_support=f"{sum(votes)}/{len(votes)}",
            controlled="yes" if controls else "no",
            residual_segments=";".join(
                f"{anchor.split('_')[0]}:{s / 1e6:.1f}-{e / 1e6:.1f}Mb({p})" for s, e, p in segs),
            terminal_bp=sum(e - s + 1 for s, e, p in segs if p != "interior"))
    return rows, readings, null
