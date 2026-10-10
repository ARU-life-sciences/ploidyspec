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
*near-allelic* when that ratio is below `ratio` (default 3), the distance is
also below `HOMEOLOG_FRAC` of the species' median homeolog distance, and the
anchor is not also near-allelic to an unrelated control chromosome in that window
(windowed.control_pairs). Shared satellites and subtelomeric repeats make a
window close to every chromosome carrying them -- ddSalPent1's chromosome ends
read 0.001 against their salicoid homeologs, below the allelic distance --
while residual tetrasomy is specific to the homeolog. Runs from before the
controls existed are read without them and flagged `controlled = no`.

The allelic level is the median over windows of the distance to the anchor's
*closest* other copy, and a chromosome's copies all use the median of their
levels, capped at ALLELE_CAP times the genome median so that a chromosome
carrying a large exchange does not raise its own yardstick. With more than two
copies one is often a divergent subgenome copy;
averaging over it (or anchoring on it) put the "allelic" level near the
homeolog distance, and three times that let most of the genome through
(ddLepDrab1, ddHypMacu1). The homeolog ceiling guards the same failure from
the other side: in a species whose homeologs are barely more diverged than its
alleles (daGleHede1: 0.03 vs 0.009) the ratio alone reads the near tail of
ordinary homeolog divergence as residual tetrasomy, while real residual
segments sit far below the homeolog level (OncMaso1: 0.0045 vs 0.048).

In residual tetrasomy every copy is close to every other, so a window also
needs the anchor's own homolog within the allelic cut. A window close to a
homeolog copy but far from its own homolog is a copy carrying its homeolog's
sequence -- a homeologous exchange, or a contig placed on the wrong homeolog
-- and is read as exchange instead. Without this an exchange in one copy
passed the pair vote from the three copies that see it (the simulated
segmental_homeology chr04), and so did homeolog-swapped contigs in OncMaso1's
reference-scaffolded HAP2: 6 of its 8 earlier residual chromosomes (54 Mb)
rested on windows whose own homolog was 0.03-0.085 away.

Null: near-allelic windows also turn up singly (shared repeats, a misplaced
contig). Runs are scanned with gaps of one window bridged; the near-allelic
labels are shuffled across all tracks of the species and the longest run
recorded, repeated `n_perm` times. A run counts only if it has more
near-allelic windows than the 95th percentile of that genome-wide maximum, and
is at least `min_segment_bp` long. In an almost homozygous genome near-allelic
windows are so rare that the shuffled maximum falls to 3-4 windows, and a 2 Mb
cluster at a chromosome end passes (laPotCris1); the null is floored at
`MIN_NULL_RUN` windows. A homeolog pair counts as residual-
tetrasomic when more than half the copies of its two numbers show such a run
(3 of 4 with two haplotypes): a real tetrasomic region shows from every copy,
but a contig scaffolded onto the wrong homeolog -- its tetrasomic sequence
aligns equally to both -- only from the two copies it touches. Each
chromosome then reports its best-covered copy.

Partners come from whole-chromosome homeolog pairs and from the homeology
map's blocks (homeology_map.py), so a chromosome whose arms have different
homeologs is tested against each; every reading is per (chromosome, partner)
pair and a chromosome sums its pairs.

Homeologous exchange, the converse reading from the same tracks: windows where
the anchor is closer to a homeolog copy than to its own chromosome's other
copies (homeolog distance < own distance / EXCHANGE_RATIO, and the own distance
beyond `ratio` times the allelic level). The anchor carries a stretch of its
homeolog's sequence: a homeologous exchange (non-reciprocal homeologous
recombination, a hallmark of allopolyploids), or a contig scaffolded onto the
wrong chromosome (masu's reference-scaffolded HAP1). Runs are tested against
their own shuffle null; no pair vote, since an exchange is usually carried by
one copy -- `exchange_support` says how many.

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
ALLELE_CAP = 1.5  # a chromosome's allelic level is at most this times the genome median
EXCHANGE_RATIO = 2.0  # homeolog copy this much closer than the anchor's own copies
HOMEOLOG_FRAC = 0.25  # near-allelic also means well below the typical homeolog distance
MIN_NULL_RUN = 8  # floor on the shuffled-run null, in windows

FIELDS = ["chrom", "partner", "anchor", "control", "n_windows", "n_near_allelic", "n_shared_with_control",
          "allele_level",
          "longest_run", "null_max_run", "segments", "covered_bp", "covered_frac",
          "n_exchange", "exchange_segments", "exchange_bp"]


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


def allele_levels(path, chrom_of=None):
    """{anchor: allelic level}: the median over its windows of the distance to its
    closest other copy; with chrom_of, every copy of a chromosome gets the median
    of those levels, so a divergent copy is measured against its chromosome's
    allelic level rather than its own distance to the rest."""
    per = defaultdict(list)
    if not os.path.exists(path):
        return {}
    by_win = defaultdict(list)
    for r, d in window_rows(path):
        by_win[(r["unit_a"], r["win_start"])].append(d)
    for (anchor, _), ds in by_win.items():
        per[anchor].append(min(ds))
    levels = {a: statistics.median(v) for a, v in per.items()}
    if not chrom_of:
        return levels
    by_chrom = defaultdict(list)
    for a, v in levels.items():
        by_chrom[chrom_of.get(a)].append(v)
    chrom_level = {c: statistics.median(v) for c, v in by_chrom.items()}
    if not chrom_level:
        return {}
    # a large exchange raises its own chromosome's median (simulated chr04,
    # half exchanged: 0.019 against 0.0035 elsewhere) and would hide itself
    cap = ALLELE_CAP * statistics.median(chrom_level.values())
    return {a: min(chrom_level[chrom_of.get(a)], cap) for a in levels}


def own_copy_distances(path):
    """{(anchor, win_start): distance to the anchor's closest own-number copy}"""
    out = {}
    if not os.path.exists(path):
        return out
    for r, d in window_rows(path):
        key = (r["unit_a"], int(r["win_start"]))
        out[key] = min(d, out.get(key, d))
    return out


def exchange_track(anchor, windows, own, allele, ratio):
    """[(start, end, exchange)]: windows closer to a homeolog copy than to the
    anchor's own copies (see module docstring)."""
    floor = ratio * max(allele, ALLELE_FLOOR)
    out = []
    for s in sorted(windows):
        d_own = own.get((anchor, s))
        d_hom = min(windows[s][1].values())
        hit = d_own is not None and d_own > floor and d_hom < d_own / EXCHANGE_RATIO
        out.append((s, windows[s][0], hit))
    return out


def near_allelic_track(windows, allele, ratio, control=None, ceiling=None, own=None):
    """[(start, end, near_allelic, shared_with_control)] along one anchor. A
    window close to the control chromosome as well is not near-allelic; the
    cut never exceeds `ceiling` (a fraction of the homeolog distance). With
    `own` ({start: distance to the anchor's own other copies}) a window whose
    own homolog is beyond the allelic cut is not near-allelic either: it is an
    exchange-like window (exchange_track), not residual tetrasomy."""
    cut = ratio * max(allele, ALLELE_FLOOR)
    own_cut = cut
    if ceiling is not None:
        cut = min(cut, ceiling)
    out = []
    for s in sorted(windows):
        near = min(windows[s][1].values()) < cut
        if near and own and own.get(s) is not None and own[s] > own_cut:
            near = False
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
    alleles = allele_levels(os.path.join(windowed_dir(outdir), "windowed_all.tsv"), chrom_of)
    controls = load_tracks(os.path.join(homeologs_dir(outdir), CONTROLS_TSV))
    control_of = {anchor: (c, w) for (anchor, c), w in controls.items()}
    ratio = thresholds.get("residual_ratio", RESIDUAL_RATIO)
    min_bp = thresholds["min_segment_bp"]
    homeolog_level = statistics.median(min(d.values()) for w in tracks.values() for _, d in w.values())
    ceiling = thresholds.get("residual_homeolog_frac", HOMEOLOG_FRAC) * homeolog_level

    own = own_copy_distances(os.path.join(windowed_dir(outdir), "windowed_all.tsv"))
    built = {}
    for (anchor, partner), windows in sorted(tracks.items()):
        if anchor not in alleles:
            continue
        built[(anchor, partner)] = near_allelic_track(
            windows, alleles[anchor], ratio, control_of.get(anchor, (None, None))[1], ceiling,
            {s: own.get((anchor, s)) for s in windows})
    lengths = [len(t) for t in built.values()]
    n_true = sum(t[2] for track in built.values() for t in track)
    null = null_max_run(lengths, n_true, thresholds.get("residual_permutations", RESIDUAL_PERMUTATIONS),
                        RESIDUAL_NULL_QUANTILE)
    null = max(null, thresholds.get("residual_min_run", MIN_NULL_RUN))

    ex_tracks = {key: exchange_track(key[0], tracks[key], own, alleles[key[0]], ratio) for key in built}
    ex_null = null_max_run([len(t) for t in ex_tracks.values()],
                           sum(t[2] for track in ex_tracks.values() for t in track),
                           thresholds.get("residual_permutations", RESIDUAL_PERMUTATIONS),
                           RESIDUAL_NULL_QUANTILE, seed=1)
    ex_null = max(ex_null, thresholds.get("residual_min_run", MIN_NULL_RUN))

    def segments(track, flags_at, threshold):
        span_start, span_end = track[0][0], track[-1][1]
        all_runs = runs([t[flags_at] for t in track])
        segs = []
        for i, j, n in all_runs:
            s, e = track[i][0], track[j][1]
            if n > threshold and e - s + 1 >= min_bp:
                segs.append((s, e, position(s, e, span_start, span_end)))
        return all_runs, segs, span_end - span_start + 1

    rows, reads_by_pair = [], defaultdict(list)
    for (anchor, partner), track in built.items():
        all_runs, segs, span = segments(track, 2, null)
        _, ex_segs, _ = segments(ex_tracks[(anchor, partner)], 2, ex_null)
        covered = sum(e - s + 1 for s, e, _ in segs)
        ex_bp = sum(e - s + 1 for s, e, _ in ex_segs)
        chrom = chrom_of[anchor]
        ctl = control_of.get(anchor, (None, None))[0]
        rows.append(dict(chrom=f"chr{chrom:02d}", partner=f"chr{partner:02d}", anchor=anchor,
                         control=f"chr{ctl:02d}" if ctl is not None else "",
                         n_windows=len(track), n_near_allelic=sum(t[2] for t in track),
                         n_shared_with_control=sum(t[3] for t in track),
                         allele_level=f"{alleles[anchor]:.5f}",
                         longest_run=max((n for _, _, n in all_runs), default=0), null_max_run=null,
                         segments=";".join(f"{s / 1e6:.2f}-{e / 1e6:.2f}Mb({p})" for s, e, p in segs),
                         covered_bp=covered, covered_frac=f"{covered / span:.3f}",
                         n_exchange=sum(t[2] for t in ex_tracks[(anchor, partner)]),
                         exchange_segments=";".join(f"{s / 1e6:.2f}-{e / 1e6:.2f}Mb({p})"
                                                    for s, e, p in ex_segs),
                         exchange_bp=ex_bp))
        reads_by_pair[tuple(sorted((chrom, partner)))].append(
            dict(chrom=chrom, partner=partner, anchor=anchor, covered=covered, frac=covered / span,
                 segs=segs, ex_bp=ex_bp, ex_segs=ex_segs))

    # residual tetrasomy is decided per homeolog pair, over the copies of both
    # numbers: real residual tetrasomy shows from every copy, a contig placed on
    # the wrong homeolog (its tetrasomic sequence aligns equally to both) from
    # only the two it touches
    readings = defaultdict(lambda: dict(residual_bp=0, residual_frac=0.0, partners=[], support=[],
                                        segments=[], terminal_bp=0, exchange_bp=0,
                                        exchange_segments=[], exchange_support=[]))
    for pair, reads in reads_by_pair.items():
        votes = [r["covered"] > 0 for r in reads]
        passes = sum(votes) > len(votes) / 2
        for chrom in sorted({r["chrom"] for r in reads}):
            mine = [r for r in reads if r["chrom"] == chrom]
            partner = mine[0]["partner"]
            rd = readings[chrom]
            rd["partners"].append(partner)
            rd["support"].append(f"chr{partner:02d}:{sum(votes)}/{len(votes)}")
            if passes:
                best = max(mine, key=lambda r: (r["covered"], r["anchor"]))
                rd["residual_bp"] += best["covered"]
                rd["residual_frac"] += best["frac"]
                rd["terminal_bp"] += sum(e - s + 1 for s, e, p in best["segs"] if p != "interior")
                rd["segments"] += [f"{best['anchor'].split('_')[0]}:{s / 1e6:.1f}-{e / 1e6:.1f}Mb({p})"
                                   + f"~chr{partner:02d}"
                                   for s, e, p in best["segs"]]
            ex = [r for r in mine if r["ex_bp"]]
            if ex:
                top = max(ex, key=lambda r: (r["ex_bp"], r["anchor"]))
                rd["exchange_bp"] += top["ex_bp"]
                rd["exchange_segments"] += [
                    f"{top['anchor'].split('_')[0]}:{s / 1e6:.1f}-{e / 1e6:.1f}Mb({p})~chr{partner:02d}"
                    for s, e, p in top["ex_segs"]]
                rd["exchange_support"].append(f"chr{partner:02d}:{len(ex)}/{len(mine)}")
    out = {}
    for chrom, rd in readings.items():
        out[chrom] = dict(
            residual_bp=rd["residual_bp"], residual_frac=min(rd["residual_frac"], 1.0),
            partner=rd["partners"][0], partners=",".join(f"chr{p:02d}" for p in rd["partners"]),
            residual_support=";".join(rd["support"]),
            controlled="yes" if controls else "no",
            residual_segments=";".join(rd["segments"]),
            terminal_bp=rd["terminal_bp"],
            exchange_bp=rd["exchange_bp"], exchange_segments=";".join(rd["exchange_segments"]),
            exchange_support=";".join(rd["exchange_support"]))
    return rows, out, null
