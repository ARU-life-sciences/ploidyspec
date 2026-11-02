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

Distances are read from one copy's windows only (-ln containment / k), so
sequence missing from the other copy reads as divergence: daPilAura1 HAP1
chr01 reads 0.015 from HAP2 chr01, HAP2 reads 0.0095 back. A pair of copies is
measured by the lower of its two directions.

In residual tetrasomy every copy is close to every other, so a window also
needs the anchor's own homolog within the same cut as the homeolog (2026-10-22;
before, within the uncapped allelic cut). In daPilAura1 the homeologs are as
close as the alleles (median 0.0115 against 0.0023-0.0197), the ceiling holds
the homeolog cut at 0.0025, and the old own cut of 3 x 0.011 let through
windows where a homeolog copy was closer than the anchor's own allele -- 79 Mb
of "residual tetrasomy" mostly carried by three of the four copies. A window close to a
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
one copy -- `exchange_support` says how many. The own copies must also show
the exchange from their side (REVERSE_FRAC): an exchange replaced sequence
they still carry, so they are far from the anchor over a comparable stretch;
a stretch merely missing from an own copy's assembly reads far one way only
(charr chr09, masu chr18), and is reported as `exchange_one_way_bp` instead.

Low-complexity windows (satellite arrays and repeat-dense stretches: zlib
compression ratio below 0.75 of the genome median) are not read for either test: homeologous centromeres share
a satellite family, which reads as near-allelic in every copy (LOW_COMPLEXITY_FRAC).

Needs `windowed` and `windowed-homeologs` outputs; returns nothing without them.
"""

import csv
import os
import random
import statistics
import zlib
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
# Satellite arrays and repeat-dense windows: a window whose sequence compresses
# to under this fraction of the genome's median compression ratio is not read.
# Position-free containment cannot tell shared repeats from shared sequence
# there: daLatClan1's chr03/chr04 "residual" segment (29-40 Mb, zlib 0.19-0.27
# against a median 0.29) read 0.0005 between homeologs, but they align at
# only 98-99% against 99.9-100% for the alleles. At 0.75 that call goes while
# daPilAura1's (homeolog identity 99.6-99.8%, as its alleles) and masu's stay.
# daGleHede1's chr14/chr17 "residual" block (16-22 Mb) is a centromeric
# satellite: k-mer containment ~1 against both homeologs (shared satellite
# family) but alignment identity 94-96%, and zlib ratio 0.06-0.13 against
# 0.32 for the rest of the chromosome. Homeologous centromeres inherit one
# satellite family, so the unrelated-chromosome control does not remove it.
LOW_COMPLEXITY_FRAC = 0.75
EXCHANGE_RATIO = 2.0  # homeolog copy this much closer than the anchor's own copies
HOMEOLOG_FRAC = 0.25  # near-allelic also means well below the typical homeolog distance
MIN_NULL_RUN = 8  # floor on the shuffled-run null, in windows
# An exchange must show from the own copy's side too: at least this fraction of
# each segment's exchange windows, as far stretches of each own copy read back
# against the anchor at the corresponding place (reverse_far_windows). Containment is read from one
# copy's windows only, so sequence missing from the other copy reads as
# divergence: charr HAP1 chr09's first 6 Mb are 0.07-0.09 from HAP2 chr09
# (51 against 62 Mb assembled), but HAP2 chr09 reads 0.0006 back.
REVERSE_FRAC = 0.5
REVERSE_MIN_RUN = 3
REVERSE_SLACK = 0.1  # fraction of the chromosome length the corresponding place may shift

FIELDS = ["chrom", "partner", "anchor", "control", "n_windows", "n_near_allelic", "n_shared_with_control",
          "allele_level",
          "longest_run", "null_max_run", "segments", "covered_bp", "covered_frac",
          "n_exchange", "exchange_segments", "exchange_bp", "exchange_reverse_windows", "exchange_one_way_bp",
          "n_low_complexity"]


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


def own_pair_tracks(path):
    """{(anchor, other copy): {win_start: distance}} from the windowed (same
    number) tracks: each anchor window's distance to each other copy."""
    out = defaultdict(dict)
    if not os.path.exists(path):
        return out
    for r, d in window_rows(path):
        out[(r["unit_a"], r["unit_b"])][int(r["win_start"])] = d
    return out


def allele_levels(pairs, chrom_of=None):
    """{anchor: allelic level} from own_pair_tracks: per pair of copies, the
    lower of the two directions' median window distance (sequence missing
    from one copy inflates only the direction read from the other; see the
    module docstring), and per anchor its closest copy. With chrom_of, every
    copy of a chromosome gets the median of those levels, so a divergent copy
    is measured against its chromosome's allelic level rather than its own
    distance to the rest."""
    med = {k: statistics.median(v.values()) for k, v in pairs.items() if v}
    per = defaultdict(list)
    for (a, b), m in med.items():
        per[a].append(min(m, med.get((b, a), m)))
    levels = {a: min(v) for a, v in per.items()}
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


def own_copy_distances(pairs):
    """{(anchor, win_start): distance to the anchor's closest own-number copy}"""
    out = {}
    for (a, _), track in pairs.items():
        for s, d in track.items():
            out[(a, s)] = min(d, out.get((a, s), d))
    return out


def reverse_far_windows(pairs, anchor, seg, anchor_len, floor, low=None, min_run=REVERSE_MIN_RUN,
                        slack=REVERSE_SLACK):
    """The fewest windows, over the anchor's own copies, that are beyond `floor`
    from the anchor read from that copy's side, in runs of at least `min_run`
    (bridged as in runs()), at the place on that copy corresponding to the
    anchor's segment `seg` = (start, end): the same fraction of its length
    +- `slack`, in either orientation. A real exchange replaced sequence the own
    copies still carry there; sequence merely missing from an own copy leaves
    that copy close to the anchor, and a far stretch elsewhere on it is another
    event (OncMaso1_wu2026 chr18: HAP1's end and HAP2's start each missing from
    the other). None when no own copy has a track back to the anchor."""
    f0, f1 = seg[0] / anchor_len, seg[1] / anchor_len
    spans = [(f0 - slack, f1 + slack), (1 - f1 - slack, 1 - f0 + slack)]
    counts = []
    for (b, a), track in pairs.items():
        if a != anchor or not track:
            continue
        starts = sorted(track)
        step = min((y - x for x, y in zip(starts, starts[1:])), default=1)
        length = starts[-1] + step
        skip = (low or {}).get(b, set())
        flags = [track[x] > floor and x not in skip for x in starts]
        best = 0
        for lo, hi in spans:
            n = 0
            for i, j, _ in runs(flags):
                if j - i + 1 < min_run:
                    continue
                n += sum(1 for k in range(i, j + 1) if flags[k] and lo <= starts[k] / length <= hi)
            best = max(best, n)
        counts.append(best)
    return min(counts) if counts else None


def compress_ratio(seq):
    return len(zlib.compress(seq, 1)) / max(len(seq), 1)


def window_complexity(outdir, windows_by_anchor):
    """{anchor: {win_start: zlib compression ratio}} for every anchor's
    windows, from its chroms/ FASTA; cached in windowed/window_complexity.tsv.
    Anchors without a FASTA are left out (no guard for them)."""
    from .kmer_tables import chrom_fasta_path

    path = os.path.join(windowed_dir(outdir), "window_complexity.tsv")
    cached = defaultdict(dict)
    if os.path.exists(path):
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                cached[r["anchor"]][int(r["win_start"])] = float(r["ratio"])
    out, new = {}, []
    for anchor, starts in windows_by_anchor.items():
        if anchor in cached and all(s in cached[anchor] for s, _ in starts):
            out[anchor] = cached[anchor]
            continue
        fa = chrom_fasta_path(outdir, anchor)
        if not os.path.exists(fa):
            continue
        with open(fa, "rb") as f:
            seq = b"".join(line.strip() for line in f if not line.startswith(b">")).upper()
        out[anchor] = {s: compress_ratio(seq[s - 1:e]) for s, e in starts}
        new += [(anchor, s, out[anchor][s]) for s, _ in starts]
    if new:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fresh = not os.path.exists(path)
        with open(path, "a", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            if fresh:
                w.writerow(["anchor", "win_start", "ratio"])
            w.writerows((a, s, f"{r:.4f}") for a, s, r in new)
    return out


def low_complexity(complexity):
    """{anchor: set of low-complexity window starts}"""
    values = [r for per in complexity.values() for r in per.values()]
    if not values:
        return {}
    cut = LOW_COMPLEXITY_FRAC * statistics.median(values)
    return {a: {s for s, r in per.items() if r < cut} for a, per in complexity.items()}


def exchange_track(anchor, windows, own, allele, ratio, low=()):
    """[(start, end, exchange)]: windows closer to a homeolog copy than to the
    anchor's own copies (see module docstring)."""
    floor = ratio * max(allele, ALLELE_FLOOR)
    out = []
    for s in sorted(windows):
        d_own = own.get((anchor, s))
        d_hom = min(windows[s][1].values())
        hit = (d_own is not None and d_own > floor and d_hom < d_own / EXCHANGE_RATIO
               and s not in low)
        out.append((s, windows[s][0], hit))
    return out


def near_allelic_track(windows, allele, ratio, control=None, ceiling=None, own=None, low=()):
    """[(start, end, near_allelic, shared_with_control)] along one anchor. A
    window close to the control chromosome as well is not near-allelic; the
    cut never exceeds `ceiling` (a fraction of the homeolog distance). With
    `own` ({start: distance to the anchor's own other copies}) a window whose
    own homolog is beyond the same cut is not near-allelic either: a homeolog
    copy closer than the anchor's own is an exchange-like window
    (exchange_track), not residual tetrasomy."""
    cut = ratio * max(allele, ALLELE_FLOOR)
    if ceiling is not None:
        cut = min(cut, ceiling)
    out = []
    for s in sorted(windows):
        near = min(windows[s][1].values()) < cut
        if near and own and own.get(s) is not None and own[s] > cut:
            near = False
        ctl = control.get(s) if control else None
        shared = bool(ctl and min(ctl[1].values()) < cut)
        out.append((s, windows[s][0], near and not shared and s not in low, near and shared))
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
    own_pairs = own_pair_tracks(os.path.join(windowed_dir(outdir), "windowed_all.tsv"))
    alleles = allele_levels(own_pairs, chrom_of)
    controls = load_tracks(os.path.join(homeologs_dir(outdir), CONTROLS_TSV))
    control_of = {anchor: (c, w) for (anchor, c), w in controls.items()}
    ratio = thresholds.get("residual_ratio", RESIDUAL_RATIO)
    min_bp = thresholds["min_segment_bp"]
    homeolog_level = statistics.median(min(d.values()) for w in tracks.values() for _, d in w.values())
    ceiling = thresholds.get("residual_homeolog_frac", HOMEOLOG_FRAC) * homeolog_level

    own = own_copy_distances(own_pairs)
    starts = defaultdict(dict)
    for (anchor, _), windows in tracks.items():
        starts[anchor].update({s: windows[s][0] for s in windows})
    low = low_complexity(window_complexity(outdir, {a: sorted(v.items()) for a, v in starts.items()}))
    built = {}
    for (anchor, partner), windows in sorted(tracks.items()):
        if anchor not in alleles:
            continue
        built[(anchor, partner)] = near_allelic_track(
            windows, alleles[anchor], ratio, control_of.get(anchor, (None, None))[1], ceiling,
            {s: own.get((anchor, s)) for s in windows}, low.get(anchor, set()))
    lengths = [len(t) for t in built.values()]
    n_true = sum(t[2] for track in built.values() for t in track)
    null = null_max_run(lengths, n_true, thresholds.get("residual_permutations", RESIDUAL_PERMUTATIONS),
                        RESIDUAL_NULL_QUANTILE)
    null = max(null, thresholds.get("residual_min_run", MIN_NULL_RUN))

    ex_tracks = {key: exchange_track(key[0], tracks[key], own, alleles[key[0]], ratio,
                                     low.get(key[0], set())) for key in built}
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

    # each exchange segment against the far stretch its own copies show back at
    # the corresponding place (REVERSE_FRAC); one-way segments are dropped
    ex_segments, reverse, one_way_bp = {}, {}, {}
    for key, t in ex_tracks.items():
        kept, counts, dropped = [], [], 0
        floor = ratio * max(alleles[key[0]], ALLELE_FLOOR)
        for seg in segments(t, 2, ex_null)[1]:
            n = sum(1 for x, _, hit in t if hit and seg[0] <= x <= seg[1])
            back = reverse_far_windows(own_pairs, key[0], seg, t[-1][1], floor, low)
            counts.append("" if back is None else str(back))
            if back is not None and back < REVERSE_FRAC * n:
                dropped += seg[1] - seg[0] + 1
            else:
                kept.append(seg)
        ex_segments[key], reverse[key], one_way_bp[key] = kept, ";".join(counts), dropped

    rows, reads_by_pair = [], defaultdict(list)
    for (anchor, partner), track in built.items():
        all_runs, segs, span = segments(track, 2, null)
        ex_segs = ex_segments[(anchor, partner)]
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
                         exchange_bp=ex_bp,
                         exchange_reverse_windows=reverse[(anchor, partner)],
                         exchange_one_way_bp=one_way_bp[(anchor, partner)],
                         n_low_complexity=sum(1 for t in track if t[0] in low.get(anchor, ()))))
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
