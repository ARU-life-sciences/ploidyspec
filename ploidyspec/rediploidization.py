"""
Rediploidization: how far, and where, a polyploid genome has returned to
diploid-like (disomic) structure. Combines three lines of evidence that the
other stages already compute, plus one new k-mer test:

1. Chromosome fusions between haplotype copies (the new test). A fusion in
   only some haplotypes splits a tetrasomic quartet into a fused and an unfused
   lineage that can no longer pair freely -- the mechanism Xie et al. 2026
   documented in the snow carps. Detected by k-mer *containment*: a fused
   scaffold contains most of the k-mers of *two* chromosome numbers, where an
   ordinary chromosome contains only background levels of any other. Tested on
   (a) placed copies much longer than their siblings (reusing the matrix stage's
   shared-k-mer counts, no new k-mer work) and (b) unplaced chromosome-length
   scaffolds whose headers matched no chromosome number (new k-mer tables, built
   only for those few sequences). A placed copy's partner must also be enriched
   over the copy's unfused siblings, which rules out shared ancient homeology
   and siblings that are only fragments.

2. Lineage structure within each chromosome number (>=3 copies): copies split
   into two groups by whole-chromosome distance (same bipartition as
   te_marker_fraction_by_lineage.tsv), read three ways -- whole-chromosome
   distance split, TE-marker split (if te-markers has run), and *where* along the
   chromosome the split holds, from the windowed track (whole chromosome vs one
   region, e.g. SchCurv1 chr19 vs chr17).

3. Ancient (WGD-derived) pairing between different chromosome numbers, from the
   homeologs stage, with distance_ratio as the resolution measure. Pairs that a
   detected fusion explains are relabelled, not counted as WGD pairs.
   Chromosomes with < 3 copies are pooled with their homeolog partner, or,
   failing that, with their reciprocal best match across a significant
   genome-wide bipartition from the structure stage.

4. Residual tetrasomy between homeologs (residual.py): stretches where a
   homeolog pair's copies are still as close as alleles, from the
   windowed-homeologs tracks. A pooled pair split into lineages on the whole
   chromosome but keeping such stretches reads as partially_resolved.

The state labels are descriptive readings of sequence data from one
individual -- not proof of inheritance mode, which needs segregation data.
Thresholds are provisional (tuned on the snow carp anchors) and exposed as
CLI flags; see OUTPUTS.md.
"""

import csv
import itertools
import math
import os
import re
import statistics
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from .common import (
    homeologs_dir,
    log,
    matrix_dir,
    rediploidization_dir,
    run_logex_batch,
    subgenomes_dir,
    sum_hist_distinct,
    window_rows,
    windowed_dir,
)
from .kmer_tables import build_one, ktab_prefix_path, load_sequences
from .homeology_map import BLOCKS_TSV, read_summary as read_map_summary
from .residual import FIELDS as RESIDUAL_FIELDS, residual_tetrasomy
from .subgenome_report import bipartition_by_distance, load_whole_chrom_distances

DEFAULT_LONG_RATIO = 1.4
DEFAULT_ORPHAN_MIN_FRAC = 0.5
DEFAULT_CONTAINMENT_Z = 10.0
DEFAULT_SIBLING_EXCESS = 2.5
DEFAULT_PARTITION_Z = 10.0
PARTITION_MIN_GROUP_FRAC = 0.25
MIN_FRAC_OF_TOP = 0.25
DEFAULT_DIST_SPLIT = 1.25
DEFAULT_TE_SPLIT = 2.0
DEFAULT_WINDOW_SPLIT = 1.25
DEFAULT_MIN_SEGMENT_BP = 2_000_000


# --- fusion detection -------------------------------------------------------


def robust_background(values):
    """(median, scaled MAD) of a list of containment values; MAD floored so a
    near-constant background can't make every small deviation look significant."""
    med = statistics.median(values)
    mad = statistics.median(abs(v - med) for v in values) * 1.4826
    return med, max(mad, 1e-3)


def containment_components(containments, background, z_min):
    """
    Chromosome numbers whose k-mers are contained in a scaffold far above the
    background rate. containments: {chrom: containment of that chromosome's
    best-matching copy in the scaffold}. Returns [(chrom, containment, z)],
    strongest first. Requires a robust z >= z_min, at least twice the
    background median, and at least MIN_FRAC_OF_TOP of the strongest hit. The
    last one matters when the background is near zero (little shared repeat
    content): the MAD then sits at its floor and repeat k-mers shared by every
    chromosome would otherwise clear the z test everywhere (caught on
    simulated genomes, where the shared-TE background is ~0.02 vs ~0.58 for
    the real fusion partner).
    """
    med, mad = background
    top = max(containments.values(), default=0.0)
    hits = []
    for chrom, c in containments.items():
        z = (c - med) / mad
        if z >= z_min and c >= 2 * med and c >= MIN_FRAC_OF_TOP * top:
            hits.append((chrom, c, z))
    return sorted(hits, key=lambda h: -h[1])


def long_copy_outliers(units, long_ratio):
    """Placed copies at >= long_ratio x the median length of their same-chromosome
    siblings -- the length signature of a fusion. Needs >=1 sibling."""
    by_chrom = defaultdict(list)
    for u in units:
        by_chrom[u["chrom"]].append(u)
    out = []
    for chrom, copies in by_chrom.items():
        for u in copies:
            sibs = [int(s["length"]) for s in copies if s["unit_id"] != u["unit_id"]]
            if sibs and int(u["length"]) >= long_ratio * statistics.median(sibs):
                out.append(u)
    return out


def load_pair_counts(outdir):
    """{(unit_a, unit_b): (shared, kmers_a, kmers_b, chosen_k)} from whole_chrom_pairs.tsv,
    keyed both ways round."""
    path = os.path.join(matrix_dir(outdir), "whole_chrom_pairs.tsv")
    counts = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            s, na, nb = int(r["shared"]), int(r["kmers_a"]), int(r["kmers_b"])
            k = int(r["chosen_k"])
            counts[(r["unit_a"], r["unit_b"])] = (s, na, nb, k)
            counts[(r["unit_b"], r["unit_a"])] = (s, nb, na, k)
    return counts


def placed_background(units, counts, exclude):
    """Containment of every copy in every *different-chromosome* copy (smaller in
    larger), skipping units in `exclude` -- the species' no-fusion baseline."""
    vals = []
    ids = [u["unit_id"] for u in units if u["unit_id"] not in exclude]
    chrom = {u["unit_id"]: u["chrom"] for u in units}
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            if chrom[a] == chrom[b] or (a, b) not in counts:
                continue
            s, na, nb, _ = counts[(a, b)]
            vals.append(s / max(min(na, nb), 1))
    return vals


def haps_by_chrom(units):
    out = defaultdict(set)
    for u in units:
        out[u["chrom"]].add(u["hap"])
    return out


def fusion_row(scaffold, hap, length, origin, own_chrom, hits, chrom_haps, k):
    """One fusions.tsv row. A fusion needs >=2 component chromosome numbers; it's
    'supported' when every non-own component is missing from this haplotype
    (the partner was absorbed, not duplicated)."""
    components = [h[0] for h in hits]
    if own_chrom is not None and own_chrom not in components:
        components = [own_chrom] + components
    partners = [c for c in components if c != own_chrom]
    if len(components) < 2:
        status = "no_fusion_signal"
    elif hap and all(hap not in chrom_haps.get(c, set()) for c in partners):
        status = "fusion"
    else:
        status = "candidate_partner_present"
    return dict(
        scaffold=scaffold,
        hap=hap,
        length=length,
        origin=origin,
        components="+".join(f"chr{c:02d}" for c in sorted(components)),
        containment=",".join(f"chr{c:02d}:{v:.3f}" for c, v, _ in hits),
        lineage_divergence=",".join(f"chr{c:02d}:{lineage_divergence(v, k):.4f}" for c, v, _ in hits),
        max_z=f"{max((z for _, _, z in hits), default=0):.1f}",
        status=status,
    )


def lineage_divergence(containment, k):
    """Mash-style per-base divergence between the fused lineage and the unfused
    copies of a component chromosome, from k-mer containment: -ln(c)/k. Ranks
    fusions by how long the fused lineage has been diverging (SchYoun1: chr19+22
    ~0.053, chr11+14 ~0.015, matching the order of Xie et al. 2026's waves)."""
    return -math.log(containment) / k if containment > 0 else float("inf")


def own_fraction(counts, a, b):
    """Fraction of a's k-mers shared with b. Normalised by a's own k-mer count, so
    a fragmentary copy reads the same per-k-mer rate as a complete one."""
    s, na, _, _ = counts[(a, b)]
    return s / max(na, 1)


def sibling_excess(units, counts, u, partner_chrom):
    """How enriched partner_chrom is in copy u relative to u's least-enriched
    sibling copy (best partner copy each time). A fused copy carries the partner,
    an unfused sibling only shares background or ancient-homeolog k-mers:
    SchCurv1 HAP3_chr19 has chr22 at 3.4x its unfused siblings, while
    drMyrSpic1 HAP1_chr13 (a full copy next to a 15 Mb fragment) has its
    homeologs chr03/chr09 at 0.9x. None when u has no sibling to compare."""
    partners = [v["unit_id"] for v in units if v["chrom"] == partner_chrom]

    def best(a):
        vals = [own_fraction(counts, a, p) for p in partners if (a, p) in counts]
        return max(vals) if vals else None

    sibs = [best(v["unit_id"]) for v in units
            if v["chrom"] == u["chrom"] and v["unit_id"] != u["unit_id"]]
    sibs = [s for s in sibs if s is not None]
    own = best(u["unit_id"])
    if own is None or not sibs:
        return None
    return own / max(min(sibs), 1e-9)


def detect_placed_fusions(units, counts, long_ratio, z_min, min_excess=DEFAULT_SIBLING_EXCESS):
    long_units = long_copy_outliers(units, long_ratio)
    if not long_units:
        return []
    exclude = {u["unit_id"] for u in long_units}
    bg_vals = placed_background(units, counts, exclude)
    if len(bg_vals) < 10:
        return []
    background = robust_background(bg_vals)
    chrom_haps = haps_by_chrom(units)
    k = max(c[3] for c in counts.values())
    rows = []
    for u in long_units:
        best = {}
        for v in units:
            if v["chrom"] == u["chrom"] or (v["unit_id"], u["unit_id"]) not in counts:
                continue
            s, nv, _, _ = counts[(v["unit_id"], u["unit_id"])]
            best[v["chrom"]] = max(best.get(v["chrom"], 0.0), s / max(nv, 1))
        hits = containment_components(best, background, z_min)
        excess = {c: sibling_excess(units, counts, u, c) for c, _, _ in hits}
        kept = [h for h in hits if excess[h[0]] is None or excess[h[0]] >= min_excess]
        row = fusion_row(u["unit_id"], u["hap"], int(u["length"]), "placed_long_copy",
                         u["chrom"], kept, chrom_haps, k)
        row["sibling_excess"] = ",".join(
            f"chr{c:02d}:{excess[c]:.2f}" for c, _, _ in hits if excess[c] is not None)
        rows.append(row)
    return rows


def orphan_candidates(outdir, units, min_len, orphan_min_frac):
    """Unplaced sequences long enough to be whole (possibly fused) chromosomes."""
    path = os.path.join(outdir, "unplaced.tsv")
    if not os.path.exists(path):
        return []
    median_len = statistics.median(int(u["length"]) for u in units)
    cutoff = max(min_len, orphan_min_frac * median_len)
    hap_of_source = defaultdict(set)
    for u in units:
        hap_of_source[u["source"]].add(u["hap"])
    out = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["reason"] == "below-min-len" or int(r["length"]) < cutoff:
                continue
            haps = hap_of_source.get(r["source"], set())
            out.append(
                dict(
                    unit_id="ORPHAN_" + re.sub(r"[^A-Za-z0-9]+", "_", r["seq_id"]).strip("_"),
                    seq_id=r["seq_id"],
                    source=r["source"],
                    length=r["length"],
                    hap=next(iter(haps)) if len(haps) == 1 else "",
                )
            )
    return out


def orphan_containments(outdir, orphan, units, k, tools, threads):
    """{chrom: best containment of that chromosome's copies in the orphan scaffold},
    cached in rediploidization/orphan_containment.tsv so reruns skip the Logex work."""
    samtools_bin, fastk_bin, logex_bin, histex_bin = tools
    cache = os.path.join(rediploidization_dir(outdir), "orphan_containment.tsv")
    cached = defaultdict(dict)
    if os.path.exists(cache):
        with open(cache) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                cached[r["orphan"]][r["unit_id"]] = float(r["containment"])
    per_unit = cached.get(orphan["unit_id"])
    if not per_unit:
        build_one(samtools_bin, fastk_bin, orphan, k, outdir)
        for u in units:  # rebuilds tables removed by --cleanup; no-op otherwise
            build_one(samtools_bin, fastk_bin, u, k, outdir)
        o_prefix = ktab_prefix_path(outdir, orphan["unit_id"], k)
        per_unit = {}

        def do_batch(i, batch):
            prefixes = [o_prefix] + [ktab_prefix_path(outdir, u["unit_id"], k) for u in batch]
            pairs = [(0, j + 1) for j in range(len(batch))]
            tmp = os.path.join(outdir, "tmp_rediploidization", f"{orphan['unit_id']}_{i}")
            shared = run_logex_batch(logex_bin, histex_bin, prefixes, pairs, tmp)
            return {
                u["unit_id"]: shared[(0, j + 1)] / max(sum_hist_distinct(histex_bin, prefixes[j + 1]), 1)
                for j, u in enumerate(batch)
            }

        chunk = 7
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [ex.submit(do_batch, i, units[i : i + chunk]) for i in range(0, len(units), chunk)]
            for fut in as_completed(futs):
                per_unit.update(fut.result())
        new_file = not os.path.exists(cache)
        with open(cache, "a", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            if new_file:
                w.writerow(["orphan", "unit_id", "k", "containment"])
            for uid, c in per_unit.items():
                w.writerow([orphan["unit_id"], uid, k, f"{c:.6f}"])
    chrom_of = {u["unit_id"]: u["chrom"] for u in units}
    best = {}
    for uid, c in per_unit.items():
        if uid in chrom_of:
            best[chrom_of[uid]] = max(best.get(chrom_of[uid], 0.0), c)
    return best


def detect_orphan_fusions(outdir, units, orphans, k, tools, threads, z_min):
    chrom_haps = haps_by_chrom(units)
    rows = []
    for o in orphans:
        log(f"  orphan scaffold {o['seq_id']} ({int(o['length']) / 1e6:.1f} Mb): k-mer containment vs {len(units)} copies")
        best = orphan_containments(outdir, o, units, k, tools, threads)
        if len(best) < 3:
            continue
        hits = containment_components(best, robust_background(list(best.values())), z_min)
        rows.append(
            fusion_row(o["seq_id"], o["hap"], int(o["length"]), "unplaced_scaffold", None, hits, chrom_haps, k)
        )
    return rows


# --- lineage structure along each chromosome ---------------------------------


def load_windowed(outdir):
    """{chrom_label: {anchor_unit: {win_start: (win_end, {other_unit: distance})}}}
    -- each copy's windows along its own coordinates, with the distance of each
    window to every other copy (see windowed.py)."""
    path = os.path.join(windowed_dir(outdir), "windowed_all.tsv")
    out = defaultdict(lambda: defaultdict(dict))
    if not os.path.exists(path):
        return out
    for r, d in window_rows(path):
        # older runs wrote a bare `chrom` number instead of the `group` label
        label = r["group"] if "group" in r else f"chr{int(r['chrom']):02d}"
        entry = out[label][r["unit_a"]].setdefault(int(r["win_start"]), (int(r["win_end"]), {}))
        entry[1][r["unit_b"]] = d
    return out


def window_split_track(windows, anchor, group_a):
    """[(start, end, cross/within ratio)] along one anchor copy, for a fixed
    bipartition: mean distance from the anchor's window to the other group's
    copies over mean distance to its own group's other copies. Empty when the
    anchor has no other copy in its own group (the lone copy of a 1-vs-rest
    split)."""
    own = set(group_a) if anchor in group_a else None
    track = []
    for start in sorted(windows):
        end, dists = windows[start]
        if own is None:
            own_group = {u for u in dists} - set(group_a)
        else:
            own_group = own
        within = [d for u, d in dists.items() if u in own_group and u != anchor]
        cross = [d for u, d in dists.items() if u not in own_group]
        if within and cross:
            w = statistics.mean(within)
            track.append((start, end, statistics.mean(cross) / w if w > 0 else float("inf")))
    return track


def split_segments(track, threshold, min_segment_bp, max_gap=1):
    """Contiguous runs of split windows (ratio >= threshold), bridging gaps of up
    to max_gap windows, kept if >= min_segment_bp long. Returns [(start, end)]."""
    segments = []
    run_start = run_end = None
    gap = 0
    for start, end, ratio in track:
        if ratio >= threshold:
            if run_start is None:
                run_start = start
            run_end = end
            gap = 0
        elif run_start is not None:
            gap += 1
            if gap > max_gap:
                segments.append((run_start, run_end))
                run_start = None
                gap = 0
    if run_start is not None:
        segments.append((run_start, run_end))
    return [(s, e) for s, e in segments if e - s + 1 >= min_segment_bp]


def load_te_split(outdir):
    path = os.path.join(subgenomes_dir(outdir), "te_marker_fraction_by_lineage.tsv")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {
            r["chrom"]: float(r["split_ratio"])
            for r in csv.DictReader(f, delimiter="\t")
            if r["split_ratio"]
        }


def copy_state(n_copies, balanced, dist_split, te_split, extent, thresholds):
    """
    Within-chromosome-number state, from the three lineage readings.
    resolved_lineages: balanced split, strong whole-chromosome distance split, and
      the split holds along most of the chromosome (e.g. SchCurv1 chr19).
    partially_resolved: >=2 of {distance split, TE split, windowed segment}
      (e.g. SchCurv1 chr17 -- regional split, strong TE signal).
    candidate: exactly one line of evidence.
    one_divergent_copy: the split isolates one copy from the rest. Either one
      odd haplotype (assembly quality -- check whether the same haplotype recurs,
      see most_frequent_outlier_hap) or a genuinely divergent genome copy (AAAB-
      like; every split of a triploid looks like this). Not two lineages.
      Needs chromosome-wide evidence (distance or TE split): a windowed segment
      alone isolating one copy is ordinary haplotype structure in a polysomic
      genome -- one of four copies carrying a divergent block somewhere
      (drLytSali1, tetrasomic by classical genetics, had 12/15 chromosomes
      called one_divergent_copy on regional segments alone).
    """
    if n_copies < 3:
        return "not_assessable"
    support = [
        dist_split is not None and dist_split >= thresholds["dist_split"],
        te_split is not None and te_split >= thresholds["te_split"],
        extent in ("whole", "regional"),
    ]
    if not balanced:
        return "one_divergent_copy" if support[0] or support[1] else "tetrasomic_like"
    if dist_split is not None and dist_split >= 2 * thresholds["dist_split"] and extent == "whole":
        return "resolved_lineages"
    n = sum(support)
    return {0: "tetrasomic_like", 1: "candidate"}.get(n, "partially_resolved")


def split_lineages(ids, dist):
    """2-way split of copies by whole-chromosome distance. Returns
    (group_a, group_b, balanced, dist_split) or None when < 3 copies. dist_split
    is mean cross-group / mean within-group distance; balanced = both groups
    have >= 2 copies (a 1-vs-rest split isolates one odd copy, not a lineage)."""
    split = bipartition_by_distance(sorted(ids), dist)
    if not split:
        return None
    ga, gb = split
    within = [dist[(x, y)] for g in (ga, gb) for i, x in enumerate(g) for y in g[i + 1 :]]
    cross = [dist[(x, y)] for x in ga for y in gb]
    ratio = None
    if within and statistics.mean(within) > 0:
        ratio = statistics.mean(cross) / statistics.mean(within)
    return ga, gb, min(len(ga), len(gb)) >= 2, ratio


def pooled_state(balanced, dist_split, thresholds):
    """
    State for a homeolog-pooled group: the copies of two chromosome numbers that
    the homeologs stage paired, read as one group. This is how a two-haplotype
    assembly of a tetraploid is assessed -- each number has only 2 copies, but a
    number plus its homeolog has 4. Only whole-chromosome distance is available
    across different numbers (windowed and TE-marker tracks compare same-number
    copies), so one line of evidence decides:
    resolved_lineages: the two numbers' copies form separate lineages,
      dist_split >= 2x threshold (allopolyploid subgenomes, or a long-diploidized
      autopolyploid).
    candidate: dist_split between 1x and 2x threshold.
    tetrasomic_like: the homeolog is as close as the same-number copy --
      the copies are interchangeable across the two numbers.
    """
    if dist_split is None:
        return "not_assessable"
    t = thresholds["dist_split"]
    if not balanced:
        return "one_divergent_copy" if dist_split >= t else "tetrasomic_like"
    if dist_split >= 2 * t:
        return "resolved_lineages"
    return "candidate" if dist_split >= t else "tetrasomic_like"


WINDOW_NULL_QUANTILE = 0.95


def windowed_reading(windows_by_anchor, group_a, thresholds):
    """The windowed split for one grouping of a chromosome's copies, read along
    every copy that has another copy in its own group. Returns the median copy's
    (covered fraction, split-window fraction, anchor, segments), with the
    split-window fraction as the median over copies, or None if no copy has a
    track."""
    reads = []
    for anchor, windows in sorted(windows_by_anchor.items()):
        track = window_split_track(windows, anchor, group_a)
        if not track:
            continue
        n_split = sum(1 for _, _, r in track if r >= thresholds["window_split"])
        segs = split_segments(track, thresholds["window_split"], thresholds["min_segment_bp"])
        span = track[-1][1] - track[0][0] + 1
        covered = sum(e - s + 1 for s, e in segs) / span
        reads.append((covered, n_split / len(track), anchor, segs))
    if not reads:
        return None
    # the median anchor: a real lineage split is seen from every copy
    covered, _, anchor, segs = sorted(reads)[(len(reads) - 1) // 2]
    return covered, statistics.median(r[1] for r in reads), anchor, segs


def crossing_groupings(ids, group_a):
    """Alternative groupings of the same copies that cut across group_a, for the
    windowed null. A balanced split is compared with the other balanced splits;
    a split isolating one copy, with the splits isolating each other copy.
    Groupings that keep the observed groups together would inherit a real
    split's signal (any grouping that keeps a tight pair together looks split)."""
    ids = sorted(ids)
    small = group_a if len(group_a) <= len(ids) - len(group_a) else [x for x in ids if x not in group_a]
    if len(small) == 1:
        return [[x] for x in ids if x != small[0]]
    out = []
    for g in itertools.combinations(ids, len(small)):
        g = list(g)
        if set(g) in (set(small), set(ids) - set(small)) or ids[0] not in g:
            continue
        out.append(g)
    return out


def quantile(values, q):
    v = sorted(values)
    return v[int(q * (len(v) - 1))] if v else 0.0


def chromosome_lineages(units, dist, windowed, te_split, thresholds):
    """
    Per chromosome number: the distance split of its copies, the TE split, and
    where along the chromosome the split holds.

    The windowed extent is tested against a null before it counts. Window
    distances carry real noise, so a run of split windows can appear by chance
    (daBudDavi1: single 2-3 Mb segments on chromosomes whose whole-chromosome
    split is ~1.0). The same scan is run for groupings that cut across the
    observed one (crossing_groupings); their covered fractions, pooled over the
    genome, give the species' noise level. A split counts as `whole`/`regional`
    only when its covered fraction exceeds both that null's 95th percentile and
    every crossing grouping of its own chromosome. `extent_raw` keeps the
    untested reading.
    """
    by_chrom = defaultdict(list)
    for u in units:
        by_chrom[u["chrom"]].append(u["unit_id"])
    rows, null = {}, []
    for chrom, ids in sorted(by_chrom.items()):
        label = f"chr{chrom:02d}"
        row = dict(n_copies=len(ids), group_a="", group_b="", dist_split=None,
                   te_split=te_split.get(label), window_split_frac=None,
                   extent="", extent_raw="", segments="", balanced=False,
                   covered=None, own_null=None)
        split = split_lineages(ids, dist)
        if split:
            ga, gb, balanced, ratio = split
            row.update(group_a=",".join(ga), group_b=",".join(gb), balanced=balanced,
                       dist_split=ratio)
            reading = windowed_reading(windowed.get(label, {}), ga, thresholds)
            if reading:
                covered, frac, anchor, segs = reading
                alts = [windowed_reading(windowed.get(label, {}), g, thresholds)
                        for g in crossing_groupings(ids, ga)]
                alts = [a[0] for a in alts if a]
                null.extend(alts)
                row.update(window_split_frac=frac, covered=covered, own_null=max(alts, default=0.0),
                           extent_raw="whole" if covered >= 0.5 else ("regional" if segs else "none"))
                row["segments"] = ";".join(
                    f"{anchor.split('_')[0]}:{s / 1e6:.1f}-{e / 1e6:.1f}Mb" for s, e in segs)
        rows[chrom] = row
    q = quantile(null, WINDOW_NULL_QUANTILE)
    for row in rows.values():
        row["window_null"] = q
        if row["covered"] is None:
            continue
        passes = row["covered"] > q and row["covered"] > row["own_null"]
        row["extent"] = row["extent_raw"] if passes else "none"
    return rows


# --- ancient pairing ----------------------------------------------------------


def load_ancient_pairs(outdir):
    """{chrom: (partner_chrom, distance_ratio)} for FDR-accepted homeolog pairs."""
    path = os.path.join(homeologs_dir(outdir), "ploidy_ancestry_summary.tsv")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("homeolog_partner"):
                ratio = r.get("distance_ratio")
                out[int(r["chrom"].replace("chr", ""))] = (
                    int(r["homeolog_partner"].replace("chr", "")),
                    float(ratio) if ratio else None,
                )
    return out


def load_chrom_distances(outdir):
    """{(chrom_a, chrom_b): distance} between chromosome numbers (both orders),
    from the homeologs stage's ranked candidates."""
    path = os.path.join(homeologs_dir(outdir), "homeolog_candidates_ranked.tsv")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            a, b = int(r["chrom_a"].replace("chr", "")), int(r["chrom_b"].replace("chr", ""))
            out[(a, b)] = out[(b, a)] = float(r["distance"])
    return out


def reciprocal_matches(group_a, group_b, cdist):
    """{chrom: partner} for chromosomes that are each other's closest match in the
    opposite group (both directions recorded)."""
    def nearest(c, others):
        cands = [(cdist[(c, o)], o) for o in others if (c, o) in cdist]
        return min(cands)[1] if cands else None

    out = {}
    for a in group_a:
        b = nearest(a, group_b)
        if b is not None and nearest(b, group_a) == a:
            out[a], out[b] = b, a
    return out


def load_partition_pairs(outdir, z_min):
    """
    Partners from a diffuse genome-wide bipartition, for species where the
    homeologs stage accepted no pairs but the structure stage found a
    significant 2-group split of chromosome numbers (daInuConz1, dmRanRepe1,
    dcCerAlpi1). Uses the k=2 partition when its z >= z_min and the smaller
    group holds >= PARTITION_MIN_GROUP_FRAC of the chromosomes -- a 2-vs-rest
    split is one tight pair, not two subgenomes. Each chromosome is paired with
    its reciprocal best match in the other group. Returns ({chrom: partner}, z).
    """
    path = os.path.join(outdir, "structure", "genome_partition.tsv")
    if not os.path.exists(path):
        return {}, None
    with open(path) as f:
        row = next((r for r in csv.DictReader(f, delimiter="\t") if r["k"] == "2"), None)
    if row is None or float(row["z_score"]) < z_min:
        return {}, None
    groups = [[int(c.strip().replace("chr", "")) for c in g.split(",")]
              for g in row["groups"].split("|")]
    if len(groups) != 2 or min(map(len, groups)) < PARTITION_MIN_GROUP_FRAC * sum(map(len, groups)):
        return {}, None
    return reciprocal_matches(groups[0], groups[1], load_chrom_distances(outdir)), float(row["z_score"])


# --- driver -------------------------------------------------------------------


def _fmt(v, nd=2):
    return "" if v is None else f"{v:.{nd}f}"


def compute_rediploidization(outdir, k_values, min_len, thresholds, tools_fn, threads):
    seq_tsv = os.path.join(outdir, "sequences.tsv")
    units = load_sequences(seq_tsv)
    for u in units:
        u["chrom"] = int(u["chrom"])
    counts = load_pair_counts(outdir)
    dist = load_whole_chrom_distances(outdir)

    log("rediploidization: fusion detection (placed long copies)")
    fusions = detect_placed_fusions(units, counts, thresholds["long_ratio"], thresholds["containment_z"],
                                    thresholds.get("sibling_excess", DEFAULT_SIBLING_EXCESS))
    orphans = orphan_candidates(outdir, units, min_len, thresholds["orphan_min_frac"])
    if orphans:
        k = max((c[3] for c in counts.values()), default=max(k_values))
        log(f"rediploidization: {len(orphans)} unplaced chromosome-length scaffold(s), k={k}")
        fusions += detect_orphan_fusions(outdir, units, orphans, k, tools_fn(), threads, thresholds["containment_z"])

    n_tested = len(fusions)
    fusions = [f for f in fusions if f["status"] != "no_fusion_signal"]
    log(f"rediploidization: {n_tested} long/unplaced scaffold(s) tested, {len(fusions)} with fusion signal")

    fused_in = defaultdict(lambda: defaultdict(set))  # chrom -> partner -> haps
    for f in fusions:
        if f["status"] != "fusion":
            continue
        comps = [int(c[3:]) for c in f["components"].split("+")]
        for c in comps:
            for p in comps:
                if p != c:
                    fused_in[c][p].add(f["hap"])

    log("rediploidization: lineage structure per chromosome")
    lineages = chromosome_lineages(units, dist, load_windowed(outdir), load_te_split(outdir), thresholds)
    ancient = load_ancient_pairs(outdir)
    partition, partition_z = {}, None
    if not ancient:
        partition, partition_z = load_partition_pairs(
            outdir, thresholds.get("partition_z", DEFAULT_PARTITION_Z))
        if partition:
            log(f"rediploidization: no homeolog pairs; pooling {len(partition)} chromosomes by "
                f"reciprocal match across the k=2 genome partition (z={partition_z:.1f})")
    all_haps = sorted({u["hap"] for u in units})
    chrom_haps = haps_by_chrom(units)
    residual_rows, residual, residual_null = residual_tetrasomy(
        outdir, {u["unit_id"]: u["chrom"] for u in units}, thresholds)
    if residual:
        log(f"rediploidization: residual tetrasomy in {sum(1 for r in residual.values() if r['residual_bp'])}"
            f" of {len(residual)} homeolog-paired chromosomes (null run length {residual_null} windows)")

    ids_by_chrom = defaultdict(list)
    for u in units:
        ids_by_chrom[u["chrom"]].append(u["unit_id"])
    block_partners = load_block_partners(outdir)

    rows = []
    for chrom, lin in sorted(lineages.items()):
        fusion_desc = "; ".join(
            f"with chr{p:02d} in {','.join(sorted(h))}" for p, h in sorted(fused_in[chrom].items())
        )
        partner = ancient.get(chrom)
        if partner is None:
            ancient_state, partner_str, ratio = "unpaired", "", None
        else:
            partner_str, ratio = f"chr{partner[0]:02d}", partner[1]
            ancient_state = "fusion_partner" if partner[0] in fused_in[chrom] else "paired"
        pooled_with = ""
        if fused_in[chrom]:
            state, basis = "fusion_lineages", "fusion"
        elif lin["n_copies"] >= 3:
            state, basis = copy_state(lin["n_copies"], lin["balanced"], lin["dist_split"],
                                      lin["te_split"], lin["extent"], thresholds), "copies"
        else:
            if ancient_state == "paired" and partner[0] in ids_by_chrom:
                mate, pool_basis = partner[0], "homeolog_pool"
            elif partition.get(chrom) in ids_by_chrom:
                mate, pool_basis = partition[chrom], "partition_pool"
            else:
                mate = None
            pooled = split_lineages(ids_by_chrom[chrom] + ids_by_chrom[mate], dist) if mate else None
            if pooled:
                ga, gb, balanced, pooled_ratio = pooled
                lin = dict(lin, group_a=",".join(ga), group_b=",".join(gb),
                           dist_split=pooled_ratio)
                state, basis = pooled_state(balanced, pooled_ratio, thresholds), pool_basis
                pooled_with = f"chr{mate:02d}"
                if state == "resolved_lineages" and residual.get(chrom, {}).get("residual_bp"):
                    state = "partially_resolved"
            else:
                state, basis = "not_assessable", "none"
        outlier = ""
        if state == "one_divergent_copy":
            ga, gb = lin["group_a"].split(","), lin["group_b"].split(",")
            outlier = (ga if len(ga) == 1 else gb)[0].split("_")[0]
        rows.append(
            dict(
                chrom=f"chr{chrom:02d}",
                n_copies=lin["n_copies"],
                haps_missing=",".join(h for h in all_haps if h not in chrom_haps[chrom]),
                fusion=fusion_desc,
                group_a=lin["group_a"],
                group_b=lin["group_b"],
                dist_split=_fmt(lin["dist_split"]),
                te_split=_fmt(lin["te_split"]),
                window_split_frac=_fmt(lin["window_split_frac"]),
                split_extent=lin["extent"],
                split_extent_raw=lin.get("extent_raw", ""),
                window_covered=_fmt(lin.get("covered")),
                window_null_q95=_fmt(lin.get("window_null")),
                split_segments=lin["segments"],
                copy_state=state,
                state_basis=basis,
                pooled_with=pooled_with,
                outlier_hap=outlier,
                ancient_partner=partner_str,
                distance_ratio=_fmt(ratio),
                ancient_state=ancient_state,
                residual_bp=residual[chrom]["residual_bp"] if chrom in residual else "",
                residual_frac=_fmt(residual[chrom]["residual_frac"], 3) if chrom in residual else "",
                residual_terminal_bp=residual[chrom]["terminal_bp"] if chrom in residual else "",
                residual_support=residual.get(chrom, {}).get("residual_support", ""),
                residual_segments=residual.get(chrom, {}).get("residual_segments", ""),
                residual_controlled=residual.get(chrom, {}).get("controlled", ""),
                residual_partners=residual.get(chrom, {}).get("partners", ""),
                exchange_bp=residual[chrom]["exchange_bp"] if chrom in residual else "",
                exchange_support=residual.get(chrom, {}).get("exchange_support", ""),
                exchange_segments=residual.get(chrom, {}).get("exchange_segments", ""),
                block_partners=block_partners.get(chrom, ""),
            )
        )

    rdir = rediploidization_dir(outdir)
    write_tsv(os.path.join(rdir, "fusions.tsv"), fusions,
              ["scaffold", "hap", "length", "origin", "components", "containment",
               "sibling_excess", "lineage_divergence", "max_z", "status"])
    write_tsv(os.path.join(rdir, "rediploidization_by_chrom.tsv"), rows, list(rows[0]) if rows else [])
    if residual_rows:
        write_tsv(os.path.join(rdir, "residual_tetrasomy.tsv"), residual_rows, RESIDUAL_FIELDS)
    summary = summarize(rows, fusions, partition_z, residual_null, read_map_summary(outdir))
    write_tsv(os.path.join(rdir, "rediploidization_summary.tsv"), summary, ["metric", "value"])
    log(f"wrote fusions.tsv, rediploidization_by_chrom.tsv, rediploidization_summary.tsv in {rdir}")
    return rows, fusions, summary


def load_block_partners(outdir):
    """{chrom: "chr05:12.3Mb,chr09:4.0Mb"} from the homeology map's blocks (the
    reference haplotype's copy of each chromosome)."""
    path = os.path.join(homeologs_dir(outdir), BLOCKS_TSV)
    if not os.path.exists(path):
        return {}
    bp = defaultdict(lambda: defaultdict(int))
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            bp[int(r["chrom"][3:])][r["partner"]] += int(r["length"])
    return {c: ",".join(f"{p}:{v / 1e6:.1f}Mb" for p, v in sorted(d.items(), key=lambda x: -x[1]))
            for c, d in bp.items()}


def summarize(rows, fusions, partition_z=None, residual_null=None, homeology=None):
    n = len(rows)
    out = [dict(metric="n_chromosome_numbers", value=n)]
    for basis in ("copies", "homeolog_pool", "partition_pool", "fusion", "none"):
        out.append(dict(metric=f"state_basis:{basis}",
                        value=sum(1 for r in rows if r.get("state_basis") == basis)))
    for state in ("tetrasomic_like", "candidate", "partially_resolved", "resolved_lineages",
                  "fusion_lineages", "one_divergent_copy", "not_assessable"):
        k = sum(1 for r in rows if r["copy_state"] == state)
        out.append(dict(metric=f"copy_state:{state}", value=k))
    outliers = Counter(r["outlier_hap"] for r in rows if r.get("outlier_hap"))
    if outliers:
        hap, k = outliers.most_common(1)[0]
        # one haplotype isolated on most chromosomes points at that haplotype
        # (assembly quality, or a divergent copy) rather than per-chromosome biology
        out.append(dict(metric="most_frequent_outlier_hap", value=f"{hap} ({k} chromosomes)"))
    confirmed = [f for f in fusions if f["status"] == "fusion"]
    out.append(dict(metric="n_distinct_fusions", value=len({f["components"] for f in confirmed})))
    out.append(dict(metric="n_fused_scaffolds", value=len(confirmed)))
    out.append(dict(metric="n_fusion_candidates_partner_present",
                    value=sum(1 for f in fusions if f["status"] == "candidate_partner_present")))
    if partition_z is not None:
        out.append(dict(metric="partition_pool_z", value=f"{partition_z:.2f}"))
    paired = [r for r in rows if r["ancient_state"] == "paired"]
    out.append(dict(metric="ancient_paired_chromosomes", value=len(paired)))
    ratios = [float(r["distance_ratio"]) for r in paired if r["distance_ratio"]]
    if ratios:
        out.append(dict(metric="ancient_distance_ratio_median", value=f"{statistics.median(ratios):.2f}"))
        if len(ratios) >= 2 and statistics.mean(ratios) > 0:
            out.append(dict(metric="ancient_distance_ratio_cv",
                            value=f"{statistics.stdev(ratios) / statistics.mean(ratios):.2f}"))
    tested = [r for r in rows if r.get("residual_bp") not in ("", None)]
    if tested:
        hit = [r for r in tested if r["residual_bp"]]
        bp = sum(r["residual_bp"] for r in hit)
        out += [dict(metric="residual_tested_chromosomes", value=len(tested)),
                dict(metric="residual_chromosomes", value=len(hit)),
                dict(metric="residual_bp", value=bp),
                dict(metric="residual_terminal_frac",
                     value=f"{sum(r['residual_terminal_bp'] for r in hit) / bp:.2f}" if bp else ""),
                dict(metric="residual_null_run_windows", value=residual_null),
                dict(metric="residual_controlled", value=tested[0].get("residual_controlled", ""))]
        ex = [r for r in tested if r.get("exchange_bp")]
        out += [dict(metric="exchange_chromosomes", value=len(ex)),
                dict(metric="exchange_bp", value=sum(r["exchange_bp"] for r in ex)),
                dict(metric="exchange_list", value=";".join(
                    f"{r['chrom']}({r['exchange_support']})" for r in ex))]
    # segmental homeology from the genome-wide map (homeology_map.py)
    for key in ("duplicated_frac", "n_blocks", "n_partner_pairs", "multi_partner_chromosomes",
                "multi_partner_list", "block_dist_median", "block_dist_cv"):
        if homeology and key in homeology:
            out.append(dict(metric=f"map_{key}", value=homeology[key]))
    return out


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
