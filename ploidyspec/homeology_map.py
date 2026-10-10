"""
Genome-wide segmental homeology: which chromosome each stretch of each
chromosome is duplicated on.

`homeologs` pairs whole chromosomes. That misses homeology that runs per arm or
per block -- after an old duplication the karyotype is reshuffled by fusions
and fissions, and a chromosome's two arms can each have a different homeolog
(salmonids: masu and Arctic charr pass only 5-6 whole-chromosome pairs, their
homeologs 0.10 apart against an unrelated background of 0.116-0.124). Here
every window of every chromosome of one haplotype is compared with every
other chromosome of that haplotype (position-free, as in windowed.py), and the
closest one is recorded.

k is 23 (the matrix stage's tables). At the windowed stages' k=15 about 40%
of a window's k-mers turn up by chance in any 100 Mb chromosome, and salmonid
homeologs barely stand out (masu 0.049 vs 0.057 unrelated); at k=23 chance
hits are negligible and what is shared is homeology or repeats.

Per window: the closest chromosome, the median over all others (background)
and a robust z = (median - closest) / (1.4826 MAD). A window is homeolog-like
when z >= MAP_Z. Single windows can pass by chance or through a repeat;
homeology shows as the *same* closest chromosome over consecutive windows.
Blocks are runs of windows with the same closest chromosome (homeolog-like,
gaps of up to MAX_GAP windows bridged), counted when longer than the 95th
percentile of the longest such run with the labels shuffled along each
chromosome (and at least MIN_BLOCK_WINDOWS and `min_segment_bp`).

Sensitivity: exact k-mer matching sees homeology up to roughly 8-10%
divergence. Glechoma's homeologs (0.042) stand far out of the background
(closest/median 0.49 per window) and 78% of its genome is in blocks; the
salmonid duplication is too old for most of the genome (closest/median 0.90,
6-10% in blocks for masu and Arctic charr, 1 Mb windows no better), so only
its least-diverged regions are mapped. An unmapped region is "no homeolog
detected", not "no homeolog".

Readings for rediploidization (rediploidization.py, summary.py):
Homeology is reciprocal: if a stretch of chr A is closest to chr B, a stretch
of B is closest to A. Blocks without a block back (`reciprocal = no`) are kept
in the table but not read: in daGleHede1 two short ones (0.007-0.014, a third
of the homeolog distance) are shared repeats, and the diploid ddMalSylv1's
seven blocks (5% of the genome) all lack a partner block.

- duplicated fraction: share of the genome in reciprocal homeolog blocks;
- partners per chromosome: a chromosome whose blocks lie on two or more other
  chromosomes is a fusion/fission or translocation since the duplication;
- block divergence spread: blocks resolved at different times differ in
  divergence (asynchronous rediploidization);
- partner pairs feed windowed-homeologs, so residual tetrasomy and homeologous
  exchange are tested per block pair, not only per whole-chromosome pair.
"""

import csv
import math
import os
import random
import statistics
from collections import defaultdict

from .common import homeologs_dir, log, window_rows
from .kmer_tables import load_sequences

MAP_K = 23
MAP_Z = 3.0
MIN_BLOCK_WINDOWS = 6
MAX_GAP = 2
NULL_PERMUTATIONS = 100
NULL_QUANTILE = 0.95
MIN_TARGETS = 3
# a pair of chromosome numbers is tested for residual tetrasomy / exchange when
# their blocks together cover at least this much
PAIR_MIN_BP = 2_000_000

WINDOWS_TSV = "windowed_homeology_map.tsv"  # raw per-window rows (gitignored)
MAP_TSV = "homeology_map.tsv"
BLOCKS_TSV = "homeolog_blocks.tsv"
SUMMARY_TSV = "homeology_map_summary.tsv"

MAP_FIELDS = ["anchor", "chrom", "win_start", "win_end", "best", "best_dist", "second", "second_dist",
              "median_dist", "z", "homeolog_like"]
BLOCK_FIELDS = ["anchor", "chrom", "partner", "start", "end", "length", "n_windows", "median_dist",
                "median_z", "position", "reciprocal"]


def reference_hap(units):
    """The haplotype with the most chromosome-scale units (ties: first by name)."""
    counts = defaultdict(int)
    for u in units:
        counts[u["hap"]] += 1
    return min(counts, key=lambda h: (-counts[h], h)) if counts else None


def compute_homeology_map(seq_tsv, outdir, samtools_bin, fastk_bin, logex_bin, histex_bin,
                          window, step, threads, min_segment_bp, k=MAP_K):
    """All-vs-all windowed comparison within one haplotype, then blocks.
    Returns the summary dict, or None when there are too few chromosomes."""
    from .windowed import compute_windowed_groups

    units = load_sequences(seq_tsv)
    ref = reference_hap(units)
    group = [u for u in units if u["hap"] == ref]
    if len(group) <= MIN_TARGETS:
        log(f"homeology map: only {len(group)} chromosomes in {ref} -- skipping")
        return None
    log(f"homeology map: {len(group)} chromosomes of {ref}, each window against every other "
        f"chromosome (k={k})")
    compute_windowed_groups({"map": group}, homeologs_dir(outdir), samtools_bin, fastk_bin,
                            logex_bin, histex_bin, k, window, step, threads, overview_title="",
                            all_tsv_name=WINDOWS_TSV, overview_png_name="", only_cross_chrom=True,
                            species_outdir=outdir, plots=False)
    return analyse_map(outdir, min_segment_bp)


def load_map_windows(path):
    """{anchor: {win_start: (win_end, chrom_a, {chrom_b: distance})}}"""
    out = defaultdict(dict)
    for r, d in window_rows(path):
        entry = out[r["unit_a"]].setdefault(int(r["win_start"]),
                                            (int(r["win_end"]), int(r["chrom_a"]), {}))
        entry[2][int(r["chrom_b"])] = d
    return out


def score_window(dists):
    """(best, best_dist, second, second_dist, median, z) for one window's
    {chrom: distance}; None with too few targets."""
    if len(dists) < MIN_TARGETS:
        return None
    ranked = sorted(dists.items(), key=lambda x: x[1])
    values = [d for _, d in ranked]
    med = statistics.median(values)
    mad = statistics.median(abs(v - med) for v in values)
    spread = max(1.4826 * mad, 1e-4)
    (best, bd), (second, sd) = ranked[0], ranked[1]
    return best, bd, second, sd, med, (med - bd) / spread


def label_runs(labels, max_gap=MAX_GAP):
    """Runs of one repeated non-None label, bridging up to max_gap windows of
    anything else. Returns [(label, first_index, last_index, n_labelled)]."""
    out = []
    i, n = 0, len(labels)
    while i < n:
        lab = labels[i]
        if lab is None:
            i += 1
            continue
        first = last = i
        count, gap, j = 1, 0, i + 1
        while j < n and gap <= max_gap:
            if labels[j] == lab:
                last, count, gap = j, count + 1, 0
            else:
                gap += 1
            j += 1
        out.append((lab, first, last, count))
        i = last + 1
    return out


def null_block_run(tracks, n_perm=NULL_PERMUTATIONS, quantile=NULL_QUANTILE, seed=0):
    """95th percentile of the genome-wide longest same-partner run with the
    window labels shuffled over the whole genome. Shuffling within a chromosome
    would keep its composition: one homeologous end to end would keep a run as
    long as itself."""
    rng = random.Random(seed)
    pooled = [lab for labels in tracks for lab in labels]
    lengths = [len(labels) for labels in tracks]
    maxima = []
    for _ in range(n_perm):
        rng.shuffle(pooled)
        best, offset = 0, 0
        for n in lengths:
            best = max([best] + [r[3] for r in label_runs(pooled[offset:offset + n])])
            offset += n
        maxima.append(best)
    maxima.sort()
    return maxima[int(quantile * (len(maxima) - 1))] if maxima else 0


def position(start, end, span_start, span_end, frac=0.1):
    lo = span_start + frac * (span_end - span_start)
    hi = span_end - frac * (span_end - span_start)
    if start <= lo and end >= hi:
        return "whole"
    if start <= lo:
        return "start"
    if end >= hi:
        return "end"
    return "interior"


def analyse_map(outdir, min_segment_bp, z_min=MAP_Z):
    hdir = homeologs_dir(outdir)
    path = os.path.join(hdir, WINDOWS_TSV)
    if not os.path.exists(path):
        return None
    windows = load_map_windows(path)
    map_rows, tracks, per_anchor = [], [], {}
    for anchor in sorted(windows):
        wins = windows[anchor]
        scored = []
        for s in sorted(wins):
            e, chrom, dists = wins[s]
            sc = score_window(dists)
            if sc is None:
                continue
            best, bd, second, sd, med, z = sc
            like = z >= z_min
            scored.append((s, e, chrom, best if like else None, bd, med, z))
            map_rows.append(dict(anchor=anchor, chrom=f"chr{chrom:02d}", win_start=s, win_end=e,
                                 best=f"chr{best:02d}", best_dist=f"{bd:.5f}", second=f"chr{second:02d}",
                                 second_dist=f"{sd:.5f}", median_dist=f"{med:.5f}", z=f"{z:.2f}",
                                 homeolog_like="yes" if like else "no"))
        if scored:
            per_anchor[anchor] = scored
            tracks.append([x[3] for x in scored])
    null = null_block_run(tracks)
    min_run = max(MIN_BLOCK_WINDOWS, null + 1)

    blocks = []
    for anchor, scored in per_anchor.items():
        chrom = scored[0][2]
        span_start, span_end = scored[0][0], scored[-1][1]
        for lab, i, j, n in label_runs([x[3] for x in scored]):
            s, e = scored[i][0], scored[j][1]
            if n < min_run or e - s + 1 < min_segment_bp:
                continue
            inside = [x for x in scored[i:j + 1] if x[3] == lab]
            blocks.append(dict(anchor=anchor, chrom=f"chr{chrom:02d}", partner=f"chr{lab:02d}",
                               start=s, end=e, length=e - s + 1, n_windows=n,
                               median_dist=f"{statistics.median(x[4] for x in inside):.5f}",
                               median_z=f"{statistics.median(x[6] for x in inside):.2f}",
                               position=position(s, e, span_start, span_end)))
    partners_of = defaultdict(set)
    for b in blocks:
        partners_of[b["chrom"]].add(b["partner"])
    for b in blocks:
        b["reciprocal"] = "yes" if b["chrom"] in partners_of.get(b["partner"], ()) else "no"

    write_tsv(os.path.join(hdir, MAP_TSV), map_rows, MAP_FIELDS)
    write_tsv(os.path.join(hdir, BLOCKS_TSV), blocks, BLOCK_FIELDS)
    summary = summarize_map(per_anchor, blocks, null, min_run)
    write_tsv(os.path.join(hdir, SUMMARY_TSV),
              [dict(metric=k, value=v) for k, v in summary.items()], ["metric", "value"])
    plot_map(hdir, per_anchor, blocks)
    log(f"homeology map: {summary['n_blocks']} reciprocal homeolog block(s) "
        f"({summary['n_blocks_unmatched']} unmatched) covering {summary['duplicated_frac']} of "
        f"the genome; {summary['n_partner_pairs']} chromosome pair(s); "
        f"{summary['multi_partner_chromosomes']} chromosome(s) with blocks on 2+ others")
    return summary


def summarize_map(per_anchor, blocks, null, min_run):
    total = sum(x[1] - x[0] + 1 for scored in per_anchor.values() for x in scored)
    unmatched = [b for b in blocks if b["reciprocal"] != "yes"]
    blocks = [b for b in blocks if b["reciprocal"] == "yes"]
    covered = sum(b["length"] for b in blocks)
    pair_bp = defaultdict(int)
    for b in blocks:
        pair_bp[tuple(sorted((b["chrom"], b["partner"])))] += b["length"]
    pairs = {p for p, bp in pair_bp.items() if bp >= PAIR_MIN_BP}
    partners = defaultdict(set)
    for b in blocks:
        if tuple(sorted((b["chrom"], b["partner"]))) in pairs:
            partners[b["chrom"]].add(b["partner"])
    multi = sorted(c for c, p in partners.items() if len(p) >= 2)
    dists = [float(b["median_dist"]) for b in blocks]
    # length-weighted, so many short blocks do not dominate
    cv = ""
    if len(dists) >= 3:
        w = [b["length"] for b in blocks]
        mean = sum(d * x for d, x in zip(dists, w)) / sum(w)
        var = sum(x * (d - mean) ** 2 for d, x in zip(dists, w)) / sum(w)
        cv = f"{math.sqrt(var) / mean:.3f}" if mean > 0 else ""
    return dict(
        n_chromosomes=len(per_anchor),
        n_blocks=len(blocks),
        n_blocks_unmatched=len(unmatched),
        duplicated_frac=f"{covered / total:.3f}" if total else "0.000",
        duplicated_bp=covered,
        n_partner_pairs=len(pairs),
        partner_pairs=";".join(f"{a}x{b}" for a, b in sorted(pairs)),
        multi_partner_chromosomes=len(multi),
        multi_partner_list=";".join(f"{c}:{'+'.join(sorted(partners[c]))}" for c in multi),
        block_dist_median=f"{statistics.median(dists):.4f}" if dists else "",
        block_dist_cv=cv,
        null_run_windows=null,
        min_block_windows=min_run,
    )


def partner_pairs(outdir):
    """[(chrom_a, chrom_b)] ints from the map summary, for windowed-homeologs."""
    path = os.path.join(homeologs_dir(outdir), SUMMARY_TSV)
    if not os.path.exists(path):
        return []
    with open(path) as f:
        s = {r["metric"]: r["value"] for r in csv.DictReader(f, delimiter="\t")}
    out = []
    for p in filter(None, s.get("partner_pairs", "").split(";")):
        a, b = p.split("x")
        out.append((int(a[3:]), int(b[3:])))
    return out


def read_summary(outdir):
    path = os.path.join(homeologs_dir(outdir), SUMMARY_TSV)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {r["metric"]: r["value"] for r in csv.DictReader(f, delimiter="\t")}


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def plot_map(hdir, per_anchor, blocks):
    """Each chromosome a bar, painted by the chromosome its blocks lie on."""
    if not per_anchor:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    anchors = sorted(per_anchor, key=lambda a: per_anchor[a][0][2])
    chroms = sorted({x[2] for a in anchors for x in per_anchor[a][:1]})
    cmap = plt.get_cmap("tab20")
    color = {f"chr{c:02d}": cmap(i % 20) for i, c in enumerate(chroms)}
    fig, ax = plt.subplots(figsize=(10, max(2.5, 0.28 * len(anchors) + 1)))
    for y, a in enumerate(anchors):
        end = per_anchor[a][-1][1]
        ax.barh(y, end / 1e6, color="0.88", height=0.7)
        for b in (b for b in blocks if b["anchor"] == a):
            matched = b.get("reciprocal") == "yes"
            ax.barh(y, b["length"] / 1e6, left=b["start"] / 1e6, color=color.get(b["partner"], "k"),
                    height=0.7, alpha=1.0 if matched else 0.45, hatch=None if matched else "////",
                    edgecolor="white" if not matched else None, linewidth=0)
            if b["length"] / 1e6 > 0.08 * end / 1e6:
                ax.text((b["start"] + b["length"] / 2) / 1e6, y, b["partner"][3:], ha="center",
                        va="center", fontsize=6)
    ax.set_yticks(range(len(anchors)))
    ax.set_yticklabels([f"chr{per_anchor[a][0][2]:02d}" for a in anchors], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("position (Mb)")
    ax.set_title("Homeolog blocks: each chromosome painted by the chromosome its duplicate lies on\n"
                 "(hatched: no block back from the partner, not counted; grey: none detected)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(hdir, "homeology_map.png"), dpi=150)
    plt.close(fig)
