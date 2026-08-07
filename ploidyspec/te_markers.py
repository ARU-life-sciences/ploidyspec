"""
Differential fossil-TE k-mer markers for subgenome resolution/phasing
(Jaron/Cerca method: github.com/KamilSJaron/k-mer-approaches-for-biodiversity-genomics).

Isolates k-mers that are both high-copy (repetitive/TE-like) and differentially
represented between two haplotype copies of the same chromosome -- a much more
targeted signal for telling apart allopolyploid subgenomes than bulk k-mer
Jaccard distance, since it directly fingerprints independent repeat-family
expansion in each parental lineage before hybridization.
"""

import csv
import os
import shutil
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from .common import log, run, subgenomes_dir
from .kmer_tables import build_one, ktab_prefix_path, load_sequences
from .windowed import build_window_ktab, make_windows

DEFAULT_MARKER_K = 13
DEFAULT_MIN_COUNT = 100
DEFAULT_MIN_RATIO = 2.0


def dump_high_copy_kmers(tabex_bin, ktab_prefix, min_count):
    """{kmer: count} for every k-mer at or above min_count, via `Tabex -A -t<min_count>`."""
    result = run([tabex_bin, "-A", f"-t{min_count}", ktab_prefix])
    counts = {}
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            counts[parts[0]] = int(parts[1])
    return counts


def differential_markers(counts_a, counts_b, min_ratio):
    """
    Classify every high-copy k-mer seen on either side as an 'a' or 'b' marker
    when its count on one side is at least min_ratio times the other (a side
    absent from one dump is treated as count 0; the ratio's denominator uses
    max(other_count, 1) so a lone hit against zero still requires the winning
    side to actually clear min_ratio, not just be nonzero).

    Returns (markers, n_markers_a, n_markers_b):
      markers: {kmer: (count_a, count_b, assigned)}, assigned in {"a", "b"}.
    """
    markers = {}
    n_a = 0
    n_b = 0
    for kmer in set(counts_a) | set(counts_b):
        ca = counts_a.get(kmer, 0)
        cb = counts_b.get(kmer, 0)
        if ca >= min_ratio * max(cb, 1):
            markers[kmer] = (ca, cb, "a")
            n_a += 1
        elif cb >= min_ratio * max(ca, 1):
            markers[kmer] = (ca, cb, "b")
            n_b += 1
    return markers, n_a, n_b


def classify_window(a_hits, b_hits, min_ratio):
    """
    Classify a window by raw marker-hit counts, NOT hits normalized by each
    side's total marker-pool size. Both a_hits and b_hits were counted
    against the SAME window (a fair, equal-denominator comparison as-is),
    but the two marker pools are usually very different sizes globally (e.g.
    3612 vs 1926 markers) -- dividing each side by its own pool size before
    comparing introduces a denominator asymmetry that systematically favors
    whichever pool happens to be smaller, independent of true enrichment.
    (Caught empirically: windows where a_hits clearly exceeded b_hits were
    still being classified "b" once normalized by each side's pool size.)
    """
    if a_hits == 0 and b_hits == 0:
        return "none"
    if a_hits >= min_ratio * max(b_hits, 1):
        return "a"
    if b_hits >= min_ratio * max(a_hits, 1):
        return "b"
    return "ambiguous"


def compute_te_markers(
    seq_tsv,
    outdir,
    samtools_bin,
    fastk_bin,
    tabex_bin,
    marker_k,
    min_count,
    min_ratio,
    threads,
):
    """
    For every chromosome number with >=2 haplotype-copy units: build (if
    missing) each unit's marker_k .ktab, dump its high-copy k-mers once per
    unit (not once per pair -- a unit is shared across every pair it's in),
    then classify differential markers for every pair within the group.
    """
    units = load_sequences(seq_tsv)
    groups = defaultdict(list)
    for u in units:
        groups[int(u["chrom"])].append(u)

    n_multi_hap = sum(1 for g in groups.values() if len(g) >= 2)
    if n_multi_hap == 0:
        log(
            f"no chromosome number has >=2 haplotype-copy units (found "
            f"{len(groups)} chromosome numbers, all single-copy) -- nothing to "
            f"compare, te-markers needs at least 2 haplotype-scale assemblies "
            f"of the same individual; writing an empty te_markers_summary.tsv"
        )

    results = []
    for chrom in sorted(groups):
        group = sorted(groups[chrom], key=lambda u: u["hap"])
        if len(group) < 2:
            continue

        log(
            f"chr{chrom:02d}: building marker-k={marker_k} tables for {len(group)} units"
        )
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [
                ex.submit(build_one, samtools_bin, fastk_bin, u, marker_k, outdir)
                for u in group
            ]
            for fut in as_completed(futs):
                fut.result()

        log(f"chr{chrom:02d}: dumping k-mers with count >= {min_count} per unit")
        dumps = {}
        for u in group:
            prefix = ktab_prefix_path(outdir, u["unit_id"], marker_k)
            dumps[u["unit_id"]] = dump_high_copy_kmers(tabex_bin, prefix, min_count)
            log(f"  {u['unit_id']}: {len(dumps[u['unit_id']])} k-mers >= {min_count}x")

        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                ua, ub = group[i]["unit_id"], group[j]["unit_id"]
                markers, n_a, n_b = differential_markers(
                    dumps[ua], dumps[ub], min_ratio
                )
                write_markers_tsv(outdir, ua, ub, markers)
                log(f"  {ua} x {ub}: {n_a} markers favoring {ua}, {n_b} favoring {ub}")
                results.append(
                    dict(
                        unit_a=ua,
                        unit_b=ub,
                        chrom=chrom,
                        n_markers_a=n_a,
                        n_markers_b=n_b,
                        n_highcopy_a=len(dumps[ua]),
                        n_highcopy_b=len(dumps[ub]),
                    )
                )

    write_summary_tsv(outdir, results)
    return results


def write_markers_tsv(outdir, unit_a, unit_b, markers):
    path = os.path.join(subgenomes_dir(outdir), f"te_markers_{unit_a}x{unit_b}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["kmer", "count_a", "count_b", "assigned"])
        for kmer, (ca, cb, assigned) in sorted(
            markers.items(), key=lambda kv: -max(kv[1][0], kv[1][1])
        ):
            w.writerow([kmer, ca, cb, assigned])


def write_summary_tsv(outdir, results):
    path = os.path.join(subgenomes_dir(outdir), "te_markers_summary.tsv")
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
            ]
        )
        for r in results:
            w.writerow(
                [
                    r["unit_a"],
                    r["unit_b"],
                    f"chr{r['chrom']:02d}",
                    r["n_markers_a"],
                    r["n_markers_b"],
                    r["n_highcopy_a"],
                    r["n_highcopy_b"],
                ]
            )


def load_markers_tsv(outdir, unit_a, unit_b):
    """{kmer} sets for each side, from a te_markers_<a>x<b>.tsv written earlier."""
    path = os.path.join(subgenomes_dir(outdir), f"te_markers_{unit_a}x{unit_b}.tsv")
    a_markers, b_markers = set(), set()
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            (a_markers if row["assigned"] == "a" else b_markers).add(row["kmer"])
    return a_markers, b_markers


def paint_unit_windows(
    samtools_bin,
    fastk_bin,
    tabex_bin,
    unit,
    marker_k,
    window,
    step,
    a_markers,
    b_markers,
    min_ratio,
    outdir,
    threads,
):
    """
    Tile one unit's own chromosome into windows and, for each, check how much
    of its k-mer content matches the A-marker vs B-marker sets (both derived
    from a whole-chromosome comparison against its haplotype partner). A
    haplotype copy's windows should mostly match "its own" side; a window
    that instead matches the *other* side is a candidate homeologous-exchange
    or introgression breakpoint, not just noise -- that's the phasing signal.
    """
    minlen = int(unit["length"])
    windows = make_windows(minlen, window, step)
    tmp_root = os.path.join(outdir, "tmp_te_windowed", unit["unit_id"])
    os.makedirs(tmp_root, exist_ok=True)

    def do_window(wi, start, end):
        tmp_dir = os.path.join(tmp_root, f"w{wi}")
        os.makedirs(tmp_dir, exist_ok=True)
        try:
            prefix = build_window_ktab(
                samtools_bin,
                fastk_bin,
                marker_k,
                unit["source"],
                unit["seq_id"],
                start,
                end,
                tmp_dir,
                unit["unit_id"],
            )
            kmers = set(dump_high_copy_kmers(tabex_bin, prefix, min_count=1))
            a_hits = len(kmers & a_markers)
            b_hits = len(kmers & b_markers)
            a_frac = a_hits / len(a_markers) if a_markers else 0.0
            b_frac = b_hits / len(b_markers) if b_markers else 0.0
            assigned = classify_window(a_hits, b_hits, min_ratio)
            return dict(
                win_start=start,
                win_end=end,
                a_hits=a_hits,
                b_hits=b_hits,
                a_frac=a_frac,
                b_frac=b_frac,
                assigned=assigned,
            )
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    rows = []
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = [
            ex.submit(do_window, wi, s, e) for wi, (s, e) in enumerate(windows)
        ]
        for fut in as_completed(futs):
            rows.append(fut.result())
    shutil.rmtree(tmp_root, ignore_errors=True)
    return sorted(rows, key=lambda r: r["win_start"])


def compute_te_markers_windowed(
    seq_tsv,
    outdir,
    samtools_bin,
    fastk_bin,
    tabex_bin,
    marker_k,
    min_ratio,
    window,
    step,
    threads,
):
    units = {u["unit_id"]: u for u in load_sequences(seq_tsv)}
    summary_path = os.path.join(subgenomes_dir(outdir), "te_markers_summary.tsv")
    if not os.path.exists(summary_path):
        raise SystemExit(
            f"{summary_path} not found -- run the `te-markers` stage first"
        )

    with open(summary_path) as f:
        pairs = list(csv.DictReader(f, delimiter="\t"))

    for row in pairs:
        ua, ub = row["unit_a"], row["unit_b"]
        a_markers, b_markers = load_markers_tsv(outdir, ua, ub)
        if not a_markers and not b_markers:
            log(f"{ua} x {ub}: no markers -- skipping windowed painting")
            continue

        log(
            f"{ua} x {ub}: painting windows against {len(a_markers)} A-markers, "
            f"{len(b_markers)} B-markers"
        )
        pair_rows = []
        for unit_id in (ua, ub):
            rows = paint_unit_windows(
                samtools_bin,
                fastk_bin,
                tabex_bin,
                units[unit_id],
                marker_k,
                window,
                step,
                a_markers,
                b_markers,
                min_ratio,
                outdir,
                threads,
            )
            for r in rows:
                r["unit"] = unit_id
            pair_rows.extend(rows)
            n_a = sum(1 for r in rows if r["assigned"] == "a")
            n_b = sum(1 for r in rows if r["assigned"] == "b")
            log(
                f"  {unit_id}: {n_a}/{len(rows)} windows assigned 'a', "
                f"{n_b}/{len(rows)} assigned 'b'"
            )

        write_windowed_tsv(outdir, ua, ub, pair_rows)
        plot_painted_track(outdir, ua, ub, pair_rows)


def write_windowed_tsv(outdir, unit_a, unit_b, rows):
    path = os.path.join(
        subgenomes_dir(outdir), f"te_markers_windowed_{unit_a}x{unit_b}.tsv"
    )
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["unit", "win_start", "win_end", "a_hits", "b_hits", "a_frac", "b_frac", "assigned"]
        )
        for r in sorted(rows, key=lambda r: (r["unit"], r["win_start"])):
            w.writerow(
                [
                    r["unit"],
                    r["win_start"],
                    r["win_end"],
                    r["a_hits"],
                    r["b_hits"],
                    f"{r['a_frac']:.4f}",
                    f"{r['b_frac']:.4f}",
                    r["assigned"],
                ]
            )


def plot_painted_track(outdir, unit_a, unit_b, rows):
    """Painted-ideogram-style plot: one horizontal color-coded strip per unit,
    x = chromosome position, color = which subgenome's marker set each window
    matches (grey = ambiguous/none) -- a visual break in an otherwise-uniform
    strip is a candidate homeologous-exchange/introgression region."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    color_of = {"a": 0, "b": 1, "ambiguous": 2, "none": 3}
    cmap = matplotlib.colors.ListedColormap(
        ["steelblue", "crimson", "0.75", "0.9"]
    )

    units = sorted(set(r["unit"] for r in rows))
    fig, axes = plt.subplots(
        len(units), 1, figsize=(10, 1.2 * len(units)), squeeze=False
    )
    maxend = max(r["win_end"] for r in rows)
    for ax, unit_id in zip(axes[:, 0], units):
        sub = sorted(
            (r for r in rows if r["unit"] == unit_id), key=lambda r: r["win_start"]
        )
        codes = np.array([[color_of[r["assigned"]] for r in sub]])
        ax.imshow(
            codes,
            cmap=cmap,
            vmin=0,
            vmax=3,
            aspect="auto",
            extent=[0, maxend / 1e6, 0, 1],
        )
        ax.set_yticks([])
        ax.set_ylabel(unit_id, rotation=0, ha="right", va="center", fontsize=8)
    axes[-1, 0].set_xlabel("position (Mb)")
    handles = [
        matplotlib.patches.Patch(color="steelblue", label=f"matches {unit_a}"),
        matplotlib.patches.Patch(color="crimson", label=f"matches {unit_b}"),
        matplotlib.patches.Patch(color="0.75", label="ambiguous"),
        matplotlib.patches.Patch(color="0.9", label="no markers"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7)
    fig.suptitle(f"{unit_a} x {unit_b}: windowed fossil-TE marker assignment")
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(
        os.path.join(
            subgenomes_dir(outdir), f"te_markers_windowed_{unit_a}x{unit_b}.png"
        ),
        dpi=150,
    )
    plt.close(fig)
