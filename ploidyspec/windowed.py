import csv
import math
import os
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from .common import (
    homeologs_dir,
    log,
    run,
    windowed_dir,
)
from .kmer_tables import build_one, chrom_fasta_path, ktab_prefix_path, load_sequences

FIELDNAMES = [
    "group",
    "win_start",
    "win_end",
    "unit_a",
    "unit_b",
    "hap_a",
    "hap_b",
    "chrom_a",
    "chrom_b",
    "kmers_a",
    "shared",
    "containment",
    "distance",
]

# Measured FastK -p memory: ~13 bytes per bp of the table's sequence (723 MB
# for a 56 Mb chromosome), so a job needs threads x this x longest chromosome.
PROFILE_BYTES_PER_BP = 14

# Profex -z prints "  <start> -   <end> (<count>)" runs, zeros omitted. This awk
# program sums run lengths into fixed bins of `g` k-mer positions, so the
# per-position profile (one line per run, ~10M lines for a 50 Mb chromosome)
# never reaches Python.
_BIN_AWK = r"""
$2 == "-" && $3 ~ /^[0-9]+$/ {
    s = $1; e = $3
    while (s <= e) {
        b = int(s / g); be = (b + 1) * g - 1; t = (e < be) ? e : be
        n[b] += t - s + 1; s = t + 1
    }
}
END { for (b in n) print b "\t" n[b] }
"""


def make_windows(minlen, window, step):
    windows = []
    start = 1
    min_tail = max(1000, window // 10)
    while start <= minlen:
        end = min(start + window - 1, minlen)
        if end - start + 1 >= min_tail:
            windows.append((start, end))
        start += step
    return windows


def build_window_ktab(
    samtools_bin, fastk_bin, k, source, seq_id, start, end, tmp_dir, unit_id
):
    region = f"{seq_id}:{start}-{end}"
    fa_path = os.path.join(tmp_dir, f"{unit_id}.fa")
    with open(fa_path, "w") as out:
        run([samtools_bin, "faidx", source, region], stdout=out)
    ktab_prefix = os.path.join(tmp_dir, unit_id)
    # -P scopes FastK's block-sort scratch files to this window's own tmp_dir: the same unit_id
    # recurs across many concurrently-processed windows, and FastK names scratch files from the
    # basename only, so sharing $TMPDIR across windows causes cross-window filename collisions.
    run(
        [fastk_bin, f"-k{k}", "-t1", "-T1", f"-N{ktab_prefix}", f"-P{tmp_dir}", fa_path]
    )
    return ktab_prefix


def find_profex(fastk_bin):
    sibling = os.path.join(os.path.dirname(fastk_bin), "Profex")
    if os.path.exists(sibling):
        return sibling
    found = shutil.which("Profex")
    if not found:
        raise SystemExit("Profex (part of FastK) not found next to FastK or on PATH")
    return found


def profile_bins(fastk_bin, profex_bin, k, fasta, table_prefix, tmp_dir, name, bin_size):
    """
    {bin: number of k-mer positions of `fasta` whose k-mer is present in the
    table}, bins of `bin_size` positions along the sequence. Uses FastK's
    profile mode against another sequence's table (-p:table), so it needs no
    alignment and no shared coordinates.
    """
    root = os.path.join(tmp_dir, name)
    run([fastk_bin, f"-k{k}", "-T1", f"-p:{table_prefix}", f"-N{root}", f"-P{tmp_dir}", fasta])
    profex = subprocess.Popen([profex_bin, "-z", root, "1"], stdout=subprocess.PIPE)
    out = subprocess.run(["awk", "-v", f"g={bin_size}", _BIN_AWK], stdin=profex.stdout,
                         capture_output=True, text=True, check=True).stdout
    profex.stdout.close()
    if profex.wait() != 0:
        raise RuntimeError(f"Profex failed on {root}")
    for f in os.listdir(tmp_dir):
        if f.startswith((name + ".", "." + name + ".")):
            os.remove(os.path.join(tmp_dir, f))
    bins = {}
    for line in out.splitlines():
        b, n = line.split("\t")
        bins[int(b)] = int(n)
    return bins


def window_counts(bins, windows, bin_size):
    """Sum per-bin counts into windows given as 1-based (start, end) bp."""
    out = []
    for s, e in windows:
        lo, hi = (s - 1) // bin_size, (e - 1) // bin_size
        out.append(sum(bins.get(b, 0) for b in range(lo, hi + 1)))
    return out


def containment_distance(shared, total, k):
    """Mash-style containment distance -ln(c)/k; 1.0 when nothing is shared."""
    if total <= 0:
        return None, None
    c = shared / total
    return c, (-math.log(c) / k if c > 0 else 1.0)


def compute_windowed_groups(
    labeled_groups,
    outdir,
    samtools_bin,
    fastk_bin,
    logex_bin,
    histex_bin,
    k,
    window,
    step,
    threads,
    overview_title,
    all_tsv_name,
    overview_png_name,
    only_cross_chrom=False,
    species_outdir=None,
):
    """
    Core windowed-comparison loop. For each label -> group of units, every unit
    is tiled along its *own* coordinates, and each window's k-mers are looked up
    in every other unit's whole-sequence k-mer table: containment c = share of the
    window's k-mer positions present anywhere in the other copy, distance =
    -ln(c)/k. Rows are directional (positions are along unit_a).

    This replaced comparing windows at equal coordinates across copies, which
    assumed the assemblies were collinear and fell out of register after the
    first indel or gap larger than a window: a pair 0.004 apart over the whole
    chromosome (ddEmpNigr1 chr03) read ~0.99 Jaccard distance in most windows.

    Valid k-mer positions per window come from profiling each unit against its
    own table, so assembly gaps (Ns) reduce the denominator rather than reading
    as divergence. Memory per concurrent FastK profile scales with the table
    size (~14 bytes per distinct k-mer of the other copy).
    """
    species_outdir = species_outdir or os.path.dirname(outdir.rstrip("/"))
    profex_bin = find_profex(fastk_bin)
    bin_size = math.gcd(window, step)
    rows_by_label = {}
    for label in sorted(labeled_groups):
        group = sorted(labeled_groups[label], key=lambda u: (u["hap"], int(u["chrom"])))
        if len(group) < 2:
            log(f"{label}: only {len(group)} unit present, skipping windowed comparison")
            continue
        for u in group:  # whole-unit FASTA and k-mer table at the window k
            build_one(samtools_bin, fastk_bin, u, k, species_outdir)

        pairs = [(a, b) for a in group for b in group
                 if a is not b and not (only_cross_chrom and a["chrom"] == b["chrom"])]
        peak_gb = min(threads, len(group) ** 2) * PROFILE_BYTES_PER_BP * max(
            int(u["length"]) for u in group) / 1e9
        log(f"{label}: {len(group)} units, {len(pairs)} directional comparisons "
            f"({window}bp windows, step {step}bp; ~{peak_gb:.1f} GB peak with {threads} threads)")
        tmp_root = os.path.join(outdir, "tmp_windowed", label)
        os.makedirs(tmp_root, exist_ok=True)

        def profile(a, b):
            tmp = os.path.join(tmp_root, f"{a['unit_id']}__{b['unit_id']}")
            os.makedirs(tmp, exist_ok=True)
            bins = profile_bins(fastk_bin, profex_bin, k, chrom_fasta_path(species_outdir, a["unit_id"]),
                                ktab_prefix_path(species_outdir, b["unit_id"], k), tmp, "p", bin_size)
            shutil.rmtree(tmp, ignore_errors=True)
            return a["unit_id"], b["unit_id"], bins

        jobs = [(u, u) for u in group] + pairs  # self-profile = valid positions
        bins = {}
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [ex.submit(profile, a, b) for a, b in jobs]
            for i, fut in enumerate(as_completed(futs), 1):
                ua, ub, bb = fut.result()
                bins[(ua, ub)] = bb
                if i % max(1, len(jobs) // 5) == 0 or i == len(jobs):
                    log(f"  {label}: [{i}/{len(jobs)}] profiles done")
        shutil.rmtree(tmp_root, ignore_errors=True)

        label_rows = []
        for a, b in pairs:
            windows = make_windows(int(a["length"]) - k + 1, window, step)
            valid = window_counts(bins[(a["unit_id"], a["unit_id"])], windows, bin_size)
            shared = window_counts(bins[(a["unit_id"], b["unit_id"])], windows, bin_size)
            for (s, e), n, m in zip(windows, valid, shared):
                c, d = containment_distance(m, n, k)
                if c is None:
                    continue
                label_rows.append(dict(
                    group=label, win_start=s, win_end=e,
                    unit_a=a["unit_id"], unit_b=b["unit_id"], hap_a=a["hap"], hap_b=b["hap"],
                    chrom_a=int(a["chrom"]), chrom_b=int(b["chrom"]),
                    kmers_a=n, shared=m, containment=c, distance=d,
                ))

        write_group_tsv(outdir, label, label_rows)
        plot_group(outdir, label, label_rows)
        plot_group_heatmap(outdir, label, label_rows)
        rows_by_label[label] = label_rows

    write_all_windows_tsv(outdir, rows_by_label, all_tsv_name)
    if rows_by_label:
        plot_overview(outdir, rows_by_label, overview_title, overview_png_name)
        heatmap_png_name = overview_png_name.replace(".png", "_heatmap.png")
        plot_overview_heatmap(outdir, rows_by_label, overview_title, heatmap_png_name)
    return rows_by_label


def compute_windowed(
    seq_tsv,
    outdir,
    samtools_bin,
    fastk_bin,
    logex_bin,
    histex_bin,
    k,
    window,
    step,
    threads,
):
    units = load_sequences(seq_tsv)
    groups = defaultdict(list)
    for u in units:
        groups[int(u["chrom"])].append(u)
    labeled_groups = {f"chr{c:02d}": g for c, g in groups.items()}
    return compute_windowed_groups(
        labeled_groups,
        windowed_dir(outdir),
        samtools_bin,
        fastk_bin,
        logex_bin,
        histex_bin,
        k,
        window,
        step,
        threads,
        overview_title="Genome-wide windowed haplotype k-mer divergence (relative chromosome position)",
        all_tsv_name="windowed_all.tsv",
        overview_png_name="windowed_genome_overview.png",
        species_outdir=outdir,
    )


def compute_windowed_homeologs(
    seq_tsv,
    outdir,
    homeolog_pairs,
    samtools_bin,
    fastk_bin,
    logex_bin,
    histex_bin,
    k,
    window,
    step,
    threads,
):
    """homeolog_pairs: list of (chrom_a:int, chrom_b:int) candidate ancestral chromosome pairs."""
    units = load_sequences(seq_tsv)
    groups = defaultdict(list)
    for u in units:
        groups[int(u["chrom"])].append(u)

    labeled_groups = {}
    for a, b in homeolog_pairs:
        label = f"chr{a:02d}x{b:02d}"
        labeled_groups[label] = groups.get(a, []) + groups.get(b, [])

    return compute_windowed_groups(
        labeled_groups,
        homeologs_dir(outdir),
        samtools_bin,
        fastk_bin,
        logex_bin,
        histex_bin,
        k,
        window,
        step,
        threads,
        overview_title="Windowed k-mer divergence between candidate ancestral (paleopolyploid) chromosome pairs",
        all_tsv_name="windowed_homeologs_all.tsv",
        overview_png_name="windowed_homeologs_overview.png",
        only_cross_chrom=True,
        species_outdir=outdir,
    )


def write_group_tsv(outdir, label, rows):
    path = os.path.join(outdir, f"windowed_{label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["win_start"], r["hap_a"], r["hap_b"])):
            w.writerow(
                {
                    k: (f"{v:.6f}" if k in ("containment", "distance") else v)
                    for k, v in r.items()
                }
            )


def write_all_windows_tsv(outdir, rows_by_label, filename):
    path = os.path.join(outdir, filename)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader()
        for label in sorted(rows_by_label):
            for r in sorted(
                rows_by_label[label],
                key=lambda r: (r["win_start"], r["hap_a"], r["hap_b"]),
            ):
                w.writerow(
                    {
                        k: (f"{v:.6f}" if k in ("containment", "distance") else v)
                        for k, v in r.items()
                    }
                )


def _pair_key_and_label(r):
    """
    Grouping key for a plotted series. Same-chromosome-number rows (the standard
    per-chromosome windowed stage) key on haplotype label alone, e.g. "HAP1 vs HAP2"
    -- consistent and comparable across every chromosome subplot. Cross-chromosome
    rows (the homeolog-pair stage) key on the specific unit pair instead: with only
    2 chromosome numbers per group but potentially several haplotype copies, two
    different unit pairs (e.g. HAP1_chr01-HAP2_chr05 and HAP2_chr01-HAP1_chr05) can
    share the same (hap_a, hap_b) label while being genuinely different comparisons,
    so collapsing them onto one line would silently interleave unrelated series.
    """
    if r["chrom_a"] == r["chrom_b"]:
        return (r["hap_a"], r["hap_b"]), f"{r['hap_a']} vs {r['hap_b']}"
    return (r["unit_a"], r["unit_b"]), f"{r['unit_a']} vs {r['unit_b']}"


def plot_group(outdir, label, rows):
    if not rows:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keyed = {}
    for r in rows:
        key, lbl = _pair_key_and_label(r)
        keyed.setdefault(key, (lbl, []))[1].append(r)

    fig, ax = plt.subplots(figsize=(10, 4))
    for key in sorted(keyed):
        lbl, sub = keyed[key]
        sub = sorted(sub, key=lambda r: r["win_start"])
        xs = [(r["win_start"] + r["win_end"]) / 2 / 1e6 for r in sub]
        ys = [r["distance"] for r in sub]
        ax.plot(xs, ys, marker="o", markersize=2, linewidth=1, label=lbl)
    ax.set_xlabel("position (Mb)")
    ax.set_ylabel("windowed k-mer distance (-ln c / k)")
    ax.set_title(f"{label} windowed divergence")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"windowed_{label}.png"), dpi=150)
    plt.close(fig)


def _pair_matrix(rows):
    """rows -> (row_labels, windows, 2D array of distance, NaN where a pair has no
    value at a given window). Each row is one directional pair along its unit_a's
    own coordinates; window starts are a shared grid, so copies of different
    lengths simply end at different columns."""
    import numpy as np

    keyed = {}
    for r in rows:
        key, lbl = _pair_key_and_label(r)
        keyed.setdefault(key, (lbl, []))[1].append(r)
    pair_keys = sorted(keyed)
    windows = sorted(set((r["win_start"], r["win_end"]) for r in rows))
    win_index = {w: i for i, w in enumerate(windows)}

    mat = np.full((len(pair_keys), len(windows)), np.nan)
    row_labels = []
    for i, key in enumerate(pair_keys):
        lbl, sub = keyed[key]
        row_labels.append(lbl)
        for r in sub:
            j = win_index.get((r["win_start"], r["win_end"]))
            if j is not None:
                mat[i, j] = r["distance"]
    return row_labels, windows, mat


def plot_group_heatmap(outdir, label, rows):
    """One row per haplotype pair, one column per window, color = jaccard distance
    -- makes mosaic/patchy divergence along a chromosome visible at a glance, which
    the line-plot equivalent (plot_group) does not do well once there are more than
    a couple of pairs."""
    if not rows:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    row_labels, windows, mat = _pair_matrix(rows)
    finite = mat[np.isfinite(mat)]
    vmin, vmax = np.percentile(finite, [2, 98]) if finite.size else (0.0, 1.0)
    if vmin == vmax:
        vmax = vmin + 1e-6

    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("0.85")
    maxend = windows[-1][1] if windows else 1

    fig, ax = plt.subplots(figsize=(10, max(1.5, 0.35 * len(row_labels) + 0.5)))
    im = ax.imshow(
        mat,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
        extent=[0, maxend / 1e6, len(row_labels), 0],
    )
    ax.set_yticks([i + 0.5 for i in range(len(row_labels))])
    ax.set_yticklabels(row_labels, fontsize=7)
    ax.set_xlabel("position (Mb)")
    ax.set_title(f"{label} windowed divergence (2nd-98th pct: {vmin:.3f}-{vmax:.3f})")
    fig.colorbar(im, ax=ax, shrink=0.7, label="distance (-ln c / k)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"windowed_{label}_heatmap.png"), dpi=150)
    plt.close(fig)


def plot_overview_heatmap(outdir, rows_by_label, title, filename):
    """Same idea as plot_group_heatmap but one subplot per chromosome (or homeolog
    pair) group, all sharing one color scale so patterns are comparable across the
    whole genome, not just within one chromosome."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    labels = sorted(rows_by_label)
    all_rows = [r for rows in rows_by_label.values() for r in rows]
    if not all_rows:
        return
    all_dist = np.array([r["distance"] for r in all_rows])
    vmin, vmax = np.percentile(all_dist, [2, 98])
    if vmin == vmax:
        vmax = vmin + 1e-6

    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("0.85")

    ncols = min(4, len(labels))
    nrows_grid = math.ceil(len(labels) / ncols)
    max_pairs = max(
        len(set(_pair_key_and_label(r)[0] for r in rows_by_label[l])) for l in labels
    )
    fig, axes = plt.subplots(
        nrows_grid,
        ncols,
        figsize=(4 * ncols, max(1.8, 0.3 * max_pairs) * nrows_grid),
        squeeze=False,
    )
    axes = axes.flatten()

    im = None
    for ax, label in zip(axes, labels):
        row_labels, windows, mat = _pair_matrix(rows_by_label[label])
        maxend = windows[-1][1] if windows else 1
        im = ax.imshow(
            mat,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
            extent=[0, maxend / 1e6, len(row_labels), 0],
        )
        ax.set_yticks([i + 0.5 for i in range(len(row_labels))])
        ax.set_yticklabels(row_labels, fontsize=5)
        ax.set_title(label, fontsize=9)

    for ax in axes[len(labels) :]:
        ax.axis("off")

    fig.suptitle(f"{title}\n(2nd-98th pct: {vmin:.3f}-{vmax:.3f})")
    fig.tight_layout(rect=[0, 0, 0.93, 0.90])
    fig.colorbar(
        im, ax=axes[: len(labels)].tolist(), shrink=0.6, label="distance (-ln c / k)"
    )
    fig.savefig(os.path.join(outdir, filename), dpi=150)
    plt.close(fig)


def plot_overview(outdir, rows_by_label, title, filename):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = sorted(rows_by_label)
    all_rows = [r for rows in rows_by_label.values() for r in rows]
    all_keys = sorted(
        set(_pair_key_and_label(r) for r in all_rows), key=lambda kl: kl[1]
    )
    cmap = plt.get_cmap("tab10")
    color_of = {key: cmap(i % 10) for i, (key, _) in enumerate(all_keys)}
    label_of = dict(all_keys)
    # same-chromosome mode (standard windowed stage): a color means the same haplotype
    # pair on every subplot, so one shared legend at the bottom is meaningful. Cross-
    # chromosome mode (homeolog pairs): colors aren't comparable across subplots (each
    # subplot compares a different pair of chromosome numbers), so each gets its own.
    same_chrom_mode = all(r["chrom_a"] == r["chrom_b"] for r in all_rows)

    ncols = min(4, len(labels))
    nrows = math.ceil(len(labels) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 2.5 * nrows), sharey=True, squeeze=False
    )
    axes = axes.flatten()

    for ax, label in zip(axes, labels):
        rows = rows_by_label[label]
        maxend = max(r["win_end"] for r in rows)
        keyed = {}
        for r in rows:
            key, _ = _pair_key_and_label(r)
            keyed.setdefault(key, []).append(r)
        for key, sub in sorted(keyed.items()):
            lbl = label_of[key]
            sub = sorted(sub, key=lambda r: r["win_start"])
            xs = [100 * (r["win_start"] + r["win_end"]) / 2 / maxend for r in sub]
            ys = [r["distance"] for r in sub]
            ax.plot(
                xs,
                ys,
                linewidth=1,
                color=color_of[key],
                label=lbl if not same_chrom_mode else None,
            )
        ax.set_title(label, fontsize=9)
        ax.set_ylim(bottom=0)
        if not same_chrom_mode:
            ax.legend(fontsize=5, loc="lower left")

    for ax in axes[len(labels) :]:
        ax.axis("off")

    fig.suptitle(title)
    if same_chrom_mode:
        handles = [
            plt.Line2D([0], [0], color=color_of[key], label=lbl)
            for key, lbl in all_keys
        ]
        fig.legend(
            handles=handles, loc="lower center", ncol=min(len(all_keys), 6), fontsize=7
        )
        fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    else:
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(outdir, filename), dpi=150)
    plt.close(fig)
