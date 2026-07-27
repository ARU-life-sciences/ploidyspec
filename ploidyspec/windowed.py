import csv
import math
import os
import shutil
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from .common import group_pair_batches, log, run, run_logex_batch, sum_hist_distinct
from .kmer_tables import load_sequences

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
    "kmers_b",
    "shared",
    "union",
    "jaccard_distance",
]


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


def build_window_ktab(samtools_bin, fastk_bin, k, source, seq_id, start, end, tmp_dir, unit_id):
    region = f"{seq_id}:{start}-{end}"
    fa_path = os.path.join(tmp_dir, f"{unit_id}.fa")
    with open(fa_path, "w") as out:
        run([samtools_bin, "faidx", source, region], stdout=out)
    ktab_prefix = os.path.join(tmp_dir, unit_id)
    # -P scopes FastK's block-sort scratch files to this window's own tmp_dir: the same unit_id
    # recurs across many concurrently-processed windows, and FastK names scratch files from the
    # basename only, so sharing $TMPDIR across windows causes cross-window filename collisions.
    run([fastk_bin, f"-k{k}", "-t1", "-T1", f"-N{ktab_prefix}", f"-P{tmp_dir}", fa_path])
    return ktab_prefix


def process_window(label, win_idx, start, end, group, samtools_bin, fastk_bin, logex_bin, histex_bin, k, base_tmp):
    tmp_dir = os.path.join(base_tmp, f"w{win_idx}")
    os.makedirs(tmp_dir, exist_ok=True)
    try:
        prefixes = []
        totals = []
        for u in group:
            p = build_window_ktab(samtools_bin, fastk_bin, k, u["source"], u["seq_id"], start, end, tmp_dir, u["unit_id"])
            prefixes.append(p)
            totals.append(sum_hist_distinct(histex_bin, p))

        n = len(group)
        shared = {}
        for idxs, local_pairs in group_pair_batches(n, letter_cap=8):
            source_prefixes = [prefixes[i] for i in idxs]
            batch_tmp = os.path.join(tmp_dir, "lx")
            res = run_logex_batch(logex_bin, histex_bin, source_prefixes, local_pairs, batch_tmp)
            for (a, b), count in res.items():
                shared[(idxs[a], idxs[b])] = count

        rows = []
        for (i, j), count in shared.items():
            union = totals[i] + totals[j] - count
            dist = 1.0 - count / union if union > 0 else 0.0
            rows.append(
                dict(
                    group=label,
                    win_start=start,
                    win_end=end,
                    unit_a=group[i]["unit_id"],
                    unit_b=group[j]["unit_id"],
                    hap_a=group[i]["hap"],
                    hap_b=group[j]["hap"],
                    chrom_a=int(group[i]["chrom"]),
                    chrom_b=int(group[j]["chrom"]),
                    kmers_a=totals[i],
                    kmers_b=totals[j],
                    shared=count,
                    union=union,
                    jaccard_distance=dist,
                )
            )
        return rows
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def compute_windowed_groups(
    labeled_groups, outdir, samtools_bin, fastk_bin, logex_bin, histex_bin, k, window, step, threads,
    overview_title, all_tsv_name, overview_png_name, only_cross_chrom=False,
):
    """
    Core windowed-comparison loop: for each label -> list-of-units group, tile the
    shortest member into windows and compute pairwise k-mer Jaccard distance at
    every window between every pair of units in that group. Used both for same-
    chromosome-number haplotype comparisons and for cross-chromosome homeolog pairs.
    """
    rows_by_label = {}
    for label in sorted(labeled_groups):
        group = sorted(labeled_groups[label], key=lambda u: (u["hap"], int(u["chrom"])))
        if len(group) < 2:
            log(f"{label}: only {len(group)} unit present, skipping windowed comparison")
            continue

        lengths = [int(u["length"]) for u in group]
        minlen = min(lengths)
        if max(lengths) > 0 and (max(lengths) - minlen) / max(lengths) > 0.2:
            log(
                f"{label}: member lengths vary by >20% (min={minlen}, max={max(lengths)}); "
                f"windows are bounded by the shortest member, assuming rough colinearity"
            )

        windows = make_windows(minlen, window, step)
        log(f"{label}: {len(group)} units x {len(windows)} windows ({window}bp, step {step}bp)")

        tmp_root = os.path.join(outdir, "tmp_windowed", label)
        os.makedirs(tmp_root, exist_ok=True)
        label_rows = []
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [
                ex.submit(process_window, label, wi, s, e, group, samtools_bin, fastk_bin, logex_bin, histex_bin, k, tmp_root)
                for wi, (s, e) in enumerate(windows)
            ]
            done = 0
            report_every = max(1, len(windows) // 10)
            for fut in as_completed(futs):
                label_rows.extend(fut.result())
                done += 1
                if done % report_every == 0 or done == len(windows):
                    log(f"  {label}: [{done}/{len(windows)}] windows done")
        shutil.rmtree(tmp_root, ignore_errors=True)

        if only_cross_chrom:
            # drop same-chromosome-number (already-known haplotype homolog) pairs -- those are
            # covered by the `windowed` stage; here we only want the cross-number ancestral signal
            label_rows = [r for r in label_rows if r["chrom_a"] != r["chrom_b"]]

        write_group_tsv(outdir, label, label_rows)
        plot_group(outdir, label, label_rows)
        rows_by_label[label] = label_rows

    write_all_windows_tsv(outdir, rows_by_label, all_tsv_name)
    if rows_by_label:
        plot_overview(outdir, rows_by_label, overview_title, overview_png_name)
    return rows_by_label


def compute_windowed(seq_tsv, outdir, samtools_bin, fastk_bin, logex_bin, histex_bin, k, window, step, threads):
    units = load_sequences(seq_tsv)
    groups = defaultdict(list)
    for u in units:
        groups[int(u["chrom"])].append(u)
    labeled_groups = {f"chr{c:02d}": g for c, g in groups.items()}
    return compute_windowed_groups(
        labeled_groups, outdir, samtools_bin, fastk_bin, logex_bin, histex_bin, k, window, step, threads,
        overview_title="Genome-wide windowed haplotype k-mer divergence (relative chromosome position)",
        all_tsv_name="windowed_all.tsv",
        overview_png_name="windowed_genome_overview.png",
    )


def compute_windowed_homeologs(seq_tsv, outdir, homeolog_pairs, samtools_bin, fastk_bin, logex_bin, histex_bin, k, window, step, threads):
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
        labeled_groups, outdir, samtools_bin, fastk_bin, logex_bin, histex_bin, k, window, step, threads,
        overview_title="Windowed k-mer divergence between candidate ancestral (paleopolyploid) chromosome pairs",
        all_tsv_name="windowed_homeologs_all.tsv",
        overview_png_name="windowed_homeologs_overview.png",
        only_cross_chrom=True,
    )


def write_group_tsv(outdir, label, rows):
    path = os.path.join(outdir, f"windowed_{label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["win_start"], r["hap_a"], r["hap_b"])):
            w.writerow({k: (f"{v:.6f}" if k == "jaccard_distance" else v) for k, v in r.items()})


def write_all_windows_tsv(outdir, rows_by_label, filename):
    path = os.path.join(outdir, filename)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader()
        for label in sorted(rows_by_label):
            for r in sorted(rows_by_label[label], key=lambda r: (r["win_start"], r["hap_a"], r["hap_b"])):
                w.writerow({k: (f"{v:.6f}" if k == "jaccard_distance" else v) for k, v in r.items()})


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
        ys = [r["jaccard_distance"] for r in sub]
        ax.plot(xs, ys, marker="o", markersize=2, linewidth=1, label=lbl)
    ax.set_xlabel("position (Mb)")
    ax.set_ylabel("windowed k-mer Jaccard distance")
    ax.set_title(f"{label} windowed divergence")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"windowed_{label}.png"), dpi=150)
    plt.close(fig)


def plot_overview(outdir, rows_by_label, title, filename):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = sorted(rows_by_label)
    all_rows = [r for rows in rows_by_label.values() for r in rows]
    all_keys = sorted(set(_pair_key_and_label(r) for r in all_rows), key=lambda kl: kl[1])
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
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 2.5 * nrows), sharey=True, squeeze=False)
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
            ys = [r["jaccard_distance"] for r in sub]
            ax.plot(xs, ys, linewidth=1, color=color_of[key], label=lbl if not same_chrom_mode else None)
        ax.set_title(label, fontsize=9)
        ax.set_ylim(0, 1)
        if not same_chrom_mode:
            ax.legend(fontsize=5, loc="lower left")

    for ax in axes[len(labels):]:
        ax.axis("off")

    fig.suptitle(title)
    if same_chrom_mode:
        handles = [plt.Line2D([0], [0], color=color_of[key], label=lbl) for key, lbl in all_keys]
        fig.legend(handles=handles, loc="lower center", ncol=min(len(all_keys), 6), fontsize=7)
        fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    else:
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(outdir, filename), dpi=150)
    plt.close(fig)
