import csv
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

from .common import group_pair_batches, log, run_logex_batch, sum_hist_distinct
from .kmer_tables import ktab_prefix_path, load_sequences
from .mash import summarize_pair_across_k


def whole_chromosome_totals(units, outdir, histex_bin, k):
    totals = {}
    for u in units:
        prefix = ktab_prefix_path(outdir, u["unit_id"], k)
        totals[u["unit_id"]] = sum_hist_distinct(histex_bin, prefix)
    return totals


def _pairwise_shared_for_k(units, outdir, logex_bin, histex_bin, threads, k):
    """Whole-chromosome pairwise shared-kmer counts + per-unit totals at one k."""
    n = len(units)
    ids = [u["unit_id"] for u in units]
    prefixes = [ktab_prefix_path(outdir, uid, k) for uid in ids]

    totals = whole_chromosome_totals(units, outdir, histex_bin, k)

    batches = group_pair_batches(n, letter_cap=8)
    log(f"  k={k}: batched into {len(batches)} Logex calls")

    shared = np.zeros((n, n), dtype=np.int64)
    tmp_root = os.path.join(outdir, "tmp_matrix", f"k{k}")
    os.makedirs(tmp_root, exist_ok=True)

    def do_batch(batch_idx, idxs, local_pairs):
        source_prefixes = [prefixes[i] for i in idxs]
        tmp_dir = os.path.join(tmp_root, f"batch{batch_idx}")
        result = run_logex_batch(
            logex_bin, histex_bin, source_prefixes, local_pairs, tmp_dir
        )
        return idxs, result

    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = [
            ex.submit(do_batch, bi, idxs, pairs)
            for bi, (idxs, pairs) in enumerate(batches)
        ]
        done = 0
        for fut in as_completed(futs):
            idxs, result = fut.result()
            for (a, b), count in result.items():
                gi, gj = idxs[a], idxs[b]
                shared[gi, gj] = count
                shared[gj, gi] = count
            done += 1
            if done % max(1, len(batches) // 20) == 0 or done == len(batches):
                log(f"  k={k}: [{done}/{len(batches)}] batches done")
    if os.path.isdir(tmp_root) and not os.listdir(tmp_root):
        os.rmdir(tmp_root)

    total_arr = np.array([totals[uid] for uid in ids], dtype=np.int64)
    return ids, shared, total_arr


def compute_matrix(seq_tsv, outdir, logex_bin, histex_bin, threads, k_values):
    """
    All-vs-all whole-chromosome k-mer distance, swept across k_values and combined
    into a Mash-corrected consensus distance per pair (see mash.py). Pairs whose
    shared-kmer count doesn't clear the chance-collision floor at any swept k are
    reported as resolution_limited rather than given a fabricated precise distance.
    """
    units = load_sequences(seq_tsv)
    n = len(units)
    if n < 2:
        raise SystemExit(
            f"only {n} chromosome-scale unit(s) in {seq_tsv} -- nothing to compare. "
            f"Check --chrom-regex/--hap-regex against this species' FASTA headers "
            f"(see unplaced.tsv for why sequences were excluded)."
        )

    log(
        f"computing whole-chromosome pairwise k-mer intersections for {n} units "
        f"({n * (n - 1) // 2} pairs) across k={','.join(str(k) for k in k_values)}"
    )

    per_k = {}
    ids = None
    for k in k_values:
        k_ids, shared, total_arr = _pairwise_shared_for_k(
            units, outdir, logex_bin, histex_bin, threads, k
        )
        if ids is None:
            ids = k_ids
        per_k[k] = dict(shared=shared, total=total_arr)

    distance = np.zeros((n, n))
    summaries = {}
    n_resolution_limited = 0
    for i in range(n):
        for j in range(i + 1, n):
            stats = {
                k: (
                    int(per_k[k]["shared"][i, j]),
                    int(per_k[k]["total"][i]),
                    int(per_k[k]["total"][j]),
                )
                for k in k_values
            }
            summary = summarize_pair_across_k(k_values, stats)
            summaries[(i, j)] = summary
            d = summary["distance"] if summary["distance"] is not None else 1.0
            distance[i, j] = distance[j, i] = d
            if summary["resolution_limited"]:
                n_resolution_limited += 1
    np.fill_diagonal(distance, 0.0)

    total_pairs = n * (n - 1) // 2
    log(
        f"  {n_resolution_limited}/{total_pairs} pairs resolution-limited at every swept k "
        f"(no k in {k_values} cleared the chance-collision noise floor)"
    )

    write_outputs(outdir, ids, units, k_values, distance, summaries)
    return ids, distance, summaries


def write_outputs(outdir, ids, units, k_values, distance, summaries):
    n = len(ids)
    k_sweep_str = ",".join(str(k) for k in k_values)

    with open(
        os.path.join(outdir, "whole_chrom_distance_matrix.csv"), "w", newline=""
    ) as f:
        w = csv.writer(f)
        w.writerow([""] + ids)
        for i, uid in enumerate(ids):
            w.writerow([uid] + [f"{distance[i, j]:.6f}" for j in range(n)])

    with open(os.path.join(outdir, "whole_chrom_pairs.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "unit_a",
                "unit_b",
                "k_sweep",
                "chosen_k",
                "kmers_a",
                "kmers_b",
                "shared",
                "union",
                "jaccard_at_chosen_k",
                "containment_at_chosen_k",
                "distance",
                "containment_p",
                "resolution_limited",
                "k_consistency_spread",
            ]
        )
        for i in range(n):
            for j in range(i + 1, n):
                s = summaries[(i, j)]
                ref_k = s["chosen_k"] if s["chosen_k"] is not None else max(k_values)
                pk = s["per_k"][ref_k]
                union = pk["n1"] + pk["n2"] - pk["shared"]
                w.writerow(
                    [
                        ids[i],
                        ids[j],
                        k_sweep_str,
                        ref_k,
                        pk["n1"],
                        pk["n2"],
                        pk["shared"],
                        union,
                        f"{pk['jaccard']:.6f}",
                        f"{pk['containment']:.6f}",
                        "" if s["distance"] is None else f"{s['distance']:.6f}",
                        "" if s["containment_p"] is None else f"{s['containment_p']:.6f}",
                        s["resolution_limited"],
                        (
                            ""
                            if s["k_consistency_spread"] is None
                            else f"{s['k_consistency_spread']:.6f}"
                        ),
                    ]
                )

    with open(os.path.join(outdir, "whole_chrom_multi_k.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "unit_a",
                "unit_b",
                "k",
                "shared",
                "kmers_a",
                "kmers_b",
                "jaccard",
                "containment",
                "distance",
                "resolution_limited",
                "z",
            ]
        )
        for i in range(n):
            for j in range(i + 1, n):
                s = summaries[(i, j)]
                for k in k_values:
                    pk = s["per_k"][k]
                    w.writerow(
                        [
                            ids[i],
                            ids[j],
                            k,
                            pk["shared"],
                            pk["n1"],
                            pk["n2"],
                            f"{pk['jaccard']:.6f}",
                            f"{pk['containment']:.6f}",
                            "" if pk["distance"] is None else f"{pk['distance']:.6f}",
                            pk["resolution_limited"],
                            "" if pk["z"] is None else f"{pk['z']:.4f}",
                        ]
                    )

    plot_heatmap(outdir, ids, units, distance)
    plot_k_resolution_diagnostic(outdir, k_values, summaries)


def average_linkage_order(distance):
    """Pure-numpy agglomerative (average-linkage) leaf ordering -- no scipy dependency."""
    n = distance.shape[0]
    if n <= 2:
        return list(range(n))

    dist = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[(i, j)] = float(distance[i, j])

    def get_d(a, b):
        return dist[(a, b)] if a < b else dist[(b, a)]

    size = {i: 1 for i in range(n)}
    children = {}
    active = list(range(n))
    next_id = n
    while len(active) > 1:
        best = None
        for ii in range(len(active)):
            for jj in range(ii + 1, len(active)):
                a, b = active[ii], active[jj]
                d = get_d(a, b)
                if best is None or d < best[0]:
                    best = (d, a, b)
        _, a, b = best
        na, nb = size[a], size[b]
        new_id = next_id
        next_id += 1
        for c in active:
            if c == a or c == b:
                continue
            newd = (na * get_d(a, c) + nb * get_d(b, c)) / (na + nb)
            key = (new_id, c) if new_id < c else (c, new_id)
            dist[key] = newd
        size[new_id] = na + nb
        children[new_id] = (a, b)
        active = [c for c in active if c not in (a, b)] + [new_id]

    order = []

    def visit(node):
        if node < n:
            order.append(node)
        else:
            left, right = children[node]
            visit(left)
            visit(right)

    visit(active[0])
    return order


def plot_heatmap(outdir, ids, units, distance):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(ids)
    order = average_linkage_order(distance)

    ordered_ids = [ids[i] for i in order]
    ordered_mat = distance[np.ix_(order, order)]

    fig_w = max(6, n * 0.25)
    fig, ax = plt.subplots(figsize=(fig_w, fig_w))
    im = ax.imshow(ordered_mat, cmap="viridis_r", vmin=0, vmax=1)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(ordered_ids, rotation=90, fontsize=max(3, 8 - n // 20))
    ax.set_yticklabels(ordered_ids, fontsize=max(3, 8 - n // 20))
    ax.set_title(
        "Whole-chromosome Mash-corrected k-mer divergence (multi-k consensus)"
    )
    fig.colorbar(im, ax=ax, shrink=0.7, label="distance")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "whole_chrom_distance_heatmap.png"), dpi=150)
    plt.close(fig)


def plot_k_resolution_diagnostic(outdir, k_values, summaries):
    """Per-pair distance-vs-k, marking where pairs cross into resolution-limited
    territory -- the empirical version of "how far back can this k-mer sweep see"."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ks_sorted = sorted(k_values)
    fig, ax = plt.subplots(figsize=(7, 5))
    for s in summaries.values():
        xs, ys, colors = [], [], []
        for k in ks_sorted:
            pk = s["per_k"][k]
            if pk["distance"] is None:
                continue
            xs.append(k)
            ys.append(pk["distance"])
            colors.append("crimson" if pk["resolution_limited"] else "steelblue")
        if len(xs) < 2:
            continue
        ax.plot(xs, ys, color="0.75", linewidth=0.6, alpha=0.5, zorder=1)
        ax.scatter(xs, ys, c=colors, s=6, alpha=0.6, zorder=2)

    ax.set_xlabel("k")
    ax.set_ylabel("Mash-corrected distance")
    ax.set_title(
        "Per-pair distance vs k -- red = below the chance-collision noise floor at that k"
    )
    ax.set_xticks(ks_sorted)
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "k_resolution_diagnostic.png"), dpi=150)
    plt.close(fig)
