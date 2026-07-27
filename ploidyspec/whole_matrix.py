import csv
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

from .common import group_pair_batches, log, run_logex_batch, sum_hist_distinct
from .kmer_tables import ktab_prefix_path, load_sequences


def whole_chromosome_totals(units, outdir, histex_bin):
    totals = {}
    for u in units:
        prefix = ktab_prefix_path(outdir, u["unit_id"])
        totals[u["unit_id"]] = sum_hist_distinct(histex_bin, prefix)
    return totals


def compute_matrix(seq_tsv, outdir, logex_bin, histex_bin, threads):
    units = load_sequences(seq_tsv)
    n = len(units)
    ids = [u["unit_id"] for u in units]
    prefixes = [ktab_prefix_path(outdir, uid) for uid in ids]

    log(f"computing whole-chromosome pairwise k-mer intersections for {n} units ({n * (n - 1) // 2} pairs)")
    totals = whole_chromosome_totals(units, outdir, histex_bin)

    batches = group_pair_batches(n, letter_cap=8)
    log(f"  batched into {len(batches)} Logex calls")

    shared = np.zeros((n, n), dtype=np.int64)
    tmp_root = os.path.join(outdir, "tmp_matrix")
    os.makedirs(tmp_root, exist_ok=True)

    def do_batch(batch_idx, idxs, local_pairs):
        source_prefixes = [prefixes[i] for i in idxs]
        tmp_dir = os.path.join(tmp_root, f"batch{batch_idx}")
        result = run_logex_batch(logex_bin, histex_bin, source_prefixes, local_pairs, tmp_dir)
        return idxs, result

    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = [ex.submit(do_batch, bi, idxs, pairs) for bi, (idxs, pairs) in enumerate(batches)]
        done = 0
        for fut in as_completed(futs):
            idxs, result = fut.result()
            for (a, b), count in result.items():
                gi, gj = idxs[a], idxs[b]
                shared[gi, gj] = count
                shared[gj, gi] = count
            done += 1
            if done % max(1, len(batches) // 20) == 0 or done == len(batches):
                log(f"  [{done}/{len(batches)}] batches done")
    os.rmdir(tmp_root) if not os.listdir(tmp_root) else None

    total_arr = np.array([totals[uid] for uid in ids], dtype=np.int64)
    union = total_arr[:, None] + total_arr[None, :] - shared
    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard = np.where(union > 0, shared / union, 0.0)
        min_total = np.minimum(total_arr[:, None], total_arr[None, :])
        containment = np.where(min_total > 0, shared / min_total, 0.0)
    np.fill_diagonal(jaccard, 1.0)
    np.fill_diagonal(containment, 1.0)
    distance = 1.0 - jaccard
    np.fill_diagonal(distance, 0.0)

    write_outputs(outdir, ids, units, total_arr, shared, distance, containment)
    return ids, distance, containment, total_arr


def write_outputs(outdir, ids, units, total_arr, shared, distance, containment):
    n = len(ids)

    with open(os.path.join(outdir, "whole_chrom_distance_matrix.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + ids)
        for i, uid in enumerate(ids):
            w.writerow([uid] + [f"{distance[i, j]:.6f}" for j in range(n)])

    with open(os.path.join(outdir, "whole_chrom_pairs.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["unit_a", "unit_b", "kmers_a", "kmers_b", "shared", "union", "jaccard_distance", "containment"])
        for i in range(n):
            for j in range(i + 1, n):
                union = int(total_arr[i] + total_arr[j] - shared[i, j])
                w.writerow(
                    [
                        ids[i],
                        ids[j],
                        int(total_arr[i]),
                        int(total_arr[j]),
                        int(shared[i, j]),
                        union,
                        f"{distance[i, j]:.6f}",
                        f"{containment[i, j]:.6f}",
                    ]
                )

    plot_heatmap(outdir, ids, units, distance)


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
    ax.set_title("Whole-chromosome k-mer Jaccard distance (1 - |A∩B| / |A∪B|)")
    fig.colorbar(im, ax=ax, shrink=0.7, label="distance")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "whole_chrom_distance_heatmap.png"), dpi=150)
    plt.close(fig)
