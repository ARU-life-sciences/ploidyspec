#!/usr/bin/env python3
"""
Prototype: genome-wide subgenome-specific repeat k-mers (SubPhaser-style),
with a shuffle null, homeologous sets across chromosome numbers, and
likelihood window painting.

Reads a finished ploidyspec run (sequences.tsv, chroms/<unit>.fa,
homeologs/homeolog_pairs.tsv) and writes <outdir>/te_subgenomes/.

1. Homeologous sets: all copies of one chromosome number, joined with the
   copies of its homeolog partners (homeolog_pairs.tsv), so the stage works
   whether homeologs share a number (wheat A/B as 1A/1B) or not (masu 1/3).
2. Differential k-mers per set: canonical k-mers whose copies split into a
   high group (each >= --min-high occurrences) and a low group at least
   --ratio times lower in density (per Mb; floor 1 occurrence), at the
   largest gap. --ratio defaults to 10 (near-absence), not 2.
3. Consistency: k-mers differential in >= 2 sets. Copies are clustered on
   the phi correlation of their high/low calls over shared k-mers; a k-mer is
   consistent when all its high copies, in every set, fall in one cluster that
   holds none of its low copies. Subgenome-specific repeat history predicts
   the same k-mers mark chr01-A, chr02-A, ...; independent per-set
   differences (autopolyploid haplotypes) do not.
4. Null: per k-mer and per set, the high/low calls are shuffled among the
   set's copies (keeps every per-set split, breaks cross-set linkage), then
   re-clustered and re-scored. z and empirical p per k.
5. Painting: each consistent k-mer marks its cluster. Windows of every copy
   are classified by Poisson likelihood of their marker counts under each
   cluster's rates (rates measured on the copies of that cluster), so marker
   pool sizes cancel. Runs of >= --min-run windows of a foreign cluster are
   candidate exchanges (or phasing/scaffolding switches).

Caveat: markers are learned from the same copies they paint; a short foreign
segment barely moves its copy's genome-wide calls, but painting is not
leave-one-out.
"""

import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

LUT = np.full(256, 4, np.uint8)
for i, c in enumerate(b"ACGT"):
    LUT[c] = i
    LUT[c + 32] = i
SENTINEL = np.uint32(0xFFFFFFFF)
CHUNK = 20_000_000

G = {}  # worker globals (set before fork)


def log(msg):
    print(time.strftime("[%H:%M:%S] ") + msg, file=sys.stderr, flush=True)


def read_fasta_codes(path, seq_id):
    raw = np.fromfile(path, dtype=np.uint8)
    nl = np.flatnonzero(raw == 10)
    header = bytes(raw[1:nl[0]]).split()[0].decode()
    if header != seq_id:
        raise SystemExit(f"{path}: header {header} != {seq_id} (stale chroms cache?)")
    body = raw[nl[0] + 1:]
    body = body[(body != 10) & (body != 13)]
    return LUT[body]


def canonical_codes(b, k):
    """canonical k-mer code per position (uint32, k <= 16); SENTINEL where the k-mer has N."""
    n = len(b) - k + 1
    out = np.empty(max(n, 0), np.uint32)
    for s in range(0, n, CHUNK):
        e = min(n, s + CHUNK)
        m = e - s
        fw = np.zeros(m, np.uint32)
        rc = np.zeros(m, np.uint32)
        bad = np.zeros(m, bool)
        for j in range(k):
            x = b[s + j:s + j + m]
            bad |= x > 3
            xv = (x & 3).astype(np.uint32)
            fw = (fw << np.uint32(2)) | xv
            rc |= (np.uint32(3) - xv) << np.uint32(2 * j)
        c = np.minimum(fw, rc)
        c[bad] = SENTINEL
        out[s:e] = c
    return out


def sorted_counts(codes):
    c = codes[codes != SENTINEL]
    c.sort()
    if len(c) == 0:
        return c, np.zeros(0, np.int64)
    starts = np.concatenate(([0], np.flatnonzero(np.diff(c)) + 1))
    counts = np.diff(np.concatenate((starts, [len(c)])))
    return c[starts], counts


def unit_codes(u):
    b = read_fasta_codes(u["fa"], u["seq_id"])
    return canonical_codes(b, G["k"])


def job_candidates(u):
    uniq, cnt = sorted_counts(unit_codes(u))
    return uniq[cnt >= G["min_high"]]


def job_lookup(u):
    codes = unit_codes(u)
    if u.get("ivals"):
        keep = np.zeros(len(codes), bool)
        for a, b in u["ivals"]:
            keep[a:max(a, b - G["k"] + 1)] = True
        codes = codes[keep]
    uniq, cnt = sorted_counts(codes)
    cand = G["cand"]
    pos = np.searchsorted(uniq, cand)
    pos[pos >= len(uniq)] = 0
    hit = uniq[pos] == cand
    out = np.zeros(len(cand), np.uint32)
    out[hit] = cnt[pos[hit]]
    return out


def job_paint(u):
    codes = unit_codes(u)
    allm, gid, W, ng = G["allm"], G["gid"], G["window"], G["ngroups"]
    pos = np.searchsorted(allm, codes)
    pos[pos >= len(allm)] = 0
    hit = (allm[pos] == codes) & (codes != SENTINEL)
    where = np.flatnonzero(hit)
    nwin = len(codes) // W + 1
    counts = np.bincount((where // W) * ng + gid[pos[where]], minlength=nwin * ng)
    return counts.reshape(nwin, ng)


def load_units(outdir):
    with open(os.path.join(outdir, "sequences.tsv")) as f:
        units = list(csv.DictReader(f, delimiter="\t"))
    for u in units:
        u["chrom"] = int(u["chrom"])
        u["length"] = int(u["length"])
        u["fa"] = os.path.join(outdir, "chroms", u["unit_id"] + ".fa")
    return units


def homeologous_sets(outdir, units):
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u in units:
        find(u["chrom"])
    pairs = os.path.join(outdir, "homeologs", "homeolog_pairs.tsv")
    n_pairs = 0
    if os.path.exists(pairs):
        with open(pairs) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                a, b = int(r["chrom_a"].replace("chr", "")), int(r["chrom_b"].replace("chr", ""))
                if a in parent and b in parent:
                    parent[find(a)] = find(b)
                    n_pairs += 1
    groups = defaultdict(list)
    for i, u in enumerate(units):
        groups[find(u["chrom"])].append(i)
    sets = [sorted(v) for v in groups.values() if len(v) >= 2]
    return sorted(sets, key=lambda s: min(units[i]["chrom"] for i in s)), n_pairs


def differential(counts, lengths_mb, sets, min_high, ratio):
    """int8 calls (units x kmers): +1 high, -1 low, 0 not differential in that unit's set."""
    n, C = counts.shape
    X = np.zeros((n, C), np.int8)
    nsets = np.zeros(C, np.int16)
    for s in sets:
        idx = np.array(s)
        raw = counts[idx].astype(np.float32)
        L = lengths_mb[idx][:, None].astype(np.float32)
        dens = np.maximum(raw, 1.0) / L  # floor of one occurrence
        order = np.argsort(dens, axis=0)
        sd = np.take_along_axis(dens, order, axis=0)
        gaps = sd[1:] / sd[:-1]
        p = np.argmax(gaps, axis=0)
        best = gaps[p, np.arange(C)]
        thresh = sd[p + 1, np.arange(C)]
        high = dens >= thresh[None, :]
        ok = (best >= ratio) & np.all(np.where(high, raw >= min_high, True), axis=0)
        calls = np.where(high, 1, -1).astype(np.int8)
        calls[:, ~ok] = 0
        X[idx] = calls
        nsets += ok
    return X, nsets


def phi_matrix(X):
    a = (X == 1).astype(np.float32)
    b = (X == -1).astype(np.float32)
    n11, n00 = a @ a.T, b @ b.T
    n10, n01 = a @ b.T, b @ a.T
    den = np.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = (n11 * n00 - n10 * n01) / den
    return np.nan_to_num(phi, nan=0.0)


def cluster(phi, k):
    d = np.clip(1.0 - phi, 0.0, 2.0)
    np.fill_diagonal(d, 0.0)
    d = (d + d.T) / 2
    Z = linkage(squareform(d, checks=False), method="average")
    return fcluster(Z, k, criterion="maxclust") - 1


def consistent_mask(X, labels):
    a = (X == 1)
    b = (X == -1)
    high_tot = a.sum(0)
    out = np.zeros(X.shape[1], bool)
    owner = np.full(X.shape[1], -1, np.int16)
    for c in np.unique(labels):
        inC = labels == c
        hin = a[inC].sum(0)
        lin = b[inC].sum(0)
        ok = (hin == high_tot) & (high_tot > 0) & (lin == 0)
        owner[ok & ~out] = c
        out |= ok
    return out, owner


def shuffle_within_sets(X, sets, rng):
    Y = X.copy()
    for s in sets:
        idx = np.array(s)
        block = X[idx]
        perm = np.argsort(rng.random(block.shape), axis=0)
        Y[idx] = np.take_along_axis(block, perm, axis=0)
    return Y


def load_distance(outdir):
    path = os.path.join(outdir, "matrix", "whole_chrom_distance_matrix.csv")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        rows = list(csv.reader(f))
    names = rows[0][1:]
    return {(r[0], names[j]): float(v) for r in rows[1:] for j, v in enumerate(r[1:]) if v}


def distance_agreement(sets, units, sub, labels, D):
    """per chromosome number with >= 3 copies that the TE clusters split: does cutting
    those copies' whole-chromosome distances into the same number of groups give the
    same partition? (Within one number only: across homeologs distance separates the
    numbers, not the subgenomes.)"""
    if D is None:
        return []
    pos = {i: r for r, i in enumerate(sub)}
    by_num = defaultdict(list)
    for i in sub:
        by_num[units[i]["chrom"]].append(i)
    groups = [(f"chr{num:02d}", s) for num, s in sorted(by_num.items()) if len(s) >= 3]
    # homeologs on different numbers with <= 2 copies each: compare across the set
    for s in sets:
        nums = {units[i]["chrom"] for i in s}
        if len(nums) >= 2 and all(len(by_num[n]) <= 2 for n in nums):
            groups.append(("+".join(f"chr{n:02d}" for n in sorted(nums)), s))
    out = []
    for si, s in groups:
        ids = [units[i]["unit_id"] for i in s]
        te = [labels[pos[i]] for i in s]
        ng = len(set(te))
        if ng < 2:
            continue
        try:
            d = np.array([[0.0 if a == b else D[(a, b)] for b in ids] for a in ids])
        except KeyError:
            continue
        d = (d + d.T) / 2
        dl = fcluster(linkage(squareform(d, checks=False), method="average"), ng, criterion="maxclust")
        part = lambda lab: sorted(sorted(x for x, l in zip(ids, lab) if l == g) for g in set(lab))
        out.append((si, part(te) == part(list(dl)), ";".join(",".join(g) for g in part(te)),
                    ";".join(",".join(g) for g in part(list(dl)))))
    return out


def length_confound(sets, units, sub, labels):
    """Sets spanning >= 2 chromosome numbers whose TE clusters split by number (each
    number in one cluster, numbers in different clusters): which cluster takes the
    shorter number? If one cluster takes it in nearly every such set, the clusters may
    be sorting chromosomes by size (repeat density per Mb differs with size), not by
    subgenome. Returns (n_sets, n_shorter_in_majority_cluster, median length ratio, rows)."""
    pos = {i: r for r, i in enumerate(sub)}
    rows = []
    for s in sets:
        by = defaultdict(list)
        for i in s:
            by[units[i]["chrom"]].append(i)
        if len(by) != 2:
            continue
        cl = {n: {labels[pos[i]] for i in ii} for n, ii in by.items()}
        if any(len(c) != 1 for c in cl.values()):
            continue
        (n1, c1), (n2, c2) = [(n, next(iter(c))) for n, c in cl.items()]
        if c1 == c2:
            continue
        L = {n: np.median([units[i]["length"] for i in by[n]]) for n in by}
        short, long_ = sorted(by, key=lambda n: L[n])
        rows.append((f"chr{short:02d}", f"chr{long_:02d}", L[long_] / L[short],
                     int(next(iter(cl[short]))), int(next(iter(cl[long_])))))
    if not rows:
        return 0, 0, float("nan"), rows
    counts = defaultdict(int)
    for r in rows:
        counts[r[3]] += 1
    return len(rows), max(counts.values()), float(np.median([r[2] for r in rows])), rows


def block_intervals(outdir, units, sets):
    """For copies in sets spanning >= 2 chromosome numbers: 0-based intervals of the
    homeology-map blocks between their number and the set's other numbers. Blocks are
    given on one anchor copy; other copies get them scaled by length (alleles are
    near-colinear). Copies with no block keep the whole chromosome."""
    path = os.path.join(outdir, "homeologs", "homeolog_blocks.tsv")
    if not os.path.exists(path):
        return 0
    blocks = defaultdict(list)  # (chrom, partner) -> [(start, end, anchor_len)]
    L = {u["unit_id"]: u["length"] for u in units}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            c, p = int(r["chrom"].replace("chr", "")), int(r["partner"].replace("chr", ""))
            blocks[(c, p)].append((int(r["start"]) - 1, int(r["end"]), L.get(r["anchor"])))
    n = 0
    for st in sets:
        nums = {units[i]["chrom"] for i in st}
        if len(nums) < 2:
            continue
        for i in st:
            u = units[i]
            iv = []
            for p in nums - {u["chrom"]}:
                for a, b, al in blocks.get((u["chrom"], p), []):
                    f = u["length"] / al if al else 1.0
                    iv.append((int(a * f), min(u["length"], int(b * f))))
            if not iv:
                continue
            iv.sort()
            merged = [list(iv[0])]
            for a, b in iv[1:]:
                if a <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            u["ivals"] = [tuple(x) for x in merged]
            u["block_bp"] = sum(b - a for a, b in merged)
            n += 1
    return n


def runs(labels_seq, home, min_run):
    """(start_idx, end_idx, cluster) of runs of >= min_run windows assigned to a foreign cluster."""
    out, i, n = [], 0, len(labels_seq)
    while i < n:
        c = labels_seq[i]
        if c >= 0 and c != home:
            j = i
            while j + 1 < n and labels_seq[j + 1] in (c, -2):  # tolerate 'none' windows inside a run
                j += 1
            while labels_seq[j] != c:
                j -= 1
            if sum(1 for x in labels_seq[i:j + 1] if x == c) >= min_run:
                out.append((i, j, c))
            i = j + 1
        else:
            i += 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--min-high", type=int, default=20)
    ap.add_argument("--ratio", type=float, default=10.0)
    ap.add_argument("--max-kmers", type=int, default=300_000, help="multi-set k-mers kept for clustering/null")
    ap.add_argument("--null", type=int, default=200)
    ap.add_argument("--max-k", type=int, default=4)
    ap.add_argument("--window", type=int, default=500_000)
    ap.add_argument("--llr", type=float, default=5.0, help="nats between best and second cluster")
    ap.add_argument("--min-run", type=int, default=3)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--blocks", action="store_true",
                    help="count copies in cross-number sets only over their homeology-map blocks "
                         "(writes te_subgenomes_blocks/)")
    ap.add_argument("--no-paint", action="store_true")
    ap.add_argument("--k-tolerance", type=float, default=0.9,
                    help="k = the smallest k whose consistent fraction is >= this times the best")
    ap.add_argument("--min-consistent", type=float, default=0.15,
                    help="provisional call threshold on the consistent fraction (panel: allo >= 0.17, auto <= 0.08)")
    ap.add_argument("--min-agree", type=float, default=2 / 3,
                    help="fraction of informative chromosomes whose TE split must match the distance split")
    ap.add_argument("--confound-frac", type=float, default=0.8)
    ap.add_argument("--confound-ratio", type=float, default=1.5,
                    help="length-confounded when one cluster takes the shorter number in >= --confound-frac "
                         "of split cross-number sets and their median length ratio is >= this")
    args = ap.parse_args()
    if args.k > 16:
        raise SystemExit("--k must be <= 16 (uint32 codes)")

    od = os.path.join(args.outdir, "te_subgenomes_blocks" if args.blocks else "te_subgenomes")
    os.makedirs(od, exist_ok=True)
    units = load_units(args.outdir)
    sets, n_pairs = homeologous_sets(args.outdir, units)
    if args.blocks:
        nb = block_intervals(args.outdir, units, sets)
        log(f"block-restricted counting for {nb} copies in cross-number sets")
    lengths_mb = np.array([u.get("block_bp", u["length"]) / 1e6 for u in units])
    log(f"{len(units)} units, {len(sets)} homeologous sets (sizes {sorted(len(s) for s in sets)}), {n_pairs} homeolog pairs used")
    with open(os.path.join(od, "sets.tsv"), "w") as f:
        f.write("set\tunits\n")
        for si, s in enumerate(sets):
            f.write(f"{si}\t{','.join(units[i]['unit_id'] for i in s)}\n")

    in_set = sorted({i for s in sets for i in s})
    G.update(k=args.k, min_high=args.min_high)

    cache = os.path.join(od, f"counts_k{args.k}_h{args.min_high}.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        cand, counts = z["cand"], z["counts"]
        log(f"loaded {cache}: {len(cand)} candidate k-mers")
    else:
        whole = os.path.join(args.outdir, "te_subgenomes", f"counts_k{args.k}_h{args.min_high}.npz")
        if args.blocks and os.path.exists(whole):
            cand = np.load(whole)["cand"]
            log(f"candidates from {whole}")
        else:
            log(f"pass 1: k-mers with >= {args.min_high} occurrences per copy")
            with Pool(args.procs) as p:
                parts = p.map(job_candidates, [units[i] for i in in_set], chunksize=1)
            cand = np.unique(np.concatenate(parts))
        log(f"{len(cand)} candidate k-mers; pass 2: their counts in every copy")
        G["cand"] = cand
        with Pool(args.procs) as p:
            rows = p.map(job_lookup, [units[i] for i in in_set], chunksize=1)
        counts = np.zeros((len(units), len(cand)), np.uint32)
        for i, r in zip(in_set, rows):
            counts[i] = r
        np.savez(cache, cand=cand, counts=counts)

    X, nsets = differential(counts, lengths_mb, sets, args.min_high, args.ratio)
    multi = np.flatnonzero(nsets >= 2)
    log(f"differential in >= 1 set: {(nsets >= 1).sum()}, in >= 2 sets: {len(multi)}")
    rng = np.random.default_rng(args.seed)
    if len(multi) > args.max_kmers:
        multi = np.sort(rng.choice(multi, args.max_kmers, replace=False))
    Xm = X[:, multi]
    with open(os.path.join(od, "per_set.tsv"), "w") as f:
        f.write("set\tn_units\tn_differential\n")
        for si, s in enumerate(sets):
            f.write(f"{si}\t{len(s)}\t{int((X[s[0]] != 0).sum())}\n")

    sub = np.array(in_set)
    phi = phi_matrix(Xm[sub])
    with open(os.path.join(od, "phi.tsv"), "w") as f:
        f.write("unit\t" + "\t".join(units[i]["unit_id"] for i in sub) + "\n")
        for r, i in enumerate(sub):
            f.write(units[i]["unit_id"] + "\t" + "\t".join(f"{v:.3f}" for v in phi[r]) + "\n")

    max_set = max(len(s) for s in sets) if sets else 0
    ks = list(range(2, max(2, min(args.max_k, max_set)) + 1))
    sets_local = [[int(np.searchsorted(sub, i)) for i in s] for s in sets]
    results = []
    for k in ks:
        labels = cluster(phi, k)
        cons, owner = consistent_mask(Xm[sub], labels)
        obs = cons.mean() if len(cons) else 0.0
        null = []
        for _ in range(args.null):
            Y = shuffle_within_sets(Xm[sub], sets_local, rng)
            lab = cluster(phi_matrix(Y), k)
            null.append(consistent_mask(Y, lab)[0].mean())
        null = np.array(null)
        sd = null.std()
        z = (obs - null.mean()) / sd if sd > 0 else float("nan")
        p = (1 + (null >= obs).sum()) / (1 + len(null))
        # does the grouping just follow the haplotype files?
        hap_groups = defaultdict(set)
        for r, i in enumerate(sub):
            hap_groups[labels[r]].add(units[i]["hap"])
        hap_pure = all(len(v) == 1 for v in hap_groups.values())
        results.append(dict(k=k, labels=labels, owner=owner, cons=cons, obs=obs, null_mean=null.mean(),
                            null_sd=null.std(), z=z, p=p, hap_pure=hap_pure,
                            sizes=np.bincount(labels, minlength=k).tolist()))
        results[-1]["null_max"] = float(null.max())
        log(f"k={k}: consistent {obs:.4f} vs null {null.mean():.4f}±{null.std():.4f} (max {null.max():.4f}) p={p:.3g} "
            f"sizes={results[-1]['sizes']} clusters_follow_hap_files={hap_pure}")

    if not results:
        log("no sets with >= 2 copies; nothing to cluster")
        return
    # z grows with k as the null tightens; take the smallest k that keeps most of the consistency
    top = max(r["obs"] for r in results)
    best = next(r for r in results if r["obs"] >= args.k_tolerance * top)
    D = load_distance(args.outdir)
    for r in results:
        agr = distance_agreement(sets, units, sub, r["labels"], D)
        r["agree"] = (sum(a[1] for a in agr), len(agr))
        n_ls, n_short, ratio, lrows = length_confound(sets, units, sub, r["labels"])
        r["length"] = (n_ls, n_short, ratio)
        r["confounded"] = n_ls >= 3 and n_short >= args.confound_frac * n_ls and ratio >= args.confound_ratio
        ok_agree = r["agree"][1] == 0 or r["agree"][0] >= args.min_agree * r["agree"][1]
        if r["obs"] < args.min_consistent:
            r["call"] = "none"
        elif r["obs"] <= r["null_max"]:
            r["call"] = "none_within_null"
        elif r["confounded"]:
            r["call"] = "none_length_confounded"
        elif not ok_agree:
            r["call"] = "none_disagrees_with_distance"
        else:
            r["call"] = "subgenome_specific_repeats"
        if r is best:
            best_agr, best_lrows = agr, lrows
    with open(os.path.join(od, "summary.tsv"), "w") as f:
        f.write("k\tchosen\tn_multiset_kmers\tconsistent_frac\tnull_mean\tnull_max\tp\tcluster_sizes\t"
                "clusters_follow_hap_files\tchroms_agreeing_with_distance\tlength_split_sets\t"
                "shorter_in_one_cluster\tmedian_length_ratio\tcall\n")
        for r in results:
            call = r["call"]
            f.write(f"{r['k']}\t{'yes' if r is best else 'no'}\t{Xm.shape[1]}\t{r['obs']:.4f}\t{r['null_mean']:.4f}\t"
                    f"{r['null_max']:.4f}\t{r['p']:.3g}\t{','.join(map(str, r['sizes']))}\t{r['hap_pure']}\t"
                    f"{r['agree'][0]}/{r['agree'][1]}\t{r['length'][0]}\t{r['length'][1]}\t"
                    f"{r['length'][2]:.2f}\t{call}\n")
    with open(os.path.join(od, "distance_agreement.tsv"), "w") as f:
        f.write("chrom\tagrees\tte_partition\tdistance_partition\n")
        for row in best_agr:
            f.write("\t".join(map(str, row)) + "\n")
    with open(os.path.join(od, "length_check.tsv"), "w") as f:
        f.write("shorter\tlonger\tlength_ratio\tshorter_cluster\tlonger_cluster\n")
        for a, b, rr, c1, c2 in best_lrows:
            f.write(f"{a}\t{b}\t{rr:.2f}\tS{c1 + 1}\tS{c2 + 1}\n")
    log(f"call: {best['call']}")
    log(f"chosen k={best['k']} (consistent {best['obs']:.4f}); chromosomes agreeing with distance {best['agree'][0]}/{best['agree'][1]}")
    k = best["k"]
    labels = best["labels"]
    with open(os.path.join(od, "clusters.tsv"), "w") as f:
        f.write("unit\thap\tchrom\tset\tcluster\n")
        set_of = {i: si for si, s in enumerate(sets) for i in s}
        for r, i in enumerate(sub):
            u = units[i]
            f.write(f"{u['unit_id']}\t{u['hap']}\tchr{u['chrom']:02d}\t{set_of[i]}\tS{labels[r] + 1}\n")

    # markers: consistent k-mers (from all multi-set k-mers, not just the subsample)
    full_multi = np.flatnonzero(nsets >= 2)
    cons_all, owner_all = consistent_mask(X[sub][:, full_multi], labels)
    mk = full_multi[cons_all]
    mown = owner_all[cons_all]
    log(f"subgenome markers (k={k}): " + ", ".join(f"S{c + 1}={int((mown == c).sum())}" for c in range(k)))
    if len(mk) == 0 or args.no_paint:
        return
    order = np.argsort(cand[mk])
    G.update(allm=cand[mk][order], gid=mown[order].astype(np.int64), window=args.window, ngroups=k)
    with Pool(args.procs) as p:
        painted = p.map(job_paint, [units[i] for i in sub], chunksize=1)

    # rates[g, h]: markers of g per bp in copies of cluster h
    tot = np.zeros((k, k))
    bp = np.zeros(k)
    for r, i in enumerate(sub):
        tot[:, labels[r]] += painted[r].sum(0)
        bp[labels[r]] += units[i]["length"]
    rates = (tot + 1.0) / bp[None, :]
    with open(os.path.join(od, "rates.tsv"), "w") as f:
        f.write("marker_set\t" + "\t".join(f"in_S{h + 1}_per_Mb" for h in range(k)) + "\n")
        for g in range(k):
            f.write(f"S{g + 1}\t" + "\t".join(f"{rates[g, h] * 1e6:.2f}" for h in range(k)) + "\n")

    W = args.window
    seg_rows, copy_rows = [], []
    with open(os.path.join(od, "windows.tsv"), "w") as fw:
        fw.write("unit\tstart\tend\t" + "\t".join(f"n_S{g + 1}" for g in range(k)) + "\tassigned\tllr\n")
        for r, i in enumerate(sub):
            u = units[i]
            home = labels[r]
            seq_labels = []
            for wi, cnt in enumerate(painted[r]):
                start = wi * W
                end = min(u["length"], start + W)
                if end <= start:
                    continue
                lam = rates * (end - start)  # (g, h)
                ll = (cnt[:, None] * np.log(lam) - lam).sum(0)
                o = np.argsort(ll)[::-1]
                llr = ll[o[0]] - ll[o[1]]
                if cnt.sum() == 0:
                    a, lab = "none", -2
                elif llr >= args.llr:
                    a, lab = f"S{o[0] + 1}", int(o[0])
                else:
                    a, lab = "ambiguous", -1
                seq_labels.append(lab)
                fw.write(f"{u['unit_id']}\t{start + 1}\t{end}\t" + "\t".join(map(str, cnt)) + f"\t{a}\t{llr:.1f}\n")
            n = len(seq_labels)
            nh = sum(1 for x in seq_labels if x == home)
            nf = sum(1 for x in seq_labels if x >= 0 and x != home)
            segs = runs(seq_labels, home, args.min_run)
            for s0, s1, c in segs:
                seg_rows.append((u["unit_id"], f"S{home + 1}", f"S{c + 1}", s0 * W, min(u["length"], (s1 + 1) * W)))
            copy_rows.append((u["unit_id"], u["hap"], f"chr{u['chrom']:02d}", f"S{home + 1}", n, nh / n if n else 0,
                              nf / n if n else 0, len(segs)))
    with open(os.path.join(od, "copies.tsv"), "w") as f:
        f.write("unit\thap\tchrom\tcluster\tn_windows\thome_frac\tforeign_frac\tn_foreign_runs\n")
        for row in copy_rows:
            f.write("\t".join(f"{x:.3f}" if isinstance(x, float) else str(x) for x in row) + "\n")
    with open(os.path.join(od, "foreign_segments.tsv"), "w") as f:
        f.write("unit\thome\tforeign\tstart\tend\n")
        for row in seg_rows:
            f.write("\t".join(map(str, row)) + "\n")
    hf = [c[5] for c in copy_rows]
    log(f"painting: median home fraction {np.median(hf):.3f}; {len(seg_rows)} foreign runs >= {args.min_run} windows")


if __name__ == "__main__":
    main()
