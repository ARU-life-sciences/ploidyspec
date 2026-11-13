#!/usr/bin/env python3
"""One row per species from scripts/te_subgenomes.py output: the chosen call of
the whole-chromosome and the blocks-only run, their scores, and which
chromosomes went to which cluster (the candidate subgenomes).

Columns
  species                 ToLID (results/<species>)
  call                    whole-chromosome call (te_subgenomes/summary.tsv, chosen k)
  k                       number of clusters chosen
  consistent_frac         share of multi-set differential k-mers that split copies consistently
  null_mean, null_max     same after shuffling calls among each set's copies
  excess_over_null        consistent_frac - null_max (> 0 = above every shuffle)
  p                       empirical p from the shuffles
  cluster_sizes           copies per cluster
  sets_agreeing           homeolog sets whose TE split matches the distance lineages
  blocks_call, blocks_consistent_frac, blocks_null_max, blocks_sets_agreeing
                          the same from te_subgenomes_blocks (homeology-map blocks only);
                          "n/a" when the blocks run restricted nothing (same k-mers and score as
                          the whole run: no set spans two numbers, or no homeology blocks)
  n_copies, n_numbers     copies (hap x chromosome) and chromosome numbers clustered
  split_numbers           chromosome numbers whose copies fall in different clusters
  painted_foreign_frac    mean share of windows painted as another cluster (whole run)
  clusters                membership, e.g. "S1: chr01,chr04 | S2: chr02,chr20[HAP1]";
                          a number is suffixed with [haps] when its copies are split

usage: scripts/te_subgenomes_panel.py [results_dir] [-o results/te_subgenomes_panel.tsv]
"""
import argparse
import csv
import glob
import os
from collections import defaultdict


def read(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def chosen(d):
    rows = read(os.path.join(d, "summary.tsv")) or []
    return next((r for r in rows if r["chosen"] == "yes"), None)


def membership(d):
    rows = read(os.path.join(d, "clusters.tsv")) or []
    by_num = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_num[r["chrom"]][r["cluster"]].append(r["hap"])
    out = defaultdict(list)
    for num in sorted(by_num):
        cl = by_num[num]
        for c, haps in cl.items():
            out[c].append(num if len(cl) == 1 else f"{num}[{'+'.join(sorted(haps))}]")
    split = [n for n in sorted(by_num) if len(by_num[n]) > 1]
    text = " | ".join(f"{c}: {','.join(out[c])}" for c in sorted(out))
    return len(rows), len(by_num), split, text


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?", default="results")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(a.results, "te_subgenomes_panel.tsv")
    cols = ["species", "call", "k", "consistent_frac", "null_mean", "null_max", "excess_over_null", "p",
            "cluster_sizes", "sets_agreeing", "blocks_call", "blocks_consistent_frac", "blocks_null_max",
            "blocks_sets_agreeing", "n_copies", "n_numbers", "split_numbers", "painted_foreign_frac", "clusters"]
    rows = []
    for d in sorted(glob.glob(os.path.join(a.results, "*", "te_subgenomes")), key=str.lower):
        sp = os.path.basename(os.path.dirname(d))
        w = chosen(d)
        if w is None:
            continue
        b = chosen(d + "_blocks") or {}
        if b and b.get("n_multiset_kmers") == w["n_multiset_kmers"] and b.get("consistent_frac") == w["consistent_frac"]:
            b = dict(call="n/a", consistent_frac="", null_max="", chroms_agreeing_with_distance="")
        cf, nm = num(w["consistent_frac"]), num(w["null_max"])
        n_copies, n_numbers, split, text = membership(d)
        copies = read(os.path.join(d, "copies.tsv")) or []
        ff = [num(r["foreign_frac"]) for r in copies if num(r["foreign_frac"]) is not None]
        rows.append(dict(
            species=sp, call=w["call"], k=w["k"], consistent_frac=w["consistent_frac"],
            null_mean=w["null_mean"], null_max=w["null_max"],
            excess_over_null=f"{cf - nm:.4f}" if cf is not None and nm is not None else "",
            p=w["p"], cluster_sizes=w["cluster_sizes"], sets_agreeing=w["chroms_agreeing_with_distance"],
            blocks_call=b.get("call", ""), blocks_consistent_frac=b.get("consistent_frac", ""),
            blocks_null_max=b.get("null_max", ""), blocks_sets_agreeing=b.get("chroms_agreeing_with_distance", ""),
            n_copies=n_copies, n_numbers=n_numbers, split_numbers=",".join(split),
            painted_foreign_frac=f"{sum(ff) / len(ff):.3f}" if ff else "", clusters=text))
    with open(out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)
    print(f"{out}: {len(rows)} species, "
          f"{sum(r['call'] == 'subgenome_specific_repeats' for r in rows)} whole-chromosome positives")


if __name__ == "__main__":
    main()
