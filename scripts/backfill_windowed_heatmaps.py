#!/usr/bin/env python3
"""
One-off backfill: regenerate the new windowed-heatmap PNGs
(windowed_<label>_heatmap.png, windowed_genome_overview_heatmap.png,
windowed_homeologs_overview_heatmap.png) for species that already have
windowed/windowed-homeologs data on disk, without rerunning FastK.

Handles two windowed_all.tsv schemas found in the existing panel:
- current: has `group`/`chrom_a`/`chrom_b` columns.
- pre-migration (daBudDavi1, daGalBore1, daGleHede1's windowed/ only): has a
  single `chrom` column and predates the windowed-homeologs feature, so it's
  always same-chromosome-number data -- chrom_a == chrom_b == chrom.
"""
import csv
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.common import homeologs_dir, windowed_dir
from ploidyspec.windowed import plot_group_heatmap, plot_overview_heatmap


def load_rows(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)
    if not rows:
        return []
    old_schema = "chrom" in fieldnames and "chrom_a" not in fieldnames
    out = []
    for r in rows:
        r["win_start"] = int(r["win_start"])
        r["win_end"] = int(r["win_end"])
        r["jaccard_distance"] = float(r["jaccard_distance"])
        if old_schema:
            chrom = int(r["chrom"])
            r["chrom_a"] = chrom
            r["chrom_b"] = chrom
            r["group"] = f"chr{chrom:02d}"
        else:
            r["chrom_a"] = int(r["chrom_a"])
            r["chrom_b"] = int(r["chrom_b"])
        out.append(r)
    return out


def group_rows(rows):
    rows_by_label = defaultdict(list)
    for r in rows:
        rows_by_label[r["group"]].append(r)
    return dict(rows_by_label)


def backfill_one(outdir, all_tsv_path, target_dir, overview_title, overview_png_name):
    rows = load_rows(all_tsv_path)
    if not rows:
        return 0
    rows_by_label = group_rows(rows)
    for label, label_rows in rows_by_label.items():
        plot_group_heatmap(target_dir, label, label_rows)
    plot_overview_heatmap(target_dir, rows_by_label, overview_title, overview_png_name)
    return len(rows_by_label)


def main():
    results_root = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"
    )
    n_species = 0
    for outdir in sorted(glob.glob(os.path.join(results_root, "*"))):
        if not os.path.isdir(outdir):
            continue
        species = os.path.basename(outdir)
        touched = False

        wdir = windowed_dir(outdir)
        w_all = os.path.join(wdir, "windowed_all.tsv")
        if os.path.exists(w_all):
            n = backfill_one(
                outdir,
                w_all,
                wdir,
                "Genome-wide windowed haplotype k-mer divergence (relative chromosome position)",
                "windowed_genome_overview_heatmap.png",
            )
            if n:
                print(f"{species}: windowed/ -- {n} chromosome groups")
                touched = True

        hdir = homeologs_dir(outdir)
        h_all = os.path.join(hdir, "windowed_homeologs_all.tsv")
        if os.path.exists(h_all):
            n = backfill_one(
                outdir,
                h_all,
                hdir,
                "Windowed k-mer divergence between candidate ancestral (paleopolyploid) chromosome pairs",
                "windowed_homeologs_overview_heatmap.png",
            )
            if n:
                print(f"{species}: homeologs/ -- {n} ancestral-pair groups")
                touched = True

        if touched:
            n_species += 1

    print(f"\nbackfilled {n_species} species")


if __name__ == "__main__":
    main()
