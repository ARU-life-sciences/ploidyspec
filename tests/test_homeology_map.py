import csv
import os
import random
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.homeology_map import (BLOCKS_TSV, WINDOWS_TSV, analyse_map, label_runs,
                                      partner_pairs, read_summary)
from ploidyspec.windowed import FIELDNAMES

W = 250_000
N_WIN = 40


def build(d, homeology, n_chrom=8, seed=0, low=(), near=()):
    """homeology: {(anchor_chrom, window_index): partner_chrom}, symmetric
    applied by the caller. Homeolog windows 0.10 from their partner, the rest
    0.30 +- noise. Anchors in `low` sit 0.03 from their partner, in `near` 0.01."""
    rng = random.Random(seed)
    os.makedirs(os.path.join(d, "homeologs"))
    with open(os.path.join(d, "homeologs", WINDOWS_TSV), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader()
        for a in range(1, n_chrom + 1):
            for i in range(N_WIN):
                for b in range(1, n_chrom + 1):
                    if a == b:
                        continue
                    close = 0.01 if a in near else 0.03 if a in low else 0.10
                    dist = close if homeology.get((a, i)) == b else 0.30 + rng.uniform(-0.02, 0.02)
                    w.writerow(dict(group="map", win_start=i * W + 1, win_end=(i + 1) * W,
                                    unit_a=f"HAP1_chr{a:02d}", unit_b=f"HAP1_chr{b:02d}", hap_a="HAP1",
                                    hap_b="HAP1", chrom_a=a, chrom_b=b, kmers_a=W, shared=0,
                                    containment=0, distance=dist))


class TestLabelRuns(unittest.TestCase):
    def test_runs_bridge_short_gaps_of_other_labels(self):
        labels = [2, 2, None, 2, 3, 3, 3, None, None, None, 3]
        self.assertEqual([(r[0], r[3]) for r in label_runs(labels)], [(2, 3), (3, 3), (3, 1)])


class TestHomeologyMap(unittest.TestCase):
    def test_arm_level_homeology_and_a_whole_chromosome_pair(self):
        # chr01's first half is homeologous to chr02, its second half to chr03
        # (a fusion since the duplication); chr04 and chr05 are a whole pair
        h = {}
        for i in range(N_WIN):
            p = 2 if i < N_WIN // 2 else 3
            h[(1, i)] = p
            h[(p, i % (N_WIN // 2) + (0 if p == 2 else N_WIN // 2))] = 1
            h[(4, i)], h[(5, i)] = 5, 4
        with tempfile.TemporaryDirectory() as d:
            build(d, h)
            s = analyse_map(d, min_segment_bp=1_000_000)
            self.assertEqual(s["multi_partner_list"], "chr01:chr02+chr03")
            self.assertEqual(sorted(partner_pairs(d)), [(1, 2), (1, 3), (4, 5)])
            self.assertEqual(read_summary(d)["n_partner_pairs"], "3")
            with open(os.path.join(d, "homeologs", BLOCKS_TSV)) as f:
                blocks = list(csv.DictReader(f, delimiter="\t"))
            whole = [b for b in blocks if b["chrom"] == "chr04"]
            self.assertEqual([b["position"] for b in whole], ["whole"])
            self.assertTrue(all(b["reciprocal"] == "yes" for b in blocks))

    def test_low_divergence_pair_is_lore_like(self):
        # five whole pairs; chr09/chr10 diverged far less than the rest
        h = {}
        for a in range(1, 11, 2):
            for i in range(N_WIN):
                h[(a, i)], h[(a + 1, i)] = a + 1, a
        with tempfile.TemporaryDirectory() as d:
            build(d, h, n_chrom=10, low=(9, 10))
            s = analyse_map(d, min_segment_bp=1_000_000)
            self.assertEqual(s["n_partner_pairs"], 5)
            self.assertEqual(s["n_lore_blocks"], 2)
            self.assertEqual(s["lore_pairs"], "chr09xchr10")
            self.assertEqual(s["lore_frac"], "0.200")
            with open(os.path.join(d, "homeologs", BLOCKS_TSV)) as f:
                low = {b["chrom"] for b in csv.DictReader(f, delimiter="\t") if b["low_divergence"] == "yes"}
            self.assertEqual(low, {"chr09", "chr10"})

    def test_near_allelic_pair_is_not_lore_like(self):
        # chr09/chr10 as close as alleles (0.1x the median): residual territory
        h = {}
        for a in range(1, 11, 2):
            for i in range(N_WIN):
                h[(a, i)], h[(a + 1, i)] = a + 1, a
        with tempfile.TemporaryDirectory() as d:
            build(d, h, n_chrom=10, near=(9, 10))
            s = analyse_map(d, min_segment_bp=1_000_000)
            self.assertEqual(s["n_lore_blocks"], 0)
            self.assertEqual(s["n_near_allelic_blocks"], 2)

    def test_uniform_divergence_has_no_lore_blocks(self):
        h = {}
        for a in range(1, 11, 2):
            for i in range(N_WIN):
                h[(a, i)], h[(a + 1, i)] = a + 1, a
        with tempfile.TemporaryDirectory() as d:
            build(d, h, n_chrom=10)
            self.assertEqual(analyse_map(d, min_segment_bp=1_000_000)["n_lore_blocks"], 0)

    def test_cross_subgenome_map_skips_own_set(self):
        # chr01-04 and chr05-08 are subgenomes: every window is closest to its
        # own set (shared repeats, 0.05), its homeolog across sets at 0.10
        h = {(a, i): a + 4 for a in range(1, 5) for i in range(N_WIN)}
        h.update({(a + 4, i): a for a in range(1, 5) for i in range(N_WIN)})
        with tempfile.TemporaryDirectory() as d:
            build(d, h)
            path = os.path.join(d, "homeologs", WINDOWS_TSV)
            with open(path) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            for r in rows:
                if (int(r["chrom_a"]) <= 4) == (int(r["chrom_b"]) <= 4):
                    r["distance"] = "0.05"
            with open(path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
                w.writeheader()
                w.writerows(rows)
            self.assertEqual(analyse_map(d, min_segment_bp=1_000_000)["n_partner_pairs"], 0)
            with open(os.path.join(d, "homeologs", "subgenome_partition.tsv"), "w") as f:
                f.write("metric\tvalue\npair_test\tcross_subgenome\n"
                        "set_a\tchr01,chr02,chr03,chr04\nset_b\tchr05,chr06,chr07,chr08\n")
            analyse_map(d, min_segment_bp=1_000_000)
            self.assertEqual(sorted(partner_pairs(d)), [(1, 5), (2, 6), (3, 7), (4, 8)])

    def test_no_homeology_gives_no_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            build(d, {})
            s = analyse_map(d, min_segment_bp=1_000_000)
            self.assertEqual(s["n_blocks"], 0)
            self.assertEqual(s["duplicated_frac"], "0.000")


if __name__ == "__main__":
    unittest.main()
