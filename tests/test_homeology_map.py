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


def build(d, homeology, n_chrom=8, seed=0):
    """homeology: {(anchor_chrom, window_index): partner_chrom}, symmetric
    applied by the caller. Homeolog windows 0.10 from their partner, the rest
    0.30 +- noise."""
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
                    dist = 0.10 if homeology.get((a, i)) == b else 0.30 + rng.uniform(-0.02, 0.02)
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

    def test_no_homeology_gives_no_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            build(d, {})
            s = analyse_map(d, min_segment_bp=1_000_000)
            self.assertEqual(s["n_blocks"], 0)
            self.assertEqual(s["duplicated_frac"], "0.000")


if __name__ == "__main__":
    unittest.main()
