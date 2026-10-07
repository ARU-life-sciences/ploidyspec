import os
import random
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.structure import (
    agglomerative_merge_sequence,
    null_distribution,
    pair_synchrony,
    separation_ratio,
)


def paired_distances(n_pairs, within=0.05, between=0.5):
    """chr01..chrNN where (1,2), (3,4), ... are close pairs."""
    chroms = [f"chr{i:02d}" for i in range(1, 2 * n_pairs + 1)]
    dist = {}
    for a in chroms:
        for b in chroms:
            if a != b:
                ia, ib = int(a[3:]), int(b[3:])
                dist[(a, b)] = within if (ia + 1) // 2 == (ib + 1) // 2 else between
    return dist, chroms


class TestPartitions(unittest.TestCase):
    def test_merge_sequence_recovers_pairs(self):
        dist, chroms = paired_distances(4)
        clusters = dict(agglomerative_merge_sequence(dist, chroms))[4]
        self.assertEqual(sorted(sorted(c) for c in clusters),
                         [["chr01", "chr02"], ["chr03", "chr04"], ["chr05", "chr06"], ["chr07", "chr08"]])

    def test_true_partition_beats_random_ones(self):
        dist, chroms = paired_distances(4)
        true = [["chr01", "chr02"], ["chr03", "chr04"], ["chr05", "chr06"], ["chr07", "chr08"]]
        null = null_distribution(dist, chroms, [2, 2, 2, 2], random.Random(1), n=50)
        self.assertGreater(separation_ratio(dist, true), max(null))


class TestPairSynchrony(unittest.TestCase):
    def test_uniform_pair_depth_has_low_cv(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "homeologs"))
            with open(os.path.join(d, "homeologs", "homeolog_pairs.tsv"), "w") as f:
                f.write("chrom_a\tchrom_b\tmean_distance\n")
                f.write("chr01\tchr02\t0.10\nchr03\tchr04\t0.10\nchr05\tchr06\t0.40\n")
            row = pair_synchrony(d, "sp")
        self.assertEqual(row["n_accepted_pairs"], 3)
        self.assertEqual(row["pair_depth_median"], "0.1000")
        self.assertAlmostEqual(float(row["pair_depth_cv"]), 0.7071, places=3)

    def test_no_pairs_gives_blank_values(self):
        with tempfile.TemporaryDirectory() as d:
            row = pair_synchrony(d, "sp")
        self.assertEqual((row["n_accepted_pairs"], row["pair_depth_cv"]), (0, ""))


if __name__ == "__main__":
    unittest.main()
