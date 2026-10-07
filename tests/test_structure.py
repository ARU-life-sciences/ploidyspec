import os
import random
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.structure import (
    agglomerative_merge_sequence,
    canonical_partition,
    metric_distance_ratio_cv,
    null_distribution,
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
    def test_canonical_partition_ignores_side_order(self):
        a = canonical_partition(["HAP1_chr01", "HAP2_chr01"], ["HAP3_chr01", "HAP4_chr01"])
        b = canonical_partition(["HAP4_chr02", "HAP3_chr02"], ["HAP2_chr02", "HAP1_chr02"])
        self.assertEqual(a, b)

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


class TestDistanceRatioCv(unittest.TestCase):
    def test_dedupes_pairs_listed_from_both_sides(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "homeologs"))
            with open(os.path.join(d, "homeologs", "ploidy_ancestry_summary.tsv"), "w") as f:
                f.write("chrom\thomeolog_partner\tdistance_ratio\n")
                f.write("chr01\tchr02\t2.0\nchr02\tchr01\t2.0\nchr03\tchr04\t4.0\nchr04\tchr03\t4.0\n")
            cv, n = metric_distance_ratio_cv(d)
        self.assertEqual(n, 2)
        self.assertAlmostEqual(cv, 1 / 3)


if __name__ == "__main__":
    unittest.main()
