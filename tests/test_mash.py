import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec import mash


class TestJaccardToDistance(unittest.TestCase):
    def test_identical_sequences(self):
        self.assertEqual(mash.jaccard_to_distance(1.0, 15), 0.0)

    def test_zero_or_negative_overlap_is_undefined(self):
        self.assertIsNone(mash.jaccard_to_distance(0.0, 15))
        self.assertIsNone(mash.jaccard_to_distance(-0.1, 15))

    def test_monotonic_in_jaccard(self):
        d_high_j = mash.jaccard_to_distance(0.9, 15)
        d_low_j = mash.jaccard_to_distance(0.1, 15)
        self.assertLess(d_high_j, d_low_j)

    def test_larger_k_gives_smaller_distance_for_same_jaccard(self):
        # the same observed Jaccard implies *less* true divergence if it was
        # measured with a larger, more specific k
        d_k11 = mash.jaccard_to_distance(0.5, 11)
        d_k25 = mash.jaccard_to_distance(0.5, 25)
        self.assertGreater(d_k11, d_k25)


class TestContainmentToP(unittest.TestCase):
    def test_full_containment(self):
        self.assertEqual(mash.containment_to_p(1.0, 15), 0.0)

    def test_zero_containment_undefined(self):
        self.assertIsNone(mash.containment_to_p(0.0, 15))

    def test_monotonic(self):
        p_high_c = mash.containment_to_p(0.9, 15)
        p_low_c = mash.containment_to_p(0.1, 15)
        self.assertLess(p_high_c, p_low_c)


class TestExpectedChanceShared(unittest.TestCase):
    def test_known_value(self):
        # canonical k-mer space at k=2 is 4**2 / 2 = 8
        expected = mash.expected_chance_shared(100, 100, 2)
        self.assertAlmostEqual(expected, 100 * 100 / 8)

    def test_shrinks_as_k_grows(self):
        floor_small_k = mash.expected_chance_shared(1000, 1000, 5)
        floor_large_k = mash.expected_chance_shared(1000, 1000, 20)
        self.assertGreater(floor_small_k, floor_large_k)


class TestResolutionZ(unittest.TestCase):
    def test_far_above_floor_is_resolved(self):
        # k=10 gives a tiny chance floor (~1.9) for n1=n2=1000; a real signal of
        # 500 shared k-mers is nowhere near noise
        self.assertFalse(
            mash.is_resolution_limited(shared=500, n1=1000, n2=1000, k=10)
        )

    def test_shared_at_the_floor_is_limited(self):
        n1 = n2 = 50
        expected = mash.expected_chance_shared(n1, n2, 4)
        shared = round(expected)
        self.assertTrue(mash.is_resolution_limited(shared=shared, n1=n1, n2=n2, k=4))

    def test_zero_inputs_are_limited(self):
        self.assertTrue(mash.is_resolution_limited(shared=0, n1=0, n2=0, k=15))


class TestSummarizePairAcrossK(unittest.TestCase):
    def test_picks_largest_non_limited_k(self):
        # same-ish absolute shared count at two k's, but the chance-collision
        # floor at k=8 is enormous relative to n1*n2 at this scale (guaranteed
        # saturation), while at k=20 it's negligible -- only k=20 should clear
        n1 = n2 = 100_000
        shared = 40_000
        floor_k8 = mash.expected_chance_shared(n1, n2, 8)
        floor_k20 = mash.expected_chance_shared(n1, n2, 20)
        self.assertGreater(floor_k8, shared)
        self.assertLess(floor_k20, shared)

        stats = {8: (shared, n1, n2), 20: (shared, n1, n2)}
        summary = mash.summarize_pair_across_k([8, 20], stats)

        self.assertEqual(summary["chosen_k"], 20)
        self.assertFalse(summary["resolution_limited"])
        self.assertTrue(summary["per_k"][8]["resolution_limited"])
        self.assertFalse(summary["per_k"][20]["resolution_limited"])
        expected_jaccard = shared / (n1 + n2 - shared)
        self.assertAlmostEqual(
            summary["distance"], mash.jaccard_to_distance(expected_jaccard, 20)
        )

    def test_all_k_limited_reports_unresolved_not_a_fake_number(self):
        n1 = n2 = 100
        stats = {k: (0, n1, n2) for k in [11, 15, 19]}
        summary = mash.summarize_pair_across_k([11, 15, 19], stats)
        self.assertIsNone(summary["chosen_k"])
        self.assertIsNone(summary["distance"])
        self.assertTrue(summary["resolution_limited"])

    def test_k_consistency_spread_zero_for_single_valid_k(self):
        n1 = n2 = 100_000
        summary = mash.summarize_pair_across_k([20], {20: (40_000, n1, n2)})
        self.assertEqual(summary["k_consistency_spread"], 0.0)

    def test_distance_clipped_to_unit_interval_in_headline_value(self):
        # a very low but nonzero jaccard at large k can push the raw Mash
        # distance above 1 -- the headline `distance` is clipped for
        # plotting/matrix purposes, but `raw_distance` keeps the true value
        n1 = n2 = 10_000_000
        shared = 10  # a handful of real matches, comfortably above the tiny floor at k=25
        summary = mash.summarize_pair_across_k([25], {25: (shared, n1, n2)})
        self.assertFalse(summary["resolution_limited"])
        self.assertLessEqual(summary["distance"], 1.0)
        self.assertGreaterEqual(summary["raw_distance"], summary["distance"])


if __name__ == "__main__":
    unittest.main()
