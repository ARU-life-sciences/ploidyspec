import itertools
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.homeologs import bh_qvalues, detect_homeolog_pairs, empirical_pvalues


class TestBhQvalues(unittest.TestCase):
    def test_matches_hand_computed_example(self):
        # classic textbook BH example: p = [0.01, 0.02, 0.03, 0.04, 0.05], n=5
        # q_(k) = min_{k'>=k} p_(k')*n/k'
        # raw p*n/k: [0.05, 0.05, 0.05, 0.05, 0.05] -- all equal here by construction
        p = {"a": 0.01, "b": 0.02, "c": 0.03, "d": 0.04, "e": 0.05}
        q = bh_qvalues(p)
        for key in p:
            self.assertAlmostEqual(q[key], 0.05, places=10)

    def test_monotone_non_decreasing_with_p(self):
        p = {"a": 0.001, "b": 0.2, "c": 0.01, "d": 0.5}
        q = bh_qvalues(p)
        ordered = sorted(p, key=lambda k: p[k])
        q_ordered = [q[k] for k in ordered]
        self.assertEqual(q_ordered, sorted(q_ordered))

    def test_single_hypothesis_q_equals_p(self):
        q = bh_qvalues({"a": 0.03})
        self.assertAlmostEqual(q["a"], 0.03)

    def test_never_exceeds_one(self):
        p = {"a": 0.9, "b": 0.95, "c": 0.99}
        q = bh_qvalues(p)
        self.assertTrue(all(v <= 1.0 for v in q.values()))


class TestEmpiricalPvalues(unittest.TestCase):
    def test_lower_distance_gets_lower_pvalue(self):
        pair_dist = {(1, 2): 0.05, (1, 3): 0.11, (2, 3): 0.12, (1, 4): 0.10}
        p, z = empirical_pvalues(pair_dist)
        self.assertLess(p[(1, 2)], p[(1, 3)])
        self.assertLess(p[(1, 2)], min(v for k, v in p.items() if k != (1, 2)))
        # z-score should be finite and move the same direction as p
        self.assertLess(z[(1, 2)], z[(1, 3)])

    def test_not_bounded_below_by_one_over_n(self):
        # regression test: a naive rank/n p-value can never go below 1/n, which
        # (combined with BH's n/rank multiplier) makes q_(1) = 1 always,
        # regardless of the data -- confirmed by hand: rank/n fed into BH gives
        # q_(k) = (k/n)*(n/k) = 1 for every k. A real signal must be able to
        # produce a p-value well under 1/n.
        random.seed(0)
        chrom_nums = list(range(1, 19))
        pair_dist = {}
        true_pair = (1, 5)
        for i, j in itertools.combinations(chrom_nums, 2):
            pair_dist[(i, j)] = (
                random.uniform(0.04, 0.05)
                if (i, j) == true_pair
                else random.uniform(0.09, 0.13)
            )
        n = len(pair_dist)
        p, z = empirical_pvalues(pair_dist)
        self.assertLess(p[true_pair], 1.0 / n)

    def test_z_score_finite_even_when_pvalue_underflows(self):
        # regression test: on real data (daGleHede1, 153 cross-chromosome
        # pairs) the background stdev was tight enough that homeolog pairs
        # sat ~20 stdevs below the mean, underflowing p to exactly 0.0 in
        # float64. z must stay finite and informative in that regime, not
        # also collapse to a fixed/degenerate value.
        random.seed(2)
        chrom_nums = list(range(1, 19))
        pair_dist = {
            (i, j): random.gauss(0.116, 0.0033)
            for i, j in itertools.combinations(chrom_nums, 2)
        }
        pair_dist[(1, 2)] = 0.047  # far below the tight background, like the real case
        p, z = empirical_pvalues(pair_dist)
        import math

        self.assertTrue(math.isfinite(z[(1, 2)]))
        self.assertLess(z[(1, 2)], -5)
        self.assertEqual(p[(1, 2)], 0.0)


class TestDetectHomeologPairs(unittest.TestCase):
    def test_recovers_planted_signal_with_fdr_control(self):
        random.seed(0)
        chrom_nums = list(range(1, 19))
        true_pairs = [
            (1, 5), (2, 3), (4, 10), (6, 9), (7, 8),
            (11, 15), (12, 13), (14, 17), (16, 18),
        ]
        pair_dist = {}
        for i, j in itertools.combinations(chrom_nums, 2):
            pair_dist[(i, j)] = (
                random.uniform(0.04, 0.06)
                if (i, j) in true_pairs
                else random.uniform(0.08, 0.14)
            )
        accepted, unmatched, background = detect_homeolog_pairs(
            pair_dist, chrom_nums, fdr_alpha=0.05
        )
        accepted_pairs = {(i, j) for i, j, *_ in accepted}
        self.assertEqual(accepted_pairs, set(true_pairs))
        self.assertEqual(unmatched, [])

    def test_pure_noise_rejects_almost_everything(self):
        random.seed(1)
        chrom_nums = list(range(1, 19))
        pair_dist = {
            (i, j): random.gauss(0.11, 0.015)
            for i, j in itertools.combinations(chrom_nums, 2)
        }
        accepted, unmatched, background = detect_homeolog_pairs(
            pair_dist, chrom_nums, fdr_alpha=0.05
        )
        # FDR control doesn't guarantee zero false positives, just that their
        # expected proportion is bounded -- with 153 pure-null tests at
        # alpha=0.05 a handful of accidental discoveries is plausible, but
        # accepting most/all chromosomes would mean the test isn't controlling
        # anything.
        self.assertLess(len(accepted), 5)

    def test_each_chromosome_gets_at_most_one_partner(self):
        random.seed(0)
        chrom_nums = list(range(1, 19))
        true_pairs = [
            (1, 5), (2, 3), (4, 10), (6, 9), (7, 8),
            (11, 15), (12, 13), (14, 17), (16, 18),
        ]
        pair_dist = {}
        for i, j in itertools.combinations(chrom_nums, 2):
            pair_dist[(i, j)] = (
                random.uniform(0.04, 0.06)
                if (i, j) in true_pairs
                else random.uniform(0.08, 0.14)
            )
        accepted, _, _ = detect_homeolog_pairs(pair_dist, chrom_nums, fdr_alpha=0.05)
        seen = set()
        for i, j, *_ in accepted:
            self.assertNotIn(i, seen)
            self.assertNotIn(j, seen)
            seen.add(i)
            seen.add(j)


if __name__ == "__main__":
    unittest.main()
