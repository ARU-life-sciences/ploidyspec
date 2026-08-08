import itertools
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.homeologs import (
    bh_qvalues,
    build_ploidy_ancestry_rows,
    detect_homeolog_pairs,
    empirical_pvalues,
    own_chrom_distances,
    rank_candidates,
)


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
        accepted, unmatched, background, p_values, q_values, z_scores = (
            detect_homeolog_pairs(pair_dist, chrom_nums, fdr_alpha=0.05)
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
        accepted, unmatched, background, p_values, q_values, z_scores = (
            detect_homeolog_pairs(pair_dist, chrom_nums, fdr_alpha=0.05)
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
        accepted, _, _, _, _, _ = detect_homeolog_pairs(
            pair_dist, chrom_nums, fdr_alpha=0.05
        )
        seen = set()
        for i, j, *_ in accepted:
            self.assertNotIn(i, seen)
            self.assertNotIn(j, seen)
            seen.add(i)
            seen.add(j)


class TestOwnChromDistances(unittest.TestCase):
    def test_mean_for_multi_copy_chromosome(self):
        groups = {1: ["a", "b", "c"]}
        mat = {
            "a": {"a": 0.0, "b": 0.02, "c": 0.04},
            "b": {"a": 0.02, "b": 0.0, "c": 0.06},
            "c": {"a": 0.04, "b": 0.06, "c": 0.0},
        }
        result = own_chrom_distances(groups, mat)
        # mean of (a,b)=0.02, (a,c)=0.04, (b,c)=0.06 -> 0.04
        self.assertAlmostEqual(result[1], 0.04)

    def test_none_for_singleton_chromosome(self):
        groups = {1: ["a"]}
        mat = {"a": {"a": 0.0}}
        result = own_chrom_distances(groups, mat)
        self.assertIsNone(result[1])


class TestRankCandidates(unittest.TestCase):
    def test_includes_all_pairs_sorted_ascending_with_accepted_flag(self):
        pair_dist = {(1, 2): 0.05, (1, 3): 0.02, (2, 3): 0.09}
        p_values, z_scores = empirical_pvalues(pair_dist)
        q_values = bh_qvalues(p_values)
        accepted_keys = {(1, 3)}  # only the closest pair accepted
        ranked = rank_candidates(pair_dist, p_values, q_values, z_scores, accepted_keys)
        self.assertEqual(len(ranked), 3)
        self.assertEqual([r["distance"] for r in ranked], sorted(pair_dist.values()))
        # accepted flag correct even for the rejected/"bad" candidates -- this
        # is the whole point, a real but FDR-rejected pair must still show up
        accepted_flags = {(r["chrom_a"], r["chrom_b"]): r["accepted"] for r in ranked}
        self.assertTrue(accepted_flags[(1, 3)])
        self.assertFalse(accepted_flags[(1, 2)])
        self.assertFalse(accepted_flags[(2, 3)])


class TestBuildPloidyAncestryRows(unittest.TestCase):
    def test_ratio_computed_when_both_sides_have_own_distance(self):
        groups = {1: ["a", "b"], 2: ["c", "d"]}
        own_dist = {1: 0.01, 2: 0.02}
        pair_dist = {(1, 2): 0.06}
        accepted = [(1, 2, 0.06, 0.001, 0.01, -5.0)]
        q_values = {(1, 2): 0.01}
        z_scores = {(1, 2): -5.0}
        rows = build_ploidy_ancestry_rows(
            groups, own_dist, pair_dist, accepted, q_values, z_scores
        )
        by_chrom = {r["chrom"]: r for r in rows}
        # ratio = 0.06 / mean(0.01, 0.02) = 0.06/0.015 = 4.0
        self.assertAlmostEqual(by_chrom[1]["distance_ratio"], 4.0)
        self.assertAlmostEqual(by_chrom[2]["distance_ratio"], 4.0)
        self.assertEqual(by_chrom[1]["homeolog_partner"], 2)
        self.assertEqual(by_chrom[2]["homeolog_partner"], 1)

    def test_no_partner_leaves_homeolog_fields_blank(self):
        groups = {1: ["a", "b"], 2: ["c", "d"], 3: ["e", "f"]}
        own_dist = {1: 0.01, 2: 0.02, 3: 0.015}
        pair_dist = {(1, 2): 0.06, (1, 3): 0.08, (2, 3): 0.09}
        accepted = [(1, 2, 0.06, 0.001, 0.01, -5.0)]
        q_values = {(1, 2): 0.01, (1, 3): 0.5, (2, 3): 0.6}
        z_scores = {(1, 2): -5.0, (1, 3): -1.0, (2, 3): -0.5}
        rows = build_ploidy_ancestry_rows(
            groups, own_dist, pair_dist, accepted, q_values, z_scores
        )
        by_chrom = {r["chrom"]: r for r in rows}
        self.assertIsNone(by_chrom[3]["homeolog_partner"])
        self.assertIsNone(by_chrom[3]["distance_ratio"])
        self.assertIsNone(by_chrom[3]["homeolog_distance"])

    def test_fewer_than_two_copies_leaves_own_distance_and_ratio_none(self):
        groups = {1: ["a"], 2: ["b", "c"]}
        own_dist = {1: None, 2: 0.02}
        pair_dist = {(1, 2): 0.06}
        accepted = [(1, 2, 0.06, 0.001, 0.01, -5.0)]
        q_values = {(1, 2): 0.01}
        z_scores = {(1, 2): -5.0}
        rows = build_ploidy_ancestry_rows(
            groups, own_dist, pair_dist, accepted, q_values, z_scores
        )
        by_chrom = {r["chrom"]: r for r in rows}
        self.assertIsNone(by_chrom[1]["own_mean_distance"])
        # ratio can't be computed since chrom 1 has no own_mean_distance
        self.assertIsNone(by_chrom[1]["distance_ratio"])
        # but homeolog_distance/z/q are still reported -- only the ratio
        # needs both sides' own distance
        self.assertIsNotNone(by_chrom[1]["homeolog_distance"])


if __name__ == "__main__":
    unittest.main()
