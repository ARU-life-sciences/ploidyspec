import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.te_markers import classify_window, differential_markers


class TestDifferentialMarkers(unittest.TestCase):
    def test_symmetric_counts_produce_no_markers(self):
        counts_a = {"AAAA": 150, "CCCC": 200}
        counts_b = {"AAAA": 150, "CCCC": 200}
        markers, n_a, n_b = differential_markers(counts_a, counts_b, min_ratio=2.0)
        self.assertEqual(markers, {})
        self.assertEqual((n_a, n_b), (0, 0))

    def test_clear_one_sided_kmer_is_a_marker(self):
        counts_a = {"AAAA": 500}
        counts_b = {"AAAA": 100}
        markers, n_a, n_b = differential_markers(counts_a, counts_b, min_ratio=2.0)
        self.assertEqual(markers["AAAA"], (500, 100, "a"))
        self.assertEqual((n_a, n_b), (1, 0))

    def test_clear_other_sided_kmer_is_b_marker(self):
        counts_a = {"TTTT": 80}
        counts_b = {"TTTT": 400}
        markers, n_a, n_b = differential_markers(counts_a, counts_b, min_ratio=2.0)
        self.assertEqual(markers["TTTT"], (80, 400, "b"))
        self.assertEqual((n_a, n_b), (0, 1))

    def test_below_ratio_is_not_a_marker(self):
        # 150 / 100 = 1.5, below min_ratio=2.0 -- shouldn't be classified
        counts_a = {"GGGG": 150}
        counts_b = {"GGGG": 100}
        markers, n_a, n_b = differential_markers(counts_a, counts_b, min_ratio=2.0)
        self.assertNotIn("GGGG", markers)
        self.assertEqual((n_a, n_b), (0, 0))

    def test_exactly_at_ratio_boundary_is_a_marker(self):
        # 200 / 100 == 2.0 exactly -- boundary is inclusive (>=)
        counts_a = {"GGGG": 200}
        counts_b = {"GGGG": 100}
        markers, n_a, n_b = differential_markers(counts_a, counts_b, min_ratio=2.0)
        self.assertEqual(markers["GGGG"], (200, 100, "a"))
        self.assertEqual((n_a, n_b), (1, 0))

    def test_kmer_absent_from_one_side_treated_as_zero(self):
        counts_a = {"ACGT": 300}
        counts_b = {}
        markers, n_a, n_b = differential_markers(counts_a, counts_b, min_ratio=2.0)
        self.assertEqual(markers["ACGT"], (300, 0, "a"))
        self.assertEqual((n_a, n_b), (1, 0))

    def test_low_count_absent_side_does_not_trivially_qualify(self):
        # denominator uses max(other, 1), so a count of 1 vs absent (0) needs
        # 1 >= min_ratio * 1 -- fails for min_ratio=2.0, so no marker
        counts_a = {"CGCG": 1}
        counts_b = {}
        markers, n_a, n_b = differential_markers(counts_a, counts_b, min_ratio=2.0)
        self.assertNotIn("CGCG", markers)
        self.assertEqual((n_a, n_b), (0, 0))


class TestClassifyWindow(unittest.TestCase):
    def test_no_hits_is_none(self):
        self.assertEqual(classify_window(0, 0, min_ratio=2.0), "none")

    def test_clear_a_enrichment(self):
        self.assertEqual(classify_window(434, 3, min_ratio=2.0), "a")

    def test_clear_b_enrichment(self):
        self.assertEqual(classify_window(3, 434, min_ratio=2.0), "b")

    def test_mild_imbalance_is_ambiguous_not_misclassified(self):
        # regression test: a_hits (230) > b_hits (173) but well below
        # min_ratio -- must not be misclassified as "b". This is the exact
        # case that exposed the pool-size-normalization bug: fraction-based
        # classification (hits / each side's total marker-pool size, pools of
        # very different sizes e.g. 3612 vs 1926) flipped this to "b" even
        # though the raw counts clearly favor "a". Raw-count classification
        # must treat it as ambiguous, and must never flip the direction.
        self.assertEqual(classify_window(230, 173, min_ratio=2.0), "ambiguous")

    def test_mild_imbalance_other_direction_is_ambiguous(self):
        self.assertEqual(classify_window(270, 293, min_ratio=2.0), "ambiguous")

    def test_boundary_at_exactly_min_ratio(self):
        self.assertEqual(classify_window(200, 100, min_ratio=2.0), "a")

    def test_one_sided_absence_of_other(self):
        self.assertEqual(classify_window(50, 0, min_ratio=2.0), "a")
        self.assertEqual(classify_window(0, 50, min_ratio=2.0), "b")

    def test_low_count_absence_does_not_trivially_qualify(self):
        # 1 vs 0: denominator uses max(other,1), so 1 >= 2.0*1 is False
        self.assertEqual(classify_window(1, 0, min_ratio=2.0), "ambiguous")


if __name__ == "__main__":
    unittest.main()
