import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.windowed import _pair_matrix


def _row(chrom_a, chrom_b, hap_a, hap_b, win_start, win_end, dist):
    return dict(
        chrom_a=chrom_a,
        chrom_b=chrom_b,
        hap_a=hap_a,
        hap_b=hap_b,
        unit_a=f"{hap_a}_chr{chrom_a:02d}",
        unit_b=f"{hap_b}_chr{chrom_b:02d}",
        win_start=win_start,
        win_end=win_end,
        jaccard_distance=dist,
    )


class TestPairMatrix(unittest.TestCase):
    def test_shape_and_row_order(self):
        rows = [
            _row(1, 1, "HAP1", "HAP2", 1, 100, 0.1),
            _row(1, 1, "HAP1", "HAP3", 1, 100, 0.2),
            _row(1, 1, "HAP1", "HAP2", 101, 200, 0.3),
            _row(1, 1, "HAP1", "HAP3", 101, 200, 0.4),
        ]
        row_labels, windows, mat = _pair_matrix(rows)
        self.assertEqual(row_labels, ["HAP1 vs HAP2", "HAP1 vs HAP3"])
        self.assertEqual(windows, [(1, 100), (101, 200)])
        self.assertEqual(mat.shape, (2, 2))
        self.assertEqual(list(mat[0]), [0.1, 0.3])
        self.assertEqual(list(mat[1]), [0.2, 0.4])

    def test_missing_window_for_a_pair_is_nan(self):
        import math

        rows = [
            _row(1, 1, "HAP1", "HAP2", 1, 100, 0.1),
            _row(1, 1, "HAP1", "HAP2", 101, 200, 0.2),
            _row(1, 1, "HAP1", "HAP3", 1, 100, 0.5),
            # HAP1 vs HAP3 has no row for window (101, 200)
        ]
        row_labels, windows, mat = _pair_matrix(rows)
        hap3_row = row_labels.index("HAP1 vs HAP3")
        self.assertTrue(math.isnan(mat[hap3_row][1]))

    def test_cross_chromosome_rows_key_on_unit_not_haplotype(self):
        # windowed-homeologs mode: same (hap_a, hap_b) label can recur across
        # genuinely different unit pairs, so grouping must use the unit ids.
        rows = [
            _row(1, 5, "HAP1", "HAP2", 1, 100, 0.1),
            _row(2, 5, "HAP1", "HAP2", 1, 100, 0.2),
        ]
        row_labels, windows, mat = _pair_matrix(rows)
        self.assertEqual(len(row_labels), 2)
        self.assertNotEqual(row_labels[0], row_labels[1])


if __name__ == "__main__":
    unittest.main()
