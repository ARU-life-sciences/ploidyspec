import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.subgenome_report import (
    summarize_windowed_assignment,
    te_marker_fraction,
)


class TestTeMarkerFraction(unittest.TestCase):
    def test_zero_when_no_high_copy_kmers(self):
        self.assertEqual(te_marker_fraction(0, 0, 0, 0), 0.0)

    def test_bounded_between_zero_and_one(self):
        # markers are a subset of high-copy kmers by construction
        frac = te_marker_fraction(3612, 1926, 8666, 6715)
        self.assertGreaterEqual(frac, 0.0)
        self.assertLessEqual(frac, 1.0)

    def test_all_shared_gives_zero(self):
        # no differential markers at all -- identical repeat content
        self.assertEqual(te_marker_fraction(0, 0, 5000, 5000), 0.0)

    def test_all_differential_gives_one(self):
        # every high-copy kmer on both sides is a marker -- completely
        # non-overlapping repeat content
        self.assertEqual(te_marker_fraction(5000, 5000, 5000, 5000), 1.0)

    def test_wheat_scale_gives_much_higher_fraction_than_glechoma_scale(self):
        # sanity check using this session's real orders of magnitude: wheat's
        # A/B subgenomes should score far higher than daGleHede1's haplotype
        # copies, since that's the entire point of the index
        glechoma_like = te_marker_fraction(3612, 1926, 8666, 6715)
        wheat_like = te_marker_fraction(331916, 497489, 600000, 750000)
        self.assertGreater(wheat_like, glechoma_like)


class TestSummarizeWindowedAssignment(unittest.TestCase):
    def _write_windowed_tsv(self, path, rows):
        import csv

        with open(path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(
                ["unit", "win_start", "win_end", "a_hits", "b_hits", "a_frac", "b_frac", "assigned"]
            )
            for r in rows:
                w.writerow(r)

    def test_self_matching_windows_not_anomalous(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "te_markers_windowed_UA_chr01xUB_chr01.tsv")
            self._write_windowed_tsv(
                path,
                [
                    ["UA_chr01", 1, 100, 5, 0, "0.5", "0.0", "a"],
                    ["UB_chr01", 1, 100, 0, 5, "0.0", "0.5", "b"],
                ],
            )
            rows, anomalous = summarize_windowed_assignment(path, "UA_chr01", "UB_chr01")
            self.assertEqual(anomalous, [])
            by_unit = {r["unit"]: r for r in rows}
            self.assertEqual(by_unit["UA_chr01"]["n_self"], 1)
            self.assertEqual(by_unit["UA_chr01"]["n_other"], 0)

    def test_cross_matching_window_is_anomalous(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "te_markers_windowed_UA_chr01xUB_chr01.tsv")
            self._write_windowed_tsv(
                path,
                [
                    ["UA_chr01", 1, 100, 0, 5, "0.0", "0.5", "b"],  # UA's window matches B -- anomalous
                    ["UB_chr01", 1, 100, 0, 5, "0.0", "0.5", "b"],  # UB's window matches B -- self, fine
                ],
            )
            rows, anomalous = summarize_windowed_assignment(path, "UA_chr01", "UB_chr01")
            self.assertEqual(len(anomalous), 1)
            self.assertEqual(anomalous[0]["unit"], "UA_chr01")
            self.assertEqual(anomalous[0]["matched_side"], "b")

    def test_ambiguous_and_none_not_counted_as_anomalous(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "te_markers_windowed_UA_chr01xUB_chr01.tsv")
            self._write_windowed_tsv(
                path,
                [
                    ["UA_chr01", 1, 100, 2, 2, "0.2", "0.2", "ambiguous"],
                    ["UA_chr01", 101, 200, 0, 0, "0.0", "0.0", "none"],
                ],
            )
            rows, anomalous = summarize_windowed_assignment(path, "UA_chr01", "UB_chr01")
            self.assertEqual(anomalous, [])
            by_unit = {r["unit"]: r for r in rows}
            self.assertEqual(by_unit["UA_chr01"]["n_ambiguous"], 1)
            self.assertEqual(by_unit["UA_chr01"]["n_none"], 1)
            self.assertEqual(by_unit["UA_chr01"]["total_windows"], 2)


if __name__ == "__main__":
    unittest.main()
