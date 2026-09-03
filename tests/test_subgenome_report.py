import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.common import matrix_dir
from ploidyspec.subgenome_report import (
    bipartition_by_distance,
    compute_lineage_te_fractions,
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


class TestBipartitionByDistance(unittest.TestCase):
    def test_fewer_than_three_units_returns_none(self):
        dist = {("A", "B"): 0.1, ("B", "A"): 0.1}
        self.assertIsNone(bipartition_by_distance(["A", "B"], dist))

    def test_recovers_a_clean_two_two_split(self):
        # A,B close to each other; C,D close to each other; both pairs far apart --
        # the classic fused/unfused (or subgenome) signature.
        units = ["A", "B", "C", "D"]
        pairs = {
            ("A", "B"): 0.01,
            ("C", "D"): 0.01,
            ("A", "C"): 0.5,
            ("A", "D"): 0.5,
            ("B", "C"): 0.5,
            ("B", "D"): 0.5,
        }
        dist = {}
        for (u, v), d in pairs.items():
            dist[(u, v)] = d
            dist[(v, u)] = d
        group_a, group_b = bipartition_by_distance(units, dist)
        self.assertEqual({frozenset(group_a), frozenset(group_b)}, {frozenset(["A", "B"]), frozenset(["C", "D"])})

    def test_missing_distances_returns_none(self):
        self.assertIsNone(bipartition_by_distance(["A", "B", "C"], {}))


class TestComputeLineageTeFractions(unittest.TestCase):
    def _write_matrix_csv(self, outdir, units, dist):
        os.makedirs(matrix_dir(outdir), exist_ok=True)
        path = os.path.join(matrix_dir(outdir), "whole_chrom_distance_matrix.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([""] + units)
            for u in units:
                w.writerow([u] + [f"{dist.get((u, v), 0.0):.6f}" for v in units])

    def test_within_vs_cross_lineage_split_recovered(self):
        # Mirrors the real SchCurv1 chr19 case: HAP1/HAP2 (unfused) close to
        # each other, HAP3/HAP4 (fused) close to each other, both groups far
        # apart -- and te_marker_fraction should show the same split.
        units = ["HAP1_chr19", "HAP2_chr19", "HAP3_chr19", "HAP4_chr19"]
        close_pairs = {("HAP1_chr19", "HAP2_chr19"), ("HAP3_chr19", "HAP4_chr19")}
        dist = {}
        for i, u in enumerate(units):
            for v in units[i + 1 :]:
                d = 0.02 if (u, v) in close_pairs else 0.3
                dist[(u, v)] = d
                dist[(v, u)] = d

        with tempfile.TemporaryDirectory() as td:
            self._write_matrix_csv(td, units, dist)
            index_rows = [
                dict(unit_a="HAP1_chr19", unit_b="HAP2_chr19", chrom="chr19", te_marker_fraction=0.14),
                dict(unit_a="HAP3_chr19", unit_b="HAP4_chr19", chrom="chr19", te_marker_fraction=0.07),
                dict(unit_a="HAP1_chr19", unit_b="HAP3_chr19", chrom="chr19", te_marker_fraction=0.59),
                dict(unit_a="HAP1_chr19", unit_b="HAP4_chr19", chrom="chr19", te_marker_fraction=0.57),
                dict(unit_a="HAP2_chr19", unit_b="HAP3_chr19", chrom="chr19", te_marker_fraction=0.60),
                dict(unit_a="HAP2_chr19", unit_b="HAP4_chr19", chrom="chr19", te_marker_fraction=0.61),
            ]
            rows = compute_lineage_te_fractions(td, index_rows)
            self.assertEqual(len(rows), 1)
            r = rows[0]
            self.assertEqual(r["n_within"], 2)
            self.assertEqual(r["n_cross"], 4)
            self.assertAlmostEqual(r["within_te_marker_fraction"], (0.14 + 0.07) / 2, places=6)
            self.assertAlmostEqual(r["cross_te_marker_fraction"], (0.59 + 0.57 + 0.60 + 0.61) / 4, places=6)
            self.assertGreater(r["split_ratio"], 4.0)

    def test_two_copy_chromosome_skipped(self):
        units = ["HAP1_chr01", "HAP2_chr01"]
        with tempfile.TemporaryDirectory() as td:
            self._write_matrix_csv(td, units, {("HAP1_chr01", "HAP2_chr01"): 0.1, ("HAP2_chr01", "HAP1_chr01"): 0.1})
            index_rows = [
                dict(unit_a="HAP1_chr01", unit_b="HAP2_chr01", chrom="chr01", te_marker_fraction=0.2),
            ]
            rows = compute_lineage_te_fractions(td, index_rows)
            self.assertEqual(rows, [])

    def test_no_real_structure_gives_split_ratio_near_one(self):
        # all four copies roughly equidistant -- no genuine lineage split
        units = ["HAP1_chr01", "HAP2_chr01", "HAP3_chr01", "HAP4_chr01"]
        dist = {}
        for i, u in enumerate(units):
            for v in units[i + 1 :]:
                dist[(u, v)] = 0.1
                dist[(v, u)] = 0.1
        with tempfile.TemporaryDirectory() as td:
            self._write_matrix_csv(td, units, dist)
            index_rows = [
                dict(unit_a="HAP1_chr01", unit_b="HAP2_chr01", chrom="chr01", te_marker_fraction=0.15),
                dict(unit_a="HAP1_chr01", unit_b="HAP3_chr01", chrom="chr01", te_marker_fraction=0.16),
                dict(unit_a="HAP1_chr01", unit_b="HAP4_chr01", chrom="chr01", te_marker_fraction=0.14),
                dict(unit_a="HAP2_chr01", unit_b="HAP3_chr01", chrom="chr01", te_marker_fraction=0.15),
                dict(unit_a="HAP2_chr01", unit_b="HAP4_chr01", chrom="chr01", te_marker_fraction=0.16),
                dict(unit_a="HAP3_chr01", unit_b="HAP4_chr01", chrom="chr01", te_marker_fraction=0.15),
            ]
            rows = compute_lineage_te_fractions(td, index_rows)
            self.assertEqual(len(rows), 1)
            self.assertAlmostEqual(rows[0]["split_ratio"], 1.0, delta=0.1)


if __name__ == "__main__":
    unittest.main()
