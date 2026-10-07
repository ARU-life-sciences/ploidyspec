import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.rediploidization import (
    chromosome_lineages,
    containment_components,
    copy_state,
    fusion_row,
    lineage_divergence,
    long_copy_outliers,
    pooled_state,
    robust_background,
    split_segments,
    summarize,
    window_split_track,
)

THRESHOLDS = dict(dist_split=1.25, te_split=2.0, window_split=1.25, min_segment_bp=2_000_000)


def unit(uid, hap, chrom, length):
    return dict(unit_id=uid, hap=hap, chrom=chrom, length=str(int(length)))


class TestFusionDetection(unittest.TestCase):
    def test_long_copy_outlier_needs_clear_length_excess(self):
        units = [
            unit("HAP1_chr19", "HAP1", 19, 38e6),
            unit("HAP2_chr19", "HAP2", 19, 40e6),
            unit("HAP3_chr19", "HAP3", 19, 69e6),
            unit("HAP1_chr05", "HAP1", 5, 50e6),
            unit("HAP2_chr05", "HAP2", 5, 52e6),
        ]
        out = [u["unit_id"] for u in long_copy_outliers(units, 1.4)]
        self.assertEqual(out, ["HAP3_chr19"])

    def test_containment_components_picks_only_far_outliers(self):
        background = robust_background([0.06, 0.063, 0.065, 0.07, 0.058, 0.062])
        hits = containment_components({22: 0.45, 5: 0.08, 17: 0.07}, background, 10)
        self.assertEqual([h[0] for h in hits], [22])

    def test_near_zero_background_does_not_promote_shared_repeats(self):
        # random-sequence background: MAD at its floor, shared-TE containment ~0.02
        background = robust_background([0.0, 0.001, 0.0005, 0.0, 0.001])
        hits = containment_components({4: 0.58, 1: 0.029, 7: 0.021}, background, 10)
        self.assertEqual([h[0] for h in hits], [4])

    def test_fusion_requires_partner_missing_from_haplotype(self):
        hits = [(22, 0.45, 40.0)]
        absorbed = fusion_row("HAP3_chr19", "HAP3", 69e6, "placed_long_copy", 19, hits,
                              {19: {"HAP1", "HAP3"}, 22: {"HAP1"}}, 23)
        self.assertEqual(absorbed["status"], "fusion")
        self.assertEqual(absorbed["components"], "chr19+chr22")
        duplicated = fusion_row("HAP3_chr19", "HAP3", 69e6, "placed_long_copy", 19, hits,
                                {19: {"HAP1", "HAP3"}, 22: {"HAP1", "HAP3"}}, 23)
        self.assertEqual(duplicated["status"], "candidate_partner_present")

    def test_orphan_with_two_components_is_a_fusion(self):
        hits = [(4, 0.5, 30.0), (15, 0.48, 29.0)]
        row = fusion_row("Sy_Chr04_15_M2", "HAP2", 59e6, "unplaced_scaffold", None, hits,
                         {4: {"HAP1"}, 15: {"HAP1"}}, 23)
        self.assertEqual((row["components"], row["status"]), ("chr04+chr15", "fusion"))

    def test_lineage_divergence_orders_by_containment(self):
        self.assertAlmostEqual(lineage_divergence(1.0, 23), 0.0)
        self.assertGreater(lineage_divergence(0.29, 23), lineage_divergence(0.71, 23))

    def test_orphan_with_one_component_is_not_a_fusion(self):
        row = fusion_row("scaffold_7", "HAP1", 5e6, "unplaced_scaffold", None,
                         [(4, 0.3, 15.0)], {4: {"HAP2"}}, 23)
        self.assertEqual(row["status"], "no_fusion_signal")


class TestSplitSegments(unittest.TestCase):
    def track(self, pattern, win=250_000):
        return [(i * win + 1, (i + 1) * win, 2.0 if c == "#" else 1.0)
                for i, c in enumerate(pattern)]

    def test_bridges_single_window_gap_and_drops_short_runs(self):
        segs = split_segments(self.track("####.#####......##"), 1.25, 2_000_000)
        self.assertEqual(segs, [(1, 2_500_000)])

    def test_no_split_windows_gives_no_segments(self):
        self.assertEqual(split_segments(self.track("........"), 1.25, 1), [])

    def test_window_track_uses_bipartition(self):
        windows = {1: (250_000, {("A", "B"): 0.1, ("C", "D"): 0.1, ("A", "C"): 0.5,
                                 ("A", "D"): 0.5, ("B", "C"): 0.5, ("B", "D"): 0.5})}
        (_, _, ratio), = window_split_track(windows, ["A", "B"])
        self.assertAlmostEqual(ratio, 5.0)


class TestCopyState(unittest.TestCase):
    def test_fusion_like_whole_chromosome_split_is_resolved(self):
        self.assertEqual(copy_state(4, True, 8.1, 5.6, "whole", THRESHOLDS), "resolved_lineages")

    def test_regional_split_with_te_support_is_partial(self):
        self.assertEqual(copy_state(4, True, 1.21, 3.2, "regional", THRESHOLDS), "partially_resolved")

    def test_single_line_of_evidence_is_candidate(self):
        self.assertEqual(copy_state(4, True, 1.08, 1.17, "regional", THRESHOLDS), "candidate")

    def test_no_evidence_is_tetrasomic_like(self):
        self.assertEqual(copy_state(4, True, 1.05, 0.9, "none", THRESHOLDS), "tetrasomic_like")

    def test_unbalanced_split_is_one_divergent_copy(self):
        self.assertEqual(copy_state(4, False, 3.45, 1.7, "none", THRESHOLDS), "one_divergent_copy")

    def test_two_copies_not_assessable(self):
        self.assertEqual(copy_state(2, False, None, None, "", THRESHOLDS), "not_assessable")


class TestPooledState(unittest.TestCase):
    def test_separate_lineages_across_homeologs(self):
        self.assertEqual(pooled_state(True, 3.5, THRESHOLDS), "resolved_lineages")

    def test_interchangeable_copies_across_homeologs(self):
        self.assertEqual(pooled_state(True, 1.05, THRESHOLDS), "tetrasomic_like")

    def test_intermediate_is_candidate(self):
        self.assertEqual(pooled_state(True, 1.6, THRESHOLDS), "candidate")

    def test_one_copy_apart(self):
        self.assertEqual(pooled_state(False, 1.45, THRESHOLDS), "one_divergent_copy")


class TestChromosomeLineages(unittest.TestCase):
    def test_two_lineages_recovered_from_distances(self):
        units = [unit(u, u[:4], 19, 40e6) for u in ("HAP1_chr19", "HAP2_chr19", "HAP3_chr19", "HAP4_chr19")]
        close = {("HAP1_chr19", "HAP2_chr19"), ("HAP3_chr19", "HAP4_chr19")}
        dist = {}
        ids = [u["unit_id"] for u in units]
        for a in ids:
            for b in ids:
                if a != b:
                    dist[(a, b)] = 0.006 if (a, b) in close or (b, a) in close else 0.054
        row = chromosome_lineages(units, dist, {}, {}, THRESHOLDS)[19]
        self.assertTrue(row["balanced"])
        self.assertAlmostEqual(row["dist_split"], 9.0)
        self.assertEqual({row["group_a"], row["group_b"]},
                         {"HAP1_chr19,HAP2_chr19", "HAP3_chr19,HAP4_chr19"})


class TestSummary(unittest.TestCase):
    def test_counts_states_and_fusions(self):
        rows = [
            dict(copy_state="tetrasomic_like", ancient_state="paired", distance_ratio="7.59"),
            dict(copy_state="fusion_lineages", ancient_state="fusion_partner", distance_ratio="3.27"),
        ]
        summary = {r["metric"]: r["value"] for r in summarize(rows, [dict(status="fusion", components="chr19+chr22"),
                                                         dict(status="fusion", components="chr19+chr22")])}
        self.assertEqual(summary["copy_state:fusion_lineages"], 1)
        self.assertEqual(summary["n_distinct_fusions"], 1)
        self.assertEqual(summary["n_fused_scaffolds"], 2)
        self.assertEqual(summary["ancient_paired_chromosomes"], 1)


if __name__ == "__main__":
    unittest.main()
