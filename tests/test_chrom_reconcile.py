import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.chrom_reconcile import detect_relabeling, reconcile_chrom_labels


def unit(unit_id, hap, chrom, source):
    return dict(
        unit_id=unit_id,
        hap=hap,
        chrom=chrom,
        seq_id=unit_id,
        length="1000000",
        source=source,
        desc="",
    )


def symmetric_matrix(n, pairs, default=0.10):
    """pairs: {(i, j): distance}. Builds a full n x n symmetric list-of-lists,
    filling anything unspecified with `default` (and 0.0 on the diagonal)."""
    m = [[default] * n for _ in range(n)]
    for i in range(n):
        m[i][i] = 0.0
    for (i, j), d in pairs.items():
        m[i][j] = d
        m[j][i] = d
    return m


class TestDetectRelabeling(unittest.TestCase):
    def test_two_way_swap_detected_and_corrected(self):
        # mirrors drLytSali1 chr01<->chr02: HAP1's file numbered these two
        # chromosomes in the opposite order from HAP2's file. fileB is listed
        # first (as it would be first in the manifest) so it's picked as the
        # reference -- only fileA's (HAP1's) labels get corrected.
        units = [
            unit("HAP2_chr01", "HAP2", "1", "fileB"),  # reference, correct
            unit("HAP1_chr01", "HAP1", "1", "fileA"),  # truly chr02
            unit("HAP2_chr02", "HAP2", "2", "fileB"),  # reference, correct
            unit("HAP1_chr02", "HAP1", "2", "fileA"),  # truly chr01
        ]
        # indices: 0=HAP2_chr01, 1=HAP1_chr01, 2=HAP2_chr02, 3=HAP1_chr02
        distance = symmetric_matrix(
            4,
            {
                (1, 0): 0.10,  # HAP1_chr01 x declared partner HAP2_chr01 -- far
                (1, 2): 0.02,  # HAP1_chr01 x true partner HAP2_chr02 -- close
                (3, 2): 0.10,  # HAP1_chr02 x declared partner HAP2_chr02 -- far
                (3, 0): 0.02,  # HAP1_chr02 x true partner HAP2_chr01 -- close
            },
        )
        corrections, ambiguous = detect_relabeling(units, distance)
        self.assertEqual(len(corrections), 2)
        self.assertEqual(ambiguous, [])
        by_unit = {c["unit_id"]: c for c in corrections}
        self.assertEqual(by_unit["HAP1_chr01"]["new_chrom"], "2")
        self.assertEqual(by_unit["HAP1_chr02"]["new_chrom"], "1")

    def test_majority_file_chosen_as_reference(self):
        # mirrors drLytSali1's real structure: one "primary" file contributes 1
        # unit/chrom (HAP1), one "alt contigs" file resolves into 3 haplotigs/chrom
        # (HAP2/HAP3/HAP4) -- the 3-unit file should win reference by unit count
        # regardless of listing order, and HAP1 is what gets corrected.
        chroms = [4, 5, 6, 7, 8, 9]
        true_target = {4: 9, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8}
        units = []
        for c in chroms:
            # HAP1 (minority file) listed first, HAP2/3/4 (majority file) after --
            # majority-by-count should still win even though it's not first-listed.
            units.append(unit(f"HAP1_chr{c:02d}", "HAP1", str(c), "fileA"))
            units.append(unit(f"HAP2_chr{c:02d}", "HAP2", str(c), "fileB"))
            units.append(unit(f"HAP3_chr{c:02d}", "HAP3", str(c), "fileB"))
            units.append(unit(f"HAP4_chr{c:02d}", "HAP4", str(c), "fileB"))
        index_of = {u["unit_id"]: i for i, u in enumerate(units)}

        pairs = {}
        # fileB's own trio is already mutually self-consistent for each chrom
        for c in chroms:
            a, b, cc = (
                index_of[f"HAP2_chr{c:02d}"],
                index_of[f"HAP3_chr{c:02d}"],
                index_of[f"HAP4_chr{c:02d}"],
            )
            pairs[(a, b)] = 0.02
            pairs[(a, cc)] = 0.02
            pairs[(b, cc)] = 0.02
        # HAP1_chr{c} truly matches true_target[c]'s HAP2/3/4 trio
        for c in chroms:
            i = index_of[f"HAP1_chr{c:02d}"]
            tc = true_target[c]
            for h in ("HAP2", "HAP3", "HAP4"):
                pairs[(i, index_of[f"{h}_chr{tc:02d}"])] = 0.02

        distance = symmetric_matrix(len(units), pairs)
        corrections, ambiguous = detect_relabeling(units, distance)
        self.assertEqual(len(corrections), 6)
        self.assertEqual(ambiguous, [])
        by_unit = {c["unit_id"]: c["new_chrom"] for c in corrections}
        for c in chroms:
            self.assertEqual(by_unit[f"HAP1_chr{c:02d}"], str(true_target[c]))
            # the reference file's own units are never touched
            self.assertNotIn(f"HAP2_chr{c:02d}", by_unit)
            self.assertNotIn(f"HAP3_chr{c:02d}", by_unit)
            self.assertNotIn(f"HAP4_chr{c:02d}", by_unit)

    def test_mid_range_ratio_is_ambiguous_not_corrected(self):
        # ratio ~2x: between the ambiguous floor (1.5) and the correction
        # threshold (3.0, default) -- flagged for review, not auto-applied.
        units = [
            unit("HAP2_chr04", "HAP2", "4", "fileB"),
            unit("HAP1_chr04", "HAP1", "4", "fileA"),
            unit("HAP2_chr05", "HAP2", "5", "fileB"),
            unit("HAP1_chr05", "HAP1", "5", "fileA"),
        ]
        distance = symmetric_matrix(
            4,
            {
                (1, 0): 0.020,  # declared partner
                (1, 2): 0.010,  # alternative -- 2x closer, not decisive enough
                (3, 2): 0.020,
                (3, 0): 0.010,
            },
        )
        corrections, ambiguous = detect_relabeling(units, distance)
        self.assertEqual(corrections, [])
        self.assertEqual(len(ambiguous), 2)

    def test_near_tie_produces_no_candidate_at_all(self):
        # mirrors the real lpElePalu1 case: ratio ~1.03x, genuinely not a
        # mismatch -- below even the ambiguous floor, no evidence of a problem.
        units = [
            unit("HAP2_chr04", "HAP2", "4", "fileB"),
            unit("HAP1_chr04", "HAP1", "4", "fileA"),
            unit("HAP2_chr05", "HAP2", "5", "fileB"),
            unit("HAP1_chr05", "HAP1", "5", "fileA"),
        ]
        distance = symmetric_matrix(
            4,
            {
                (1, 0): 0.0176,
                (1, 2): 0.0170,
                (3, 2): 0.0176,
                (3, 0): 0.0170,
            },
        )
        corrections, ambiguous = detect_relabeling(units, distance)
        self.assertEqual(corrections, [])
        self.assertEqual(ambiguous, [])

    def test_fully_consistent_panel_is_idempotent(self):
        units = [
            unit("HAP2_chr01", "HAP2", "1", "fileB"),
            unit("HAP1_chr01", "HAP1", "1", "fileA"),
            unit("HAP2_chr02", "HAP2", "2", "fileB"),
            unit("HAP1_chr02", "HAP1", "2", "fileA"),
            unit("HAP2_chr03", "HAP2", "3", "fileB"),
            unit("HAP1_chr03", "HAP1", "3", "fileA"),
        ]
        distance = symmetric_matrix(
            6,
            {
                (0, 1): 0.02,
                (2, 3): 0.02,
                (4, 5): 0.02,
            },
            default=0.10,
        )
        corrections, ambiguous = detect_relabeling(units, distance)
        self.assertEqual(corrections, [])
        self.assertEqual(ambiguous, [])
        # re-running detection on the same (already-clean) state changes nothing
        corrections2, ambiguous2 = detect_relabeling(units, distance)
        self.assertEqual(corrections2, [])
        self.assertEqual(ambiguous2, [])

    def test_correction_then_rerun_finds_nothing_left(self):
        units = [
            unit("HAP2_chr01", "HAP2", "1", "fileB"),
            unit("HAP1_chr01", "HAP1", "1", "fileA"),
            unit("HAP2_chr02", "HAP2", "2", "fileB"),
            unit("HAP1_chr02", "HAP1", "2", "fileA"),
        ]
        distance = symmetric_matrix(
            4,
            {
                (1, 0): 0.10,
                (1, 2): 0.02,
                (3, 2): 0.10,
                (3, 0): 0.02,
            },
        )
        corrections, _ = detect_relabeling(units, distance)
        self.assertEqual(len(corrections), 2)
        for c in corrections:
            units[c["index"]]["chrom"] = c["new_chrom"]
        corrections2, ambiguous2 = detect_relabeling(units, distance)
        self.assertEqual(corrections2, [])
        self.assertEqual(ambiguous2, [])

    def test_two_unit_species_does_not_crash(self):
        units = [
            unit("HAP2_chr01", "HAP2", "1", "fileB"),
            unit("HAP1_chr01", "HAP1", "1", "fileA"),
        ]
        distance = symmetric_matrix(2, {(0, 1): 0.01})
        corrections, ambiguous = detect_relabeling(units, distance)
        self.assertEqual(corrections, [])
        self.assertEqual(ambiguous, [])

    def test_move_into_occupied_non_moving_slot_is_rejected_not_applied(self):
        # regression: real lpElePalu1 run crashed because HAP2_chr12 was
        # "corrected" to chr11 while the genuinely-correct HAP2_chr11 stayed
        # put -- two units both ending up named HAP2_chr11 (a KeyError two
        # steps downstream in write_ploidy_and_homology_reports). The target
        # slot being occupied by a unit that isn't itself relocating must
        # block the move rather than create a duplicate unit_id.
        units = [
            unit("HAP1_chr10", "HAP1", "10", "fileA"),  # reference
            unit("HAP1_chr11", "HAP1", "11", "fileA"),  # reference
            unit("HAP1_chr12", "HAP1", "12", "fileA"),  # reference
            unit("HAP2_chr10", "HAP2", "10", "fileB"),  # correctly matches chr10
            unit("HAP2_chr11", "HAP2", "11", "fileB"),  # correctly matches chr11 -- stays
            unit("HAP2_chr12", "HAP2", "12", "fileB"),  # mislabeled, true match is chr11
        ]
        distance = symmetric_matrix(
            6,
            {
                (3, 0): 0.02,  # HAP2_chr10 x HAP1_chr10 -- correct
                (4, 1): 0.02,  # HAP2_chr11 x HAP1_chr11 -- correct
                (5, 2): 0.10,  # HAP2_chr12 x declared partner HAP1_chr12 -- far
                (5, 1): 0.02,  # HAP2_chr12 x HAP1_chr11 -- closer, but chr11 is taken
            },
        )
        corrections, ambiguous = detect_relabeling(units, distance)
        self.assertEqual(corrections, [])
        self.assertEqual(len(ambiguous), 1)
        self.assertEqual(ambiguous[0]["unit_id"], "HAP2_chr12")

        # applying this must never produce two units with the same unit_id
        with tempfile.TemporaryDirectory() as td:
            seq_tsv = os.path.join(td, "sequences.tsv")
            with open(seq_tsv, "w", newline="") as f:
                w = csv.DictWriter(
                    f,
                    fieldnames=[
                        "unit_id",
                        "hap",
                        "chrom",
                        "seq_id",
                        "length",
                        "source",
                        "desc",
                    ],
                    delimiter="\t",
                )
                w.writeheader()
                for u in units:
                    w.writerow(u)
            reconcile_chrom_labels(td, seq_tsv, units, distance)
            unit_ids = [u["unit_id"] for u in units]
            self.assertEqual(len(unit_ids), len(set(unit_ids)))


class TestReconcileChromLabelsIO(unittest.TestCase):
    def test_applies_correction_renames_files_and_logs(self):
        units = [
            unit("HAP2_chr01", "HAP2", "1", "fileB"),
            unit("HAP1_chr01", "HAP1", "1", "fileA"),
            unit("HAP2_chr02", "HAP2", "2", "fileB"),
            unit("HAP1_chr02", "HAP1", "2", "fileA"),
        ]
        distance = symmetric_matrix(
            4,
            {
                (1, 0): 0.10,
                (1, 2): 0.02,
                (3, 2): 0.10,
                (3, 0): 0.02,
            },
        )
        with tempfile.TemporaryDirectory() as td:
            os.makedirs(os.path.join(td, "chroms"))
            os.makedirs(os.path.join(td, "ktabs_k13"))
            for u in units:
                open(os.path.join(td, "chroms", u["unit_id"] + ".fa"), "w").close()
                open(
                    os.path.join(td, "ktabs_k13", u["unit_id"] + ".ktab"), "w"
                ).close()
                open(
                    os.path.join(td, "ktabs_k13", u["unit_id"] + ".hist"), "w"
                ).close()
                open(
                    os.path.join(td, "ktabs_k13", "." + u["unit_id"] + ".ktab.1"),
                    "w",
                ).close()

            seq_tsv = os.path.join(td, "sequences.tsv")
            with open(seq_tsv, "w", newline="") as f:
                w = csv.DictWriter(
                    f,
                    fieldnames=[
                        "unit_id",
                        "hap",
                        "chrom",
                        "seq_id",
                        "length",
                        "source",
                        "desc",
                    ],
                    delimiter="\t",
                )
                w.writeheader()
                for u in units:
                    w.writerow(u)

            corrections, ambiguous = reconcile_chrom_labels(
                td, seq_tsv, units, distance
            )
            self.assertEqual(len(corrections), 2)

            # unit_ids swapped in place (HAP1's file is the one that moved)
            self.assertEqual(units[1]["unit_id"], "HAP1_chr02")
            self.assertEqual(units[3]["unit_id"], "HAP1_chr01")

            # renamed files exist on disk under the new unit_id
            self.assertTrue(
                os.path.exists(os.path.join(td, "chroms", "HAP1_chr02.fa"))
            )
            self.assertTrue(
                os.path.exists(os.path.join(td, "ktabs_k13", "HAP1_chr02.ktab"))
            )
            self.assertTrue(
                os.path.exists(os.path.join(td, "ktabs_k13", "HAP1_chr02.hist"))
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(td, "ktabs_k13", ".HAP1_chr02.ktab.1")
                )
            )
            # it's a swap: HAP1_chr01.fa also still exists, just now belonging
            # to the unit that used to be HAP1_chr02 -- no file was lost, and no
            # leftover .reconcile_tmp intermediate remains
            self.assertTrue(
                os.path.exists(os.path.join(td, "chroms", "HAP1_chr01.fa"))
            )
            self.assertEqual(
                sorted(os.listdir(os.path.join(td, "chroms"))),
                ["HAP1_chr01.fa", "HAP1_chr02.fa", "HAP2_chr01.fa", "HAP2_chr02.fa"],
            )
            # the reference file's own files are untouched
            self.assertTrue(
                os.path.exists(os.path.join(td, "chroms", "HAP2_chr01.fa"))
            )

            # sequences.tsv rewritten with corrected chrom/unit_id
            with open(seq_tsv) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            unit_ids = {r["unit_id"] for r in rows}
            self.assertIn("HAP1_chr02", unit_ids)
            self.assertIn("HAP1_chr01", unit_ids)

            # corrections log written and auditable
            log_path = os.path.join(td, "matrix", "chrom_label_corrections.tsv")
            self.assertTrue(os.path.exists(log_path))
            with open(log_path) as f:
                log_rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(len(log_rows), 2)
            self.assertTrue(all(r["status"] == "corrected" for r in log_rows))

    def test_no_corrections_still_writes_empty_log(self):
        units = [
            unit("HAP2_chr01", "HAP2", "1", "fileB"),
            unit("HAP1_chr01", "HAP1", "1", "fileA"),
        ]
        distance = symmetric_matrix(2, {(0, 1): 0.01})
        with tempfile.TemporaryDirectory() as td:
            seq_tsv = os.path.join(td, "sequences.tsv")
            with open(seq_tsv, "w", newline="") as f:
                w = csv.DictWriter(
                    f,
                    fieldnames=[
                        "unit_id",
                        "hap",
                        "chrom",
                        "seq_id",
                        "length",
                        "source",
                        "desc",
                    ],
                    delimiter="\t",
                )
                w.writeheader()
                for u in units:
                    w.writerow(u)

            corrections, ambiguous = reconcile_chrom_labels(
                td, seq_tsv, units, distance
            )
            self.assertEqual(corrections, [])
            self.assertEqual(ambiguous, [])
            log_path = os.path.join(td, "matrix", "chrom_label_corrections.tsv")
            self.assertTrue(os.path.exists(log_path))
            with open(log_path) as f:
                log_rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(log_rows, [])


if __name__ == "__main__":
    unittest.main()
