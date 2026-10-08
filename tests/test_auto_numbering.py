import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.chrom_reconcile import detect_relabeling
from ploidyspec.manifest import AUTO_OFFSET, auto_number, chromosome_scale


def seq(hap, seq_id, length, source):
    return dict(hap=hap, seq_id=seq_id, length=length, source=source, desc="")


class TestChromosomeScale(unittest.TestCase):
    def test_drops_small_unplaced_contigs(self):
        # masu salmon HAP2: 33 chromosomes of 18-127 Mb and unplaced 1-5 Mb contigs
        lengths = {f"group{i}": int(127e6 - i * 3.3e6) for i in range(33)}
        lengths.update({"h2tg190": 5_000_000, "h2tg419": 3_700_000, "h2tg136": 1_400_000})
        kept = chromosome_scale(lengths)
        self.assertEqual(len(kept), 33)
        self.assertNotIn("h2tg190", kept)

    def test_all_similar_lengths_are_kept(self):
        self.assertEqual(len(chromosome_scale({f"c{i}": 30_000_000 + i for i in range(10)})), 10)


class TestAutoNumber(unittest.TestCase):
    def candidates(self):
        a = [seq("HAP1", f"a{i}", L, "h1.fa") for i, L in enumerate([50e6, 40e6, 30e6])]
        b = [seq("HAP2", f"b{i}", L, "h2.fa") for i, L in enumerate([31e6, 52e6, 39e6, 2e5])]
        return a + b

    def test_reference_by_length_others_provisional(self):
        with tempfile.TemporaryDirectory() as d:
            rows, excluded = auto_number(self.candidates(), d)
        by = {r["seq_id"]: r for r in rows}
        # HAP2 has more sequences >= min_len, but its 0.2 Mb one is not chromosome-scale;
        # HAP1 and HAP2 tie at three, so the first in the manifest (HAP1) is the reference
        self.assertEqual([by[f"a{i}"]["chrom"] for i in range(3)], [1, 2, 3])
        self.assertEqual(sorted(by[f"b{i}"]["chrom"] for i in range(3)),
                         [AUTO_OFFSET, AUTO_OFFSET + 1, AUTO_OFFSET + 2])
        self.assertEqual(by["b1"]["chrom"], AUTO_OFFSET)  # longest first
        self.assertEqual([e[0] for e in excluded], ["b3"])
        self.assertTrue(all(r["numbering"] == "auto" for r in rows))

    def test_resume_reuses_assigned_numbers(self):
        with tempfile.TemporaryDirectory() as d:
            rows, _ = auto_number(self.candidates(), d)
            # the matrix stage matched HAP2's provisional numbers to 3, 1, 2
            final = {"b0": 3, "b1": 1, "b2": 2}
            with open(os.path.join(d, "sequences.tsv"), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["unit_id", "hap", "chrom", "seq_id", "source", "numbering"],
                                   delimiter="\t")
                w.writeheader()
                for r in rows:
                    chrom = final.get(r["seq_id"], r["chrom"])
                    w.writerow(dict(unit_id="", hap=r["hap"], chrom=chrom, seq_id=r["seq_id"],
                                    source=r["source"], numbering="auto"))
            again, _ = auto_number(self.candidates(), d)
        self.assertEqual({r["seq_id"]: r["chrom"] for r in again if r["hap"] == "HAP2"}, final)


class TestAutoReconcile(unittest.TestCase):
    def unit(self, uid, hap, chrom, source):
        return dict(unit_id=uid, hap=hap, chrom=str(chrom), seq_id=uid, length="1",
                    source=source, desc="", numbering="auto")

    def test_provisional_numbers_take_their_match_even_near_homeologs(self):
        # reference HAP1 chr1, chr2 are homeologs (0.05 apart); HAP2's two
        # sequences are allelic copies (0.01) but only ~1.2x closer to their
        # true chromosome than to its homeolog -- below the 3x needed to
        # relabel a named unit, accepted for an auto-numbered one
        units = [self.unit("HAP1_chr01", "HAP1", 1, "h1"), self.unit("HAP1_chr02", "HAP1", 2, "h1"),
                 self.unit("HAP2_chr1001", "HAP2", AUTO_OFFSET, "h2"),
                 self.unit("HAP2_chr1002", "HAP2", AUTO_OFFSET + 1, "h2")]
        d = [[0, 0.05, 0.060, 0.050],
             [0.05, 0, 0.050, 0.060],
             [0.060, 0.050, 0, 0.05],
             [0.050, 0.060, 0.05, 0]]
        corrections, ambiguous = detect_relabeling(units, d)
        self.assertEqual(ambiguous, [])
        self.assertEqual({c["unit_id"]: c["new_chrom"] for c in corrections},
                         {"HAP2_chr1001": "2", "HAP2_chr1002": "1"})

    def test_unmatched_extra_sequence_keeps_its_provisional_number(self):
        units = [self.unit("HAP1_chr01", "HAP1", 1, "h1"),
                 self.unit("HAP2_chr1001", "HAP2", AUTO_OFFSET, "h2"),
                 self.unit("HAP2_chr1002", "HAP2", AUTO_OFFSET + 1, "h2")]
        d = [[0, 0.01, 0.12], [0.01, 0, 0.12], [0.12, 0.12, 0]]
        corrections, _ = detect_relabeling(units, d)
        self.assertEqual([(c["unit_id"], c["new_chrom"]) for c in corrections], [("HAP2_chr1001", "1")])


if __name__ == "__main__":
    unittest.main()


class TestReferenceChoice(unittest.TestCase):
    def test_contig_level_haplotype_is_never_the_reference(self):
        contigs = [seq("HAP1", f"c{i}", 5e6 + i * 1e5, "h1.fa") for i in range(40)]
        chroms = [seq("HAP2", f"g{i}", 60e6 - i * 2e6, "h2.fa") for i in range(10)]
        with tempfile.TemporaryDirectory() as d:
            rows, _ = auto_number(contigs + chroms, d)
        self.assertEqual({r["chrom"] for r in rows if r["hap"] == "HAP2"}, set(range(1, 11)))


class TestLimitationNotes(unittest.TestCase):
    def test_contig_level_haplotype_is_flagged(self):
        from ploidyspec.manifest import contig_level_notes
        rows = [dict(source="h1.fa", hap="HAP1")] * 33 + [dict(source="h2.fa", hap="HAP2")] * 300
        notes = contig_level_notes(rows)
        self.assertEqual(len(notes), 1)
        self.assertIn("HAP2", notes[0])
        self.assertIn("scaffold_by_reference.py", notes[0])

    def test_balanced_haplotypes_are_not_flagged(self):
        from ploidyspec.manifest import contig_level_notes
        rows = [dict(source="h1.fa", hap="HAP1")] * 33 + [dict(source="h2.fa", hap="HAP2")] * 34
        self.assertEqual(contig_level_notes(rows), [])
