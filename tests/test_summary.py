import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.summary import summarize_species


def write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def species(d, copies, states, pairs=(), partition=(), fusions=(), outlier=""):
    n = len(states)
    write(os.path.join(d, "matrix", "ploidy_summary.tsv"), ["chrom", "n_haplotype_copies", "haplotype_labels"],
          [[f"chr{i:02d}", copies, ""] for i in range(1, n + 1)])
    write(os.path.join(d, "matrix", "whole_chrom_pairs.tsv"), ["unit_a", "unit_b", "distance"],
          [["HAP1_chr01", "HAP2_chr01", "0.01"], ["HAP1_chr01", "HAP1_chr02", "0.12"]])
    write(os.path.join(d, "homeologs", "homeolog_pairs.tsv"), ["chrom_a", "chrom_b", "mean_distance"], pairs)
    write(os.path.join(d, "structure", "genome_partition.tsv"),
          ["species", "k", "group_sizes", "separation_ratio", "null_mean_ratio", "z_score", "groups"], partition)
    write(os.path.join(d, "rediploidization", "rediploidization_by_chrom.tsv"),
          ["chrom", "copy_state", "state_basis", "split_extent"],
          [[f"chr{i:02d}", s, "copies" if copies >= 3 else "homeolog_pool", "none"] for i, s in enumerate(states, 1)])
    write(os.path.join(d, "rediploidization", "rediploidization_summary.tsv"), ["metric", "value"],
          [["most_frequent_outlier_hap", outlier]] if outlier else [])
    write(os.path.join(d, "rediploidization", "fusions.tsv"), ["scaffold", "components", "status"], fusions)


def answers(d):
    return {r["question"]: r for r in summarize_species(d)}


class TestSummary(unittest.TestCase):
    def test_tetrasomic_autotetraploid(self):
        with tempfile.TemporaryDirectory() as d:
            species(d, 4, ["tetrasomic_like"] * 9 + ["candidate"])
            a = answers(d)
            self.assertIn("4 haplotype copies", a["ploidy"]["answer"])
            self.assertIn("Auto-like", a["origin_like_structure"]["answer"])
            self.assertEqual(a["origin_like_structure"]["confidence"], "high")
            self.assertIn("No sign", a["rediploidization"]["answer"])

    def test_fusions_mean_rediploidization(self):
        with tempfile.TemporaryDirectory() as d:
            species(d, 4, ["fusion_lineages"] * 2 + ["tetrasomic_like"] * 8,
                    fusions=[["HAP3_chr19", "chr19+chr22", "fusion"]])
            self.assertEqual(answers(d)["rediploidization"]["confidence"], "high")

    def test_one_haplotype_fusion_in_a_diploid_is_not_rediploidization(self):
        # xgMonCant1: HAP1_chr23 contains chr23 and chr24, which makes them the
        # only "homeolog pair"; nothing else is duplicated
        with tempfile.TemporaryDirectory() as d:
            species(d, 2, ["fusion_lineages"] * 2 + ["not_assessable"] * 8,
                    pairs=[["chr01", "chr02", "0.07"]],
                    fusions=[["HAP1_chr01", "chr01+chr02", "fusion"]])
            a = answers(d)
            self.assertTrue(a["rediploidization"]["answer"].startswith("Not read as rediploidization"))
            self.assertTrue(a["origin_like_structure"]["answer"].startswith("Not assessable"))
            self.assertIn("no older duplicated sets", a["ploidy"]["answer"])

    def test_two_haplotype_allotetraploid_from_pairs(self):
        with tempfile.TemporaryDirectory() as d:
            pairs = [[f"chr{i:02d}", f"chr{i + 4:02d}", "0.05"] for i in range(1, 5)]
            part = [["sp", "4", "2,2,2,2", "2.0", "1.0", "20.0", "x"]]
            species(d, 2, ["resolved_lineages"] * 8, pairs=pairs, partition=part)
            a = answers(d)
            self.assertIn("sets of 2", a["ploidy"]["answer"])
            self.assertIn("x ≈ 4", a["ploidy"]["answer"])
            self.assertIn("Allo-like", a["origin_like_structure"]["answer"])
            self.assertIn("rediploidized autopolyploid", a["origin_like_structure"]["answer"])

    def test_same_odd_haplotype_reads_as_assembly(self):
        with tempfile.TemporaryDirectory() as d:
            species(d, 4, ["one_divergent_copy"] * 8, outlier="HAP1 (8 chromosomes)")
            self.assertIn("assembly", answers(d)["origin_like_structure"]["answer"])

    def test_diploid_is_not_assessable(self):
        with tempfile.TemporaryDirectory() as d:
            species(d, 2, ["not_assessable"] * 10)
            a = answers(d)
            self.assertEqual(a["origin_like_structure"]["confidence"], "not assessable")
            self.assertIn("no older duplicated sets", a["ploidy"]["answer"])


if __name__ == "__main__":
    unittest.main()
