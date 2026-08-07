import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.whole_matrix import write_ploidy_and_homology_reports


class TestPloidyAndHomologyReports(unittest.TestCase):
    def _run(self, units, summaries):
        ids = [u["unit_id"] for u in units]
        with tempfile.TemporaryDirectory() as td:
            write_ploidy_and_homology_reports(td, ids, units, summaries)
            with open(os.path.join(td, "matrix", "ploidy_summary.tsv")) as f:
                ploidy_rows = list(csv.DictReader(f, delimiter="\t"))
            with open(os.path.join(td, "matrix", "homologous_chromosomes.tsv")) as f:
                homolog_rows = list(csv.DictReader(f, delimiter="\t"))
        return ploidy_rows, homolog_rows

    def test_detects_mixed_ploidy_across_chromosomes(self):
        units = [
            dict(unit_id="HAP1_chr01", hap="HAP1", chrom=1),
            dict(unit_id="HAP2_chr01", hap="HAP2", chrom=1),
            dict(unit_id="HAP3_chr01", hap="HAP3", chrom=1),
            dict(unit_id="HAP4_chr01", hap="HAP4", chrom=1),
            dict(unit_id="HAP1_chr02", hap="HAP1", chrom=2),
            dict(unit_id="HAP2_chr02", hap="HAP2", chrom=2),
        ]
        ids = [u["unit_id"] for u in units]
        n = len(ids)
        summaries = {
            (i, j): dict(distance=0.01, resolution_limited=False)
            for i in range(n)
            for j in range(i + 1, n)
        }
        ploidy_rows, homolog_rows = self._run(units, summaries)

        by_chrom = {r["chrom"]: r for r in ploidy_rows}
        self.assertEqual(by_chrom["chr01"]["n_haplotype_copies"], "4")
        self.assertEqual(by_chrom["chr02"]["n_haplotype_copies"], "2")
        self.assertEqual(by_chrom["chr01"]["haplotype_labels"], "HAP1,HAP2,HAP3,HAP4")

        # chr01 has C(4,2)=6 homologous pairs, chr02 has C(2,2)... =1
        chr01_pairs = [r for r in homolog_rows if r["chrom"] == "chr01"]
        chr02_pairs = [r for r in homolog_rows if r["chrom"] == "chr02"]
        self.assertEqual(len(chr01_pairs), 6)
        self.assertEqual(len(chr02_pairs), 1)

    def test_homologous_pairs_never_cross_chromosome_numbers(self):
        units = [
            dict(unit_id="HAP1_chr01", hap="HAP1", chrom=1),
            dict(unit_id="HAP2_chr01", hap="HAP2", chrom=1),
            dict(unit_id="HAP1_chr02", hap="HAP1", chrom=2),
            dict(unit_id="HAP2_chr02", hap="HAP2", chrom=2),
        ]
        ids = [u["unit_id"] for u in units]
        n = len(ids)
        summaries = {
            (i, j): dict(distance=0.05, resolution_limited=False)
            for i in range(n)
            for j in range(i + 1, n)
        }
        _, homolog_rows = self._run(units, summaries)
        for r in homolog_rows:
            self.assertTrue(r["unit_a"].endswith(r["chrom"]))
            self.assertTrue(r["unit_b"].endswith(r["chrom"]))

    def test_distance_and_resolution_limited_carried_through(self):
        units = [
            dict(unit_id="HAP1_chr01", hap="HAP1", chrom=1),
            dict(unit_id="HAP2_chr01", hap="HAP2", chrom=1),
        ]
        ids = [u["unit_id"] for u in units]
        summaries = {(0, 1): dict(distance=0.0099, resolution_limited=False)}
        _, homolog_rows = self._run(units, summaries)
        self.assertEqual(len(homolog_rows), 1)
        self.assertAlmostEqual(float(homolog_rows[0]["distance"]), 0.0099)
        self.assertEqual(homolog_rows[0]["resolution_limited"], "False")

    def test_unresolved_distance_written_as_empty(self):
        units = [
            dict(unit_id="HAP1_chr01", hap="HAP1", chrom=1),
            dict(unit_id="HAP2_chr01", hap="HAP2", chrom=1),
        ]
        summaries = {(0, 1): dict(distance=None, resolution_limited=True)}
        _, homolog_rows = self._run(units, summaries)
        self.assertEqual(homolog_rows[0]["distance"], "")
        self.assertEqual(homolog_rows[0]["resolution_limited"], "True")


if __name__ == "__main__":
    unittest.main()
