import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.common import cleanup_intermediates


class TestCleanup(unittest.TestCase):
    def test_removes_intermediates_keeps_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            for sub in ("ktabs_k15", "ktabs_k23", "chroms", "tmp_matrix",
                        "windowed/tmp_windowed", "matrix", "rediploidization"):
                os.makedirs(os.path.join(d, sub))
            for keep in ("matrix/whole_chrom_distance_matrix.csv", "sequences.tsv",
                         "rediploidization/fusions.tsv"):
                open(os.path.join(d, keep), "w").close()
            open(os.path.join(d, "ktabs_k15", "HAP1_chr01.ktab"), "w").close()

            removed = cleanup_intermediates(d)

            self.assertEqual(len(removed), 5)
            left = sorted(os.path.relpath(os.path.join(r, f), d)
                          for r, _, files in os.walk(d) for f in files)
            self.assertEqual(left, ["matrix/whole_chrom_distance_matrix.csv",
                                    "rediploidization/fusions.tsv", "sequences.tsv"])


if __name__ == "__main__":
    unittest.main()


class TestWindowedHomeologsWithoutPairs(unittest.TestCase):
    def test_no_pairs_is_skipped_not_fatal(self):
        # dcCerAlpi1: `all --with-windowed-homeologs` stopped at this stage
        # because the species has no accepted homeolog pairs
        from ploidyspec.cli import main
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "sequences.tsv"), "w").write("unit_id\thap\tchrom\tseq_id\tlength\tsource\tdesc\n")
            os.makedirs(os.path.join(d, "homeologs"))
            open(os.path.join(d, "homeologs", "homeolog_pairs.tsv"), "w").write("chrom_a\tchrom_b\tmean_distance\n")
            manifest = os.path.join(d, "m.tsv")
            open(manifest, "w").write("x.fa\tHAP1\n")
            main(["windowed-homeologs", "--manifest", manifest, "--outdir", d])
