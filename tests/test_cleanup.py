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
