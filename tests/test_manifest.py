import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.manifest import match_chrom, read_manifest


class TestReadManifest(unittest.TestCase):
    def test_relative_paths_resolve_against_manifest_dir(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "manifest.tsv")
            with open(path, "w") as f:
                f.write("# fasta_path\thap_label\n")
                f.write("sp.HAP1.fa\tHAP1\n")
                f.write("/abs/sp.HAP2.fa\tHAP2\n")
            rows = read_manifest(path)
        self.assertEqual(rows, [(os.path.join(d, "sp.HAP1.fa"), "HAP1"),
                                ("/abs/sp.HAP2.fa", "HAP2")])


class TestMatchChrom(unittest.TestCase):
    def test_alternation_reads_the_group_that_matched(self):
        pat = [r"chromosome:?\s*(\d+)\b|SUPER[_-](\d+)_HAP\d+\b|SUPER[_-](\d+)\b(?!_)"]
        self.assertEqual(match_chrom("OZ282501.1 Buddleja davidii genome assembly, chromosome: 1", pat), 1)
        self.assertEqual(match_chrom("SUPER_7_HAP2", pat), 7)
        self.assertEqual(match_chrom("SUPER_12", pat), 12)
        self.assertIsNone(match_chrom("scaffold_3", pat))


if __name__ == "__main__":
    unittest.main()
