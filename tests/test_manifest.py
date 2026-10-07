import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.manifest import read_manifest


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


if __name__ == "__main__":
    unittest.main()
