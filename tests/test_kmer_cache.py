import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.kmer_tables import (
    cached_fasta_matches,
    chrom_fasta_path,
    remove_unit_cache,
    unit_signature,
)


def unit(seq_id="SUPER_1_HAP1", length=150, source="/data/asm.fa.gz"):
    return dict(unit_id="HAP1_chr01", seq_id=seq_id, length=str(length), source=source)


class TestCachedFasta(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.tmp, "chroms"))
        self.fa = chrom_fasta_path(self.tmp, "HAP1_chr01")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def write_fasta(self, seq_id, length, width=60):
        seq = "A" * length
        with open(self.fa, "w") as f:
            f.write(f">{seq_id}\n")
            for i in range(0, length, width):
                f.write(seq[i : i + width] + "\n")

    def test_legacy_file_matching_name_and_length_is_reused(self):
        self.write_fasta("SUPER_1_HAP1", 150)
        self.assertTrue(cached_fasta_matches(self.fa, unit()))

    def test_legacy_file_with_other_width_is_reused(self):
        self.write_fasta("SUPER_1_HAP1", 150, width=80)
        self.assertTrue(cached_fasta_matches(self.fa, unit()))

    def test_renumbered_chromosome_is_stale(self):
        # drLytSali1: HAP1_chr01.fa held SUPER_2_HAP1 after renumbering
        self.write_fasta("SUPER_2_HAP1", 150)
        self.assertFalse(cached_fasta_matches(self.fa, unit()))

    def test_same_name_different_length_is_stale(self):
        self.write_fasta("SUPER_1_HAP1", 140)
        self.assertFalse(cached_fasta_matches(self.fa, unit()))

    def test_signature_decides_when_present(self):
        self.write_fasta("SUPER_1_HAP1", 150)
        with open(self.fa + ".src", "w") as f:
            f.write(unit_signature(unit(source="/data/old.fa.gz")))
        self.assertFalse(cached_fasta_matches(self.fa, unit()))
        with open(self.fa + ".src", "w") as f:
            f.write(unit_signature(unit()))
        self.assertTrue(cached_fasta_matches(self.fa, unit()))

    def test_remove_unit_cache_clears_tables_at_every_k(self):
        self.write_fasta("SUPER_1_HAP1", 150)
        keep = []
        for k in (13, 23):
            kdir = os.path.join(self.tmp, f"ktabs_k{k}")
            os.makedirs(kdir)
            for name in ("HAP1_chr01.ktab", ".HAP1_chr01.ktab.1", "HAP1_chr01.hist",
                         "HAP1_chr010.ktab", "HAP2_chr01.ktab"):
                open(os.path.join(kdir, name), "w").close()
            keep += [os.path.join(kdir, "HAP1_chr010.ktab"), os.path.join(kdir, "HAP2_chr01.ktab")]
        remove_unit_cache(self.tmp, "HAP1_chr01")
        self.assertFalse(os.path.exists(self.fa))
        remaining = sorted(os.path.join(d, n) for d in (os.path.join(self.tmp, "ktabs_k13"),
                                                         os.path.join(self.tmp, "ktabs_k23"))
                           for n in os.listdir(d))
        self.assertEqual(remaining, sorted(keep))


if __name__ == "__main__":
    unittest.main()
