import csv
import gzip
import os
import stat
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.cli import write_fasta_manifest
from ploidyspec.manifest import gzip_kind, prepare, read_manifest, usable_fasta


def write_fasta(path, records):
    """records: [(seq_id, description, length)]; writes the FASTA and its .fai
    (so prepare needs no samtools)."""
    with open(path, "w") as f, open(path + ".fai", "w") as fai:
        for seq_id, desc, length in records:
            f.write(f">{seq_id} {desc}\n")
            offset = f.tell()
            f.write("A" * length + "\n")
            fai.write(f"{seq_id}\t{length}\t{offset}\t{length}\t{length + 1}\n")


class TestUsableFasta(unittest.TestCase):
    def test_plain_gzip_is_staged_in_outdir(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "hap1.fa.gz")
            with gzip.open(src, "wt") as f:
                f.write(">s1\nACGT\n")
            self.assertEqual(gzip_kind(src), "gzip")
            out = usable_fasta(src, os.path.join(d, "out"), "samtools")
            self.assertTrue(out.startswith(os.path.join(d, "out", "inputs")))
            self.assertIn(gzip_kind(out), ("bgzf", None))  # bgzip if found, else plain

    def test_read_only_directory_is_linked(self):
        with tempfile.TemporaryDirectory() as d:
            ro = os.path.join(d, "ro")
            os.makedirs(ro)
            src = os.path.join(ro, "hap1.fa")
            with open(src, "w") as f:
                f.write(">s1\nACGT\n")
            os.chmod(ro, stat.S_IRUSR | stat.S_IXUSR)
            try:
                if os.access(ro, os.W_OK):
                    self.skipTest("running as a user who can write anywhere")
                out = usable_fasta(src, os.path.join(d, "out"), "samtools")
                self.assertTrue(os.path.islink(out))
                self.assertEqual(os.path.realpath(out), os.path.realpath(src))
            finally:
                os.chmod(ro, stat.S_IRWXU)

    def test_indexed_or_writable_fasta_is_used_in_place(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "hap1.fa")
            with open(src, "w") as f:
                f.write(">s1\nACGT\n")
            self.assertEqual(usable_fasta(src, os.path.join(d, "out"), "samtools"), src)


class TestFastaManifest(unittest.TestCase):
    def test_one_haplotype_per_file_in_order(self):
        with tempfile.TemporaryDirectory() as d:
            fas = [os.path.join(d, f"h{i}.fa") for i in (1, 2)]
            for fa in fas:
                open(fa, "w").close()
            rows = read_manifest(write_fasta_manifest(fas, d))
        self.assertEqual(rows, [(fas[0], "HAP1"), (fas[1], "HAP2")])


class TestDetectNaming(unittest.TestCase):
    """xgLitLitt3: hap1 headers carry 'chromosome: N', hap2's scaffolds only
    'contig: RL_2'. detect numbers both automatically instead of dropping hap2."""

    def test_unnamed_haplotype_switches_to_auto(self):
        with tempfile.TemporaryDirectory() as d:
            h1, h2 = os.path.join(d, "h1.fa"), os.path.join(d, "h2.fa")
            write_fasta(h1, [(f"OZ{i}", f"genome assembly, chromosome: {i}", 3000 - i * 100) for i in (1, 2, 3)])
            write_fasta(h2, [(f"CA{i}", f"genome assembly, contig: RL_{i}", 3000 - i * 100) for i in (1, 2, 3)])
            manifest = write_fasta_manifest([h1, h2], d)
            prepare(manifest, os.path.join(d, "out"), "samtools", 1000, None, None)
            with open(os.path.join(d, "out", "sequences.tsv")) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(sorted({r["hap"] for r in rows}), ["HAP1", "HAP2"])
        self.assertTrue(all(r["numbering"] == "auto" for r in rows))

    def test_named_haplotypes_keep_their_names(self):
        with tempfile.TemporaryDirectory() as d:
            h1, h2 = os.path.join(d, "h1.fa"), os.path.join(d, "h2.fa")
            write_fasta(h1, [(f"SUPER_{i}", "", 2000) for i in (1, 2)])
            write_fasta(h2, [(f"SUPER_{i}_HAP2", "", 2000) for i in (1, 2)])
            prepare(write_fasta_manifest([h1, h2], d), os.path.join(d, "out"), "samtools", 1000, None, None)
            with open(os.path.join(d, "out", "sequences.tsv")) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(sorted(r["unit_id"] for r in rows),
                         ["HAP1_chr01", "HAP1_chr02", "HAP2_chr01", "HAP2_chr02"])
        self.assertTrue(all(r["numbering"] == "names" for r in rows))


if __name__ == "__main__":
    unittest.main()
