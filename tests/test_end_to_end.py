"""
End-to-end: simulate genomes with a known answer, run `ploidyspec all` on each,
check the pipeline recovers it. Needs samtools + FastK on PATH and takes a few
minutes, so it only runs when PLOIDYSPEC_E2E=1:

    PLOIDYSPEC_E2E=1 python3 -m unittest tests.test_end_to_end
"""

import csv
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.cli import main
from ploidyspec.simulate import SCENARIOS, simulate

TOOLS = ("samtools", "FastK", "Logex", "Histex", "Tabex")
ENABLED = os.environ.get("PLOIDYSPEC_E2E") == "1" and all(shutil.which(t) for t in TOOLS)
N_CHROM = 7
SIM_PARAMS = dict(n_chrom=N_CHROM, chrom_len=800_000)
RUN_ARGS = ["--k", "15,23", "--min-len", "400000", "--window", "100000",
            "--min-segment-bp", "200000", "--with-te-markers", "--threads", "4"]


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


@unittest.skipUnless(ENABLED, "set PLOIDYSPEC_E2E=1 with samtools + FastK on PATH")
class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="ploidyspec_e2e_")
        simulate("all", cls.tmp, seed=7, params=SIM_PARAMS)
        for name in SCENARIOS:
            sdir = os.path.join(cls.tmp, name)
            main(["all", "--manifest", os.path.join(sdir, "manifest.tsv"),
                  "--outdir", os.path.join(sdir, "out")] + RUN_ARGS)

    @classmethod
    def tearDownClass(cls):
        if os.environ.get("PLOIDYSPEC_E2E_KEEP") != "1":
            shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, scenario, *parts):
        return os.path.join(self.tmp, scenario, "out", *parts)

    def states(self, scenario):
        rows = read_tsv(self.out(scenario, "rediploidization", "rediploidization_by_chrom.tsv"))
        return {r["chrom"]: r["copy_state"] for r in rows}

    def expected_states(self, scenario):
        return {r["chrom"]: r["copy_state"]
                for r in read_tsv(os.path.join(self.tmp, scenario, "truth.tsv"))}

    def test_copy_counts(self):
        for scenario, copies in (("diploid", "2"), ("autotetraploid", "4")):
            rows = read_tsv(self.out(scenario, "matrix", "ploidy_summary.tsv"))
            self.assertEqual({r["n_haplotype_copies"] for r in rows}, {copies}, scenario)

    def test_states_match_truth(self):
        for scenario in ("diploid", "autotetraploid", "autotetraploid_2hap", "allotetraploid"):
            self.assertEqual(self.states(scenario), self.expected_states(scenario), scenario)

    def test_rediploidized_states(self):
        self.assertEqual(self.states("rediploidized"), self.expected_states("rediploidized"))

    def test_rediploidized_fusions(self):
        rows = read_tsv(self.out("rediploidized", "rediploidization", "fusions.tsv"))
        found = {(r["hap"], r["components"]) for r in rows if r["status"] == "fusion"}
        expected = {(h, c) for h in ("HAP3", "HAP4") for c in ("chr03+chr04", "chr05+chr06")}
        self.assertEqual(found, expected)

    def test_homeolog_pairs_in_two_haplotype_tetraploids(self):
        expected = {(f"chr{i:02d}", f"chr{i + N_CHROM:02d}") for i in range(1, N_CHROM + 1)}
        for scenario in ("allotetraploid", "autotetraploid_2hap"):
            rows = read_tsv(self.out(scenario, "homeologs", "homeolog_pairs.tsv"))
            self.assertEqual({(r["chrom_a"], r["chrom_b"]) for r in rows}, expected, scenario)

    def test_no_fusions_without_rediploidization(self):
        for scenario in ("diploid", "autotetraploid", "autotetraploid_2hap", "allotetraploid"):
            rows = read_tsv(self.out(scenario, "rediploidization", "fusions.tsv"))
            self.assertEqual([r for r in rows if r["status"] == "fusion"], [], scenario)


if __name__ == "__main__":
    unittest.main()
