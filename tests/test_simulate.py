import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.simulate import Simulator, build_scenario

SMALL = dict(n_chrom=7, chrom_len=20_000, te_copies=5, te_len=200)


class TestSimulator(unittest.TestCase):
    def test_mutate_hits_requested_rate_and_always_changes_base(self):
        sim = Simulator(1)
        seq = sim.random_seq(200_000)
        out = sim.mutate(seq, 0.02)
        self.assertAlmostEqual(np.mean(out != seq), 0.02, delta=0.002)

    def test_mutate_respects_region(self):
        sim = Simulator(1)
        seq = sim.random_seq(10_000)
        out = sim.mutate(seq, 0.5, 0, 4_000)
        self.assertTrue(np.array_equal(out[4_000:], seq[4_000:]))
        self.assertGreater(np.mean(out[:4_000] != seq[:4_000]), 0.4)

    def test_te_insertion_adds_copies(self):
        sim = Simulator(1, te_copy_div=0.0)
        fam = sim.te_family()
        seq = sim.insert_copies(sim.random_seq(10_000), fam, 7)
        self.assertEqual(len(seq), 10_000 + 7 * len(fam))


class TestScenarios(unittest.TestCase):
    def test_diploid_and_autotetraploid_haplotype_counts(self):
        for name, n_hap in (("diploid", 2), ("autotetraploid", 4)):
            haps, truth = build_scenario(name, 1, SMALL)
            self.assertEqual(len(haps), n_hap)
            self.assertEqual(len(truth), SMALL["n_chrom"])

    def test_allotetraploid_numbers_subgenomes_separately(self):
        haps, truth = build_scenario("allotetraploid", 1, SMALL)
        self.assertEqual(len(haps["HAP1"]), 2 * SMALL["n_chrom"])
        self.assertEqual({r["chrom"]: r["homeolog"] for r in truth}["chr01"], "chr08")

    def test_autotetraploid_2hap_numbers_copies_twice(self):
        haps, truth = build_scenario("autotetraploid_2hap", 1, SMALL)
        self.assertEqual(len(haps), 2)
        self.assertEqual(len(haps["HAP1"]), 2 * SMALL["n_chrom"])
        self.assertEqual({r["copy_state"] for r in truth}, {"tetrasomic_like"})

    def test_rediploidized_layout(self):
        haps, truth = build_scenario("rediploidized", 1, SMALL)
        names = [r[0] for r in haps["HAP3"]]
        self.assertIn("HAP3_fused_5_6", names)
        self.assertNotIn("HAP3_chr4", names)
        fused = {r[0]: r for r in haps["HAP3"]}
        self.assertEqual(fused["HAP3_fused_5_6"][1], "")  # unplaced: no chromosome tag
        self.assertGreater(len(fused["HAP3_chr3"][2]), 1.8 * len(haps["HAP1"][2][2]))
        states = {r["chrom"]: r["copy_state"] for r in truth}
        self.assertEqual(states["chr01"], "resolved_lineages")
        self.assertEqual(states["chr07"], "tetrasomic_like")

    def test_same_seed_same_genome(self):
        a, _ = build_scenario("diploid", 3, SMALL)
        b, _ = build_scenario("diploid", 3, SMALL)
        self.assertTrue(np.array_equal(a["HAP1"][0][2], b["HAP1"][0][2]))


if __name__ == "__main__":
    unittest.main()
