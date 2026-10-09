import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.residual import allele_levels, near_allelic_track, null_max_run, residual_tetrasomy, runs
from ploidyspec.windowed import FIELDNAMES, control_pairs

W = 250_000


class TestRuns(unittest.TestCase):
    def test_bridges_single_gaps_only(self):
        flags = [1, 1, 0, 1, 0, 0, 1]
        self.assertEqual(runs([bool(f) for f in flags]), [(0, 3, 3), (6, 6, 1)])

    def test_null_grows_with_density(self):
        sparse = null_max_run([100] * 10, 20, 100, 0.95)
        dense = null_max_run([100] * 10, 300, 100, 0.95)
        self.assertLess(sparse, dense)
        self.assertLessEqual(sparse, 3)


class TestControl(unittest.TestCase):
    def test_window_close_to_the_control_is_not_near_allelic(self):
        windows = {0: (W - 1, {"B": 0.001}), W: (2 * W - 1, {"B": 0.002}), 2 * W: (3 * W - 1, {"B": 0.05})}
        control = {0: (W - 1, {"C": 0.001}), W: (2 * W - 1, {"C": 0.06})}
        track = near_allelic_track(windows, 0.003, 3.0, control)
        self.assertEqual([t[2] for t in track], [False, True, False])
        self.assertEqual([t[3] for t in track], [True, False, False])

    def test_ceiling_caps_an_inflated_allelic_level(self):
        windows = {0: (W - 1, {"HAP1_chr02": 0.04})}
        self.assertTrue(near_allelic_track(windows, 0.02, 3.0)[0][2])
        self.assertFalse(near_allelic_track(windows, 0.02, 3.0, ceiling=0.0125)[0][2])

    def test_control_pairs_cover_every_paired_chromosome(self):
        ctl = control_pairs([(1, 8), (2, 9), (3, 10)], range(1, 11))
        self.assertEqual(ctl, [(1, 9), (2, 10), (3, 8)])
        anchors = {a for a, _ in ctl} | {c for _, c in ctl}
        self.assertEqual(anchors, {1, 2, 3, 8, 9, 10})
        self.assertEqual(control_pairs([(1, 2)], [1, 2, 3]), [(1, 3), (2, 3)])


def write_windows(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def row(a, b, start, d):
    return dict(group="g", win_start=start, win_end=start + W - 1, unit_a=a, unit_b=b,
                hap_a=a.split("_")[0], hap_b=b.split("_")[0], chrom_a=int(a[-2:]), chrom_b=int(b[-2:]),
                kmers_a=W, shared=0, containment=0, distance=d)


class TestAlleleLevels(unittest.TestCase):
    """Four copies of chr01, HAP4 a divergent subgenome copy: every copy is
    measured at the close copies' allelic level, not one inflated by HAP4."""

    def test_divergent_copy_does_not_inflate_the_level(self):
        units = [f"HAP{h}_chr01" for h in (1, 2, 3, 4)]
        rows = []
        for a in units:
            for b in units:
                if a != b:
                    d = 0.04 if "HAP4" in (a, b) else 0.003
                    rows += [row(a, b, i * W, d) for i in range(5)]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "windowed_all.tsv")
            write_windows(path, rows)
            levels = allele_levels(path, {u: 1 for u in units})
        self.assertEqual(set(levels.values()), {0.003})


class TestPairVote(unittest.TestCase):
    """chr01/chr02 homeologs, 40 windows each, both haplotypes. The last 12
    windows are near-allelic from all four copies (residual tetrasomy) on
    chr01/chr02; a misplaced contig makes 12 windows near-allelic from only
    two copies of chr03/chr04 -- not enough."""

    def build(self, d, exchange=False):
        units = [f"HAP{h}_chr{c:02d}" for c in (1, 2, 3, 4) for h in (1, 2)]
        partner = {1: 2, 2: 1, 3: 4, 4: 3}
        homeo, allelic, ctl = [], [], []
        for u in units:
            h, c = u.split("_")
            c = int(c[3:])
            other_hap = "HAP2" if h == "HAP1" else "HAP1"
            for i in range(40):
                s = i * W
                near = (c in (1, 2) and i >= 28) or (u in ("HAP1_chr03", "HAP2_chr04") and i >= 28)
                for hb in ("HAP1", "HAP2"):
                    d_h = 0.004 if near else 0.05
                    if exchange and u == "HAP1_chr03" and i < 12:
                        d_h = 0.003
                    homeo.append(row(u, f"{hb}_chr{partner[c]:02d}", s, d_h))
                # with exchange: HAP1_chr03's first 12 windows carry chr04 sequence
                swapped = exchange and u == "HAP1_chr03" and i < 12
                allelic.append(row(u, f"{other_hap}_chr{c:02d}", s, 0.05 if swapped else 0.003))
                ctl.append(row(u, f"HAP1_chr{(c % 4) + 1:02d}", s, 0.06))
        os.makedirs(os.path.join(d, "homeologs"))
        os.makedirs(os.path.join(d, "windowed"))
        write_windows(os.path.join(d, "homeologs", "windowed_homeologs_all.tsv"), homeo)
        write_windows(os.path.join(d, "homeologs", "windowed_homeolog_controls.tsv"), ctl)
        write_windows(os.path.join(d, "windowed", "windowed_all.tsv"), allelic)
        return {u: int(u[-2:]) for u in units}

    def test_residual_needs_most_copies_of_the_pair(self):
        with tempfile.TemporaryDirectory() as d:
            chrom_of = self.build(d)
            rows, readings, null = residual_tetrasomy(d, chrom_of, dict(min_segment_bp=1_000_000))
        self.assertLess(null, 12)
        self.assertEqual(readings[1]["residual_bp"], 12 * W)
        self.assertEqual(readings[1]["residual_support"], "chr02:4/4")
        self.assertIn("(end)", readings[1]["residual_segments"])
        self.assertEqual(readings[3]["residual_bp"], 0)
        self.assertEqual(readings[3]["residual_support"], "chr04:2/4")
        self.assertEqual(readings[1]["controlled"], "yes")


    def test_homeologous_exchange_is_read_from_one_copy(self):
        with tempfile.TemporaryDirectory() as d:
            chrom_of = self.build(d, exchange=True)
            rows, readings, _ = residual_tetrasomy(d, chrom_of, dict(min_segment_bp=1_000_000))
        self.assertEqual(readings[3]["exchange_bp"], 12 * W)
        self.assertEqual(readings[3]["exchange_support"], "chr04:1/2")
        self.assertIn("HAP1:0.0-3.0Mb(start)~chr04", readings[3]["exchange_segments"])
        self.assertEqual(readings[1]["exchange_bp"], 0)
        # the carrier's exchanged windows are not near-allelic (own homolog far)
        track = {r["anchor"]: r for r in rows if r["partner"] == "chr04"}["HAP1_chr03"]
        self.assertEqual(track["n_near_allelic"], 12)  # only its last 12 (misplaced-contig) windows


if __name__ == "__main__":
    unittest.main()
