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
from ploidyspec.panel import build_panel
from ploidyspec.simulate import SCENARIOS, simulate

TOOLS = ("samtools", "FastK", "Logex", "Histex", "Tabex", "Profex")
ENABLED = os.environ.get("PLOIDYSPEC_E2E") == "1" and all(shutil.which(t) for t in TOOLS)
N_CHROM = 7
SIM_PARAMS = dict(n_chrom=N_CHROM, chrom_len=800_000)
RUN_ARGS = ["--k", "15,23", "--min-len", "400000", "--window", "100000",
            # segment floors scaled to the 1.5 Mb simulated chromosomes (a 6-window residual stretch)
            "--min-segment-bp", "200000", "--residual-min-run", "3",
            "--with-te-markers", "--with-windowed-homeologs",
            "--threads", "4"]


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
            main(["all", "--manifest", os.path.join(cls.tmp, name, "manifest.tsv"),
                  "--outdir", os.path.join(cls.tmp, "results", name)] + RUN_ARGS)
        build_panel(os.path.join(cls.tmp, "results"), os.path.join(cls.tmp, "panel"))

    @classmethod
    def tearDownClass(cls):
        if os.environ.get("PLOIDYSPEC_E2E_KEEP") != "1":
            shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, scenario, *parts):
        return os.path.join(self.tmp, "results", scenario, *parts)

    def states(self, scenario):
        rows = read_tsv(self.out(scenario, "rediploidization", "rediploidization_by_chrom.tsv"))
        return {r["chrom"]: r["copy_state"] for r in rows}

    def expected_states(self, scenario):
        return {r["chrom"]: r["copy_state"]
                for r in read_tsv(os.path.join(self.tmp, scenario, "truth.tsv"))}

    def test_copy_counts(self):
        for scenario, copies in (("diploid", "2"), ("autotetraploid", "4"), ("mislabelled", "4")):
            rows = read_tsv(self.out(scenario, "matrix", "ploidy_summary.tsv"))
            self.assertEqual({r["n_haplotype_copies"] for r in rows}, {copies}, scenario)

    def test_states_match_truth(self):
        for scenario in ("diploid", "autotetraploid", "autotetraploid_2hap", "allotetraploid",
                         "mislabelled", "residual_tetrasomy", "allotriploid_one_file"):
            self.assertEqual(self.states(scenario), self.expected_states(scenario), scenario)

    def test_second_set_in_one_file_becomes_a_third_copy(self):
        seqs = read_tsv(self.out("allotriploid_one_file", "sequences.tsv"))
        self.assertEqual({r["hap"] for r in seqs}, {"HAP1", "HAP1B", "HAP2"})
        self.assertEqual({r["chrom"] for r in seqs}, {str(c) for c in range(1, N_CHROM + 1)})
        rows = {r["question"]: r for r in read_tsv(self.out("allotriploid_one_file", "summary.tsv"))}
        self.assertTrue(rows["origin_like_structure"]["answer"].startswith("AAB-like"))
        self.assertIn("odd copy number", rows["ploidy"]["answer"])

    def test_homeology_map_finds_arm_level_partners(self):
        summary = {r["metric"]: r["value"] for r in
                   read_tsv(self.out("segmental_homeology", "homeologs", "homeology_map_summary.tsv"))}
        self.assertIn("chr01:chr02+chr03", summary["multi_partner_list"])
        for scenario in ("diploid", "autotetraploid"):
            blocks = read_tsv(self.out(scenario, "homeologs", "homeolog_blocks.tsv"))
            self.assertEqual(blocks, [], scenario)

    def test_homeologous_exchange_found_in_one_copy(self):
        rows = {r["chrom"]: r for r in
                read_tsv(self.out("segmental_homeology", "rediploidization", "rediploidization_by_chrom.tsv"))}
        self.assertTrue(int(rows["chr04"]["exchange_bp"] or 0) > 0)
        self.assertIn("HAP1:", rows["chr04"]["exchange_segments"])
        self.assertIn("~chr05", rows["chr04"]["exchange_segments"])
        for c in ("chr02", "chr03", "chr06", "chr07"):
            self.assertEqual(rows[c]["exchange_bp"] or "0", "0", c)

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

    def test_mislabelled_chromosomes_are_relabelled(self):
        rows = read_tsv(self.out("mislabelled", "matrix", "chrom_label_corrections.tsv"))
        fixed = {(r["unit_id"], r["new_chrom"]) for r in rows if r["status"] == "corrected"}
        expected = {(f"HAP4_chr{c:02d}", str(c % N_CHROM + 1)) for c in range(1, N_CHROM + 1)}
        expected |= {("HAP3_chr05", "6"), ("HAP3_chr06", "5")}
        self.assertEqual(fixed, expected)

    def test_window_tracks_follow_whole_chromosome_distance(self):
        # every copy but HAP1 carries a deletion and an inversion, so equal-
        # coordinate windows would fall out of register halfway along; the
        # closest pair's windows must stay close all along the chromosome
        rows = read_tsv(self.out("autotetraploid", "windowed", "windowed_all.tsv"))
        far = [float(r["distance"]) for r in rows if float(r["distance"]) > 0.05]
        self.assertLess(len(far) / len(rows), 0.05)

    def test_no_fusions_without_rediploidization(self):
        for scenario in ("diploid", "autotetraploid", "autotetraploid_2hap", "allotetraploid"):
            rows = read_tsv(self.out(scenario, "rediploidization", "fusions.tsv"))
            self.assertEqual([r for r in rows if r["status"] == "fusion"], [], scenario)

    def test_panel_tables(self):
        summary = {r["species"]: r for r in
                   read_tsv(os.path.join(self.tmp, "panel", "panel_summary.tsv"))}
        self.assertEqual(set(summary), set(SCENARIOS))
        self.assertEqual(summary["autotetraploid"]["copies_per_chromosome"], "4")
        self.assertEqual(summary["allotetraploid"]["n_resolved_lineages"], str(2 * N_CHROM))
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "panel", "supplementary",
                                                    "te_markers_panel.tsv")))
        redip = {r["species"]: r for r in
                 read_tsv(os.path.join(self.tmp, "panel", "rediploidization_panel.tsv"))}
        self.assertEqual(redip["rediploidized"]["n_distinct_fusions"], "2")

    def residual(self, scenario):
        rows = read_tsv(self.out(scenario, "rediploidization", "rediploidization_by_chrom.tsv"))
        return {r["chrom"]: r for r in rows if r.get("residual_bp") not in ("", None)}

    def test_residual_tetrasomy_found_where_simulated(self):
        res = self.residual("residual_tetrasomy")
        self.assertEqual(len(res), 2 * N_CHROM)  # every pooled pair tested
        hit = {c for c, r in res.items() if int(r["residual_bp"])}
        self.assertEqual(hit, {"chr01", "chr02", f"chr{1 + N_CHROM:02d}", f"chr{2 + N_CHROM:02d}"})
        # the undiverged stretch is the chromosome's last 40%, not the satellite at its start
        for c in hit:
            self.assertNotIn("(start)", res[c]["residual_segments"], c)
            self.assertEqual(res[c]["residual_controlled"], "yes")

    def test_shared_satellite_is_excluded_by_the_control(self):
        rows = read_tsv(self.out("residual_tetrasomy", "rediploidization", "residual_tetrasomy.tsv"))
        self.assertTrue(all(int(r["n_shared_with_control"]) >= 1 for r in rows))

    def test_no_residual_tetrasomy_in_diverged_pairs(self):
        res = self.residual("allotetraploid")
        self.assertEqual(len(res), 2 * N_CHROM)
        self.assertEqual({c for c, r in res.items() if int(r["residual_bp"])}, set())

    def test_residual_tetrasomy_summary(self):
        rows = {r["question"]: r for r in read_tsv(self.out("residual_tetrasomy", "summary.tsv"))}
        self.assertTrue(rows["rediploidization"]["answer"].startswith("Partly"))
        self.assertTrue(rows["origin_like_structure"]["answer"].startswith("Auto-like origin"))

    def test_structure_outputs(self):
        for scenario in SCENARIOS:
            self.assertTrue(os.path.exists(self.out(scenario, "structure", "pair_synchrony.tsv")))
        rows = read_tsv(self.out("allotetraploid", "structure", "genome_partition.tsv"))
        # 7 homeolog pairs: the chromosome numbers factor into pairs
        self.assertIn(str(N_CHROM), {r["k"] for r in rows})


if __name__ == "__main__":
    unittest.main()


@unittest.skipUnless(ENABLED, "set PLOIDYSPEC_E2E=1 with samtools + FastK on PATH")
class TestRerunOnChangedAssembly(unittest.TestCase):
    """A re-run on a changed assembly in the same output directory must rebuild
    the cached per-unit FASTA, not reuse it by name (drLytSali1, 2026-09)."""

    def test_changed_sequence_is_re_extracted(self):
        from ploidyspec.simulate import build_scenario, write_fasta

        tmp = tempfile.mkdtemp(prefix="ploidyspec_rerun_")
        try:
            params = dict(n_chrom=2, chrom_len=60_000, te_copies=5, te_len=200)
            haps, _ = build_scenario("diploid", 1, params)
            manifest = os.path.join(tmp, "manifest.tsv")
            with open(manifest, "w") as m:
                for hap, recs in haps.items():
                    path = os.path.join(tmp, f"{hap}.fa")
                    write_fasta(path, recs)
                    m.write(f"{path}\t{hap}\n")
            out = os.path.join(tmp, "out")
            stage = ["--manifest", manifest, "--outdir", out, "--min-len", "10000", "--k", "15"]
            main(["prepare"] + stage)
            main(["kmers"] + stage)
            # new assembly: HAP2's chr1 and chr2 sequences exchanged under the same names
            recs = haps["HAP2"]
            write_fasta(os.path.join(tmp, "HAP2.fa"),
                        [(recs[0][0], recs[0][1], recs[1][2]), (recs[1][0], recs[1][1], recs[0][2])])
            for f in os.listdir(tmp):
                if f.startswith("HAP2.fa."):
                    os.remove(os.path.join(tmp, f))
            main(["prepare"] + stage)
            main(["kmers"] + stage)
            with open(os.path.join(out, "chroms", "HAP2_chr01.fa")) as f:
                f.readline()
                cached = "".join(line.strip() for line in f)
            from ploidyspec.simulate import to_text
            self.assertEqual(cached, to_text(recs[1][2]))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@unittest.skipUnless(ENABLED, "set PLOIDYSPEC_E2E=1 with samtools + FastK on PATH")
class TestAutoNumberingEndToEnd(unittest.TestCase):
    """Assemblies with no chromosome names at all (ROADMAP 1.7): headers are
    replaced by arbitrary contig names and the sequences shuffled. ploidyspec
    must find the chromosome-scale sequences, number them consistently across
    haplotypes, and reach the same states as with names."""

    SCENARIOS = ("allotetraploid", "autotetraploid")

    @classmethod
    def setUpClass(cls):
        import random
        cls.tmp = tempfile.mkdtemp(prefix="ploidyspec_auto_")
        for name in cls.SCENARIOS:
            simulate(name, os.path.join(cls.tmp, "named"), seed=11,
                     params=dict(SIM_PARAMS, n_chrom=5))
        cls.truth_chrom = {}
        rng = random.Random(3)
        for name in cls.SCENARIOS:
            src = os.path.join(cls.tmp, "named", name)
            dst = os.path.join(cls.tmp, "unnamed", name)
            os.makedirs(dst)
            manifest = []
            for line in open(os.path.join(src, "manifest.tsv")):
                if line.startswith("#"):
                    continue
                fasta, hap = line.rstrip("\n").split("\t")
                fasta = os.path.join(src, fasta) if not os.path.isabs(fasta) else fasta
                recs = open(fasta).read().split(">")[1:]
                rng.shuffle(recs)
                out = os.path.join(dst, f"{hap}.fa")
                with open(out, "w") as f:
                    for i, rec in enumerate(recs):
                        header, body = rec.split("\n", 1)
                        new = f"ctg{rng.randrange(10**6):06d}_{i}"
                        cls.truth_chrom[(name, hap, new)] = int(header.split("chromosome: ")[1])
                        f.write(f">{new}\n{body}")
                manifest.append(f"{out}\t{hap}\n")
            with open(os.path.join(dst, "manifest.tsv"), "w") as f:
                f.writelines(manifest)
            main(["all", "--manifest", os.path.join(dst, "manifest.tsv"),
                  "--outdir", os.path.join(dst, "out")] + RUN_ARGS)

    @classmethod
    def tearDownClass(cls):
        if os.environ.get("PLOIDYSPEC_E2E_KEEP") != "1":
            shutil.rmtree(cls.tmp, ignore_errors=True)

    def units(self, name):
        return read_tsv(os.path.join(self.tmp, "unnamed", name, "out", "sequences.tsv"))

    def test_each_number_holds_one_original_chromosome(self):
        n = 5
        for name in self.SCENARIOS:
            groups = {}
            for u in self.units(name):
                orig = self.truth_chrom[(name, u["hap"], u["seq_id"])]
                # in the allotetraploid, i and i+n are homeologs from different
                # progenitors: they must stay distinct numbers
                groups.setdefault(u["chrom"], set()).add(orig)
            self.assertTrue(all(len(g) == 1 for g in groups.values()), (name, groups))
            self.assertEqual(len(groups), 2 * n if name == "allotetraploid" else n, name)

    def test_every_copy_was_matched(self):
        for name in self.SCENARIOS:
            self.assertTrue(all(int(u["chrom"]) < 1000 for u in self.units(name)), name)

    def test_states_match_named_run(self):
        expected = {"allotetraploid": {"resolved_lineages"}, "autotetraploid": {"tetrasomic_like"}}
        for name in self.SCENARIOS:
            rows = read_tsv(os.path.join(self.tmp, "unnamed", name, "out", "rediploidization",
                                         "rediploidization_by_chrom.tsv"))
            self.assertEqual({r["copy_state"] for r in rows}, expected[name], name)
