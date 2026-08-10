import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.report import RANKED_CANDIDATES_LIMIT, render_report_html


def empty_data(species="testSpecies1"):
    return {
        "species": species,
        "matrix": {
            "ploidy_summary": None,
            "homologous_chromosomes": None,
            "heatmap": None,
            "heatmap_contrast": None,
            "k_resolution": None,
        },
        "windowed": {"overview": None},
        "homeologs": {
            "homeolog_pairs": None,
            "ploidy_ancestry_summary": None,
            "ranked_candidates": None,
            "homeolog_pairs_plot": None,
        },
        "subgenomes": None,
    }


class TestRenderReportHtml(unittest.TestCase):
    def test_fully_populated_species_has_all_sections(self):
        data = empty_data()
        data["matrix"]["ploidy_summary"] = [
            {"chrom": "chr01", "n_haplotype_copies": "2", "haplotype_labels": "HAP1,HAP2"}
        ]
        data["matrix"]["heatmap"] = "data:image/png;base64,AAAA"
        data["windowed"]["overview"] = "data:image/png;base64,BBBB"
        data["homeologs"]["homeolog_pairs"] = [
            {"chrom_a": "chr01", "chrom_b": "chr02", "mean_distance": "0.05"}
        ]
        data["subgenomes"] = {
            "auto_allo_index": [{"chrom": "chr01", "te_marker_fraction": "0.22"}],
            "windows_summary": [{"unit": "HAP1_chr01", "pct_self": "90.0"}],
            "windowed_plots": ["data:image/png;base64,CCCC"],
        }

        out = render_report_html(data)
        self.assertIn("<html>", out)
        self.assertIn("testSpecies1", out)
        self.assertIn("Whole-chromosome matrix", out)
        self.assertIn("data:image/png;base64,AAAA", out)
        self.assertIn("Windowed divergence", out)
        self.assertIn("data:image/png;base64,BBBB", out)
        self.assertIn("Ancient homeolog pairing", out)
        self.assertIn("Subgenomes / auto-allo index", out)
        self.assertIn("0.22", out)
        self.assertIn("data:image/png;base64,CCCC", out)
        self.assertNotIn("Not run yet", out)
        self.assertNotIn("te-markers not run", out)

    def test_partial_species_degrades_gracefully(self):
        # mirrors most of the panel before te-markers was ever run: only
        # matrix/ has real data, everything else is absent
        data = empty_data()
        data["matrix"]["ploidy_summary"] = [
            {"chrom": "chr01", "n_haplotype_copies": "2", "haplotype_labels": "HAP1,HAP2"}
        ]

        out = render_report_html(data)
        # no exception, valid-looking document, and the missing stages say so
        self.assertIn("<html>", out)
        self.assertIn("Whole-chromosome matrix", out)
        self.assertIn("chr01", out)  # the populated ploidy_summary table rendered
        self.assertIn("Windowed divergence", out)
        self.assertIn("Not run yet", out)
        self.assertIn("Ancient homeolog pairing", out)
        self.assertIn("Subgenomes / auto-allo index", out)
        self.assertIn("te-markers not run", out)

    def test_completely_empty_species_does_not_crash(self):
        # prepare-only, nothing else has run at all
        out = render_report_html(empty_data())
        self.assertIn("<html>", out)
        self.assertIn("Not run yet", out)
        self.assertIn("te-markers not run", out)

    def test_ranked_candidates_truncated_with_note(self):
        data = empty_data()
        data["homeologs"]["ranked_candidates"] = [
            {
                "chrom_a": f"chr{i:02d}",
                "chrom_b": f"chr{i+1:02d}",
                "distance": f"{0.01 * i:.4f}",
                "z_score": "-5.0",
                "p_value": "0.001",
                "q_value": "0.01",
                "accepted": "False",
            }
            for i in range(30)
        ]
        out = render_report_html(data)
        self.assertIn(f"showing top {RANKED_CANDIDATES_LIMIT} of 30 rows", out)
        # only the first RANKED_CANDIDATES_LIMIT rows' chrom_a values appear as table cells
        self.assertIn("chr00", out)
        self.assertIn(f"chr{RANKED_CANDIDATES_LIMIT - 1:02d}", out)
        self.assertNotIn(f"chr{RANKED_CANDIDATES_LIMIT + 1:02d}", out)


if __name__ == "__main__":
    unittest.main()
