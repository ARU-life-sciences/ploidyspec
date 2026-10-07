# Retired panel tables

Last outputs of metrics retired on 2026-10-07 (ROADMAP Phase 2, option B),
kept because INTERPRETATION.md quotes them:

- `auto_allo_spectrum.tsv`: partition_consistency, run length / flip_rate,
  windowed CVs, pair_depth_cv, distance_ratio_cv, combined_allo_score.
- `poly_space_*`: the PCA over those metrics.

Several values are known to be wrong: `ddLepDrab1` (uncorrected HAP4
chromosome numbering), `ddHesMatr1` (HAP1 chr01/chr02 swap), `drLytSali1`
(previous assembly), and every windowed-derived column (equal-coordinate
windows). See INTERPRETATION.md, "Phase 2 hand-check: corrections". The
current tables are `meta/panel_summary.tsv` and `meta/supplementary/`.
