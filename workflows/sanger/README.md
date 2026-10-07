# Sanger farm workflow

How this project's own panel is run on the Sanger LSF farm. None of this is
needed to use ploidyspec elsewhere: it wraps the `ploidyspec` command with
`module load` and `bsub`, and points at Tree of Life data paths.

- `build_manifests.py` — writes `manifests/*.tsv` and `jobs.tsv` from `meta/meta.tsv`.
- `stage_curated.sh` — re-bgzips read-only curated assemblies into `data/` so
  `samtools faidx` can index them.
- `submit_all.sh` / `run_one.sh` — one LSF job per species running `ploidyspec all`
  (plus optional stages) with the panel's settings.
- `verify_species.sh`, `verify_te_markers.sh` — re-run subsets of stages.
- `rediploidization_status.sh` — progress of a batch submitted with job IDs in
  `logs/redip_panel_jobs.txt` (see ROADMAP.md's progress log).
- `archive/` — one-off batch and repair scripts kept for the record.

Run them from anywhere; each changes to the repository root first.
