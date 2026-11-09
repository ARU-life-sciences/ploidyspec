#!/usr/bin/env python3
"""Weekly scan report: read each screened release's k-mer outputs, update
scan/state.tsv, and write scan/candidates.tsv (every screened release) and
scan/reports/<date>.md (this run's releases, then all candidates so far).

Readings, from the panel's own readers (ploidyspec.panel):
- copies per chromosome number (matrix): more copies than haplotype files
  means each haplotype holds more than one set (daBudDavi1, wheat);
- homeolog pairs (structure, FDR-accepted) and their median distance;
- the subgenome partition (homeologs/subgenome_partition.tsv).

Verdicts, calibrated on the panel:
- `assembled polyploid`: >= 3 haplotypes, or more copies per number than files;
- `WGD signal`: >= 3 homeolog pairs, covering at least PAIR_FRAC of the
  chromosome numbers or with an accepted subgenome partition (a partition
  with no pairs behind it, bSomMol2, is not enough). A whole-genome
  duplication pairs most numbers (0.5 when every number has one partner; wheat
  0.33-0.5); the birds and Lepidoptera of the 9 November run paired 0.03-0.23,
  from fissions, microchromosomes and shared repeats. The pairs' median
  distance places it: young < 0.075 (confirmed allopolyploids
  0.05-0.07, Ranunculus 0.071), intermediate 0.075-0.11 (salmonid 4R 0.098-0.105,
  Potamogeton 0.077-0.098), old > 0.11;
- `weak signal`: >= 3 pairs, fewer than PAIR_FRAC of the numbers, no partition;
- `metadata only`: a metadata flag (weekly.flags) without k-mer support;
- `no signal`. One or two pairs is background: gar, bowfin and the diploid
  molluscs each show one.
Usage: report.py DATE
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from ploidyspec.panel import modal_copies, species_structure  # noqa: E402
from weekly import SCAN, STATE_FIELDS, flags, read, write  # noqa: E402

MIN_PAIRS = 3
PAIR_FRAC = 0.25
YOUNG, OLD = 0.075, 0.11
FIELDS = ["tolid", "species", "verdict", "age", "n_haps", "copies_per_chromosome", "n_chromosome_numbers",
          "n_pairs", "pair_fraction", "pair_depth_median", "subgenome_partition", "subgenome_z", "metadata_flags",
          "screen_mode", "screened", "status"]
ORDER = ["assembled polyploid", "WGD signal", "weak signal", "metadata only", "no signal", "failed"]


def kv(path):
    rows = read(path)
    return {r["metric"]: r["value"] for r in rows if "metric" in r}


def reading(tolid, cat, meta, st):
    d = os.path.join(SCAN, "results", tolid)
    n_haps = int(cat.get("n_haps") or 2)
    row = dict(tolid=tolid, species=cat.get("species", ""), n_haps=n_haps, screened=st["date"],
               status=st["status"], metadata_flags="; ".join(flags(meta)))
    if st["status"] != "done":
        row["verdict"] = "failed" if st["status"] == "failed" else "pending"
        return row
    mode_path = os.path.join(d, "screen.mode")
    mode = open(mode_path).read().strip() if os.path.exists(mode_path) else ""
    row["screen_mode"] = mode
    if mode == "hap1":
        n_haps = 1  # hap2 was not chromosome-scale: copies are read from hap1 alone
    copies, n_chrom = modal_copies(d)
    sync, _ = species_structure(tolid, d)
    sub = kv(os.path.join(d, "homeologs", "subgenome_partition.tsv"))
    n_pairs = int(sync.get("n_accepted_pairs") or 0)
    depth = sync.get("pair_depth_median") or ""
    accepted = sub.get("accepted") == "yes"
    frac = n_pairs / n_chrom if n_chrom else 0.0
    row.update(copies_per_chromosome=copies, n_chromosome_numbers=n_chrom, n_pairs=n_pairs,
               pair_fraction=f"{frac:.2f}",
               pair_depth_median=depth, subgenome_partition="yes" if accepted else "no",
               subgenome_z=sub.get("z_score", ""))
    if n_haps > 2 or (copies and int(copies) > n_haps):
        row["verdict"] = "assembled polyploid"
    elif n_pairs >= MIN_PAIRS and (frac >= PAIR_FRAC or accepted):
        row["verdict"] = "WGD signal"
    elif n_pairs >= MIN_PAIRS:
        row["verdict"] = "weak signal"
    elif row["metadata_flags"]:
        row["verdict"] = "metadata only"
    else:
        row["verdict"] = "no signal"
    if depth and n_pairs >= MIN_PAIRS:
        x = float(depth)
        row["age"] = "young" if x < YOUNG else "old" if x > OLD else "intermediate"
    return row


def table(rows):
    cols = ["tolid", "species", "verdict", "age", "n_haps", "screen_mode", "copies_per_chromosome", "n_pairs",
            "pair_fraction", "pair_depth_median", "subgenome_partition", "metadata_flags"]
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return "\n".join(out)


def main():
    date = sys.argv[1]
    catalog = {r["tolid"]: r for r in read(os.path.join(SCAN, "catalog.tsv"))}
    meta = {r["tolid"]: r for r in read(os.path.join(SCAN, "metadata.tsv"))}
    state_path = os.path.join(SCAN, "state.tsv")
    state = read(state_path)
    for st in state:
        if st["status"] == "submitted":
            d = os.path.join(SCAN, "results", st["tolid"])
            if os.path.exists(os.path.join(d, "screen.done")):
                st["status"] = "done"
            elif os.path.exists(os.path.join(d, "screen.failed")):
                st["status"] = "failed"
    write(state_path, state, STATE_FIELDS)

    rows = [reading(st["tolid"], catalog.get(st["tolid"], {}), meta.get(st["tolid"], {}), st) for st in state]
    rank = {v: i for i, v in enumerate(ORDER)}
    key = lambda r: (rank.get(r.get("verdict"), 9), -(r.get("n_pairs") or 0))  # noqa: E731
    rows.sort(key=key)
    write(os.path.join(SCAN, "candidates.tsv"), rows, FIELDS)

    this = [r for r in rows if r["screened"] == date]
    hits = [r for r in rows if r.get("verdict") in ("assembled polyploid", "WGD signal")]
    unscreened = [m for t, m in meta.items() if flags(m) and t not in {s["tolid"] for s in state}]
    counts = {v: sum(1 for r in this if r.get("verdict") == v) for v in ORDER + ["pending"]}
    lines = [f"# Polyploidy scan, {date}", "",
             f"{len(catalog)} ToL releases with two or more haplotypes; {len(state)} screened so far, "
             f"{len(this)} this run (" + ", ".join(f"{n} {v}" for v, n in counts.items() if n) + ").", "",
             "Verdicts are leads, not calls: run the full pipeline before reading anything into one. "
             "One or two homeolog pairs is background (gar, bowfin and diploid molluscs each show one); "
             f"so is a few pairs among many chromosome numbers (`weak signal`: < {PAIR_FRAC} of numbers paired, "
             "as in birds and Lepidoptera).", "",
             "## This run", "", table(this) if this else "Nothing screened.", "",
             f"## All k-mer candidates so far ({len(hits)})", "", table(hits) if hits else "None yet.", "",
             f"## Metadata flags not yet screened ({len(unscreened)})", ""]
    lines += [f"- {m['tolid']} ({m['species']}): {'; '.join(flags(m))}" for m in unscreened[:50]]
    if len(unscreened) > 50:
        lines.append(f"- ... and {len(unscreened) - 50} more (scan/metadata.tsv)")
    path = os.path.join(SCAN, "reports", f"{date}.md")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(path)


if __name__ == "__main__":
    main()
