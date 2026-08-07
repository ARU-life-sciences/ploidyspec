"""
Cross-source-file chromosome-number reconciliation.

`chrom` numbers are extracted independently per source FASTA file (see
`manifest.py`'s `match_chrom`), purely from regex matches against that file's own
headers -- there is no cross-file consistency check. When a species has multiple
source files (the common case), there is no guarantee that e.g. `SUPER_4` in file A
and `SUPER_4` in file B are the same physical chromosome; they are independently
assigned local numbers that can simply collide. Every downstream stage that groups
units by `chrom` (matrix's ploidy/homology reports, homeologs.py, windowed.py,
te_markers.py, subgenome_report.py) silently trusts this grouping.

This module detects and corrects that using the already-computed whole-genome
distance matrix (chrom-agnostic by construction). One source file is designated the
*reference* (the one contributing the most units -- ties broken by earliest
appearance in the manifest/`units` list); only non-reference files' labels are ever
corrected, always against the reference file's groups. This is a deliberate
asymmetry, not an oversight: scoring every unit against "its own group vs the best
alternative group" independently is unsound when a mismatch is a clean swap between
two groups -- both sides of the swap then look equally wrong to themselves, and
"correcting" both directions at once can merge two genuinely different chromosomes
under one label (caught by a regression test in tests/test_chrom_reconcile.py after
it happened once during development). Picking one side as ground truth and only
ever moving the other side avoids that failure mode entirely, and for the one real
case with a structural majority (a "primary" file vs an "alternate contigs" file
that resolved into several haplotigs per locus), it also happens to pick the file
that's actually self-consistent.

A unit whose distance to its declared same-chrom reference-group members is
decisively worse than its distance to some other reference chrom-group is almost
certainly mislabeled, not genuinely divergent -- confirmed empirically (mismatch
ratios of 4.6x-95.7x across affected species, versus ~1.03x for the one borderline
case that turned out not to be a real mismatch).
"""

import csv
import os
from collections import Counter, defaultdict

from .common import log, matrix_dir

AMBIGUOUS_RATIO_FLOOR = 1.5
DEFAULT_RATIO_THRESHOLD = 3.0


def choose_reference_source(units):
    """The source file most likely to have internally-consistent chrom numbering:
    whichever contributes the most units (ties broken by earliest first-appearance,
    i.e. manifest order) -- see the module docstring for why this side is never
    itself corrected."""
    counts = Counter(u["source"] for u in units)
    first_seen = {}
    for i, u in enumerate(units):
        first_seen.setdefault(u["source"], i)
    return max(counts, key=lambda s: (counts[s], -first_seen[s]))


def mean_distance(i, group_indices, distance):
    """Mean distance from unit i to every unit in group_indices. None if the
    group is empty."""
    if not group_indices:
        return None
    return sum(distance[i][j] for j in group_indices) / len(group_indices)


def detect_relabeling(units, distance, ratio_threshold=DEFAULT_RATIO_THRESHOLD):
    """
    Designates a reference source file (see choose_reference_source), then for
    every unit from a DIFFERENT source file, compares its mean distance to its own
    declared chrom's reference-group members against every other reference
    chrom-group's members. If the best alternative is decisively closer
    (ratio > ratio_threshold), it's a confident mislabeling candidate.

    Candidates are validated as a batch per non-reference source file, not
    accepted independently: a source file's set of qualifying moves is only
    safe to apply together if the moves are injective (no two old chroms
    target the same new chrom) AND every target chrom is either not currently
    used by that source file at all, or is itself being vacated by the same
    batch (a swap or cycle). Accepting a move whose target is occupied by a
    unit that ISN'T also relocating would create two units sharing one
    unit_id -- caught by a regression test after it crashed a real run
    (lpElePalu1: HAP2_chr12 tried to become HAP2_chr11 while the real
    HAP2_chr11 stayed put, a KeyError two steps downstream in
    write_ploidy_and_homology_reports). If a source's batch fails this check,
    every candidate in it is demoted to ambiguous -- no partial application.

    Returns (corrections, ambiguous):
      corrections: list of dicts (index, unit_id, old_chrom, new_chrom, own_dist,
        alt_dist, ratio) for confidently-resolved mismatches.
      ambiguous: same shape, for candidates with AMBIGUOUS_RATIO_FLOOR <= ratio <
        ratio_threshold, or that lost a slot-claiming conflict -- not auto-corrected,
        logged for manual review.
    """
    ref_source = choose_reference_source(units)
    ref_groups = defaultdict(list)
    for i, u in enumerate(units):
        if u["source"] == ref_source:
            ref_groups[u["chrom"]].append(i)

    candidates = []
    for i, u in enumerate(units):
        if u["source"] == ref_source:
            continue
        own_chrom = u["chrom"]
        own_dist = mean_distance(i, ref_groups.get(own_chrom, []), distance)
        if own_dist is None:
            continue
        best_alt_chrom = None
        best_alt_dist = None
        for chrom, idxs in ref_groups.items():
            if chrom == own_chrom:
                continue
            d = mean_distance(i, idxs, distance)
            if d is None:
                continue
            if best_alt_dist is None or d < best_alt_dist:
                best_alt_dist = d
                best_alt_chrom = chrom
        if best_alt_dist is None:
            continue
        ratio = own_dist / best_alt_dist if best_alt_dist > 0 else float("inf")
        if ratio >= AMBIGUOUS_RATIO_FLOOR:
            candidates.append(
                dict(
                    index=i,
                    unit_id=u["unit_id"],
                    old_chrom=own_chrom,
                    new_chrom=best_alt_chrom,
                    own_dist=own_dist,
                    alt_dist=best_alt_dist,
                    ratio=ratio,
                )
            )

    qualifying = [c for c in candidates if c["ratio"] >= ratio_threshold]
    ambiguous = [c for c in candidates if c["ratio"] < ratio_threshold]

    by_source = defaultdict(list)
    for c in qualifying:
        by_source[units[c["index"]]["source"]].append(c)

    corrections = []
    for src, cands in by_source.items():
        moving_chroms = {c["old_chrom"] for c in cands}
        targets = [c["new_chrom"] for c in cands]
        injective = len(set(targets)) == len(targets)
        occupied_chroms = {
            u["chrom"] for u in units if u["source"] == src
        }
        closed = all(
            t in moving_chroms or t not in occupied_chroms for t in targets
        )
        if injective and closed:
            corrections.extend(cands)
        else:
            ambiguous.extend(cands)

    return corrections, ambiguous


def _file_belongs_to_unit(fname, unit_id):
    stripped = fname[1:] if fname.startswith(".") else fname
    return stripped == unit_id or stripped.startswith(unit_id + ".")


def _unit_dirs(outdir):
    dirs = [os.path.join(outdir, "chroms")]
    if os.path.isdir(outdir):
        dirs += [
            os.path.join(outdir, d)
            for d in os.listdir(outdir)
            if d.startswith("ktabs_k")
        ]
    return [d for d in dirs if os.path.isdir(d)]


def _rename_unit_files_batch(outdir, renames):
    """Renames every already-built file (chrom fasta, .ktab/.hist and hidden
    FastK shard files across every already-built k) for a batch of
    (old_unit_id, new_unit_id) corrections. Cheap os.rename, no recomputation --
    this is what makes correction safe to apply after expensive k-mer tables
    already exist.

    Two-phase (old -> unique temp -> new): a swap between two units (chr01<->
    chr02, the common case) means unit A's new name is unit B's old name and
    vice versa -- renaming directly in one pass would have the first rename
    clobber the file the second rename still needs to read."""
    dirs = _unit_dirs(outdir)
    temp_moves = []  # (temp_path, final_path)
    for old_unit_id, new_unit_id in renames:
        for d in dirs:
            for fname in os.listdir(d):
                if _file_belongs_to_unit(fname, old_unit_id):
                    new_fname = fname.replace(old_unit_id, new_unit_id, 1)
                    old_path = os.path.join(d, fname)
                    tmp_path = old_path + ".reconcile_tmp"
                    os.rename(old_path, tmp_path)
                    temp_moves.append((tmp_path, os.path.join(d, new_fname)))
    for tmp_path, final_path in temp_moves:
        os.rename(tmp_path, final_path)


def _rewrite_sequences_tsv(seq_tsv, units):
    fieldnames = ["unit_id", "hap", "chrom", "seq_id", "length", "source", "desc"]
    tmp = seq_tsv + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for u in units:
            w.writerow({k: u.get(k, "") for k in fieldnames})
    os.replace(tmp, seq_tsv)


def _write_corrections_log(outdir, corrections, ambiguous):
    path = os.path.join(matrix_dir(outdir), "chrom_label_corrections.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "status",
                "unit_id",
                "old_chrom",
                "new_chrom",
                "own_group_dist",
                "alt_group_dist",
                "ratio",
            ]
        )
        for c in sorted(corrections, key=lambda c: c["unit_id"]):
            w.writerow(
                [
                    "corrected",
                    c["unit_id"],
                    c["old_chrom"],
                    c["new_chrom"],
                    f"{c['own_dist']:.6f}",
                    f"{c['alt_dist']:.6f}",
                    f"{c['ratio']:.2f}",
                ]
            )
        for c in sorted(ambiguous, key=lambda c: c["unit_id"]):
            w.writerow(
                [
                    "ambiguous",
                    c["unit_id"],
                    c["old_chrom"],
                    c["new_chrom"],
                    f"{c['own_dist']:.6f}",
                    f"{c['alt_dist']:.6f}",
                    f"{c['ratio']:.2f}",
                ]
            )


def reconcile_chrom_labels(
    outdir, seq_tsv, units, distance, ratio_threshold=DEFAULT_RATIO_THRESHOLD
):
    """
    Detects and applies cross-source-file chromosome-number corrections in place:
    mutates `units`' chrom/unit_id for confidently-corrected entries, renames the
    corresponding on-disk .ktab/.hist/chrom-fasta files to match, and rewrites
    `seq_tsv`. Always writes matrix/chrom_label_corrections.tsv (possibly empty)
    for auditability -- corrections are never silent. Returns (corrections,
    ambiguous) as applied/found.
    """
    corrections, ambiguous = detect_relabeling(units, distance, ratio_threshold)

    if corrections:
        log(
            f"chrom-label reconciliation: correcting {len(corrections)} mislabeled "
            f"unit(s) (cross-source-file chromosome-numbering mismatch)"
        )
        for c in corrections:
            log(
                f"  {c['unit_id']}: declared chr{c['old_chrom']} -> chr{c['new_chrom']} "
                f"(own-group dist={c['own_dist']:.4f}, true-group dist={c['alt_dist']:.4f}, "
                f"ratio={c['ratio']:.1f}x)"
            )
    if ambiguous:
        log(
            f"chrom-label reconciliation: {len(ambiguous)} ambiguous case(s) NOT "
            f"auto-corrected (see matrix/chrom_label_corrections.tsv) -- manual "
            f"review recommended"
        )

    renames = []
    for c in corrections:
        i = c["index"]
        old_unit_id = units[i]["unit_id"]
        new_chrom = c["new_chrom"]
        new_unit_id = f"{units[i]['hap']}_chr{int(new_chrom):02d}"
        renames.append((old_unit_id, new_unit_id))
        units[i]["chrom"] = new_chrom
        units[i]["unit_id"] = new_unit_id
    if renames:
        _rename_unit_files_batch(outdir, renames)

    _write_corrections_log(outdir, corrections, ambiguous)
    if corrections:
        _rewrite_sequences_tsv(seq_tsv, units)

    return corrections, ambiguous
