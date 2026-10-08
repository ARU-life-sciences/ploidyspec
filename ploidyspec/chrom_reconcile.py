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
from .manifest import AUTO_OFFSET

AMBIGUOUS_RATIO_FLOOR = 1.5
DEFAULT_RATIO_THRESHOLD = 3.0
# Two units that each match the other's declared chrom corroborate each other,
# so a mutual swap needs less per-unit evidence. ddHesMatr1's HAP1 chr01/chr02
# swap (ratios 3.1 and 2.3) sits below the single-unit threshold because its
# allele distances are only 2-4x below its between-chromosome distance.
MUTUAL_SWAP_RATIO = 2.0


def is_auto(u):
    return u.get("numbering") == "auto"


def group_key(u):
    """Units relabelled together: one source file, or for automatically numbered
    units one haplotype (a file can hold several when labels are AUTO)."""
    return f"{u['source']}|{u['hap']}" if is_auto(u) else u["source"]


def choose_reference_source(units):
    """The source file most likely to have internally-consistent chrom numbering:
    whichever contributes the most units (ties broken by earliest first-appearance,
    i.e. manifest order) -- see the module docstring for why this side is never
    itself corrected."""
    auto_ref = {group_key(u) for u in units if is_auto(u) and int(u["chrom"]) < AUTO_OFFSET}
    if auto_ref:
        # automatic numbering: the haplotype prepare numbered 1..n is the reference
        return next(iter(auto_ref))
    counts = Counter(group_key(u) for u in units)
    first_seen = {}
    for i, u in enumerate(units):
        first_seen.setdefault(group_key(u), i)
    return max(counts, key=lambda s: (counts[s], -first_seen[s]))


def mean_distance(i, group_indices, distance):
    """Mean distance from unit i to every unit in group_indices. None if the
    group is empty."""
    if not group_indices:
        return None
    return sum(distance[i][j] for j in group_indices) / len(group_indices)


def match_to_reference(units, idxs, ref_groups, distance):
    """
    One-to-one assignment of the units `idxs` (one source file) to reference
    chrom groups: greedy on mean distance, closest pairs first, with a unit's
    own declared chrom winning ties. Returns {unit index: chrom}; units left
    over when chroms run out are absent.
    """
    pairs = []
    for i in idxs:
        for chrom, members in ref_groups.items():
            d = mean_distance(i, members, distance)
            if d is not None:
                pairs.append((d, chrom != units[i]["chrom"], i, chrom))
    pairs.sort(key=lambda p: (p[0], p[1]))
    assigned, taken = {}, set()
    for _, _, i, chrom in pairs:
        if i not in assigned and chrom not in taken:
            assigned[i] = chrom
            taken.add(chrom)
    return assigned


def detect_relabeling(units, distance, ratio_threshold=DEFAULT_RATIO_THRESHOLD):
    """
    Designates a reference source file (see choose_reference_source), then for
    every unit from a DIFFERENT source file, compares its mean distance to its own
    declared chrom's reference-group members against every other reference
    chrom-group's members. If the best alternative is decisively closer
    (ratio > ratio_threshold), it's a confident mislabeling candidate.

    Each flagged source file's units are then matched one-to-one to the
    reference chrom groups (match_to_reference), so two units can't both claim
    the same chrom: ddLepDrab1's HAP4 numbered 14 of its 16 chromosomes
    differently, and with independent nearest-group choices HAP4_chr07 (0.040
    to chr06, 0.041 to chr03) and HAP4_chr08 (0.028 to chr06) both picked chr06,
    which failed the batch and left the whole haplotype mislabelled. A move is
    applied when the matched chrom is decisively closer than the declared one
    (ratio >= ratio_threshold).

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
        if group_key(u) == ref_source:
            ref_groups[u["chrom"]].append(i)

    candidates = []
    for i, u in enumerate(units):
        if group_key(u) == ref_source:
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

    flagged = {c["index"]: c for c in candidates}
    by_source = defaultdict(list)
    for i, u in enumerate(units):
        if group_key(u) != ref_source:
            by_source[group_key(u)].append(i)

    corrections, ambiguous = [], []
    for src, idxs in by_source.items():
        auto = all(is_auto(units[i]) for i in idxs)
        if not auto and not any(i in flagged for i in idxs):
            continue
        assigned = match_to_reference(units, idxs, ref_groups, distance)
        moves, held = [], []
        for i in idxs:
            own = units[i]["chrom"]
            new = assigned.get(i)
            if new is None or new == own:
                if i in flagged:
                    held.append(flagged[i])
                continue
            own_dist = mean_distance(i, ref_groups.get(own, []), distance)
            new_dist = mean_distance(i, ref_groups[new], distance)
            if own_dist is None:  # provisional auto number: no reference group
                own_dist, ratio = float("nan"), float("inf")
            else:
                ratio = own_dist / new_dist if new_dist > 0 else float("inf")
            move = dict(index=i, unit_id=units[i]["unit_id"], old_chrom=own, new_chrom=new,
                        own_dist=own_dist, alt_dist=new_dist, ratio=ratio)
            # automatically numbered units take their match whatever the ratio
            (moves if auto or ratio >= ratio_threshold else held).append(move)
        # mutual swaps: a held unit whose matched slot belongs to a unit matched
        # into its own slot moves too, if both clear MUTUAL_SWAP_RATIO
        by_old = {c["old_chrom"]: c for c in moves + held if assigned.get(c["index"]) == c["new_chrom"]}
        for c in list(held):
            partner = by_old.get(c["new_chrom"])
            if (assigned.get(c["index"]) == c["new_chrom"] and partner is not None
                    and partner is not c and partner["new_chrom"] == c["old_chrom"]
                    and min(c["ratio"], partner["ratio"]) >= MUTUAL_SWAP_RATIO):
                for x in (c, partner):
                    if x in held:
                        held.remove(x)
                        moves.append(x)
        # a held move below the floor is only the matching's leftover, not evidence
        held = [c for c in held if c["ratio"] >= AMBIGUOUS_RATIO_FLOOR]
        moving_chroms = {c["old_chrom"] for c in moves}
        occupied_chroms = {units[i]["chrom"] for i in idxs}
        closed = all(c["new_chrom"] in moving_chroms or c["new_chrom"] not in occupied_chroms
                     for c in moves)
        if closed:
            corrections.extend(moves)
            ambiguous.extend(held)
        else:
            ambiguous.extend(moves + held)

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
    fieldnames = ["unit_id", "hap", "chrom", "seq_id", "length", "source", "desc", "numbering"]
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

    assigned = [c for c in corrections if is_auto(units[c["index"]])]
    if assigned:
        log(f"automatic chromosome numbering: {len(assigned)} unit(s) given the number of "
            f"their closest reference chromosome (one-to-one k-mer matching)")
        corrections = [c for c in corrections if not is_auto(units[c["index"]])]
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

    moved = {c["index"] for c in assigned}
    unmatched = [u["unit_id"] for i, u in enumerate(units)
                 if is_auto(u) and int(u["chrom"]) >= AUTO_OFFSET and i not in moved]
    if unmatched:
        log(f"automatic chromosome numbering: {len(unmatched)} sequence(s) matched no reference "
            f"chromosome and keep provisional numbers ({', '.join(unmatched[:5])}"
            f"{', ...' if len(unmatched) > 5 else ''}): extra or fragmentary sequence, or a "
            f"chromosome the reference haplotype lacks")
    corrections = corrections + assigned
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
