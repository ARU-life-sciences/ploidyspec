"""
Plain-language answers to the four questions for one species, read from the
outputs of the other stages (no new k-mer work):

1. Ploidy: haplotype copies per chromosome number, then whether the
   chromosome numbers themselves form duplicated sets (homeolog groups), and
   how deep those sets are.
2. Origin-like structure: are copies interchangeable (tetrasomic-like) or in
   separate lineages (disomic-like), and how consistently across chromosomes.
3. TE markers: does repeat history track the lineages (within-genome contrast).
4. Rediploidization: fusions, chromosomes split into lineages among
   tetrasomic-like ones, regional splits, residual tetrasomy between
   homeologs, asynchronous ancient pairs.

Each answer has a confidence (high / medium / low / not assessable) and the
evidence it rests on. The rules follow docs/guide (Chapter 6). They describe
structure in one individual's assembly: not inheritance, not origin.

Writes results/<species>/summary.tsv and summary.md.
"""

import csv
import os
import statistics
from collections import Counter

from .common import log

FIELDS = ["species", "question", "answer", "confidence", "evidence"]
SPLIT_STATES = ("resolved_lineages", "fusion_lineages", "partially_resolved")

# A homeolog group partition is trusted as "sets per haplotype" at this z, and
# only when its groups are small (homeolog groups, not subgenomes).
SET_PARTITION_Z = 10.0
# Ancient sets are flagged as possibly paleopolyploid when their distance is
# this close to the unrelated-chromosome background.
PALEO_RATIO = 0.8
# Residual tetrasomy on this many paired chromosome numbers reads as an
# auto-like origin (one could be a homeologous exchange).
RESIDUAL_MIN_CHROMS = 2


def _read(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def gather(outdir):
    """Everything the four answers need, from the stage outputs."""
    from .panel import copy_divergence, modal_copies, te_markers_row  # avoid import cycle

    species = os.path.basename(os.path.normpath(outdir))
    copies, n_numbers = modal_copies(outdir)
    allele, cross = copy_divergence(outdir)
    pairs = _read(os.path.join(outdir, "homeologs", "homeolog_pairs.tsv")) or []
    depths = [float(r["mean_distance"]) for r in pairs if r.get("mean_distance")]
    sync = (_read(os.path.join(outdir, "structure", "pair_synchrony.tsv")) or [{}])[0]
    parts = _read(os.path.join(outdir, "structure", "genome_partition.tsv")) or []
    summary = {r["metric"]: r["value"] for r in
               _read(os.path.join(outdir, "rediploidization", "rediploidization_summary.tsv")) or []}
    by_chrom = _read(os.path.join(outdir, "rediploidization", "rediploidization_by_chrom.tsv")) or []
    fusions = [r for r in _read(os.path.join(outdir, "rediploidization", "fusions.tsv")) or []
               if r.get("status") == "fusion"]
    te = te_markers_row(species, outdir)
    lineage = _read(os.path.join(outdir, "subgenomes", "te_marker_fraction_by_lineage.tsv")) or []
    return dict(
        species=species,
        copies=int(copies) if str(copies).isdigit() else 0,
        n_numbers=n_numbers,
        allele=allele, cross=cross,
        n_pairs=len(pairs),
        pair_depth=statistics.median(depths) if depths else None,
        pair_cv=_num(sync.get("pair_depth_cv")),
        parts=parts,
        summary=summary,
        by_chrom=by_chrom,
        fusions=fusions,
        te=te,
        te_splits=[_num(r["split_ratio"]) for r in lineage if _num(r.get("split_ratio")) is not None],
    )


def duplicated_sets(g):
    """(sets per haplotype, base number x, basis, confidence)."""
    n = g["n_numbers"] or 0
    small = [r for r in g["parts"]
             if int(r["k"]) >= n / 4 and _num(r["z_score"]) and _num(r["z_score"]) >= SET_PARTITION_Z]
    if small:
        best = max(small, key=lambda r: _num(r["z_score"]))
        sizes = [int(s) for s in best["group_sizes"].split(",")]
        modal, count = Counter(sizes).most_common(1)[0]
        if modal >= 2:
            uniform = count / len(sizes)
            conf = "high" if uniform >= 0.8 and _num(best["z_score"]) >= 15 else "medium"
            return modal, int(best["k"]), (f"{best['k']} groups of mostly {modal} chromosome numbers "
                                           f"(genome partition z = {_num(best['z_score']):.1f})"), conf
    if n and 2 * g["n_pairs"] / n >= 0.5:
        frac = 2 * g["n_pairs"] / n
        return 2, round(n / 2), f"{g['n_pairs']} homeolog pairs covering {frac:.0%} of numbers", \
            "high" if frac >= 0.8 else "medium"
    k2 = [r for r in g["parts"] if r["k"] == "2" and _num(r["z_score"]) and _num(r["z_score"]) >= SET_PARTITION_Z]
    if k2:
        sizes = [int(s) for s in k2[0]["group_sizes"].split(",")]
        if min(sizes) >= n / 4:
            return 2, round(n / 2), (f"no accepted pairs, but a balanced two-group partition "
                       f"({k2[0]['group_sizes'].replace(',', ' vs ')}, z = {_num(k2[0]['z_score']):.1f})"), "medium"
    return 1, n, "no duplicated sets of chromosome numbers detected", "medium"


def answer_ploidy(g):
    c, n = g["copies"], g["n_numbers"]
    sets, x, basis, conf = duplicated_sets(g)
    if c <= 1:
        return ("Not assessable from this assembly: only one usable haplotype "
                f"({n} chromosome numbers)."
                + (f" The numbers do form duplicated sets ({basis})." if sets > 1 else ""),
                "not assessable", basis)
    text = f"{c} haplotype copies of each of {n} chromosome numbers"
    evidence = [f"copies from the assembly ({c} per number)"]
    if sets > 1:
        text += (f"; the numbers form sets of {sets} ({basis.split(' (')[0]}), so x ≈ {x} "
                 f"and up to {c * sets}x relative to x")
        evidence.append(basis)
        if g["pair_depth"] and g["cross"]:
            ratio = g["pair_depth"] / g["cross"]
            evidence.append(f"homeolog distance {g['pair_depth']:.3f} vs unrelated {g['cross']:.3f}")
            if ratio >= PALEO_RATIO:
                text += "; the sets are nearly as diverged as unrelated chromosomes, so they may be an old (paleo) duplication"
    else:
        text += "; no older duplicated sets detected"
    text += "."
    # copy count is set by the assembly; the set count is inferred
    confidence = conf if sets > 1 else "high"
    return text, confidence, "; ".join(evidence)


def residual(g):
    """(tested, with residual tetrasomy, bp, terminal fraction, evidence, confidence) or None
    when windowed-homeologs did not run."""
    s = g["summary"]
    if "residual_tested_chromosomes" not in s:
        return None
    tested, hit, bp = int(s["residual_tested_chromosomes"]), int(s["residual_chromosomes"]), int(s["residual_bp"])
    term = _num(s.get("residual_terminal_frac"))
    controlled = s.get("residual_controlled") == "yes"
    ev = (f"residual tetrasomy on {hit}/{tested} homeolog-paired chromosome numbers"
          + (f", {bp / 1e6:.1f} Mb, {term:.0%} of it terminal" if hit else "")
          + ("" if controlled else " (no unrelated-chromosome control; shared repeats not excluded)"))
    return tested, hit, bp, term, ev, "medium" if controlled else "low"


def _state_counts(g):
    states = Counter(r["copy_state"] for r in g["by_chrom"])
    bases = Counter(r.get("state_basis", "") for r in g["by_chrom"])
    return states, bases


def odd_haplotype(g, onediv):
    """The haplotype that is the odd copy on >= 75% of one-copy-apart chromosomes, or None."""
    hap = g["summary"].get("most_frequent_outlier_hap", "")
    if "(" not in hap or not onediv:
        return None
    return hap if int(hap.split("(")[1].split()[0]) >= 0.75 * onediv else None


def one_copy_apart(g, onediv, ev):
    hap = odd_haplotype(g, onediv)
    if hap:
        return (f"One copy apart on most chromosomes, nearly always the same haplotype ({hap}): "
                "most likely an assembly or phasing problem in that haplotype.", "medium", ev)
    return ("One copy apart on most chromosomes, changing haplotype between chromosomes: "
            "AAAB-like (one divergent genome copy).", "medium", ev)


def answer_structure(g):
    c = g["copies"]
    states, bases = _state_counts(g)
    n = sum(states.values())
    assessable = n - states.get("not_assessable", 0)
    if c <= 1 or assessable == 0:
        if c == 2 and g["n_pairs"] == 0 and not states.get("resolved_lineages"):
            return ("Not assessable: two haplotypes and no duplicated sets to compare. "
                    "Diploid-like as assembled.", "not assessable", "no homeolog pairs, no genome partition")
        return ("Not assessable from this assembly.", "not assessable", "fewer than three copies and nothing to pool")
    split = sum(states.get(s, 0) for s in SPLIT_STATES)
    tet = states.get("tetrasomic_like", 0)
    onediv = states.get("one_divergent_copy", 0)
    frac = lambda k: k / assessable  # noqa: E731
    ev = (f"{assessable} chromosome numbers assessed: {split} split into lineages, {tet} tetrasomic-like, "
          f"{onediv} one divergent copy, {states.get('candidate', 0)} candidate")
    pooled = bases.get("homeolog_pool", 0) + bases.get("partition_pool", 0)
    if c == 2 or pooled >= assessable / 2:
        # two-haplotype assemblies are read through pooling (guide 6.2): resolved
        # only shows the sets do not interchange; tetrasomic-like is the informative case
        via_partition = bases.get("partition_pool", 0) > bases.get("homeolog_pool", 0)
        if onediv >= assessable / 2:
            return one_copy_apart(g, onediv, ev)
        if frac(tet) >= 0.5:
            return ("Auto-like: the duplicated sets are about as close as alleles, so the copies "
                    "look interchangeable across chromosome numbers.", "medium", ev)
        res = residual(g)
        if res and res[1] >= RESIDUAL_MIN_CHROMS:
            return ("Auto-like origin, mostly rediploidized: the duplicated sets have separated along "
                    f"most of their length, but {res[1]} of {res[0]} paired chromosome numbers keep "
                    "stretches where the homeologs are as close as alleles (residual tetrasomy). "
                    "That is expected after an autopolyploidization; recent homeologous exchanges "
                    "in an allopolyploid can leave similar stretches.", res[5], f"{ev}; {res[4]}")
        if frac(split) >= 0.6:
            conf = "medium" if via_partition else ("high" if (g["pair_cv"] or 1) <= 0.06 else "medium")
            why = ("a genome-wide partition (no accepted pairs)" if via_partition
                   else f"homeolog pairs (pair-depth CV {g['pair_cv']:.2f})" if g["pair_cv"] is not None
                   else "homeolog pairs")
            return ("Allo-like (disomic-like): the duplicated sets do not interchange, from "
                    f"{why}. This is also what a long-rediploidized autopolyploid looks like.", conf, ev)
        return ("Mixed: some duplicated sets look interchangeable, others separate.", "medium", ev)
    if onediv >= assessable / 2:
        return one_copy_apart(g, onediv, ev)
    if frac(split) >= 0.6:
        return ("Allo-like (disomic-like): most chromosomes split into two lineages. A long-"
                "rediploidized autopolyploid would look the same.",
                "high" if frac(split) >= 0.8 else "medium", ev)
    if frac(tet) >= 0.6 and not g["fusions"]:
        return ("Auto-like (tetrasomic-like): on most chromosomes all copies are about equally close.",
                "high" if frac(tet) >= 0.8 else "medium", ev)
    return ("Mixed: tetrasomic-like and split chromosomes in one genome (see rediploidization).",
            "medium", ev)


def answer_te(g):
    if not g["te"]:
        return "TE markers not run for this species.", "not assessable", ""
    abs_frac = g["te"].get("te_marker_fraction_mean") or ""
    if g["te_splits"]:
        strong = sum(1 for s in g["te_splits"] if s >= 2)
        med = statistics.median(g["te_splits"])
        ev = (f"TE split >= 2 on {strong}/{len(g['te_splits'])} chromosomes (median {med:.2f}); "
              f"absolute te_marker_fraction {abs_frac} (reference only)")
        if strong >= max(1, len(g["te_splits"]) // 4):
            onediv = sum(1 for r in g["by_chrom"] if r["copy_state"] == "one_divergent_copy")
            if odd_haplotype(g, onediv):
                return ("Repeat content differs between the split copies on several chromosomes, "
                        "but the split isolates the same odd haplotype throughout, so this may "
                        "reflect that haplotype's assembly rather than separate TE histories.",
                        "low", ev)
            return ("Repeat history tracks the lineages on several chromosomes: separate TE "
                    "histories in the split copies.", "medium", ev)
        if strong:
            return (f"Repeat history tracks the lineages on {strong} chromosome(s) only.", "low", ev)
        return ("Repeat history does not track the copies' lineages (no within-genome TE contrast).",
                "medium", ev)
    return ("Only the absolute marker fraction is available (two copies per number); it does not "
            "separate diploid, auto- and allopolyploid, so no reading is made from it.",
            "not assessable", f"te_marker_fraction mean {abs_frac}")


def answer_rediploidization(g):
    states, bases = _state_counts(g)
    c = g["copies"]
    fus = sorted({f["components"] for f in g["fusions"]})
    split = sum(states.get(s, 0) for s in SPLIT_STATES)
    tet = states.get("tetrasomic_like", 0)
    regional = sum(1 for r in g["by_chrom"] if r.get("split_extent") == "regional"
                   and r["copy_state"] in ("partially_resolved", "candidate"))
    ev = []
    if fus:
        ev.append(f"{len(fus)} fusion(s): {', '.join(fus)}")
    if c >= 3:
        ev.append(f"{split} split / {tet} tetrasomic-like chromosomes; {regional} regional split(s)")
    res = residual(g)
    if res:
        ev.append(res[4])
    if g["pair_cv"] is not None:
        ev.append(f"pair-depth CV {g['pair_cv']:.2f}")
    ev = "; ".join(ev)
    if fus:
        return (f"Yes: {len(fus)} chromosome fusion(s) between copies, which separate fused and "
                "unfused lineages (the snow carp mechanism).", "high", ev)
    if c >= 3 and split and tet:
        if split >= max(2, 0.1 * (split + tet)):
            return ("Likely under way: some chromosomes have split into lineages while others stay "
                    "tetrasomic-like.", "medium", ev)
        return (f"Little sign: {split} chromosome(s) partly split among tetrasomic-like ones.",
                "low", ev)
    if c >= 3 and tet and not split:
        return ("No sign yet: copies interchangeable throughout.", "medium", ev)
    if c >= 3 and split and not tet:
        return ("Split throughout: either rediploidization is complete or the genome was "
                "allo-like from the start; these data cannot tell which.", "medium", ev)
    if c < 3 and res and res[1]:
        where = f", {res[3]:.0%} of it at chromosome ends" if res[3] is not None else ""
        return (f"Partly: {res[1]} of {res[0]} homeolog-paired chromosome numbers keep stretches of "
                f"residual tetrasomy ({res[2] / 1e6:.1f} Mb{where}) where the homeologs are still as "
                "close as alleles; elsewhere they have diverged.", res[5], ev)
    if c < 3 and res:
        return ("Diverged throughout the paired chromosomes: no stretch where homeologs are as close "
                "as alleles. Rediploidization is complete, or the genome was allo-like from the "
                "start.", res[5], ev)
    if g["pair_cv"] is not None and g["pair_cv"] >= 0.15:
        return ("Possibly: ancient pairs diverged to different depths (asynchronous resolution).",
                "low", ev)
    if c == 2 and g["n_pairs"]:
        return ("Not tested along the chromosomes (run --with-windowed-homeologs for the residual-"
                "tetrasomy test); ancient pairs are synchronous.", "low", ev)
    return ("Not assessable from this assembly.", "not assessable", ev)


QUESTIONS = [
    ("ploidy", answer_ploidy),
    ("origin_like_structure", answer_structure),
    ("te_markers", answer_te),
    ("rediploidization", answer_rediploidization),
]


def summarize_species(outdir):
    g = gather(outdir)
    rows = []
    for q, fn in QUESTIONS:
        text, conf, ev = fn(g)
        rows.append(dict(species=g["species"], question=q, answer=text, confidence=conf, evidence=ev))
    return rows


def to_markdown(rows):
    species = rows[0]["species"] if rows else ""
    titles = dict(ploidy="Ploidy", origin_like_structure="Auto-like or allo-like?",
                  te_markers="TE markers", rediploidization="Rediploidization")
    out = [f"# {species}", ""]
    for r in rows:
        out += [f"## {titles[r['question']]}", "", r["answer"], "",
                f"*Confidence:* {r['confidence']}" + (f"  \n*Evidence:* {r['evidence']}" if r["evidence"] else ""), ""]
    return "\n".join(out)


def compute_summary(outdir):
    rows = summarize_species(outdir)
    with open(os.path.join(outdir, "summary.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(outdir, "summary.md"), "w") as f:
        f.write(to_markdown(rows) + "\n")
    log(f"wrote summary.tsv and summary.md in {outdir}")
    return rows
