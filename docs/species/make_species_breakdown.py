"""
Write docs/species/species_breakdown.md: a detailed per-species breakdown of
every result directory -- assembly, copies and duplicated sets, each
rediploidization signal, the four summary answers and caveats -- plus a panel
table. Run from the repository root after the species have been (re)run:

    python3 docs/species/make_species_breakdown.py [--results results] [--json out.json]

Reads only result files and meta/literature_ploidy.tsv; nothing is computed
beyond counting and formatting. `--json` also writes the same content as data
(for an HTML view).
"""

import argparse
import csv
import json
import os
import statistics
import sys
from collections import Counter, OrderedDict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from ploidyspec.panel import copy_divergence  # noqa: E402

OUT = "docs/species/species_breakdown.md"

# species with no meta/literature_ploidy.tsv row: (name, what is known, source)
EXTRA = {
    "SchCurv1": ("Schizothorax curvilabiatus", "4n = 98, autotetraploid snow carp", "Xie et al. 2026"),
    "SchYoun1": ("Schizopygopsis younghusbandi", "4n = 90, autotetraploid snow carp", "Xie et al. 2026"),
    "OncMaso1": ("Oncorhynchus masou", "salmonid whole-genome duplication (Ss4R, ~90 Mya), autotetraploid "
                 "origin; residual tetrasomy mapped in Atlantic salmon", "Lien et al. 2016"),
    "fSalAlp3": ("Salvelinus alpinus", "salmonid whole-genome duplication (Ss4R)", "Lien et al. 2016"),
    "fCorLav1": ("Coregonus lavaretus", "salmonid whole-genome duplication (Ss4R)", "Lien et al. 2016"),
    "fLepOcu1": ("Lepisosteus oculatus", "no duplication since vertebrate 2R; outgroup to the teleost "
                 "duplication (negative control)", "Braasch et al. 2016"),
    "fAmiCal2": ("Amia calva", "no duplication since vertebrate 2R; outgroup to the teleost "
                 "duplication (negative control)", "Thompson et al. 2021"),
    "icStrMela3": ("Strophosoma melanogrammum", "triploid (AAB) by the ToL k-mer ploidy plot",
                   "ToL genomic_data k31 ploidy plot"),
    "xgCepNemo3": ("Cepaea nemoralis", "selected as a WGD candidate (land snail)", "candidate list"),
    "xgHygCinc1": ("Hygromia cinctella", "selected as a WGD candidate (land snail)", "candidate list"),
    "xgMonCant1": ("Monacha cantiana", "selected as a WGD candidate (land snail)", "candidate list"),
    "xgDauRufa1": ("Daudebardia rufa", "selected as a WGD candidate (land snail)", "candidate list"),
    "xgPomEleg1": ("Pomatias elegans", "selected as a WGD candidate (land snail)", "candidate list"),
    "xgStaPalu1": ("Stagnicola palustris", "selected as a WGD candidate (freshwater snail)", "candidate list"),
    "xgLitLitt3": ("Littorina littorea", "selected as a WGD candidate (periwinkle)", "candidate list"),
    "lpTriTurg1_A": ("Triticum turgidum, A subgenome", "2n = 4x = 28 (AABB); A run alone", "known"),
    "lpTriTurg1_B": ("Triticum turgidum, B subgenome", "2n = 4x = 28 (AABB); B run alone", "known"),
    "lpTriTurg1_AB": ("Triticum turgidum", "2n = 4x = 28 (AABB), allotetraploid", "known"),
}

# known issues that the result files cannot show by themselves
NOTES = {
    "OncMaso1": "HAP2 was contig-level and was scaffolded on HAP1 (workflows/prep/scaffold_by_reference.py); "
                "positions along HAP2 are borrowed from HAP1, and exchange candidates can be scaffolding "
                "placements.",
    "xgMonCant1": "HAP1 carries chr23+chr24 fused, HAP2 has them apart: a fusion polymorphism or a "
                  "scaffolding join; HAP2 was numbered automatically.",
    "lpElePalu1": "Holocentric chromosomes (Eleocharis): fission and fusion are frequent, so chromosome "
                  "numbers need not reflect ancestral units.",
    "icStrMela3": "HAP1 holds two chromosome sets (A as chr1-10, B as chr11-20); the matrix stage splits B "
                  "out as HAP1B. HAP1 chr09 sits between A and B (0.030 to HAP2 chr09, 0.027 to B).",
    "fCorLav1": "HAP2's NCBI chromosome numbers differ from HAP1's on 30 of 40 chromosomes (cycles); "
                "renumbered by the matrix stage. HAP2 chr33 has no allele in HAP1.",
}

QUESTIONS = OrderedDict(ploidy="Ploidy", origin_like_structure="Auto- or allo-like",
                        te_markers="TE markers", rediploidization="Rediploidization")
STATES = ("tetrasomic_like", "candidate", "partially_resolved", "resolved_lineages", "fusion_lineages",
          "one_divergent_copy", "not_assessable")


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def metrics(path):
    return {r["metric"]: r["value"] for r in read_tsv(path)}


def literature():
    out = {}
    for r in read_tsv("meta/literature_ploidy.tsv"):
        bits = [x for x in (r.get("literature_ploidy_level"), r.get("literature_chrom_number"),
                            r.get("literature_origin") and f"origin: {r['literature_origin']}") if x]
        out[r["species"]] = (r.get("scientific_name", ""), "; ".join(bits),
                             r.get("key_citations", ""))
    out.update(EXTRA)
    return out


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mb(bp):
    v = num(bp)
    return f"{v / 1e6:.1f} Mb" if v is not None else ""


def species_record(sp, d, lit):
    seqs = read_tsv(os.path.join(d, "sequences.tsv"))
    unplaced = read_tsv(os.path.join(d, "unplaced.tsv"))
    corr = read_tsv(os.path.join(d, "matrix", "chrom_label_corrections.tsv"))
    ploidy = read_tsv(os.path.join(d, "matrix", "ploidy_summary.tsv"))
    pairs = read_tsv(os.path.join(d, "homeologs", "homeolog_pairs.tsv"))
    hmap = metrics(os.path.join(d, "homeologs", "homeology_map_summary.tsv"))
    rsum = metrics(os.path.join(d, "rediploidization", "rediploidization_summary.tsv"))
    by_chrom = read_tsv(os.path.join(d, "rediploidization", "rediploidization_by_chrom.tsv"))
    fusions = [f for f in read_tsv(os.path.join(d, "rediploidization", "fusions.tsv"))
               if f.get("status") == "fusion"]
    answers = {r["question"]: r for r in read_tsv(os.path.join(d, "summary.tsv"))}
    sync = (read_tsv(os.path.join(d, "structure", "pair_synchrony.tsv")) or [{}])[0]
    allele, cross = copy_divergence(d)

    haps = Counter(r["hap"] for r in seqs)
    copies = Counter(r["n_haplotype_copies"] for r in ploidy)
    states = Counter(r["copy_state"] for r in by_chrom)
    status = Counter(r["status"] for r in corr)
    name, known, source = lit.get(sp, ("", "", ""))

    residual = [dict(chrom=r["chrom"], bp=int(r["residual_bp"]), support=r.get("residual_support", ""),
                     segments=r.get("residual_segments", ""))
                for r in by_chrom if (num(r.get("residual_bp")) or 0) > 0]
    exchange = [dict(chrom=r["chrom"], bp=int(r["exchange_bp"]), support=r.get("exchange_support", ""),
                     segments=r.get("exchange_segments", ""))
                for r in by_chrom if (num(r.get("exchange_bp")) or 0) > 0]
    regional = [dict(chrom=r["chrom"], state=r["copy_state"], segments=r.get("split_segments", ""))
                for r in by_chrom if r.get("split_extent") == "regional"]

    caveats = []
    if len(haps) <= 1:
        caveats.append("one usable haplotype: copy-level readings are not assessable")
    if status.get("ambiguous"):
        caveats.append(f"{status['ambiguous']} chromosome-numbering mismatch(es) left ambiguous (not applied)")
    if status.get("extra_set"):
        caveats.append("a haplotype file held a second chromosome set; split out as "
                       + ", ".join(sorted({r.get("new_hap", "") for r in corr if r["status"] == "extra_set"})))
    term = num(rsum.get("residual_terminal_frac"))
    if residual and term is not None and term >= 0.9:
        caveats.append("residual tetrasomy only at chromosome ends: homeolog-specific subtelomeric repeats "
                       "are not excluded")
    if len(copies) > 1:
        caveats.append("copy number varies between chromosomes: "
                       + ", ".join(f"{k} copies x{v}" for k, v in sorted(copies.items())))
    if sp in NOTES:
        caveats.append(NOTES[sp])

    return OrderedDict(
        species=sp, name=name, known=known, source=source,
        haplotypes=", ".join(f"{h} ({n})" for h, n in sorted(haps.items())),
        n_numbers=len(ploidy), copies=copies.most_common(1)[0][0] if copies else "",
        units=len(seqs), unplaced=len(unplaced),
        corrections=", ".join(f"{k} {v}" for k, v in sorted(status.items())) or "none",
        allele=allele, cross=cross,
        n_pairs=len(pairs),
        pair_depth=statistics.median([float(r["mean_distance"]) for r in pairs]) if pairs else None,
        pair_cv=sync.get("pair_depth_cv", ""),
        map=hmap,
        states=OrderedDict((s, states.get(s, 0)) for s in STATES if states.get(s)),
        fusions=sorted({f"{f['components']} ({f['hap']})" for f in fusions}),
        regional=regional,
        residual=residual, residual_summary=dict(
            tested=rsum.get("residual_tested_chromosomes", ""), bp=rsum.get("residual_bp", ""),
            terminal=rsum.get("residual_terminal_frac", ""), controlled=rsum.get("residual_controlled", "")),
        exchange=exchange,
        answers=OrderedDict((q, dict(answer=answers.get(q, {}).get("answer", ""),
                                     confidence=answers.get(q, {}).get("confidence", ""),
                                     evidence=answers.get(q, {}).get("evidence", "")))
                            for q in QUESTIONS),
        caveats=caveats,
    )


def short(text, n=90):
    text = text.split(":")[0] if ":" in text[:40] else text
    return text if len(text) <= n else text[: n - 1] + "…"


def fmt(x, nd=4):
    return f"{x:.{nd}f}" if isinstance(x, float) else (x or "–")


def markdown(records):
    out = ["# Species breakdown", "",
           "Generated by `docs/species/make_species_breakdown.py` from `results/` -- do not edit by hand.",
           "Every reading is from one individual's haplotype assemblies; states describe the data, not "
           "inheritance (see OUTPUTS.md).", "",
           "## Panel", "",
           "| Species | Name | Copies × numbers | Duplicated (map) | Origin-like | Rediploidization |",
           "|---|---|---|---|---|---|"]
    for r in records:
        a = r["answers"]
        dup = r["map"].get("duplicated_frac", "")
        out.append(f"| [{r['species']}](#{r['species'].lower()}) | {r['name']} | {r['copies']} × {r['n_numbers']} | "
                   f"{float(dup):.0%} | " if dup else
                   f"| [{r['species']}](#{r['species'].lower()}) | {r['name']} | {r['copies']} × {r['n_numbers']} | – | ")
        out[-1] += (f"{short(a['origin_like_structure']['answer'])} ({a['origin_like_structure']['confidence']}) | "
                    f"{short(a['rediploidization']['answer'])} ({a['rediploidization']['confidence']}) |")
    for r in records:
        m = r["map"]
        out += ["", f"## {r['species']}", "", f"*{r['name']}*" if r["name"] else "", ""]
        if r["known"]:
            out.append(f"**Known:** {r['known']}" + (f" ({r['source']})" if r["source"] else ""))
            out.append("")
        out += ["**Assembly**", "",
                f"- Haplotypes (chromosome-scale units): {r['haplotypes'] or '–'}; {r['n_numbers']} chromosome "
                f"numbers, {r['copies']} copies each (modal); {r['unplaced']} sequences unplaced",
                f"- Numbering changes by the matrix stage: {r['corrections']}",
                "", "**Copies and duplicated sets**", "",
                f"- Same-number (allelic) distance {fmt(r['allele'])}; different-number distance {fmt(r['cross'])}",
                f"- Whole-chromosome homeolog pairs: {r['n_pairs']}"
                + (f", median distance {r['pair_depth']:.3f}, pair-depth CV {r['pair_cv'] or '–'}"
                   if r["pair_depth"] else "")]
        if m:
            out.append(f"- Homeology map: {float(m.get('duplicated_frac', 0)):.0%} of the genome in "
                       f"{m.get('n_blocks', 0)} homeolog blocks; {m.get('n_partner_pairs', 0)} chromosome pairs; "
                       f"block divergence median {m.get('block_dist_median') or '–'}, CV {m.get('block_dist_cv') or '–'}")
            if m.get("multi_partner_list"):
                out.append(f"- Chromosomes with blocks on two or more others: "
                           f"{m['multi_partner_list'].replace(';', ', ')}")
        else:
            out.append("- Homeology map: not run")
        out += ["", "**Rediploidization signals**", ""]
        if r["states"]:
            out.append("- Per-chromosome states: " + ", ".join(f"{k.replace('_', ' ')} {v}"
                                                             for k, v in r["states"].items()))
        out.append(f"- Fusions between copies: {', '.join(r['fusions']) or 'none'}")
        if r["regional"]:
            out.append("- Regional lineage splits: " + "; ".join(
                f"{x['chrom']} ({x['state'].replace('_', ' ')}: {x['segments']})" for x in r["regional"]))
        rs = r["residual_summary"]
        if rs["tested"]:
            out.append(f"- Residual tetrasomy: {len(r['residual'])} of {rs['tested']} tested chromosomes, "
                       f"{mb(rs['bp']) or '0 Mb'}"
                       + (f", {float(rs['terminal']):.0%} at chromosome ends" if rs["terminal"] else "")
                       + ("" if rs["controlled"] == "yes" else " (no unrelated-chromosome control)"))
            for x in r["residual"]:
                out.append(f"  - {x['chrom']}: {mb(x['bp'])}, support {x['support']}; {x['segments']}")
            out.append(f"- Homeologous-exchange candidates: {len(r['exchange']) or 'none'}")
            for x in r["exchange"]:
                out.append(f"  - {x['chrom']}: {mb(x['bp'])}, copies {x['support']}; {x['segments']}")
        else:
            out.append("- Residual tetrasomy / exchange: not tested (no homeolog pairs or blocks)")
        out += ["", "**Summary answers**", "", "| Question | Answer | Confidence |", "|---|---|---|"]
        for q, title in QUESTIONS.items():
            a = r["answers"][q]
            out.append(f"| {title} | {a['answer'].replace('|', '/')} | {a['confidence']} |")
        out += ["", "<details><summary>Evidence</summary>", ""]
        for q, title in QUESTIONS.items():
            out.append(f"- {title}: {r['answers'][q]['evidence'] or '–'}")
        out += ["", "</details>"]
        if r["caveats"]:
            out += ["", "**Caveats**", ""] + [f"- {c}" for c in r["caveats"]]
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--json")
    args = ap.parse_args()
    lit = literature()
    records = []
    for sp in sorted(os.listdir(args.results), key=str.lower):
        d = os.path.join(args.results, sp)
        if os.path.exists(os.path.join(d, "summary.tsv")):
            records.append(species_record(sp, d, lit))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        f.write(markdown(records))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(records, f, indent=1, default=str)
    print(f"wrote {args.out} ({len(records)} species)")


if __name__ == "__main__":
    main()
