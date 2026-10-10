"""
Write docs/guide/handbook_species.typ (the handbook's species chapter) from the
result directories, using the same per-species records as
docs/species/make_species_breakdown.py. Run from the repository root:

    python3 docs/guide/make_handbook_species.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "species"))
sys.path.insert(0, os.path.join(HERE, "..", ".."))

from make_species_breakdown import QUESTIONS, literature, species_record  # noqa: E402

OUT = os.path.join(HERE, "handbook_species.typ")

GROUPS = [
    ("Two-haplotype genomes with duplicated sets (allo-like anchors and relatives)",
     ["daGleHede1", "drTriRepe1", "drSorDevo1", "ddSalPent1", "drMyrVert1", "drMyrSpic1", "drRosSpin1",
      "daLatSqua1", "daLatClan1", "daSenVulg1", "laPotCris1", "laPotLuce1", "laPotNata1", "laPotNodo1",
      "laPotPerf1", "daPilAura1", "llColAutu1", "lpElePalu1"]),
    ("Cryptic or partial duplication in two-haplotype genomes",
     ["daInuConz1", "dmRanRepe1", "dcCerAlpi1", "daSonOler1", "dcHonPepl1"]),
    ("Genomes assembled with three or more copies",
     ["SchCurv1", "SchYoun1", "drLytSali1", "daBudDavi1", "daGalBore1", "ddLepDrab1", "ddEmpNigr1",
      "drAriEdul1", "ddHesMatr1", "ddHypMacu1", "ddHypPerf1", "ddSalTria1", "icStrMela3"]),
    ("Salmonids", ["OncMaso1", "fSalAlp3", "fCorLav1"]),
    ("Diploid-like and negative controls",
     ["ddMalSylv1", "fLepOcu1", "fAmiCal2", "xgCepNemo3", "xgHygCinc1", "xgDauRufa1", "xgLitLitt3",
      "xgMonCant1", "xgPomEleg1", "xgStaPalu1", "lpTriTurg1_A", "lpTriTurg1_B"]),
    ("One usable haplotype", ["daEupConf1", "ddAraThal4", "ddPopNigr1", "ddSalCine1", "drIngLaur1",
                              "drTriDubi3", "drUrtDioi1", "lpTriTurg1_AB"]),
]


def s(text):
    """A Typst string literal: shown verbatim, no markup."""
    return '"' + str(text).replace("\\", "\\\\").replace('"', '\\"') + '"'


def mb(bp):
    try:
        return f"{float(bp) / 1e6:.1f} Mb"
    except (TypeError, ValueError):
        return ""


def short(text, n=70):
    head = text.split(":")[0] if ":" in text[:45] else text.split(".")[0]
    return head if len(head) <= n else head[: n - 1] + "…"


def facts(r):
    m = r["map"]
    rows = [
        ("Assembly", f"{r['haplotypes'] or '–'}; {r['n_numbers']} chromosome numbers, {r['copies']} copies "
                     f"each; {r['unplaced']} sequences unplaced"),
        ("Numbering fixes", r["corrections"]),
        ("Distances", f"same number {r['allele']:.4f}, different numbers {r['cross']:.4f}"
         if isinstance(r["allele"], float) and isinstance(r["cross"], float) else "–"),
        ("Homeolog pairs", f"{r['n_pairs']}" + (f" (median {r['pair_depth']:.3f}, CV {r['pair_cv'] or '–'})"
                                              if r["pair_depth"] else "")),
    ]
    if m:
        rows.append(("Homeology map",
                     f"{float(m.get('duplicated_frac', 0)):.0%} of the genome in {m.get('n_blocks', 0)} blocks; "
                     f"{m.get('n_partner_pairs', 0)} chromosome pairs; block CV {m.get('block_dist_cv') or '–'}"
                     + (f"; two or more partners: {m['multi_partner_list'].replace(';', ', ')}"
                        if m.get("multi_partner_list") else "")))
    if r["states"]:
        rows.append(("States", ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in r["states"].items())))
    rows.append(("Fusions", ", ".join(r["fusions"]) or "none"))
    rs = r["residual_summary"]
    if rs["tested"]:
        res = (f"{len(r['residual'])} of {rs['tested']} tested, {mb(rs['bp']) or '0 Mb'}"
               + (f", {float(rs['terminal']):.0%} terminal" if rs["terminal"] else ""))
        if r["residual"]:
            res += "; " + "; ".join(f"{x['chrom']} {mb(x['bp'])} ({x['support']})" for x in r["residual"][:6])
            if len(r["residual"]) > 6:
                res += f"; +{len(r['residual']) - 6} more"
        rows.append(("Residual tetrasomy", res))
        ex = "; ".join(f"{x['chrom']} {mb(x['bp'])} ({x['support']})" for x in r["exchange"][:6]) or "none"
        if len(r["exchange"]) > 6:
            ex += f"; +{len(r['exchange']) - 6} more"
        rows.append(("Exchange candidates", ex))
    return rows


def card(r):
    out = [f"#spcard({s(r['species'])}, {s(r['name'])}, {s(r['known'])}, {s(r['source'])}, ("]
    for k, v in facts(r):
        out.append(f"  ({s(k)}, {s(v)}),")
    out.append("), (")
    for q, title in QUESTIONS.items():
        a = r["answers"][q]
        out.append(f"  ({s(title)}, {s(a['confidence'] or '–')}, {s(a['answer'] or '–')}, {s(a['evidence'] or '')}),")
    out.append("), (")
    for c in r["caveats"]:
        out.append(f"  {s(c)},")
    out.append("))")
    return "\n".join(out)


def main():
    lit = literature()
    records = {}
    for sp in os.listdir("results"):
        d = os.path.join("results", sp)
        if os.path.exists(os.path.join(d, "summary.tsv")):
            records[sp] = species_record(sp, d, lit)
    placed = {sp for _, members in GROUPS for sp in members}
    rest = sorted(set(records) - placed, key=str.lower)
    groups = GROUPS + ([("Other", rest)] if rest else [])

    out = ["// generated by docs/guide/make_handbook_species.py -- do not edit by hand",
           '#import "handbook_lib.typ": *', ""]
    out.append("#table(columns: (1.05fr, 1.45fr, 0.62fr, 0.62fr, 1.55fr, 1.55fr),")
    out.append("  table.header([*Species*], [*Name*], [*Copies × nos.*], [*Dup. (map)*], "
               "[*Origin-like structure*], [*Rediploidization*]),")
    for title, members in groups:
        members = [m for m in members if m in records]
        if not members:
            continue
        out.append(f"  table.cell(colspan: 6, fill: soft, text(weight: \"bold\", {s(title)})),")
        for sp in members:
            r = records[sp]
            dup = r["map"].get("duplicated_frac", "")
            a = r["answers"]
            out.append(
                f"  raw({s(sp)}), emph({s(r['name'])}), {s(str(r['copies']) + ' × ' + str(r['n_numbers']))}, "
                f"{s(f'{float(dup):.0%}' if dup else '–')}, "
                f"{s(short(a['origin_like_structure']['answer']) + ' (' + a['origin_like_structure']['confidence'] + ')')}, "
                f"{s(short(a['rediploidization']['answer']) + ' (' + a['rediploidization']['confidence'] + ')')},")
    out.append(")")
    out.append("")
    for title, members in groups:
        members = [m for m in members if m in records]
        if not members:
            continue
        out.append(f"== {title}")
        out.append("")
        for sp in members:
            out.append(card(records[sp]))
            out.append("")
    with open(OUT, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {OUT} ({len(records)} species)")


if __name__ == "__main__":
    main()
