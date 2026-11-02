#!/usr/bin/env python3
"""Metadata flags for a ToL release, before any k-mer work: chromosome counts
and span from the insdc assembly report, BUSCO duplication from the release's
analysis directory, and karyotype / genome size / ploidy records from GoaT
(goat.genomehubs.org). Each flag compares the assembly with what is recorded
for the species or its genus:

- karyotype_ratio: GoaT's directly recorded haploid number over the assembled
  chromosomes per haplotype. daBudDavi1 assembles 19 per haplotype against a
  recorded n = 38 (a 4x whose two haplotypes each hold one of two subgenome
  pairs' worth): 2.0.
- genus_ratio: assembled chromosomes per haplotype over the genus's modal
  haploid number (a 4x on the genus's base number reads ~2), when at least
  GENUS_MIN_RECORDS species of the genus have one. Weak in groups with free
  fission and fusion (holocentric Lepidoptera) and in birds' microchromosomes.
- busco_dup near 100% means every BUSCO is in two copies within one haplotype:
  both subgenomes assembled into each haplotype (Cardamine pratensis 93%).
- size_ratio: GoaT's directly measured genome size over the haplotype span.
  C-values of polyploids are often holoploid (all sets), so ~2 is a flag.
- busco_dup: % duplicated BUSCOs in hap1, most specific lineage.

Usage: metadata.py CATALOG.tsv TOLIDS_FILE OUT.tsv [CACHE.json]
"""
import csv
import glob
import json
import os
import sys
import time
import urllib.parse
import urllib.request

GOAT = "https://goat.genomehubs.org/api/v2/search"
GENUS_MIN_RECORDS = 3  # a genus mode from fewer species' records is not used
FIELDS = ["tolid", "species", "n_haps", "chroms_hap1", "chroms_hap2", "span_hap1_mb",
          "busco_lineage", "busco_dup", "goat_haploid", "goat_ploidy", "goat_genome_size_mb",
          "genus_haploid_mode", "karyotype_ratio", "genus_ratio", "size_ratio"]


def report_counts(fasta):
    """(chromosomes, span bp) from the insdc assembly report next to the FASTA"""
    reports = glob.glob(os.path.join(os.path.dirname(fasta), "*.assembly_report.txt"))
    if not reports:
        return None, None
    n = span = 0
    with open(reports[0]) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > 8 and p[1] == "assembled-molecule" and p[3] == "Chromosome":
                n += 1
                span += int(p[8]) if p[8].isdigit() else 0
    return n, span


def busco_dup(analysis):
    best = None
    for path in glob.glob(os.path.join(analysis, "busco", "*", "*.short_summary.json")):
        lineage = os.path.basename(os.path.dirname(path))
        if lineage.startswith(("archaea", "bacteria")):
            continue
        try:
            with open(path) as f:
                res = json.load(f)["results"]
        except (OSError, ValueError, KeyError):
            continue
        n = res.get("n_markers", 0)
        if best is None or n > best[0]:
            best = (n, lineage, res.get("Multi copy percentage"))
    return (best[1], best[2]) if best else ("", None)


def goat(query, fields, estimates=False):
    params = dict(query=query, result="taxon", fields=",".join(fields), taxonomy="ncbi",
                  includeEstimates="true" if estimates else "false")
    url = GOAT + "?" + urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                res = json.load(r).get("results", [])
            return res[0]["result"]["fields"] if res else {}
        except Exception:
            time.sleep(5 * (attempt + 1))
    return None


def genus_taxid(taxon_id):
    """The NCBI taxid of the taxon's genus, from GoaT's lineage: genus names are
    not unique across kingdoms (Linaria, Oenanthe: plants and birds)."""
    url = GOAT.replace("/search", "/record") + "?" + urllib.parse.urlencode(
        dict(recordId=taxon_id, result="taxon", taxonomy="ncbi"))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                rec = json.load(r)["records"][0]["record"]
            return next((x["taxon_id"] for x in rec.get("lineage", []) if x.get("taxon_rank") == "genus"), None)
        except Exception:
            time.sleep(5 * (attempt + 1))
    return None


def direct(fields, name):
    v = (fields or {}).get(name) or {}
    src = v.get("aggregation_source")
    src = src if isinstance(src, list) else [src]
    return v.get("value") if "direct" in src else None


def ratio(a, b):
    return f"{a / b:.2f}" if a and b else ""


def main():
    catalog, tolids, out = sys.argv[1:4]
    cache_path = sys.argv[4] if len(sys.argv) > 4 else None
    cache = {}
    if cache_path and os.path.exists(cache_path):
        with open(cache_path) as f:
            cache = json.load(f)
    want = {line.strip() for line in open(tolids) if line.strip()}
    rows = [r for r in csv.DictReader(open(catalog), delimiter="\t") if r["tolid"] in want]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for r in rows:
            fastas = r["fastas"].split(";")
            (c1, span1), (c2, _) = report_counts(fastas[0]), report_counts(fastas[1])
            lineage, dup = busco_dup(r["hap1_analysis"])
            key = f"taxon:{r['taxon_id']}"
            if key not in cache and r["taxon_id"]:
                sp = goat(f"tax_eq({r['taxon_id']})", ["haploid_number", "ploidy", "genome_size"])
                if sp is not None:
                    cache[key] = dict(haploid=direct(sp, "haploid_number"), ploidy=direct(sp, "ploidy"),
                                      genome_size=direct(sp, "genome_size"))
            gkey = f"genus_of:{r['taxon_id']}"
            if gkey not in cache and r["taxon_id"]:
                gid = genus_taxid(r["taxon_id"])
                g = goat(f"tax_eq({gid})", ["haploid_number"]) if gid else None
                if g is not None:
                    h = g.get("haploid_number") or {}
                    count = h.get("count") or 0
                    cache[gkey] = dict(mode=(h.get("mode") or h.get("value")) if count >= GENUS_MIN_RECORDS else None)
            sp, g = cache.get(key, {}), cache.get(gkey, {})
            gs = sp.get("genome_size")
            w.writerow(dict(
                tolid=r["tolid"], species=r["species"], n_haps=r["n_haps"],
                chroms_hap1="" if c1 is None else c1, chroms_hap2="" if c2 is None else c2,
                span_hap1_mb=f"{span1 / 1e6:.0f}" if span1 else "",
                busco_lineage=lineage, busco_dup="" if dup is None else dup,
                goat_haploid=sp.get("haploid") or "", goat_ploidy=sp.get("ploidy") or "",
                goat_genome_size_mb=f"{gs / 1e6:.0f}" if gs else "",
                genus_haploid_mode=g.get("mode") or "",
                karyotype_ratio=ratio(sp.get("haploid"), c1),
                genus_ratio=ratio(c1, g.get("mode")),
                size_ratio=ratio(gs, span1)))
    if cache_path:
        with open(cache_path + ".part", "w") as f:
            json.dump(cache, f)
        os.replace(cache_path + ".part", cache_path)


if __name__ == "__main__":
    main()
