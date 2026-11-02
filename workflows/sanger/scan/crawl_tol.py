#!/usr/bin/env python3
"""List every Tree of Life release with two or more haplotype assemblies.

Walks /lustre/scratch122/tol/data/<6 hex levels>/<Species>/assembly/release/
and keeps, per ToLID, the latest version that has a hap1 and at least a hap2
(`<tolid>.hap<N>.<version>` symlinks to the GCA directory). One row per
release: tolid, species, taxon id, version, number of haplotypes, the insdc
FASTA of each (bgzipped and indexed by the ToL pipeline, so samtools reads
them in place), their release time (FASTA mtime) and the analysis directory
of hap1 (BUSCO). Usage: crawl_tol.py OUT.tsv [ROOT]
"""
import json
import os
import re
import sys
import time

ROOT = "/lustre/scratch122/tol/data"
HAP = re.compile(r"^(?P<tolid>[a-zA-Z]+\d+)\.hap(?P<hap>\d+)\.(?P<ver>\d+)$")
FIELDS = ["tolid", "species", "taxon_id", "version", "n_haps", "released", "fastas", "hap1_analysis"]


def species_dirs(root):
    level = [root]
    for _ in range(6):
        nxt = []
        for d in level:
            try:
                nxt += [e.path for e in os.scandir(d) if e.is_dir(follow_symlinks=False) and len(e.name) == 1]
            except OSError:
                pass
        level = nxt
    for d in level:
        try:
            yield from (e.path for e in os.scandir(d) if e.is_dir(follow_symlinks=False))
        except OSError:
            pass


def insdc_fasta(release, link):
    d = os.path.join(release, link, "insdc")
    try:
        names = os.listdir(d)
    except OSError:
        return None
    fa = [n for n in names if n.endswith((".fa.gz", ".fasta.gz")) and n.startswith("GCA_")]
    return os.path.join(d, sorted(fa)[0]) if fa else None


def releases(species_dir):
    release = os.path.join(species_dir, "assembly", "release")
    try:
        names = os.listdir(release)
    except OSError:
        return
    by = {}
    for n in names:
        m = HAP.match(n)
        if m:
            by.setdefault((m["tolid"], int(m["ver"])), {})[int(m["hap"])] = n
    latest = {}
    for (tolid, ver), haps in by.items():
        if 1 in haps and 2 in haps and ver > latest.get(tolid, (0, None))[0]:
            latest[tolid] = (ver, haps)
    if not latest:
        return
    try:
        with open(os.path.join(species_dir, "speciesops.json")) as f:
            ops = json.load(f)
    except (OSError, ValueError):
        ops = {}
    for tolid, (ver, haps) in sorted(latest.items()):
        fastas = [insdc_fasta(release, haps[h]) for h in sorted(haps)]
        if not all(fastas):
            continue
        yield dict(tolid=tolid, species=ops.get("species_name", os.path.basename(species_dir).replace("_", " ")),
                   taxon_id=ops.get("taxon_id", ""), version=ver, n_haps=len(haps),
                   released=time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(fastas[0]))),
                   fastas=";".join(fastas),
                   hap1_analysis=os.path.join(species_dir, "analysis", haps[1]))


def main():
    out = sys.argv[1]
    root = sys.argv[2] if len(sys.argv) > 2 else ROOT
    n = 0
    with open(out + ".part", "w") as f:
        f.write("\t".join(FIELDS) + "\n")
        for sd in species_dirs(root):
            for r in releases(sd):
                f.write("\t".join(str(r[k]) for k in FIELDS) + "\n")
                n += 1
    os.replace(out + ".part", out)
    print(f"{n} releases with >= 2 haplotypes", file=sys.stderr)


if __name__ == "__main__":
    main()
