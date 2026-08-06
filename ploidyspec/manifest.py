import csv
import os
import re
import subprocess

from .common import log, run

DEFAULT_CHROM_REGEXES = [
    r"chromosome:?\s*(\d+)\b",
    r"SUPER[_-](\d+)\b(?!_)",
]
DEFAULT_HAP_REGEX = r"HAP(\d+)[_-]SUPER"


def read_fai(samtools_bin, fasta):
    fai = fasta + ".fai"
    if not os.path.exists(fai):
        log(f"indexing {fasta}")
        run([samtools_bin, "faidx", fasta])
    lengths = {}
    with open(fai) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def read_descriptions(fasta):
    """One streaming decompress pass to pull full FASTA deflines (id -> description)."""
    descs = {}
    if fasta.endswith(".gz"):
        proc = subprocess.Popen(
            ["gzip", "-dc", fasta], stdout=subprocess.PIPE, text=True
        )
        stream = proc.stdout
    else:
        proc = None
        stream = open(fasta)
    try:
        for line in stream:
            if line.startswith(">"):
                seq_id, _, desc = line[1:].rstrip("\n").partition(" ")
                descs[seq_id] = desc
    finally:
        if proc:
            proc.stdout.close()
            proc.wait()
        else:
            stream.close()
    return descs


def match_chrom(desc, regexes):
    for pat in regexes:
        m = re.search(pat, desc, re.IGNORECASE)
        if m:
            return int(m.group(1))
    return None


def match_hap(desc, hap_regex):
    m = re.search(hap_regex, desc, re.IGNORECASE)
    return f"HAP{m.group(1)}" if m else None


def read_manifest(manifest_path):
    rows = []
    with open(manifest_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise SystemExit(
                    f"manifest line malformed (need <fasta>\\t<hap_label|AUTO>): {line!r}"
                )
            rows.append((parts[0], parts[1]))
    return rows


def prepare(manifest_path, outdir, samtools_bin, min_len, chrom_regexes, hap_regex):
    chrom_regexes = chrom_regexes or DEFAULT_CHROM_REGEXES
    hap_regex = hap_regex or DEFAULT_HAP_REGEX

    rows = []
    unplaced = []
    for fasta, hap_label in read_manifest(manifest_path):
        fasta = os.path.abspath(fasta)
        if not os.path.exists(fasta):
            raise SystemExit(f"manifest references missing file: {fasta}")
        lengths = read_fai(samtools_bin, fasta)
        descs = read_descriptions(fasta)
        for seq_id, length in lengths.items():
            desc = descs.get(seq_id, "")
            if length < min_len:
                unplaced.append((seq_id, fasta, length, "below-min-len"))
                continue
            if hap_label.upper() == "AUTO":
                # "curated" assemblies use bare headers (e.g. ">SUPER_1", no
                # space) so desc comes back empty -- fall back to seq_id,
                # which is where release-vs-curated conventions put the tag.
                hap = match_hap(desc, hap_regex) or match_hap(seq_id, hap_regex)
                if hap is None:
                    unplaced.append((seq_id, fasta, length, "no-hap-match"))
                    continue
            else:
                hap = hap_label
            chrom = match_chrom(desc, chrom_regexes) or match_chrom(
                seq_id, chrom_regexes
            )
            if chrom is None:
                unplaced.append((seq_id, fasta, length, "no-chrom-match"))
                continue
            unit_id = f"{hap}_chr{chrom:02d}"
            rows.append(
                dict(
                    unit_id=unit_id,
                    hap=hap,
                    chrom=chrom,
                    seq_id=seq_id,
                    length=length,
                    source=fasta,
                    desc=desc,
                )
            )

    seen = {}
    for r in rows:
        if r["unit_id"] in seen and seen[r["unit_id"]] != r["seq_id"]:
            raise SystemExit(
                f"duplicate unit_id {r['unit_id']!r}: both {seen[r['unit_id']]!r} and {r['seq_id']!r} "
                f"map to it -- refine --chrom-regex/--hap-regex"
            )
        seen[r["unit_id"]] = r["seq_id"]

    rows.sort(key=lambda r: (r["hap"], r["chrom"]))
    os.makedirs(outdir, exist_ok=True)
    seq_tsv = os.path.join(outdir, "sequences.tsv")
    fieldnames = ["unit_id", "hap", "chrom", "seq_id", "length", "source", "desc"]
    with open(seq_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if unplaced:
        with open(os.path.join(outdir, "unplaced.tsv"), "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["seq_id", "source", "length", "reason"])
            w.writerows(unplaced)

    haps = sorted(set(r["hap"] for r in rows))
    chroms = sorted(set(r["chrom"] for r in rows))
    log(
        f"prepared {len(rows)} chromosome-scale units "
        f"({len(haps)} haplotypes: {', '.join(haps)}; {len(chroms)} chromosome numbers); "
        f"{len(unplaced)} sequences excluded -> {os.path.join(outdir, 'unplaced.tsv')}"
    )
    return seq_tsv
