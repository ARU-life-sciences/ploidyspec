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
# Unlocalised pieces of a chromosome carry its number in NCBI-style headers
# ("chromosome 22 SUPER_22_unloc_1"); they are never the chromosome itself, and
# one over --min-len would collide with it (fLepOcu1).
UNLOCALISED = re.compile(r"unloc|unlocali[sz]ed|_random\b", re.I)


def is_unlocalised(seq_id, desc):
    return bool(UNLOCALISED.search(seq_id) or UNLOCALISED.search(desc))

# Automatic chromosome numbering (--chrom-naming auto, or detect with no named
# sequences). The reference haplotype is numbered 1..n by length; every other
# haplotype gets provisional numbers from AUTO_OFFSET, which the matrix stage
# replaces by one-to-one k-mer matching to the reference (chrom_reconcile).
AUTO_OFFSET = 1001
# A sequence is chromosome-scale, in auto mode, if it is at least this fraction
# of the median length of the haplotype's larger sequences (the top half of
# those >= --min-len). Masu salmon HAP2: 33 chromosomes of 18-127 Mb kept,
# unplaced contigs of 1-5 Mb dropped.
AUTO_MIN_FRAC = 0.1
CHROM_NAMING = ("detect", "names", "auto")
# A haplotype with this many times the median number of chromosome-scale
# sequences of the others is probably contig-level, not chromosome-scale.
CONTIG_LEVEL_FACTOR = 2.0


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
    """[(fasta_path, hap_label)]. Relative FASTA paths are resolved against the
    manifest's own directory, so a manifest and its FASTAs can be moved together."""
    base = os.path.dirname(os.path.abspath(manifest_path))
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
            fasta = os.path.expanduser(parts[0])
            rows.append((fasta if os.path.isabs(fasta) else os.path.join(base, fasta), parts[1]))
    return rows


def chromosome_scale(lengths, min_frac=AUTO_MIN_FRAC):
    """From {seq_id: length} (already >= min_len), the ids long enough to be
    chromosomes: >= min_frac x the median length of the longer half."""
    if not lengths:
        return set()
    ordered = sorted(lengths.values(), reverse=True)
    top = ordered[: max(1, (len(ordered) + 1) // 2)]
    cutoff = min_frac * top[len(top) // 2]
    return {s for s, n in lengths.items() if n >= cutoff}


def previous_auto_numbers(outdir):
    """{(source, seq_id): chrom} from an earlier auto-numbered sequences.tsv, so a
    resumed run keeps the numbers the matrix stage already assigned (and the
    k-mer tables cached under them)."""
    path = os.path.join(outdir, "sequences.tsv")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {(r["source"], r["seq_id"]): int(r["chrom"]) for r in csv.DictReader(f, delimiter="\t")
                if r.get("numbering") == "auto"}


def auto_number(candidates, outdir, min_frac=AUTO_MIN_FRAC):
    """
    Number chromosome-scale sequences without names. candidates: [dict(hap,
    seq_id, length, source, desc)] for sequences >= min_len. Returns (rows,
    excluded) where excluded are [(seq_id, source, length, reason)].
    """
    groups = {}
    for c in candidates:
        groups.setdefault((c["source"], c["hap"]), []).append(c)
    rows, excluded = [], []
    kept = {}
    for key, seqs in groups.items():
        ok = chromosome_scale({s["seq_id"]: s["length"] for s in seqs}, min_frac)
        kept[key] = sorted((s for s in seqs if s["seq_id"] in ok), key=lambda s: -s["length"])
        excluded += [(s["seq_id"], s["source"], s["length"], "auto-not-chromosome-scale")
                     for s in seqs if s["seq_id"] not in ok]
    if not kept:
        return rows, excluded
    # reference: most chromosome-scale sequences, ties to manifest order, among
    # haplotypes that do not look contig-level (far more sequences than the
    # others): masu salmon HAP1, 430 contigs against HAP2's 33 chromosomes
    order = list(groups)

    def contig_level(k):
        others = sorted(len(kept[o]) for o in order if o != k)
        return bool(others) and len(kept[k]) > CONTIG_LEVEL_FACTOR * others[len(others) // 2]

    eligible = [k for k in order if not contig_level(k)] or order
    ref = max(eligible, key=lambda k: (len(kept[k]), -order.index(k)))
    previous = previous_auto_numbers(outdir)
    reuse = all((s["source"], s["seq_id"]) in previous for seqs in kept.values() for s in seqs)
    for key in order:
        for rank, s in enumerate(kept[key]):
            if reuse:
                chrom = previous[(s["source"], s["seq_id"])]
            else:
                chrom = rank + 1 if key == ref else AUTO_OFFSET + rank
            rows.append(dict(s, chrom=chrom, unit_id=f"{s['hap']}_chr{chrom:02d}", numbering="auto"))
    log(f"auto chromosome numbering: reference {ref[1]} ({len(kept[ref])} chromosome-scale "
        f"sequences, numbered by length); other haplotypes numbered by k-mer matching in the "
        f"matrix stage" + (" (reusing the numbering of the existing sequences.tsv)" if reuse else ""))
    return rows, excluded



ROADMAP_17 = "not yet handled by ploidyspec (ROADMAP 1.7)"


def contig_level_notes(rows):
    """Warnings for haplotypes that look contig-level: far more chromosome-scale
    units than the other haplotypes, or none while the others have some."""
    counts = {}
    for r in rows:
        counts[(r["source"], r["hap"])] = counts.get((r["source"], r["hap"]), 0) + 1
    notes = []
    if len(counts) < 2:
        return notes
    for key, n in counts.items():
        others = sorted(v for k, v in counts.items() if k != key)
        typical = others[len(others) // 2]
        if typical and n > CONTIG_LEVEL_FACTOR * typical:
            notes.append(
                f"{key[1]} ({os.path.basename(key[0])}) has {n} chromosome-scale sequences against "
                f"{typical} in the other haplotypes: it looks contig-level. Contig-level haplotypes "
                f"are {ROADMAP_17}; place its contigs on a scaffolded haplotype of the same "
                f"individual first with workflows/prep/scaffold_by_reference.py.")
    return notes


def prepare(manifest_path, outdir, samtools_bin, min_len, chrom_regexes, hap_regex,
            chrom_naming="detect"):
    """
    chrom_naming: "names" reads chromosome numbers from headers
    (--chrom-regex); "auto" ignores names and numbers by length and k-mer
    matching (auto_number); "detect" uses names when any sequence has one and
    falls back to auto otherwise.
    """
    chrom_regexes = chrom_regexes or DEFAULT_CHROM_REGEXES
    hap_regex = hap_regex or DEFAULT_HAP_REGEX

    rows = []
    unplaced = []
    candidates = []  # >= min_len with a haplotype, for auto numbering
    notes = []
    for fasta, hap_label in read_manifest(manifest_path):
        if not os.path.exists(fasta):
            raise SystemExit(f"manifest references missing file: {fasta}")
        lengths = read_fai(samtools_bin, fasta)
        descs = read_descriptions(fasta)
        n_before = len(candidates)
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
            if is_unlocalised(seq_id, desc):
                unplaced.append((seq_id, fasta, length, "unlocalised"))
                continue
            candidates.append(dict(hap=hap, seq_id=seq_id, length=length, source=fasta, desc=desc))
            chrom = None if chrom_naming == "auto" else (
                match_chrom(desc, chrom_regexes) or match_chrom(seq_id, chrom_regexes))
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
                    numbering="names",
                )
            )
        if hap_label.upper() == "AUTO" and len(candidates) == n_before and any(
                n >= min_len for n in lengths.values()):
            notes.append(
                f"{os.path.basename(fasta)} is labelled AUTO but no header matched the haplotype "
                f"pattern (--hap-regex), so none of its sequences were used. Separating several "
                f"haplotypes inside one file without names is {ROADMAP_17}: split the file by "
                f"haplotype and label each in the manifest, or pass --hap-regex.")

    if chrom_naming == "auto" or (chrom_naming == "detect" and not rows and candidates):
        if chrom_naming == "detect":
            log("no sequence header matched a chromosome-number pattern -- numbering "
                "chromosomes automatically (--chrom-naming auto)")
        rows, excluded = auto_number(candidates, outdir)
        unplaced = [u for u in unplaced if u[3] != "no-chrom-match"] + excluded

    notes += contig_level_notes(rows)
    named = {(r["source"], r["hap"]) for r in rows}
    if rows and any(r["numbering"] == "names" for r in rows):
        for key in {(c["source"], c["hap"]) for c in candidates} - named:
            notes.append(
                f"{key[1]} ({os.path.basename(key[0])}) has no sequence with a chromosome number while "
                f"the other haplotypes do. If its headers simply use another naming scheme, pass "
                f"--chrom-regex or --chrom-naming auto. If it is contig-level, that is {ROADMAP_17}: "
                f"place its contigs on a scaffolded haplotype first with "
                f"workflows/prep/scaffold_by_reference.py.")
    for note in notes:
        log(f"NOTE: {note}")

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
    fieldnames = ["unit_id", "hap", "chrom", "seq_id", "length", "source", "desc", "numbering"]
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
