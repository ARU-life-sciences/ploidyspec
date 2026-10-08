"""
Place the contigs of an unscaffolded haplotype onto the chromosomes of a
scaffolded one from the same individual, so both can go through ploidyspec
with the same chromosome numbers.

Input: the contig FASTA (indexed), and a PAF of those contigs aligned to the
scaffolded haplotype (minimap2 -x asm5 -c). Each contig of at least
--min-contig bp is assigned to the reference chromosome that takes most of its
aligned bases, if that chromosome holds >= --min-share of them and the
alignments cover >= --min-cover of the contig. Contigs on a chromosome are
ordered by the weighted median reference position of their alignments and
oriented by the majority strand, then joined with runs of N into one
pseudo-chromosome per reference chromosome.

Output: a FASTA with one pseudo-chromosome per reference chromosome (same name,
`chromosome: N` in the description when the name ends in a number), the
unplaced contigs unchanged, and a TSV of every contig's placement.

The order and orientation come from the other haplotype, so positions along
these pseudo-chromosomes are reference-guided. ploidyspec's whole-chromosome
distances do not depend on order, and its windows are position-free, so the
effect is limited to where along the chromosome a windowed segment is reported.
"""

import argparse
import re
import subprocess
from collections import defaultdict

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fai(path):
    with open(path) as f:
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in f}


def read_paf(path, min_mapq):
    """{query: {target: [(aligned_bp, target_mid, strand)]}}"""
    hits = defaultdict(lambda: defaultdict(list))
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if int(p[11]) < min_mapq:
                continue
            q, strand, t = p[0], p[4], p[5]
            ts, te, nmatch = int(p[7]), int(p[8]), int(p[9])
            hits[q][t].append((nmatch, (ts + te) / 2, strand))
    return hits


def weighted_median(pairs):
    pairs = sorted(pairs, key=lambda p: p[1])
    half = sum(w for w, _ in pairs) / 2
    acc = 0
    for w, v in pairs:
        acc += w
        if acc >= half:
            return v
    return pairs[-1][1]


def place(qlen, hits, targets, min_contig, min_share, min_cover):
    placements, rows = defaultdict(list), []
    for q, length in qlen.items():
        by_t = {t: hs for t, hs in hits.get(q, {}).items() if t in targets}
        aligned = {t: sum(h[0] for h in hs) for t, hs in by_t.items()}
        total = sum(aligned.values())
        status, target, share, cover, pos, strand = "unplaced", "", 0.0, total / length, None, "+"
        if length >= min_contig and total:
            target = max(aligned, key=aligned.get)
            share = aligned[target] / total
            if share >= min_share and cover >= min_cover:
                hs = by_t[target]
                pos = weighted_median([(h[0], h[1]) for h in hs])
                plus = sum(h[0] for h in hs if h[2] == "+")
                strand = "+" if plus >= aligned[target] / 2 else "-"
                status = "placed"
                placements[target].append((pos, q, strand))
        rows.append(dict(contig=q, length=length, status=status, target=target,
                         share=f"{share:.3f}", cover=f"{min(cover, 1):.3f}",
                         position=f"{pos:.0f}" if pos is not None else "", strand=strand))
    return placements, rows


def fetch(fasta, name):
    out = subprocess.run(["samtools", "faidx", fasta, name], capture_output=True, text=True,
                         check=True).stdout
    return "".join(out.split("\n")[1:])


def write_fasta(f, name, desc, seq, width=80):
    f.write(f">{name} {desc}\n" if desc else f">{name}\n")
    for i in range(0, len(seq), width):
        f.write(seq[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contigs", required=True, help="contig FASTA (samtools faidx indexed)")
    ap.add_argument("--reference-fai", required=True, help=".fai of the scaffolded haplotype")
    ap.add_argument("--paf", required=True, help="minimap2 PAF of contigs vs scaffolded haplotype")
    ap.add_argument("--targets", default=r"^group\d+$",
                    help="regex for reference sequences that are chromosomes (default ^group\\d+$)")
    ap.add_argument("--min-ref-len", type=int, default=1_000_000)
    ap.add_argument("--min-contig", type=int, default=100_000)
    ap.add_argument("--min-share", type=float, default=0.5)
    ap.add_argument("--min-cover", type=float, default=0.2)
    ap.add_argument("--min-mapq", type=int, default=5)
    ap.add_argument("--gap", type=int, default=100)
    ap.add_argument("--out", required=True, help="output FASTA")
    ap.add_argument("--placements", required=True, help="output TSV of contig placements")
    a = ap.parse_args()

    ref = read_fai(a.reference_fai)
    targets = {t for t, n in ref.items() if re.search(a.targets, t) and n >= a.min_ref_len}
    qlen = read_fai(a.contigs + ".fai")
    placements, rows = place(qlen, read_paf(a.paf, a.min_mapq), targets,
                             a.min_contig, a.min_share, a.min_cover)

    with open(a.placements, "w") as f:
        cols = list(rows[0])
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    placed = set()
    with open(a.out, "w") as f:
        for t in sorted(targets, key=lambda x: [int(n) if n.isdigit() else n for n in re.split(r"(\d+)", x)]):
            parts = []
            for _, q, strand in sorted(placements.get(t, [])):
                seq = fetch(a.contigs, q)
                parts.append(seq if strand == "+" else seq.translate(COMPLEMENT)[::-1])
                placed.add(q)
            if not parts:
                continue
            num = re.search(r"(\d+)$", t)
            write_fasta(f, t, f"chromosome: {num.group(1)}" if num else "", ("N" * a.gap).join(parts))
        for q in qlen:
            if q not in placed:
                write_fasta(f, q, "", fetch(a.contigs, q))
    n_placed = len(placed)
    bp = sum(qlen[q] for q in placed)
    print(f"placed {n_placed} contigs ({bp / 1e9:.2f} Gb, {bp / sum(qlen.values()):.0%}) "
          f"on {len(placements)} of {len(targets)} reference chromosomes")


if __name__ == "__main__":
    main()
