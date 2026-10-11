"""
Alignment dotplot of two chromosome copies, for checking a ploidyspec
reading by hand (a residual-tetrasomy block, an exchange, a fusion, a
breakpoint). Aligns query against target with minimap2 (-x asm20: tolerates
homeolog-level divergence), writes the PAF and a dotplot coloured by
identity, and prints identity by target position bins.

Usage: dotplot_pair.py TARGET.fa QUERY.fa OUT_PREFIX [--minimap2 PATH] [--threads N]
Needs minimap2 on PATH (module load minimap2) and matplotlib.
"""

import argparse
import os
import subprocess
from collections import defaultdict


def run_minimap2(minimap2, target, query, paf, threads, approx=False):
    """approx: no base-level alignment (-c); identity is then matches / block
    length. Needed for very large repetitive chromosomes (ddHesMatr1, ~850 Mb),
    where base-level alignment exceeded 64 GB."""
    with open(paf, "w") as out:
        subprocess.run([minimap2, "-x", "asm20", *([] if approx else ["-c"]), "-t", str(threads),
                        "--secondary=no", target, query], stdout=out, check=True)


def read_paf(paf, min_len=5000):
    rows = []
    with open(paf) as f:
        for line in f:
            p = line.split("\t")
            qlen, qs, qe, strand = int(p[1]), int(p[2]), int(p[3]), p[4]
            tlen, ts, te = int(p[6]), int(p[7]), int(p[8])
            matches, block = int(p[9]), int(p[10])
            if block < min_len:
                continue
            de = next((float(x[5:]) for x in p[12:] if x.startswith("de:f:")), None)
            ident = 1 - de if de is not None else matches / block
            rows.append((qs, qe, ts, te, strand, ident, qlen, tlen))
    return rows


def identity_profile(rows, tlen, bin_size):
    """Best alignment identity covering each target bin."""
    best = defaultdict(float)
    for qs, qe, ts, te, strand, ident, _, _ in rows:
        for b in range(ts // bin_size, te // bin_size + 1):
            best[b] = max(best[b], ident)
    return [(b * bin_size, best.get(b)) for b in range(tlen // bin_size + 1)]


def plot(rows, out_png, tname, qname):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 6))
    cmap = plt.get_cmap("viridis")
    for qs, qe, ts, te, strand, ident, qlen, tlen in rows:
        y = (qs, qe) if strand == "+" else (qe, qs)
        c = cmap(min(1.0, max(0.0, (ident - 0.90) / 0.10)))
        ax.plot((ts / 1e6, te / 1e6), (y[0] / 1e6, y[1] / 1e6), color=c, linewidth=1.2)
    if rows:
        ax.set_xlim(0, rows[0][7] / 1e6)
        ax.set_ylim(0, rows[0][6] / 1e6)
    ax.set_xlabel(f"{tname} (Mb)")
    ax.set_ylabel(f"{qname} (Mb)")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.90, 1.0))
    fig.colorbar(sm, ax=ax, label="alignment identity")
    ax.set_title(f"{qname} vs {tname} (minimap2 asm20, blocks >= 5 kb)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target")
    ap.add_argument("query")
    ap.add_argument("out_prefix")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--bin", type=int, default=1_000_000)
    ap.add_argument("--approx", action="store_true", help="no base-level alignment (large genomes)")
    args = ap.parse_args()
    tname = os.path.basename(args.target).split(".")[0]
    qname = os.path.basename(args.query).split(".")[0]
    paf = args.out_prefix + ".paf"
    if not os.path.exists(paf):
        run_minimap2(args.minimap2, args.target, args.query, paf, args.threads, args.approx)
    rows = read_paf(paf)
    plot(rows, args.out_prefix + ".png", tname, qname)
    if not rows:
        print(f"{qname} vs {tname}: no alignment >= 5 kb")
        return
    tlen = rows[0][7]
    aligned = sum(te - ts for _, _, ts, te, _, _, _, _ in rows)
    print(f"{qname} vs {tname}: {len(rows)} blocks, {aligned / 1e6:.1f} Mb aligned of {tlen / 1e6:.1f} Mb target")
    print("target_bin_Mb\tbest_identity")
    for start, ident in identity_profile(rows, tlen, args.bin):
        print(f"{start / 1e6:.0f}\t{ident:.4f}" if ident else f"{start / 1e6:.0f}\t-")


if __name__ == "__main__":
    main()
