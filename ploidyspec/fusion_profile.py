"""
Divergence along fused chromosomes.

Fusion detection (rediploidization.py) says which copies are fused and how
far the fused lineage has diverged overall (`lineage_divergence`). Xie et al.
(2026) found more: in snow carps divergence between fused and unfused copies
is highest at the fusion site and spreads outwards along the arms, because
recombination with the unfused copies stops first next to the junction. This
reads that profile.

Each fused sequence (a long placed copy, or an unplaced scaffold) is windowed
along its own coordinates (position-free, windowed.py) against:
- every unfused copy of each component chromosome: per window the closest one
  gives the window's component and `d_unfused`, the fused-vs-unfused lineage
  divergence there;
- the other fused copies of the same fusion, if any (`d_fused`): the
  within-lineage distance, as a baseline.

The junction is the window boundary that best separates the two components'
windows (most windows on the right side of it). Per fusion: the median
`d_unfused` in the windows nearest the junction (the inner 20% of each arm)
and at the arm ends (the outer 20%), and their ratio `junction_gradient`:
above 1 means divergence is highest at the fusion site, as Xie et al.
describe; about 1 means it is uniform (an old fusion whose divergence has
spread along the whole chromosome, or the simulated fusions, diverged evenly).

Writes rediploidization/fusion_profile.tsv (one row per window) and
fusion_profile_<scaffold>.png, and returns per-fusion summaries for
fusions.tsv.
"""

import os
import shutil
import statistics

from .common import log
from .kmer_tables import build_one, chrom_fasta_path, ktab_prefix_path

PROFILE_FIELDS = ["scaffold", "win_start", "win_end", "component", "d_unfused", "d_fused"]
ARM_FRAC = 0.2


def window_distances(fastk_bin, profex_bin, k, outdir, fused, targets, window, tmp_root):
    """{target_unit_id: [(start, end, distance)]} along the fused sequence."""
    from .windowed import containment_distance, make_windows, profile_bins, window_counts

    fasta = chrom_fasta_path(outdir, fused["unit_id"])
    os.makedirs(tmp_root, exist_ok=True)
    windows = make_windows(int(fused["length"]) - k + 1, window, window)
    valid = window_counts(profile_bins(fastk_bin, profex_bin, k, fasta,
                                       ktab_prefix_path(outdir, fused["unit_id"], k), tmp_root, "self",
                                       window), windows, window)
    out = {}
    for t in targets:
        shared = window_counts(profile_bins(fastk_bin, profex_bin, k, fasta,
                                            ktab_prefix_path(outdir, t["unit_id"], k), tmp_root,
                                            "t", window), windows, window)
        out[t["unit_id"]] = [(s, e, containment_distance(m, n, k)[1])
                             for (s, e), n, m in zip(windows, valid, shared) if n > 0]
    shutil.rmtree(tmp_root, ignore_errors=True)
    return out


def junction(labels):
    """Index i such that labels[:i] are mostly one component and labels[i:]
    the other; returns (i, left_component, right_component, agreement)."""
    comps = sorted({c for c in labels if c is not None})
    if len(comps) < 2 or len(labels) < 4:
        return None
    best = None
    for a in comps:
        for b in comps:
            if a == b:
                continue
            for i in range(1, len(labels)):
                agree = sum(c == a for c in labels[:i]) + sum(c == b for c in labels[i:])
                if best is None or agree > best[3]:
                    best = (i, a, b, agree)
    i, a, b, agree = best
    return i, a, b, agree / len(labels)


def profile_fusion(fused, components, unfused, fused_sibs, tools, k, window, outdir):
    """Window rows and a summary dict for one fused sequence."""
    samtools_bin, fastk_bin, _, _ = tools
    from .windowed import find_profex

    profex_bin = find_profex(fastk_bin)
    for u in [fused] + unfused + fused_sibs:
        build_one(samtools_bin, fastk_bin, u, k, outdir)
    tmp = os.path.join(outdir, "tmp_rediploidization", f"profile_{fused['unit_id']}")
    dist = window_distances(fastk_bin, profex_bin, k, outdir, fused, unfused + fused_sibs, window, tmp)
    chrom_of = {u["unit_id"]: int(u["chrom"]) for u in unfused}
    starts = sorted({s for track in dist.values() for s, _, _ in track})
    by_target = {t: {s: (e, d) for s, e, d in track} for t, track in dist.items()}
    rows = []
    for s in starts:
        per_comp = {}
        for t, track in by_target.items():
            if t in chrom_of and s in track:
                c = chrom_of[t]
                per_comp[c] = min(per_comp.get(c, 1.0), track[s][1])
        if not per_comp:
            continue
        comp = min(per_comp, key=per_comp.get)
        d_f = [by_target[u["unit_id"]][s][1] for u in fused_sibs if s in by_target.get(u["unit_id"], {})]
        end = next(track[s][0] for track in by_target.values() if s in track)
        rows.append(dict(scaffold=fused["seq_id"] if fused["unit_id"].startswith("ORPHAN_") else fused["unit_id"],
                         win_start=s, win_end=end, component=f"chr{comp:02d}",
                         d_unfused=per_comp[comp], d_fused=min(d_f) if d_f else None))
    summary = dict(junction_mb="", junction_agreement="", near_junction_div="", arm_end_div="",
                   junction_gradient="", fused_within_div="")
    j = junction([int(r["component"][3:]) for r in rows])
    if j:
        i, a, b, agree = j
        left, right = rows[:i], rows[i:]
        n_l, n_r = max(1, int(ARM_FRAC * len(left))), max(1, int(ARM_FRAC * len(right)))
        near = [r["d_unfused"] for r in left[-n_l:] + right[:n_r]]
        far = [r["d_unfused"] for r in left[:n_l] + right[-n_r:]]
        near_m, far_m = statistics.median(near), statistics.median(far)
        within = [r["d_fused"] for r in rows if r["d_fused"] is not None]
        summary = dict(junction_mb=f"{rows[i]['win_start'] / 1e6:.1f}", junction_agreement=f"{agree:.2f}",
                       near_junction_div=f"{near_m:.4f}", arm_end_div=f"{far_m:.4f}",
                       junction_gradient=f"{near_m / far_m:.2f}" if far_m > 0 else "",
                       fused_within_div=f"{statistics.median(within):.4f}" if within else "")
    return rows, summary


def plot_profile(path, rows, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 3.2))
    comps = sorted({r["component"] for r in rows})
    for c in comps:
        xs = [(r["win_start"] + r["win_end"]) / 2e6 for r in rows if r["component"] == c]
        ys = [r["d_unfused"] for r in rows if r["component"] == c]
        ax.plot(xs, ys, ".", markersize=3, label=f"vs unfused {c}")
    fx = [(r["win_start"] + r["win_end"]) / 2e6 for r in rows if r["d_fused"] is not None]
    fy = [r["d_fused"] for r in rows if r["d_fused"] is not None]
    if fx:
        ax.plot(fx, fy, ".", color="0.6", markersize=3, label="vs other fused copy")
    ax.set_xlabel("position along the fused sequence (Mb)")
    ax.set_ylabel("window distance (-ln c / k)")
    ax.set_ylim(bottom=0)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fusion_profiles(outdir, rdir, units, fusions, fused_units, tools, k, window):
    """fused_units: {scaffold name as in fusions.tsv: unit dict}. Adds the
    profile summary to each fusion row in place and writes the window table."""
    fused_ids = {u["unit_id"] for u in fused_units.values()}
    confirmed = [f for f in fusions if f["status"] == "fusion" and f["scaffold"] in fused_units]
    all_rows = []
    for f in confirmed:
        comps = [int(c[3:]) for c in f["components"].split("+")]
        fused = fused_units[f["scaffold"]]
        unfused = [u for u in units if int(u["chrom"]) in comps and u["unit_id"] not in fused_ids]
        sibs = [fused_units[g["scaffold"]] for g in confirmed
                if g is not f and g["components"] == f["components"]]
        if not unfused:
            continue
        log(f"  fusion profile: {f['scaffold']} ({f['components']}) against {len(unfused)} unfused and "
            f"{len(sibs)} fused cop(ies)")
        rows, summary = profile_fusion(fused, comps, unfused, sibs, tools, k, window, outdir)
        f.update(summary)
        all_rows += rows
        if rows:
            safe = "".join(ch if ch.isalnum() else "_" for ch in f["scaffold"])
            plot_profile(os.path.join(rdir, f"fusion_profile_{safe}.png"), rows,
                         f"{f['scaffold']} ({f['components']}): fused vs unfused lineage along the fusion"
                         + (f"; junction {summary['junction_mb']} Mb, gradient {summary['junction_gradient']}"
                            if summary["junction_mb"] else ""))
    return all_rows
