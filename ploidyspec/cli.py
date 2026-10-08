#!/usr/bin/env python3
"""
ploidyspec: k-mer based subgenome/haplotype divergence profiling for polyploid genomes.

Pipeline:
  prepare   -> parse a manifest of haplotype assemblies into per-chromosome units (sequences.tsv)
  kmers     -> FastK canonical k-mer table per chromosome-scale unit
  matrix    -> all-vs-all whole-chromosome k-mer Jaccard distance matrix (Logex + Histex)
  homeologs -> detect candidate ancestral (paleopolyploid) chromosome pairs from the
               whole-chromosome matrix, e.g. chr01<->chr05 retained homeology from an
               old whole-genome duplication, distinct from same-numbered haplotype copies
  windowed  -> sliding-window k-mer divergence between haplotype copies of the same chromosome
  report    -> self-contained HTML report for one species, from whatever stages have run
  all       -> run prepare -> kmers -> matrix -> homeologs -> windowed -> report in order
               (cheap stages only by default; --with-te-markers/--with-te-markers-windowed/
               --with-windowed-homeologs opt into the more expensive stages below)

  windowed-homeologs -> sliding-window k-mer divergence between candidate ancestral pairs

  te-markers         -> differential high-copy (fossil-TE) k-mer markers between same-chromosome
                         haplotype copies -- subgenome resolver for recent allopolyploids
  te-markers-windowed -> paints each haplotype copy's windows by which subgenome's markers
                         they match -- the phaser
  subgenome-report   -> consolidates te-markers output into a continuous auto<->allo index
                         and subgenome-assignment summary
  rediploidization   -> chromosome fusions between haplotype copies, lineage structure along
                         each chromosome, and a per-chromosome rediploidization state
                         (run by `all`; uses te-markers output if present)

  structure          -> inheritance-mode metrics (partition consistency, divergence-depth
                         CVs) and diffuse genome-wide chromosome partitions (run by `all`)
  panel              -> cross-species tables and poly-space PCA from a results directory
  simulate           -> synthetic haplotype assemblies with a known answer (diploid,
                         autotetraploid x2 layouts, allotetraploid, rediploidized)
"""

import argparse
import csv
import os
import sys

from .common import (
    cleanup_intermediates,
    find_tool,
    homeologs_dir,
    log,
    matrix_dir,
    subgenomes_dir,
    windowed_dir,
)
from .manifest import DEFAULT_CHROM_REGEXES, DEFAULT_HAP_REGEX, prepare
from .kmer_tables import build_all
from .whole_matrix import compute_matrix
from .windowed import compute_windowed, compute_windowed_homeologs
from .homeologs import DEFAULT_MIN_EFFECT, run as run_homeolog_detection
from .rediploidization import (
    DEFAULT_CONTAINMENT_Z,
    DEFAULT_DIST_SPLIT,
    DEFAULT_LONG_RATIO,
    DEFAULT_MIN_SEGMENT_BP,
    DEFAULT_ORPHAN_MIN_FRAC,
    DEFAULT_PARTITION_Z,
    DEFAULT_SIBLING_EXCESS,
    DEFAULT_TE_SPLIT,
    DEFAULT_WINDOW_SPLIT,
    compute_rediploidization,
)
from .report import generate_report
from .summary import compute_summary
from .simulate import DEFAULTS as SIM_DEFAULTS, SCENARIOS, simulate
from .structure import compute_structure
from .panel import build_panel
from .subgenome_report import compute_subgenome_report
from .te_markers import (
    DEFAULT_MARKER_K,
    DEFAULT_MIN_COUNT,
    DEFAULT_MIN_RATIO,
    compute_te_markers,
    compute_te_markers_windowed,
)


def parse_k_list(s):
    try:
        values = sorted({int(x) for x in s.split(",")})
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--k must be a comma-separated list of ints, got {s!r}"
        )
    if not values:
        raise argparse.ArgumentTypeError("--k must contain at least one value")
    return values


def add_common_args(p):
    p.add_argument(
        "--manifest",
        required=True,
        help="TSV: <fasta_path>\\t<hap_label|AUTO> per line",
    )
    p.add_argument(
        "--outdir", required=True, help="output directory (created if missing)"
    )
    p.add_argument(
        "--k",
        type=parse_k_list,
        default=[15, 23],
        help="k-mer size, or comma-separated list for a multi-k sweep (default 15,23). "
        "The distance comes from the largest k that clears the chance-collision floor; "
        "on the 45-species panel that was k=23 for every pair, so 15 is a fallback for "
        "very divergent pairs. The panel itself ran 11,13,15,17,19,23 (identical "
        "distances, ~3x the k-mer counting). "
        "Used by kmers/matrix/homeologs -- matrix/homeologs combine the swept k's into "
        "a Mash-corrected consensus distance per pair (see ploidyspec/mash.py), flagging "
        "pairs whose shared-kmer count never clears the chance-collision noise floor as "
        "resolution_limited instead of reporting a number that's indistinguishable from "
        "noise. windowed/windowed-homeologs ignore this and use --window-k instead.",
    )
    p.add_argument(
        "--min-len",
        type=int,
        default=1_000_000,
        help="minimum sequence length to treat as chromosome-scale (default 1e6)",
    )
    p.add_argument(
        "--chrom-regex",
        action="append",
        default=None,
        help="regex with one capturing group for the chromosome number, tried against the FASTA "
        "description in order given; may be passed multiple times (default: built-in patterns "
        "for 'chromosome: N' and 'SUPER_N' style headers)",
    )
    p.add_argument(
        "--hap-regex",
        default=None,
        help=f"regex with one capturing group for the haplotype number, used only for manifest rows "
        f"marked AUTO (default: {DEFAULT_HAP_REGEX!r})",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count() or 4,
        help="parallel workers (default: all cores)",
    )
    p.add_argument(
        "--samtools", default=None, help="path to samtools binary or its containing dir"
    )
    p.add_argument(
        "--fastk-dir", default=None, help="dir containing FastK/Logex/Histex binaries"
    )


def add_window_k_arg(p):
    p.add_argument(
        "--window-k",
        type=int,
        default=15,
        help="k-mer size for the windowed track (default 15). Decoupled from --k: "
        "windowed stages build their own per-window k-mer tables independently and "
        "don't currently support a multi-k sweep or Mash correction.",
    )


def add_marker_args(p):
    p.add_argument(
        "--marker-k",
        type=int,
        default=DEFAULT_MARKER_K,
        help=f"k-mer size for differential fossil-TE marker extraction (default "
        f"{DEFAULT_MARKER_K}, matching the Jaron/Cerca method). Decoupled from --k "
        f"and --window-k: builds its own .ktab per unit if not already present.",
    )
    p.add_argument(
        "--min-count",
        type=int,
        default=DEFAULT_MIN_COUNT,
        help=f"minimum within-chromosome k-mer count to treat as high-copy/repetitive "
        f"(default {DEFAULT_MIN_COUNT}, matching the tutorial's validated choice). "
        f"Check `Histex -h` on a unit's .hist if unsure this falls in a real high-copy "
        f"tail for a given species.",
    )
    p.add_argument(
        "--min-ratio",
        type=float,
        default=DEFAULT_MIN_RATIO,
        help=f"minimum count ratio between haplotype copies for a high-copy k-mer to "
        f"be called a differential marker (default {DEFAULT_MIN_RATIO}).",
    )


def add_rediploidization_args(p):
    p.add_argument("--long-ratio", type=float, default=DEFAULT_LONG_RATIO,
                   help="a placed copy at >= this x its siblings' median length is tested as a "
                   f"candidate fusion (default {DEFAULT_LONG_RATIO})")
    p.add_argument("--orphan-min-frac", type=float, default=DEFAULT_ORPHAN_MIN_FRAC,
                   help="unplaced scaffolds >= this x the median chromosome length are tested as "
                   f"candidate fused chromosomes (default {DEFAULT_ORPHAN_MIN_FRAC})")
    p.add_argument("--containment-z", type=float, default=DEFAULT_CONTAINMENT_Z,
                   help="robust z above the species' background k-mer containment for a chromosome "
                   f"to count as a fusion component (default {DEFAULT_CONTAINMENT_Z})")
    p.add_argument("--sibling-excess", type=float, default=DEFAULT_SIBLING_EXCESS,
                   help="a long copy's fusion partner must be this many times more enriched in it "
                   "than in its least-enriched sibling copy, per k-mer (default "
                   f"{DEFAULT_SIBLING_EXCESS}; rules out ancient homeology and fragmentary siblings)")
    p.add_argument("--partition-z", type=float, default=DEFAULT_PARTITION_Z,
                   help="with no homeolog pairs, pool < 3-copy chromosomes across the structure "
                   "stage's k=2 genome partition when its z >= this "
                   f"(default {DEFAULT_PARTITION_Z})")
    p.add_argument("--dist-split", type=float, default=DEFAULT_DIST_SPLIT,
                   help="cross/within whole-chromosome distance ratio counted as lineage-split "
                   f"evidence (default {DEFAULT_DIST_SPLIT}; 2x this plus a whole-chromosome "
                   "windowed split = resolved_lineages)")
    p.add_argument("--te-split", type=float, default=DEFAULT_TE_SPLIT,
                   help=f"TE-marker split_ratio counted as lineage-split evidence (default {DEFAULT_TE_SPLIT})")
    p.add_argument("--window-split", type=float, default=DEFAULT_WINDOW_SPLIT,
                   help="per-window cross/within distance ratio counted as a split window "
                   f"(default {DEFAULT_WINDOW_SPLIT})")
    p.add_argument("--min-segment-bp", type=int, default=DEFAULT_MIN_SEGMENT_BP,
                   help="minimum length of a contiguous split region to count "
                   f"(default {DEFAULT_MIN_SEGMENT_BP})")


def resolve_tools(args):
    samtools_bin = find_tool(args.samtools, "samtools*", "samtools")
    fastk_bin = find_tool(args.fastk_dir, "FASTK*", "FastK")
    logex_bin = find_tool(args.fastk_dir, "FASTK*", "Logex")
    histex_bin = find_tool(args.fastk_dir, "FASTK*", "Histex")
    tabex_bin = find_tool(args.fastk_dir, "FASTK*", "Tabex")
    # FastK shells out to sibling tools (e.g. Fastrm) by bare name on error-cleanup paths.
    for d in {os.path.dirname(p) for p in (samtools_bin, fastk_bin)}:
        if d not in os.environ.get("PATH", "").split(os.pathsep):
            os.environ["PATH"] = d + os.pathsep + os.environ.get("PATH", "")
    return samtools_bin, fastk_bin, logex_bin, histex_bin, tabex_bin


def cmd_prepare(args):
    samtools_bin, _, _, _, _ = resolve_tools(args)
    prepare(
        args.manifest,
        args.outdir,
        samtools_bin,
        args.min_len,
        args.chrom_regex,
        args.hap_regex,
    )


def cmd_kmers(args):
    samtools_bin, fastk_bin, _, _, _ = resolve_tools(args)
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        raise SystemExit(f"{seq_tsv} not found -- run the `prepare` stage first")
    for k in args.k:
        build_all(seq_tsv, args.outdir, samtools_bin, fastk_bin, k, args.threads)


def cmd_matrix(args):
    _, _, logex_bin, histex_bin, _ = resolve_tools(args)
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        raise SystemExit(
            f"{seq_tsv} not found -- run the `prepare` and `kmers` stages first"
        )
    compute_matrix(seq_tsv, args.outdir, logex_bin, histex_bin, args.threads, args.k)
    mdir = matrix_dir(args.outdir)
    log(
        f"wrote {os.path.join(mdir, 'whole_chrom_distance_matrix.csv')}, "
        f"{os.path.join(mdir, 'whole_chrom_distance_heatmap.png')}, "
        f"ploidy_summary.tsv and homologous_chromosomes.tsv in {mdir}"
    )


def cmd_windowed(args):
    samtools_bin, fastk_bin, logex_bin, histex_bin, _ = resolve_tools(args)
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        raise SystemExit(
            f"{seq_tsv} not found -- run the `prepare` and `kmers` stages first"
        )
    step = args.step or args.window
    compute_windowed(
        seq_tsv,
        args.outdir,
        samtools_bin,
        fastk_bin,
        logex_bin,
        histex_bin,
        args.window_k,
        args.window,
        step,
        args.threads,
    )
    log(
        f"wrote per-chromosome windowed_chrNN.tsv/.png, windowed_all.tsv and "
        f"windowed_genome_overview.png in {windowed_dir(args.outdir)}"
    )


def cmd_homeologs(args):
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    matrix_csv = os.path.join(matrix_dir(args.outdir), "whole_chrom_distance_matrix.csv")
    if not os.path.exists(seq_tsv) or not os.path.exists(matrix_csv):
        raise SystemExit(
            f"need sequences.tsv and {matrix_csv} -- run `prepare`, `kmers`, `matrix` first"
        )
    run_homeolog_detection(seq_tsv, args.outdir, args.fdr_alpha, args.min_effect)
    log(
        f"wrote {os.path.join(homeologs_dir(args.outdir), 'homeolog_pairs.tsv')} and homeolog_pairs.png"
    )


def cmd_windowed_homeologs(args):
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    pairs_tsv = os.path.join(homeologs_dir(args.outdir), "homeolog_pairs.tsv")
    if not os.path.exists(seq_tsv):
        raise SystemExit(
            f"{seq_tsv} not found -- run the `prepare` and `kmers` stages first"
        )
    if not os.path.exists(pairs_tsv):
        raise SystemExit(f"{pairs_tsv} not found -- run the `homeologs` stage first")
    with open(pairs_tsv) as f:
        pairs = [
            (
                int(row["chrom_a"].replace("chr", "")),
                int(row["chrom_b"].replace("chr", "")),
            )
            for row in csv.DictReader(f, delimiter="\t")
        ]
    if not pairs:
        # not an error: many genomes have no accepted homeolog pairs, and `all`
        # must carry on to the later stages (dcCerAlpi1 stopped here)
        log(f"{pairs_tsv} has no accepted homeolog pairs -- skipping windowed-homeologs")
        return
    samtools_bin, fastk_bin, logex_bin, histex_bin, _ = resolve_tools(args)
    step = args.step or args.window
    compute_windowed_homeologs(
        seq_tsv,
        args.outdir,
        pairs,
        samtools_bin,
        fastk_bin,
        logex_bin,
        histex_bin,
        args.window_k,
        args.window,
        step,
        args.threads,
    )
    log(
        f"wrote windowed_chrAAxBB.tsv/.png per pair, windowed_homeologs_all.tsv and "
        f"windowed_homeologs_overview.png in {homeologs_dir(args.outdir)}"
    )


def cmd_te_markers(args):
    samtools_bin, fastk_bin, _, _, tabex_bin = resolve_tools(args)
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        raise SystemExit(
            f"{seq_tsv} not found -- run the `prepare` and `kmers` stages first"
        )
    compute_te_markers(
        seq_tsv,
        args.outdir,
        samtools_bin,
        fastk_bin,
        tabex_bin,
        args.marker_k,
        args.min_count,
        args.min_ratio,
        args.threads,
    )
    log(
        f"wrote te_markers_<unit_a>x<unit_b>.tsv per same-chromosome haplotype pair "
        f"and te_markers_summary.tsv in {subgenomes_dir(args.outdir)}"
    )


def cmd_te_markers_windowed(args):
    samtools_bin, fastk_bin, _, _, tabex_bin = resolve_tools(args)
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        raise SystemExit(
            f"{seq_tsv} not found -- run the `prepare` and `kmers` stages first"
        )
    step = args.step or args.window
    compute_te_markers_windowed(
        seq_tsv,
        args.outdir,
        samtools_bin,
        fastk_bin,
        tabex_bin,
        args.marker_k,
        args.min_ratio,
        args.window,
        step,
        args.threads,
    )
    log(
        f"wrote te_markers_windowed_<unit_a>x<unit_b>.tsv/.png per pair in "
        f"{subgenomes_dir(args.outdir)}"
    )


def cmd_subgenome_report(args):
    compute_subgenome_report(args.outdir)


def cmd_rediploidization(args):
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    pairs_tsv = os.path.join(matrix_dir(args.outdir), "whole_chrom_pairs.tsv")
    if not os.path.exists(seq_tsv) or not os.path.exists(pairs_tsv):
        raise SystemExit(
            f"need sequences.tsv and {pairs_tsv} -- run `prepare`, `kmers`, `matrix` first"
        )
    thresholds = dict(
        long_ratio=args.long_ratio,
        orphan_min_frac=args.orphan_min_frac,
        containment_z=args.containment_z,
        sibling_excess=args.sibling_excess,
        partition_z=args.partition_z,
        dist_split=args.dist_split,
        te_split=args.te_split,
        window_split=args.window_split,
        min_segment_bp=args.min_segment_bp,
    )

    def tools():
        samtools_bin, fastk_bin, logex_bin, histex_bin, _ = resolve_tools(args)
        return samtools_bin, fastk_bin, logex_bin, histex_bin

    compute_rediploidization(
        args.outdir, args.k, args.min_len, thresholds, tools, args.threads
    )


def cmd_structure(args):
    # every metric degrades to blank when its input (homeolog pairs, windowed
    # tracks) is missing, so this never blocks `all` on a species with too few
    # chromosomes for the homeolog search
    compute_structure(args.outdir)


def cmd_summary(args):
    compute_summary(args.outdir)


def cmd_panel(args):
    build_panel(args.results, args.outdir, args.categories)


def cmd_cleanup(args):
    removed = cleanup_intermediates(args.outdir)
    log(f"removed {len(removed)} intermediate dir(s) from {args.outdir}"
        + "".join(f"\n  {p}" for p in removed))


def cmd_simulate(args):
    params = dict(n_chrom=args.n_chrom, chrom_len=args.chrom_len)
    for d in simulate(args.scenario, args.outdir, args.seed, params):
        log(f"wrote {d}/manifest.tsv, truth.tsv and one FASTA per haplotype")


def cmd_report(args):
    path = generate_report(args.outdir)
    log(f"wrote {path}")


def cmd_all(args):
    cmd_prepare(args)
    cmd_kmers(args)
    cmd_matrix(args)
    cmd_homeologs(args)
    cmd_windowed(args)
    if args.with_windowed_homeologs:
        cmd_windowed_homeologs(args)
    cmd_structure(args)
    run_te_markers = args.with_te_markers or args.with_te_markers_windowed
    if run_te_markers:
        cmd_te_markers(args)
    if args.with_te_markers_windowed:
        cmd_te_markers_windowed(args)
    if run_te_markers:
        cmd_subgenome_report(args)
    cmd_rediploidization(args)
    cmd_summary(args)
    cmd_report(args)
    if args.cleanup:
        cmd_cleanup(args)


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="ploidyspec",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="stage", required=True)

    p_prep = sub.add_parser(
        "prepare", help="build sequences.tsv from a manifest of haplotype assemblies"
    )
    add_common_args(p_prep)
    p_prep.set_defaults(func=cmd_prepare)

    p_kmers = sub.add_parser("kmers", help="build per-chromosome FastK k-mer tables")
    add_common_args(p_kmers)
    p_kmers.set_defaults(func=cmd_kmers)

    p_matrix = sub.add_parser(
        "matrix", help="whole-chromosome all-vs-all k-mer distance matrix"
    )
    add_common_args(p_matrix)
    p_matrix.set_defaults(func=cmd_matrix)

    p_win = sub.add_parser(
        "windowed",
        help="sliding-window k-mer divergence between haplotype copies of each chromosome",
    )
    add_common_args(p_win)
    add_window_k_arg(p_win)
    p_win.add_argument(
        "--window", type=int, default=250_000, help="window size in bp (default 250000)"
    )
    p_win.add_argument(
        "--step",
        type=int,
        default=None,
        help="step size in bp (default: same as --window, i.e. tumbling windows)",
    )
    p_win.set_defaults(func=cmd_windowed)

    p_all = sub.add_parser(
        "all",
        help="run prepare -> kmers -> matrix -> homeologs -> windowed -> structure -> "
        "rediploidization -> summary -> report in sequence (see --with-* flags for the opt-in stages)",
    )
    add_common_args(p_all)
    add_window_k_arg(p_all)
    add_marker_args(p_all)
    add_rediploidization_args(p_all)
    p_all.add_argument(
        "--window", type=int, default=250_000, help="window size in bp (default 250000)"
    )
    p_all.add_argument(
        "--step",
        type=int,
        default=None,
        help="step size in bp (default: same as --window, i.e. tumbling windows)",
    )
    p_all.add_argument(
        "--fdr-alpha",
        type=float,
        default=0.05,
        help="Benjamini-Hochberg FDR threshold for the homeologs stage (default 0.05)",
    )
    p_all.add_argument(
        "--min-effect",
        type=float,
        default=DEFAULT_MIN_EFFECT,
        help="a homeolog pair must also be at least this fraction closer than the median "
        f"cross-chromosome distance (default {DEFAULT_MIN_EFFECT}; guards against "
        "significant-but-tiny differences when the background is very tight)",
    )
    p_all.add_argument(
        "--with-windowed-homeologs",
        action="store_true",
        help="also run windowed-homeologs after homeologs (moderate cost, scales with "
        "how many ancestral pairs are found)",
    )
    p_all.add_argument(
        "--with-te-markers",
        action="store_true",
        help="also run te-markers + subgenome-report (opt-in: moderate cost, reuses "
        "existing k-mer tables where possible)",
    )
    p_all.add_argument(
        "--with-te-markers-windowed",
        action="store_true",
        help="also run te-markers-windowed (implies --with-te-markers) + "
        "subgenome-report (opt-in: the most expensive stage -- builds a fresh k-mer "
        "table per window)",
    )
    p_all.add_argument(
        "--cleanup",
        action="store_true",
        help="delete k-mer tables, per-chromosome FASTAs and temp dirs after the run "
        "(keeps every table and plot; re-running a stage later rebuilds what it needs)",
    )
    p_all.set_defaults(func=cmd_all)

    p_homeo = sub.add_parser(
        "homeologs",
        help="detect candidate ancestral (paleopolyploid) chromosome pairs from the whole-chromosome matrix",
    )
    add_common_args(p_homeo)
    p_homeo.add_argument(
        "--fdr-alpha",
        type=float,
        default=0.05,
        help="Benjamini-Hochberg FDR threshold for accepting a candidate ancestral "
        "chromosome pair (default 0.05). Empirical p-value per pair = fraction of "
        "all cross-chromosome-number distances that are as small or smaller.",
    )
    p_homeo.add_argument(
        "--min-effect",
        type=float,
        default=DEFAULT_MIN_EFFECT,
        help="a homeolog pair must also be at least this fraction closer than the median "
        f"cross-chromosome distance (default {DEFAULT_MIN_EFFECT}; guards against "
        "significant-but-tiny differences when the background is very tight)",
    )
    p_homeo.set_defaults(func=cmd_homeologs)

    p_win_homeo = sub.add_parser(
        "windowed-homeologs",
        help="sliding-window k-mer divergence between candidate ancestral chromosome pairs found by `homeologs`",
    )
    add_common_args(p_win_homeo)
    add_window_k_arg(p_win_homeo)
    p_win_homeo.add_argument(
        "--window", type=int, default=250_000, help="window size in bp (default 250000)"
    )
    p_win_homeo.add_argument(
        "--step",
        type=int,
        default=None,
        help="step size in bp (default: same as --window)",
    )
    p_win_homeo.set_defaults(func=cmd_windowed_homeologs)

    p_te = sub.add_parser(
        "te-markers",
        help="differential high-copy (fossil-TE) k-mer markers between same-chromosome "
        "haplotype copies -- subgenome resolver for recent allopolyploids",
    )
    add_common_args(p_te)
    add_marker_args(p_te)
    p_te.set_defaults(func=cmd_te_markers)

    p_te_win = sub.add_parser(
        "te-markers-windowed",
        help="paint each haplotype copy's windows by which subgenome's fossil-TE "
        "markers they match -- the phaser; run `te-markers` first",
    )
    add_common_args(p_te_win)
    add_marker_args(p_te_win)
    p_te_win.add_argument(
        "--window", type=int, default=250_000, help="window size in bp (default 250000)"
    )
    p_te_win.add_argument(
        "--step",
        type=int,
        default=None,
        help="step size in bp (default: same as --window, i.e. tumbling windows)",
    )
    p_te_win.set_defaults(func=cmd_te_markers_windowed)

    p_report = sub.add_parser(
        "report",
        help="self-contained HTML report for one species (embeds whatever plots/tables "
        "the completed stages produced; degrades gracefully if some stages haven't run)",
    )
    p_report.add_argument(
        "--outdir", required=True, help="species results directory to report on"
    )
    p_report.set_defaults(func=cmd_report)

    p_subgenome = sub.add_parser(
        "subgenome-report",
        help="consolidate te-markers/te-markers-windowed output into a continuous "
        "auto<->allo index and a subgenome-assignment summary (pure aggregation, "
        "no new FastK/Tabex calls; run `te-markers`/`te-markers-windowed` first)",
    )
    add_common_args(p_subgenome)
    p_subgenome.set_defaults(func=cmd_subgenome_report)

    p_redip = sub.add_parser(
        "rediploidization",
        help="chromosome fusions, per-chromosome lineage structure and rediploidization "
        "state (run after matrix/homeologs/windowed; uses te-markers output if present)",
    )
    add_common_args(p_redip)
    add_rediploidization_args(p_redip)
    p_redip.set_defaults(func=cmd_rediploidization)

    p_struct = sub.add_parser(
        "structure",
        help="genome-wide chromosome partitions and homeolog pair synchrony "
        "(from matrix/homeologs output, no new k-mer work)",
    )
    p_struct.add_argument("--outdir", required=True, help="species results directory")
    p_struct.set_defaults(func=cmd_structure)

    p_summ = sub.add_parser(
        "summary",
        help="plain-language answers to the four questions (ploidy, auto/allo-like "
        "structure, TE markers, rediploidization) with confidence, from the other stages' output",
    )
    p_summ.add_argument("--outdir", required=True, help="species results directory")
    p_summ.set_defaults(func=cmd_summary)

    p_panel = sub.add_parser(
        "panel",
        help="cross-species tables (auto/allo spectrum, genome partitions, rediploidization) "
        "and the poly-space PCA, from a directory of species results",
    )
    p_panel.add_argument("--results", required=True,
                         help="directory with one results subdirectory per species")
    p_panel.add_argument("--outdir", required=True, help="where to write the panel tables")
    p_panel.add_argument("--categories", default=None,
                         help="optional TSV (species, category) used only to colour the PCA plot")
    p_panel.set_defaults(func=cmd_panel)

    p_clean = sub.add_parser(
        "cleanup",
        help="delete a species' k-mer tables, per-chromosome FASTAs and temp dirs "
        "(the bulk of its disk use), keeping all outputs",
    )
    p_clean.add_argument("--outdir", required=True, help="species results directory")
    p_clean.set_defaults(func=cmd_cleanup)

    p_sim = sub.add_parser(
        "simulate",
        help="write synthetic haplotype assemblies with a known answer (manifest.tsv + "
        "truth.tsv per scenario) for testing and calibration",
    )
    p_sim.add_argument("--scenario", choices=SCENARIOS + ("all",), default="all")
    p_sim.add_argument("--outdir", required=True, help="one subdirectory per scenario")
    p_sim.add_argument("--seed", type=int, default=1)
    p_sim.add_argument("--n-chrom", type=int, default=SIM_DEFAULTS["n_chrom"],
                       help=f"chromosomes per genome (default {SIM_DEFAULTS['n_chrom']}; "
                       "allotetraploid has 2x this)")
    p_sim.add_argument("--chrom-len", type=int, default=SIM_DEFAULTS["chrom_len"],
                       help=f"chromosome length in bp (default {SIM_DEFAULTS['chrom_len']})")
    p_sim.set_defaults(func=cmd_simulate)

    args = parser.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
