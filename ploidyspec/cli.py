#!/usr/bin/env python3
"""
ploidyspec: k-mer based subgenome/haplotype divergence profiling for polyploid genomes.

Pipeline:
  prepare   -> parse a manifest of haplotype assemblies into per-chromosome units (sequences.tsv)
  kmers     -> FastK canonical k-mer table per chromosome-scale unit
  matrix    -> all-vs-all whole-chromosome k-mer Jaccard distance matrix (Logex + Histex)
  windowed  -> sliding-window k-mer divergence between haplotype copies of the same chromosome
  all       -> run the four stages above in order

  homeologs         -> detect candidate ancestral (paleopolyploid) chromosome pairs from the
                        whole-chromosome matrix, e.g. chr01<->chr05 retained homeology from an
                        old whole-genome duplication, distinct from same-numbered haplotype copies
  windowed-homeologs -> sliding-window k-mer divergence between those candidate ancestral pairs

  te-markers         -> differential high-copy (fossil-TE) k-mer markers between same-chromosome
                         haplotype copies -- subgenome resolver for recent allopolyploids
"""

import argparse
import csv
import os
import sys

from .common import find_tool, log
from .manifest import DEFAULT_CHROM_REGEXES, DEFAULT_HAP_REGEX, prepare
from .kmer_tables import build_all
from .whole_matrix import compute_matrix
from .windowed import compute_windowed, compute_windowed_homeologs
from .homeologs import run as run_homeolog_detection
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
        default=[15],
        help="k-mer size, or comma-separated list for a multi-k sweep (default 15). "
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
    log(
        f"wrote {os.path.join(args.outdir, 'whole_chrom_distance_matrix.csv')} and "
        f"{os.path.join(args.outdir, 'whole_chrom_distance_heatmap.png')}"
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
        f"wrote per-chromosome windowed_chrNN.tsv/.png, windowed_all.tsv and windowed_genome_overview.png in {args.outdir}"
    )


def cmd_homeologs(args):
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    matrix_csv = os.path.join(args.outdir, "whole_chrom_distance_matrix.csv")
    if not os.path.exists(seq_tsv) or not os.path.exists(matrix_csv):
        raise SystemExit(
            f"need sequences.tsv and whole_chrom_distance_matrix.csv -- run `prepare`, `kmers`, `matrix` first"
        )
    run_homeolog_detection(seq_tsv, args.outdir)
    log(
        f"wrote {os.path.join(args.outdir, 'homeolog_pairs.tsv')} and homeolog_pairs.png"
    )


def cmd_windowed_homeologs(args):
    samtools_bin, fastk_bin, logex_bin, histex_bin, _ = resolve_tools(args)
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    pairs_tsv = os.path.join(args.outdir, "homeolog_pairs.tsv")
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
        raise SystemExit(f"{pairs_tsv} has no candidate pairs -- nothing to do")
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
        f"wrote windowed_chrAAxBB.tsv/.png per pair, windowed_homeologs_all.tsv and windowed_homeologs_overview.png in {args.outdir}"
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
        f"and te_markers_summary.tsv in {args.outdir}"
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
        f"wrote te_markers_windowed_<unit_a>x<unit_b>.tsv/.png per pair in {args.outdir}"
    )


def cmd_all(args):
    cmd_prepare(args)
    cmd_kmers(args)
    cmd_matrix(args)
    cmd_windowed(args)


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
        "all", help="run prepare -> kmers -> matrix -> windowed in sequence"
    )
    add_common_args(p_all)
    add_window_k_arg(p_all)
    p_all.add_argument(
        "--window", type=int, default=250_000, help="window size in bp (default 250000)"
    )
    p_all.add_argument(
        "--step",
        type=int,
        default=None,
        help="step size in bp (default: same as --window, i.e. tumbling windows)",
    )
    p_all.set_defaults(func=cmd_all)

    p_homeo = sub.add_parser(
        "homeologs",
        help="detect candidate ancestral (paleopolyploid) chromosome pairs from the whole-chromosome matrix",
    )
    add_common_args(p_homeo)
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

    args = parser.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
