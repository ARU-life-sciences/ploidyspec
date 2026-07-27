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


def add_common_args(p):
    p.add_argument(
        "--manifest",
        required=True,
        help="TSV: <fasta_path>\\t<hap_label|AUTO> per line",
    )
    p.add_argument(
        "--outdir", required=True, help="output directory (created if missing)"
    )
    p.add_argument("--k", type=int, default=15, help="k-mer size (default 15)")
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


def resolve_tools(args):
    samtools_bin = find_tool(args.samtools, "samtools*", "samtools")
    fastk_bin = find_tool(args.fastk_dir, "FASTK*", "FastK")
    logex_bin = find_tool(args.fastk_dir, "FASTK*", "Logex")
    histex_bin = find_tool(args.fastk_dir, "FASTK*", "Histex")
    # FastK shells out to sibling tools (e.g. Fastrm) by bare name on error-cleanup paths.
    for d in {os.path.dirname(p) for p in (samtools_bin, fastk_bin)}:
        if d not in os.environ.get("PATH", "").split(os.pathsep):
            os.environ["PATH"] = d + os.pathsep + os.environ.get("PATH", "")
    return samtools_bin, fastk_bin, logex_bin, histex_bin


def cmd_prepare(args):
    samtools_bin, _, _, _ = resolve_tools(args)
    prepare(
        args.manifest,
        args.outdir,
        samtools_bin,
        args.min_len,
        args.chrom_regex,
        args.hap_regex,
    )


def cmd_kmers(args):
    samtools_bin, fastk_bin, _, _ = resolve_tools(args)
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        raise SystemExit(f"{seq_tsv} not found -- run the `prepare` stage first")
    build_all(seq_tsv, args.outdir, samtools_bin, fastk_bin, args.k, args.threads)


def cmd_matrix(args):
    _, _, logex_bin, histex_bin = resolve_tools(args)
    seq_tsv = os.path.join(args.outdir, "sequences.tsv")
    if not os.path.exists(seq_tsv):
        raise SystemExit(
            f"{seq_tsv} not found -- run the `prepare` and `kmers` stages first"
        )
    compute_matrix(seq_tsv, args.outdir, logex_bin, histex_bin, args.threads)
    log(
        f"wrote {os.path.join(args.outdir, 'whole_chrom_distance_matrix.csv')} and "
        f"{os.path.join(args.outdir, 'whole_chrom_distance_heatmap.png')}"
    )


def cmd_windowed(args):
    samtools_bin, fastk_bin, logex_bin, histex_bin = resolve_tools(args)
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
        args.k,
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
    samtools_bin, fastk_bin, logex_bin, histex_bin = resolve_tools(args)
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
        args.k,
        args.window,
        step,
        args.threads,
    )
    log(
        f"wrote windowed_chrAAxBB.tsv/.png per pair, windowed_homeologs_all.tsv and windowed_homeologs_overview.png in {args.outdir}"
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

    args = parser.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
