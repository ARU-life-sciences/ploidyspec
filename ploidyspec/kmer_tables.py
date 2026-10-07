import csv
import glob
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from .common import log, run


def load_sequences(seq_tsv):
    with open(seq_tsv) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def chrom_fasta_path(outdir, unit_id):
    return os.path.join(outdir, "chroms", unit_id + ".fa")


def ktab_prefix_path(outdir, unit_id, k):
    return os.path.join(outdir, f"ktabs_k{k}", unit_id)


def unit_signature(unit):
    """Where a unit's cached FASTA came from: source file (with its size and
    modification time, so a rewritten assembly is re-extracted even when names
    and lengths are unchanged), sequence name and length."""
    try:
        st = os.stat(unit["source"])
        stamp = f"{st.st_size}\t{st.st_mtime_ns}"
    except OSError:
        stamp = "missing"
    return f"{unit['source']}\t{stamp}\t{unit['seq_id']}\t{int(unit['length'])}\n"


def cached_fasta_matches(chrom_fa, unit):
    """
    Whether an existing per-unit FASTA was extracted from the sequence `unit`
    now describes. Cached files are reused by name (HAP1_chr01.fa), so without
    this check a re-run on a new assembly, or after chromosome renumbering,
    silently keeps the old sequence: drLytSali1's 2026-09 re-run computed every
    distance from its previous assembly this way. Checked against the
    `.src` signature written at extraction; files from before signatures
    existed are checked by header name and by length (file size given the line
    width).
    """
    sig_path = chrom_fa + ".src"
    if os.path.exists(sig_path):
        with open(sig_path) as f:
            return f.read() == unit_signature(unit)
    with open(chrom_fa) as f:
        header = f.readline()
        width = len(f.readline().rstrip("\n"))
    if not header.startswith(">") or header[1:].split()[0] != unit["seq_id"] or width == 0:
        return False
    length = int(unit["length"])
    return os.path.getsize(chrom_fa) == len(header) + length + (length + width - 1) // width


def remove_unit_cache(outdir, unit_id):
    """Delete a unit's FASTA and its k-mer tables at every k."""
    chrom_fa = chrom_fasta_path(outdir, unit_id)
    for path in (chrom_fa, chrom_fa + ".fai", chrom_fa + ".src"):
        if os.path.exists(path):
            os.remove(path)
    for kdir in glob.glob(os.path.join(outdir, "ktabs_k*")):
        for name in os.listdir(kdir):
            if name.startswith((unit_id + ".", "." + unit_id + ".")):
                os.remove(os.path.join(kdir, name))


def build_one(samtools_bin, fastk_bin, unit, k, outdir):
    chrom_fa = chrom_fasta_path(outdir, unit["unit_id"])
    ktab_prefix = ktab_prefix_path(outdir, unit["unit_id"], k)
    os.makedirs(os.path.dirname(chrom_fa), exist_ok=True)
    os.makedirs(os.path.dirname(ktab_prefix), exist_ok=True)

    if os.path.exists(chrom_fa) and not cached_fasta_matches(chrom_fa, unit):
        log(f"  {unit['unit_id']}: cached sequence is not {unit['seq_id']} from "
            f"{os.path.basename(unit['source'])} -- rebuilding its FASTA and k-mer tables; "
            "re-run the stages after `kmers` for this output directory")
        remove_unit_cache(outdir, unit["unit_id"])

    if not os.path.exists(chrom_fa):
        tmp = chrom_fa + ".tmp"
        with open(tmp, "w") as out:
            run([samtools_bin, "faidx", unit["source"], unit["seq_id"]], stdout=out)
        os.replace(tmp, chrom_fa)
        with open(chrom_fa + ".src", "w") as f:
            f.write(unit_signature(unit))

    if not os.path.exists(ktab_prefix + ".ktab"):
        run(
            [
                fastk_bin,
                f"-k{k}",
                "-t1",
                "-T1",
                f"-N{ktab_prefix}",
                f"-P{os.path.dirname(ktab_prefix)}",
                chrom_fa,
            ]
        )

    return unit["unit_id"]


def build_all(seq_tsv, outdir, samtools_bin, fastk_bin, k, threads):
    units = load_sequences(seq_tsv)
    log(
        f"building k={k} k-mer tables for {len(units)} chromosome-scale units ({threads} workers)"
    )
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = {
            ex.submit(build_one, samtools_bin, fastk_bin, u, k, outdir): u["unit_id"]
            for u in units
        }
        done = 0
        for fut in as_completed(futs):
            uid = futs[fut]
            fut.result()
            done += 1
            log(f"  [{done}/{len(units)}] {uid}")
    return units
