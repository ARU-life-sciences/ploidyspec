import csv
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from .common import log, run


def load_sequences(seq_tsv):
    with open(seq_tsv) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def chrom_fasta_path(outdir, unit_id):
    return os.path.join(outdir, "chroms", unit_id + ".fa")


def ktab_prefix_path(outdir, unit_id):
    return os.path.join(outdir, "ktabs", unit_id)


def build_one(samtools_bin, fastk_bin, unit, k, outdir):
    chrom_fa = chrom_fasta_path(outdir, unit["unit_id"])
    ktab_prefix = ktab_prefix_path(outdir, unit["unit_id"])
    os.makedirs(os.path.dirname(chrom_fa), exist_ok=True)
    os.makedirs(os.path.dirname(ktab_prefix), exist_ok=True)

    if not os.path.exists(chrom_fa):
        tmp = chrom_fa + ".tmp"
        with open(tmp, "w") as out:
            run([samtools_bin, "faidx", unit["source"], unit["seq_id"]], stdout=out)
        os.replace(tmp, chrom_fa)

    if not os.path.exists(ktab_prefix + ".ktab"):
        run([fastk_bin, f"-k{k}", "-t1", "-T1", f"-N{ktab_prefix}", f"-P{os.path.dirname(ktab_prefix)}", chrom_fa])

    return unit["unit_id"]


def build_all(seq_tsv, outdir, samtools_bin, fastk_bin, k, threads):
    units = load_sequences(seq_tsv)
    log(f"building k={k} k-mer tables for {len(units)} chromosome-scale units ({threads} workers)")
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = {ex.submit(build_one, samtools_bin, fastk_bin, u, k, outdir): u["unit_id"] for u in units}
        done = 0
        for fut in as_completed(futs):
            uid = futs[fut]
            fut.result()
            done += 1
            log(f"  [{done}/{len(units)}] {uid}")
    return units
