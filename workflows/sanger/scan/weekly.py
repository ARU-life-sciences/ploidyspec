#!/usr/bin/env python3
"""Weekly polyploidy scan of Tree of Life releases (report only).

1. Crawl the ToL data tree for releases with >= 2 haplotypes (crawl_tol.py).
2. Metadata flags for every release not seen before (metadata.py, GoaT).
3. Submit the k-mer screen (screen_one.sh) for up to --max releases: new
   releases first, then the backlog with metadata flags, then the rest,
   newest first. Releases already in results/ (the panel) are skipped.
4. A report job, after every screen job has ended (report.py), writes
   scan/reports/<date>.md and scan/candidates.tsv.
5. Unless --once, resubmit itself for the same time next week.

State lives in scan/ (gitignored): catalog.tsv, metadata.tsv, state.tsv
(tolid, version, status, date, job), results/<tolid>/, logs/, reports/.
Run on the farm: bsub -q transfer ... python3 workflows/sanger/scan/weekly.py
"""
import argparse
import csv
import datetime
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SCAN = os.path.join(REPO, "scan")
STATE_FIELDS = ["tolid", "version", "status", "date", "job"]
FLAG_CUT = 1.8  # karyotype / genus / size ratios at or above this are flags
BUSCO_DUP_CUT = 10.0


def read(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path, rows, fields):
    with open(path + ".part", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(path + ".part", path)


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def flags(m):
    out = []
    if int(m.get("n_haps") or 2) > 2:
        out.append(f"{m['n_haps']} haplotypes")
    for col, label in (("karyotype_ratio", "recorded n"), ("genus_ratio", "genus n"), ("size_ratio", "C-value")):
        v = num(m.get(col))
        if v is not None and v >= FLAG_CUT:
            out.append(f"{label} {v:.1f}x")
    v = num(m.get("busco_dup"))
    if v is not None and v >= BUSCO_DUP_CUT:
        out.append(f"BUSCO D {v:.0f}%")
    if num(m.get("goat_ploidy")) and num(m["goat_ploidy"]) > 2:
        out.append(f"GoaT ploidy {m['goat_ploidy']}")
    return out


def max_chrom_len(fastas):
    best = 0
    for fa in fastas:
        fai = fa + ".fai"
        if os.path.exists(fai):
            with open(fai) as f:
                for line in f:
                    p = line.split("\t")
                    if len(p) > 1 and p[1].isdigit():
                        best = max(best, int(p[1]))
    return best


def bsub(args, cmd):
    out = subprocess.run(["bsub", *args, cmd], capture_output=True, text=True, check=True).stdout
    m = re.search(r"Job <(\d+)>", out)
    return m.group(1) if m else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=int(os.environ.get("SCAN_MAX", 60)),
                    help="k-mer screens submitted per run (default 60)")
    ap.add_argument("--next-max", type=int, default=int(os.environ.get("SCAN_MAX", 60)),
                    help="--max for the resubmitted runs (default 60); lets a first run clear the backlog")
    ap.add_argument("--queue", default=os.environ.get("SCAN_QUEUE", "transfer"))
    ap.add_argument("--big-queue", default=os.environ.get("SCAN_BIG_QUEUE", "basement"),
                    help="queue for screens needing more than --queue-mem MB (transfer hosts have 42 GB)")
    ap.add_argument("--queue-mem", type=int, default=30000)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--once", action="store_true", help="do not resubmit for next week")
    ap.add_argument("--dry-run", action="store_true", help="crawl, metadata and selection only")
    a = ap.parse_args()

    today = datetime.date.today().isoformat()
    for sub in ("results", "logs", "reports"):
        os.makedirs(os.path.join(SCAN, sub), exist_ok=True)
    catalog_path = os.path.join(SCAN, "catalog.tsv")
    previous = {(r["tolid"], r["version"]) for r in read(catalog_path)}
    subprocess.run([sys.executable, os.path.join(HERE, "crawl_tol.py"), catalog_path], check=True)
    catalog = read(catalog_path)

    state_path = os.path.join(SCAN, "state.tsv")
    state = {(r["tolid"], r["version"]): r for r in read(state_path)}
    panel = set(os.listdir(os.path.join(REPO, "results")))
    todo = [r for r in catalog if (r["tolid"], r["version"]) not in state and r["tolid"] not in panel]

    meta_path = os.path.join(SCAN, "metadata.tsv")
    meta = {r["tolid"]: r for r in read(meta_path)}
    missing = [r["tolid"] for r in todo if r["tolid"] not in meta]
    if missing:
        ids = os.path.join(SCAN, "logs", f"metadata_{today}.ids")
        with open(ids, "w") as f:
            f.write("\n".join(missing) + "\n")
        new_meta = os.path.join(SCAN, "logs", f"metadata_{today}.tsv")
        subprocess.run([sys.executable, os.path.join(HERE, "metadata.py"), catalog_path, ids, new_meta,
                        os.path.join(SCAN, "goat_cache.json")], check=True)
        rows = read(new_meta)
        meta.update({r["tolid"]: r for r in rows})
        write(meta_path, list(meta.values()), list(rows[0]) if rows else [])

    def priority(r):
        new = bool(previous) and (r["tolid"], r["version"]) not in previous
        return (not new, not flags(meta.get(r["tolid"], {})), -int(r["released"].replace("-", "") or 0))

    batch = sorted(todo, key=priority)[:a.max]
    print(f"{len(catalog)} releases, {len(todo)} not yet screened, {len(batch)} submitted this run")
    if a.dry_run:
        for r in batch:
            print(r["tolid"], r["species"], r["released"], "; ".join(flags(meta.get(r["tolid"], {}))))
        return

    group = f"scan.{today}"
    for r in batch:
        d = os.path.join(SCAN, "results", r["tolid"])
        os.makedirs(d, exist_ok=True)
        fastas = r["fastas"].split(";")
        with open(os.path.join(d, "manifest.tsv"), "w") as f:
            f.write("# fasta_path\thap_label\n")
            f.writelines(f"{fa}\tHAP{i}\n" for i, fa in enumerate(fastas, 1))
        mem = max(8000, int(40 * max_chrom_len(fastas) / 1e6) + 4000)
        queue = a.queue if mem <= a.queue_mem else a.big_queue
        job = bsub(["-q", queue, "-n", str(a.threads), "-M", str(mem),
                    "-R", f"span[hosts=1] select[mem>={mem}] rusage[mem={mem}]",
                    "-o", os.path.join(SCAN, "logs", f"{r['tolid']}.%J.out"),
                    "-e", os.path.join(SCAN, "logs", f"{r['tolid']}.%J.err"),
                    "-J", f"{group}.{r['tolid']}"],
                   f"bash {os.path.join(HERE, 'screen_one.sh')} {d} {a.threads}")
        state[(r["tolid"], r["version"])] = dict(tolid=r["tolid"], version=r["version"], status="submitted",
                                                 date=today, job=job)
    write(state_path, list(state.values()), STATE_FIELDS)

    report = f"{sys.executable} {os.path.join(HERE, 'report.py')} {today}"
    dep = ["-w", f"ended({group}.*)"] if batch else []
    bsub(["-q", a.queue, "-n", "1", "-M", "4000", "-R", "select[mem>=4000] rusage[mem=4000]", *dep,
          "-o", os.path.join(SCAN, "logs", f"report.{today}.%J.out"), "-J", f"scanreport.{today}"], report)

    if not a.once:
        nxt = datetime.datetime.now() + datetime.timedelta(days=7)
        cmd = (f"{sys.executable} {os.path.abspath(__file__)} --max {a.next_max} --next-max {a.next_max}"
               f" --queue {a.queue} --big-queue {a.big_queue} --queue-mem {a.queue_mem}")
        bsub(["-q", a.queue, "-n", "1", "-M", "4000", "-R", "select[mem>=4000] rusage[mem=4000]",
              "-b", nxt.strftime("%Y:%m:%d:%H:%M"),
              "-o", os.path.join(SCAN, "logs", "weekly.%J.out"), "-J", "scanweekly"], cmd)
        print("next run", nxt.strftime("%Y-%m-%d %H:%M"))


if __name__ == "__main__":
    main()
