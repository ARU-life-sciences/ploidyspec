import glob
import itertools
import os
import shutil
import subprocess
import sys
import time


def log(msg):
    ts = time.strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", file=sys.stderr, flush=True)


def repo_root():
    # this file lives at <repo_root>/ploidyspec/common.py
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find_tool(explicit, dir_glob, binary_name):
    """Locate a binary: explicit path/dir > sibling dir matching dir_glob > PATH."""
    if explicit:
        cand = explicit if os.path.basename(explicit) == binary_name else os.path.join(explicit, binary_name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return os.path.abspath(cand)
        raise SystemExit(f"{binary_name} not found at {explicit!r}")
    for d in sorted(glob.glob(os.path.join(repo_root(), dir_glob))):
        cand = os.path.join(d, binary_name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return os.path.abspath(cand)
    which = shutil.which(binary_name)
    if which:
        return which
    raise SystemExit(
        f"could not locate {binary_name} (looked in {dir_glob}/ and $PATH); "
        f"pass an explicit path via the CLI flags"
    )


def run(cmd, stdout=None, check=True):
    """Run a command with no shell involved. Returns CompletedProcess (stdout captured as text unless redirected)."""
    kwargs = dict(stderr=subprocess.PIPE, text=True)
    if stdout is not None:
        kwargs["stdout"] = stdout
    else:
        kwargs["stdout"] = subprocess.PIPE
    result = subprocess.run(cmd, **kwargs)
    if check and result.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n{result.stderr}")
    return result


def sum_hist_distinct(histex_bin, hist_prefix):
    """Sum column 2 of `Histex -A <prefix>` = total distinct k-mers represented by a .hist file."""
    result = run([histex_bin, "-A", hist_prefix])
    total = 0
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            total += int(parts[1])
    return total


def group_pair_batches(n, letter_cap=8):
    """
    Partition all C(n,2) unordered pairs among n items into Logex-sized batches.

    Each batch is (global_indices, local_pairs): global_indices is the ordered list of
    item indices to pass as Logex sources (mapped positionally to letters A,B,C...),
    local_pairs is a list of (local_i, local_j) index pairs *within that source list*
    to intersect. If n <= letter_cap, a single batch covering every pair is returned
    (one Logex process computes the whole group at once). Otherwise the pairs are
    chunked so each batch uses at most letter_cap source tables.
    """
    if n <= letter_cap:
        return [(list(range(n)), list(itertools.combinations(range(n), 2)))]
    batches = []
    chunk_size = letter_cap - 1
    for i in range(n):
        targets = list(range(i + 1, n))
        for k in range(0, len(targets), chunk_size):
            chunk = targets[k : k + chunk_size]
            if not chunk:
                continue
            idxs = [i] + chunk
            local_pairs = [(0, pos + 1) for pos in range(len(chunk))]
            batches.append((idxs, local_pairs))
    return batches


def run_logex_batch(logex_bin, histex_bin, source_prefixes, local_pairs, tmp_dir):
    """Run one Logex call computing several pairwise intersections at once, return {(a,b): shared_count}."""
    os.makedirs(tmp_dir, exist_ok=True)
    letters = [chr(ord("A") + k) for k in range(len(source_prefixes))]
    exprs = []
    names = []
    for a, b in local_pairs:
        name = os.path.join(tmp_dir, f"{a}_{b}")
        exprs.append(f"{name}={letters[a]} &. {letters[b]}")
        names.append(name)
    cmd = [logex_bin, "-H"] + exprs + source_prefixes
    run(cmd)
    shared = {}
    for (a, b), name in zip(local_pairs, names):
        shared[(a, b)] = sum_hist_distinct(histex_bin, name)
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return shared
