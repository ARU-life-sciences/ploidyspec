#!/usr/bin/env python3
"""Generate manifests/*.tsv from meta/meta.tsv (real farm paths) for the ploidyspec batch run.

Column layout of meta/meta.tsv (no header row): note, species, dir, path1[, path2[, path3, path4]]
Paths are listed in haplotype order (hap1, hap2, [hap3, hap4] or hap1, alternate_haplotype).
"""
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
META = os.path.join(REPO, "meta", "meta.tsv")
MANIFEST_DIR = os.path.join(REPO, "manifests")

# species with existing hand-curated manifests -- leave hap_label scheme untouched,
# only refresh the fasta paths to the real farm locations.
EXISTING = {"daBudDavi1", "daGalBore1", "daGleHede1"}

# daSenVisc1's two "haplotype" paths are both cataloged as hap2 (GCA_965637505.1 and
# .2 -- same accession root, different assembly versions) -- not a real hap1/hap2 pair.
# Needs a human to confirm the correct hap1 path before this can run.
SKIP = {"daSenVisc1"}

# durum wheat (AABB tetraploid): every "haplotype" file contains BOTH subgenomes for
# every chromosome (SUPER_1A/SUPER_1B/...), so it's split into two independent runs,
# one per subgenome, each comparing hap1 vs hap2 within that subgenome only.
SPLIT_SUBGENOME = {"lpTriTurg1": ["A", "B"]}

# Default chrom-regex is [r"chromosome:?\s*(\d+)\b", r"SUPER[_-](\d+)\b(?!_)"]. Several
# species suffix the chromosome number with a haplotype tag the *other* way round
# (e.g. "SUPER_1_HAP2" instead of "HAP2_SUPER_1") -- the bare SUPER pattern's
# (?!_) lookahead deliberately excludes that (it also excludes "_unloc_" junk, which
# we still want excluded), so add a pattern that specifically matches the
# "SUPER_<n>_HAP<m>" case ahead of the bare one. Applied uniformly since it's a
# strict superset of the default and doesn't false-positive on any sampled species.
GLOBAL_CHROM_REGEX = [
    r"chromosome:?\s*(\d+)\b",
    r"SUPER[_-](\d+)_HAP\d+\b",
    r"SUPER[_-](\d+)\b(?!_)",
]

# Default hap-regex only catches "HAP<n>_SUPER" (prefix order). Header sampling turned
# up several species whose "hap2" (or hap3/hap4) file is *not* actually single-haplotype
# -- it's a grab-bag of several haplotype copies distinguished only by a per-sequence
# "..._HAP<n>" or "HAP<n>_..." tag (ddEmpNigr1 and ddHypPerf1's hap2 files each turned
# out to contain both HAP2- and HAP3-tagged contigs; trusting the file's position label
# caused a duplicate-unit_id crash for both). Wherever a file's header convention embeds
# a haplotype tag at all, use AUTO (which self-sorts by the tag, whether the file turns
# out single- or mixed-haplotype) instead of trusting file position -- position is only
# used where the header carries no haplotype tag whatsoever (so there's nothing to
# auto-detect from, e.g. plain "SUPER_N" or bare "SUPER_N" with no hap suffix).
UNIFIED_HAP_REGEX = r"(?:SUPER_\d+_)?HAP(\d+)(?:_SUPER)?"

# species -> 1-based positions (within that species' path list) whose header embeds a
# haplotype tag and should therefore use AUTO instead of an explicit HAP<position> label.
# (daBudDavi1/daGalBore1 already use AUTO for position 2 in their existing curated
# manifests, handled separately via EXISTING.)
AUTO_POSITIONS = {
    "ddEmpNigr1": {2},
    "ddHypPerf1": {2},
    "ddHypMacu1": {1, 2, 3, 4},
    "drLytSali1": {1, 2},
    "lpElePalu1": {1, 2},
}


def resolve(path):
    if os.path.exists(path):
        return path
    if path.endswith(".fasta.gz"):
        alt = path[: -len(".fasta.gz")] + ".fa.gz"
    elif path.endswith(".fa.gz"):
        alt = path[: -len(".fa.gz")] + ".fasta.gz"
    else:
        alt = path
    if os.path.exists(alt):
        return alt
    raise SystemExit(f"cannot resolve path: {path}")


def restage_if_curated(species, path):
    """Curated assemblies are plain-gzip (not bgzip) in read-only dirs -- samtools
    faidx can't index them in place. build_manifests assumes scripts/stage_curated.sh
    has already re-bgzipped them into data/<species>/<basename>."""
    if "/assembly/curated/" in path:
        staged = os.path.join(REPO, "data", species, os.path.basename(path))
        return staged
    return path


def main():
    os.makedirs(MANIFEST_DIR, exist_ok=True)
    jobs = []  # (outdir_species_id, manifest_path, extra_chrom_regex_args or None)

    with open(META) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            note = parts[0]
            species = parts[1]
            raw_paths = [p for p in parts[3:] if p.strip()]
            paths = [restage_if_curated(species, resolve(p)) for p in raw_paths]

            if species in SKIP:
                print(f"SKIP {species}: {note or 'flagged for manual review'}")
                continue

            if species in EXISTING:
                manifest_path = os.path.join(MANIFEST_DIR, f"{species}.tsv")
                with open(manifest_path) as mf:
                    existing_labels = [
                        l.rstrip("\n").split("\t")[1]
                        for l in mf
                        if l.strip() and not l.startswith("#")
                    ]
                with open(manifest_path, "w") as mf:
                    mf.write("# fasta_path\thap_label\n")
                    for p, label in zip(paths, existing_labels):
                        mf.write(f"{p}\t{label}\n")
                jobs.append((species, manifest_path, GLOBAL_CHROM_REGEX, None))
                print(f"UPDATED {species}: {len(paths)} paths (labels kept: {existing_labels})")
                continue

            if species in SPLIT_SUBGENOME:
                for genome in SPLIT_SUBGENOME[species]:
                    sub_id = f"{species}_{genome}"
                    manifest_path = os.path.join(MANIFEST_DIR, f"{sub_id}.tsv")
                    with open(manifest_path, "w") as mf:
                        mf.write("# fasta_path\thap_label\n")
                        for i, p in enumerate(paths, start=1):
                            mf.write(f"{p}\tHAP{i}\n")
                    extra_regex = [
                        rf"chromosome:?\s*(\d+){genome}\b",
                        rf"SUPER[_-](\d+){genome}\b",
                    ]
                    jobs.append((sub_id, manifest_path, extra_regex, None))
                    print(f"NEW {sub_id}: {len(paths)} paths, subgenome {genome}")
                continue

            auto_positions = AUTO_POSITIONS.get(species, set())
            manifest_path = os.path.join(MANIFEST_DIR, f"{species}.tsv")
            with open(manifest_path, "w") as mf:
                mf.write("# fasta_path\thap_label\n")
                for i, p in enumerate(paths, start=1):
                    label = "AUTO" if i in auto_positions else f"HAP{i}"
                    mf.write(f"{p}\t{label}\n")
            hap_regex = UNIFIED_HAP_REGEX if auto_positions else None
            jobs.append((species, manifest_path, GLOBAL_CHROM_REGEX, hap_regex))
            print(f"NEW {species}: {len(paths)} paths" + (f" (AUTO positions: {sorted(auto_positions)})" if auto_positions else ""))

    with open(os.path.join(REPO, "scripts", "jobs.tsv"), "w") as jf:
        for outdir_id, manifest_path, extra_regex, hap_regex in jobs:
            regex_field = "|".join(extra_regex) if extra_regex else ""
            jf.write(f"{outdir_id}\t{manifest_path}\t{regex_field}\t{hap_regex or ''}\n")
    print(f"\n{len(jobs)} jobs written to scripts/jobs.tsv")


if __name__ == "__main__":
    main()
