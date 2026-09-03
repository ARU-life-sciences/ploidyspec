#!/usr/bin/env python3
"""
One-off setup for the 2026-09 Darwin assembly batch: 14 species with new
hap1/hap2(+) Darwin Tree of Life assemblies, mixing `release`-stage (direct
use, already valid bgzip) and `curated`-stage (needs zcat|bgzip re-staging,
same reason as scripts/stage_curated.sh) sources per haplotype. Writes
manifests/*.tsv directly rather than going through build_manifests.py's
meta.tsv-driven flow, since several species mix sources per-haplotype in a
way that flow doesn't model.
"""
import gzip
import os
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")
MANIFESTS = os.path.join(REPO, "manifests")

GLOBAL_CHROM_REGEX = [
    r"chromosome:?\s*(\d+)\b",
    r"SUPER[_-](\d+)_HAP\d+\b",
    r"SUPER[_-](\d+)\b(?!_)",
]
UNIFIED_HAP_REGEX = r"(?:SUPER_\d+_)?HAP(\d+)(?:_SUPER)?"
AUTO_SPECIES = {"drLytSali1", "ddHesMatr1"}
# ddHypMacu1/lpElePalu1: their new Darwin `release` assemblies split cleanly
# to one haplotype per file with INSDC "chromosome: N" headers (no HAP tag
# to auto-detect at all), unlike their old `curated`-stage files -- switched
# from AUTO to positional after a real prepare failure (0 units parsed,
# AUTO couldn't find any HAP tag in "chromosome: N" headers).
# ddHesMatr1: the opposite discovery -- its curated hap2 file is a grab-bag
# mixing HAP2/HAP3/HAP4 tags together (same pattern as ddEmpNigr1/ddHypPerf1
# below), caught via a duplicate-unit_id crash under positional labeling.

# (source_type, path) per haplotype. "release" paths used directly (already
# valid bgzip, confirmed elsewhere in this project). "curated" paths staged
# into data/<species>/ via zcat|bgzip (matches stage_curated.sh; these are
# plain gzip, not bgzip, so samtools faidx can't index them in place).
SPECIES = {
    "dcCerAlpi1": [
        ("release", "/lustre/scratch122/tol/data/a/2/6/9/1/4/Cerastium_alpinum/assembly/release/GCA_985607225.1/insdc/GCA_985607225.1.fa.gz"),
        ("release", "/lustre/scratch122/tol/data/a/2/6/9/1/4/Cerastium_alpinum/assembly/release/GCA_985607205.1/insdc/GCA_985607205.1.fa.gz"),
    ],
    "ddHypMacu1": [
        ("release", "/lustre/scratch122/tol/data/f/3/f/3/4/3/Hypericum_maculatum/assembly/release/GCA_986264625.1/insdc/GCA_986264625.1.fa.gz"),
        ("release", "/lustre/scratch122/tol/data/f/3/f/3/4/3/Hypericum_maculatum/assembly/release/GCA_986264605.1/insdc/GCA_986264605.1.fa.gz"),
        ("release", "/lustre/scratch122/tol/data/f/3/f/3/4/3/Hypericum_maculatum/assembly/release/GCA_986264615.1/insdc/GCA_986264615.1.fa.gz"),
        ("release", "/lustre/scratch122/tol/data/f/3/f/3/4/3/Hypericum_maculatum/assembly/release/GCA_986264645.1/insdc/GCA_986264645.1.fa.gz"),
    ],
    "drSorDevo1": [
        ("release", "/lustre/scratch122/tol/data/1/1/c/2/0/9/Sorbus_devoniensis/assembly/release/GCA_986449145.1/insdc/GCA_986449145.1.fa.gz"),
        ("release", "/lustre/scratch122/tol/data/1/1/c/2/0/9/Sorbus_devoniensis/assembly/release/GCA_986449215.1/insdc/GCA_986449215.1.fa.gz"),
    ],
    "lpElePalu1": [
        ("release", "/lustre/scratch122/tol/data/e/e/e/d/3/9/Eleocharis_palustris/assembly/release/GCA_986264755.1/insdc/GCA_986264755.1.fa.gz"),
        ("release", "/lustre/scratch122/tol/data/e/e/e/d/3/9/Eleocharis_palustris/assembly/release/GCA_986264795.1/insdc/GCA_986264795.1.fa.gz"),
    ],
    "daSenVulg1": [
        ("curated", "/lustre/scratch122/tol/data/f/f/7/c/3/e/Senecio_vulgaris/assembly/curated/daSenVulg1.1/daSenVulg1.hap1.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/f/f/7/c/3/e/Senecio_vulgaris/assembly/curated/daSenVulg1.1/daSenVulg1.hap2.1.primary.fa.gz"),
    ],
    "drMyrSpic1": [
        ("curated", "/lustre/scratch122/tol/data/8/5/c/0/0/1/Myriophyllum_spicatum/assembly/curated/drMyrSpic1.1/drMyrSpic1.hap1.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/8/5/c/0/0/1/Myriophyllum_spicatum/assembly/curated/drMyrSpic1.1/drMyrSpic1.hap2.1.primary.fa.gz"),
    ],
    "drMyrVert1": [
        ("curated", "/lustre/scratch122/tol/data/b/3/d/c/b/3/Myriophyllum_verticillatum/assembly/curated/drMyrVert1.1/drMyrVert1.hap1.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/b/3/d/c/b/3/Myriophyllum_verticillatum/assembly/curated/drMyrVert1.1/drMyrVert1.hap2.1.primary.fa.gz"),
    ],
    "daPilAura1": [
        ("curated", "/lustre/scratch122/tol/data/6/3/3/2/5/4/Pilosella_aurantiaca/assembly/curated/daPilAura1.1/daPilAura1.hap1.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/6/3/3/2/5/4/Pilosella_aurantiaca/assembly/curated/daPilAura1.1/daPilAura1.hap2.1.primary.fa.gz"),
    ],
    "ddSalTria1": [
        ("curated", "/lustre/scratch122/tol/data/0/6/b/6/1/4/Salix_triandra/assembly/curated/ddSalTria1.1/ddSalTria1.hap1.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/0/6/b/6/1/4/Salix_triandra/assembly/curated/ddSalTria1.1/ddSalTria1.hap2.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/0/6/b/6/1/4/Salix_triandra/assembly/curated/ddSalTria1.1/ddSalTria1.hap3.1.primary.fa.gz"),
    ],
    "ddLepDrab1": [
        ("curated", "/lustre/scratch122/tol/data/4/5/a/e/4/1/Lepidium_draba/assembly/curated/ddLepDrab1.1/ddLepDrab1.hap1.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/4/5/a/e/4/1/Lepidium_draba/assembly/curated/ddLepDrab1.1/ddLepDrab1.hap2.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/4/5/a/e/4/1/Lepidium_draba/assembly/curated/ddLepDrab1.1/ddLepDrab1.hap3.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/4/5/a/e/4/1/Lepidium_draba/assembly/curated/ddLepDrab1.1/ddLepDrab1.hap4.1.primary.fa.gz"),
    ],
    "ddHesMatr1": [
        ("curated", "/lustre/scratch122/tol/data/7/9/4/f/a/9/Hesperis_matronalis/assembly/curated/ddHesMatr1.1/ddHesMatr1.hap1.1.primary.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/7/9/4/f/a/9/Hesperis_matronalis/assembly/curated/ddHesMatr1.1/ddHesMatr1.hap2.1.primary.fa.gz"),
    ],
    "dmRanRepe1": [
        ("release", "/lustre/scratch122/tol/data/3/9/6/a/7/e/Ranunculus_repens/assembly/release/GCA_986449105.1/insdc/GCA_986449105.1.fa.gz"),
        ("curated", "/lustre/scratch122/tol/data/3/9/6/a/7/e/Ranunculus_repens/assembly/curated/dmRanRepe1.1/dmRanRepe1.hap2.1.primary.fa.gz"),
    ],
    "daLatClan1": [
        ("release", "/lustre/scratch122/tol/data/8/6/2/6/9/0/Lathraea_clandestina/assembly/release/GCA_986341785.1/insdc/GCA_986341785.1.fa.gz"),
        ("curated_plain", "/lustre/scratch122/tol/data/8/6/2/6/9/0/Lathraea_clandestina/assembly/curated/daLatClan1.1/daLatClan1.hap2.1.primary.curated.fa"),
    ],
    "drLytSali1": [
        ("curated", "/lustre/scratch122/tol/data/7/3/1/d/7/f/Lythrum_salicaria/assembly/curated/drLytSali1.1/drLytSali1.hap1.1.primary.fa.gz"),
        ("release", "/lustre/scratch122/tol/data/7/3/1/d/7/f/Lythrum_salicaria/assembly/release/GCA_986280985.1/insdc/GCA_986280985.1.fa.gz"),
    ],
}


def stage_curated(species, src, plain=False):
    """zcat|bgzip (or bgzip directly if the source is already plain-text) a
    curated fasta into data/<species>/, matching stage_curated.sh's method."""
    os.makedirs(os.path.join(DATA, species), exist_ok=True)
    base = os.path.basename(src)
    dest = os.path.join(DATA, species, base if base.endswith(".gz") else base + ".gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  skip (exists): {dest}")
        return dest
    tmp = dest + ".tmp"
    print(f"  staging: {src} -> {dest}")
    if plain:
        with open(src, "rb") as fin, open(tmp, "wb") as fout:
            subprocess.run(["bgzip", "-@4", "-c"], stdin=fin, stdout=fout, check=True)
    else:
        with open(src, "rb") as fin, open(tmp, "wb") as fout:
            zcat = subprocess.Popen(["zcat", src], stdout=subprocess.PIPE)
            bgzip = subprocess.Popen(["bgzip", "-@4", "-c"], stdin=zcat.stdout, stdout=fout)
            zcat.stdout.close()
            bgzip.communicate()
            if bgzip.returncode != 0:
                raise SystemExit(f"bgzip failed for {src}")
    os.rename(tmp, dest)
    subprocess.run(["samtools", "faidx", dest], check=True)
    return dest


def main():
    os.makedirs(MANIFESTS, exist_ok=True)
    jobs = []
    for species, haps in SPECIES.items():
        print(f"=== {species} ===")
        resolved = []
        for kind, path in haps:
            if kind == "release":
                resolved.append(path)
            elif kind == "curated":
                resolved.append(stage_curated(species, path, plain=False))
            elif kind == "curated_plain":
                resolved.append(stage_curated(species, path, plain=True))
            else:
                raise SystemExit(f"unknown kind {kind}")

        label_mode = "AUTO" if species in AUTO_SPECIES else "positional"
        manifest_path = os.path.join(MANIFESTS, f"{species}.tsv")
        with open(manifest_path, "w") as f:
            f.write("# fasta_path\thap_label\n")
            for i, p in enumerate(resolved, start=1):
                label = "AUTO" if label_mode == "AUTO" else f"HAP{i}"
                f.write(f"{p}\t{label}\n")
        print(f"  wrote {manifest_path} ({len(resolved)} haplotypes, {label_mode})")

        hap_regex = UNIFIED_HAP_REGEX if label_mode == "AUTO" else ""
        jobs.append((species, manifest_path, "|".join(GLOBAL_CHROM_REGEX), hap_regex))

    # merge into jobs.tsv: replace existing rows for these species, keep the rest
    jobs_path = os.path.join(REPO, "scripts", "jobs.tsv")
    existing = []
    species_set = set(SPECIES)
    if os.path.exists(jobs_path):
        with open(jobs_path) as f:
            for line in f:
                if line.split("\t")[0] not in species_set:
                    existing.append(line.rstrip("\n"))
    with open(jobs_path, "w") as f:
        for line in existing:
            f.write(line + "\n")
        for species, manifest_path, regex_field, hap_regex in jobs:
            f.write(f"{species}\t{manifest_path}\t{regex_field}\t{hap_regex}\n")
    print(f"\nupdated {jobs_path}: {len(jobs)} species (re)written")


if __name__ == "__main__":
    main()
