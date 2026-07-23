#!/usr/bin/env python3
"""
compare_to_ucsc.py

Compares our pipeline's CE (04-phastcons/regions) and CNEE
(05-cnees/phastcons/bed) predictions, for both the num_seqs and ns real runs,
against the translated UCSC mm39 35-way conserved-elements track, per
chromosome. Writes one flat TSV of overlap statistics.

bedtools is called by absolute path (bedtools-env conda env) rather than via
`conda activate`, since this repo's sessions have twice hit a real bug where
`conda activate` doesn't reliably win on $PATH.

Usage: python3 compare_to_ucsc.py <ucsc_genbank_bed> <output_tsv> [chrom ...]
       (chrom args restrict to a subset, for quick sanity-check runs; omit
       for the full sweep)
"""

import csv
import os
import subprocess
import sys
import tempfile

BEDTOOLS = "/n/home07/gthomas/miniconda3/envs/bedtools-env/bin/bedtools"

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

RUNS = {
    "num_seqs": os.path.join(REPO_ROOT, "data", "hamsters", "workflow-tests", "test-full-split-aln"),
    "ns": os.path.join(REPO_ROOT, "data", "hamsters", "workflow-tests", "test-full-gc"),
}

CHROM_GROUPS = {
    "autosomes": [f"CM0009{n:02d}.3" for n in range(94, 100)] + [f"CM0010{n:02d}.3" for n in range(0, 13)],
    "xy": ["CM001013.3", "CM001014.3"],
}


def track_bed_path(run_dir, track, group, chrom):
    if track == "ce":
        return os.path.join(run_dir, "04-phastcons", "regions", group, f"{chrom}.bed")
    return os.path.join(run_dir, "05-cnees", "phastcons", "bed", group, f"{chrom}.cnees.bed4")


def cds_bed_path(run_dir, group, chrom):
    # extract_cds_bed_chr's output path was flattened mid-session (3-level nested
    # -> flat) - test-full-split-aln postdates that fix, test-full-gc predates it,
    # so check both rather than assuming one shape.
    flat = os.path.join(run_dir, "05-cnees", "phastcons", "bed", group, f"{chrom}.cds.bed")
    if os.path.isfile(flat):
        return flat
    nested = os.path.join(run_dir, "05-cnees", "phastcons", "bed", group, chrom, f"{chrom}.cds.bed")
    if os.path.isfile(nested):
        return nested
    return None


def chrom_subset(bed_path, chrom, out_path):
    with open(out_path, "w") as out:
        with open(bed_path) as f:
            for line in f:
                if line.startswith(chrom + "\t"):
                    out.write(line)


def bp_total(bed_path):
    n = 0
    bp = 0
    with open(bed_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            n += 1
            bp += int(p[2]) - int(p[1])
    return n, bp


def run_jaccard(a_path, b_path):
    result = subprocess.run(
        [BEDTOOLS, "jaccard", "-a", a_path, "-b", b_path],
        capture_output=True, text=True, check=True,
    )
    lines = result.stdout.strip().split("\n")
    header = lines[0].split("\t")
    vals = lines[1].split("\t")
    return dict(zip(header, vals))


def run_intersect_bp(a_path, b_path):
    result = subprocess.run(
        [BEDTOOLS, "intersect", "-a", a_path, "-b", b_path],
        capture_output=True, text=True, check=True,
    )
    bp = 0
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        p = line.split("\t")
        bp += int(p[2]) - int(p[1])
    return bp


def main():
    ucsc_bed = sys.argv[1]
    out_tsv = sys.argv[2]
    chrom_filter = set(sys.argv[3:]) if len(sys.argv) > 3 else None

    rows = []
    with tempfile.TemporaryDirectory() as tmpdir:
        ucsc_by_chrom = {}
        for group, chroms in CHROM_GROUPS.items():
            for chrom in chroms:
                if chrom_filter and chrom not in chrom_filter:
                    continue
                ucsc_chrom_bed = os.path.join(tmpdir, f"ucsc.{chrom}.bed")
                chrom_subset(ucsc_bed, chrom, ucsc_chrom_bed)
                ucsc_by_chrom[chrom] = ucsc_chrom_bed

        for run_name, run_dir in RUNS.items():
            for track in ("ce", "cnee"):
                for group, chroms in CHROM_GROUPS.items():
                    for chrom in chroms:
                        if chrom_filter and chrom not in chrom_filter:
                            continue
                        our_bed = track_bed_path(run_dir, track, group, chrom)
                        ucsc_chrom_bed = ucsc_by_chrom[chrom]
                        if not os.path.isfile(our_bed):
                            print(f"SKIP missing: {our_bed}", file=sys.stderr)
                            continue

                        # CE beds (04-phastcons/regions) have phastCons's own truncated
                        # chrom name (e.g. "CM000994", no version suffix) rather than the
                        # full accession every other output uses - cnees_from_conserved_chr
                        # already has to normalize this internally when reading the same
                        # file, so do the same here rather than assume it matches `chrom`.
                        our_normalized = os.path.join(tmpdir, f"{run_name}.{track}.{chrom}.normalized.bed")
                        with open(our_normalized, "w") as out_f, open(our_bed) as in_f:
                            for line in in_f:
                                parts = line.rstrip("\n").split("\t")
                                if not parts or not parts[0]:
                                    continue
                                parts[0] = chrom
                                out_f.write("\t".join(parts) + "\n")

                        our_sorted = os.path.join(tmpdir, f"{run_name}.{track}.{chrom}.sorted.bed")
                        subprocess.run(
                            [BEDTOOLS, "sort", "-i", our_normalized],
                            stdout=open(our_sorted, "w"), check=True,
                        )
                        ucsc_sorted = os.path.join(tmpdir, f"ucsc.{chrom}.sorted.bed")
                        subprocess.run(
                            [BEDTOOLS, "sort", "-i", ucsc_chrom_bed],
                            stdout=open(ucsc_sorted, "w"), check=True,
                        )

                        n_ours, bp_ours = bp_total(our_sorted)
                        n_ucsc, bp_ucsc = bp_total(ucsc_sorted)
                        jac = run_jaccard(our_sorted, ucsc_sorted)
                        bp_intersect = run_intersect_bp(our_sorted, ucsc_sorted)

                        rows.append({
                            "our_run": run_name,
                            "our_track": track,
                            "group": group,
                            "chrom": chrom,
                            "n_ours": n_ours,
                            "bp_ours": bp_ours,
                            "n_ucsc": n_ucsc,
                            "bp_ucsc": bp_ucsc,
                            "bp_intersect": bp_intersect,
                            "jaccard": jac["jaccard"],
                            "pct_ours_in_ucsc": round(100 * bp_intersect / bp_ours, 2) if bp_ours else None,
                            "pct_ucsc_in_ours": round(100 * bp_intersect / bp_ucsc, 2) if bp_ucsc else None,
                        })
                        print(f"{run_name} {track} {chrom}: n_ours={n_ours} n_ucsc={n_ucsc} jaccard={jac['jaccard']}", file=sys.stderr)

                        # Fairer CNEE comparison: our CNEEs already have CDS subtracted,
                        # but UCSC's track doesn't - subtract our own CDS bed (same GFF,
                        # already computed by extract_cds_bed_chr) from the UCSC bed
                        # before comparing, so both sides exclude coding regions.
                        if track == "cnee":
                            cds_bed = cds_bed_path(run_dir, group, chrom)
                            if cds_bed is None:
                                print(f"SKIP cnee_vs_nocds_ucsc (no cds bed): {run_name} {chrom}", file=sys.stderr)
                                continue
                            cds_sorted = os.path.join(tmpdir, f"{run_name}.{chrom}.cds.sorted.bed")
                            subprocess.run(
                                [BEDTOOLS, "sort", "-i", cds_bed],
                                stdout=open(cds_sorted, "w"), check=True,
                            )
                            ucsc_nocds = os.path.join(tmpdir, f"ucsc.{chrom}.nocds.bed")
                            subprocess.run(
                                [BEDTOOLS, "subtract", "-a", ucsc_sorted, "-b", cds_sorted],
                                stdout=open(ucsc_nocds, "w"), check=True,
                            )
                            n_ucsc_nc, bp_ucsc_nc = bp_total(ucsc_nocds)
                            jac_nc = run_jaccard(our_sorted, ucsc_nocds)
                            bp_intersect_nc = run_intersect_bp(our_sorted, ucsc_nocds)
                            rows.append({
                                "our_run": run_name,
                                "our_track": "cnee_vs_nocds_ucsc",
                                "group": group,
                                "chrom": chrom,
                                "n_ours": n_ours,
                                "bp_ours": bp_ours,
                                "n_ucsc": n_ucsc_nc,
                                "bp_ucsc": bp_ucsc_nc,
                                "bp_intersect": bp_intersect_nc,
                                "jaccard": jac_nc["jaccard"],
                                "pct_ours_in_ucsc": round(100 * bp_intersect_nc / bp_ours, 2) if bp_ours else None,
                                "pct_ucsc_in_ours": round(100 * bp_intersect_nc / bp_ucsc_nc, 2) if bp_ucsc_nc else None,
                            })
                            print(f"{run_name} cnee_vs_nocds_ucsc {chrom}: n_ours={n_ours} n_ucsc_nocds={n_ucsc_nc} jaccard={jac_nc['jaccard']}", file=sys.stderr)

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "our_run", "our_track", "group", "chrom", "n_ours", "bp_ours",
            "n_ucsc", "bp_ucsc", "bp_intersect", "jaccard",
            "pct_ours_in_ucsc", "pct_ucsc_in_ours",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
