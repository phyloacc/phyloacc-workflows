#!/usr/bin/env python3
"""
compare_to_atac.py

How often do our CEs/CNEEs (both real runs) and UCSC's raw mm39 35-way CEs
overlap an independent ATAC-seq peak set (skin, late anagen,
GSM4156597_skin.late.anagen.peaks.mm39.bed, translated to GenBank accessions
the same way as the UCSC conservation track), per chromosome. Writes one flat
TSV of overlap statistics. bedtools called by absolute path - see
compare_to_ucsc.py for why.

Usage: python3 compare_to_atac.py <atac_genbank_bed> <ucsc_genbank_bed> <output_tsv> [chrom ...]
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


def chrom_subset(bed_path, chrom, out_path):
    with open(out_path, "w") as out:
        with open(bed_path) as f:
            for line in f:
                if line.startswith(chrom + "\t"):
                    out.write(line)


def normalize_chrom(bed_path, chrom, out_path):
    # Same CE chrom-truncation workaround as compare_to_ucsc.py.
    with open(out_path, "w") as out_f, open(bed_path) as in_f:
        for line in in_f:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            parts[0] = chrom
            out_f.write("\t".join(parts) + "\n")


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


def run_intersect_count(a_path, b_path):
    # -u: report each element in a_path at most once, if it overlaps b_path at all -
    # this is "how many elements overlap", not "how many overlapping pairs".
    result = subprocess.run(
        [BEDTOOLS, "intersect", "-a", a_path, "-b", b_path, "-u"],
        capture_output=True, text=True, check=True,
    )
    return len([line for line in result.stdout.strip().split("\n") if line])


def run_jaccard(a_path, b_path):
    result = subprocess.run(
        [BEDTOOLS, "jaccard", "-a", a_path, "-b", b_path],
        capture_output=True, text=True, check=True,
    )
    lines = result.stdout.strip().split("\n")
    header = lines[0].split("\t")
    vals = lines[1].split("\t")
    return dict(zip(header, vals))


def sort_bed(path, tmpdir, label):
    out = os.path.join(tmpdir, f"{label}.sorted.bed")
    subprocess.run([BEDTOOLS, "sort", "-i", path], stdout=open(out, "w"), check=True)
    return out


def main():
    atac_bed = sys.argv[1]
    ucsc_bed = sys.argv[2]
    out_tsv = sys.argv[3]
    chrom_filter = set(sys.argv[4:]) if len(sys.argv) > 4 else None

    rows = []
    with tempfile.TemporaryDirectory() as tmpdir:
        atac_by_chrom = {}
        ucsc_by_chrom = {}
        for group, chroms in CHROM_GROUPS.items():
            for chrom in chroms:
                if chrom_filter and chrom not in chrom_filter:
                    continue
                atac_chrom_bed = os.path.join(tmpdir, f"atac.{chrom}.bed")
                chrom_subset(atac_bed, chrom, atac_chrom_bed)
                if bp_total(atac_chrom_bed)[0] == 0:
                    atac_by_chrom[chrom] = None
                    print(f"NOTE: no ATAC peaks on {chrom} (expected for Y/M)", file=sys.stderr)
                else:
                    atac_by_chrom[chrom] = sort_bed(atac_chrom_bed, tmpdir, f"atac.{chrom}")

                ucsc_chrom_bed = os.path.join(tmpdir, f"ucsc.{chrom}.bed")
                chrom_subset(ucsc_bed, chrom, ucsc_chrom_bed)
                ucsc_by_chrom[chrom] = sort_bed(ucsc_chrom_bed, tmpdir, f"ucsc.{chrom}")

        def add_row(source, group, chrom, elem_bed):
            atac_sorted = atac_by_chrom[chrom]
            n_elem, bp_elem = bp_total(elem_bed)
            if atac_sorted is None:
                n_atac, bp_atac, bp_intersect, n_elem_overlap, jaccard = 0, 0, 0, 0, 0.0
            else:
                n_atac, bp_atac = bp_total(atac_sorted)
                bp_intersect = run_intersect_bp(elem_bed, atac_sorted)
                n_elem_overlap = run_intersect_count(elem_bed, atac_sorted)
                jaccard = float(run_jaccard(elem_bed, atac_sorted)["jaccard"])
            rows.append({
                "source": source,
                "group": group,
                "chrom": chrom,
                "n_elements": n_elem,
                "bp_elements": bp_elem,
                "n_atac_peaks": n_atac,
                "bp_atac_peaks": bp_atac,
                "n_elements_overlapping": n_elem_overlap,
                "pct_elements_overlapping": round(100 * n_elem_overlap / n_elem, 2) if n_elem else None,
                "bp_intersect": bp_intersect,
                "jaccard": jaccard,
                "pct_elements_bp_in_atac": round(100 * bp_intersect / bp_elem, 2) if bp_elem else None,
                "pct_atac_bp_in_elements": round(100 * bp_intersect / bp_atac, 2) if bp_atac else None,
            })
            print(f"{source} {chrom}: n_elements={n_elem} n_atac={n_atac} n_overlapping={n_elem_overlap} jaccard={jaccard} pct_elements_in_atac={rows[-1]['pct_elements_bp_in_atac']}", file=sys.stderr)

        for run_name, run_dir in RUNS.items():
            for track in ("ce", "cnee"):
                for group, chroms in CHROM_GROUPS.items():
                    for chrom in chroms:
                        if chrom_filter and chrom not in chrom_filter:
                            continue
                        our_bed = track_bed_path(run_dir, track, group, chrom)
                        if not os.path.isfile(our_bed):
                            print(f"SKIP missing: {our_bed}", file=sys.stderr)
                            continue
                        normalized = os.path.join(tmpdir, f"{run_name}.{track}.{chrom}.normalized.bed")
                        normalize_chrom(our_bed, chrom, normalized)
                        elem_sorted = sort_bed(normalized, tmpdir, f"{run_name}.{track}.{chrom}")
                        add_row(f"{run_name}_{track}", group, chrom, elem_sorted)

        for group, chroms in CHROM_GROUPS.items():
            for chrom in chroms:
                if chrom_filter and chrom not in chrom_filter:
                    continue
                add_row("ucsc_ce", group, chrom, ucsc_by_chrom[chrom])

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "source", "group", "chrom", "n_elements", "bp_elements",
            "n_atac_peaks", "bp_atac_peaks", "n_elements_overlapping",
            "pct_elements_overlapping", "bp_intersect", "jaccard",
            "pct_elements_bp_in_atac", "pct_atac_bp_in_elements",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
