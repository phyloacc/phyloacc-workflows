#!/usr/bin/env python3
# Flattens Snakemake benchmark TSVs from a real (non-dry-run) workflow-tests
# run into one long-format TSV for downstream analysis in R. Path shape is
# logs/benchmarks/<rule>/<group>/<chrom>.txt (per-chromosome rules),
# logs/benchmarks/<rule>/<group>.txt (per-group rules), or
# logs/benchmarks/<rule>/run.txt (single whole-genome rules).

import os
import sys
import csv

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <benchmarks_dir> <output_tsv>", file=sys.stderr)
        sys.exit(1)

    bench_dir, out_path = sys.argv[1], sys.argv[2]
    rows = []

    for rule in sorted(os.listdir(bench_dir)):
        rule_dir = os.path.join(bench_dir, rule)
        if not os.path.isdir(rule_dir):
            continue
        for root, _, files in os.walk(rule_dir):
            for fname in sorted(files):
                if not fname.endswith(".txt"):
                    continue
                fpath = os.path.join(root, fname)
                rel = os.path.relpath(fpath, rule_dir)
                parts = rel.replace(".txt", "").split(os.sep)

                if len(parts) == 1:
                    if parts[0] == "run":
                        group, chrom = "", ""
                    else:
                        group, chrom = parts[0], ""
                elif len(parts) == 2:
                    group, chrom = parts[0], parts[1]
                else:
                    group, chrom = parts[0], os.sep.join(parts[1:])

                with open(fpath) as fh:
                    lines = fh.read().strip().split("\n")
                if len(lines) < 2:
                    continue
                header = lines[0].split("\t")
                for line in lines[1:]:
                    vals = line.split("\t")
                    d = dict(zip(header, vals))
                    rows.append({
                        "rule": rule,
                        "group": group,
                        "chrom": chrom,
                        "wall_s": d.get("s", ""),
                        "max_rss_mb": d.get("max_rss", ""),
                        "cpu_time_s": d.get("cpu_time", ""),
                        "mean_load_pct": d.get("mean_load", ""),
                    })

    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=[
            "rule", "group", "chrom", "wall_s", "max_rss_mb", "cpu_time_s", "mean_load_pct"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_path}")

if __name__ == "__main__":
    main()
