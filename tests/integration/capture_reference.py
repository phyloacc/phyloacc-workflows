#!/usr/bin/env python3
"""
capture_reference.py

Runs the tier-3 fixture pipeline for real and writes its filter-summary
numbers into tests/integration/data/reference_summary.json - the "silver
standard" tests/integration/test_pipeline_e2e.py softly (warning, not
failing) compares future runs against. Re-run this deliberately after an
intentional change that's expected to shift these numbers (e.g. a real
algorithm change to CE/CNEE filtering) to re-baseline; not run automatically
by the test suite itself.

Usage: python3 tests/integration/capture_reference.py
"""

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

import yaml

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SNAKEFILE = os.path.join(REPO_ROOT, "Snakefile")
FIXTURE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/data"
FIXTURE_CONFIG = os.path.join(FIXTURE_DIR, "config.yaml")
REFERENCE_SUMMARY = os.path.join(FIXTURE_DIR, "reference_summary.json")

GROUP = "group1"
CHROM = "CM000994.3"

REQUIRED_TOOLS = ["mafutils", "phyloFit", "phastCons"]

# See test_pipeline_e2e.py's SNAKEMAKE_EXE comment: `python -m snakemake` makes
# lib/common.py's pipelineSetup() misclassify this as a Snakemake "worker"
# re-invocation rather than the real top-level run.
SNAKEMAKE_EXE = os.path.join(os.path.dirname(sys.executable), "snakemake")


def main():
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        print(f"Missing required tool(s) on PATH: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    with tempfile.TemporaryDirectory() as tmp_dir:
        with open(FIXTURE_CONFIG) as f:
            config = yaml.safe_load(f)
        config["output_dir"] = os.path.join(tmp_dir, "out")
        config["tmp_dir"] = os.path.join(tmp_dir, "tmp")

        config_path = os.path.join(tmp_dir, "config.yaml")
        with open(config_path, "w") as f:
            yaml.dump(config, f)

        result = subprocess.run(
            [SNAKEMAKE_EXE, "-j", "1", "-s", SNAKEFILE, "--configfile", config_path],
            cwd=REPO_ROOT,
        )
        if result.returncode != 0:
            print("Pipeline run failed - not writing a reference from a failed run.", file=sys.stderr)
            sys.exit(1)

        summary_path = os.path.join(
            config["output_dir"], "05-cnees", "phastcons", "summary", GROUP, f"{CHROM}.cnees-filter-summary.tsv"
        )
        with open(summary_path) as f:
            values = {row["metric"]: int(row["value"]) for row in csv.DictReader(f, delimiter="\t")}

    with open(REFERENCE_SUMMARY, "w") as f:
        json.dump(values, f, indent=2)
        f.write("\n")

    print(f"Wrote {REFERENCE_SUMMARY}:")
    for k, v in values.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
