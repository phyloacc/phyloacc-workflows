#############################################################################
# Tier-3 test: run the real pipeline (real mafutils/phyloFit/phastCons calls,
# not mocked) end-to-end against a tiny real-data fixture - a genuine 260kb
# window (CM000994.3:3,030,000-3,290,000) extracted from the real hamster
# MAF/GFF, plus the real (unmodified) tree file. Nothing here is fabricated;
# see tests/integration/README.md and tests/integration/data/config.yaml for
# how the fixture was built.
#
# Needs the real conda env (mafutils/phyloFit/phastCons actually installed) -
# skips cleanly if they're not on PATH, rather than failing confusingly.
#
# Assertions are structural/invariant-based, not exact-value and not a
# golden-file comparison against a fixed run - see the "silver standard"
# section below for the one place a soft, non-failing comparison is made.
#############################################################################

import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import warnings

import pytest
import yaml

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SNAKEFILE = os.path.join(REPO_ROOT, "Snakefile")
FIXTURE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/data"
FIXTURE_CONFIG = os.path.join(FIXTURE_DIR, "config.yaml")
REFERENCE_SUMMARY = os.path.join(FIXTURE_DIR, "reference_summary.json")

# Invoke the actual `snakemake` console-script entry point, not `python -m snakemake` -
# see tests/test_dag_validation.py's identical SNAKEMAKE_EXE comment for why: `-m`
# always sets sys.argv[0] to snakemake's own __main__.py, which lib/common.py's
# pipelineSetup() misreads as a Snakemake-internal "worker" re-invocation rather than
# the real top-level run.
SNAKEMAKE_EXE = os.path.join(os.path.dirname(sys.executable), "snakemake")

REQUIRED_TOOLS = ["mafutils", "phyloFit", "phastCons"]

GROUP = "group1"
CHROM = "CM000994.3"

# Relative change beyond this fraction triggers a soft warning, not a failure -
# see the plan file's "silver standard" section for why this isn't a hard assert.
DRIFT_WARNING_THRESHOLD = 0.5


def _missing_tools():
    return [t for t in REQUIRED_TOOLS if shutil.which(t) is None]


@pytest.fixture(scope="module")
def pipeline_run(tmp_path_factory):
    missing = _missing_tools()
    if missing:
        pytest.skip(f"Real tool(s) not on PATH, skipping tier-3 integration test: {', '.join(missing)}")

    tmp_path = tmp_path_factory.mktemp("tier3")
    with open(FIXTURE_CONFIG) as f:
        config = yaml.safe_load(f)
    config["output_dir"] = str(tmp_path / "out")
    config["tmp_dir"] = str(tmp_path / "tmp")

    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    result = subprocess.run(
        [SNAKEMAKE_EXE, "-j", "1", "-s", SNAKEFILE, "--configfile", str(config_path)],
        capture_output=True, text=True, cwd=REPO_ROOT,
    )
    return result, config["output_dir"]


def test_pipeline_completes_successfully(pipeline_run):
    result, _ = pipeline_run
    assert result.returncode == 0, result.stdout + result.stderr


def test_filter_summary_invariants(pipeline_run):
    _, output_dir = pipeline_run
    summary_path = os.path.join(
        output_dir, "05-cnees", "phastcons", "summary", GROUP, f"{CHROM}.cnees-filter-summary.tsv"
    )
    assert os.path.isfile(summary_path)
    with open(summary_path) as f:
        rows = {row["metric"]: int(row["value"]) for row in csv.DictReader(f, delimiter="\t")}

    assert rows["ces_raw"] >= 1
    assert rows["ces_merged"] <= rows["ces_raw"]
    # The CDS-drop invariant confirmed by hand against real data earlier this session:
    # dropping can only ever reduce the count, never exceed the merged CE count.
    assert rows["cnees_after_cds_drop"] <= rows["ces_merged"]
    assert rows["ces_dropped_cds_overlap"] == rows["ces_merged"] - rows["cnees_after_cds_drop"]
    # This fixture window was specifically chosen (see plan file) to contain both a real
    # CDS-overlapping CE (dropped) and real non-overlapping CEs (kept) - confirm both
    # actually happened, not just that the counts are internally consistent.
    assert rows["ces_dropped_cds_overlap"] >= 1
    assert rows["cnees_after_cds_drop"] >= 1


def test_final_cnee_bed4_lengths(pipeline_run):
    _, output_dir = pipeline_run
    bed4_path = os.path.join(output_dir, "05-cnees", "phastcons", "bed", GROUP, f"{CHROM}.cnees.bed4")
    assert os.path.isfile(bed4_path)
    rows = []
    with open(bed4_path) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, s, e, cid = line.rstrip("\n").split("\t")
            rows.append((chrom, int(s), int(e), cid))
    assert len(rows) >= 1
    for chrom, s, e, cid in rows:
        assert chrom == CHROM
        assert (e - s) > 50  # cnee_min_len_bp in the fixture config


def test_cnee_fasta_output_produced(pipeline_run):
    _, output_dir = pipeline_run
    fasta_dir = os.path.join(output_dir, "05-cnees", "phastcons", "fasta", GROUP, CHROM)
    manifest_path = os.path.join(fasta_dir, "manifest.txt")
    assert os.path.isfile(manifest_path)
    with open(manifest_path) as f:
        manifest_files = [line.strip() for line in f if line.strip()]
    assert len(manifest_files) >= 1
    for fname in manifest_files:
        fasta_path = os.path.join(fasta_dir, fname)
        assert os.path.isfile(fasta_path)
        assert os.path.getsize(fasta_path) > 0


def test_summary_report_produced(pipeline_run):
    _, output_dir = pipeline_run
    # Report is named "<config-stem>.summary_report.html" (config-file-stem prefix), so
    # match by suffix rather than a hard-coded basename.
    reports = glob.glob(os.path.join(output_dir, "*summary_report.html"))
    assert reports, "no summary report produced"
    assert os.path.getsize(reports[0]) > 0


def test_silver_standard_no_dramatic_drift(pipeline_run):
    # Soft regression check (see plan file) - compares against a committed reference of
    # this exact fixture's own past real numbers, and warns (does not fail) if something
    # drifted by more than DRIFT_WARNING_THRESHOLD. Never a hard assert - real tool output
    # legitimately isn't guaranteed stable across phastCons/mafutils versions.
    if not os.path.isfile(REFERENCE_SUMMARY):
        pytest.skip("No reference_summary.json yet - run tests/integration/capture_reference.py once to seed it.")

    _, output_dir = pipeline_run
    summary_path = os.path.join(
        output_dir, "05-cnees", "phastcons", "summary", GROUP, f"{CHROM}.cnees-filter-summary.tsv"
    )
    with open(summary_path) as f:
        current = {row["metric"]: int(row["value"]) for row in csv.DictReader(f, delimiter="\t")}
    with open(REFERENCE_SUMMARY) as f:
        reference = json.load(f)

    for key, ref_val in reference.items():
        cur_val = current.get(key)
        if cur_val is None or ref_val == 0:
            continue
        relative_change = abs(cur_val - ref_val) / ref_val
        if relative_change > DRIFT_WARNING_THRESHOLD:
            warnings.warn(
                f"tier-3 silver-standard drift: '{key}' changed from {ref_val} (reference) to "
                f"{cur_val} (this run), a {relative_change:.0%} relative change - "
                f"worth a human look (could be a real regression, or just a tool version "
                f"difference - see the plan file for why this warns rather than fails).",
                UserWarning,
            )


#############################################################################
# phyloP -> CNEE path end-to-end (workflow/cnees.smk with source=phylop).
#
# The fixture tree is shallow (hamster), so the phyloP power gate would normally
# stop this stage - we set phylop_power_override to run it anyway, and phyloP finds
# few/no conserved sites, so the phyloP CNEE set is legitimately (near-)empty. The
# point is to exercise the real machinery end-to-end: gate override, run_phylop,
# clustering, and the graceful empty short-circuit through cnees.smk - and to confirm
# it produces its source-namespaced outputs and exits 0 rather than erroring on empty.

@pytest.fixture(scope="module")
def pipeline_run_phylop(tmp_path_factory):
    missing = _missing_tools()
    if missing:
        pytest.skip(f"Real tool(s) not on PATH, skipping tier-3 phyloP test: {', '.join(missing)}")

    tmp_path = tmp_path_factory.mktemp("tier3_phylop")
    with open(FIXTURE_CONFIG) as f:
        config = yaml.safe_load(f)
    config.update({
        "run_phylofit": True,
        "run_phylop": True,
        "run_phastcons": False,
        "build_cnees": True,
        "cnee_output_format": "fasta",
        "phylop_power_override": True,  # shallow fixture: proceed despite the gate
    })
    config["output_dir"] = str(tmp_path / "out")
    config["tmp_dir"] = str(tmp_path / "tmp")
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    result = subprocess.run(
        [SNAKEMAKE_EXE, "-j", "1", "-s", SNAKEFILE, "--configfile", str(config_path)],
        capture_output=True, text=True, cwd=REPO_ROOT,
    )
    return result, config["output_dir"]


def test_phylop_pipeline_completes(pipeline_run_phylop):
    result, _ = pipeline_run_phylop
    assert result.returncode == 0, result.stdout + result.stderr


def test_phylop_power_check_report_written(pipeline_run_phylop):
    _, output_dir = pipeline_run_phylop
    report = os.path.join(output_dir, "03-phylop", "power-check", GROUP, f"{CHROM}.power.tsv")
    assert os.path.isfile(report)
    with open(report) as f:
        rows = {r["metric"]: r["value"] for r in csv.DictReader(f, delimiter="\t")}
    # Shallow tree -> the gate would fail; override let it proceed.
    assert rows["passes"] == "False"
    assert float(rows["tree_length"]) > 0


def test_phylop_cnee_outputs_namespaced_and_present(pipeline_run_phylop):
    # The phyloP-source CNEE set is produced under its own namespace (possibly empty).
    _, output_dir = pipeline_run_phylop
    base = os.path.join(output_dir, "05-cnees", "phylop")
    cnees_bed = os.path.join(base, "bed", GROUP, f"{CHROM}.cnees.bed")
    summary = os.path.join(base, "summary", GROUP, f"{CHROM}.cnees-filter-summary.tsv")
    manifest = os.path.join(base, "fasta", GROUP, CHROM, "manifest.txt")
    assert os.path.isfile(cnees_bed), "phyloP CNEE bed missing"
    assert os.path.isfile(summary), "phyloP CNEE filter-summary missing"
    assert os.path.isfile(manifest), "phyloP CNEE manifest missing"
    # No phastCons CNEE set this run (run_phastcons was False).
    assert not os.path.isdir(os.path.join(output_dir, "05-cnees", "phastcons"))
