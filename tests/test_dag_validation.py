#############################################################################
# Tier-2 tests: does the Snakemake DAG build correctly, and are bad config
# values rejected with the expected error? Runs the real `snakemake -n`
# (dry-run - no external tools, no real execution) as a subprocess against a
# minimal config, rather than importing Python functions directly like
# tests/test_intervals.py / tests/test_parsing.py do.
#
# Baseline config values (maf/ref_gff/tree_file paths, etc.) are placeholder
# strings, not real files - every check here is reachable by config VALUE
# alone. The one exception (species parsed from a real tree file's content)
# is intentionally not covered here - see the plan file for why.
#############################################################################

import os
import subprocess
import sys

import yaml

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SNAKEFILE = os.path.join(REPO_ROOT, "Snakefile")

# Invoke the actual `snakemake` console-script entry point, not `python -m snakemake`.
# lib/common.py's pipelineSetup() classifies a run as a Snakemake-internal "worker"
# re-invocation (vs. the real top-level run) by checking for "__main__.py" in
# sys.argv[0] - which `python -m snakemake` always sets, since that's how Python's -m
# mechanics work for any package. That misclassifies every dry-run here as a "worker",
# which then requires a pre-existing state file from an earlier top-level run - fine by
# accident in a repo that already has one lying around, but breaks on a fresh checkout
# (confirmed: this is exactly what failed in CI). The installed `snakemake` executable's
# sys.argv[0] never contains "__main__.py", so it's always classified correctly.
SNAKEMAKE_EXE = os.path.join(os.path.dirname(sys.executable), "snakemake")

BASELINE = {
    "maf": "/placeholder/does-not-exist.maf",
    "maf_ref_id": "Species_name",
    "maf_chr_prefix": "",
    "maf_ref_chr_joiner": ".",
    "ref_gff": "/placeholder/does-not-exist.gff",
    "tree_file": "/placeholder/does-not-exist.tre",
    "ref_chromosome_groups": {"group1": ["chr1"]},
    "split_strategy": "num_seqs",
    "cnee_output_format": "none",
    "run_phylofit": False,
    "run_phylop": False,
    "run_phastcons": False,
    "rule_resources": {"default": {"partition": "test_partition", "mem_mb": 1000, "cpus": 1, "time": 10}},
}


def run_dryrun(tmp_path, overrides):
    config = dict(BASELINE)
    config.update(overrides)
    config["output_dir"] = str(tmp_path / "out")
    config["tmp_dir"] = str(tmp_path / "tmp")

    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    result = subprocess.run(
        [SNAKEMAKE_EXE, "-n", "-s", SNAKEFILE, "--configfile", str(config_path)],
        capture_output=True, text=True, cwd=REPO_ROOT,
    )
    # Which stream Snakemake writes a parse-time error to isn't stable across versions -
    # confirmed directly: the identical error goes to stderr on snakemake 9.2.0 but
    # stdout on 9.23.1. Check both combined rather than assume one.
    return result.returncode, result.stdout + result.stderr

#############################################################################
# Snakefile-level (no run_* flags needed - always evaluated)

def test_invalid_cnee_output_format(tmp_path):
    rc, output = run_dryrun(tmp_path, {"cnee_output_format": "bogus"})
    assert rc != 0
    assert "Invalid cnee_output_format" in output

#############################################################################
# workflow/phastcons_cnees.smk (run_phastcons=True, others False)

PHASTCONS_ONLY = {"run_phastcons": True, "run_phylofit": False, "run_phylop": False}


def test_invalid_split_strategy(tmp_path):
    overrides = {**PHASTCONS_ONLY, "split_strategy": "bogus"}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "Invalid split_strategy" in output


def test_ns_split_strategy_requires_ref_fasta(tmp_path):
    overrides = {**PHASTCONS_ONLY, "split_strategy": "ns"}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "requires ref_fasta to be set" in output


def test_ref_fasta_index_must_match_ref_fasta(tmp_path):
    overrides = {
        **PHASTCONS_ONLY, "split_strategy": "ns",
        "ref_fasta": "/placeholder/genome.fasta",
        "ref_fasta_index": "/placeholder/mismatched.fai",
    }
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "ref_fasta_index must match ref_fasta" in output


def test_ns_split_strategy_requires_picard_extension(tmp_path):
    overrides = {**PHASTCONS_ONLY, "split_strategy": "ns", "ref_fasta": "/placeholder/genome.notfasta"}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "Picard-recognized FASTA extension" in output


def test_num_seqs_bounds_must_be_nonnegative(tmp_path):
    overrides = {**PHASTCONS_ONLY, "num_seqs_min_gap_bp": -1}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "num_seqs_max_for_gap, num_seqs_min_gap_bp, and num_seqs_min_keep_region_len must be >= 0" in output


def test_window_overlap_bp_must_be_nonnegative(tmp_path):
    overrides = {**PHASTCONS_ONLY, "window_overlap_bp": -1}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "window_overlap_bp must be >= 0" in output


def test_fixed_windows_requires_window_size_bp(tmp_path):
    overrides = {**PHASTCONS_ONLY, "split_strategy": "fixed_windows", "ref_fasta": "/placeholder/genome.fasta"}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "requires window_size_bp in config" in output


def test_window_size_bp_must_be_positive(tmp_path):
    overrides = {
        **PHASTCONS_ONLY, "split_strategy": "fixed_windows",
        "ref_fasta": "/placeholder/genome.fasta", "window_size_bp": 0,
    }
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "window_size_bp must be > 0" in output


def test_window_overlap_bp_must_be_smaller_than_window_size_bp(tmp_path):
    overrides = {
        **PHASTCONS_ONLY, "split_strategy": "fixed_windows", "ref_fasta": "/placeholder/genome.fasta",
        "window_size_bp": 100, "window_overlap_bp": 100,
    }
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "window_overlap_bp must be smaller than window_size_bp" in output


def test_no_chromosomes_selected(tmp_path):
    overrides = {**PHASTCONS_ONLY, "target_ref_chromosomes": ["does-not-exist"]}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "No chromosomes selected" in output


def test_invalid_rho_mode(tmp_path):
    overrides = {**PHASTCONS_ONLY, "rho_mode": "bogus"}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "Invalid rho_mode" in output


def test_invalid_global_rho_stat(tmp_path):
    overrides = {**PHASTCONS_ONLY, "global_rho_stat": "bogus"}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "Invalid global_rho_stat" in output


def test_build_cnees_requires_ref_gff(tmp_path):
    overrides = {**PHASTCONS_ONLY, "build_cnees": True, "ref_gff": ""}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "build_cnees=true requires ref_gff to be set" in output


def test_cnee_output_format_requires_build_cnees(tmp_path):
    overrides = {**PHASTCONS_ONLY, "build_cnees": False, "cnee_output_format": "maf"}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "cnee_output_format requires build_cnees=true" in output


def test_cnee_ces_merge_gap_bp_must_be_nonnegative(tmp_path):
    overrides = {**PHASTCONS_ONLY, "cnee_ces_merge_gap_bp": -1}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "cnee_ces_merge_gap_bp must be >= 0" in output


def test_cnee_min_len_bp_must_be_nonnegative(tmp_path):
    overrides = {**PHASTCONS_ONLY, "cnee_min_len_bp": -1}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "cnee_min_len_bp must be >= 0" in output


def test_cnee_density_bin_bp_must_be_positive(tmp_path):
    overrides = {**PHASTCONS_ONLY, "cnee_density_bin_bp": 0}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "cnee_density_bin_bp must be > 0" in output


def test_fasta_output_with_no_species_info_at_all(tmp_path):
    # tree_file unset and no explicit species list - reachable without any real file,
    # since TREE_FILE is falsy and the "else: raise" branch never opens anything.
    overrides = {**PHASTCONS_ONLY, "cnee_output_format": "fasta", "tree_file": ""}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "requires species information" in output


def test_fasta_output_with_unparseable_tree_file(tmp_path):
    # tree_file is set and exists, but has no parseable tip labels - a minimal valid
    # (empty) Newick tree, just the terminating semicolon.
    tree_path = tmp_path / "empty.tre"
    tree_path.write_text(";")
    overrides = {**PHASTCONS_ONLY, "cnee_output_format": "fasta", "tree_file": str(tree_path)}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "No tip labels parsed from tree_file" in output

#############################################################################
# workflow/phylofit_models.smk (run_phylofit=True, others False)

PHYLOFIT_ONLY = {"run_phylofit": True, "run_phylop": False, "run_phastcons": False}


def test_phylofit_requires_maf(tmp_path):
    overrides = {**PHYLOFIT_ONLY, "maf": ""}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "run_phylofit=true requires maf to be set" in output


def test_phylofit_requires_ref_gff(tmp_path):
    overrides = {**PHYLOFIT_ONLY, "ref_gff": ""}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "run_phylofit=true requires ref_gff to be set" in output


def test_phylofit_requires_tree_file(tmp_path):
    overrides = {**PHYLOFIT_ONLY, "tree_file": ""}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "run_phylofit=true requires tree_file to be set" in output

#############################################################################
# workflow/phylop_regions.smk (run_phylop=True, others False)

PHYLOP_ONLY = {"run_phylop": True, "run_phylofit": False, "run_phastcons": False}


def test_phylop_no_chromosomes_selected(tmp_path):
    overrides = {**PHYLOP_ONLY, "target_ref_chromosomes": ["does-not-exist"]}
    rc, output = run_dryrun(tmp_path, overrides)
    assert rc != 0
    assert "No chromosomes selected" in output


# --- conserved-region clustering method selection (cluster_conserved_sites) ---

def test_invalid_phylop_cluster_method(tmp_path):
    rc, output = run_dryrun(tmp_path, {**PHYLOP_ONLY, "phylop_cluster_method": "bogus"})
    assert rc != 0
    assert "Invalid phylop_cluster_method" in output


def test_phylop_cluster_method_hdbscan_blocked(tmp_path):
    # hdbscan is implemented but blocked in the pipeline; the error must explain why (over-calls).
    rc, output = run_dryrun(tmp_path, {**PHYLOP_ONLY, "phylop_cluster_method": "hdbscan"})
    assert rc != 0
    assert "'hdbscan' is not supported" in output
    assert "false-positive" in output


def test_windowed_window_bp_must_be_positive(tmp_path):
    rc, output = run_dryrun(tmp_path, {**PHYLOP_ONLY, "phylop_cluster_method": "windowed",
                                       "windowed_window_bp": 0})
    assert rc != 0
    assert "windowed_window_bp must be > 0" in output


def test_hmm_probability_out_of_range(tmp_path):
    rc, output = run_dryrun(tmp_path, {**PHYLOP_ONLY, "phylop_cluster_method": "hmm",
                                       "hmm_t1_1": 1.5})
    assert rc != 0
    assert "hmm_t1_1 must be strictly between 0 and 1" in output


def test_valid_phylop_cluster_method_passes_validation(tmp_path):
    # A valid method must NOT trip the cluster-method validation. (The dry-run may still fail
    # later - a full phyloP DAG needs real phyloFit inputs - which is tier 3's job; here we only
    # assert our validation didn't fire.)
    rc, output = run_dryrun(tmp_path, {**PHYLOP_ONLY, "phylop_cluster_method": "hmm"})
    assert "Invalid phylop_cluster_method" not in output
    assert "is not supported" not in output

# --- phyloP power gate config validation (workflow/phylop_regions.smk) ---

def test_invalid_phylop_power_num_sites(tmp_path):
    rc, output = run_dryrun(tmp_path, {**PHYLOP_ONLY, "phylop_power_num_sites": "lots"})
    assert rc != 0
    assert "Invalid phylop_power_num_sites" in output


def test_phylop_power_num_sites_estimate_and_int_ok(tmp_path):
    # 'estimate' (default) and a positive int must both pass validation.
    _, out_est = run_dryrun(tmp_path, {**PHYLOP_ONLY, "phylop_power_num_sites": "estimate"})
    assert "Invalid phylop_power_num_sites" not in out_est
    _, out_int = run_dryrun(tmp_path, {**PHYLOP_ONLY, "phylop_power_num_sites": 5000000})
    assert "Invalid phylop_power_num_sites" not in out_int


# No "valid config successfully builds a full DAG" test here with placeholder paths -
# reaching a genuinely complete DAG needs real maf/ref_gff/tree_file content, which is
# what the fixture-backed structural checks below (and tier 3) use instead. Tier 2's
# placeholder-path tests stay scoped to config-validation error paths.

#############################################################################
# CNEE source fan-out (workflow/cnees.smk) - uses the tier-3 real fixture files so the
# DAG actually resolves, but only dry-runs and inspects the planned target paths.

INTEGRATION_CONFIG = os.path.join(REPO_ROOT, "tests", "integration", "data", "config.yaml")


def _dryrun_integration(tmp_path, overrides):
    with open(INTEGRATION_CONFIG) as f:
        config = yaml.safe_load(f)
    if not os.path.exists(config.get("maf", "")):
        import pytest
        pytest.skip("integration fixture MAF not present")
    config.update(overrides)
    config["output_dir"] = str(tmp_path / "out")
    config["tmp_dir"] = str(tmp_path / "tmp")
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    result = subprocess.run(
        [SNAKEMAKE_EXE, "-n", "-s", SNAKEFILE, "--configfile", str(config_path)],
        capture_output=True, text=True, cwd=REPO_ROOT,
    )
    return result.returncode, result.stdout + result.stderr


def test_both_sources_build_namespaced_cnees(tmp_path):
    # run_phylop + run_phastcons + build_cnees -> a CNEE set per source, namespaced.
    rc, output = _dryrun_integration(tmp_path, {"run_phylop": True, "run_phastcons": True})
    assert rc == 0, output[-3000:]
    assert "05-cnees/phastcons/bed/group1/CM000994.3.cnees.bed" in output
    assert "05-cnees/phylop/bed/group1/CM000994.3.cnees.bed" in output
    # CDS + MAF index are shared (computed once, not per source).
    assert "05-cnees/cds/group1/CM000994.3.cds.bed" in output


def test_phylop_only_builds_phylop_cnees_without_phastcons(tmp_path):
    # phyloP-only + build_cnees still builds CNEEs (and the shared MAF index) with no phastCons.
    rc, output = _dryrun_integration(tmp_path, {"run_phylop": True, "run_phastcons": False})
    assert rc == 0, output[-3000:]
    assert "05-cnees/phylop/bed/group1/CM000994.3.cnees.bed" in output
    assert "05-cnees/phastcons/" not in output
    assert "maf_index_chr" in output


def test_phylop_power_check_in_phylop_dag(tmp_path):
    # The pre-flight power gate is wired ahead of the phyloP scan.
    rc, output = _dryrun_integration(tmp_path, {"run_phylop": True, "run_phastcons": False})
    assert rc == 0, output[-3000:]
    assert "phylop_power_check" in output
