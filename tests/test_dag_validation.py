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
        [sys.executable, "-m", "snakemake", "-n", "-s", SNAKEFILE, "--configfile", str(config_path)],
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

# No "valid config successfully builds a full DAG" test here - reaching a genuinely
# complete DAG needs real maf/ref_gff/tree_file content (not just placeholder paths),
# which starts to mean building a small synthetic dataset - that's tier 3's job, not
# tier 2's. Tier 2 stays scoped to config-validation error paths (above), which are
# all reachable with placeholder paths alone since they raise before any real file
# needs to exist or be read.
