#############################################################################
# Unit tests for lib/parsing.py - Newick tip-name extraction, phastCons
# stderr rho parsing, and CNEE FASTA duplicate-species detection, extracted
# from workflow/*.smk rules' inline `run:` blocks. Run with `pytest tests/`
# from the repo root.
#############################################################################

import math

import pytest

import lib.parsing as parsing

#############################################################################
# parse_newick_tip_names

def test_parse_newick_tip_names_simple(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A,B,(C,D));")
    assert parsing.parse_newick_tip_names(str(tree)) == ["A", "B", "C", "D"]


def test_parse_newick_tip_names_quoted_and_branch_lengths(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("('Species A':0.1,B:0.2,(C:0.3,D:0.4)E:0.5);")
    # "E" is an internal node label, not a tip - shouldn't appear in the tip list.
    assert parsing.parse_newick_tip_names(str(tree)) == ["Species A", "B", "C", "D"]


def test_parse_newick_tip_names_comments_stripped(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A[a comment],B);")
    assert parsing.parse_newick_tip_names(str(tree)) == ["A", "B"]


def test_parse_newick_tip_names_duplicates_deduped(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(A,B,A);")
    assert parsing.parse_newick_tip_names(str(tree)) == ["A", "B"]

#############################################################################
# rho_from_phastcons_stderr

def test_rho_from_phastcons_stderr_finds_match():
    text = "some phastCons banner\nestimated rho = 0.452\nDone.\n"
    assert parsing.rho_from_phastcons_stderr(text) == 0.452


def test_rho_from_phastcons_stderr_takes_last_match():
    text = "rho = 0.1\n...\nrho = 0.9\n"
    assert parsing.rho_from_phastcons_stderr(text) == 0.9


def test_rho_from_phastcons_stderr_no_match_is_nan():
    assert math.isnan(parsing.rho_from_phastcons_stderr("no rho here"))


def test_rho_from_phastcons_stderr_empty_text_is_nan():
    assert math.isnan(parsing.rho_from_phastcons_stderr(""))

#############################################################################
# summarize_rho

def test_summarize_rho_basic_stats():
    # sorted: 0.1, 0.2, 0.3, 0.4, 0.5 - median=0.3, mean=0.3, p90 index = ceil(0.9*5)-1=4 -> 0.5
    summary = parsing.summarize_rho([0.5, 0.1, 0.3, 0.2, 0.4], stat="p90")
    assert summary["n"] == 5
    assert summary["median"] == pytest.approx(0.3)
    assert summary["mean"] == pytest.approx(0.3)
    assert summary["p90"] == pytest.approx(0.5)
    assert summary["selected"] == pytest.approx(0.5)


def test_summarize_rho_selects_requested_stat():
    summary = parsing.summarize_rho([0.1, 0.2, 0.3], stat="median")
    assert summary["selected"] == summary["median"]


def test_summarize_rho_filters_nan_and_nonpositive():
    summary = parsing.summarize_rho([float("nan"), -0.5, 0.0, 0.2, 0.4])
    assert summary["n"] == 2
    assert summary["mean"] == pytest.approx(0.3)


def test_summarize_rho_empty_list():
    summary = parsing.summarize_rho([])
    assert summary["n"] == 0
    assert math.isnan(summary["selected"])

#############################################################################
# filter_duplicate_species_fasta

def test_filter_duplicate_species_fasta_no_duplicates():
    lines = [">Mus_musculus:chr1-100-200\n", "ACGT\n", ">Rattus_norvegicus:chr1-100-200\n", "ACGT\n"]
    assert parsing.filter_duplicate_species_fasta(lines) is False


def test_filter_duplicate_species_fasta_detects_duplicate():
    lines = [
        ">Mus_musculus:chr1-100-200\n", "ACGT\n",
        ">Mus_musculus:chr2-300-400\n", "ACGT\n",
    ]
    assert parsing.filter_duplicate_species_fasta(lines) is True


def test_filter_duplicate_species_fasta_non_conforming_header_falls_back_to_full_token():
    # No genus_species-shaped prefix - falls back to using the whole species token as
    # the dedup key, so two distinct oddly-named entries are not flagged as duplicates.
    lines = [">weird123:chr1-1-2\n", "ACGT\n", ">weird456:chr1-1-2\n", "ACGT\n"]
    assert parsing.filter_duplicate_species_fasta(lines) is False
