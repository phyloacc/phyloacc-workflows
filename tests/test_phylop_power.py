#############################################################################
# Tier-1 unit tests for lib/phylop_power.py - the statistical-power ceiling
# behind the phyloP pre-flight gate (workflow/phylop_regions.smk:phylop_power_check).
#
# Uses two tiny self-contained PHAST .mod fixtures (tests/data/{shallow,deep}.mod,
# same JC rate matrix, differing only in branch lengths) so the tests need no
# external tools and no gitignored real data. Pure math - no Snakemake.
#############################################################################

import math
import os

import pytest

import lib.phylop_power as PP

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
SHALLOW = os.path.join(DATA, "shallow.mod")   # total tree length 0.5 subs/site
DEEP = os.path.join(DATA, "deep.mod")         # total tree length 16 subs/site


#############################################################################
# .mod parsing

def test_parse_mod_fields():
    m = PP.parse_mod(DEEP)
    assert m["Q"].shape == (4, 4)
    assert len(m["pi"]) == 4
    assert abs(float(m["pi"].sum()) - 1.0) < 1e-9
    assert PP.count_tips(PP.parse_newick(m["tree_newick"])) == 4  # fixture has 4 taxa
    # PHAST normalizes Q to 1 expected substitution/unit branch length.
    assert abs(m["expected_rate"] - 1.0) < 1e-6


def test_parse_mod_incomplete(tmp_path):
    bad = tmp_path / "bad.mod"
    bad.write_text("BACKGROUND: 0.25 0.25 0.25 0.25\nTREE: (A:0.1,B:0.1);\n")  # no RATE_MAT
    with pytest.raises(ValueError):
        PP.parse_mod(str(bad))


def test_tree_length():
    assert abs(PP.tree_length(PP.parse_newick(PP.parse_mod(SHALLOW)["tree_newick"])) - 0.5) < 1e-9
    assert abs(PP.tree_length(PP.parse_newick(PP.parse_mod(DEEP)["tree_newick"])) - 16.0) < 1e-9


#############################################################################
# Ceiling

def test_ceiling_reports_tree_length_and_tips():
    c = PP.ceiling(PP.parse_mod(DEEP))
    assert abs(c["tree_length"] - 16.0) < 1e-9
    assert c["tips"] == 4


def test_ceiling_increases_with_depth():
    # A longer neutral tree makes a perfectly conserved column more surprising.
    assert PP.ceiling(PP.parse_mod(DEEP))["neglog10p"] > PP.ceiling(PP.parse_mod(SHALLOW))["neglog10p"]


def test_ceiling_regression_values():
    # Pin the computed ceilings so a silent change in the math is caught.
    assert PP.ceiling(PP.parse_mod(SHALLOW))["neglog10p"] == pytest.approx(0.493, abs=0.01)
    assert PP.ceiling(PP.parse_mod(DEEP))["neglog10p"] == pytest.approx(2.405, abs=0.01)


#############################################################################
# LRT -> -log10 p

def test_neglog10_from_lrt_basics():
    assert PP.neglog10_from_lrt(0.0) == 0.0
    assert PP.neglog10_from_lrt(-5.0) == 0.0
    # p = erfc(sqrt(LRT/2)); LRT=2 -> p = erfc(1) = 0.15730
    assert PP.neglog10_from_lrt(2.0) == pytest.approx(-math.log10(math.erfc(1.0)), abs=1e-9)
    # monotonic increasing in LRT
    assert PP.neglog10_from_lrt(10.0) > PP.neglog10_from_lrt(5.0)


def test_neglog10_from_lrt_large_underflow():
    # Very large LRT underflows erfc to 0; the asymptotic branch still returns a finite value.
    v = PP.neglog10_from_lrt(2000.0)
    assert math.isfinite(v) and v > 100


#############################################################################
# FDR bar + gate

def test_fdr_threshold_value_and_monotonicity():
    assert PP.fdr_neglog10_threshold(1e8, 0.05) == pytest.approx(math.log10(1e8 / 0.05), abs=1e-9)
    assert PP.fdr_neglog10_threshold(1e8, 0.05) > PP.fdr_neglog10_threshold(1e6, 0.05)


def test_fdr_threshold_validation():
    with pytest.raises(ValueError):
        PP.fdr_neglog10_threshold(0, 0.05)
    with pytest.raises(ValueError):
        PP.fdr_neglog10_threshold(1e6, 1.5)


def test_power_gate_fails_on_shallow_tree():
    g = PP.power_gate(SHALLOW, m_sites=100, alpha=0.05)
    assert g["passes"] is False
    assert g["ceiling_neglog10p"] < g["threshold_neglog10p"]


def test_power_gate_passes_when_bar_is_low_enough():
    # deep.mod ceiling ~2.405 clears the bar only for a tiny M (log10(M/0.05) <= 2.405).
    assert PP.power_gate(DEEP, m_sites=10, alpha=0.05)["passes"] is True
    assert PP.power_gate(DEEP, m_sites=1_000_000, alpha=0.05)["passes"] is False


def test_power_gate_report_keys():
    g = PP.power_gate(DEEP, m_sites=10, alpha=0.05)
    for k in ("tree_length", "tips", "ceiling_neglog10p", "threshold_neglog10p", "m_sites", "alpha", "passes"):
        assert k in g
