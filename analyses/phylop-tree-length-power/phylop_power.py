#!/usr/bin/env python3
"""
phylop_power.py - analytic backbone for the phyloP tree-length / power analysis.

The question this supports: why does per-site phyloP (LRT, CONACC) + genome-wide
BH-FDR return ~zero conserved sites on shallow-tree datasets (hamsters, birds,
turtles; total tree length ~1 subs/site) but millions on the 241-taxon Zoonomia
mammal alignment (tree length ~16)?

Core claim: the *maximum achievable* per-site conservation score is capped by the
total neutral branch length. A perfectly conserved column (zero substitutions) is
only as surprising as the number of substitutions you *expected* to see across the
tree. This module computes that ceiling directly from the real fitted neutral
models (PHAST .mod files), with no simulated data - it only evaluates each model's
own likelihood for a fully-conserved column.

phyloP's per-site score is -log10(p); the p-value comes from a one-sided LRT on a
single branch-scaling parameter rho (rho=1 neutral, rho<1 conserved). Under the
null LRT ~ chi^2_1, so the one-sided conserved p-value is p = Phi(-sqrt(LRT)).
For a fully-conserved all-b column the MLE is rho_hat -> 0 (fewer substitutions is
always more likely), giving

    LRT_max(b) = 2 * ( log pi_b - logL_neutral(all leaves = b) )

and the site-score ceiling is max over b of -log10 Phi(-sqrt(LRT_max(b))).

Everything here uses numpy only (matrix exponentials via eigendecomposition of the
reversible rate matrix), so it runs in the repo's phyloacc-workflows env with no
extra installs.
"""

import math
import re
import sys

import numpy as np

BASES = "ACGT"


# ---------------------------------------------------------------------------
# .mod parsing
# ---------------------------------------------------------------------------

def parse_mod(path):
    """Parse a PHAST .mod file into background freqs, rate matrix, and tree."""
    pi = None
    Q = None
    tree = None
    with open(path) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("BACKGROUND:"):
            pi = np.array([float(x) for x in line.split(":", 1)[1].split()])
        elif line.startswith("RATE_MAT:"):
            rows = []
            for j in range(1, 5):
                rows.append([float(x) for x in lines[i + j].split()])
            Q = np.array(rows)
            i += 4
        elif line.startswith("TREE:"):
            tree = line.split(":", 1)[1].strip()
        i += 1

    if pi is None or Q is None or tree is None:
        raise ValueError(f"Incomplete .mod file (missing BACKGROUND/RATE_MAT/TREE): {path}")

    # PHAST normalizes Q so branch lengths are in expected subs/site: -sum pi_i Q_ii = 1.
    exp_rate = -float(np.sum(pi * np.diag(Q)))
    return {"pi": pi, "Q": Q, "tree_newick": tree, "expected_rate": exp_rate}


# ---------------------------------------------------------------------------
# Newick parsing -> nested tree of (branch_length, [children]); leaves have []
# ---------------------------------------------------------------------------

class Node:
    __slots__ = ("length", "children", "name")

    def __init__(self, length=0.0, name=None):
        self.length = length
        self.children = []
        self.name = name


def parse_newick(s):
    s = s.strip().rstrip(";").strip()
    pos = 0

    def parse_clade():
        nonlocal pos
        node = Node()
        if s[pos] == "(":
            pos += 1  # consume '('
            while True:
                child = parse_clade()
                node.children.append(child)
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        # read optional name (leaf label or internal label)
        m = re.match(r"[^():,]+", s[pos:])
        if m:
            token = m.group(0)
            pos += len(token)
            # token may be "name" or "name:len" - but ':' is handled next; PHAST puts name then :len
            node.name = token
        # read optional branch length
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m = re.match(r"[0-9eE.+\-]+", s[pos:])
            node.length = float(m.group(0))
            pos += len(m.group(0))
        return node

    root = parse_clade()
    return root


def tree_length(node):
    return node.length + sum(tree_length(c) for c in node.children)


def count_tips(node):
    if not node.children:
        return 1
    return sum(count_tips(c) for c in node.children)


# ---------------------------------------------------------------------------
# Substitution model: P(t) = expm(Q t) via eigendecomposition
# ---------------------------------------------------------------------------

class SubstModel:
    def __init__(self, Q, pi):
        self.Q = Q
        self.pi = pi
        lam, V = np.linalg.eig(Q)
        self.lam = lam.real
        self.V = V.real
        self.Vinv = np.linalg.inv(V).real

    def P(self, t):
        if t <= 0:
            return np.eye(4)
        return (self.V * np.exp(self.lam * t)) @ self.Vinv


# ---------------------------------------------------------------------------
# Felsenstein pruning for a fully-conserved (all leaves = base b) column,
# with all branch lengths scaled by rho.
# ---------------------------------------------------------------------------

def loglik_conserved(node, model, base_idx, rho):
    """log P(all leaves = base_idx) with branch lengths scaled by rho."""
    partial = _prune(node, model, base_idx, rho, is_root=True)
    total = float(np.sum(model.pi * partial))
    return math.log(total) if total > 0 else -math.inf


def _prune(node, model, base_idx, rho, is_root=False):
    if not node.children:  # leaf: observed state = base_idx
        L = np.zeros(4)
        L[base_idx] = 1.0
        return L
    # internal: product over children of (P(branch) @ child_partial)
    L = np.ones(4)
    for child in node.children:
        cp = _prune(child, model, base_idx, rho)
        P = model.P(child.length * rho)
        L = L * (P @ cp)
    return L


# ---------------------------------------------------------------------------
# Conservation-score ceiling
# ---------------------------------------------------------------------------

def neglog10_from_lrt(lrt):
    """phyloP CONACC conserved p-value from the LRT.

    Empirically validated against the real phyloP binary on a fully-conserved
    probe column (see analyses note): phyloP uses the full chi^2_1 upper tail,
    p = P(chi^2_1 >= LRT) = 2*Phi(-sqrt(LRT)) = erfc(sqrt(LRT/2)) - i.e. the
    two-sided CONACC convention, NOT the one-sided Phi(-sqrt(LRT)). Using the
    one-sided form overestimates -log10(p) by log10(2) ~ 0.301.
    """
    if lrt <= 0:
        return 0.0
    p = math.erfc(math.sqrt(lrt / 2.0))  # = 2*Phi(-sqrt(LRT))
    if p <= 0:
        # underflow for very large LRT: log-domain asymptotic of the normal tail,
        # including the factor of 2 (=> +log10(2) relative to the one-sided form).
        x = math.sqrt(lrt)
        log10p = math.log10(2) - (x * x) / (2 * math.log(10)) - math.log10(x * math.sqrt(2 * math.pi))
        return -log10p
    return -math.log10(p)


def invariant_probability(model_dict, rho_scale=1.0):
    """P1 = probability a NEUTRAL single site is invariant (identical in all species),
    with the tree scaled by rho_scale. This is the exact quantity the 29-mammals paper
    (Lindblad-Toh et al. 2011) calls P1: P1 = sum_b P(all leaves = b | neutral tree).
    Computed from the real fitted model via Felsenstein pruning - no approximation.
    """
    model = SubstModel(model_dict["Q"], model_dict["pi"])
    root = parse_newick(model_dict["tree_newick"])
    return sum(math.exp(loglik_conserved(root, model, b, rho_scale)) for b in range(4))


def ceiling(model_dict, rho_scale=1.0):
    """Max achievable conserved score for this model with tree scaled by rho_scale.

    Returns dict with tree_length (of the scaled tree), LRT_max, neglog10p ceiling,
    and kappa = LRT_max / tree_length.
    """
    model = SubstModel(model_dict["Q"], model_dict["pi"])
    root = parse_newick(model_dict["tree_newick"])
    T_base = tree_length(root)
    T = T_base * rho_scale

    best = None
    for b in range(4):
        logL0 = math.log(model_dict["pi"][b])            # rho -> 0 : P(all b) = pi_b
        logL1 = loglik_conserved(root, model, b, rho_scale)  # neutral (scaled) tree
        lrt = 2.0 * (logL0 - logL1)
        nlp = neglog10_from_lrt(lrt)
        if best is None or nlp > best["neglog10p"]:
            best = {"base": BASES[b], "LRT_max": lrt, "neglog10p": nlp}

    best["tree_length"] = T
    best["tips"] = count_tips(root)
    best["kappa"] = best["LRT_max"] / T if T > 0 else float("nan")
    return best


# ---------------------------------------------------------------------------
# CLI: cross-check against the real turtle observed max score
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    mods = sys.argv[1:]
    if not mods:
        # default cross-check: the turtle model whose real phylop.bed had max
        # conserved score 1.746 (p=0.018) on N=6.32e6 sites.
        mods = ["data/turtles/small-test2/03-phylofit/shortest-scaffolds/NC_050095.1-corrected.mod"]

    print(f"{'tips':>5} {'tree_len':>9} {'LRT_max':>8} {'kappa':>6} {'ceiling_-log10p':>15}  model")
    for m in mods:
        md = parse_mod(m)
        c = ceiling(md)
        print(f"{c['tips']:>5} {c['tree_length']:>9.3f} {c['LRT_max']:>8.3f} "
              f"{c['kappa']:>6.2f} {c['neglog10p']:>15.3f}  {m}")
