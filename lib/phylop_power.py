#############################################################################
# Statistical-power ceiling for per-site phyloP conservation calling.
#
# Per-site phyloP (LRT, CONACC) has a hard detection ceiling set by the total
# neutral tree length: a perfectly conserved column (zero substitutions) is only
# as surprising as the number of substitutions expected across the tree. On a
# shallow tree even the best-possible site cannot clear a genome/chromosome-wide
# FDR threshold, so phyloP returns ~zero conserved sites by construction (not a
# bug). This module computes that ceiling directly from a fitted neutral model
# (PHAST .mod) and compares it to the FDR bar, so the pipeline can gate the
# phyloP stage before spending the (expensive) per-chromosome phyloP scan.
#
# The math is the importable, pipeline-side copy of the derivation in
# analyses/phylop-tree-length-power/phylop_power.py (kept in sync by hand; the
# analysis version stays frozen as the research artifact). numpy-only.
#############################################################################

import math
import re

import numpy as np

BASES = "ACGT"


#############################################################################
# PHAST .mod parsing

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


#############################################################################
# Newick parsing -> nested tree; leaves have no children

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
        m = re.match(r"[^():,]+", s[pos:])
        if m:
            token = m.group(0)
            pos += len(token)
            node.name = token
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m = re.match(r"[0-9eE.+\-]+", s[pos:])
            node.length = float(m.group(0))
            pos += len(m.group(0))
        return node

    return parse_clade()


def tree_length(node):
    return node.length + sum(tree_length(c) for c in node.children)


def count_tips(node):
    if not node.children:
        return 1
    return sum(count_tips(c) for c in node.children)


#############################################################################
# Substitution model: P(t) = expm(Q t) via eigendecomposition

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


def _prune(node, model, base_idx, rho, is_root=False):
    if not node.children:  # leaf: observed state = base_idx
        L = np.zeros(4)
        L[base_idx] = 1.0
        return L
    L = np.ones(4)
    for child in node.children:
        cp = _prune(child, model, base_idx, rho)
        P = model.P(child.length * rho)
        L = L * (P @ cp)
    return L


def loglik_conserved(node, model, base_idx, rho):
    """log P(all leaves = base_idx) with branch lengths scaled by rho."""
    partial = _prune(node, model, base_idx, rho, is_root=True)
    total = float(np.sum(model.pi * partial))
    return math.log(total) if total > 0 else -math.inf


#############################################################################
# Conservation-score ceiling

def neglog10_from_lrt(lrt):
    """phyloP CONACC conserved p-value from the LRT, as -log10(p).

    Empirically validated against the real phyloP binary on a fully-conserved
    probe column: phyloP uses the full chi^2_1 upper tail,
    p = P(chi^2_1 >= LRT) = 2*Phi(-sqrt(LRT)) = erfc(sqrt(LRT/2)) - i.e. the
    two-sided CONACC convention, NOT the one-sided Phi(-sqrt(LRT)).
    """
    if lrt <= 0:
        return 0.0
    p = math.erfc(math.sqrt(lrt / 2.0))  # = 2*Phi(-sqrt(LRT))
    if p <= 0:
        # underflow for very large LRT: log-domain asymptotic of the normal tail.
        x = math.sqrt(lrt)
        log10p = math.log10(2) - (x * x) / (2 * math.log(10)) - math.log10(x * math.sqrt(2 * math.pi))
        return -log10p
    return -math.log10(p)


def ceiling(model_dict, rho_scale=1.0):
    """Max achievable conserved -log10(p) for this model (tree scaled by rho_scale).

    Returns a dict with tree_length (of the scaled tree), LRT_max, neglog10p
    ceiling, base, tips, and kappa = LRT_max / tree_length.
    """
    model = SubstModel(model_dict["Q"], model_dict["pi"])
    root = parse_newick(model_dict["tree_newick"])
    T = tree_length(root) * rho_scale

    best = None
    for b in range(4):
        logL0 = math.log(model_dict["pi"][b])                 # rho -> 0 : P(all b) = pi_b
        logL1 = loglik_conserved(root, model, b, rho_scale)   # neutral (scaled) tree
        lrt = 2.0 * (logL0 - logL1)
        nlp = neglog10_from_lrt(lrt)
        if best is None or nlp > best["neglog10p"]:
            best = {"base": BASES[b], "LRT_max": lrt, "neglog10p": nlp}

    best["tree_length"] = T
    best["tips"] = count_tips(root)
    best["kappa"] = best["LRT_max"] / T if T > 0 else float("nan")
    return best


#############################################################################
# Power gate: can the best-possible conserved site clear FDR?

def fdr_neglog10_threshold(m_sites, alpha):
    """-log10 of the FDR bar the most significant site must clear.

    Benjamini-Hochberg at the top rank requires p_(1) <= alpha / M, so in
    -log10 space the site score must reach log10(M / alpha). (This is also the
    Bonferroni bar for the single most-significant test - the two coincide at
    rank 1, which is the only rank a single best-possible site can occupy.)
    """
    if m_sites <= 0:
        raise ValueError(f"m_sites must be > 0; got {m_sites}")
    if not (0.0 < alpha < 1.0):
        raise ValueError(f"alpha must be in (0, 1); got {alpha}")
    return math.log10(m_sites / alpha)


def power_gate(mod_path, m_sites, alpha):
    """Evaluate the phyloP power gate for one fitted neutral model.

    Returns a dict describing the decision; `passes` is True iff the
    best-possible conserved site can clear the FDR bar for `m_sites` tests.
    """
    model = parse_mod(mod_path)
    c = ceiling(model)
    thr = fdr_neglog10_threshold(m_sites, alpha)
    return {
        "tree_length": c["tree_length"],
        "tips": c["tips"],
        "ceiling_neglog10p": c["neglog10p"],
        "threshold_neglog10p": thr,
        "m_sites": int(m_sites),
        "alpha": alpha,
        "passes": c["neglog10p"] >= thr,
    }
