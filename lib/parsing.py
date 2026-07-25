#############################################################################
# Smaller, more heterogeneous parsing helpers shared across workflow/*.smk
# rules - Newick tip-name extraction, phastCons stderr rho parsing, and CNEE
# FASTA duplicate-species detection. Extracted from those rules' inline
# `run:` blocks so the logic is importable and unit-testable outside of a
# live Snakemake job (see tests/test_parsing.py) - a pure refactor, not a
# behavior change.
#############################################################################

import math
import re

#############################################################################

def parse_newick_tip_names(tree_path: str):
    with open(tree_path, "r", encoding="utf-8") as fp:
        data = fp.read()

    tips = []
    seen = set()
    i = 0
    expect_label = False
    n = len(data)

    while i < n:
        ch = data[i]

        if ch in "(,":
            expect_label = True
            i += 1
            continue

        if ch in " \t\r\n":
            i += 1
            continue

        if ch == "[":
            i += 1
            while i < n and data[i] != "]":
                i += 1
            if i < n:
                i += 1
            continue

        if ch == ")":
            expect_label = False
            i += 1
            continue

        if ch == ";":
            break

        if expect_label:
            if ch in "\'\"":
                quote = ch
                i += 1
                start = i
                while i < n and data[i] != quote:
                    i += 1
                label = data[start:i].strip()
                if i < n:
                    i += 1
            else:
                start = i
                while i < n and data[i] not in ":,()[];":
                    i += 1
                label = data[start:i].strip()

            if label and label not in seen:
                tips.append(label)
                seen.add(label)
            expect_label = False
            continue

        i += 1

    return tips

#############################################################################

RHO_RE = re.compile(r"rho\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")


def rho_from_phastcons_stderr(text):
    # Returns the last "rho = X.XX"-style match in a phastCons --estimate-rho stderr
    # capture, or NaN if there's no match (or the run failed before printing one).
    matches = RHO_RE.findall(text or "")
    if not matches:
        return float("nan")
    try:
        return float(matches[-1])
    except ValueError:
        return float("nan")


def summarize_rho(values, stat="p90"):
    # values: iterable of per-chunk rho estimates (may include NaN/non-positive
    # values from failed chunks - those are filtered out here, not by the caller).
    # Returns a dict with mean/median/p90 and the one selected by `stat`.
    vals = sorted(v for v in values if math.isfinite(v) and v > 0.0)
    if not vals:
        return {"mean": float("nan"), "median": float("nan"), "p90": float("nan"),
                 "selected": float("nan"), "n": 0}
    n = len(vals)
    idx = int(math.ceil(0.9 * n) - 1)
    idx = max(0, min(idx, n - 1))
    p90 = vals[idx]
    median = float(vals[n // 2]) if n % 2 == 1 else float((vals[n // 2 - 1] + vals[n // 2]) / 2.0)
    mean = float(sum(vals) / n)
    selected = {"p90": p90, "median": median, "mean": mean}.get(stat, mean)
    return {"mean": mean, "median": median, "p90": p90, "selected": selected, "n": n}

#############################################################################

SPECIES_KEY_RE = re.compile(r"^([A-Za-z]+)_([A-Za-z]+)")


def filter_duplicate_species_fasta(fasta_lines):
    # fasta_lines: iterable of lines from a CNEE FASTA file (as produced by mafutils
    # fetch -f). Returns True if two or more headers share the same genus_species key
    # (first two underscore-separated tokens of the species name before ":"), meaning
    # this alignment has a duplicated species and should be dropped.
    seen = set()
    for line in fasta_lines:
        if not line.startswith(">"):
            continue
        token = line[1:].strip().split()[0]
        species = token.split(":", 1)[0]
        m = SPECIES_KEY_RE.match(species)
        species_key = f"{m.group(1)}_{m.group(2)}" if m else species
        if species_key in seen:
            return True
        seen.add(species_key)
    return False
