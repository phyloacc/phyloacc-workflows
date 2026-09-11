#############################################################################
# Collects pipeline inputs/filtering/outputs into a single self-contained
# HTML summary report, reading only files the pipeline's own rules already
# produce (see analyses/resource-usage or the Snakefile's rule summary_report
# for what depends on what). Called from rule summary_report in Snakefile.
#
# Usage: python summary_report.py <manifest.json> <output.html>
#############################################################################

import base64
import io
import json
import math
import os
import sys
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from jinja2 import Template

# Invoked as a standalone script (python summary_report.py ...), so the repo root
# isn't on sys.path by default the way it is for the Snakefile/workflow/*.smk files -
# add it so lib.intervals is importable the same way they already do.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import lib.intervals as INTERVALS
import lib.phylop_power as PHYLOP_POWER

#############################################################################
# Small file readers

def read_lines(path):
    if not path or not os.path.isfile(path):
        return None
    with open(path) as fh:
        return [line.rstrip("\n") for line in fh if line.strip()]


def read_bed_lengths(path):
    lines = read_lines(path)
    if lines is None:
        return None
    lengths = []
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            lengths.append(int(fields[2]) - int(fields[1]))
        except ValueError:
            continue
    return lengths


def read_bed_positions(path):
    lines = read_lines(path)
    if lines is None:
        return None
    positions = []
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            positions.append((int(fields[1]), int(fields[2])))
        except ValueError:
            continue
    return positions


def chrom_pairs(chromosome_groups):
    return [(group, chrom) for group, chroms in chromosome_groups.items() for chrom in chroms]

#############################################################################
# Newick tree parsing (simplified - just tip count + total branch length,
# not the full robustness of the pipeline's own parser since it's only
# feeding a summary count here)

def parse_tree_stats(tree_str):
    tips = set()
    branch_lengths = []
    i, n = 0, len(tree_str)
    expect_label = False
    while i < n:
        ch = tree_str[i]
        if ch in "(,":
            expect_label = True
            i += 1
            continue
        if ch == ")":
            expect_label = False
            i += 1
            continue
        if ch == ";":
            break
        if ch == ":":
            i += 1
            start = i
            while i < n and tree_str[i] not in ",()[];":
                i += 1
            try:
                branch_lengths.append(float(tree_str[start:i]))
            except ValueError:
                pass
            continue
        if expect_label:
            start = i
            while i < n and tree_str[i] not in ":,()[];":
                i += 1
            label = tree_str[start:i].strip()
            if label:
                tips.add(label)
            expect_label = False
            continue
        i += 1
    return len(tips), sum(branch_lengths)

#############################################################################
# Collection functions - each returns None if its stage wasn't run/found

def read_prepost_tsv(fpath, var_name):
    # Parses the small "var / filter.cat / value" long-format summary TSVs
    # written by filter_maf_by_gap and filter_4d_sites (pre.filter/post.filter
    # rows for a given var). Returns (pre, post) as ints, or (None, None).
    if not os.path.isfile(fpath):
        return None, None
    vals = {}
    with open(fpath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3 or parts[0] != var_name:
                continue
            vals[parts[1]] = parts[2]
    return vals.get("pre.filter"), vals.get("post.filter")


def read_metric_tsv(fpath):
    # Parses the small "metric / value" TSVs written by cnees_from_conserved_chr.
    if not os.path.isfile(fpath):
        return {}
    vals = {}
    with open(fpath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2 or parts[0] == "metric":
                continue
            vals[parts[0]] = parts[1]
    return vals


def collect_maf_filtering(m):
    # Reads filter_maf_by_gap's durable summary TSV (maf_chunk_summary_dir), not
    # the manifest.txt/manifest.filtered.txt files themselves - those live under
    # the per-chromosome chunked-mafs directory, which phastcons_concat_chr
    # deletes wholesale once it's done (cleanup_chunk_intermediates), so they
    # won't survive to the end of a normal run.
    summary_dir = m["paths"].get("maf_chunk_summary_dir")
    if not summary_dir:
        return None
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        fpath = os.path.join(summary_dir, group, f"{chrom}.maf-chunk-filter-summary.tsv")
        total, kept = read_prepost_tsv(fpath, "num.chunks")
        if total is None:
            continue
        rows.append({
            "group": group, "chrom": chrom,
            "total_chunks": int(total),
            "kept_chunks": int(kept) if kept is not None else None,
        })
    return pd.DataFrame(rows) if rows else None


def collect_4d_filtering(m):
    summary_dir = m["paths"].get("neutral_summary_dir")
    threshold = m["paths"].get("filter_threshold_4d")
    prefix = m.get("maf_chr_prefix", "")
    if not summary_dir or threshold is None:
        return None
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        fpath = os.path.join(summary_dir, group, f"{prefix}{chrom}.4d-sites-filtered-{threshold}-summary.tsv")
        sites_pre, sites_post = read_prepost_tsv(fpath, "num.sites")
        if sites_pre is None:
            continue
        rows.append({
            "group": group, "chrom": chrom,
            "sites_pre": int(sites_pre),
            "sites_post": int(sites_post) if sites_post is not None else None,
        })
    return pd.DataFrame(rows) if rows else None


def collect_cnee_filtering(m, summary_dir_key="cnees_summary_dir", bed_dir_key=None):
    # Full CE/region -> CNEE filtering funnel, one row per chromosome. From
    # cnees_from_conserved_chr's durable summary: raw conserved elements -> merged (within
    # cnee_ces_merge_gap_bp) -> any element overlapping a CDS at all is dropped entirely.
    # Then the FINAL step (cnee_min_len_bp length filter, applied by cnees_to_bed4_chr) is
    # recomputed here directly from cnees.bed - the same length filter, but regardless of
    # cnee_output_format - so the funnel shows every step and its count, not just up to the
    # CDS drop. summary_dir_key/bed_dir_key select the source (phastcons vs phylop).
    summary_dir = m["paths"].get(summary_dir_key)
    if not summary_dir:
        return None
    bed_dir = m["paths"].get(bed_dir_key) if bed_dir_key else None
    min_len_bp = m["paths"].get("cnee_min_len_bp")
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        fpath = os.path.join(summary_dir, group, f"{chrom}.cnees-filter-summary.tsv")
        vals = read_metric_tsv(fpath)
        if not vals:
            continue
        after_cds = int(vals.get("cnees_after_cds_drop", 0))
        cnees_final = after_cds  # fallback if the bed isn't available
        if bed_dir:
            lengths = read_bed_lengths(os.path.join(bed_dir, group, f"{chrom}.cnees.bed"))
            if lengths is not None:
                cnees_final = sum(1 for l in lengths if (min_len_bp is None or l > min_len_bp))
        rows.append({
            "group": group, "chrom": chrom,
            "ces_raw": int(vals.get("ces_raw", 0)),
            "ces_merged": int(vals.get("ces_merged", 0)),
            "ces_dropped_cds_overlap": int(vals.get("ces_dropped_cds_overlap", 0)),
            "cnees_after_cds_drop": after_cds,
            "cnees_dropped_short": max(0, after_cds - cnees_final),
            "cnees_final": cnees_final,
        })
    return pd.DataFrame(rows) if rows else None


def collect_neutral_models(m):
    phylofit_dir = m["paths"].get("phylofit_dir")
    if not phylofit_dir:
        return None
    prefix = m.get("maf_chr_prefix", "")
    use_gc = m["flags"].get("use_gc_corrected_models")
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        if use_gc:
            fpath = os.path.join(phylofit_dir, group, f"{prefix}{chrom}-corrected.mod")
        else:
            fpath = os.path.join(phylofit_dir, group, "uncorrected-mods", f"{prefix}{chrom}.mod")
        if not os.path.isfile(fpath):
            continue
        tree_str = None
        with open(fpath) as fh:
            for line in fh:
                if line.startswith("TREE:"):
                    tree_str = line[len("TREE:"):].strip()
                    break
        n_species, total_branch_len = (None, None)
        if tree_str:
            n_species, total_branch_len = parse_tree_stats(tree_str)
        rows.append({
            "group": group, "chrom": chrom,
            "n_species": n_species,
            "total_branch_length": round(total_branch_len, 4) if total_branch_len is not None else None,
        })
    return pd.DataFrame(rows) if rows else None


def collect_gc_summary(m):
    if not m["flags"].get("use_gc_corrected_models"):
        return None
    avg_gc_file = m["paths"].get("avg_gc_file")
    if not avg_gc_file or not os.path.isfile(avg_gc_file):
        return None
    with open(avg_gc_file) as fh:
        try:
            return float(fh.read().strip())
        except ValueError:
            return None


def collect_region_lengths(m, dir_key, filename_fn, min_len_bp=None):
    base_dir = m["paths"].get(dir_key)
    if not base_dir:
        return None
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        fpath = os.path.join(base_dir, group, filename_fn(group, chrom))
        lengths = read_bed_lengths(fpath)
        if lengths is None:
            continue
        if min_len_bp is not None:
            lengths = [l for l in lengths if l > min_len_bp]
        rows.append({"group": group, "chrom": chrom, "n": len(lengths), "total_bp": sum(lengths), "lengths": lengths})
    return rows if rows else None


def collect_cnee_density(m, bin_bp, raw_dir_key, cnee_dir_key, raw_filename_fn):
    # Bins final-CNEE positions (and the raw conserved elements/regions they came from)
    # along each chromosome for the distribution plot. Source-parameterized: phastCons
    # (conserve_dir CEs + cnees_dir CNEEs) or phyloP (phylop_regions_dir + phylop_cnees_dir).
    # Final CNEEs are length-filtered by cnee_min_len_bp, matching the CNEE sections.
    maf_index_dir = m["paths"].get("maf_index_dir")
    maf_chr_prefix = m.get("maf_chr_prefix", "")
    raw_dir = m["paths"].get(raw_dir_key)
    cnee_dir = m["paths"].get(cnee_dir_key)
    min_len_bp = m["paths"].get("cnee_min_len_bp")
    if not maf_index_dir or not cnee_dir:
        return None
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        # Block index is co-located with the chromosome MAF and MAF-named (prefix + chrom).
        chrom_len = INTERVALS.read_chrom_length(os.path.join(maf_index_dir, group, f"{maf_chr_prefix}{chrom}.maf.block.idx"))
        if chrom_len is None:
            continue

        cnee_positions = read_bed_positions(os.path.join(cnee_dir, group, f"{chrom}.cnees.bed")) or []
        if min_len_bp is not None:
            cnee_positions = [(s, e) for s, e in cnee_positions if e - s > min_len_bp]

        raw_positions = None
        if raw_dir:
            raw_positions = read_bed_positions(os.path.join(raw_dir, group, raw_filename_fn(group, chrom)))

        n_bins = max(1, -(-chrom_len // bin_bp))  # ceil division, no float rounding

        def bin_counts(positions):
            counts = [0] * n_bins
            for s, _ in positions:
                counts[min(int(s // bin_bp), n_bins - 1)] += 1
            return counts

        rows.append({
            "group": group, "chrom": chrom, "chrom_len": chrom_len, "n_bins": n_bins,
            "cnee_counts": bin_counts(cnee_positions),
            "raw_counts": bin_counts(raw_positions) if raw_positions is not None else None,
        })
    return rows if rows else None

def read_count_tsv(fpath):
    # site_counts writes one "<chrom>\t<count>" row per chromosome file; sum the last
    # column so it's robust whether a file has one row or several.
    lines = read_lines(fpath)
    if lines is None:
        return None
    total = 0
    for line in lines:
        parts = line.split("\t")
        try:
            total += int(parts[-1])
        except ValueError:
            continue
    return total


def collect_power_gate(m):
    # Reads phylop_power_check's per-chromosome "metric/value" TSV: the LRT-based score
    # ceiling vs the multiple-testing bar log10(M/alpha) that decides whether per-site
    # phyloP has any power on this tree (see analyses/phylop-tree-length-power).
    power_dir = m["paths"].get("phylop_power_dir")
    if not power_dir:
        return None
    prefix = m.get("maf_chr_prefix", "")
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        vals = read_metric_tsv(os.path.join(power_dir, group, f"{prefix}{chrom}.power.tsv"))
        if not vals:
            continue

        def _num(key, cast=float, ndigits=None):
            try:
                v = cast(float(vals[key]))
            except (KeyError, ValueError):
                return None
            return round(v, ndigits) if (ndigits is not None) else v

        rows.append({
            "group": group, "chrom": chrom,
            "tips": _num("tips", int),
            "tree_length": _num("tree_length", ndigits=4),
            "ceiling": _num("ceiling_neglog10p", ndigits=3),
            "bar": _num("threshold_neglog10p", ndigits=3),
            "m_sites": _num("m_sites", int),
            "passes": str(vals.get("passes", "")).strip().lower() in ("true", "1", "yes"),
        })
    return pd.DataFrame(rows) if rows else None


def collect_site_counts(m):
    # Per-chromosome counts of FDR-significant conserved and accelerated sites.
    summary_dir = m["paths"].get("phylop_summary_dir")
    alpha = m["paths"].get("phylop_alpha")
    if not summary_dir or alpha is None:
        return None
    prefix = m.get("maf_chr_prefix", "")
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        base = os.path.join(summary_dir, group, f"{prefix}{chrom}")
        cons = read_count_tsv(f"{base}.conserved-site-counts.{alpha}.tsv")
        acc = read_count_tsv(f"{base}.accelerated-site-counts.{alpha}.tsv")
        if cons is None and acc is None:
            continue
        rows.append({"group": group, "chrom": chrom, "conserved": cons or 0, "accelerated": acc or 0})
    return pd.DataFrame(rows) if rows else None


def interval_overlap_stats(a, b):
    # a, b: iterables of (start, end). Returns bp overlap / union / Jaccard between the
    # two (merged) sets plus how many original intervals of each overlap the other set.
    def merge2(intervals):
        # merge overlapping (start, end) 2-tuples (INTERVALS.merge_intervals is for
        # (chrom, start, end) triples, so use a local 2-tuple version here)
        out = []
        for s, e in sorted(intervals):
            if out and s <= out[-1][1]:
                if e > out[-1][1]:
                    out[-1][1] = e
            else:
                out.append([s, e])
        return [(s, e) for s, e in out]

    a = [(s, e) for s, e in a if e > s]
    b = [(s, e) for s, e in b if e > s]
    ma = merge2(a)
    mb = merge2(b)
    i = j = inter = 0
    while i < len(ma) and j < len(mb):
        lo, hi = max(ma[i][0], mb[j][0]), min(ma[i][1], mb[j][1])
        if hi > lo:
            inter += hi - lo
        if ma[i][1] < mb[j][1]:
            i += 1
        else:
            j += 1
    bp_a = sum(e - s for s, e in ma)
    bp_b = sum(e - s for s, e in mb)
    union = bp_a + bp_b - inter

    def count_hits(orig, other):
        other = sorted(other)
        k = hits = 0
        for s, e in sorted(orig):
            while k < len(other) and other[k][1] <= s:
                k += 1
            if k < len(other) and other[k][0] < e:
                hits += 1
        return hits

    return {
        "n_a": len(a), "n_b": len(b),
        "n_a_hit": count_hits(a, mb), "n_b_hit": count_hits(b, ma),
        "a_bp": bp_a, "b_bp": bp_b,
        "overlap_bp": inter, "union_bp": union,
        "jaccard": (inter / union) if union else 0.0,
    }


def collect_cnee_overlap(m):
    # Concordance between the phastCons-source and phyloP-source CNEE sets (both must
    # exist). Length-filters both by cnee_min_len_bp so it matches the CNEE sections.
    pc_dir = m["paths"].get("cnees_dir")
    pp_dir = m["paths"].get("phylop_cnees_dir")
    if not pc_dir or not pp_dir:
        return None
    min_len = m["paths"].get("cnee_min_len_bp")
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        a = read_bed_positions(os.path.join(pc_dir, group, f"{chrom}.cnees.bed"))
        b = read_bed_positions(os.path.join(pp_dir, group, f"{chrom}.cnees.bed"))
        if a is None or b is None:
            continue
        if min_len is not None:
            a = [(s, e) for s, e in a if e - s > min_len]
            b = [(s, e) for s, e in b if e - s > min_len]
        st = interval_overlap_stats(a, b)
        rows.append({
            "group": group, "chrom": chrom,
            "n_phastcons": st["n_a"], "n_phylop": st["n_b"],
            "phastcons_hit": st["n_a_hit"], "phylop_hit": st["n_b_hit"],
            "phastcons_bp": st["a_bp"], "phylop_bp": st["b_bp"],
            "overlap_bp": st["overlap_bp"], "union_bp": st["union_bp"],
            "jaccard": round(st["jaccard"], 3),
        })
    return pd.DataFrame(rows) if rows else None

#############################################################################
# Plotting - each returns a base64-encoded PNG string, embedded directly
# into the HTML so the report has no external file dependencies.

def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def bar_plot(df, x_col, y_specs, title, ylabel):
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(df)), 4.5))
    x = range(len(df))
    width = 0.8 / len(y_specs)
    for i, (col, label) in enumerate(y_specs):
        ax.bar([xi + i * width for xi in x], df[col], width=width, label=label)
    ax.set_xticks([xi + width * (len(y_specs) - 1) / 2 for xi in x])
    ax.set_xticklabels(df[x_col], rotation=90, fontsize=7)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig_to_base64(fig)


def power_gate_nomogram(df, alpha, model_paths=None):
    # Detectability nomogram (see analyses/phylop-tree-length-power): a site is detectable
    # iff  ceiling(T) >= log10(M/alpha), i.e. below the boundary M_bound(T) = alpha *
    # 10^ceiling(T). The shaded band is the THEORETICAL closed form spanning kappa=1
    # (flat/JC) to kappa=1.5 (typical structure). When the fitted .mod is available
    # (model_paths), the ACTUAL fitted-model boundary for this run is drawn on top (dark
    # line, ceiling recomputed from the real model rescaled in T) - the definitive boundary
    # for these data, which resolves points that fall inside the generic band. Each
    # chromosome is a point at (its tree depth T, its M), colored by the gate result.
    tvals = [t for t in df["tree_length"].tolist() if t and t > 0]
    mvals = [m for m in df["m_sites"].tolist() if m and m > 0]
    if not tvals or not mvals:
        return None
    T_grid = np.geomspace(0.03, max(40.0, max(tvals) * 2.0), 220)
    boundary = lambda k: alpha * np.power(10.0, np.array([PHYLOP_POWER.neglog10_from_lrt(2.0 * k * t) for t in T_grid]))
    Mb1, Mb15 = boundary(1.0), boundary(1.5)
    ylo, yhi = min(min(mvals), 1e1) / 3.0, max(max(mvals), 1e8) * 3.0

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.fill_between(T_grid, ylo, Mb1, color="#c7dcef", alpha=0.7)          # detectable (any model)
    ax.fill_between(T_grid, Mb15, yhi, color="#f2c9c4", alpha=0.7)         # power-limited (any model)
    ax.fill_between(T_grid, Mb1, Mb15, color="#ece3b0", alpha=0.75)        # model-dependent band
    ax.plot(T_grid, Mb1, color="#8a6d1a", lw=1.3); ax.plot(T_grid, Mb15, color="#8a6d1a", lw=1.3)

    # fitted-model boundary/boundaries for this run (one per chromosome with a parseable .mod)
    fitted_drawn = False
    if model_paths:
        Tg2 = np.geomspace(T_grid[0], T_grid[-1], 60)
        for _, r in df.iterrows():
            mp = model_paths.get(r["chrom"])
            if not mp:
                continue
            try:
                md = PHYLOP_POWER.parse_mod(mp)
                tbase = PHYLOP_POWER.tree_length(PHYLOP_POWER.parse_newick(md["tree_newick"]))
                if tbase <= 0:
                    continue
                mb_real = alpha * np.power(10.0, np.array(
                    [PHYLOP_POWER.ceiling(md, rho_scale=t / tbase)["neglog10p"] for t in Tg2]))
            except Exception:
                continue
            ax.plot(Tg2, mb_real, color="#222222", lw=2.0, zorder=4,
                    label=("fitted model (this run)" if not fitted_drawn else None))
            fitted_drawn = True

    for _, r in df.iterrows():
        if not (r["tree_length"] and r["m_sites"]):
            continue
        passed = bool(r["passes"])
        ax.scatter([r["tree_length"]], [r["m_sites"]], s=90, zorder=5,
                   c=("#1f6ac0" if passed else "#c0392b"), edgecolor="black", linewidth=0.7,
                   marker=("o" if passed else "X"))
        ax.annotate(str(r["chrom"]), (r["tree_length"], r["m_sites"]),
                    textcoords="offset points", xytext=(7, 4), fontsize=8)
    ax.set_xlim(T_grid[0], T_grid[-1]); ax.set_ylim(ylo, yhi)
    ax.set_xlabel("neutral tree depth  T  (subs/site)")
    ax.set_ylabel("sites tested  M")
    ax.set_title("phyloP power: each chromosome vs the detectability limit")
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    handles = [
        Patch(fc="#c7dcef", label="detectable"),
        Patch(fc="#ece3b0", label="model-dependent (κ 1→1.5)"),
        Patch(fc="#f2c9c4", label="power-limited"),
    ]
    if fitted_drawn:
        handles.append(Line2D([0], [0], color="#222222", lw=2.0, label="fitted model (this run)"))
    handles += [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#1f6ac0", markeredgecolor="black", label="gate PASS"),
        Line2D([0], [0], marker="X", color="w", markerfacecolor="#c0392b", markeredgecolor="black", label="gate FAIL"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7.5, framealpha=0.92)
    fig.tight_layout()
    return fig_to_base64(fig)


def _lens_area(rA, rB, d):
    # area of intersection of two circles of radii rA, rB whose centers are distance d apart
    if d >= rA + rB:
        return 0.0
    if d <= abs(rA - rB):
        return math.pi * min(rA, rB) ** 2
    ca = max(-1.0, min(1.0, (d * d + rA * rA - rB * rB) / (2 * d * rA)))
    cb = max(-1.0, min(1.0, (d * d + rB * rB - rA * rA) / (2 * d * rB)))
    tri = max(0.0, (-d + rA + rB) * (d + rA - rB) * (d - rA + rB) * (d + rA + rB))
    return rA * rA * math.acos(ca) + rB * rB * math.acos(cb) - 0.5 * math.sqrt(tri)


def concordance_venn(phastcons_bp, phylop_bp, shared_bp):
    # Area-proportional 2-set base-pair Venn: each circle's AREA is that set's total bp and
    # the overlap lens AREA is the shared bp (solved for the center distance numerically).
    # bp overlap is symmetric so areas are unambiguous (element counts aren't - one CNEE can
    # overlap several).
    from matplotlib.patches import Circle
    rA = math.sqrt(max(phastcons_bp, 1) / math.pi)
    rB = math.sqrt(max(phylop_bp, 1) / math.pi)
    shared = max(0.0, min(shared_bp, math.pi * min(rA, rB) ** 2))
    lo, hi = abs(rA - rB) + 1e-9, rA + rB - 1e-9   # lens decreases as d grows
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _lens_area(rA, rB, mid) > shared:
            lo = mid
        else:
            hi = mid
    d = 0.5 * (lo + hi)
    xa, xb = 0.0, d

    def fmt(bp):
        return f"{bp / 1e6:.2f} Mb" if bp >= 1e6 else f"{bp / 1e3:.0f} kb"

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Circle((xa, 0), rA, facecolor="#4a90c2", alpha=0.5, edgecolor="#2a5a7a", lw=1.2))
    ax.add_patch(Circle((xb, 0), rB, facecolor="#e08a6b", alpha=0.5, edgecolor="#a0522d", lw=1.2))
    # region labels pushed to the outer crescents / overlap centre
    ax.text(xa - rA * 0.78, 0, f"phastCons only\n{fmt(max(0, phastcons_bp - shared_bp))}", ha="center", va="center", fontsize=8.5)
    ax.text(xb + rB * 0.78, 0, f"phyloP only\n{fmt(max(0, phylop_bp - shared_bp))}", ha="center", va="center", fontsize=8.5)
    if shared > 0:
        ax.text((xa + xb) / 2, 0, f"shared\n{fmt(shared_bp)}", ha="center", va="center", fontsize=8.5, weight="bold")
    span = max(rA, rB, d)
    ax.set_xlim(xa - rA - 0.18 * span, xb + rB + 0.18 * span)
    ax.set_ylim(-max(rA, rB) - 0.12 * span, max(rA, rB) + 0.12 * span)
    ax.set_title("CNEE overlap (area ∝ base pairs)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc="#4a90c2", alpha=0.5, label="phastCons CNEEs"),
                       Patch(fc="#e08a6b", alpha=0.5, label="phyloP CNEEs")],
              loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=9, frameon=False)
    return fig_to_base64(fig)


def box_plot(lengths, title):
    fig, ax = plt.subplots(figsize=(7, 1.8))
    flierprops = dict(marker=".", markersize=3, alpha=0.3)
    try:
        ax.boxplot(lengths, orientation="horizontal", flierprops=flierprops)
    except TypeError:
        ax.boxplot(lengths, vert=False, flierprops=flierprops)
    ax.set_xscale("log")
    ax.set_xlabel("Length (bp, log scale)")
    ax.set_yticks([])
    ax.set_title(title)
    fig.tight_layout()
    return fig_to_base64(fig)


def cnee_density_plot(rows, bin_bp, raw_label, source_label):
    # One combined figure, one row per chromosome. Final CNEEs (the subject) go on the
    # PRIMARY axis so they're always visible; the raw elements/regions they came from go
    # on a faint SECONDARY axis (their own scale) - otherwise the far more numerous raw
    # elements crush the CNEE bars to an invisible line.
    fig, axes = plt.subplots(nrows=len(rows), ncols=1, figsize=(10, 1.5 * len(rows)), squeeze=False)
    axes = axes[:, 0]
    cnee_handle = raw_handle = None
    for ax, row in zip(axes, rows):
        bin_starts_mb = [i * bin_bp / 1e6 for i in range(row["n_bins"])]
        if row["raw_counts"] is not None:
            ax2 = ax.twinx()
            raw_handle = ax2.bar(bin_starts_mb, row["raw_counts"], width=bin_bp / 1e6, align="edge",
                                 color="C1", alpha=0.30, zorder=1, label=f"{raw_label} (raw)")
            ax2.set_ylim(bottom=0)
            ax2.tick_params(axis="y", labelsize=6, colors="C1")
        cnee_handle = ax.bar(bin_starts_mb, row["cnee_counts"], width=bin_bp / 1e6, align="edge",
                             color="C0", zorder=3, label="CNEEs (final)")
        ax.set_ylim(bottom=0)
        ax.set_zorder(2); ax.patch.set_visible(False)  # keep CNEE bars above the twin-axis raw bars
        # Zoom the x-axis to the ALIGNED extent (first->last non-empty bin) rather than
        # 0..chrom_len, so a sub-chromosome window (e.g. a 2 Mb chunk of chr1) fills the
        # panel instead of hugging one edge; full-chromosome runs start near 0 unchanged.
        occupied = [i for i in range(row["n_bins"])
                    if row["cnee_counts"][i] or (row["raw_counts"] and row["raw_counts"][i])]
        if occupied:
            lo = occupied[0] * bin_bp / 1e6
            hi = min(row["chrom_len"], (occupied[-1] + 1) * bin_bp) / 1e6
            pad = 0.02 * (hi - lo) if hi > lo else bin_bp / 1e6
            ax.set_xlim(max(0, lo - pad), hi + pad)
        else:
            ax.set_xlim(0, row["chrom_len"] / 1e6)
        ax.set_ylabel(row["chrom"], rotation=0, ha="right", va="center", fontsize=8)
        ax.tick_params(axis="y", labelsize=6, colors="C0")
    axes[-1].set_xlabel("Position (Mb)")
    fig.suptitle(f"{source_label} CNEE density along each chromosome ({bin_bp // 1000}kb bins)\n"
                 f"blue = CNEEs (left axis) · orange = {raw_label} (right axis, own scale)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig_to_base64(fig)


def percentile(sorted_vals, pct):
    if not sorted_vals:
        return None
    k = (len(sorted_vals) - 1) * (pct / 100)
    f = int(k)
    c = min(f + 1, len(sorted_vals) - 1)
    if f == c:
        return sorted_vals[f]
    return sorted_vals[f] + (sorted_vals[c] - sorted_vals[f]) * (k - f)


def length_threshold_stats(lengths):
    # A long right tail makes a plain histogram hard to read, so alongside the
    # boxplot, quantify it directly: how many elements sit above a threshold
    # picked from the data itself (nearest round number to the 95th percentile).
    if not lengths:
        return None
    sorted_lengths = sorted(lengths)
    p95 = percentile(sorted_lengths, 95)
    threshold = max(50, round(p95 / 50) * 50)
    n_above = sum(1 for l in lengths if l > threshold)
    return {
        "threshold": threshold,
        "n_above": n_above,
        "pct_above": round(100 * n_above / len(lengths), 1),
    }

def add_total_row(df, sum_cols, pct_specs=None):
    # Appends a whole-genome TOTAL row: sums the given count columns, then
    # recomputes any percentage columns from those summed totals (summing
    # per-chromosome percentages directly wouldn't be meaningful).
    total = {"group": "TOTAL", "chrom": ""}
    for c in sum_cols:
        total[c] = df[c].sum()
    if pct_specs:
        for new_col, numer_col, denom_col in pct_specs:
            denom = total[denom_col]
            total[new_col] = round(100 * total[numer_col] / denom, 1) if denom else None
    return pd.concat([df, pd.DataFrame([total])], ignore_index=True)

#############################################################################

def region_summary(rows, label):
    if not rows:
        return None
    all_lengths = [l for r in rows for l in r["lengths"]]
    table = pd.DataFrame([{"group": r["group"], "chrom": r["chrom"], f"n_{label}": r["n"], "total_bp": r["total_bp"]} for r in rows])
    n_total = sum(r["n"] for r in rows)
    total_bp = sum(r["total_bp"] for r in rows)
    table = add_total_row(table, sum_cols=[f"n_{label}", "total_bp"])
    return {
        "table": table,
        "n_total": n_total,
        "total_bp": total_bp,
        "mean_len": round(sum(all_lengths) / len(all_lengths), 1) if all_lengths else None,
        "box": box_plot(all_lengths, f"{label} length distribution") if all_lengths else None,
        "threshold_stats": length_threshold_stats(all_lengths),
    }


def build_funnel(cnee_filter, raw_label, min_len_bp):
    # Per-source CE/region -> CNEE funnel context (stat totals + table + bar plot) showing
    # EVERY step and its count: raw -> merged (within cnee_ces_merge_gap_bp) -> after
    # CDS-overlap drop -> final (>= cnee_min_len_bp length filter). raw_label is what the
    # source's raw elements are called ("CEs" for phastCons, "phyloP regions" for phyloP).
    final_label = f"final (≥{min_len_bp}bp)" if min_len_bp is not None else "final"
    plot = bar_plot(
        cnee_filter, "chrom",
        [("ces_raw", f"raw {raw_label}"), ("ces_merged", "merged"),
         ("cnees_after_cds_drop", "−CDS overlap"), ("cnees_final", final_label)],
        f"{raw_label} → CNEE funnel (every step)", "Count",
    )
    totals = {k: int(cnee_filter[k].sum())
              for k in ("ces_raw", "ces_merged", "ces_dropped_cds_overlap",
                        "cnees_after_cds_drop", "cnees_dropped_short", "cnees_final")}
    cf = cnee_filter.copy()
    cf["pct_after_merge"] = (100 * cf["ces_merged"] / cf["ces_raw"]).round(1)
    cf["pct_kept_no_cds_overlap"] = (100 * cf["cnees_after_cds_drop"] / cf["ces_merged"]).round(1)
    cf["pct_kept_length"] = (100 * cf["cnees_final"] / cf["cnees_after_cds_drop"]).round(1)
    table = add_total_row(
        cf, sum_cols=["ces_raw", "ces_merged", "ces_dropped_cds_overlap",
                      "cnees_after_cds_drop", "cnees_dropped_short", "cnees_final"],
        pct_specs=[("pct_after_merge", "ces_merged", "ces_raw"),
                   ("pct_kept_no_cds_overlap", "cnees_after_cds_drop", "ces_merged"),
                   ("pct_kept_length", "cnees_final", "cnees_after_cds_drop")],
    )
    return {"raw_label": raw_label, "final_label": final_label, "totals": totals, "table": table, "plot": plot}


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <manifest.json> <output.html>", file=sys.stderr)
        sys.exit(1)

    manifest_path, output_path = sys.argv[1], sys.argv[2]
    with open(manifest_path) as fh:
        m = json.load(fh)

    ctx = {
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "output_dir": m["output_dir"],
        "snakemake_command": m.get("snakemake_command", ""),
        "main_inputs": sorted(m["main_inputs"].items()),
        "config_rows": sorted(m["config_display"].items()),
        "flags": m["flags"],
    }

    maf_filter = collect_maf_filtering(m)
    if maf_filter is not None:
        ctx["maf_filter_plot"] = bar_plot(
            maf_filter, "chrom", [("total_chunks", "total"), ("kept_chunks", "kept")],
            "MAF chunks kept after gap filtering", "Chunks"
        )
        ctx["maf_filter_total"] = int(maf_filter["total_chunks"].sum())
        ctx["maf_filter_kept"] = int(maf_filter["kept_chunks"].sum(skipna=True))
        maf_filter["pct_kept"] = (100 * maf_filter["kept_chunks"] / maf_filter["total_chunks"]).round(1)
        ctx["maf_filter_table"] = add_total_row(
            maf_filter, sum_cols=["total_chunks", "kept_chunks"],
            pct_specs=[("pct_kept", "kept_chunks", "total_chunks")]
        )

    sites_4d = collect_4d_filtering(m)
    if sites_4d is not None:
        plot_df = sites_4d.dropna(subset=["sites_pre", "sites_post"]).astype({"sites_pre": int, "sites_post": int}, errors="ignore")
        if not plot_df.empty:
            ctx["sites_4d_plot"] = bar_plot(
                plot_df, "chrom", [("sites_pre", "pre-filter"), ("sites_post", "post-filter")],
                "4d sites kept after missing-data filtering", "Sites"
            )
        sites_4d["pct_kept"] = (100 * sites_4d["sites_post"] / sites_4d["sites_pre"]).round(1)
        ctx["sites_4d_table"] = add_total_row(
            sites_4d, sum_cols=["sites_pre", "sites_post"],
            pct_specs=[("pct_kept", "sites_post", "sites_pre")]
        )

    # CNEE-producing sources present in this run (phastCons CEs and/or phyloP regions).
    # The funnel and the CNEE-output sections both iterate this list, so a phyloP-only
    # run shows just phyloP, a both-branches run shows the two side by side.
    cnee_sources = []
    if m["flags"].get("build_cnees"):
        if m["flags"].get("run_phastcons"):
            cnee_sources.append({"label": "phastCons", "raw_label": "CEs",
                                 "summary_key": "cnees_summary_dir", "bed_key": "cnees_dir"})
        if m["flags"].get("run_phylop"):
            cnee_sources.append({"label": "phyloP", "raw_label": "phyloP regions",
                                 "summary_key": "phylop_cnees_summary_dir", "bed_key": "phylop_cnees_dir"})

    funnels = []
    funnel_min_len = m["paths"].get("cnee_min_len_bp")
    for src in cnee_sources:
        cf = collect_cnee_filtering(m, src["summary_key"], src["bed_key"])
        if cf is not None:
            funnel = build_funnel(cf, src["raw_label"], funnel_min_len)
            funnel["label"] = src["label"]
            funnels.append(funnel)
    if funnels:
        ctx["cnee_funnels"] = funnels

    neutral_models = collect_neutral_models(m)
    if neutral_models is not None:
        ctx["neutral_models_table"] = neutral_models
    ctx["avg_gc"] = collect_gc_summary(m)

    ce_rows = collect_region_lengths(m, "conserve_dir", lambda g, c: f"{c}.bed")
    ctx["ces"] = region_summary(ce_rows, "CEs")

    # --- phyloP branch: power gate, FDR-significant sites, clustered conserved regions ---
    power_gate = collect_power_gate(m)
    if power_gate is not None:
        ctx["power_gate_all_pass"] = bool(power_gate["passes"].all())
        ctx["power_gate_n"] = int(len(power_gate))
        ctx["power_gate_n_pass"] = int(power_gate["passes"].sum())
        ctx["power_gate_enabled"] = bool(m["paths"].get("phylop_power_gate_enabled"))
        try:
            alpha_val = float(m["paths"].get("phylop_alpha"))
        except (TypeError, ValueError):
            alpha_val = 0.05
        # fitted-model .mod path per chromosome (same layout collect_neutral_models uses),
        # so the nomogram can draw the actual fitted-model boundary for this run.
        pfx = m.get("maf_chr_prefix", "")
        phylofit_dir = m["paths"].get("phylofit_dir")
        use_gc = m["flags"].get("use_gc_corrected_models")
        model_paths = {}
        if phylofit_dir:
            for _, r in power_gate.iterrows():
                if use_gc:
                    p = os.path.join(phylofit_dir, r["group"], f"{pfx}{r['chrom']}-corrected.mod")
                else:
                    p = os.path.join(phylofit_dir, r["group"], "uncorrected-mods", f"{pfx}{r['chrom']}.mod")
                if os.path.isfile(p):
                    model_paths[r["chrom"]] = p
        nomogram = power_gate_nomogram(power_gate, alpha_val, model_paths)
        if nomogram is not None:
            ctx["power_gate_plot"] = nomogram
        pg = power_gate.copy()
        pg["passes"] = pg["passes"].map({True: "yes", False: "NO"})
        ctx["power_gate_table"] = pg

    site_counts = collect_site_counts(m)
    if site_counts is not None:
        ctx["site_counts_totals"] = {
            "conserved": int(site_counts["conserved"].sum()),
            "accelerated": int(site_counts["accelerated"].sum()),
        }
        ctx["site_counts_plot"] = bar_plot(
            site_counts, "chrom", [("conserved", "conserved"), ("accelerated", "accelerated")],
            "phyloP FDR-significant sites", "Sites"
        )
        ctx["site_counts_table"] = add_total_row(site_counts, sum_cols=["conserved", "accelerated"])

    prefix = m.get("maf_chr_prefix", "")
    pp_region_rows = collect_region_lengths(m, "phylop_regions_dir", lambda g, c: f"{prefix}{c}.bed")
    ctx["phylop_regions"] = region_summary(pp_region_rows, "phyloP regions")
    ctx["phylop_cluster_method"] = m["paths"].get("phylop_cluster_method")

    # --- CNEE output sets, one per source (length-filtered by cnee_min_len_bp) ---
    cnee_outputs = []
    min_len = m["paths"].get("cnee_min_len_bp")
    for src in cnee_sources:
        rows = collect_region_lengths(m, src["bed_key"], lambda g, c: f"{c}.cnees.bed", min_len_bp=min_len)
        summ = region_summary(rows, "CNEEs")
        if summ is not None:
            cnee_outputs.append({"label": src["label"], "summary": summ})
    if cnee_outputs:
        ctx["cnee_outputs"] = cnee_outputs

    # --- phastCons vs phyloP CNEE concordance (Phase 3; both sets required) ---
    cnee_overlap = collect_cnee_overlap(m)
    if cnee_overlap is not None:
        tot_overlap = int(cnee_overlap["overlap_bp"].sum())
        tot_union = int(cnee_overlap["union_bp"].sum())
        tot_pc = int(cnee_overlap["n_phastcons"].sum())
        tot_pp = int(cnee_overlap["n_phylop"].sum())
        ctx["cnee_overlap_totals"] = {
            "jaccard": round(tot_overlap / tot_union, 3) if tot_union else 0.0,
            "pct_phastcons_w_phylop": round(100 * cnee_overlap["phastcons_hit"].sum() / tot_pc, 1) if tot_pc else None,
            "pct_phylop_w_phastcons": round(100 * cnee_overlap["phylop_hit"].sum() / tot_pp, 1) if tot_pp else None,
        }
        ctx["cnee_overlap_plot"] = concordance_venn(
            int(cnee_overlap["phastcons_bp"].sum()),
            int(cnee_overlap["phylop_bp"].sum()),
            tot_overlap,
        )
        ctx["cnee_overlap_bar"] = bar_plot(
            cnee_overlap, "chrom",
            [("n_phastcons", "phastCons"), ("n_phylop", "phyloP"), ("phastcons_hit", "phastCons ∩ phyloP")],
            "CNEE counts and overlap by chromosome", "CNEEs"
        )
        ctx["cnee_overlap_table"] = add_total_row(
            cnee_overlap, sum_cols=["n_phastcons", "n_phylop", "phastcons_hit", "phylop_hit", "overlap_bp", "union_bp"]
        )

    # --- CNEE distribution along chromosomes, one figure per source ---
    density_bin_bp = m["paths"].get("cnee_density_bin_bp") or 1_000_000
    density_figs = []
    if m["flags"].get("build_cnees"):
        if m["flags"].get("run_phastcons"):
            rows = collect_cnee_density(m, density_bin_bp, "conserve_dir", "cnees_dir", lambda g, c: f"{c}.bed")
            if rows:
                density_figs.append({"label": "phastCons",
                                     "plot": cnee_density_plot(rows, density_bin_bp, "CEs", "phastCons")})
        if m["flags"].get("run_phylop"):
            rows = collect_cnee_density(m, density_bin_bp, "phylop_regions_dir", "phylop_cnees_dir",
                                        lambda g, c: f"{prefix}{c}.bed")
            if rows:
                density_figs.append({"label": "phyloP",
                                     "plot": cnee_density_plot(rows, density_bin_bp, "regions", "phyloP")})
    if density_figs:
        ctx["cnee_density_bin_bp"] = density_bin_bp
        ctx["cnee_density_figs"] = density_figs

    template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates", "summary_report.html.j2")
    with open(template_path) as fh:
        template = Template(fh.read())

    def to_html_table(df):
        return df.to_html(index=False, na_rep="-", classes="data-table", border=0)

    html = template.render(ctx, table=to_html_table)
    with open(output_path, "w") as fh:
        fh.write(html)


if __name__ == "__main__":
    main()
