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
import os
import sys
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from jinja2 import Template

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


def collect_cnee_filtering(m):
    # Reads cnees_from_conserved_chr's durable summary: raw phastCons CEs -> merged
    # (within cnee_ces_merge_gap_bp) -> any CE overlapping a CDS at all is dropped
    # entirely (no flanking fragments kept). The final post-length-filter CNEE count
    # is computed separately in collect_region_lengths (applying cnee_min_len_bp
    # directly to cnees.bed), since that's the same filter cnees_to_bed4_chr applies
    # but works regardless of cnee_output_format.
    summary_dir = m["paths"].get("cnees_summary_dir")
    if not summary_dir:
        return None
    rows = []
    for group, chrom in chrom_pairs(m["chromosome_groups"]):
        fpath = os.path.join(summary_dir, group, f"{chrom}.cnees-filter-summary.tsv")
        vals = read_metric_tsv(fpath)
        if not vals:
            continue
        rows.append({
            "group": group, "chrom": chrom,
            "ces_raw": int(vals.get("ces_raw", 0)),
            "ces_merged": int(vals.get("ces_merged", 0)),
            "ces_dropped_cds_overlap": int(vals.get("ces_dropped_cds_overlap", 0)),
            "cnees_after_cds_drop": int(vals.get("cnees_after_cds_drop", 0)),
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

    cnee_filter = collect_cnee_filtering(m)
    if cnee_filter is not None:
        ctx["cnee_filter_plot"] = bar_plot(
            cnee_filter, "chrom",
            [("ces_raw", "raw CEs"), ("ces_merged", "merged CEs"), ("cnees_after_cds_drop", "no CDS overlap")],
            "CE -> CNEE candidate funnel (before length filter)", "Count"
        )
        ctx["cnee_filter_totals"] = {
            "ces_raw": int(cnee_filter["ces_raw"].sum()),
            "ces_merged": int(cnee_filter["ces_merged"].sum()),
            "ces_dropped_cds_overlap": int(cnee_filter["ces_dropped_cds_overlap"].sum()),
            "cnees_after_cds_drop": int(cnee_filter["cnees_after_cds_drop"].sum()),
        }
        cnee_filter["pct_after_merge"] = (100 * cnee_filter["ces_merged"] / cnee_filter["ces_raw"]).round(1)
        cnee_filter["pct_kept_no_cds_overlap"] = (100 * cnee_filter["cnees_after_cds_drop"] / cnee_filter["ces_merged"]).round(1)
        ctx["cnee_filter_table"] = add_total_row(
            cnee_filter, sum_cols=["ces_raw", "ces_merged", "ces_dropped_cds_overlap", "cnees_after_cds_drop"],
            pct_specs=[
                ("pct_after_merge", "ces_merged", "ces_raw"),
                ("pct_kept_no_cds_overlap", "cnees_after_cds_drop", "ces_merged"),
            ]
        )

    neutral_models = collect_neutral_models(m)
    if neutral_models is not None:
        ctx["neutral_models_table"] = neutral_models
    ctx["avg_gc"] = collect_gc_summary(m)

    ce_rows = collect_region_lengths(m, "conserve_dir", lambda g, c: f"{c}.bed")
    ctx["ces"] = region_summary(ce_rows, "CEs")

    cnee_rows = collect_region_lengths(
        m, "cnees_dir", lambda g, c: f"{c}.cnees.bed", min_len_bp=m["paths"].get("cnee_min_len_bp")
    )
    ctx["cnees"] = region_summary(cnee_rows, "CNEEs")

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
