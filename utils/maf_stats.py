#############################################################################
# Summarize an indexed MAF quickly at overall, species, and block levels.
#
# Gregg Thomas + Codex, March 2026
#############################################################################

"""
maf_stats.py

Compute summary statistics from a MAF file using an existing block index
produced by maf_index.py.

Outputs:
    <prefix>.overall.tsv  : one metric per line
    <prefix>.species.tsv  : one row per species
    <prefix>.block.tsv    : one row per block (optional; default enabled)
    <prefix>.dashboard.html : self-contained summary dashboard (optional)

Designed for speed:
    - Parallel by index block chunks with ProcessPoolExecutor
    - One file handle per worker
    - Optional temp-file based block-row writing to keep RAM usage low
"""

import argparse
import base64
import datetime
import gzip
import html
import io
import logging
import random
import os
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import lib.common as COMMON
import lib.loginit as LOGINIT

#############################################################################

SITE_GAP_HIST_BINS = 20

#############################################################################


def optParse():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize an indexed MAF at overall/species/block levels using "
            "parallel block processing."
        )
    )
    parser.add_argument("maf_file", help="Input MAF file (.maf or .maf.gz)")
    parser.add_argument("index_file", help="Block index from maf_index.py")
    parser.add_argument(
        "--output-prefix",
        "-o",
        default="maf_stats",
        help="Output prefix/path (default: maf_stats)",
    )
    parser.add_argument(
        "--processes",
        "-p",
        type=int,
        default=1,
        help="Number of worker processes (default: 1)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
        help="Blocks per worker task (default: 5000)",
    )
    parser.add_argument(
        "--no-block-table",
        action="store_true",
        default=False,
        help="Skip writing the per-block output table for faster/lighter runs.",
    )
    parser.add_argument(
        "--expected-species",
        default="",
        help="Comma-separated species names for exact per-block missing lists.",
    )
    parser.add_argument(
        "--expected-species-file",
        default="",
        help="File with one species name per line for exact missing lists.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--html-dashboard",
        action="store_true",
        default=False,
        help="Write an HTML dashboard with summary metrics and species plots.",
    )
    parser.add_argument(
        "--dashboard-top-species",
        type=int,
        default=25,
        help="Number of species to show in each dashboard bar plot (default: 25).",
    )
    parser.add_argument(
        "--dashboard-max-block-points",
        type=int,
        default=50000,
        help="Maximum block rows sampled for dashboard block plots (default: 50000).",
    )
    return parser.parse_args()


#############################################################################


def parseExpectedSpecies(args):
    species = set()
    if args.expected_species:
        for sp in args.expected_species.split(","):
            sp = sp.strip()
            if sp:
                species.add(sp)

    if args.expected_species_file:
        with open(args.expected_species_file, "r", encoding="utf-8") as fp:
            for line in fp:
                line = line.strip()
                if line and not line.startswith("#"):
                    species.add(line)
    return species


def parseIndex(index_file, LOG):
    entries = []
    with open(index_file, "r", encoding="utf-8") as fp:
        block_id = 0
        for line in fp:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                LOG.warning(f"Skipping malformed index line: {line.strip()}")
                continue
            try:
                block_id += 1
                entries.append(
                    (
                        block_id,
                        fields[0],  # scaffold
                        int(fields[1]),  # ref_start
                        int(fields[2]),  # ref_length
                        int(fields[3]),  # aln_length
                        int(fields[5]),  # n_seq_lines (index)
                        int(fields[6]),  # offset_start
                        int(fields[7]),  # offset_end
                    )
                )
            except ValueError:
                LOG.warning(f"Skipping malformed index line: {line.strip()}")
    return entries


def chunker(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


def speciesFromSrc(src):
    if "." in src:
        return src.split(".", 1)[0]
    return src


def parseBlock(block_text):
    species_stats = defaultdict(lambda: {"copies": 0, "gaps": 0, "nongaps": 0})
    all_gap_species = set()
    present_species = set()
    n_seq_lines = 0
    parse_errors = 0
    seqs = []

    for raw_line in block_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        fields = line.split()
        if not fields or fields[0] != "s":
            continue
        if len(fields) < 7:
            parse_errors += 1
            continue
        sp = speciesFromSrc(fields[1])
        seq = fields[6]
        gaps = seq.count("-")
        nongaps = len(seq) - gaps

        n_seq_lines += 1
        present_species.add(sp)
        seqs.append(seq)
        rec = species_stats[sp]
        rec["copies"] += 1
        rec["gaps"] += gaps
        rec["nongaps"] += nongaps

    for sp, rec in species_stats.items():
        if rec["nongaps"] == 0:
            all_gap_species.add(sp)

    total_gaps = sum(rec["gaps"] for rec in species_stats.values())
    total_cells = sum((rec["gaps"] + rec["nongaps"]) for rec in species_stats.values())
    unique_species = len(species_stats)
    duplicate_copies = sum(max(0, rec["copies"] - 1) for rec in species_stats.values())
    site_gap_hist = [0] * (SITE_GAP_HIST_BINS + 1)
    all_gap_sites = 0
    invariant_sites = 0
    variable_sites = 0
    parsimony_informative_sites = 0

    if seqs:
        for col in zip(*seqs):
            gap_count = col.count("-")
            gap_fraction = gap_count / len(col)
            gap_bin = min(SITE_GAP_HIST_BINS, int(round(gap_fraction * SITE_GAP_HIST_BINS)))
            site_gap_hist[gap_bin] += 1

            nongap_states = [c.upper() for c in col if c != "-"]
            if not nongap_states:
                all_gap_sites += 1
                continue

            state_counts = defaultdict(int)
            for state in nongap_states:
                state_counts[state] += 1

            if len(state_counts) == 1:
                invariant_sites += 1
            else:
                variable_sites += 1
                informative_states = sum(1 for count in state_counts.values() if count >= 2)
                if informative_states >= 2:
                    parsimony_informative_sites += 1

    return {
        "species_stats": species_stats,
        "present_species": present_species,
        "all_gap_species": all_gap_species,
        "n_seq_lines": n_seq_lines,
        "unique_species": unique_species,
        "duplicate_copies": duplicate_copies,
        "total_gaps": total_gaps,
        "total_cells": total_cells,
        "gap_fraction_cells": (float(total_gaps) / total_cells if total_cells else 0.0),
        "parse_errors": parse_errors,
        "site_gap_hist": site_gap_hist,
        "all_gap_sites": all_gap_sites,
        "invariant_sites": invariant_sites,
        "variable_sites": variable_sites,
        "parsimony_informative_sites": parsimony_informative_sites,
    }


def workerTask(
    task_id,
    maf_file,
    maf_compression,
    entries,
    write_block_rows,
    tmp_dir,
):
    opener = gzip.open if maf_compression == "gz" else open

    overall = {
        "total_blocks": 0,
        "total_alignment_columns": 0,
        "total_seq_lines": 0,
        "sum_unique_species_per_block": 0,
        "sum_duplicate_copies_per_block": 0,
        "sum_gap_fraction_per_block": 0.0,
        "total_gap_chars": 0,
        "total_cells": 0,
        "total_parse_errors": 0,
        "blocks_with_zero_parsed_seq_lines_but_index_nonzero": 0,
        "all_gap_sites": 0,
        "invariant_sites": 0,
        "variable_sites": 0,
        "parsimony_informative_sites": 0,
        "site_gap_hist": [0] * (SITE_GAP_HIST_BINS + 1),
    }
    species = defaultdict(
        lambda: {
            "blocks_present": 0,
            "copy_lines": 0,
            "duplicated_blocks": 0,
            "duplicate_copies_total": 0,
            "gaps_total": 0,
            "nongaps_total": 0,
            "gap_pct_block_sum": 0.0,
            "all_gap_blocks": 0,
        }
    )
    observed_species = set()

    block_tmp = None
    block_tmp_path = None
    if write_block_rows:
        block_tmp_path = os.path.join(tmp_dir, f"maf_stats.blocks.task{task_id}.tsv")
        block_tmp = open(block_tmp_path, "w", encoding="utf-8")

    with opener(maf_file, "rb") as maf_fp:
        for entry in entries:
            (
                block_id,
                scaffold,
                ref_start,
                ref_length,
                aln_length,
                n_seq_lines_idx,
                offset_start,
                offset_end,
            ) = entry

            maf_fp.seek(offset_start)
            block_bytes = maf_fp.read(offset_end - offset_start)
            block_text = block_bytes.decode("utf-8", errors="replace")
            parsed = parseBlock(block_text)

            overall["total_blocks"] += 1
            overall["total_alignment_columns"] += aln_length
            overall["total_seq_lines"] += parsed["n_seq_lines"]
            overall["sum_unique_species_per_block"] += parsed["unique_species"]
            overall["sum_duplicate_copies_per_block"] += parsed["duplicate_copies"]
            overall["sum_gap_fraction_per_block"] += parsed["gap_fraction_cells"]
            overall["total_gap_chars"] += parsed["total_gaps"]
            overall["total_cells"] += parsed["total_cells"]
            overall["total_parse_errors"] += parsed["parse_errors"]
            overall["all_gap_sites"] += parsed["all_gap_sites"]
            overall["invariant_sites"] += parsed["invariant_sites"]
            overall["variable_sites"] += parsed["variable_sites"]
            overall["parsimony_informative_sites"] += parsed["parsimony_informative_sites"]
            for i, count in enumerate(parsed["site_gap_hist"]):
                overall["site_gap_hist"][i] += count
            if n_seq_lines_idx > 0 and parsed["n_seq_lines"] == 0:
                overall["blocks_with_zero_parsed_seq_lines_but_index_nonzero"] += 1

            observed_species.update(parsed["present_species"])

            for sp, rec in parsed["species_stats"].items():
                s = species[sp]
                s["blocks_present"] += 1
                s["copy_lines"] += rec["copies"]
                if rec["copies"] > 1:
                    s["duplicated_blocks"] += 1
                    s["duplicate_copies_total"] += rec["copies"] - 1
                s["gaps_total"] += rec["gaps"]
                s["nongaps_total"] += rec["nongaps"]
                denom = rec["gaps"] + rec["nongaps"]
                s["gap_pct_block_sum"] += (float(rec["gaps"]) / denom) if denom else 0.0
                if rec["nongaps"] == 0:
                    s["all_gap_blocks"] += 1

            if block_tmp is not None:
                block_tmp.write(
                    "\t".join(
                        [
                            str(block_id),
                            scaffold,
                            str(ref_start),
                            str(ref_length),
                            str(aln_length),
                            str(n_seq_lines_idx),
                            str(parsed["n_seq_lines"]),
                            str(parsed["unique_species"]),
                            str(parsed["duplicate_copies"]),
                            str(parsed["total_gaps"]),
                            f"{parsed['gap_fraction_cells']:.8f}",
                            str(len(parsed["all_gap_species"])),
                            ",".join(sorted(parsed["all_gap_species"])) if parsed["all_gap_species"] else ".",
                            ",".join(sorted(parsed["present_species"])) if parsed["present_species"] else ".",
                        ]
                    )
                    + "\n"
                )

    if block_tmp is not None:
        block_tmp.close()

    return {
        "overall": overall,
        "species": dict(species),
        "observed_species": sorted(observed_species),
        "block_tmp_path": block_tmp_path,
    }


def mergeStats(results):
    overall = defaultdict(float)
    species = defaultdict(
        lambda: {
            "blocks_present": 0,
            "copy_lines": 0,
            "duplicated_blocks": 0,
            "duplicate_copies_total": 0,
            "gaps_total": 0,
            "nongaps_total": 0,
            "gap_pct_block_sum": 0.0,
            "all_gap_blocks": 0,
        }
    )
    observed_species = set()
    block_tmp_paths = []

    for res in results:
        for k, v in res["overall"].items():
            if k == "site_gap_hist":
                if "site_gap_hist" not in overall:
                    overall["site_gap_hist"] = [0] * (SITE_GAP_HIST_BINS + 1)
                for i, count in enumerate(v):
                    overall["site_gap_hist"][i] += count
            else:
                overall[k] += v
        for sp, rec in res["species"].items():
            s = species[sp]
            for k, v in rec.items():
                s[k] += v
        observed_species.update(res["observed_species"])
        if res["block_tmp_path"]:
            block_tmp_paths.append(res["block_tmp_path"])

    return overall, species, observed_species, block_tmp_paths


def writeOverall(out_path, overall, n_species_total):
    total_blocks = int(overall["total_blocks"])
    total_cells = int(overall["total_cells"])

    metrics = [
        ("total_blocks", total_blocks),
        ("total_species", n_species_total),
        ("total_alignment_columns", int(overall["total_alignment_columns"])),
        ("total_sequence_lines", int(overall["total_seq_lines"])),
        ("all_gap_sites", int(overall["all_gap_sites"])),
        ("invariant_sites", int(overall["invariant_sites"])),
        ("variable_sites", int(overall["variable_sites"])),
        ("parsimony_informative_sites", int(overall["parsimony_informative_sites"])),
        ("avg_species_per_block", (overall["sum_unique_species_per_block"] / total_blocks) if total_blocks else 0.0),
        ("avg_duplicate_species_per_block", (overall["sum_duplicate_copies_per_block"] / total_blocks) if total_blocks else 0.0),
        ("avg_gap_fraction_per_block", (overall["sum_gap_fraction_per_block"] / total_blocks) if total_blocks else 0.0),
        ("global_gap_fraction_cells", (overall["total_gap_chars"] / total_cells) if total_cells else 0.0),
        (
            "avg_species_coverage_per_block",
            ((overall["sum_unique_species_per_block"] / total_blocks) / n_species_total) if total_blocks and n_species_total else 0.0,
        ),
        ("total_parse_errors", int(overall["total_parse_errors"])),
        (
            "blocks_with_zero_parsed_seq_lines_but_index_nonzero",
            int(overall["blocks_with_zero_parsed_seq_lines_but_index_nonzero"]),
        ),
    ]

    with open(out_path, "w", encoding="utf-8") as fp:
        fp.write("metric\tvalue\n")
        for key, val in metrics:
            if isinstance(val, float):
                fp.write(f"{key}\t{val:.8f}\n")
            else:
                fp.write(f"{key}\t{val}\n")


def writeSpecies(out_path, species, total_blocks):
    header = [
        "species",
        "blocks_present",
        "blocks_missing",
        "presence_fraction",
        "copy_lines",
        "duplicated_blocks",
        "duplicate_copies_total",
        "gaps_total",
        "nongaps_total",
        "gap_fraction_cells",
        "avg_gaps_per_present_block",
        "avg_gap_fraction_per_present_block",
        "all_gap_blocks",
    ]
    with open(out_path, "w", encoding="utf-8") as fp:
        fp.write("\t".join(header) + "\n")
        for sp in sorted(species):
            rec = species[sp]
            blocks_present = int(rec["blocks_present"])
            blocks_missing = max(0, total_blocks - blocks_present)
            total_cells = rec["gaps_total"] + rec["nongaps_total"]
            fp.write(
                "\t".join(
                    [
                        sp,
                        str(blocks_present),
                        str(blocks_missing),
                        f"{(blocks_present / total_blocks) if total_blocks else 0.0:.8f}",
                        str(int(rec["copy_lines"])),
                        str(int(rec["duplicated_blocks"])),
                        str(int(rec["duplicate_copies_total"])),
                        str(int(rec["gaps_total"])),
                        str(int(rec["nongaps_total"])),
                        f"{(rec['gaps_total'] / total_cells) if total_cells else 0.0:.8f}",
                        f"{(rec['gaps_total'] / blocks_present) if blocks_present else 0.0:.8f}",
                        f"{(rec['gap_pct_block_sum'] / blocks_present) if blocks_present else 0.0:.8f}",
                        str(int(rec["all_gap_blocks"])),
                    ]
                )
                + "\n"
            )


def writeBlockTable(out_path, block_tmp_paths, n_species_total, expected_species):
    expected_species = set(expected_species) if expected_species else None
    header = [
        "block_id",
        "scaffold",
        "ref_start",
        "ref_length",
        "aln_length",
        "n_seq_lines_index",
        "n_seq_lines_parsed",
        "n_unique_species",
        "n_duplicate_copies",
        "total_gaps",
        "gap_fraction_cells",
        "n_all_gap_species",
        "all_gap_species",
        "missing_species_count",
        "all_gap_or_missing_count",
        "missing_species",
        "all_gap_or_missing_species",
        "present_species",
    ]
    with open(out_path, "w", encoding="utf-8") as out_fp:
        out_fp.write("\t".join(header) + "\n")

        for tmp_path in sorted(block_tmp_paths):
            with open(tmp_path, "r", encoding="utf-8") as in_fp:
                for line in in_fp:
                    fields = line.rstrip("\n").split("\t")
                    n_unique_species = int(fields[7])
                    all_gap_species = set() if fields[12] == "." else set(fields[12].split(","))
                    present_species = set() if fields[13] == "." else set(fields[13].split(","))

                    if expected_species is None:
                        missing_species_count = max(0, n_species_total - n_unique_species)
                        missing_species_str = "."
                        all_gap_or_missing_species = all_gap_species
                    else:
                        missing_species = expected_species - present_species
                        missing_species_count = len(missing_species)
                        missing_species_str = ",".join(sorted(missing_species)) if missing_species else "."
                        all_gap_or_missing_species = all_gap_species.union(missing_species)

                    all_gap_or_missing_count = missing_species_count + len(all_gap_species)
                    all_gap_or_missing_species_str = (
                        ",".join(sorted(all_gap_or_missing_species))
                        if all_gap_or_missing_species
                        else "."
                    )

                    out_fp.write(
                        "\t".join(
                            [
                                fields[0],  # block_id
                                fields[1],  # scaffold
                                fields[2],  # ref_start
                                fields[3],  # ref_length
                                fields[4],  # aln_length
                                fields[5],  # n_seq_lines_index
                                fields[6],  # n_seq_lines_parsed
                                fields[7],  # n_unique_species
                                fields[8],  # n_duplicate_copies
                                fields[9],  # total_gaps
                                fields[10],  # gap_fraction_cells
                                fields[11],  # n_all_gap_species
                                fields[12],  # all_gap_species
                                str(missing_species_count),
                                str(all_gap_or_missing_count),
                                missing_species_str,
                                all_gap_or_missing_species_str,
                                fields[13],  # present_species
                            ]
                        )
                        + "\n"
                    )

            os.remove(tmp_path)


#############################################################################


def _figureToDataUri(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=105, bbox_inches="tight")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def sampleBlockRows(block_tmp_paths, max_points):
    max_points = max(1, int(max_points))
    sample = []
    seen = 0
    rng = random.Random(1)

    for tmp_path in sorted(block_tmp_paths):
        with open(tmp_path, "r", encoding="utf-8") as fp:
            for line in fp:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue
                try:
                    row = {
                        "aln_length": int(fields[4]),
                        "n_unique_species": int(fields[7]),
                        "n_duplicate_copies": int(fields[8]),
                        "gap_fraction_cells": float(fields[10]),
                    }
                except ValueError:
                    continue

                seen += 1
                if len(sample) < max_points:
                    sample.append(row)
                else:
                    j = rng.randint(1, seen)
                    if j <= max_points:
                        sample[j - 1] = row

    return sample, seen


def writeDashboard(
    out_path,
    overall,
    species,
    species_total,
    top_n,
    max_block_points,
    block_tmp_paths,
    cmdline,
    generated_at,
    LOG,
):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        LOG.warning("Skipping dashboard: matplotlib unavailable (%s)", exc)
        return

    total_blocks = int(overall["total_blocks"])
    total_cells = int(overall["total_cells"])
    total_sites = int(overall["total_alignment_columns"])
    avg_species_per_block = (overall["sum_unique_species_per_block"] / total_blocks) if total_blocks else 0.0
    avg_gap_fraction_per_block = (overall["sum_gap_fraction_per_block"] / total_blocks) if total_blocks else 0.0
    global_gap_fraction_cells = (overall["total_gap_chars"] / total_cells) if total_cells else 0.0
    site_gap_hist = list(overall.get("site_gap_hist", [0] * (SITE_GAP_HIST_BINS + 1)))
    invariant_sites = int(overall["invariant_sites"])
    variable_sites = int(overall["variable_sites"])
    parsimony_informative_sites = int(overall["parsimony_informative_sites"])
    all_gap_sites = int(overall["all_gap_sites"])

    rows = []
    for sp, rec in species.items():
        blocks_present = int(rec["blocks_present"])
        total_sp_cells = rec["gaps_total"] + rec["nongaps_total"]
        rows.append(
            {
                "species": sp,
                "presence_fraction": (blocks_present / total_blocks) if total_blocks else 0.0,
                "gap_fraction_cells": (rec["gaps_total"] / total_sp_cells) if total_sp_cells else 0.0,
                "duplicated_blocks": int(rec["duplicated_blocks"]),
                "blocks_present": blocks_present,
            }
        )

    top_n = max(1, int(top_n))
    presence_low = sorted(rows, key=lambda r: (r["presence_fraction"], r["species"]))[:top_n]
    presence_high = sorted(rows, key=lambda r: (r["presence_fraction"], r["species"]), reverse=True)[:top_n]
    gap_sorted = sorted(rows, key=lambda r: (r["gap_fraction_cells"], r["species"]), reverse=True)[:top_n]
    dup_sorted = sorted(rows, key=lambda r: (r["duplicated_blocks"], r["species"]), reverse=True)[:top_n]
    presence_zero = sum(1 for r in rows if r["blocks_present"] == 0)
    presence_full = sum(1 for r in rows if total_blocks and r["blocks_present"] == total_blocks)
    duplicated_any = sum(1 for r in rows if r["duplicated_blocks"] > 0)
    gap_free = sum(1 for r in rows if r["gap_fraction_cells"] == 0.0)

    species_summary_items = [
        f"Lowest presence species: {presence_low[0]['species']} ({presence_low[0]['presence_fraction']:.4f})" if presence_low else "Lowest presence species: n/a",
        f"Highest presence species: {presence_high[0]['species']} ({presence_high[0]['presence_fraction']:.4f})" if presence_high else "Highest presence species: n/a",
        f"Highest gap species: {gap_sorted[0]['species']} ({gap_sorted[0]['gap_fraction_cells']:.4f})" if gap_sorted else "Highest gap species: n/a",
        f"Most duplicated species: {dup_sorted[0]['species']} ({dup_sorted[0]['duplicated_blocks']:,} duplicated blocks)" if dup_sorted else "Most duplicated species: n/a",
        f"Species present in every block: {presence_full:,} of {species_total:,}",
        f"Species absent from all blocks: {presence_zero:,} of {species_total:,}",
        f"Species with zero cell-level gaps: {gap_free:,} of {species_total:,}",
        f"Species with any duplicated block: {duplicated_any:,} of {species_total:,}",
    ]
    species_table_rows = sorted(
        rows,
        key=lambda r: (-r["presence_fraction"], r["gap_fraction_cells"], r["species"]),
    )

    def makeBarFigure(data_rows, x_key, title, color):
        fig_height = max(2.5, min(5.0, 1.2 + 0.16 * len(data_rows)))
        fig, ax = plt.subplots(figsize=(8.4, fig_height))
        labels = [r["species"] for r in data_rows]
        values = [r[x_key] for r in data_rows]
        y = list(range(len(labels)))
        ax.barh(y, values, color=color)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_title(title)
        ax.grid(axis="x", linestyle="--", alpha=0.35)
        return fig

    fig_presence_low = makeBarFigure(
        presence_low,
        "presence_fraction",
        f"Lowest Species Presence Fraction (bottom {len(presence_low)})",
        "#1f77b4",
    )
    fig_presence_high = makeBarFigure(
        presence_high,
        "presence_fraction",
        f"Highest Species Presence Fraction (top {len(presence_high)})",
        "#6a3d9a",
    )
    fig_gap = makeBarFigure(
        gap_sorted,
        "gap_fraction_cells",
        f"Highest Species Gap Fraction (top {len(gap_sorted)})",
        "#d62728",
    )
    fig_dup = makeBarFigure(
        dup_sorted,
        "duplicated_blocks",
        f"Most Duplicated Species (top {len(dup_sorted)})",
        "#2ca02c",
    )

    fig_sp_scatter, ax_sp_scatter = plt.subplots(figsize=(6.3, 4.0))
    ax_sp_scatter.scatter(
        [r["presence_fraction"] for r in rows],
        [r["gap_fraction_cells"] for r in rows],
        s=16,
        alpha=0.7,
        color="#1f78b4",
        edgecolors="none",
    )
    ax_sp_scatter.set_title("Species: Presence Fraction vs Gap Fraction")
    ax_sp_scatter.set_xlabel("Presence fraction")
    ax_sp_scatter.set_ylabel("Gap fraction (cells)")
    ax_sp_scatter.grid(linestyle="--", alpha=0.3)

    fig_site_gap, ax_site_gap = plt.subplots(figsize=(7.8, 4.6))
    bin_labels = [i / SITE_GAP_HIST_BINS for i in range(SITE_GAP_HIST_BINS + 1)]
    ax_site_gap.bar(bin_labels, site_gap_hist, width=(1.0 / SITE_GAP_HIST_BINS) * 0.9, color="#dd8452")
    ax_site_gap.set_title("Sites: Gap Fraction Distribution")
    ax_site_gap.set_xlabel("Gap fraction")
    ax_site_gap.set_ylabel("Count")
    ax_site_gap.grid(axis="y", linestyle="--", alpha=0.3)

    block_rows_sample, block_rows_seen = sampleBlockRows(block_tmp_paths, max_block_points)
    block_plot_html = []
    block_plot_notice = ""
    if block_rows_sample:
        fig_block_gap, ax_block_gap = plt.subplots(figsize=(6.3, 4.0))
        ax_block_gap.hist([r["gap_fraction_cells"] for r in block_rows_sample], bins=35, color="#e6550d")
        ax_block_gap.set_title("Blocks: Gap Fraction Distribution")
        ax_block_gap.set_xlabel("Gap fraction (cells)")
        ax_block_gap.set_ylabel("Count")
        ax_block_gap.grid(linestyle="--", alpha=0.3)
        block_plot_html.append(
            f'<div class="plot"><img alt="Blocks gap fraction distribution" src="{_figureToDataUri(fig_block_gap)}" /></div>'
        )
        plt.close(fig_block_gap)

        fig_block_species, ax_block_species = plt.subplots(figsize=(6.3, 4.0))
        ax_block_species.hist([r["n_unique_species"] for r in block_rows_sample], bins=35, color="#3182bd")
        ax_block_species.set_title("Blocks: Unique Species per Block")
        ax_block_species.set_xlabel("Unique species")
        ax_block_species.set_ylabel("Count")
        ax_block_species.grid(linestyle="--", alpha=0.3)
        block_plot_html.append(
            f'<div class="plot"><img alt="Blocks unique species distribution" src="{_figureToDataUri(fig_block_species)}" /></div>'
        )
        plt.close(fig_block_species)

        fig_block_scatter, ax_block_scatter = plt.subplots(figsize=(6.3, 4.0))
        ax_block_scatter.scatter(
            [r["aln_length"] for r in block_rows_sample],
            [r["gap_fraction_cells"] for r in block_rows_sample],
            s=10,
            alpha=0.45,
            color="#756bb1",
            edgecolors="none",
        )
        ax_block_scatter.set_title("Blocks: Alignment Length vs Gap Fraction")
        ax_block_scatter.set_xlabel("Alignment length")
        ax_block_scatter.set_ylabel("Gap fraction (cells)")
        ax_block_scatter.grid(linestyle="--", alpha=0.3)
        block_plot_html.append(
            f'<div class="plot"><img alt="Blocks alignment length versus gap fraction" src="{_figureToDataUri(fig_block_scatter)}" /></div>'
        )
        plt.close(fig_block_scatter)
    elif not block_tmp_paths:
        block_plot_notice = (
            "Block-level plots are unavailable because no block-row data were available to the dashboard. "
            "This usually means block-level rows were not produced before dashboard generation."
        )
    else:
        block_plot_notice = (
            "Block-level plots are unavailable because block-row files were present, but no usable block rows "
            "were parsed for plotting."
        )

    presence_low_img = _figureToDataUri(fig_presence_low)
    presence_high_img = _figureToDataUri(fig_presence_high)
    gap_img = _figureToDataUri(fig_gap)
    dup_img = _figureToDataUri(fig_dup)
    sp_scatter_img = _figureToDataUri(fig_sp_scatter)
    site_gap_img = _figureToDataUri(fig_site_gap)
    plt.close(fig_presence_low)
    plt.close(fig_presence_high)
    plt.close(fig_gap)
    plt.close(fig_dup)
    plt.close(fig_sp_scatter)
    plt.close(fig_site_gap)

    cards = [
        ("Total blocks", f"{total_blocks:,}"),
        ("Total species", f"{species_total:,}"),
        ("Total alignment columns", f"{int(overall['total_alignment_columns']):,}"),
        ("Total sequence lines", f"{int(overall['total_seq_lines']):,}"),
        ("Avg species per block", f"{avg_species_per_block:.4f}"),
        ("Avg gap fraction per block", f"{avg_gap_fraction_per_block:.4f}"),
        ("Global gap fraction (cells)", f"{global_gap_fraction_cells:.4f}"),
        ("Invariant sites", f"{invariant_sites:,}"),
        ("Variable sites", f"{variable_sites:,}"),
        ("Parsimony-informative sites", f"{parsimony_informative_sites:,}"),
        ("Total parse errors", f"{int(overall['total_parse_errors']):,}"),
    ]

    card_html = "\n".join(
        [
            f"<div class='card'><div class='k'>{html.escape(k)}</div><div class='v'>{html.escape(v)}</div></div>"
            for k, v in cards
        ]
    )

    species_summary_html = "".join(f"<li>{html.escape(item)}</li>" for item in species_summary_items)
    species_table_html = "\n".join(
        [
            (
                "<tr>"
                f"<td>{html.escape(r['species'])}</td>"
                f"<td>{r['blocks_present']:,}</td>"
                f"<td>{total_blocks:,}</td>"
                f"<td>{r['presence_fraction']:.4f}</td>"
                f"<td>{r['gap_fraction_cells']:.4f}</td>"
                f"<td>{r['duplicated_blocks']:,}</td>"
                "</tr>"
            )
            for r in species_table_rows
        ]
    )
    plot_sections = [
        (
            "plot-presence-low",
            "Lowest species presence fraction",
            "Species at the low end of this plot are missing from many blocks. This is useful for spotting taxa with sparse representation across the alignment set.",
            presence_low_img,
            "Lowest species presence fraction plot",
        ),
        (
            "plot-presence-high",
            "Highest species presence fraction",
            "Species at the high end of this plot are present in nearly all blocks. This gives a quick sense of the best-covered taxa in the dataset.",
            presence_high_img,
            "Highest species presence fraction plot",
        ),
        (
            "plot-gap-high",
            "Highest species gap fraction",
            "These species contribute the largest share of gap cells across all parsed alignment cells. High values can indicate fragmented or low-coverage representation.",
            gap_img,
            "Highest species gap fraction plot",
        ),
        (
            "plot-dup",
            "Most duplicated species",
            "This highlights species that appear more than once within the same block most often, which can indicate duplicated alignments or paralog-like structure.",
            dup_img,
            "Most duplicated species plot",
        ),
        (
            "plot-species-scatter",
            "Species presence versus gap fraction",
            "This compares two species-level summaries at once: how often a species appears and how gappy it is when it does appear.",
            sp_scatter_img,
            "Species presence versus gap fraction scatter plot",
        ),
        (
            "plot-site-gap",
            "Site gap-fraction distribution",
            "This shows the distribution of per-column gappiness across all alignment columns, helping distinguish mostly complete columns from heavily gapped columns.",
            site_gap_img,
            "Site gap fraction distribution",
        ),
    ]

    block_section_html = ""
    if block_rows_sample:
        block_section_html = "".join(
            [
                f"""
  <div class="section" id="plot-block-gap">
    <b>Block gap-fraction distribution</b>
    <p>This summarizes how gappy whole alignment blocks are. It helps identify whether missing data is concentrated in a few blocks or spread broadly.</p>
    {block_plot_html[0]}
  </div>
  <div class="section" id="plot-block-species">
    <b>Unique species per block</b>
    <p>This shows how many distinct species are represented in each sampled block and gives a quick view of block-level occupancy.</p>
    {block_plot_html[1]}
  </div>
  <div class="section" id="plot-block-scatter">
    <b>Alignment length versus block gap fraction</b>
    <p>This compares block length to block-level gappiness, which can help reveal whether shorter or longer blocks tend to carry more missing data.</p>
    {block_plot_html[2]}
  </div>
"""
            ]
        )
    elif block_plot_notice:
        block_section_html = f"""
  <div class="section" id="plot-block-gap">
    <b>Block-level plots unavailable</b>
    <p>{html.escape(block_plot_notice)}</p>
  </div>
"""
    plot_section_html = "".join(
        [
            f"""
  <div class="section" id="{section_id}">
    <b>{html.escape(title)}</b>
    <p>{html.escape(desc)}</p>
    <div class="plot"><img alt="{html.escape(alt)}" src="{img}" /></div>
  </div>
"""
            for section_id, title, desc, img, alt in plot_sections
        ]
    )

    dashboard_html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>MAF Stats Dashboard</title>
<style>
body {{
  margin: 0;
  font-family: "Segoe UI", Arial, sans-serif;
  color: #14202b;
  background: #f5f7fa;
}}
.page {{
  display: grid;
  grid-template-columns: 260px minmax(0, 1fr);
  gap: 24px;
  padding: 18px;
}}
h1 {{
  margin: 0 0 12px 0;
  font-size: 24px;
}}
.header-panel,
.section,
.defs,
.plot {{
  background: #ffffff;
  border: 1px solid #d4dde8;
  border-radius: 8px;
}}
.header-panel,
.section,
.defs {{
  padding: 12px 14px;
  margin-bottom: 16px;
}}
.toc {{
  position: sticky;
  top: 18px;
  align-self: start;
  background: #ffffff;
  border: 1px solid #d4dde8;
  border-radius: 8px;
  padding: 14px;
}}
.toc ul {{
  list-style: none;
  margin: 10px 0 0 0;
  padding: 0;
}}
.toc li {{
  margin: 0 0 8px 0;
}}
.toc a {{
  color: #284760;
  text-decoration: none;
}}
.toc a:hover {{
  text-decoration: underline;
}}
.header-panel pre {{
  margin: 8px 0 0 0;
  white-space: pre-wrap;
  word-break: break-word;
  font-size: 12px;
  color: #31485f;
}}
.cards {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 10px;
  margin-bottom: 16px;
}}
.card {{
  background: #ffffff;
  border: 1px solid #d4dde8;
  border-radius: 8px;
  padding: 10px 12px;
}}
.k {{ font-size: 12px; color: #47617a; }}
.v {{ font-size: 20px; font-weight: 700; margin-top: 4px; }}
.plot {{
  padding: 8px;
  width: min(100%, 860px);
  margin: 12px auto 0 auto;
}}
.plot img {{
  width: 100%;
  height: auto;
  display: block;
  margin: 0 auto;
}}
.section p {{
  max-width: 78ch;
}}
.meta {{
  font-size: 12px;
  color: #4b5f73;
  margin-top: 10px;
}}
.defs,
.section {{
  font-size: 13px;
}}
.defs ul,
.section ul {{
  margin: 8px 0 0 18px;
  padding: 0;
}}
.defs li,
.section li {{
  margin: 0 0 6px 0;
}}
table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}}
th, td {{
  border-bottom: 1px solid #dfe6ee;
  padding: 6px 8px;
  text-align: left;
}}
th {{
  background: #eef3f8;
  position: sticky;
  top: 0;
}}
.table-wrap {{
  max-height: 320px;
  overflow: auto;
  border: 1px solid #dfe6ee;
  border-radius: 6px;
  margin-top: 10px;
}}
@media (max-width: 1100px) {{
  .page {{
    grid-template-columns: 1fr;
  }}
  .toc {{
    position: static;
  }}
}}
</style>
</head>
<body>
  <div class="page">
  <div class="toc">
    <b>Contents</b>
    <ul>
      <li><a href="#top">Top</a></li>
      <li><a href="#overview">Overview metrics</a></li>
      <li><a href="#species-summary">Species summary</a></li>
      <li><a href="#species-table">Per-species summary</a></li>
      <li><a href="#plot-presence-low">Lowest presence plot</a></li>
      <li><a href="#plot-presence-high">Highest presence plot</a></li>
      <li><a href="#plot-gap-high">Highest gap plot</a></li>
      <li><a href="#plot-dup">Duplication plot</a></li>
      <li><a href="#plot-species-scatter">Species scatter</a></li>
      <li><a href="#plot-site-gap">Site gap plot</a></li>
      <li><a href="#plot-block-gap">Block plots</a></li>
      <li><a href="#notes">Metric notes</a></li>
    </ul>
  </div>
  <div class="main" id="top">
  <h1>MAF Summary Dashboard</h1>
  <div class="header-panel">
    <b>Generated</b>: {html.escape(generated_at)}<br />
    <b>Dashboard file</b>: {html.escape(out_path)}<br />
    <b>Command</b>:
    <pre>{html.escape(cmdline)}</pre>
  </div>
  <div class="cards" id="overview">
    {card_html}
  </div>
  <div class="section" id="species-summary">
    <b>Species summary</b>
    <ul>
      {species_summary_html}
    </ul>
  </div>
  <div class="section" id="species-table">
    <b>Per-species summary</b>
    <div class="table-wrap">
      <table>
        <thead>
          <tr>
            <th>Species</th>
            <th>Blocks present</th>
            <th>Total blocks</th>
            <th>Presence fraction</th>
            <th>Gap fraction (cells)</th>
            <th>Duplicated blocks</th>
          </tr>
        </thead>
        <tbody>
          {species_table_html}
        </tbody>
      </table>
    </div>
  </div>
  {plot_section_html}
  {block_section_html}
  <div class="meta">
    Generated by maf_stats.py from in-memory summary statistics.
    Block plots use reservoir sampling: {len(block_rows_sample):,} of {block_rows_seen:,} blocks.
  </div>
  <div class="defs" id="notes">
    <b>Metric notes</b>
    <ul>
      <li>Global gap fraction (cells) = total gap characters across all sequence lines and blocks, divided by total alignment cells (gaps + non-gaps) across all sequence lines and blocks.</li>
      <li>Avg species per block = average number of unique species observed per alignment block.</li>
      <li>Avg gap fraction per block = mean of each block's cell-level gap fraction, not the same as the global pooled gap fraction.</li>
      <li>Species presence fraction = fraction of blocks in which a species appears at least once.</li>
      <li>Species gap fraction (cells) = a species' total gap characters divided by its total gap + non-gap cells across all parsed blocks.</li>
      <li>Duplicated blocks = number of blocks where a species appears more than once.</li>
      <li>Invariant sites have one observed non-gap state across the column; variable sites have more than one observed non-gap state.</li>
      <li>Parsimony-informative sites have at least two non-gap states, each observed at least twice.</li>
      <li>All-gap sites have no observed non-gap state in the column.</li>
      <li>Total parse errors counts blocks that could not be parsed cleanly and were skipped or partially counted.</li>
      <li>Block plots are based on reservoir sampling of block rows when the full block table is large.</li>
    </ul>
  </div>
  </div>
  </div>
</body>
</html>
"""

    with open(out_path, "w", encoding="utf-8") as fp:
        fp.write(dashboard_html)


#############################################################################


def main():
    args = optParse()
    cmdline = " ".join(sys.argv)
    generated_at = datetime.datetime.now().astimezone().isoformat(timespec="seconds")

    out_prefix = args.output_prefix
    out_dir = os.path.dirname(out_prefix) if os.path.dirname(out_prefix) else "."
    os.makedirs(out_dir, exist_ok=True)

    log_file = out_prefix + ".log"
    logger_name = LOGINIT.configureLogging(
        log_level=args.log_level,
        log_verbosity="BOTH",
        log_filename=log_file,
        logger_name="maf_stats_logger",
        overwrite_log_file=True,
    )
    LOG = logging.getLogger(logger_name)

    LOG.info(f"maf_stats.py called as: {cmdline}")
    LOG.info("-" * 40)
    for arg, value in vars(args).items():
        LOG.info(f"{arg:22} : {value}")
    LOG.info("-" * 40)

    if args.processes < 1:
        LOG.error("--processes must be >= 1")
        sys.exit(1)
    if args.chunk_size < 1:
        LOG.error("--chunk-size must be >= 1")
        sys.exit(1)

    maf_compression = COMMON.detectCompression(args.maf_file)
    if maf_compression not in ("none", "gz"):
        LOG.error(f"Unsupported MAF compression: {maf_compression}")
        sys.exit(1)

    expected_species = parseExpectedSpecies(args)
    if expected_species:
        LOG.info(f"Loaded {len(expected_species)} expected species for exact missing-species reporting.")

    LOG.info(f"Parsing index file: {args.index_file}")
    entries = parseIndex(args.index_file, LOG)
    if not entries:
        LOG.error("No valid index entries found.")
        sys.exit(1)
    LOG.info(f"Loaded {len(entries)} indexed blocks.")

    chunks = list(chunker(entries, args.chunk_size))
    LOG.info(f"Running {len(chunks)} tasks with chunk size {args.chunk_size} across {args.processes} process(es).")

    write_block_rows = (not args.no_block_table) or args.html_dashboard
    with tempfile.TemporaryDirectory(prefix="maf_stats_tmp_", dir=out_dir) as tmp_dir:
        results = []
        with ProcessPoolExecutor(max_workers=args.processes) as executor:
            futures = []
            for task_id, chunk in enumerate(chunks, start=1):
                futures.append(
                    executor.submit(
                        workerTask,
                        task_id,
                        args.maf_file,
                        maf_compression,
                        chunk,
                        write_block_rows,
                        tmp_dir,
                    )
                )
            for f in futures:
                results.append(f.result())

        overall, species, observed_species, block_tmp_paths = mergeStats(results)

        zero_seq_blocks = int(overall["blocks_with_zero_parsed_seq_lines_but_index_nonzero"])
        if zero_seq_blocks:
            LOG.error(
                "Detected %d blocks where index reports sequence rows but parser found none. "
                "Aborting to avoid invalid summary output.",
                zero_seq_blocks,
            )
            LOG.error(
                "Check that the index matches this exact MAF file and rerun after updating to the fixed parser."
            )
            sys.exit(2)

        if expected_species:
            species_total = len(expected_species)
            # Keep species that are expected but absent in all blocks.
            for sp in expected_species:
                if sp not in species:
                    species[sp] = {
                        "blocks_present": 0,
                        "copy_lines": 0,
                        "duplicated_blocks": 0,
                        "duplicate_copies_total": 0,
                        "gaps_total": 0,
                        "nongaps_total": 0,
                        "gap_pct_block_sum": 0.0,
                        "all_gap_blocks": 0,
                    }
        else:
            species_total = len(observed_species)

        overall_file = out_prefix + ".overall.tsv"
        species_file = out_prefix + ".species.tsv"
        block_file = out_prefix + ".block.tsv"

        LOG.info(f"Writing overall stats: {overall_file}")
        writeOverall(overall_file, overall, species_total)

        LOG.info(f"Writing species stats: {species_file}")
        writeSpecies(species_file, species, int(overall["total_blocks"]))

        if args.html_dashboard:
            dashboard_file = out_prefix + ".dashboard.html"
            LOG.info(f"Writing HTML dashboard: {dashboard_file}")
            writeDashboard(
                dashboard_file,
                overall,
                species,
                species_total,
                args.dashboard_top_species,
                args.dashboard_max_block_points,
                block_tmp_paths,
                cmdline,
                generated_at,
                LOG,
            )

        if not args.no_block_table:
            LOG.info(f"Writing block stats: {block_file}")
            writeBlockTable(block_file, block_tmp_paths, species_total, expected_species)
        else:
            LOG.info("Skipping block table (--no-block-table).")

    LOG.info("Done.")


if __name__ == "__main__":
    main()
