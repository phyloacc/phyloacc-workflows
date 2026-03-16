#############################################################################

# snakemake -j 5 -s <snakefile> --configfile <config file> --dryrun

#############################################################################

import os
import sys
import glob
import logging
import traceback

import lib.common as COMMON
from lib.common import spacedOut as SO

from functools import partial


def gap_pct_other(maf_path: str) -> float:
    total_gaps = 0
    total_positions = 0
    in_block = False
    block_seqs = []

    with open(maf_path, "r") as f:
        for line in f:
            if not line.strip():
                if in_block and block_seqs:
                    ref_seq = block_seqs[0]
                    aln_len = len(ref_seq)
                    nonref = block_seqs[1:]
                    if nonref:
                        total_positions += aln_len * len(nonref)
                        for seq in nonref:
                            total_gaps += seq.count('-')
                in_block = False
                block_seqs = []
                continue

            if line.startswith('a'):
                in_block = True
                block_seqs = []
                continue

            if in_block and line.startswith('s '):
                parts = line.rstrip("\n").split()
                if len(parts) >= 7:
                    block_seqs.append(parts[6])

        # handle file without trailing blank line
        if in_block and block_seqs:
            ref_seq = block_seqs[0]
            aln_len = len(ref_seq)
            nonref = block_seqs[1:]
            if nonref:
                total_positions += aln_len * len(nonref)
                for seq in nonref:
                    total_gaps += seq.count('-')

    if total_positions == 0:
        return 1.0
    return total_gaps / total_positions

#############################################################################
# System setup

config_flag = config.get("display", False);
version_flag = config.get("version", False);
info_flag = config.get("info", False);
debug = config.get("debug", False);
#debug = True;
# A hacky way to get some custom command line arguments for the pipeline
# These just control preprocessing flags that stop the pipeline early anyways

MAIN, DRY_RUN, OUTPUT_DIR, LOG_DIR, TMPDIR, LOG_LEVEL, LOG_VERBOSITY = COMMON.pipelineSetup(config, sys.argv, version_flag, info_flag, config_flag, debug, workflow);
# Setup the pipeline, including the output directory, log directory, and tmp directory

MLOG = logging.getLogger('META')
# Setup logging if debugging

PIPELINE_DIR = os.path.dirname(workflow.snakefile);
# The directory where the snakemake file is located

ENVS_DIR = os.path.join(PIPELINE_DIR, "envs");
# The directory where the conda environments are located

UTILS_DIR = os.path.join(PIPELINE_DIR, "utils");
# The directory where the utility scripts are located

getRuleResources = partial(COMMON.getResources, config)
# This maps the function to get rule resources from the config file
# so we don't have to pass config each time we call it

#############################################################################
# Input files and output paths

PHYLOFIT_DIR = config.get(
    "Phylofit_per_chromosome",
    config.get("phylofit_dir", os.path.join(OUTPUT_DIR, "phylofit"))
)
MIN_KEEP_REGION_LEN = int(config.get("min_keep_region_len", 6))
CONSERVE_DIR = os.path.join(OUTPUT_DIR, "02-conserved-bed")
RHO_STATS_DIR = os.path.join(OUTPUT_DIR, "03-rho_stats")
CNEES_DIR = os.path.join(OUTPUT_DIR, "04-CNEES")
CNEE_MAF_DIR = os.path.join(OUTPUT_DIR, "05-CNEE_fastas")

MAF_SPLIT_BY_CHROM_DIR = config.get("maf_split_by_chromosome", os.path.join(OUTPUT_DIR, "mafs"))
# Directory containing chromosome/scaffold MAFs expected as:
# {MAF_SPLIT_BY_CHROM_DIR}/{chromosome_group}/{ref_chromosome}.maf

MAF_REF_ID = config["maf_ref_id"];
MAF_CHR_PREFIX = config["maf_chr_prefix"];
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"];
MAF_REF_PREFIX = MAF_REF_ID + MAF_REF_CHR_JOINER + MAF_CHR_PREFIX;
# MAF reference and scaffold id parsing

REF_FASTA = config["ref_fasta"];
REF_CHROMOSOME_GROUPS = config["ref_chromosome_groups"];
# Reference genome info

#############################################################################
# Other params

MIN_NS_TO_SPLIT_BY = config["min_Ns_to_split_by"];
NS_INTERVAL_FILE = os.path.join(OUTPUT_DIR, "ns-intervals.txt")

target_chrom_cfg = config.get("target_ref_chromosomes", [])
if isinstance(target_chrom_cfg, str):
    target_chrom_set = set([c.strip() for c in target_chrom_cfg.split(",") if c.strip()])
elif isinstance(target_chrom_cfg, (list, tuple, set)):
    target_chrom_set = set(target_chrom_cfg)
else:
    target_chrom_set = set()

if target_chrom_set:
    flattened_chromosome_groups = [
        (group, chromosome)
        for group, chromosome_list in REF_CHROMOSOME_GROUPS.items()
        for chromosome in chromosome_list
        if chromosome in target_chrom_set
    ]
else:
    flattened_chromosome_groups = [
        (group, chromosome)
        for group, chromosome_list in REF_CHROMOSOME_GROUPS.items()
        for chromosome in chromosome_list
    ]

if not flattened_chromosome_groups:
    raise ValueError("No chromosomes selected. Check ref_chromosome_groups/target_ref_chromosomes in config.")

REF_CHR_GROUPS_LIST, REF_CHROMOSOMES = zip(*flattened_chromosome_groups);
REF_CHR_PATHS = [ REF_CHR_GROUPS_LIST[i] + "/" + REF_CHROMOSOMES[i] for i in range(len(REF_CHROMOSOMES)) ];
# Get two lists of equal length to use for wild cards

####################

if LOG_LEVEL == "debug":
    for i in range(len(REF_CHR_GROUPS_LIST)):
        MLOG.debug(SO(f"CHR GROUP {i}", 15) + f"{REF_CHR_GROUPS_LIST[i]} : {REF_CHROMOSOMES[i]}");
    MLOG.debug("EXITING BEFORE RULES. DEBUG MODE.");
    sys.exit(0);
# Exit before running rules if in debug mode

#############################################################################
localrules: all

MIN_KEEP_REGION_LEN = int(config.get("min_keep_region_len", 6))
MAF_SPLIT_NS_DIR = os.path.join(OUTPUT_DIR, f"01-mafsplit_by_Ns_min{MIN_KEEP_REGION_LEN}")
MAX_GAP_PCT = float(config.get("max_gap_pct", 0.9))

def _as_bool(x, default=False):
    if x is None:
        return default
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)
    s = str(x).strip().lower()
    if s in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "f", "no", "n", "off", ""}:
        return False
    raise ValueError(f"Cannot interpret boolean value from '{x}'")

CLEANUP_CHUNK_INTERMEDIATES = _as_bool(config.get("cleanup_chunk_intermediates", True), True)
RHO_ONLY = _as_bool(config.get("rho_only", False), False)

# Rho mode:
# - fixed (default): do NOT estimate chunk rho; use fixed_rho (default 0.3)
# - estimate: estimate chunk rho and derive chromosome-wide rho via global_rho_stat
RHO_MODE = str(config.get("rho_mode", "fixed")).strip().lower()
estimate_rho_override = config.get("estimate_rho", None)
if estimate_rho_override is not None:
    estimate_rho_override = _as_bool(estimate_rho_override, False)
if estimate_rho_override is True:
    RHO_MODE = "estimate"
if estimate_rho_override is False and config.get("estimate_rho", None) is not None:
    RHO_MODE = "fixed"
if RHO_MODE not in {"fixed", "estimate"}:
    raise ValueError(f"Invalid rho_mode '{RHO_MODE}'. Use 'fixed' or 'estimate'.")

FIXED_RHO = float(config.get("fixed_rho", config.get("rho", 0.3)))
GLOBAL_RHO_STAT = str(config.get("global_rho_stat", "p90")).strip().lower()
if GLOBAL_RHO_STAT in {"90pct", "q90", "90", "percentile90"}:
    GLOBAL_RHO_STAT = "p90"
if GLOBAL_RHO_STAT not in {"p90", "median", "mean"}:
    raise ValueError(f"Invalid global_rho_stat '{GLOBAL_RHO_STAT}'. Use 'p90', 'median', or 'mean'.")

MAKE_CNEES = _as_bool(config.get("make_cnees", False), False)
REF_GFF = config.get("ref_gff", "")
if MAKE_CNEES and not REF_GFF:
    raise ValueError("make_cnees=true requires ref_gff to be set in config.")
MAKE_CNEE_MAFS = _as_bool(config.get("make_cnee_mafs", False), False)
CNEE_CES_MERGE_GAP_BP = int(config.get("cnee_ces_merge_gap_bp", 5))
if CNEE_CES_MERGE_GAP_BP < 0:
    raise ValueError("cnee_ces_merge_gap_bp must be >= 0.")
CNEE_MIN_LEN_BP = int(config.get("cnee_min_len_bp", 50))
if CNEE_MIN_LEN_BP < 0:
    raise ValueError("cnee_min_len_bp must be >= 0.")
CNE_EXTRACT_FORMAT = str(config.get("cne_extract_format", "fasta")).strip().lower()
if CNE_EXTRACT_FORMAT in {"fa", "fna"}:
    CNE_EXTRACT_FORMAT = "fasta"
if CNE_EXTRACT_FORMAT not in {"fasta", "maf"}:
    raise ValueError("cne_extract_format must be 'fasta' or 'maf'.")
CNE_FASTA_HEADER = str(config.get("cne_fasta_header", "species-coords-id")).strip()
KEEP_CNEE_SIDECARS = _as_bool(config.get("keep_cnee_sidecars", False), False)

ALL_BED_TARGETS = expand(
    os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}.phastcon-conserved.bed"),
    zip,
    chromosome_group=REF_CHR_GROUPS_LIST,
    ref_chromosome=REF_CHROMOSOMES
)

ALL_CNEES_TARGETS = expand(
    os.path.join(CNEES_DIR, "{chromosome_group}", "{ref_chromosome}.cnees.bed"),
    zip,
    chromosome_group=REF_CHR_GROUPS_LIST,
    ref_chromosome=REF_CHROMOSOMES
)

ALL_CNEE_MAF_TARGETS = expand(
    os.path.join(CNEE_MAF_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.txt"),
    zip,
    chromosome_group=REF_CHR_GROUPS_LIST,
    ref_chromosome=REF_CHROMOSOMES
)

rule all:
    input:
        ALL_BED_TARGETS
        + (ALL_CNEES_TARGETS if MAKE_CNEES else [])
        + (ALL_CNEE_MAF_TARGETS if MAKE_CNEE_MAFS else [])

#############################################################################



rule picard_scatter_by_ns:
    input:
        ref_fasta = REF_FASTA
    output:
        ref_interval_file = NS_INTERVAL_FILE
    params:
        min_Ns_to_split_by = MIN_NS_TO_SPLIT_BY,
        rule_name = "picard_scatter_by_ns"
    log:
        job_log = os.path.join(LOG_DIR, "picard-scatter-by-ns.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "picard-scatter-by-ns.txt")
    resources:
        **getRuleResources("picard_scatter_by_ns")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                mem = config.get("rule_resources", {}).get("picard_scatter_by_ns", {}).get("mem_mb", 4000)

                cmd = [ "picard", "ScatterIntervalsByNs",
                        f"-Xmx{mem}m",
                        "--REFERENCE", input.ref_fasta,
                        "--OUTPUT_TYPE", "ACGT",
                        "--N", str(params.min_Ns_to_split_by),
                        "--OUTPUT", output.ref_interval_file ];

                COMMON.runCommand(cmd, log_stream, log_stream, params.rule_name);
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

####################

rule ns_to_bed:
    input:
        ref_interval_file = NS_INTERVAL_FILE
    output:
        chr_bed_file = os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", "windows", "Ns", MAF_CHR_PREFIX + "{ref_chromosome}-Ns.bed")
    log:
        job_log = os.path.join(LOG_DIR, "ns-to-bed.{chromosome_group}.{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "ns-to-bed.{chromosome_group}.{ref_chromosome}.txt")
    resources:
        **getRuleResources("ns_to_bed")
    run:
        with open(input.ref_interval_file, "r") as input_intervals, \
            open(output.chr_bed_file, "w") as out_file:
                for line in input_intervals:
                    if not line.startswith("@"):
                        line = line.strip().split("\t")
                        chrom, start, end = line[0], line[1], line[2]
                        if chrom == wildcards.ref_chromosome:
                            print(f"{chrom}:{start}-{end}", file=out_file)

####################

# in config, set something like:
# MIN_KEEP_REGION_LEN = config.get("min_keep_region_len", 6)

rule filter_ns_bed_minlen:
    input:
        chr_bed_file = os.path.join(
            OUTPUT_DIR, "beds", "{chromosome_group}", "windows", "Ns",
            MAF_CHR_PREFIX + "{ref_chromosome}-Ns.bed"
        )
    output:
        chr_bed_min_file = os.path.join(
            OUTPUT_DIR, "beds", "{chromosome_group}", "windows", "Ns",
            MAF_CHR_PREFIX + "{ref_chromosome}-Ns.min" + str(MIN_KEEP_REGION_LEN) + ".bed"
        )
    params:
        minlen = MIN_KEEP_REGION_LEN
    log:
        job_log = os.path.join(LOG_DIR, "filter-ns-bed-minlen.{chromosome_group}.{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "filter-ns-bed-minlen.{chromosome_group}.{ref_chromosome}.txt")
    resources:
        **getRuleResources("filter_ns_bed_minlen")
    shell:
        r"""
        set -euo pipefail
        awk -F'[:-]' -v MIN={params.minlen} 'NF==3 {{ if (($3-$2+1) >= MIN) print $0 }}' \
          {input.chr_bed_file} > {output.chr_bed_min_file}
        """

####################

rule ns_minlen_to_bed3:
    input:
        chr_bed_min_file = rules.filter_ns_bed_minlen.output.chr_bed_min_file
    output:
        chr_bed_min_fixed = os.path.join(
            OUTPUT_DIR, "beds", "{chromosome_group}", "windows", "Ns",
            MAF_CHR_PREFIX + "{ref_chromosome}-Ns.min" + str(MIN_KEEP_REGION_LEN) + ".fixed.bed"
        )
    log:
        job_log = os.path.join(LOG_DIR, "ns-minlen-to-bed3.{chromosome_group}.{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "ns-minlen-to-bed3.{chromosome_group}.{ref_chromosome}.txt")
    resources:
        **getRuleResources("ns_minlen_to_bed3")
    shell:
        r"""
        set -euo pipefail
        awk -F'[:-]' 'NF==3 {{print $1 "\t" ($2-1) "\t" $3}}' {input.chr_bed_min_file} > {output.chr_bed_min_fixed}
        """



####################



rule maf_index_chr:
    input:
        maf = os.path.join(MAF_SPLIT_BY_CHROM_DIR, "{chromosome_group}", "{ref_chromosome}.maf")
    output:
        maf_index_block = os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}", "{ref_chromosome}.maf.block.idx"),
        maf_index_scaff = os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}", "{ref_chromosome}.maf.scaff.idx")
    params:
        script_path = os.path.join(UTILS_DIR, "maf_index.py")
    log:
        job_log = os.path.join(LOG_DIR, "maf-index-chr.{chromosome_group}.{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "maf-index-chr.{chromosome_group}.{ref_chromosome}.txt")
    resources:
        **getRuleResources("maf_index_chr")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = ["python", params.script_path, input.maf, output.maf_index_block, output.maf_index_scaff]
                COMMON.runCommand(cmd, log_stream, log_stream, "maf_index_chr",
                                  wc=f"{wildcards.chromosome_group}.{wildcards.ref_chromosome}")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise



####################

rule maf_split_by_ns_minlen:
    input:
        maf = rules.maf_index_chr.input.maf,
        maf_index_block = rules.maf_index_chr.output.maf_index_block,
        bed3 = rules.ns_minlen_to_bed3.output.chr_bed_min_fixed
    output:
        manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.txt")
    params:
        script_path = os.path.join(UTILS_DIR, "maf_fetch.py"),
        rule_name = "maf_split_by_ns_minlen"
    log:
        job_log = os.path.join(LOG_DIR, "maf-split-by-ns-minlen.{chromosome_group}.{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "maf-split-by-ns-minlen.{chromosome_group}.{ref_chromosome}.txt")
    resources:
        **getRuleResources("maf_split_by_ns_minlen")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                p = min(int(resources.cpus_per_task), 2)

                # create output directory
                split_outdir = os.path.dirname(output.manifest)
                os.makedirs(split_outdir, exist_ok=True)

                cmd = [
                    "python", params.script_path,
                    input.maf,
                    input.maf_index_block,
                    input.bed3,
                    "-o", split_outdir,
                    "-p", str(p),
                    "-m", "block"
                ]
                COMMON.runCommand(
                    cmd, log_stream, log_stream, params.rule_name,
                    wc=f"{wildcards.chromosome_group}.{wildcards.ref_chromosome}"
                )

                # Write manifest listing produced chunk MAFs (basenames)
                mafs = sorted(glob.glob(os.path.join(split_outdir, f"{wildcards.ref_chromosome}-*.maf")))
                with open(output.manifest, "w") as out:
                    for m in mafs:
                        out.write(os.path.basename(m) + "\n")

            except Exception:
                traceback.print_exc(file=log_stream)
                raise




####################

checkpoint filter_maf_by_gap:
    input:
        manifest = rules.maf_split_by_ns_minlen.output.manifest
    output:
        filtered_manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.filtered.txt")
    params:
        max_gap_pct = MAX_GAP_PCT,
        outdir = lambda wc: os.path.join(MAF_SPLIT_NS_DIR, wc.chromosome_group, wc.ref_chromosome),
        rule_name = "filter_maf_by_gap"
    log:
        job_log = os.path.join(LOG_DIR, "filter-maf-by-gap.{chromosome_group}.{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "filter-maf-by-gap.{chromosome_group}.{ref_chromosome}.txt")
    resources:
        **getRuleResources("filter_maf_by_gap")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.filtered_manifest), exist_ok=True)
                with open(input.manifest) as mf:
                    all_lines = [l.strip() for l in mf if l.strip()]

                kept = []
                for fname in all_lines:
                    maf_path = os.path.join(params.outdir, fname)
                    if not os.path.exists(maf_path):
                        log_stream.write(f"MISSING: {maf_path}\n")
                        continue
                    gap = gap_pct_other(maf_path)
                    if gap <= params.max_gap_pct:
                        kept.append(fname)

                with open(output.filtered_manifest, "w") as out:
                    for k in kept:
                        out.write(k + "\n")

                log_stream.write(f"Kept {len(kept)} of {len(all_lines)} (gap_pct_other <= {params.max_gap_pct})\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise
####################



def beds_for_chr(wc):
    # Wait until the filtering checkpoint has run for this chromosome
    ckpt = checkpoints.filter_maf_by_gap.get(
        chromosome_group=wc.chromosome_group,
        ref_chromosome=wc.ref_chromosome
    )

    filtered_manifest = ckpt.output.filtered_manifest
    outdir = os.path.dirname(filtered_manifest)

    mafs = []
    with open(filtered_manifest) as mf:
        for line in mf:
            line = line.strip()
            if line:
                mafs.append(os.path.join(outdir, line))

    chunks = [os.path.splitext(os.path.basename(m))[0] for m in mafs]

    return expand(
        os.path.join(CONSERVE_DIR, wc.chromosome_group, wc.ref_chromosome, "{chunk}.conserved.bed"),
        chunk=chunks
    )
####################

def rho_values_for_chr(wc):
    ckpt = checkpoints.filter_maf_by_gap.get(
        chromosome_group=wc.chromosome_group,
        ref_chromosome=wc.ref_chromosome
    )
    filtered_manifest = ckpt.output.filtered_manifest
    outdir = os.path.dirname(filtered_manifest)

    rho_files = []
    with open(filtered_manifest) as mf:
        for line in mf:
            line = line.strip()
            if not line:
                continue
            chunk = os.path.splitext(os.path.basename(line))[0]
            rho_files.append(
                os.path.join(CONSERVE_DIR, wc.chromosome_group, wc.ref_chromosome, f"{chunk}.rho.txt")
            )
    return rho_files

rule phastcons_estimate_rho_chunk:
    input:
        maf = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.maf"),
        mod = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "{ref_chromosome}-corrected.mod"),
    output:
        stderr_file = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.phastcons.stderr.log"),
        rho_value = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.rho.txt"),
    params:
        outdir = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}"),
        rho_prefix = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.rho", "rho"),
        rule_name = "phastcons_estimate_rho_chunk",
    log:
        job_log = os.path.join(LOG_DIR, "phastcons-estimate-rho.{chromosome_group}.{ref_chromosome}.{chunk}.log"),
    resources:
        **getRuleResources("phastcons_chunk")
    run:
        import glob, os, re, shlex, subprocess, traceback, math

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(params.outdir, exist_ok=True)

                # Count UNIQUE sequence IDs
                cmd_uniq = (
                    f"awk '$1==\"s\"{{print $2}}' {shlex.quote(input.maf)} | sort -u | wc -l"
                )
                uniq_ids = int(subprocess.check_output(["bash", "-lc", cmd_uniq], text=True).strip())

                if uniq_ids < 2:
                    log_stream.write(f"SKIP: only {uniq_ids} unique sequence id(s) in {input.maf}\n")
                    with open(output.rho_value, "w") as rf:
                        rf.write("nan\n")
                    with open(output.stderr_file, "w") as sf:
                        sf.write(f"SKIP: only {uniq_ids} unique sequence id(s)\n")
                    return

                tmp_prefix = os.path.join(params.outdir, f".{wildcards.chunk}.estimate_rho.tmp")
                cmd = [
                    "phastCons",
                    input.maf,
                    input.mod,
                    "--estimate-rho",
                    tmp_prefix,
                ]

                log_stream.write(f"Running: {' '.join(shlex.quote(c) for c in cmd)}\n")
                log_stream.flush()

                result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                stderr_text = result.stderr or ""
                with open(output.stderr_file, "w") as sf:
                    sf.write(stderr_text)
                if stderr_text:
                    log_stream.write(stderr_text)
                    log_stream.flush()

                if result.returncode != 0:
                    # Don't crash pipeline; keep stderr, but make other outputs empty
                    log_stream.write(f"FAIL: phastCons exit {result.returncode} on {input.maf}\n")
                    with open(output.rho_value, "w") as rf:
                        rf.write("nan\n")
                    return

                # Parse rho directly from phastCons stderr text.
                rho = float("nan")
                try:
                    m_all = re.findall(r"rho\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", stderr_text)
                    if m_all:
                        rho = float(m_all[-1])
                except Exception:
                    pass

                with open(output.rho_value, "w") as rf:
                    rf.write(f"{rho}\n")

                # Remove temporary files written by --estimate-rho prefix.
                for fp in glob.glob(tmp_prefix + "*"):
                    try:
                        os.remove(fp)
                    except OSError:
                        pass

            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule global_rho:
    input:
        manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.filtered.txt"),
        mod = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "{ref_chromosome}-corrected.mod")
    output:
        global_rho = (
            os.path.join(RHO_STATS_DIR, "{chromosome_group}", "{ref_chromosome}", "global_rho", "selected.txt")
            if RHO_MODE == "estimate"
            else os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "global_rho.txt")
        )
    log:
        job_log = os.path.join(LOG_DIR, "global-rho.{chromosome_group}.{ref_chromosome}.log")
    resources:
        **getRuleResources("phastcons_chunk")
    run:
        import glob, math, os, re, shlex, subprocess

        with open(log.job_log, "w") as log_stream:
            log_stream.write("START global_rho\n")
            log_stream.write(f"RHO_MODE={RHO_MODE}\n")
            log_stream.flush()

            chr_dir = os.path.join(CONSERVE_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
            rho_chr_dir = os.path.join(RHO_STATS_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
            rho_global_dir = os.path.join(rho_chr_dir, "global_rho")
            maf_dir = os.path.join(MAF_SPLIT_NS_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
            os.makedirs(chr_dir, exist_ok=True)

            with open(input.manifest) as mf:
                chunks = [os.path.splitext(line.strip())[0] for line in mf if line.strip()]

            log_stream.write(f"Found {len(chunks)} filtered chunks\n")
            log_stream.flush()

            chunk_rho = {}
            if RHO_MODE == "estimate":
                os.makedirs(rho_chr_dir, exist_ok=True)
                os.makedirs(rho_global_dir, exist_ok=True)
                # Stage 1: per-chunk rho estimation from chunk MAFs.
                for chunk in chunks:
                    maf_path = os.path.join(maf_dir, f"{chunk}.maf")
                    rho_txt = os.path.join(rho_chr_dir, f"{chunk}.rho.txt")
                    rho_stderr = os.path.join(rho_chr_dir, f"{chunk}.phastcons.stderr.log")
                    tmp_prefix = os.path.join(rho_chr_dir, f".{chunk}.estimate_rho.tmp")

                    cmd = [
                        "phastCons",
                        maf_path,
                        input.mod,
                        "--estimate-rho",
                        tmp_prefix,
                    ]
                    log_stream.write(f"RHO: {' '.join(shlex.quote(c) for c in cmd)}\n")
                    log_stream.flush()

                    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                    stderr_text = result.stderr or ""
                    with open(rho_stderr, "w") as sf:
                        sf.write(stderr_text)

                    rho_val = float("nan")
                    if result.returncode == 0:
                        try:
                            m_all = re.findall(r"rho\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", stderr_text)
                            if m_all:
                                rho_val = float(m_all[-1])
                        except Exception:
                            pass
                    else:
                        log_stream.write(f"RHO FAIL: phastCons exit {result.returncode} on {maf_path}\n")

                    with open(rho_txt, "w") as rf:
                        rf.write(f"{rho_val}\n")
                    chunk_rho[chunk] = rho_val

                    for fp in glob.glob(tmp_prefix + "*"):
                        try:
                            os.remove(fp)
                        except OSError:
                            pass

                vals = [v for v in chunk_rho.values() if math.isfinite(v) and v > 0.0]
                if not vals:
                    rho_global = float("nan")
                    log_stream.write("NO_VALID_RHO_VALUES\n")
                    rho_mean = float("nan")
                    rho_median = float("nan")
                    rho_p90 = float("nan")
                else:
                    vals.sort()
                    idx = int(math.ceil(0.9 * len(vals)) - 1)
                    idx = max(0, min(idx, len(vals) - 1))
                    rho_p90 = vals[idx]
                    rho_median = float(vals[len(vals) // 2]) if len(vals) % 2 == 1 else float((vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2.0)
                    rho_mean = float(sum(vals) / len(vals))
                    if GLOBAL_RHO_STAT == "p90":
                        rho_global = rho_p90
                    elif GLOBAL_RHO_STAT == "median":
                        rho_global = rho_median
                    else:
                        rho_global = rho_mean
                    log_stream.write(f"GLOBAL_RHO_{GLOBAL_RHO_STAT.upper()}: {rho_global} (n={len(vals)})\n")

                with open(os.path.join(rho_global_dir, "mean.txt"), "w") as f:
                    f.write(f"{rho_mean}\n")
                with open(os.path.join(rho_global_dir, "median.txt"), "w") as f:
                    f.write(f"{rho_median}\n")
                with open(os.path.join(rho_global_dir, "p90.txt"), "w") as f:
                    f.write(f"{rho_p90}\n")
                with open(os.path.join(rho_global_dir, "selected.txt"), "w") as f:
                    f.write(f"{rho_global}\n")
            else:
                rho_global = FIXED_RHO
                log_stream.write(f"USING_FIXED_RHO: {rho_global}\n")

            with open(output.global_rho, "w") as out:
                out.write(f"{rho_global}\n")

####################

rule phastcons_apply_rho_chr:
    input:
        manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.filtered.txt"),
        mod = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "{ref_chromosome}-corrected.mod"),
        global_rho = rules.global_rho.output.global_rho
    output:
        chunks_done = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "chunks.done")
    log:
        job_log = os.path.join(LOG_DIR, "phastcons-apply-rho.{chromosome_group}.{ref_chromosome}.log")
    resources:
        **getRuleResources("phastcons_chunk")
    run:
        import math, os, shlex, subprocess, traceback

        with open(log.job_log, "w") as log_stream:
            try:
                chr_dir = os.path.join(CONSERVE_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
                maf_dir = os.path.join(MAF_SPLIT_NS_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
                os.makedirs(chr_dir, exist_ok=True)

                with open(input.global_rho) as gf:
                    rho_str = gf.read().strip()
                try:
                    rho_global = float(rho_str)
                except Exception:
                    rho_global = float("nan")

                with open(input.manifest) as mf:
                    chunks = [os.path.splitext(line.strip())[0] for line in mf if line.strip()]

                for chunk in chunks:
                    maf_path = os.path.join(maf_dir, f"{chunk}.maf")
                    bed = os.path.join(chr_dir, f"{chunk}.conserved.bed")
                    wig = os.path.join(chr_dir, f"{chunk}.scores.wig")
                    err = os.path.join(chr_dir, f"{chunk}.phastcons.stderr.log")

                    if not (math.isfinite(rho_global) and rho_global > 0.0):
                        open(bed, "w").close()
                        open(wig, "w").close()
                        with open(err, "w") as ef:
                            ef.write(f"SKIP: invalid rho '{rho_str}'\n")
                        continue

                    cmd = (
                        f"phastCons {shlex.quote(maf_path)} {shlex.quote(input.mod)} "
                        f"--rho {rho_global} "
                        f"--most-conserved {shlex.quote(bed)} "
                        f"> {shlex.quote(wig)} 2> {shlex.quote(err)}"
                    )
                    log_stream.write(f"CHUNK: {cmd}\n")
                    log_stream.flush()
                    rc = subprocess.call(["bash", "-lc", cmd])
                    if rc != 0:
                        log_stream.write(f"CHUNK FAIL: phastCons exit {rc} on {maf_path}\n")
                        open(bed, "w").close()
                        open(wig, "w").close()

                with open(output.chunks_done, "w") as done:
                    done.write("ok\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule phastcons_concat_chr:
    input:
        manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.filtered.txt"),
        chunks_done = rules.phastcons_apply_rho_chr.output.chunks_done
    output:
        chr_bed = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}.phastcon-conserved.bed")
    params:
        rule_name = "phastcons_concat_chr",
        cleanup_chunk_intermediates = CLEANUP_CHUNK_INTERMEDIATES
    log:
        job_log = os.path.join(LOG_DIR, "phastcons-concat.{chromosome_group}.{ref_chromosome}.log")
    resources:
        **getRuleResources("phastcons_concat_chr")
    run:
        import os, shlex, subprocess, traceback, glob

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.chr_bed), exist_ok=True)

                tmp = output.chr_bed + ".tmp"

                # Concatenate only non-empty chunk beds discovered on disk
                chunk_dir = os.path.join(CONSERVE_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
                bed_files = sorted(glob.glob(os.path.join(chunk_dir, "*.conserved.bed")))
                log_stream.write(f"Found {len(bed_files)} chunk bed files in {chunk_dir}\n")

                with open(tmp, "w") as out:
                    for bedf in bed_files:
                        if os.path.exists(bedf) and os.path.getsize(bedf) > 0:
                            with open(bedf) as bf:
                                for line in bf:
                                    line = line.strip()
                                    if line and not line.startswith("#"):
                                        out.write(line + "\n")

                # Sort into final output (works even if tmp is empty)
                sort_cmd = f"sort -k1,1 -k2,2n -k3,3n {shlex.quote(tmp)} > {shlex.quote(output.chr_bed)}"
                log_stream.write(f"Running: {sort_cmd}\n")
                log_stream.flush()

                rc = subprocess.call(["bash", "-lc", sort_cmd])
                if rc != 0:
                    raise RuntimeError(f"Sort failed with exit code {rc}")

                try:
                    os.remove(tmp)
                except OSError:
                    pass

                # Keep only final chromosome outputs by cleaning noisy chunk-level intermediates.
                if params.cleanup_chunk_intermediates:
                    cleanup_patterns = [
                        "*.conserved.bed",
                        "*.scores.wig",
                        "*.phastcons.stderr.log",
                    ]
                    for pat in cleanup_patterns:
                        for fp in glob.glob(os.path.join(chunk_dir, pat)):
                            try:
                                if os.path.isdir(fp):
                                    import shutil
                                    shutil.rmtree(fp, ignore_errors=True)
                                else:
                                    os.remove(fp)
                            except OSError:
                                pass

                    # Remove heavy 01 split chunks once chromosome-level bed is finalized.
                    split_dir = os.path.join(MAF_SPLIT_NS_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
                    try:
                        if os.path.isdir(split_dir):
                            import shutil
                            shutil.rmtree(split_dir, ignore_errors=True)
                    except OSError:
                        pass

            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule extract_cds_bed_chr:
    input:
        ref_gff = REF_GFF
    output:
        cds_bed = os.path.join(CNEES_DIR, "{chromosome_group}", "{ref_chromosome}", "{ref_chromosome}.cds.bed")
    log:
        job_log = os.path.join(LOG_DIR, "extract-cds-bed.{chromosome_group}.{ref_chromosome}.log")
    resources:
        **getRuleResources("get_conserved_sites")
    run:
        import os
        import traceback

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.cds_bed), exist_ok=True)
                with open(input.ref_gff) as gf, open(output.cds_bed, "w") as out:
                    for line in gf:
                        if not line or line.startswith("#"):
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 5:
                            continue
                        chrom, feature, start_s, end_s = parts[0], parts[2], parts[3], parts[4]
                        if chrom != wildcards.ref_chromosome or feature.upper() != "CDS":
                            continue
                        try:
                            start = int(start_s) - 1
                            end = int(end_s)
                        except ValueError:
                            continue
                        if end > start >= 0:
                            out.write(f"{chrom}\t{start}\t{end}\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule cnees_from_conserved_chr:
    input:
        conserved_bed = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}.phastcon-conserved.bed"),
        cds_bed = rules.extract_cds_bed_chr.output.cds_bed
    output:
        cnees_bed = os.path.join(CNEES_DIR, "{chromosome_group}", "{ref_chromosome}.cnees.bed")
    log:
        job_log = os.path.join(LOG_DIR, "cnees-from-conserved.{chromosome_group}.{ref_chromosome}.log")
    resources:
        **getRuleResources("get_conserved_sites")
    run:
        import os
        import traceback

        def parse_bed3(path, normalize_to=None):
            rows = []
            with open(path) as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 3:
                        continue
                    chrom = p[0]
                    if normalize_to is not None:
                        # This rule is per-chromosome; normalize aliases like NC_085107 -> NC_085107.1
                        base_norm = normalize_to.rsplit(".", 1)[0]
                        base_chrom = chrom.rsplit(".", 1)[0]
                        if chrom == normalize_to or base_chrom == base_norm:
                            chrom = normalize_to
                        else:
                            # Skip any unexpected chromosome labels in this per-chromosome file.
                            continue
                    try:
                        s = int(p[1])
                        e = int(p[2])
                    except ValueError:
                        continue
                    if e > s:
                        rows.append((chrom, s, e))
            return rows

        def merge_intervals(intervals, max_gap_bp=0):
            if not intervals:
                return []
            intervals = sorted(intervals, key=lambda x: (x[1], x[2]))
            merged = [list(intervals[0])]
            for _, s, e in intervals[1:]:
                # Merge overlaps and near-adjacent intervals within max_gap_bp.
                if s <= merged[-1][2] + max_gap_bp:
                    if e > merged[-1][2]:
                        merged[-1][2] = e
                else:
                    merged.append([wildcards.ref_chromosome, s, e])
            return [(c, s, e) for c, s, e in merged]

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.cnees_bed), exist_ok=True)
                conserved_raw = parse_bed3(input.conserved_bed, normalize_to=wildcards.ref_chromosome)
                conserved = merge_intervals(conserved_raw, CNEE_CES_MERGE_GAP_BP)
                cds = merge_intervals(parse_bed3(input.cds_bed, normalize_to=wildcards.ref_chromosome), 0)
                conserved.sort(key=lambda x: (x[1], x[2]))
                out_rows = []
                j = 0

                log_stream.write(
                    f"Conserved raw intervals: {len(conserved_raw)}; "
                    f"merged (gap<={CNEE_CES_MERGE_GAP_BP}bp): {len(conserved)}\\n"
                )

                for chrom, s, e in conserved:
                    cur = s
                    while j < len(cds) and cds[j][2] <= s:
                        j += 1
                    k = j
                    while k < len(cds) and cds[k][1] < e:
                        _, cs, ce = cds[k]
                        if cs > cur:
                            out_rows.append((chrom, cur, min(cs, e)))
                        cur = max(cur, ce)
                        if cur >= e:
                            break
                        k += 1
                    if cur < e:
                        out_rows.append((chrom, cur, e))

                # Final pass: merge adjacent/overlapping CNEE fragments.
                out_rows = merge_intervals(out_rows)
                log_stream.write(f"CNEE intervals after merge: {len(out_rows)}\n")

                with open(output.cnees_bed, "w") as out:
                    for chrom, s, e in out_rows:
                        if e > s:
                            out.write(f"{chrom}\t{s}\t{e}\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule cnees_to_bed4_chr:
    input:
        cnees_bed = rules.cnees_from_conserved_chr.output.cnees_bed
    output:
        cnees_bed4 = os.path.join(CNEES_DIR, "{chromosome_group}", "{ref_chromosome}.cnees.bed4")
    log:
        job_log = os.path.join(LOG_DIR, "cnees-to-bed4.{chromosome_group}.{ref_chromosome}.log")
    resources:
        **getRuleResources("get_conserved_sites")
    run:
        import os
        import traceback

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.cnees_bed4), exist_ok=True)
                n = 0
                kept = 0
                dropped = 0
                with open(input.cnees_bed) as inf, open(output.cnees_bed4, "w") as outf:
                    for line in inf:
                        if not line.strip() or line.startswith("#"):
                            continue
                        p = line.rstrip("\n").split("\t")
                        if len(p) < 3:
                            continue
                        chrom = p[0]
                        try:
                            s = int(p[1])
                            e = int(p[2])
                        except ValueError:
                            continue
                        if e <= s:
                            continue
                        if (e - s) <= CNEE_MIN_LEN_BP:
                            dropped += 1
                            continue
                        n += 1
                        kept += 1
                        cid = f"{wildcards.ref_chromosome}.cnee{n:07d}"
                        outf.write(f"{chrom}\t{s}\t{e}\t{cid}\n")
                log_stream.write(
                    f"Filtered CNEEs by length > {CNEE_MIN_LEN_BP} bp: "
                    f"kept={kept}, dropped={dropped}\\n"
                )
                log_stream.write(f"Wrote {n} CNEE intervals to BED4\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule cnee_mafs_chr:
    input:
        cnees_bed4 = rules.cnees_to_bed4_chr.output.cnees_bed4,
        maf = rules.maf_index_chr.input.maf,
        maf_index_block = rules.maf_index_chr.output.maf_index_block
    output:
        manifest = os.path.join(CNEE_MAF_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.txt")
    params:
        script_path = os.path.join(UTILS_DIR, "maf_fetch.py"),
        outdir = os.path.join(CNEE_MAF_DIR, "{chromosome_group}", "{ref_chromosome}"),
        rule_name = "cnee_mafs_chr"
    log:
        job_log = os.path.join(LOG_DIR, "cnee-mafs.{chromosome_group}.{ref_chromosome}.log")
    resources:
        **getRuleResources("maf_split_by_ns_minlen")
    run:
        import glob
        import os
        import traceback

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(params.outdir, exist_ok=True)
                # Remove stale per-CNEE outputs so manifest reflects current filtering.
                stale = glob.glob(os.path.join(params.outdir, f"{wildcards.ref_chromosome}.cnee*.maf"))
                stale += glob.glob(os.path.join(params.outdir, f"{wildcards.ref_chromosome}.cnee*.fa"))
                for fp in stale:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass
                p = min(int(resources.cpus_per_task), 4)
                cmd = [
                    "python3", params.script_path,
                    input.maf,
                    input.maf_index_block,
                    input.cnees_bed4,
                    "-o", params.outdir,
                    "-p", str(p),
                    "-m", "block",
                    "-b", "id",
                ]
                if CNE_EXTRACT_FORMAT == "fasta":
                    cmd += ["-f", "-fh", CNE_FASTA_HEADER]
                COMMON.runCommand(
                    cmd, log_stream, log_stream, params.rule_name,
                    wc=f"{wildcards.chromosome_group}.{wildcards.ref_chromosome}"
                )

                ext = "fa" if CNE_EXTRACT_FORMAT == "fasta" else "maf"
                outs = sorted(glob.glob(os.path.join(params.outdir, f"{wildcards.ref_chromosome}.cnee*.{ext}")))

                if CNE_EXTRACT_FORMAT == "fasta":
                    import re

                    dropped_files = 0
                    kept_outs = []

                    for fp in outs:
                        seen = set()
                        has_dup = False
                        with open(fp, "r") as inf:
                            for line in inf:
                                if not line.startswith(">"):
                                    continue
                                token = line[1:].strip().split()[0]
                                species = token.split(":", 1)[0]
                                # Duplicate key is genus+species (first two underscore-separated tokens).
                                m = re.match(r"^([A-Za-z]+)_([A-Za-z]+)", species)
                                species_key = f"{m.group(1)}_{m.group(2)}" if m else species
                                if species_key in seen:
                                    has_dup = True
                                    break
                                seen.add(species_key)
                        if has_dup:
                            try:
                                os.remove(fp)
                            except OSError:
                                pass
                            dropped_files += 1
                        else:
                            kept_outs.append(fp)

                    outs = kept_outs
                    log_stream.write(
                        "Filtered CNEE FASTAs: "
                        f"dropped_files_with_duplicate_species={dropped_files}\n"
                    )

                with open(output.manifest, "w") as out:
                    for m in outs:
                        out.write(os.path.basename(m) + "\n")
                log_stream.write(f"Wrote manifest with {len(outs)} CNEE {CNE_EXTRACT_FORMAT.upper()} files\n")

                # Keep only essential outputs by default: manifest + extracted alignments.
                if not KEEP_CNEE_SIDECARS:
                    for side_pat in ("*.tsv", "*.log"):
                        for fp in glob.glob(os.path.join(params.outdir, side_pat)):
                            try:
                                os.remove(fp)
                            except OSError:
                                pass

            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule phastcons_chunk:
    input:
        maf = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.maf"),
        mod = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "{ref_chromosome}-corrected.mod"),
        rho_value = rules.phastcons_estimate_rho_chunk.output.rho_value,
        global_rho = rules.global_rho.output.global_rho,
    output:
        bed = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.conserved.bed"),
        wig = temp(os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.scores.wig")),
        err = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.phastcons.stderr.log"),
    params:
        outdir = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}"),
        rule_name = "phastcons_chunk",
    log:
        job_log = os.path.join(LOG_DIR, "phastcons-chunk.{chromosome_group}.{ref_chromosome}.{chunk}.log"),
    resources:
        **getRuleResources("phastcons_chunk")
    run:
        import os, shlex, subprocess, traceback, math

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(params.outdir, exist_ok=True)

                # Count UNIQUE sequence IDs
                cmd_uniq = (
                    f"awk '$1==\"s\"{{print $2}}' {shlex.quote(input.maf)} | sort -u | wc -l"
                )
                uniq_ids = int(subprocess.check_output(["bash", "-lc", cmd_uniq], text=True).strip())

                if uniq_ids < 2:
                    log_stream.write(f"SKIP: only {uniq_ids} unique sequence id(s) in {input.maf}\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    with open(output.err, "w") as ef:
                        ef.write(f"SKIP: only {uniq_ids} unique sequence id(s)\n")
                    return

                with open(input.rho_value) as rf:
                    rho_str = rf.read().strip()
                try:
                    rho_chunk = float(rho_str)
                except Exception:
                    rho_chunk = float("nan")

                with open(input.global_rho) as gf:
                    global_str = gf.read().strip()
                try:
                    rho_global = float(global_str)
                except Exception:
                    rho_global = float("nan")

                if not (rho_global > 0.0 and math.isfinite(rho_global)):
                    log_stream.write(f"SKIP: invalid global rho '{global_str}'\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    with open(output.err, "w") as ef:
                        ef.write(f"SKIP: invalid global rho '{global_str}'\n")
                    return

                if not (rho_chunk > 0.0 and math.isfinite(rho_chunk)):
                    log_stream.write(f"SKIP: invalid chunk rho '{rho_str}' for {input.maf}\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    with open(output.err, "w") as ef:
                        ef.write(f"SKIP: invalid chunk rho '{rho_str}'\n")
                    return

                if rho_chunk > rho_global:
                    log_stream.write(f"SKIP: rho_chunk {rho_chunk} > global_rho {rho_global}\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    with open(output.err, "w") as ef:
                        ef.write(f"SKIP: rho_chunk {rho_chunk} > global_rho {rho_global}\n")
                    return

                # KEEP post probs (do NOT add --no-post-probs)
                cmd = (
                    f"phastCons {shlex.quote(input.maf)} {shlex.quote(input.mod)} "
                    f"--rho {rho_global} "
                    f"--most-conserved {shlex.quote(output.bed)} "
                    f"> {shlex.quote(output.wig)} 2> {shlex.quote(output.err)}"
                )

                log_stream.write(f"Running: {cmd}\n")
                log_stream.flush()

                rc = subprocess.call(["bash", "-lc", cmd])

                if rc != 0:
                    # Don't crash pipeline; keep stderr, but make other outputs empty
                    log_stream.write(f"FAIL: phastCons exit {rc} on {input.maf}\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    return

            except Exception:
                traceback.print_exc(file=log_stream)
                raise
