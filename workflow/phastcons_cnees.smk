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


def parse_species_list(species_text: str):
    species = []
    seen = set()
    for sp in species_text.split(','):
        sp = sp.strip()
        if sp and sp not in seen:
            species.append(sp)
            seen.add(sp)
    return species


def picard_dict_path(ref_path: str) -> str:
    lower = ref_path.lower()
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if lower.endswith(suffix):
            return ref_path[:-len(suffix)] + ".dict"
    return ref_path + ".dict"


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

if config.get("__master_workflow__", False):
    setup = config["__pipeline_setup__"];
    MAIN = setup["MAIN"];
    DRY_RUN = setup["DRY_RUN"];
    OUTPUT_DIR = setup["OUTPUT_DIR"];
    LOG_DIR = setup["LOG_DIR"];
    TMPDIR = setup["TMPDIR"];
    LOG_LEVEL = setup["LOG_LEVEL"];
    LOG_VERBOSITY = setup["LOG_VERBOSITY"];
else:
    MAIN, DRY_RUN, OUTPUT_DIR, LOG_DIR, TMPDIR, LOG_LEVEL, LOG_VERBOSITY = COMMON.pipelineSetup(config, sys.argv, version_flag, info_flag, config_flag, debug, workflow);
# Setup the pipeline once in the master workflow, or locally when this module
# is run as a standalone Snakefile.

MLOG = logging.getLogger('META')
# Setup logging if debugging

MODULE_DIR = os.path.dirname(workflow.snakefile);
PIPELINE_DIR = os.path.dirname(MODULE_DIR) if os.path.basename(MODULE_DIR) == "workflow" else MODULE_DIR;
# Resolve the repository root whether this file is run directly from workflow/
# or included from the top-level Snakefile.

UTILS_DIR = os.path.join(PIPELINE_DIR, "utils");
# The directory where the utility scripts are located

getRuleResources = partial(COMMON.getResources, config)
# This maps the function to get rule resources from the config file
# so we don't have to pass config each time we call it

#############################################################################
# Input files and output paths

REFERENCE_DIR = os.path.join(OUTPUT_DIR, "00-reference")
REFERENCE_GROUP_BEDS_DIR = os.path.join(REFERENCE_DIR, "group-beds")

MAF_PREP_DIR = os.path.join(OUTPUT_DIR, "01-maf-prep")
MAF_INDEX_DIR = os.path.join(MAF_PREP_DIR, "maf-index")
NS_INTERVALS_DIR = os.path.join(MAF_PREP_DIR, "ns-intervals")
NS_INTERVALS_RAW_DIR = os.path.join(NS_INTERVALS_DIR, "raw")
NS_INTERVALS_FILTERED_DIR = os.path.join(NS_INTERVALS_DIR, "filtered")
CHUNK_BEDS_DIR = os.path.join(MAF_PREP_DIR, "chunk-beds")
CHUNKED_MAFS_DIR = os.path.join(MAF_PREP_DIR, "chunked-mafs")
MAF_CHUNK_SUMMARY_DIR = os.path.join(MAF_PREP_DIR, "summary")
# Kept outside chunked-mafs/ deliberately: phastcons_concat_chr rmtree's the whole
# per-chromosome chunked-mafs directory once done (cleanup_chunk_intermediates), so
# a durable pre/post-filter count summary needs to live somewhere that survives that.
MAF_SPLIT_BY_CHROM_DIR = COMMON.getOptionalConfigPath(
    config,
    "maf_split_chr_dir",
    os.path.join(MAF_PREP_DIR, "maf-by-chromosome"),
)

NEUTRAL_MODEL_DIR = os.path.join(OUTPUT_DIR, "02-neutral-model")
PHYLOFIT_DIR = COMMON.getOptionalConfigPath(
    config,
    "phylofit_chr_dir",
    os.path.join(NEUTRAL_MODEL_DIR, "phylofit"),
)

PHASTCONS_DIR = os.path.join(OUTPUT_DIR, "04-phastcons")
CONSERVE_DIR = os.path.join(PHASTCONS_DIR, "regions")
RHO_STATS_DIR = os.path.join(PHASTCONS_DIR, "rho")

CNEES_ROOT_DIR = os.path.join(OUTPUT_DIR, "05-cnees", "phastcons")
CNEES_DIR = os.path.join(CNEES_ROOT_DIR, "bed")
CNEES_SUMMARY_DIR = os.path.join(CNEES_ROOT_DIR, "summary")
REF_FASTA = config.get("ref_fasta") or ""
if not REF_FASTA:
    raise ValueError("run_phastcons=true requires ref_fasta to be set in config.")
REF_INDEX = COMMON.getOptionalConfigPath(
    config,
    "ref_fasta_index",
    COMMON.getOptionalConfigPath(config, "ref_genome_index", REF_FASTA + ".fai"),
)
if os.path.abspath(REF_INDEX) != os.path.abspath(REF_FASTA + ".fai"):
    raise ValueError(
        f"ref_fasta_index must match ref_fasta + '.fai' because samtools faidx writes next to the FASTA. "
        f"Got ref_fasta={REF_FASTA}, ref_fasta_index={REF_INDEX}"
    )
REF_DICT = picard_dict_path(REF_FASTA)
REF_GFF = config.get("ref_gff", "")
TREE_FILE = config.get("tree_file", "")
REF_CHROMOSOME_GROUPS = config["ref_chromosome_groups"]

def get_chr_maf_for_stage(wildcards):
    if "checkpoints" in globals() and hasattr(checkpoints, "maf_split_chr_by_group"):
        cp_output = checkpoints.maf_split_chr_by_group.get(chromosome_group=wildcards.chromosome_group).output
        group_dir = os.path.join(MAF_SPLIT_BY_CHROM_DIR, wildcards.chromosome_group)
        manifest_file = cp_output.maf_manifest
        with open(manifest_file, "r") as mf:
            expected = [line.strip() for line in mf if line.strip()]
        target = f"{wildcards.ref_chromosome}.maf"
        for fname in expected:
            if fname == target:
                return os.path.join(group_dir, fname)
        raise ValueError(f"No chromosome MAF for {wildcards.ref_chromosome} listed in {manifest_file}")
    return os.path.join(MAF_SPLIT_BY_CHROM_DIR, wildcards.chromosome_group, f"{wildcards.ref_chromosome}.maf")

# Directory containing chromosome/scaffold MAFs expected as:
# {MAF_SPLIT_BY_CHROM_DIR}/{chromosome_group}/{ref_chromosome}.maf

MAF_REF_ID = config["maf_ref_id"]
MAF_CHR_PREFIX = config["maf_chr_prefix"]
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"]
MAF_REF_PREFIX = MAF_REF_ID + MAF_REF_CHR_JOINER + MAF_CHR_PREFIX

#############################################################################
# Other params

SPLIT_STRATEGY = str(config.get("split_strategy", "ns")).strip().lower()
if SPLIT_STRATEGY not in {"ns", "fixed_windows"}:
    raise ValueError(f"Invalid split_strategy '{SPLIT_STRATEGY}'. Use 'ns' or 'fixed_windows'.")

PICARD_FASTA_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
)
if SPLIT_STRATEGY == "ns" and not REF_FASTA.lower().endswith(PICARD_FASTA_SUFFIXES):
    raise ValueError(
        "split_strategy=Ns requires ref_fasta to have a Picard-recognized FASTA extension "
        f"({', '.join(PICARD_FASTA_SUFFIXES)}). Got: {REF_FASTA}"
    )

MIN_NS_TO_SPLIT_BY = int(config.get("min_Ns_to_split_by", 100))
MIN_KEEP_REGION_LEN = int(config.get("min_keep_region_len", 6))
WINDOW_SIZE_BP = config.get("window_size_bp", None)
WINDOW_OVERLAP_BP = int(config.get("window_overlap_bp", 0))
if WINDOW_OVERLAP_BP < 0:
    raise ValueError("window_overlap_bp must be >= 0.")
if SPLIT_STRATEGY == "fixed_windows":
    if WINDOW_SIZE_BP is None:
        raise ValueError("split_strategy=fixed_windows requires window_size_bp in config.")
    WINDOW_SIZE_BP = int(WINDOW_SIZE_BP)
    if WINDOW_SIZE_BP <= 0:
        raise ValueError("window_size_bp must be > 0.")
    if WINDOW_OVERLAP_BP >= WINDOW_SIZE_BP:
        raise ValueError("window_overlap_bp must be smaller than window_size_bp.")
else:
    WINDOW_SIZE_BP = int(config.get("window_size_bp", 1000000))
WINDOW_STEP_BP = WINDOW_SIZE_BP - WINDOW_OVERLAP_BP

SPLIT_LABEL = (
    f"ns_min{MIN_KEEP_REGION_LEN}"
    if SPLIT_STRATEGY == "ns"
    else f"fixed_windows_w{WINDOW_SIZE_BP}_o{WINDOW_OVERLAP_BP}"
)

NS_INTERVAL_FILE = os.path.join(NS_INTERVALS_DIR, f"{SPLIT_LABEL}.txt")
CHUNK_BED_DIR = os.path.join(CHUNK_BEDS_DIR, SPLIT_LABEL)
MAF_SPLIT_NS_DIR = os.path.join(CHUNKED_MAFS_DIR, SPLIT_LABEL)

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
PHASTCONS_STANDALONE = not bool(config.get("__master_workflow__", False))
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

DEBUG_KEEP_INTERMEDIATES = _as_bool(config.get("debug_keep_intermediates", False), False)
CLEANUP_CHUNK_INTERMEDIATES = _as_bool(
    config.get("cleanup_chunk_intermediates", not DEBUG_KEEP_INTERMEDIATES),
    not DEBUG_KEEP_INTERMEDIATES,
)
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

MAKE_CNEES = _as_bool(config.get("build_cnees", config.get("make_cnees", True)), True)
if MAKE_CNEES and not REF_GFF:
    raise ValueError("build_cnees=true requires ref_gff to be set in config.")
CNEE_OUTPUT_FORMAT = str(config.get("cnee_output_format", "fasta")).strip().lower()
legacy_make_cnee_mafs = config.get("make_cnee_mafs", None)
legacy_cnee_extract_format = config.get("cnee_extract_format", config.get("cne_extract_format", None))
if CNEE_OUTPUT_FORMAT in {"fa", "fna"}:
    CNEE_OUTPUT_FORMAT = "fasta"
if CNEE_OUTPUT_FORMAT == "none" and legacy_make_cnee_mafs is not None and _as_bool(legacy_make_cnee_mafs, False):
    CNEE_OUTPUT_FORMAT = str(legacy_cnee_extract_format or "fasta").strip().lower()
    if CNEE_OUTPUT_FORMAT in {"fa", "fna"}:
        CNEE_OUTPUT_FORMAT = "fasta"
if CNEE_OUTPUT_FORMAT not in {"none", "fasta", "maf"}:
    raise ValueError("cnee_output_format must be 'none', 'fasta', or 'maf'.")
if CNEE_OUTPUT_FORMAT != "none" and not MAKE_CNEES:
    raise ValueError("cnee_output_format requires build_cnees=true.")
MAKE_CNEE_MAFS = CNEE_OUTPUT_FORMAT != "none"
CNEE_MAF_DIR = os.path.join(CNEES_ROOT_DIR, CNEE_OUTPUT_FORMAT if CNEE_OUTPUT_FORMAT == "fasta" else "maf")
CNEE_CES_MERGE_GAP_BP = int(config.get("cnee_ces_merge_gap_bp", 5))
if CNEE_CES_MERGE_GAP_BP < 0:
    raise ValueError("cnee_ces_merge_gap_bp must be >= 0.")
CNEE_MIN_LEN_BP = int(config.get("cnee_min_len_bp", 50))
if CNEE_MIN_LEN_BP < 0:
    raise ValueError("cnee_min_len_bp must be >= 0.")
CNEE_FASTA_HEADER = str(config.get("cnee_fasta_header", config.get("cne_fasta_header", "species-coords-id"))).strip()
CNEE_EXPECTED_SPECIES = []
if CNEE_OUTPUT_FORMAT == "fasta":
    explicit_species = parse_species_list(str(config.get("cnee_expected_species") or "").strip())
    explicit_species_file = str(config.get("cnee_expected_species_file") or "").strip()
    if explicit_species_file:
        with open(explicit_species_file, "r", encoding="utf-8") as fp:
            for line in fp:
                sp = line.strip()
                if sp and not sp.startswith("#") and sp not in explicit_species:
                    explicit_species.append(sp)
    if explicit_species:
        CNEE_EXPECTED_SPECIES = explicit_species
    elif TREE_FILE:
        CNEE_EXPECTED_SPECIES = parse_newick_tip_names(TREE_FILE)
        if not CNEE_EXPECTED_SPECIES:
            raise ValueError(f"No tip labels parsed from tree_file '{TREE_FILE}' for CNEE FASTA extraction.")
    else:
        raise ValueError(
            "cnee_output_format=fasta requires species information: set "
            "cnee_expected_species, cnee_expected_species_file, or tree_file."
        )
KEEP_CNEE_SIDECARS = _as_bool(
    config.get("keep_cnee_sidecars", DEBUG_KEEP_INTERMEDIATES),
    DEBUG_KEEP_INTERMEDIATES,
)

ALL_BED_TARGETS = expand(
    os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}.bed"),
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

if PHASTCONS_STANDALONE:
    localrules: all

    rule all:
        input:
            ALL_BED_TARGETS
            + (ALL_CNEES_TARGETS if MAKE_CNEES else [])
            + (ALL_CNEE_MAF_TARGETS if MAKE_CNEE_MAFS else [])

wildcard_constraints:
    chromosome_group = r"[^/]+",
    ref_chromosome = r"[^/]+",
    chunk = r"[^/]+"

#############################################################################

def chunk_bed_for_chr(wc):
    if SPLIT_STRATEGY == "ns":
        return os.path.join(CHUNK_BED_DIR, wc.chromosome_group, f"{wc.ref_chromosome}.bed")
    return os.path.join(CHUNK_BEDS_DIR, "fixed_windows", wc.chromosome_group, f"{wc.ref_chromosome}.bed")


if not bool(config.get("__ref_fasta_index_rule_defined__", False)):
    config["__ref_fasta_index_rule_defined__"] = True

    rule ref_fasta_index:
        input:
            ref_fasta = REF_FASTA
        output:
            ref_fasta_index = REF_INDEX
        log:
            job_log = os.path.join(LOG_DIR, "ref_fasta_index", "run.log")
        benchmark:
            os.path.join(LOG_DIR, "benchmarks", "ref_fasta_index", "run.txt")
        resources:
            **getRuleResources("ref_fasta_index")
        run:
            with open(log.job_log, "w") as log_stream:
                try:
                    cmd = ["samtools", "faidx", input.ref_fasta]
                    COMMON.runCommand(cmd, log_stream, log_stream, "ref_fasta_index")
                except Exception:
                    traceback.print_exc(file=log_stream)
                    raise


if not bool(config.get("__ref_fasta_dict_rule_defined__", False)):
    config["__ref_fasta_dict_rule_defined__"] = True

    rule ref_fasta_dict:
        input:
            ref_fasta = REF_FASTA
        output:
            ref_fasta_dict = REF_DICT
        log:
            job_log = os.path.join(LOG_DIR, "ref_fasta_dict", "run.log")
        benchmark:
            os.path.join(LOG_DIR, "benchmarks", "ref_fasta_dict", "run.txt")
        resources:
            **getRuleResources("ref_fasta_dict")
        run:
            with open(log.job_log, "w") as log_stream:
                try:
                    mem = config.get("rule_resources", {}).get("ref_fasta_dict", {}).get("mem_mb", 4000)
                    cmd = [
                        "picard",
                        "CreateSequenceDictionary",
                        f"-Xmx{mem}m",
                        "--REFERENCE", input.ref_fasta,
                        "--OUTPUT", output.ref_fasta_dict,
                    ]
                    COMMON.runCommand(cmd, log_stream, log_stream, "ref_fasta_dict")
                except Exception:
                    traceback.print_exc(file=log_stream)
                    raise


rule picard_scatter_by_ns:
    input:
        ref_fasta = REF_FASTA,
        ref_fasta_index = REF_INDEX,
        ref_fasta_dict = REF_DICT
    output:
        ref_interval_file = NS_INTERVAL_FILE
    params:
        min_Ns_to_split_by = MIN_NS_TO_SPLIT_BY,
        rule_name = "picard_scatter_by_ns"
    log:
        job_log = os.path.join(LOG_DIR, "picard_scatter_by_ns", f"{SPLIT_LABEL}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "picard_scatter_by_ns", f"{SPLIT_LABEL}.txt")
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
        chr_bed_file = os.path.join(NS_INTERVALS_RAW_DIR, "{chromosome_group}", "{ref_chromosome}.bed")
    log:
        job_log = os.path.join(LOG_DIR, "ns_to_bed", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "ns_to_bed", "{chromosome_group}", "{ref_chromosome}.txt")
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
            NS_INTERVALS_RAW_DIR, "{chromosome_group}", "{ref_chromosome}.bed"
        )
    output:
        chr_bed_min_file = os.path.join(
            NS_INTERVALS_FILTERED_DIR, "{chromosome_group}", "{ref_chromosome}.bed"
        )
    params:
        minlen = MIN_KEEP_REGION_LEN
    log:
        job_log = os.path.join(LOG_DIR, "filter_ns_bed_minlen", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "filter_ns_bed_minlen", "{chromosome_group}", "{ref_chromosome}.txt")
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
        chr_bed_min_fixed = os.path.join(CHUNK_BED_DIR, "{chromosome_group}", "{ref_chromosome}.bed")
    log:
        job_log = os.path.join(LOG_DIR, "ns_minlen_to_bed3", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "ns_minlen_to_bed3", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("ns_minlen_to_bed3")
    shell:
        r"""
        set -euo pipefail
        awk -F'[:-]' 'NF==3 {{print $1 "\t" ($2-1) "\t" $3}}' {input.chr_bed_min_file} > {output.chr_bed_min_fixed}
        """



####################

rule fixed_windows_bed:
    input:
        ref_index = REF_INDEX
    output:
        chr_bed = os.path.join(CHUNK_BEDS_DIR, "fixed_windows", "{chromosome_group}", "{ref_chromosome}.bed")
    log:
        job_log = os.path.join(LOG_DIR, "fixed_windows_bed", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "fixed_windows_bed", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("fixed_windows_bed")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.chr_bed), exist_ok=True)
                chrom_len = None
                with open(input.ref_index) as inf:
                    for line in inf:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 2:
                            continue
                        if parts[0] == wildcards.ref_chromosome:
                            chrom_len = int(parts[1])
                            break
                if chrom_len is None:
                    raise ValueError(f"Chromosome {wildcards.ref_chromosome} not found in {input.ref_index}")

                with open(output.chr_bed, "w") as out:
                    start = 0
                    while start < chrom_len:
                        end = min(start + WINDOW_SIZE_BP, chrom_len)
                        out.write(f"{wildcards.ref_chromosome}\t{start}\t{end}\n")
                        if end >= chrom_len:
                            break
                        start += WINDOW_STEP_BP

                log_stream.write(
                    f"window_size_bp={WINDOW_SIZE_BP}; window_overlap_bp={WINDOW_OVERLAP_BP}; step_bp={WINDOW_STEP_BP}\n"
                )
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule maf_index_chr:
    input:
        maf = get_chr_maf_for_stage
    output:
        maf_index_block = os.path.join(MAF_INDEX_DIR, "{chromosome_group}", "{ref_chromosome}.maf.block.idx"),
        maf_index_scaff = os.path.join(MAF_INDEX_DIR, "{chromosome_group}", "{ref_chromosome}.maf.scaff.idx")
    log:
        job_log = os.path.join(LOG_DIR, "maf_index_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "maf_index_chr", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("maf_index_chr")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = ["mafutils", "index",
                       input.maf, output.maf_index_block, output.maf_index_scaff]
                COMMON.runCommand(cmd, log_stream, log_stream, "maf_index_chr",
                                  wc=f"{wildcards.chromosome_group}.{wildcards.ref_chromosome}")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise



####################

rule maf_split_chunks:
    input:
        maf = rules.maf_index_chr.input.maf,
        maf_index_block = rules.maf_index_chr.output.maf_index_block,
        bed3 = chunk_bed_for_chr
    output:
        manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.txt")
    params:
        rule_name = "maf_split_chunks"
    log:
        job_log = os.path.join(LOG_DIR, "maf_split_chunks", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "maf_split_chunks", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("maf_split_chunks")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                p = min(int(resources.cpus_per_task), 2)

                # create output directory
                split_outdir = os.path.dirname(output.manifest)
                os.makedirs(split_outdir, exist_ok=True)

                cmd = [
                    "mafutils", "fetch",
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
        manifest = rules.maf_split_chunks.output.manifest
    output:
        filtered_manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.filtered.txt"),
        filter_summary = os.path.join(MAF_CHUNK_SUMMARY_DIR, "{chromosome_group}", "{ref_chromosome}.maf-chunk-filter-summary.tsv")
    params:
        max_gap_pct = MAX_GAP_PCT,
        outdir = lambda wc: os.path.join(MAF_SPLIT_NS_DIR, wc.chromosome_group, wc.ref_chromosome),
        rule_name = "filter_maf_by_gap"
    log:
        job_log = os.path.join(LOG_DIR, "filter_maf_by_gap", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "filter_maf_by_gap", "{chromosome_group}", "{ref_chromosome}.txt")
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

                # Written here (rather than just relying on manifest.txt/manifest.filtered.txt)
                # because phastcons_concat_chr rmtree's the whole chunked-mafs directory for
                # this chromosome once it's done with it - this summary lives elsewhere
                # (MAF_CHUNK_SUMMARY_DIR) specifically so it survives that cleanup.
                os.makedirs(os.path.dirname(output.filter_summary), exist_ok=True)
                with open(output.filter_summary, "w") as sf:
                    sf.write("var\tfilter.cat\tvalue\n")
                    sf.write(f"num.chunks\tpre.filter\t{len(all_lines)}\n")
                    sf.write(f"num.chunks\tpost.filter\t{len(kept)}\n")

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

# rho_values_for_chr() and rule phastcons_estimate_rho_chunk moved to
# workflow/stash.smk - dead code, orphaned since global_rho inlined this work.

####################

rule global_rho:
    input:
        manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.filtered.txt"),
        mod = PHYLOFIT_ACTIVE_MODEL_PATH
    output:
        global_rho = (
            os.path.join(RHO_STATS_DIR, "{chromosome_group}", "{ref_chromosome}", "global_rho", "selected.txt")
            if RHO_MODE == "estimate"
            else os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "global_rho.txt")
        )
    log:
        job_log = os.path.join(LOG_DIR, "global_rho", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "global_rho", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **(getRuleResources("global_rho") if RHO_MODE == "estimate" else getRuleResources("default"))
    run:
        import glob, math, os, re, shlex, subprocess
        from concurrent.futures import ThreadPoolExecutor, as_completed

        with open(log.job_log, "w") as log_stream:
            log_stream.write("START global_rho\n")
            log_stream.write(f"RHO_MODE={RHO_MODE}\n")
            log_stream.flush()

            chr_dir = os.path.join(CONSERVE_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
            rho_chr_dir = os.path.join(RHO_STATS_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
            rho_global_dir = os.path.join(rho_chr_dir, "global_rho")
            maf_dir = os.path.join(MAF_SPLIT_NS_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
            rho_log_dir = os.path.join(LOG_DIR, "global_rho", wildcards.chromosome_group, wildcards.ref_chromosome)
            os.makedirs(chr_dir, exist_ok=True)
            os.makedirs(rho_log_dir, exist_ok=True)

            with open(input.manifest) as mf:
                chunks = [os.path.splitext(line.strip())[0] for line in mf if line.strip()]

            log_stream.write(f"Found {len(chunks)} filtered chunks\n")
            log_stream.flush()

            chunk_rho = {}
            if RHO_MODE == "estimate":
                os.makedirs(rho_chr_dir, exist_ok=True)
                os.makedirs(rho_global_dir, exist_ok=True)
                os.makedirs(resources.tmpdir, exist_ok=True)

                # Stage 1: per-chunk rho estimation from chunk MAFs, run concurrently
                # (each --estimate-rho call is an independent, comparatively expensive
                # EM fit, unlike the fixed-rho apply in run_phastcons_chr). Per-chunk
                # scratch files live on resources.tmpdir (node-local scratch), not the
                # shared output filesystem -- with 100k+ chunks per chromosome, the
                # small-file traffic from phastCons's own --estimate-rho working files
                # was bottlenecking on shared-filesystem metadata contention rather
                # than CPU, which is why run_phastcons_chr (3 known files per chunk,
                # no globbing) doesn't hit this even at the same chunk counts.
                def estimate_chunk_rho(chunk):
                    maf_path = os.path.join(maf_dir, f"{chunk}.maf")
                    tmp_prefix = os.path.join(
                        resources.tmpdir,
                        f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}-{chunk}.estimate_rho.tmp",
                    )

                    cmd = [
                        "phastCons",
                        maf_path,
                        input.mod,
                        "--estimate-rho",
                        tmp_prefix,
                    ]

                    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                    stderr_text = result.stderr or ""

                    if result.returncode != 0:
                        rho_stderr = os.path.join(rho_log_dir, f"{chunk}.stderr.log")
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

                    for fp in glob.glob(tmp_prefix + "*"):
                        try:
                            os.remove(fp)
                        except OSError:
                            pass

                    if DEBUG_KEEP_INTERMEDIATES:
                        rho_txt = os.path.join(rho_chr_dir, f"{chunk}.rho.txt")
                        with open(rho_txt, "w") as rf:
                            rf.write(f"{rho_val}\n")

                    return chunk, maf_path, result.returncode, rho_val

                max_workers = max(1, min(int(resources.cpus_per_task), len(chunks) if chunks else 1))
                log_stream.write(
                    f"Estimating rho for {len(chunks)} chunks with parallel_workers={max_workers} "
                    f"(cpus_per_task={int(resources.cpus_per_task)})\n"
                )
                log_stream.flush()

                future_to_chunk = {}
                with ThreadPoolExecutor(max_workers=max_workers) as pool:
                    for chunk in chunks:
                        maf_path = os.path.join(maf_dir, f"{chunk}.maf")
                        log_stream.write(f"RHO SUBMIT: {maf_path}\n")
                        future_to_chunk[pool.submit(estimate_chunk_rho, chunk)] = chunk
                    log_stream.flush()

                    for fut in as_completed(future_to_chunk):
                        chunk_name, maf_path, rc, rho_val = fut.result()
                        if rc != 0:
                            log_stream.write(f"RHO FAIL: phastCons exit {rc} on {maf_path}\n")
                        else:
                            log_stream.write(f"RHO DONE: {chunk_name} rho={rho_val}\n")
                        chunk_rho[chunk_name] = rho_val
                        log_stream.flush()

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

rule run_phastcons_chr:
    input:
        manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.filtered.txt"),
        mod = PHYLOFIT_ACTIVE_MODEL_PATH,
        global_rho = rules.global_rho.output.global_rho
    output:
        chunks_done = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "chunks.done")
    log:
        job_log = os.path.join(LOG_DIR, "run_phastcons_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "run_phastcons_chr", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("phastcons_per_chunk")
    run:
        import math, os, shlex, subprocess, traceback
        from concurrent.futures import ThreadPoolExecutor, as_completed

        with open(log.job_log, "w") as log_stream:
            try:
                chr_dir = os.path.join(CONSERVE_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
                maf_dir = os.path.join(MAF_SPLIT_NS_DIR, wildcards.chromosome_group, wildcards.ref_chromosome)
                chunk_log_dir = os.path.join(
                    LOG_DIR,
                    "run_phastcons_chr",
                    wildcards.chromosome_group,
                    wildcards.ref_chromosome,
                )
                os.makedirs(chr_dir, exist_ok=True)
                os.makedirs(chunk_log_dir, exist_ok=True)

                with open(input.global_rho) as gf:
                    rho_str = gf.read().strip()
                try:
                    rho_global = float(rho_str)
                except Exception:
                    rho_global = float("nan")

                with open(input.manifest) as mf:
                    chunks = [os.path.splitext(line.strip())[0] for line in mf if line.strip()]

                max_workers = max(1, min(int(resources.cpus_per_task), len(chunks) if chunks else 1))
                log_stream.write(
                    f"Applying phastCons to {len(chunks)} chunks with parallel_workers={max_workers} "
                    f"(cpus_per_task={int(resources.cpus_per_task)})\n"
                )
                log_stream.flush()

                def run_chunk(chunk):
                    maf_path = os.path.join(maf_dir, f"{chunk}.maf")
                    bed = os.path.join(chr_dir, f"{chunk}.conserved.bed")
                    wig = os.path.join(chr_dir, f"{chunk}.scores.wig")
                    err = os.path.join(chunk_log_dir, f"{chunk}.stderr.log")

                    if not (math.isfinite(rho_global) and rho_global > 0.0):
                        open(bed, "w").close()
                        open(wig, "w").close()
                        with open(err, "w") as ef:
                            ef.write(f"SKIP: invalid rho '{rho_str}'\n")
                        return chunk, maf_path, "skip-invalid-rho", 0

                    cmd = [
                        "phastCons",
                        maf_path,
                        input.mod,
                        "--rho",
                        str(rho_global),
                        "--most-conserved",
                        bed,
                    ]
                    with open(wig, "w") as wig_stream, open(err, "w") as err_stream:
                        result = subprocess.run(cmd, stdout=wig_stream, stderr=err_stream, text=True)
                    if result.returncode != 0:
                        return chunk, maf_path, "fail", result.returncode
                    return chunk, maf_path, "ok", result.returncode

                failed_chunks = []

                if chunks:
                    if max_workers == 1:
                        for chunk in chunks:
                            maf_path = os.path.join(maf_dir, f"{chunk}.maf")
                            log_stream.write(f"CHUNK: {maf_path}\n")
                            log_stream.flush()
                            chunk_name, maf_path, status, rc = run_chunk(chunk)
                            if status == "fail":
                                log_stream.write(f"CHUNK FAIL: phastCons exit {rc} on {maf_path}\n")
                                failed_chunks.append((chunk_name, maf_path, rc))
                            else:
                                log_stream.write(f"CHUNK DONE: {chunk_name} status={status}\n")
                            log_stream.flush()
                    else:
                        future_to_chunk = {}
                        with ThreadPoolExecutor(max_workers=max_workers) as pool:
                            for chunk in chunks:
                                maf_path = os.path.join(maf_dir, f"{chunk}.maf")
                                log_stream.write(f"CHUNK SUBMIT: {maf_path}\n")
                                future_to_chunk[pool.submit(run_chunk, chunk)] = chunk
                            log_stream.flush()
                            for fut in as_completed(future_to_chunk):
                                chunk_name, maf_path, status, rc = fut.result()
                                if status == "fail":
                                    log_stream.write(f"CHUNK FAIL: phastCons exit {rc} on {maf_path}\n")
                                    failed_chunks.append((chunk_name, maf_path, rc))
                                else:
                                    log_stream.write(f"CHUNK DONE: {chunk_name} status={status}\n")
                                log_stream.flush()

                if failed_chunks:
                    first_chunk, first_maf, first_rc = failed_chunks[0]
                    raise RuntimeError(
                        f"phastCons failed on {len(failed_chunks)} chunk(s); "
                        f"first failure: {first_chunk} ({first_maf}) exit={first_rc}"
                    )

                with open(output.chunks_done, "w") as done:
                    done.write("ok\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule phastcons_concat_chr:
    input:
        manifest = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.filtered.txt"),
        chunks_done = rules.run_phastcons_chr.output.chunks_done
    output:
        chr_bed = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}.bed")
    params:
        rule_name = "phastcons_concat_chr",
        cleanup_chunk_intermediates = CLEANUP_CHUNK_INTERMEDIATES
    log:
        job_log = os.path.join(LOG_DIR, "phastcons_concat_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "phastcons_concat_chr", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("phastcons_concat_chr")
    run:
        import glob
        import os
        import subprocess
        import traceback

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

                sort_cmd = ["sort", "-k1,1", "-k2,2n", "-k3,3n", tmp]
                log_stream.write(f"Running: {' '.join(sort_cmd)} > {output.chr_bed}\n")
                log_stream.flush()

                with open(output.chr_bed, "w") as out_stream:
                    result = subprocess.run(sort_cmd, stdout=out_stream, text=True)
                if result.returncode != 0:
                    raise RuntimeError(f"Sort failed with exit code {result.returncode}")

                try:
                    os.remove(tmp)
                except OSError:
                    pass

                # Keep only final chromosome outputs by cleaning noisy chunk-level intermediates.
                if params.cleanup_chunk_intermediates:
                    cleanup_patterns = [
                        "*.conserved.bed",
                        "*.scores.wig",
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
        job_log = os.path.join(LOG_DIR, "extract_cds_bed_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "extract_cds_bed_chr", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("extract_cds_bed_chr")
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
        conserved_bed = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}.bed"),
        cds_bed = rules.extract_cds_bed_chr.output.cds_bed
    output:
        cnees_bed = os.path.join(CNEES_DIR, "{chromosome_group}", "{ref_chromosome}.cnees.bed"),
        filter_summary = os.path.join(CNEES_SUMMARY_DIR, "{chromosome_group}", "{ref_chromosome}.cnees-filter-summary.tsv")
    log:
        job_log = os.path.join(LOG_DIR, "cnees_from_conserved_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "cnees_from_conserved_chr", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("cnees_from_conserved_chr")
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

                os.makedirs(os.path.dirname(output.filter_summary), exist_ok=True)
                with open(output.filter_summary, "w") as sf:
                    sf.write("metric\tvalue\n")
                    sf.write(f"ces_raw\t{len(conserved_raw)}\n")
                    sf.write(f"ces_merged\t{len(conserved)}\n")
                    sf.write(f"cnees_after_cds_subtract\t{len(out_rows)}\n")
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
        job_log = os.path.join(LOG_DIR, "cnees_to_bed4_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "cnees_to_bed4_chr", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("cnees_to_bed4_chr")
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

rule cnee_alignments_chr:
    input:
        cnees_bed4 = rules.cnees_to_bed4_chr.output.cnees_bed4,
        maf = rules.maf_index_chr.input.maf,
        maf_index_block = rules.maf_index_chr.output.maf_index_block
    output:
        manifest = os.path.join(CNEE_MAF_DIR, "{chromosome_group}", "{ref_chromosome}", "manifest.txt")
    params:
        outdir = os.path.join(CNEE_MAF_DIR, "{chromosome_group}", "{ref_chromosome}"),
        rule_name = "cnee_alignments_chr"
    log:
        job_log = os.path.join(LOG_DIR, "cnee_alignments_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "cnee_alignments_chr", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("cnee_alignments_chr")
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
                    "mafutils", "fetch",
                    input.maf,
                    input.maf_index_block,
                    input.cnees_bed4,
                    "-o", params.outdir,
                    "-p", str(p),
                    "-m", "block",
                    "-b", "id",
                ]
                if CNEE_OUTPUT_FORMAT == "fasta":
                    cmd += ["-f", "-fh", CNEE_FASTA_HEADER]
                    if CNEE_EXPECTED_SPECIES:
                        cmd += ["--expected-species", ",".join(CNEE_EXPECTED_SPECIES)]
                COMMON.runCommand(
                    cmd, log_stream, log_stream, params.rule_name,
                    wc=f"{wildcards.chromosome_group}.{wildcards.ref_chromosome}"
                )

                ext = "fa" if CNEE_OUTPUT_FORMAT == "fasta" else "maf"
                outs = sorted(glob.glob(os.path.join(params.outdir, f"{wildcards.ref_chromosome}.cnee*.{ext}")))

                if CNEE_OUTPUT_FORMAT == "fasta":
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
                log_stream.write(f"Wrote manifest with {len(outs)} CNEE {CNEE_OUTPUT_FORMAT.upper()} files\n")

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

# rule phastcons_chunk moved to workflow/stash.smk - dead code, orphaned since
# run_phastcons_chr inlined this per-chunk phastCons scoring work.
