#############################################################################
# phyloP site calling and naive conserved-region clustering
#############################################################################

import os
import sys
import logging
import traceback

import lib.common as COMMON
from lib.common import spacedOut as SO
import lib.intervals as INTERVALS
# lib.clustering is imported at runtime inside the cluster_conserved_sites run block, not
# here: it pulls in numpy, and the Snakefile must parse for DAG validation in a numpy-free
# environment. Importing it at the point of use keeps parsing dependency-free and still
# fails fast (clear ImportError) if numpy is ever missing when the rule actually runs.

from functools import partial

#############################################################################
# System setup

config_flag = config.get("display", False)
version_flag = config.get("version", False)
info_flag = config.get("info", False)
debug = config.get("debug", False)

if config.get("__master_workflow__", False):
    setup = config["__pipeline_setup__"]
    MAIN = setup["MAIN"]
    DRY_RUN = setup["DRY_RUN"]
    OUTPUT_DIR = setup["OUTPUT_DIR"]
    LOG_DIR = setup["LOG_DIR"]
    TMPDIR = setup["TMPDIR"]
    LOG_LEVEL = setup["LOG_LEVEL"]
    LOG_VERBOSITY = setup["LOG_VERBOSITY"]
else:
    MAIN, DRY_RUN, OUTPUT_DIR, LOG_DIR, TMPDIR, LOG_LEVEL, LOG_VERBOSITY = COMMON.pipelineSetup(
        config, sys.argv, version_flag, info_flag, config_flag, debug, workflow
    )

MLOG = logging.getLogger("META")
MODULE_DIR = os.path.dirname(workflow.snakefile)
PIPELINE_DIR = os.path.dirname(MODULE_DIR) if os.path.basename(MODULE_DIR) == "workflow" else MODULE_DIR
getRuleResources = partial(COMMON.getResources, config)

#############################################################################
# Inputs and outputs

MAF_PREP_DIR = os.path.join(OUTPUT_DIR, "01-maf-prep")
NEUTRAL_MODEL_DIR = os.path.join(OUTPUT_DIR, "02-neutral-model")
PHYLOP_STAGE_DIR = os.path.join(OUTPUT_DIR, "03-phylop")

MAF_SPLIT_BY_CHROM_DIR = COMMON.getOptionalConfigPath(
    config,
    "maf_split_chr_dir",
    os.path.join(MAF_PREP_DIR, "maf-by-chromosome"),
)
PHYLOFIT_DIR = COMMON.getOptionalConfigPath(
    config,
    "phylofit_chr_dir",
    os.path.join(NEUTRAL_MODEL_DIR, "phylofit"),
)

MAF_REF_ID = config["maf_ref_id"]
# Chromosome-name reconciliation. The config lists the CORE chromosome id (both
# prefixes stripped); the MAF and GFF names derive from it. maf_prefix/gff_prefix
# default "". maf_chr_prefix is kept as a legacy alias for maf_prefix.
MAF_PREFIX = str(config.get("maf_prefix", config.get("maf_chr_prefix", "")))
GFF_PREFIX = str(config.get("gff_prefix", ""))
MAF_CHR_PREFIX = MAF_PREFIX  # legacy alias: existing path/filename code uses MAF_CHR_PREFIX
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"]

def maf_chrom(core):
    """Chromosome name as it appears in the MAF (used for all bed contents / mafutils args)."""
    return f"{MAF_PREFIX}{core}"

REF_CHROMOSOME_GROUPS = config["ref_chromosome_groups"]


def get_chr_maf_for_stage(wildcards):
    if "checkpoints" in globals() and hasattr(checkpoints, "maf_split_chr_by_group"):
        cp_output = checkpoints.maf_split_chr_by_group.get(chromosome_group=wildcards.chromosome_group).output
        group_dir = os.path.join(MAF_SPLIT_BY_CHROM_DIR, wildcards.chromosome_group)
        manifest_file = cp_output.maf_manifest
        with open(manifest_file, "r") as mf:
            expected = [line.strip() for line in mf if line.strip()]
        target = f"{MAF_CHR_PREFIX}{wildcards.ref_chromosome}.maf"
        for fname in expected:
            if fname == target:
                return os.path.join(group_dir, fname)
        raise ValueError(f"No chromosome MAF for {wildcards.ref_chromosome} listed in {manifest_file}")
    return os.path.join(MAF_SPLIT_BY_CHROM_DIR, wildcards.chromosome_group, f"{MAF_CHR_PREFIX}{wildcards.ref_chromosome}.maf")

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

REF_CHR_GROUPS_LIST, REF_CHROMOSOMES = zip(*flattened_chromosome_groups)

PHYLOP_ALPHA = str(config.get("phylop_alpha", config.get("alpha", "0.05")))
NAIVE_CLUSTER_MAX_GAP_BP = int(config.get("naive_merge_gap_bp", config.get("naive_cluster_max_gap_bp", 20)))
NAIVE_MIN_REGION_SITES = int(config.get("naive_min_region_sites", 5))
NAIVE_MIN_REGION_LEN_BP = int(config.get("naive_min_region_len_bp", 20))


def _pbool(x, default=False):
    if x is None:
        return default
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y", "on"}


# --- pre-flight phyloP statistical-power gate ---
# Per-site phyloP has a hard detection ceiling set by total neutral tree length:
# on a shallow tree even a perfectly conserved site cannot clear a per-chromosome
# FDR threshold, so phyloP returns ~0 conserved sites by construction. This gate
# computes that ceiling from each chromosome's fitted .mod (lib.phylop_power) and
# fails early (before the expensive phyloP scan) when it can't clear the FDR bar,
# unless phylop_power_override is set. See analyses/phylop-tree-length-power.
PHYLOP_POWER_GATE = _pbool(config.get("phylop_power_gate", True), True)
PHYLOP_POWER_OVERRIDE = _pbool(config.get("phylop_power_override", False), False)
# Number of sites tested (M) for the FDR bar log10(M/alpha). "estimate" (default)
# sums the reference block lengths in the MAF block index for this chromosome - the
# actual number of aligned/scored bases, which equals the chromosome length for a
# whole-chromosome MAF but correctly reports the true (smaller) count for a sliced
# MAF. Falls back to the per-chromosome MAF srcSize (chromosome length), then
# ref_fasta's .fai, then a constant. An integer sets M directly.
PHYLOP_POWER_NUM_SITES = config.get("phylop_power_num_sites", "estimate")
if isinstance(PHYLOP_POWER_NUM_SITES, bool) or not (
    (isinstance(PHYLOP_POWER_NUM_SITES, int))
    or (isinstance(PHYLOP_POWER_NUM_SITES, str)
        and (PHYLOP_POWER_NUM_SITES.strip().lower() == "estimate" or PHYLOP_POWER_NUM_SITES.strip().isdigit()))
):
    raise ValueError(
        f"Invalid phylop_power_num_sites '{PHYLOP_POWER_NUM_SITES}'. Use 'estimate' or a positive integer."
    )
if isinstance(PHYLOP_POWER_NUM_SITES, int) and PHYLOP_POWER_NUM_SITES <= 0:
    raise ValueError("phylop_power_num_sites must be > 0.")
PHYLOP_POWER_FALLBACK_M = int(config.get("phylop_power_fallback_num_sites", 100_000_000))
_POWER_REF_FASTA = str(config.get("ref_fasta") or "").strip()
_POWER_REF_FAI = str(config.get("ref_fasta_index") or "").strip() or (
    _POWER_REF_FASTA + ".fai" if _POWER_REF_FASTA else "")
# Whole-MAF block index (maf_index writes it as <maf>.block.idx). Built early in every
# run - transitively upstream of the gate via maf_split/phyloFit - so it is present by
# gate time; read opportunistically (fall back to srcSize if absent).
_POWER_MAF_BLOCK_IDX = (str(config.get("maf")) + ".block.idx") if config.get("maf") else ""


def _sum_reflen_from_block_index(chrom):
    """Sum of reference block lengths (aligned/scored bases) for `chrom` from the MAF
    block index, or None if unavailable. The index names chromosomes as they appear in
    the MAF (with maf_chr_prefix), so match against the prefixed name."""
    if not _POWER_MAF_BLOCK_IDX or not os.path.exists(_POWER_MAF_BLOCK_IDX):
        return None
    target = f"{MAF_CHR_PREFIX}{chrom}"
    total = 0
    found = False
    try:
        with open(_POWER_MAF_BLOCK_IDX) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                parts = line.split("\t")
                if len(parts) >= 3 and parts[0] == target:
                    try:
                        total += int(parts[2])
                        found = True
                    except ValueError:
                        pass
    except OSError:
        return None
    return total if found else None


def _chrom_len_from_fai(chrom):
    """Reference length of `chrom` from a .fai index, or None if unavailable."""
    if not _POWER_REF_FAI or not os.path.exists(_POWER_REF_FAI):
        return None
    with open(_POWER_REF_FAI) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] == chrom:
                try:
                    return int(parts[1])
                except ValueError:
                    return None
    return None


def _srcsize_from_maf(maf_path):
    """Reference chromosome length from the MAF's first reference s-line (srcSize,
    field 6). O(1): reads only up to the reference sequence line of the first block."""
    ref_prefix = f"{MAF_REF_ID}{MAF_REF_CHR_JOINER}"
    try:
        with open(maf_path) as fh:
            for i, line in enumerate(fh):
                if line[:2] in ("s ", "s\t"):
                    parts = line.split()
                    if len(parts) >= 6 and parts[1].startswith(ref_prefix):
                        try:
                            return int(parts[5])
                        except ValueError:
                            return None
                if i > 100000:  # the reference appears in the first block; avoid scanning a huge MAF on misconfig
                    break
    except OSError:
        return None
    return None


def _resolve_power_m_sites(chrom, maf_path=None):
    v = PHYLOP_POWER_NUM_SITES
    if isinstance(v, int) and not isinstance(v, bool):
        m = v
    elif isinstance(v, str) and v.strip().isdigit():
        m = int(v.strip())
    elif isinstance(v, str) and v.strip().lower() == "estimate":
        # aligned/scored bases: block-index sum(reflen) -> MAF srcSize -> .fai -> constant.
        # sum(reflen) is the true count (correct for sliced MAFs, where srcSize overcounts).
        m = (_sum_reflen_from_block_index(chrom)
             or (_srcsize_from_maf(maf_path) if maf_path else None)
             or _chrom_len_from_fai(chrom) or PHYLOP_POWER_FALLBACK_M)
    else:
        raise ValueError(f"Invalid phylop_power_num_sites '{v}'. Use 'estimate' or a positive integer.")
    if m <= 0:
        raise ValueError(f"phylop_power_num_sites resolved to {m}; must be > 0.")
    return m

# --- conserved-region clustering method (turns significant conserved sites into regions) ---
# Default 'hmm' (best phastCons concordance in analyses/phylop-clustering-comparison). Options:
# 'hmm', 'windowed', 'gap_merge'. 'hdbscan' is implemented in lib/clustering.py but blocked here.
PHYLOP_CLUSTER_METHOD = str(config.get("phylop_cluster_method", "hmm")).strip().lower()

if PHYLOP_CLUSTER_METHOD == "hdbscan":
    raise ValueError(
        "phylop_cluster_method: 'hdbscan' is not supported. On real 241-mammal data it over-calls "
        "conserved sequence badly (~569 kb of false-positive conserved bp vs phastCons's 182 kb; "
        "precision ~0.22), and it can't be tuned out of it. Use 'hmm' (default), 'windowed', or "
        "'gap_merge'. See analyses/phylop-clustering-comparison for the comparison."
    )
_VALID_CLUSTER_METHODS = ("hmm", "windowed", "gap_merge")
if PHYLOP_CLUSTER_METHOD not in _VALID_CLUSTER_METHODS:
    raise ValueError(
        f"Invalid phylop_cluster_method '{PHYLOP_CLUSTER_METHOD}'. "
        f"One of: {', '.join(_VALID_CLUSTER_METHODS)}."
    )

# windowed params
WINDOWED_WINDOW_BP = int(config.get("windowed_window_bp", 20))
WINDOWED_MIN_SITES = int(config.get("windowed_min_sites_per_window", 5))
if PHYLOP_CLUSTER_METHOD == "windowed":
    if WINDOWED_WINDOW_BP <= 0:
        raise ValueError("windowed_window_bp must be > 0")
    if WINDOWED_MIN_SITES < 1:
        raise ValueError("windowed_min_sites_per_window must be >= 1")

# hmm params (2-state online HMM; defaults from analyses/phylop-clustering-comparison)
HMM_T0_0 = float(config.get("hmm_t0_0", 0.9))
HMM_T1_1 = float(config.get("hmm_t1_1", 0.99))
HMM_E0_0 = float(config.get("hmm_e0_0", 0.8))
HMM_E1_1 = float(config.get("hmm_e1_1", 0.5))
HMM_S0 = float(config.get("hmm_s0", 0.9))
HMM_MIN_LEN = int(config.get("hmm_min_len", 20))
HMM_MAX_LEN = int(config.get("hmm_max_len", 100000))
if PHYLOP_CLUSTER_METHOD == "hmm":
    for _pname, _pval in (("hmm_t0_0", HMM_T0_0), ("hmm_t1_1", HMM_T1_1), ("hmm_e0_0", HMM_E0_0),
                          ("hmm_e1_1", HMM_E1_1), ("hmm_s0", HMM_S0)):
        if not (0.0 < _pval < 1.0):
            raise ValueError(f"{_pname} must be strictly between 0 and 1; got {_pval}")
    if not (0 <= HMM_MIN_LEN < HMM_MAX_LEN):
        raise ValueError(f"require 0 <= hmm_min_len < hmm_max_len; got {HMM_MIN_LEN}, {HMM_MAX_LEN}")

PHYLOP_SITES_DIR = os.path.join(PHYLOP_STAGE_DIR, "sites")
PHYLOP_REGIONS_DIR = os.path.join(PHYLOP_STAGE_DIR, "regions")
PHYLOP_SUMMARY_DIR = os.path.join(PHYLOP_STAGE_DIR, "summary")

if LOG_LEVEL == "debug":
    for i in range(len(REF_CHR_GROUPS_LIST)):
        MLOG.debug(SO(f"CHR GROUP {i}", 15) + f"{REF_CHR_GROUPS_LIST[i]} : {REF_CHROMOSOMES[i]}")
    MLOG.debug("EXITING BEFORE RULES. DEBUG MODE.")
    sys.exit(0)

PHYLOP_STANDALONE = not bool(config.get("__master_workflow__", False))

if PHYLOP_STANDALONE:
    localrules: all

    rule all:
        input:
            expand(
                os.path.join(PHYLOP_SUMMARY_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.conserved-site-counts." + PHYLOP_ALPHA + ".tsv"),
                zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES
            ),
            expand(
                os.path.join(PHYLOP_SUMMARY_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.accelerated-site-counts." + PHYLOP_ALPHA + ".tsv"),
                zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES
            ),
            expand(
                os.path.join(PHYLOP_REGIONS_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.bed"),
                zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES
            )

#############################################################################

rule phylop_power_check:
    input:
        mod = PHYLOFIT_ACTIVE_MODEL_PATH,
        chromosome_maf = get_chr_maf_for_stage
    output:
        report = os.path.join(PHYLOP_STAGE_DIR, "power-check", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.power.tsv")
    log:
        job_log = os.path.join(LOG_DIR, "phylop_power_check", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "phylop_power_check", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("phylop_power_check")
    run:
        import lib.phylop_power as PP
        try:
            from snakemake.exceptions import WorkflowError
        except Exception:
            WorkflowError = RuntimeError

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.report), exist_ok=True)
                m_sites = _resolve_power_m_sites(wildcards.ref_chromosome, input.chromosome_maf)
                alpha = float(PHYLOP_ALPHA)
                g = PP.power_gate(input.mod, m_sites, alpha)

                with open(output.report, "w") as rep:
                    rep.write("metric\tvalue\n")
                    for k in ("tree_length", "tips", "ceiling_neglog10p", "threshold_neglog10p",
                              "m_sites", "alpha", "passes"):
                        rep.write(f"{k}\t{g[k]}\n")
                    rep.write("check\tceiling_neglog10p >= log10(m_sites / alpha)\n")
                    rep.write(f"gate_enabled\t{PHYLOP_POWER_GATE}\n")
                    rep.write(f"override\t{PHYLOP_POWER_OVERRIDE}\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

            log_stream.write(
                f"phyloP power check {wildcards.ref_chromosome}: tree_length={g['tree_length']:.2f} subs/site; "
                f"ceiling(-log10 p)={g['ceiling_neglog10p']:.2f} vs log10(M/alpha)="
                f"log10({g['m_sites']}/{alpha})={g['threshold_neglog10p']:.2f}; passes={g['passes']}\n"
            )
            gate_failed = PHYLOP_POWER_GATE and not g["passes"]
            fail_msg = None
            if gate_failed:
                fail_msg = (
                    f"phyloP power gate FAILED for '{wildcards.ref_chromosome}'.\n"
                    f"  gate passes iff  ceiling(-log10 p)  >=  log10(M / alpha):\n"
                    f"    LHS  ceiling = {g['ceiling_neglog10p']:.2f}   <- from neutral tree length "
                    f"{g['tree_length']:.2f} subs/site ({g['tips']} tips)\n"
                    f"    RHS  FDR bar = log10({g['m_sites']} / {alpha}) = {g['threshold_neglog10p']:.2f}\n"
                    f"  {g['ceiling_neglog10p']:.2f} < {g['threshold_neglog10p']:.2f}  =>  FAIL "
                    f"(tree too shallow for per-site phyloP).\n"
                    f"  Set 'phylop_power_override: true' to run anyway, or use phastCons on shallow "
                    f"trees (see analyses/phylop-tree-length-power)."
                )
                # Write the full explanation to the log too (not only the console error),
                # so it's captured per-rule (e.g. under SLURM the console error lands elsewhere).
                log_stream.write(fail_msg + "\n")
                if PHYLOP_POWER_OVERRIDE:
                    log_stream.write("phylop_power_override is set: proceeding despite the failed gate.\n")

        # Raised outside the try/except (after the log is closed) so a failed gate is a clean
        # WorkflowError (no Python traceback); same text as written to the log above.
        if gate_failed and not PHYLOP_POWER_OVERRIDE:
            raise WorkflowError(fail_msg)

####################

rule run_phylop:
    input:
        adj_mod_file = PHYLOFIT_ACTIVE_MODEL_PATH,
        chromosome_maf = get_chr_maf_for_stage,
        power_ok = rules.phylop_power_check.output.report
    output:
        phylop_wig_file = os.path.join(PHYLOP_SITES_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.wig")
    log:
        job_log = os.path.join(LOG_DIR, "run_phylop", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "run_phylop", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("run_phylop")
    run:
        with open(log.job_log, "w") as log_stream, open(output.phylop_wig_file, "w") as out_stream:
            try:
                os.makedirs(os.path.dirname(output.phylop_wig_file), exist_ok=True)
                cmd = [
                    "phyloP", "--method", "LRT",
                    "--mode", "CONACC",
                    "--wig-scores",
                    "-i", "MAF",
                    input.adj_mod_file,
                    input.chromosome_maf,
                ]
                COMMON.runCommand(cmd, log_stream, out_stream, "run_phylop", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule convert_wig_to_bed:
    input:
        phylop_wig_file = rules.run_phylop.output.phylop_wig_file
    output:
        phylop_bed_file = os.path.join(PHYLOP_SITES_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.bed")
    params:
        script_path = os.path.join(PIPELINE_DIR, "utils", "convert_wig_to_bed.awk")
    log:
        job_log = os.path.join(LOG_DIR, "convert_wig_to_bed", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "convert_wig_to_bed", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("convert_wig_to_bed")
    run:
        with open(log.job_log, "w") as log_stream, open(output.phylop_bed_file, "w") as out_stream:
            try:
                os.makedirs(os.path.dirname(output.phylop_bed_file), exist_ok=True)
                cmd = ["awk", "-f", params.script_path, input.phylop_wig_file]
                COMMON.runCommand(cmd, log_stream, out_stream, "convert_wig_to_bed", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule adjust_pvals:
    input:
        phylop_bed_file = rules.convert_wig_to_bed.output.phylop_bed_file
    output:
        phylop_bed_file_fdr = os.path.join(PHYLOP_SITES_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + PHYLOP_ALPHA + ".bed")
    params:
        alpha = PHYLOP_ALPHA,
        tmp_dir = TMPDIR,
        script_path = os.path.join(PIPELINE_DIR, "utils", "adjust_pvals.sh")
    log:
        job_log = os.path.join(LOG_DIR, "adjust_pvals", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "adjust_pvals", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("adjust_pvals")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = [
                    "bash", params.script_path,
                    input.phylop_bed_file,
                    output.phylop_bed_file_fdr,
                    params.tmp_dir,
                    str(resources.cpus_per_task),
                    str(params.alpha),
                ]
                COMMON.runCommand(cmd, log_stream, log_stream, "adjust_pvals", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule get_conserved_sites:
    input:
        phylop_bed_file_fdr = rules.adjust_pvals.output.phylop_bed_file_fdr
    output:
        conserved_sites_bed = os.path.join(PHYLOP_SITES_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + PHYLOP_ALPHA + ".conserved.bed")
    params:
        script_path = os.path.join(PIPELINE_DIR, "utils", "get_conserved_sites.awk")
    log:
        job_log = os.path.join(LOG_DIR, "get_conserved_sites", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "get_conserved_sites", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("get_conserved_sites")
    run:
        with open(log.job_log, "w") as log_stream, open(output.conserved_sites_bed, "w") as out_stream:
            try:
                cmd = ["awk", "-f", params.script_path, input.phylop_bed_file_fdr]
                COMMON.runCommand(cmd, log_stream, out_stream, "get_conserved_sites", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule get_accelerated_sites:
    input:
        phylop_bed_file_fdr = rules.adjust_pvals.output.phylop_bed_file_fdr
    output:
        accelerated_sites_bed = os.path.join(PHYLOP_SITES_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + PHYLOP_ALPHA + ".accelerated.bed")
    params:
        script_path = os.path.join(PIPELINE_DIR, "utils", "get_accelerated_sites.awk")
    log:
        job_log = os.path.join(LOG_DIR, "get_accelerated_sites", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "get_accelerated_sites", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("get_accelerated_sites")
    run:
        with open(log.job_log, "w") as log_stream, open(output.accelerated_sites_bed, "w") as out_stream:
            try:
                cmd = ["awk", "-f", params.script_path, input.phylop_bed_file_fdr]
                COMMON.runCommand(cmd, log_stream, out_stream, "get_accelerated_sites", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule site_counts:
    input:
        conserved_sites_bed = rules.get_conserved_sites.output.conserved_sites_bed,
        accelerated_sites_bed = rules.get_accelerated_sites.output.accelerated_sites_bed
    output:
        conserved_site_counts = os.path.join(PHYLOP_SUMMARY_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.conserved-site-counts." + PHYLOP_ALPHA + ".tsv"),
        accelerated_site_counts = os.path.join(PHYLOP_SUMMARY_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.accelerated-site-counts." + PHYLOP_ALPHA + ".tsv")
    log:
        job_log = os.path.join(LOG_DIR, "site_counts", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "site_counts", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("site_counts")
    run:
        os.makedirs(os.path.dirname(output.conserved_site_counts), exist_ok=True)
        with open(log.job_log, "w") as log_stream, open(output.conserved_site_counts, "w") as out_stream:
            try:
                cnt = sum(1 for _ in open(input.conserved_sites_bed))
                out_stream.write(f"{wildcards.ref_chromosome}\t{cnt}\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

        with open(log.job_log, "a") as log_stream, open(output.accelerated_site_counts, "w") as out_stream:
            try:
                cnt = sum(1 for _ in open(input.accelerated_sites_bed))
                out_stream.write(f"{wildcards.ref_chromosome}\t{cnt}\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule cluster_conserved_sites:
    input:
        conserved_sites_bed = rules.get_conserved_sites.output.conserved_sites_bed
    output:
        conserved_regions_bed = os.path.join(PHYLOP_REGIONS_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.bed")
    log:
        job_log = os.path.join(LOG_DIR, "cluster_conserved_sites", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "cluster_conserved_sites", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("cluster_conserved_sites")
    run:
        # Cluster FDR-significant conserved sites into regions using phylop_cluster_method.
        # Output BED4+count: chrom, start, end, <method>_cluster_NNNNNNN, n_sites.
        with open(log.job_log, "w") as log_stream:
            try:
                import lib.clustering as CLUSTER  # runtime-only (numpy-backed); see import note near top
                os.makedirs(os.path.dirname(output.conserved_regions_bed), exist_ok=True)
                sites = []
                with open(input.conserved_sites_bed) as inf:
                    for line in inf:
                        if not line.strip() or line.startswith("#"):
                            continue
                        _c, start_s, end_s = line.rstrip("\n").split("\t")[:3]
                        sites.append((int(start_s), int(end_s)))

                chrom = maf_chrom(wildcards.ref_chromosome)  # regions feed mafutils; use the MAF name
                # length = last conserved-site end; there are no conserved regions past it, and
                # windowed/hmm only need to scan up to there (no chromosome-length lookup needed).
                length = max((e for _s, e in sites), default=0)

                if PHYLOP_CLUSTER_METHOD == "gap_merge":
                    clusters, total = INTERVALS.cluster_sites(
                        [(chrom, s, e) for s, e in sites],
                        NAIVE_CLUSTER_MAX_GAP_BP, NAIVE_MIN_REGION_SITES, NAIVE_MIN_REGION_LEN_BP)
                    regions = [(s, e, n) for (_c, s, e, n, _idx) in clusters]
                    params_desc = (f"gap={NAIVE_CLUSTER_MAX_GAP_BP}, min_sites={NAIVE_MIN_REGION_SITES}, "
                                   f"min_len={NAIVE_MIN_REGION_LEN_BP}, pre_filter_regions={total}")
                elif PHYLOP_CLUSTER_METHOD == "windowed":
                    regions = CLUSTER.windowed(sites, length, WINDOWED_WINDOW_BP, WINDOWED_MIN_SITES)
                    params_desc = f"window_bp={WINDOWED_WINDOW_BP}, min_sites_per_window={WINDOWED_MIN_SITES}"
                elif PHYLOP_CLUSTER_METHOD == "hmm":
                    regions = CLUSTER.hmm(sites, length, t0_0=HMM_T0_0, t1_1=HMM_T1_1, e0_0=HMM_E0_0,
                                          e1_1=HMM_E1_1, s0=HMM_S0, min_len=HMM_MIN_LEN, max_len=HMM_MAX_LEN)
                    params_desc = (f"t0_0={HMM_T0_0}, t1_1={HMM_T1_1}, e0_0={HMM_E0_0}, e1_1={HMM_E1_1}, "
                                   f"s0={HMM_S0}, min_len={HMM_MIN_LEN}, max_len={HMM_MAX_LEN}")
                else:  # already validated at parse time, but guard anyway
                    raise ValueError(f"Unhandled phylop_cluster_method '{PHYLOP_CLUSTER_METHOD}'")

                with open(output.conserved_regions_bed, "w") as out:
                    for idx, (start, end, count) in enumerate(regions, start=1):
                        out.write(f"{chrom}\t{start}\t{end}\t{PHYLOP_CLUSTER_METHOD}_cluster_{idx:07d}\t{count}\n")

                log_stream.write(
                    f"method={PHYLOP_CLUSTER_METHOD}; {params_desc}; "
                    f"sites={len(sites)}; regions={len(regions)}\n"
                )
                if not sites:
                    log_stream.write(
                        "No FDR-significant conserved sites on this chromosome; wrote an empty region "
                        "set. Downstream CNEE building skips gracefully. If unexpected, check the "
                        "phyloP power gate / tree length (see analyses/phylop-tree-length-power).\n"
                    )
            except Exception:
                traceback.print_exc(file=log_stream)
                raise
