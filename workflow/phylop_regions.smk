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
MAF_CHR_PREFIX = config["maf_chr_prefix"]
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"]

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

rule run_phylop:
    input:
        adj_mod_file = PHYLOFIT_ACTIVE_MODEL_PATH,
        chromosome_maf = get_chr_maf_for_stage
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

rule naive_cluster_conserved_sites:
    input:
        conserved_sites_bed = rules.get_conserved_sites.output.conserved_sites_bed
    output:
        conserved_regions_bed = os.path.join(PHYLOP_REGIONS_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.bed")
    log:
        job_log = os.path.join(LOG_DIR, "naive_cluster_conserved_sites", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "naive_cluster_conserved_sites", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("naive_cluster_conserved_sites")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.conserved_regions_bed), exist_ok=True)
                sites = []
                with open(input.conserved_sites_bed) as inf:
                    for line in inf:
                        if not line.strip() or line.startswith("#"):
                            continue
                        chrom, start_s, end_s = line.rstrip("\n").split("\t")[:3]
                        sites.append((chrom, int(start_s), int(end_s)))

                clusters, total_regions = INTERVALS.cluster_sites(
                    sites, NAIVE_CLUSTER_MAX_GAP_BP, NAIVE_MIN_REGION_SITES, NAIVE_MIN_REGION_LEN_BP
                )

                kept = len(clusters)
                dropped = total_regions - kept
                with open(output.conserved_regions_bed, "w") as out:
                    for chrom, start, end, count, idx in clusters:
                        out.write(
                            f"{chrom}\t{start}\t{end}\t"
                            f"naive_cluster_{idx:07d}\t{count}\n"
                        )

                log_stream.write(
                    f"naive_merge_gap_bp={NAIVE_CLUSTER_MAX_GAP_BP}; min_sites={NAIVE_MIN_REGION_SITES}; "
                    f"min_len_bp={NAIVE_MIN_REGION_LEN_BP}; kept={kept}; dropped={dropped}\n"
                )
            except Exception:
                traceback.print_exc(file=log_stream)
                raise
