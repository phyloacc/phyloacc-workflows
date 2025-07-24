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

PREFIX = config["prefix"];
# The prefix for all output files

getRuleResources = partial(COMMON.getResources, config)
# This maps the function to get rule resources from the config file
# so we don't have to pass config each time we call it

#############################################################################
# Input files and output paths

MOD_DIR = os.path.join(OUTPUT_DIR, "03-phylofit")

MAF_PATH = config["maf"];
MAF_FILE = os.path.basename(MAF_PATH);
MAF_DIR = os.path.dirname(MAF_PATH);
# MAF info

MAF_INDEX_BLOCK = os.path.join(MAF_DIR, MAF_FILE + ".block.idx");
MAF_INDEX_SCAFF = os.path.join(MAF_DIR, MAF_FILE + ".scaffold.idx");
# This will be created

MAF_REF_ID = config["maf_ref_id"];
MAF_CHR_PREFIX = config["maf_chr_prefix"];
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"];
MAF_REF_PREFIX = MAF_REF_ID + MAF_REF_CHR_JOINER + MAF_CHR_PREFIX;
# MAF reference and scaffold id parsing

REF_FASTA = config["ref_fasta"];
REF_FASTA_INDEX = config["ref_fasta_index"];
REF_CHROMOSOME_GROUPS = config["ref_chromosome_groups"];
# Reference genome info

#############################################################################
# Other params

MIN_NS_TO_SPLIT_BY = config["min_Ns_to_split_by"];

flattened_chromosome_groups = [(group, chromosome) for group, chromosome_list in REF_CHROMOSOME_GROUPS.items() for chromosome in chromosome_list];

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

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", "windows", "Ns", MAF_CHR_PREFIX + "{ref_chromosome}-Ns.bed"),
                zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)

#############################################################################

rule picard_scatter_by_ns:
    input:
        ref_fasta = REF_FASTA
    output:
        ref_interval_file = os.path.join(OUTPUT_DIR, os.path.basename(REF_FASTA) + ".N-intervals")
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
        ref_interval_file = os.path.join(OUTPUT_DIR, os.path.basename(REF_FASTA) + ".N-intervals")
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