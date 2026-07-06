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

REF_INDEX = config["ref_genome_index"];
REF_CHROMOSOME_GROUPS = config["ref_chromosome_groups"];
# Reference genome info

#############################################################################
# Other params

WINDOW_SIZE = config.get("window_size");

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
        expand(os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", "windows", str(WINDOW_SIZE), MAF_CHR_PREFIX + "{ref_chromosome}-" + str(WINDOW_SIZE) + ".bed"),
                zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)

        # expand(os.path.join(LOG_DIR, "{chromosome_group}" + "-" + "{ref_chromosome}" + "-" + str(WINDOW_SIZE) + ".done"),
        #         zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)

# rule all:
#     input:
#         expand(os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", "windows", MAF_CHR_PREFIX + "{ref_chromosome}-" + str(WINDOW_SIZE) + ".bed"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)

        #expand(os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.bed"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),

        #os.path.join(OUTPUT_DIR, os.path.basename(REF_INDEX) + ".bed")
        #expand(os.path.join(OUTPUT_DIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)

#############################################################################

rule split_group_beds:
    input:
        chr_group_bed = os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}.bed")
    output:
        chr_bed = os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", "scaffolds", MAF_CHR_PREFIX + "{ref_chromosome}.bed")
    params:
        rule_name = "split_group_beds"
    log:
        job_log = os.path.join(LOG_DIR, "split-group-beds", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.log")
    resources:
        **getRuleResources("split_group_beds")
    run:
        with open(log.job_log, "w") as log_stream, \
            open(output.chr_bed, "w") as out_stream:
            try:
                cmd = [ "grep", f"^{wildcards.ref_chromosome}", input.chr_group_bed ];

                COMMON.runCommand(cmd, log_stream, out_stream, params.rule_name, wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

####################

checkpoint generate_window_beds:
    input:
        chr_bed = os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", "scaffolds", MAF_CHR_PREFIX + "{ref_chromosome}.bed")
    output:
        chr_window_bed = os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", "windows", str(WINDOW_SIZE), MAF_CHR_PREFIX + "{ref_chromosome}-" + str(WINDOW_SIZE) + ".bed")
    params:
        window_size = WINDOW_SIZE,
        rule_name = "generate_window_beds"
    log:
        job_log = os.path.join(LOG_DIR, "generate_window_beds", str(WINDOW_SIZE), "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-windows.log")
    resources:
        **getRuleResources("generate_window_beds")
    run:
        with open(log.job_log, "w") as log_stream, \
            open(output.chr_window_bed, "w") as out_stream:
            try:
                cmd = [ "bedtools", "makewindows",
                        "-b", input.chr_bed,
                        "-w", str(params.window_size) ];

                COMMON.runCommand(cmd, log_stream, out_stream, params.rule_name, wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

####################

# checkpoint fetch_window_mafs:
#     input:
#         chr_window_bed = os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}", "windows", str(WINDOW_SIZE), MAF_CHR_PREFIX + "{ref_chromosome}-" + str(WINDOW_SIZE) + ".bed"),
#         maf_file = MAF_PATH,
#         maf_index = MAF_INDEX_BLOCK
#     output:
#         window_mafs_dir = directory(os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}", "windows", str(WINDOW_SIZE), "{ref_chromosome}"))
#     params:
#         script_path = os.path.join(UTILS_DIR, "maf_fetch.py"),
#         rule_name = "fetch_window_mafs"
#     log:
#         job_log = os.path.join(LOG_DIR, "fetch_window_mafs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.log")
#     resources:
#         **getRuleResources("fetch_window_mafs")
#     run:
#         import os
#         os.makedirs(output.window_mafs_dir, exist_ok=True)
        
#         with open(log.job_log, "w") as log_stream:
#             try:
#                 cpus = config.get("rule_resources", {}).get("fetch_window_mafs", {}).get("cpus", 1)
#                 cmd = [
#                     "python", params.script_path,
#                     "--basename", "coords",
#                     "--processes", str(cpus),
#                     "--output", output.window_mafs_dir,
#                     input.maf_file,
#                     input.maf_index,
#                     input.chr_window_bed
#                 ]

#                 # maf_fetch creates files directly in the output directory
#                 COMMON.runCommand(cmd, log_stream, log_stream, params.rule_name,
#                                 wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}")
                                
#             except Exception as e:
#                 traceback.print_exc(file=log_stream)
#                 raise

# ####################

# def getWindowMAFs(wildcards):
#     checkpoint_output = checkpoints.fetch_window_mafs.get(**wildcards)
    
#     bed_file = os.path.join(OUTPUT_DIR, "beds", wildcards.chromosome_group, "windows", 
#                             MAF_CHR_PREFIX + wildcards.ref_chromosome + "-" + str(WINDOW_SIZE) + ".bed")
    
#     window_mafs = []
#     with open(bed_file, 'r') as f:
#         for line in f:
#             if line.strip():
#                 parts = line.strip().split('\t')
#                 chrom, start, end = parts[0], parts[1], parts[2]
#                 window_id = f"{chrom}-{start}-{end}"
                
#                 # Fix the path to match the checkpoint output structure
#                 maf_file = os.path.join(OUTPUT_DIR, "mafs", "windows", str(WINDOW_SIZE),
#                                         wildcards.chromosome_group, 
#                                         wildcards.ref_chromosome, 
#                                         f"{window_id}.maf")
#                 window_mafs.append(maf_file)
    
#     return window_mafs

# rule finalize_chromosome_windows:
#     input:
#         getWindowMAFs
#     output:
#         # Move done flag OUTSIDE the directory output
#         done_flag = os.path.join(LOG_DIR, "{chromosome_group}"+ "-" + "{ref_chromosome}" + "-" + str(WINDOW_SIZE) + ".done")
#     shell:
#         "touch {output.done_flag}"
