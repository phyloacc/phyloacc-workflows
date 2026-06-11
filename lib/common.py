#############################################################################
# Functions to help with the phyloacc snakemake pipeline
#
# Gregg Thomas, April 2025
#############################################################################

import sys
import os
import shutil
import re
import requests
import subprocess
import logging
from datetime import datetime
import yaml
import contextlib
import traceback

#############################################################################

meta_logger = logging.getLogger('META')
# Get the logger for the cactuslib module

#############################################################################

# Example of your ColoredFormatter
class ColoredFormatter(logging.Formatter):
    RESET = "\033[0m"
    COLOR_MAP = {
        logging.DEBUG: "\033[35m",      # Purple (Magenta) for DEBUG
        logging.INFO: "\033[36m",       # Cyan for INFO
        logging.WARNING: "\033[33m",    # Yellow for WARNING
        logging.ERROR: "\033[31m",      # Red for ERROR
        logging.CRITICAL: "\033[1;31m"    # Bold red for CRITICAL
    }

    def format(self, record):
        color = self.COLOR_MAP.get(record.levelno, self.RESET)
        message = super().format(record)
        return f"{color}{message}{self.RESET}"

# Define the filter that will block records flagged as file_only.
def no_file_only(record):
    return not getattr(record, 'file_only', False)

# Assume meta_logger is created earlier:
meta_logger = logging.getLogger('META')

def configureLogging(log_filename: str, log_level: str, log_verbosity: str) -> None:
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    # Check if logger has handlers already to avoid duplicate handlers
    if not meta_logger.hasHandlers():
        if log_verbosity in ["BOTH", "FILE"]:
            handler_file = logging.FileHandler(log_filename)
            handler_file.setFormatter(logging.Formatter(fmt=log_format, datefmt=date_format))
            meta_logger.addHandler(handler_file)
        # Add file handler if specified

        if log_verbosity in ["BOTH", "SCREEN"]:
            # Create a handler for DEBUG and INFO messages that prints to stdout
            handler_stdout = logging.StreamHandler(sys.stdout)
            handler_stdout.setFormatter(ColoredFormatter(fmt=log_format, datefmt=date_format))
            # Only allow messages with level less than WARNING (DEBUG, INFO)
            handler_stdout.addFilter(lambda record: record.levelno < logging.WARNING)
            # And add our extra filter to skip file-only records
            handler_stdout.addFilter(no_file_only)
            meta_logger.addHandler(handler_stdout)

            # Create a separate handler for WARNING and above that prints to stderr
            handler_stderr = logging.StreamHandler(sys.stderr)
            handler_stderr.setFormatter(ColoredFormatter(fmt=log_format, datefmt=date_format))
            # Only allow messages with level WARNING or higher
            handler_stderr.addFilter(lambda record: record.levelno >= logging.WARNING)
            # Also add our filter here to skip file-only records
            handler_stderr.addFilter(no_file_only)
            meta_logger.addHandler(handler_stderr)

        if log_verbosity not in ["BOTH", "FILE", "SCREEN"]:
            raise ValueError("Invalid log verbosity: " + log_verbosity + ". Choose from BOTH, FILE, SCREEN.")
    
    # Set the logging level
    log_level = log_level.upper()
    level_map = {
        'CRITICAL': logging.CRITICAL,
        'ERROR': logging.ERROR,
        'WARNING': logging.WARNING,
        'INFO': logging.INFO,
        'DEBUG': logging.DEBUG,
        'NOTSET': logging.NOTSET
    }
    try:
        meta_logger.setLevel(level_map[log_level])
        meta_logger.debug(f"Logging level set to {log_level}")
    except KeyError:
        raise ValueError(f"Invalid logging level: {log_level}. Choose from {list(level_map.keys())}")


#############################################################################

def getInfo(version_flag, info_flag, args):
    color = "\033[36m"; # cyan
    reset_color = "\033[0m";
    # Some color codes for printing

    info_path = os.path.join(os.path.dirname(__file__), "info.yaml")
    with open(info_path, "r") as file:
        info = yaml.safe_load(file);
    # Read the meta info from the info.yaml file

    if version_flag:
        print(f"\n{color}Snakemake conserved elements pipeline version {info['version']} released on {info['releasedate-patch']}{reset_color}");
        sys.exit();
    # If the version flag is set, print the version and exit

    if info_flag:
        snakefile = "";
        if '-s' in args:
            idx = args.index('-s');
            if idx + 1 < len(args):
                snakefile = args[idx + 1];
        elif '--snakefile' in args:
            idx = args.index('--snakefile');
            if idx + 1 < len(args):
                snakefile = args[idx + 1];
        # Try to get the path to the snakefile

        snakefile_path = os.path.abspath(snakefile);
        mod_timestamp = os.path.getmtime(snakefile_path);
        mod_datetime = datetime.fromtimestamp(mod_timestamp);
        # Get the last modified date of the snakefile

        print(f"\n{color}---Snakemake conserved elements pipeline---{reset_color}")
        for key, value in info.items():
            if value:
                print(f"{color}{key}: {value}{reset_color}");

                if key == "latest-commit-date" and snakefile:
                    print(f"{color}{os.path.basename(snakefile)} last modified: {mod_datetime.date()}{reset_color}");
                # Print the last modified date of the snakefile if it exists
        sys.exit();
    # If the info flag is set, print the meta info and exit

#############################################################################

def pipelineSetup(config, args, version_flag, info_flag, config_flag, debug, workflow):
    main_flag = True;
    if "__main__.py" in args[0]:
        main_flag = False;
    # Whether the pipeline is being run as a main script or not 

    if main_flag and version_flag or info_flag:
        info = getInfo(version_flag, info_flag, args);
    # Print version or info if specified, which will terminate the program early

    dry_run_flag = False;
    if any([arg in args for arg in ["--dry-run", "--dryrun", "-n"]]):
        dry_run_flag = True;
    # Whether the pipeline is running in dry-run mode

    log_level = "info";
    if any([arg in args for arg in ["--rulegraph", "--dag"]]):
        log_level = "notset";
    if debug or config_flag:
        log_level = "debug";
    # Set the log level based on the arguments

    #output_dir_cfg = os.path.abspath(config["output_dir"]);
    output_dir = os.path.abspath(config["output_dir"]);
    # if dry_run_flag:
    #     import atexit, shutil
    #     output_dir = os.path.join("/tmp/", "cactus-smk-dryrun");
    #     atexit.register(lambda: shutil.rmtree(output_dir, ignore_errors=True));
    log_dir = os.path.join(output_dir, "logs");
    # The output directory where all the files and logs are stored

    if main_flag and not version_flag or not info_flag:
        outdir_log_msg, outdir_err_flag = createOutputDirs(output_dir, log_dir, dry_run_flag);
    # Create the output directories if they don't exist

    log_verbosity = "both"; # "screen", "file", "both"
    if dry_run_flag:
        log_verbosity = "screen";
    # Set the log verbosity based on the arguments

    log_filename = os.path.join(log_dir, f"{log_level}.log"); # Log file name if log_verbosity is "file" or "both"
    configureLogging(log_filename, log_level.upper(), log_verbosity.upper());
    meta_logger = logging.getLogger('META')
    # Set up the logger

    if config_flag or debug:
        pad = 30;
        meta_logger.debug(spacedOut("Config file", pad) + os.path.abspath(workflow.configfiles[0])); 
        meta_logger.debug("---");   
        for key, value in config.items():
            meta_logger.debug(spacedOut(key, pad) + str(value));
        meta_logger.debug("=" * 80);
        if config_flag:
            sys.exit();

    if main_flag:
        meta_logger.info(f"MAIN call: {' '.join(args)}");
        # Log the command that was run to start the pipeline

        if outdir_err_flag:
            meta_logger.error(outdir_log_msg);
            sys.exit(1);
        else:
            meta_logger.info(outdir_log_msg);
        # Log the command that was run to create the output directory
        target_jobs = None;
    else:
        target_jobs = None
        if "--target-jobs" in args:
            target_jobs_index = args.index("--target-jobs")
            if target_jobs_index + 1 < len(args):
                target_jobs = args[target_jobs_index + 1];
                if target_jobs[-1] == ":":
                    target_jobs = target_jobs[:-1];
        meta_logger.info(f"RULE {target_jobs} call: {' '.join(args)}", extra={'file_only': True});
        # Log the command that was run for a rule

    tmp_dir = config["tmp_dir"];
    if not os.path.exists(tmp_dir):
        if main_flag:
            meta_logger.info(f"Creating temporary directory at {tmp_dir}");
        os.makedirs(tmp_dir);
    # A directory with lots of space to use for temporary files generated by the cactus-align command

    return main_flag, dry_run_flag, output_dir, log_dir, tmp_dir, log_level, log_verbosity

#############################################################################

def createOutputDirs(outdir, logdir, dry_run, outdir_cfg=""):
    err_flag = False;

    outdir_exists = os.path.exists(outdir);

    # if dry_run:
    #     msg = f"Output directory {outdir_cfg} will be created.";
    # else:
    if not outdir_exists:
        os.makedirs(outdir);
        msg = f"Created output directory: {outdir}";
    # elif outdir_exists and not overwite_output_dir:
    #     msg = f"Output directory already exists: {outdir}. Remove the directory or set overwrite_output_dir to True in your config file to use it anyway and potentially overwrite files from previous runs.";
    #     err_flag = True;
    else:
        msg = f"Output directory {outdir} already exists. Continuing.";
    # Logging for the output directory        
    # Make the output directory if it doesn't exist

    if not os.path.exists(logdir):
        os.makedirs(logdir);
    # Make the log directory if it doesn't exist. This has to be done before the rest
    # so the logging works

    return msg, err_flag;

#############################################################################

def spacedOut(string, totlen, sep="."):
# Properly adds spaces to the end of a message to make it a given length
    spaces = sep * (totlen - len(string));
    if len(string) > totlen:
        spaces += sep * 4;
    return string + spaces;

def printWrite(string, stream):
    if "- INFO -" in string:
        color = "\033[36m"; # cyan
    elif "- ERROR -" in string:
        color = "\033[31m"; # red
    reset_color = "\033[0m";
    color_string = color + string + reset_color;
    # Format the string with color codes

    print(color_string, flush=True);
    stream.write(string + "\n");
    stream.flush();
# For logging in runCommand(), print the string and write it to the file stream

def writeFlush(string, stream):
    stream.write(string + "\n");
    stream.flush();
# For logging in runCommand(), write the string to the file stream

#############################################################################

def fmtDT():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')

def fmtDTLog():
    return datetime.now().strftime('%Y%m%d%H%M%S')

#############################################################################

def writeBedManifest(bed_file, manifest_file):
        expected_files = []  # Build this list from the input BED file for the group.
        # For example, read the BED file and determine for each line the expected output filename.
        with open(bed_file, "r") as bed:
            for line in bed:
                fields = line.strip().split()

                expected_file = fields[3] + ".maf";
                # Assuming the 4th field is the name of the file to be created

                expected_files.append(expected_file);

        # Write the manifest:
        with open(manifest_file, "w") as mf:
            for fname in expected_files:
                mf.write(fname + "\n");

#############################################################################

def hasSnakemakeOption(option, args=None):
    # Detect a Snakemake CLI option in either "--flag value" or "--flag=value" form.

    if args is None:
        args = sys.argv

    return any(arg == option or arg.startswith(option + "=") for arg in args)

def hasAutoPartitionSelection(args=None):
    # The slurm executor plugin can auto-select a partition when this option is provided.

    return hasSnakemakeOption("--slurm-partition-config", args)

#############################################################################

def getResources(config, rule_name, keys=("partition", "mem_mb", "cpus", "time")):
# Return dict of all requested resource keys for a rule (with fallback to defaults).
    
    slurm_resource_map = { "partition" : "slurm_partition", "mem_mb" : "mem_mb", 
                            "cpus" : "cpus_per_task", "time" : "runtime" };
    # Because I use slightly different resource names from what snakemake does for slurm

    rule_resources = {};
    for resource in keys:
        resource_value = getResource(config, rule_name, resource);
        if resource == "partition" and resource_value is None:
            continue;
        rule_resources[slurm_resource_map[resource]] = resource_value;

    # for key, value in rule_resources.items():
    #     meta_logger.info(f"Rule {rule_name} resource '{key}' set to {value}");

    return rule_resources

def getResource(config, rule_name, resource):
    # Get a specific resource value from the Snakemake config.yaml.
    rule_val = config.get("rule_resources", {}).get(rule_name, {}).get(resource)
    default_val = config.get("rule_resources", {}).get("default", {}).get(resource)

    if rule_val is not None:
        return rule_val
    elif default_val is not None:
        return default_val
    elif resource == "partition" and hasAutoPartitionSelection():
        return None
    else:
        meta_logger.error(f"Missing resource '{resource}' for rule '{rule_name}' and no default set.");
        raise ValueError();

#############################################################################    

def runCommand(cmd, log_stream, out_stream, rule, wc=""):

    if wc:
        wc = " ~ " + wc;
    # If a wild card is specified, add a hyphen so it is formatted
    # nicer in the print statements

    cmd_str = " ".join(cmd);
    # Join the command list into a string for printing

    printWrite(f"{fmtDT()} - RULE {rule}{wc} - INFO - Running command: {cmd_str}", log_stream);
    writeFlush("-" * 20 + " COMMAND LOG BEGIN " + "-" * 20 + "\n", log_stream);
    proc = subprocess.run(cmd, stdout=out_stream, stderr=log_stream);
    writeFlush("-" * 20 + "  COMMAND LOG END  " + "-" * 20 + "\n", log_stream);
    # Run the command and write the output to the log file

    rcode = proc.returncode;
    printWrite(f"{fmtDT()} - RULE {rule}{wc} - INFO - Command finished with return code: {rcode}", log_stream);
    # Get the return code

    if rcode != 0:
        printWrite(f"{fmtDT()} - RULE {rule}{wc} - ERROR - Command failed: {cmd_str}", log_stream);
        raise RuntimeError(f"Command failed with return code {rcode}: {cmd_str}");
    # If the command failed, raise an exception

#############################################################################
