import os
import sys

import lib.common as COMMON


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

run_phylofit = _as_bool(config.get("run_phylofit", config.get("run_neutral_model", True)), True)
run_phylop = _as_bool(config.get("run_phylop", True), True)
run_phastcons = _as_bool(config.get("run_phastcons", True), True)
build_cnees = _as_bool(config.get("build_cnees", config.get("make_cnees", True)), True)
cnee_output_format = str(config.get("cnee_output_format", "fasta")).strip().lower()
legacy_make_cnee_mafs = config.get("make_cnee_mafs", None)
legacy_cnee_extract_format = config.get("cnee_extract_format", config.get("cne_extract_format", None))
if cnee_output_format in {"fa", "fna"}:
    cnee_output_format = "fasta"
if cnee_output_format == "none" and legacy_make_cnee_mafs is not None and _as_bool(legacy_make_cnee_mafs, False):
    cnee_output_format = str(legacy_cnee_extract_format or "fasta").strip().lower()
    if cnee_output_format in {"fa", "fna"}:
        cnee_output_format = "fasta"
if cnee_output_format not in {"none", "fasta", "maf"}:
    raise ValueError(f"Invalid cnee_output_format '{cnee_output_format}'. Use 'none', 'fasta', or 'maf'.")
config["build_cnees"] = build_cnees
config["cnee_output_format"] = cnee_output_format

config["__master_workflow__"] = True
config_flag = config.get("display", False)
version_flag = config.get("version", False)
info_flag = config.get("info", False)
debug = config.get("debug", False)

MAIN, DRY_RUN, OUTPUT_DIR, LOG_DIR, TMPDIR, LOG_LEVEL, LOG_VERBOSITY = COMMON.pipelineSetup(
    config, sys.argv, version_flag, info_flag, config_flag, debug, workflow
)

NEUTRAL_MODEL_DIR = os.path.join(OUTPUT_DIR, "02-neutral-model")
PHYLOP_STAGE_DIR = os.path.join(OUTPUT_DIR, "03-phylop")

config["__pipeline_setup__"] = {
    "MAIN": MAIN,
    "DRY_RUN": DRY_RUN,
    "OUTPUT_DIR": OUTPUT_DIR,
    "LOG_DIR": LOG_DIR,
    "TMPDIR": TMPDIR,
    "LOG_LEVEL": LOG_LEVEL,
    "LOG_VERBOSITY": LOG_VERBOSITY,
}

if run_phylofit:
    include: "workflow/phylofit_models.smk"
if run_phylop:
    include: "workflow/phylop_regions.smk"
if run_phastcons:
    include: "workflow/phastcons_cnees.smk"

localrules: all

NEUTRAL_MODEL_TARGETS = []
PHYLOP_SITE_TARGETS = []
PHYLOP_REGION_TARGETS = []
ALL_BED_TARGETS = globals().get("ALL_BED_TARGETS", []) if run_phastcons else []
ALL_CNEES_TARGETS = globals().get("ALL_CNEES_TARGETS", []) if run_phastcons else []
ALL_CNEE_MAF_TARGETS = globals().get("ALL_CNEE_MAF_TARGETS", []) if run_phastcons else []

if run_phylofit:
    NEUTRAL_MODEL_TARGETS = expand(
        PHYLOFIT_ACTIVE_MODEL_PATH,
        zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES
    )

if run_phylop:
    PHYLOP_REGION_TARGETS = expand(
        os.path.join(PHYLOP_REGIONS_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.bed"),
        zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES
    )
    PHYLOP_SITE_TARGETS = expand(
        os.path.join(PHYLOP_STAGE_DIR, "summary", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.conserved-site-counts." + PHYLOP_ALPHA + ".tsv"),
        zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES
    ) + expand(
        os.path.join(PHYLOP_STAGE_DIR, "summary", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.accelerated-site-counts." + PHYLOP_ALPHA + ".tsv"),
        zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES
    )

ALL_TARGETS = (
    (NEUTRAL_MODEL_TARGETS if run_phylofit else [])
    + (PHYLOP_SITE_TARGETS if run_phylop else [])
    + (PHYLOP_REGION_TARGETS if run_phylop else [])
    + (ALL_BED_TARGETS if run_phastcons else [])
    + (ALL_CNEES_TARGETS if run_phastcons and build_cnees else [])
    + (ALL_CNEE_MAF_TARGETS if run_phastcons and build_cnees and cnee_output_format != "none" else [])
)

PIPELINE_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))
UTILS_DIR = os.path.join(PIPELINE_DIR, "utils")

SUMMARY_REPORT_PATH = os.path.join(OUTPUT_DIR, "summary_report.html")

rule all:
    input:
        ALL_TARGETS + [SUMMARY_REPORT_PATH]

rule summary_report:
    input:
        targets = ALL_TARGETS
        # Depending on the same conditional target lists rule all does is enough:
        # every intermediate file this report reads (filter_maf_by_gap manifests,
        # filter_4d_sites summaries, site_counts tsvs, etc.) is already a required
        # upstream dependency of these targets in the DAG, so it's guaranteed to
        # exist by the time this rule runs without needing to be listed separately.
    output:
        report = SUMMARY_REPORT_PATH
    params:
        script_path = os.path.join(UTILS_DIR, "summary_report.py"),
        manifest_path = os.path.join(LOG_DIR, "summary_report", "manifest.json"),
        manifest = {
            "output_dir": OUTPUT_DIR,
            "snakemake_command": config.get("__top_level_command__", ""),
            "main_inputs": {
                "output_dir": OUTPUT_DIR,
                "maf": config.get("maf"),
                "maf_ref_id": config.get("maf_ref_id"),
                "ref_fasta": config.get("ref_fasta"),
                "ref_gff": config.get("ref_gff"),
                "tree_file": config.get("tree_file"),
                "sample_file": config.get("sample_file"),
                "n_chromosome_groups": len(config["ref_chromosome_groups"]),
                "n_chromosomes": sum(len(v) for v in config["ref_chromosome_groups"].values()),
                "run_phylofit": run_phylofit,
                "run_phastcons": run_phastcons,
                "build_cnees": build_cnees,
                "cnee_output_format": cnee_output_format,
                "use_gc_corrected_models": USE_GC_CORRECTED_MODELS if run_phylofit else None,
                "rho_mode": RHO_MODE if run_phastcons else None,
            },
            "config_display": {
                k: str(v) for k, v in config.items()
                if not k.startswith("__") and k != "rule_resources"
            },
            "flags": {
                "run_phylofit": run_phylofit,
                "run_phastcons": run_phastcons,
                "run_phylop": run_phylop,
                "build_cnees": build_cnees,
                "cnee_output_format": cnee_output_format,
                "use_gc_corrected_models": USE_GC_CORRECTED_MODELS if run_phylofit else None,
                "rho_mode": RHO_MODE if run_phastcons else None,
            },
            "chromosome_groups": config["ref_chromosome_groups"],
            "maf_chr_prefix": config.get("maf_chr_prefix", ""),
            "paths": {
                "phylofit_dir": PHYLOFIT_DIR if run_phylofit else None,
                "neutral_summary_dir": NEUTRAL_SUMMARY_DIR if run_phylofit else None,
                "filter_threshold_4d": SEQ_THRESHOLD_4D if run_phylofit else None,
                "avg_gc_file": AVG_GC_FILE if run_phylofit else None,
                "maf_chunk_summary_dir": MAF_CHUNK_SUMMARY_DIR if run_phastcons else None,
                "maf_index_dir": MAF_INDEX_DIR if run_phastcons else None,
                "conserve_dir": CONSERVE_DIR if run_phastcons else None,
                "cnees_dir": CNEES_DIR if (run_phastcons and build_cnees) else None,
                "cnees_summary_dir": CNEES_SUMMARY_DIR if (run_phastcons and build_cnees) else None,
                "cnee_min_len_bp": CNEE_MIN_LEN_BP if (run_phastcons and build_cnees) else None,
                "cnee_density_bin_bp": CNEE_DENSITY_BIN_BP if (run_phastcons and build_cnees) else None,
            },
        }
    log:
        job_log = os.path.join(LOG_DIR, "summary_report", "run.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "summary_report", "run.txt")
    resources:
        **COMMON.getResources(config, "summary_report")
    run:
        import json
        import traceback

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(params.manifest_path), exist_ok=True)
                with open(params.manifest_path, "w") as mf:
                    json.dump(params.manifest, mf, indent=2, default=str)

                cmd = ["python", params.script_path, params.manifest_path, output.report]
                COMMON.runCommand(cmd, log_stream, log_stream, "summary_report")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise
