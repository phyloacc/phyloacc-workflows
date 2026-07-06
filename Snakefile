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

rule all:
    input:
        (NEUTRAL_MODEL_TARGETS if run_phylofit else [])
        + (PHYLOP_SITE_TARGETS if run_phylop else [])
        + (PHYLOP_REGION_TARGETS if run_phylop else [])
        + (ALL_BED_TARGETS if run_phastcons else [])
        + (ALL_CNEES_TARGETS if run_phastcons and build_cnees else [])
        + (ALL_CNEE_MAF_TARGETS if run_phastcons and build_cnees and cnee_output_format != "none" else [])
