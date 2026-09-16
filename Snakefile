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

# --- align-only mode: inferred from presence of align_input_dir (no explicit switch) ---
# If align_input_dir is set, bypass the entire MAF-based prediction pipeline and only align
# a directory of user-supplied unaligned per-element FASTAs (workflow/align_elements.smk).
align_input_dir = str(config.get("align_input_dir", "")).strip()
align_only = bool(align_input_dir)
if align_only:
    if str(config.get("maf", "")).strip():
        raise ValueError(
            "Ambiguous mode: set either 'maf' (predict pipeline) OR 'align_input_dir' "
            "(align-only), not both."
        )
    run_phylofit = run_phylop = run_phastcons = build_cnees = False
    config["build_cnees"] = False

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
# Shared CNEE stage. Included before phastcons_cnees.smk because it hosts the
# per-chromosome MAF index the phastCons chunking rules depend on, and it builds
# CNEEs source-namespaced (phastcons + phylop). Needed whenever phastCons runs
# (for the MAF index) or CNEEs are built from any enabled source.
include_cnees = run_phastcons or (build_cnees and run_phylop)
if include_cnees:
    include: "workflow/cnees.smk"
if run_phastcons:
    include: "workflow/phastcons_cnees.smk"
if align_only:
    include: "workflow/align_elements.smk"

localrules: all

NEUTRAL_MODEL_TARGETS = []
PHYLOP_SITE_TARGETS = []
PHYLOP_REGION_TARGETS = []
ALL_BED_TARGETS = globals().get("ALL_BED_TARGETS", []) if run_phastcons else []
# CNEE target lists come from workflow/cnees.smk (source-namespaced, fanned out
# over all active conservation sources).
ALL_CNEES_TARGETS = globals().get("ALL_CNEES_TARGETS", []) if include_cnees else []
ALL_CNEE_MAF_TARGETS = globals().get("ALL_CNEE_MAF_TARGETS", []) if include_cnees else []

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
    + (ALL_CNEES_TARGETS if build_cnees else [])
    + (ALL_CNEE_MAF_TARGETS if build_cnees and cnee_output_format != "none" else [])
)

# In align-only mode the sole target is the aligned-elements manifest from
# workflow/align_elements.smk; none of the MAF-based sublists above apply.
ALL_ALIGN_TARGETS = globals().get("ALL_ALIGN_TARGETS", []) if align_only else []
if align_only:
    ALL_TARGETS = list(ALL_ALIGN_TARGETS)

PIPELINE_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))
UTILS_DIR = os.path.join(PIPELINE_DIR, "utils")

# Prefix the report with the config-file stem (e.g. mammals-2mb-chunk.yaml ->
# mammals-2mb-chunk.summary_report.html) so reports from different runs are
# self-identifying; fall back to a bare name when run via --config (no file).
_report_configfiles = list(getattr(workflow, "configfiles", []) or [])
_report_stem = os.path.splitext(os.path.basename(str(_report_configfiles[0])))[0] if _report_configfiles else ""
SUMMARY_REPORT_PATH = os.path.join(
    OUTPUT_DIR, f"{_report_stem}.summary_report.html" if _report_stem else "summary_report.html"
)

rule all:
    input:
        ALL_TARGETS + ([] if align_only else [SUMMARY_REPORT_PATH])

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
                "n_chromosome_groups": len(config.get("ref_chromosome_groups", {})),
                "n_chromosomes": sum(len(v) for v in config.get("ref_chromosome_groups", {}).values()),
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
            "chromosome_groups": config.get("ref_chromosome_groups", {}),
            # Resolve the canonical MAF chromosome prefix from the new maf_prefix key,
            # falling back to the legacy maf_chr_prefix alias - the summary collectors
            # build MAF-named filenames ("chr1.bed") from this, so an unset value on a
            # prefixed config would silently miss every per-chromosome file.
            "maf_chr_prefix": config.get("maf_prefix", config.get("maf_chr_prefix", "")),
            "paths": {
                "phylofit_dir": PHYLOFIT_DIR if run_phylofit else None,
                "neutral_summary_dir": NEUTRAL_SUMMARY_DIR if run_phylofit else None,
                "filter_threshold_4d": SEQ_THRESHOLD_4D if run_phylofit else None,
                "avg_gc_file": AVG_GC_FILE if run_phylofit else None,
                "maf_chunk_summary_dir": MAF_CHUNK_SUMMARY_DIR if run_phastcons else None,
                # Per-chromosome block index is now co-located with the chromosome MAF
                # (maf_index_chr), named <prefix><chrom>.maf.block.idx - not the old maf-index/ dir.
                "maf_index_dir": MAF_SPLIT_BY_CHROM_DIR if (run_phastcons or (run_phylop and build_cnees)) else None,
                "conserve_dir": CONSERVE_DIR if run_phastcons else None,
                "cnees_dir": CNEES_DIR if (run_phastcons and build_cnees) else None,
                "cnees_summary_dir": CNEES_SUMMARY_DIR if (run_phastcons and build_cnees) else None,
                "cnee_min_len_bp": CNEE_MIN_LEN_BP if (build_cnees and (run_phastcons or run_phylop)) else None,
                "cnee_density_bin_bp": CNEE_DENSITY_BIN_BP if (run_phastcons and build_cnees) else None,
                # phyloP branch (per-site LRT -> clustered regions -> phyloP-source CNEEs)
                "phylop_power_dir": os.path.join(PHYLOP_STAGE_DIR, "power-check") if run_phylop else None,
                "phylop_summary_dir": PHYLOP_SUMMARY_DIR if run_phylop else None,
                "phylop_regions_dir": PHYLOP_REGIONS_DIR if run_phylop else None,
                "phylop_alpha": PHYLOP_ALPHA if run_phylop else None,
                "phylop_cluster_method": PHYLOP_CLUSTER_METHOD if run_phylop else None,
                "phylop_power_gate_enabled": PHYLOP_POWER_GATE if run_phylop else None,
                "phylop_cnees_dir": os.path.join(OUTPUT_DIR, "05-cnees", "phylop", "bed") if (run_phylop and build_cnees) else None,
                "phylop_cnees_summary_dir": os.path.join(OUTPUT_DIR, "05-cnees", "phylop", "summary") if (run_phylop and build_cnees) else None,
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
