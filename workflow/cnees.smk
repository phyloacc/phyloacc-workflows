#############################################################################
# Source-agnostic CNEE stage.
#
# Turns a set of conserved elements (from ANY conservation source) into CNEEs:
#   merge within cnee_ces_merge_gap_bp -> drop CDS-overlapping elements ->
#   length-filter/assign IDs -> extract per-CNEE alignments (mafutils).
#
# The transformation is a pure function of (conserved-elements bed, CDS bed,
# MAF); the only per-source difference is which conserved-elements bed feeds it.
# So the rules carry a `{source}` wildcard and fan out over whichever sources
# are enabled:
#   - phastcons : phastCons conserved elements   (04-phastcons/regions)
#   - phylop    : phyloP-clustered conserved regions (03-phylop/regions)
# Outputs are namespaced by source: 05-cnees/{source}/... . The `phastcons`
# paths are identical to the pre-split layout, so downstream consumers (e.g. the
# summary report) are unaffected.
#
# This file is include-only (always run under the master Snakefile); it also
# hosts the shared per-chromosome MAF index (maf_index_chr), which the phastCons
# chunking rules depend on too - so it is included before phastcons_cnees.smk.
#############################################################################

import os
import glob
import traceback

from functools import partial

import lib.common as COMMON
import lib.intervals as INTERVALS
from lib.parsing import parse_newick_tip_names

_setup = config["__pipeline_setup__"]
OUTPUT_DIR = _setup["OUTPUT_DIR"]
LOG_DIR = _setup["LOG_DIR"]

getRuleResources = partial(COMMON.getResources, config)

#############################################################################
# Config / paths

# Config lists the CORE chromosome id; MAF/GFF names derive via maf_prefix/gff_prefix
# (maf_chr_prefix = legacy alias for maf_prefix). All bed contents / mafutils args use
# the MAF name; the GFF is read with the GFF name and its CDS relabeled to the MAF name.
MAF_PREFIX = str(config.get("maf_prefix", config.get("maf_chr_prefix", "")))
GFF_PREFIX = str(config.get("gff_prefix", ""))
MAF_CHR_PREFIX = MAF_PREFIX  # legacy alias
TREE_FILE = config.get("tree_file", "")
REF_GFF = config.get("ref_gff", "")

def maf_chrom(core):
    """Chromosome name as it appears in the MAF (for bed contents / mafutils args)."""
    return f"{MAF_PREFIX}{core}"

def gff_chrom(core):
    """Chromosome name as it appears in the GFF (for GFF reads only)."""
    return f"{GFF_PREFIX}{core}"

MAF_PREP_DIR = os.path.join(OUTPUT_DIR, "01-maf-prep")
MAF_INDEX_DIR = os.path.join(MAF_PREP_DIR, "maf-index")
MAF_SPLIT_BY_CHROM_DIR = COMMON.getOptionalConfigPath(
    config, "maf_split_chr_dir", os.path.join(MAF_PREP_DIR, "maf-by-chromosome"),
)

# Conserved-source region directories (mirror the producing modules' layouts).
PHYLOP_REGIONS_DIR = os.path.join(OUTPUT_DIR, "03-phylop", "regions")
PHASTCONS_REGIONS_DIR = os.path.join(OUTPUT_DIR, "04-phastcons", "regions")

CNEES_STAGE_DIR = os.path.join(OUTPUT_DIR, "05-cnees")
CDS_DIR = os.path.join(CNEES_STAGE_DIR, "cds")   # source-independent, built once

# Which conservation sources are active this run (independent stage toggles).
CNEE_SOURCES = []
if run_phastcons:
    CNEE_SOURCES.append("phastcons")
if run_phylop:
    CNEE_SOURCES.append("phylop")

#############################################################################
# CNEE parameters (build_cnees / cnee_output_format already normalized into
# config by the Snakefile).

MAKE_CNEES = _as_bool(config.get("build_cnees", True), True)
CNEE_OUTPUT_FORMAT = str(config.get("cnee_output_format", "fasta")).strip().lower()
MAKE_CNEE_MAFS = CNEE_OUTPUT_FORMAT != "none"

CNEE_CES_MERGE_GAP_BP = int(config.get("cnee_ces_merge_gap_bp", 5))
if CNEE_CES_MERGE_GAP_BP < 0:
    raise ValueError("cnee_ces_merge_gap_bp must be >= 0.")
CNEE_MIN_LEN_BP = int(config.get("cnee_min_len_bp", 50))
if CNEE_MIN_LEN_BP < 0:
    raise ValueError("cnee_min_len_bp must be >= 0.")
CNEE_FASTA_HEADER = str(config.get("cnee_fasta_header", config.get("cne_fasta_header", "species-coords-id"))).strip()

# How mafutils collapses a species that appears more than once in a block (paralogous /
# both-strand alignments): 'most-seq' keeps the copy with the most non-gap bases,
# 'none' emits every copy. Deep alignments (e.g. 241-way) routinely duplicate species,
# so 'none' leaves the downstream duplicate-species filter to discard nearly every CNEE;
# 'most-seq' is the right default for PhyloAcc-ready one-row-per-species alignments.
CNEE_FASTA_DEDUPE = str(config.get("cnee_fasta_dedupe", "most-seq")).strip().lower() or "most-seq"
if CNEE_FASTA_DEDUPE not in ("none", "most-seq"):
    raise ValueError("cnee_fasta_dedupe must be one of: none, most-seq.")

CNEE_EXPECTED_SPECIES = []
if MAKE_CNEES and CNEE_OUTPUT_FORMAT == "fasta":
    _explicit = [s.strip() for s in str(config.get("cnee_expected_species") or "").replace(",", " ").split() if s.strip()]
    _species_file = str(config.get("cnee_expected_species_file") or "").strip()
    if _species_file:
        with open(_species_file, "r", encoding="utf-8") as fp:
            for line in fp:
                sp = line.strip()
                if sp and not sp.startswith("#") and sp not in _explicit:
                    _explicit.append(sp)
    if _explicit:
        CNEE_EXPECTED_SPECIES = _explicit
    elif TREE_FILE:
        CNEE_EXPECTED_SPECIES = parse_newick_tip_names(TREE_FILE)
        if not CNEE_EXPECTED_SPECIES:
            raise ValueError(f"No tip labels parsed from tree_file '{TREE_FILE}' for CNEE FASTA extraction.")
    else:
        raise ValueError(
            "cnee_output_format=fasta requires species information: set "
            "cnee_expected_species, cnee_expected_species_file, or tree_file."
        )

if MAKE_CNEES and not REF_GFF:
    raise ValueError("build_cnees=true requires ref_gff to be set in config (CDS coordinates for CNEE calling).")

KEEP_CNEE_SIDECARS = _as_bool(config.get("keep_cnee_sidecars", False), False)

#############################################################################
# Chromosome list (recomputed here so this file is order-independent).

_target_chrom_cfg = config.get("target_ref_chromosomes", [])
if isinstance(_target_chrom_cfg, str):
    _target_chrom_set = set(c.strip() for c in _target_chrom_cfg.split(",") if c.strip())
elif isinstance(_target_chrom_cfg, (list, tuple, set)):
    _target_chrom_set = set(_target_chrom_cfg)
else:
    _target_chrom_set = set()

_flattened = [
    (group, chromosome)
    for group, chromosome_list in config["ref_chromosome_groups"].items()
    for chromosome in chromosome_list
    if not _target_chrom_set or chromosome in _target_chrom_set
]
if not _flattened:
    raise ValueError("No chromosomes selected. Check ref_chromosome_groups/target_ref_chromosomes in config.")
_CNEE_CHR_GROUPS, _CNEE_CHROMS = zip(*_flattened)

#############################################################################
# Source -> conserved-elements bed resolver.

def _conserved_bed_for_source(wildcards):
    if wildcards.source == "phastcons":
        return os.path.join(PHASTCONS_REGIONS_DIR, wildcards.chromosome_group, f"{wildcards.ref_chromosome}.bed")
    if wildcards.source == "phylop":
        return os.path.join(PHYLOP_REGIONS_DIR, wildcards.chromosome_group, f"{MAF_CHR_PREFIX}{wildcards.ref_chromosome}.bed")
    raise ValueError(f"Unknown CNEE source '{wildcards.source}'")

# Per-chromosome MAF resolver for the shared index. Copied from the phastCons
# stage verbatim to preserve maf_index_chr behavior exactly (note: this matches
# the per-chromosome MAF files carry maf_chr_prefix ("chr1.maf"), matching how
# maf_split_chr_by_group names them and how the phyloP stage looks them up.
def _cnee_chr_maf(wildcards):
    fname_target = f"{MAF_CHR_PREFIX}{wildcards.ref_chromosome}.maf"
    if "checkpoints" in globals() and hasattr(checkpoints, "maf_split_chr_by_group"):
        cp_output = checkpoints.maf_split_chr_by_group.get(chromosome_group=wildcards.chromosome_group).output
        group_dir = os.path.join(MAF_SPLIT_BY_CHROM_DIR, wildcards.chromosome_group)
        with open(cp_output.maf_manifest, "r") as mf:
            expected = [line.strip() for line in mf if line.strip()]
        for fname in expected:
            if fname == fname_target:
                return os.path.join(group_dir, fname)
        raise ValueError(f"No chromosome MAF for {wildcards.ref_chromosome} listed in {cp_output.maf_manifest}")
    return os.path.join(MAF_SPLIT_BY_CHROM_DIR, wildcards.chromosome_group, fname_target)

#############################################################################
# Target lists (fan out over active sources x chromosomes).

def _expand_over_sources(path_template):
    out = []
    for src in CNEE_SOURCES:
        out += expand(path_template, zip,
                      source=[src] * len(_CNEE_CHROMS),
                      chromosome_group=list(_CNEE_CHR_GROUPS),
                      ref_chromosome=list(_CNEE_CHROMS))
    return out

ALL_CNEES_TARGETS = _expand_over_sources(
    os.path.join(CNEES_STAGE_DIR, "{source}", "bed", "{chromosome_group}", "{ref_chromosome}.cnees.bed")
)
ALL_CNEE_MAF_TARGETS = _expand_over_sources(
    os.path.join(CNEES_STAGE_DIR, "{source}", CNEE_OUTPUT_FORMAT if CNEE_OUTPUT_FORMAT == "fasta" else "maf",
                 "{chromosome_group}", "{ref_chromosome}", "manifest.txt")
)

wildcard_constraints:
    source = r"phastcons|phylop"

#############################################################################

rule maf_index_chr:
    input:
        maf = _cnee_chr_maf
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

rule extract_cds_bed_chr:
    input:
        ref_gff = REF_GFF
    output:
        cds_bed = os.path.join(CDS_DIR, "{chromosome_group}", "{ref_chromosome}.cds.bed")
    log:
        job_log = os.path.join(LOG_DIR, "extract_cds_bed_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "extract_cds_bed_chr", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("extract_cds_bed_chr")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.cds_bed), exist_ok=True)
                # Match CDS rows by the GFF's chromosome name, but write the CDS bed with the
                # MAF name so it lines up with the (MAF-named) conserved elements at drop time.
                mchrom = maf_chrom(wildcards.ref_chromosome)
                with open(input.ref_gff) as gf:
                    cds_rows = INTERVALS.gff_to_cds_bed(gf, gff_chrom(wildcards.ref_chromosome))
                with open(output.cds_bed, "w") as out:
                    for _chrom, start, end in cds_rows:
                        out.write(f"{mchrom}\t{start}\t{end}\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule cnees_from_conserved_chr:
    input:
        conserved_bed = _conserved_bed_for_source,
        cds_bed = rules.extract_cds_bed_chr.output.cds_bed
    output:
        cnees_bed = os.path.join(CNEES_STAGE_DIR, "{source}", "bed", "{chromosome_group}", "{ref_chromosome}.cnees.bed"),
        filter_summary = os.path.join(CNEES_STAGE_DIR, "{source}", "summary", "{chromosome_group}", "{ref_chromosome}.cnees-filter-summary.tsv")
    log:
        job_log = os.path.join(LOG_DIR, "cnees_from_conserved_chr", "{source}", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "cnees_from_conserved_chr", "{source}", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("cnees_from_conserved_chr")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.cnees_bed), exist_ok=True)
                mchrom = maf_chrom(wildcards.ref_chromosome)  # CNEEs are MAF-coordinate; use the MAF name
                conserved_raw = INTERVALS.parse_bed3(input.conserved_bed, normalize_to=mchrom)

                # Graceful skip-and-succeed: an empty conserved-elements bed (e.g. phyloP
                # found no FDR-significant sites on a shallow tree) yields an empty, clearly
                # labeled CNEE set rather than an error or downstream churn.
                if not conserved_raw:
                    log_stream.write(
                        f"No conserved elements for {wildcards.source}/{wildcards.ref_chromosome}; "
                        "writing empty CNEE set. For phyloP this usually means no FDR-significant "
                        "conserved sites - see analyses/phylop-tree-length-power.\n"
                    )

                conserved = INTERVALS.merge_intervals(conserved_raw, CNEE_CES_MERGE_GAP_BP)
                cds = INTERVALS.merge_intervals(
                    INTERVALS.parse_bed3(input.cds_bed, normalize_to=mchrom), 0
                )
                log_stream.write(
                    f"Conserved raw intervals: {len(conserved_raw)}; "
                    f"merged (gap<={CNEE_CES_MERGE_GAP_BP}bp): {len(conserved)}\n"
                )

                out_rows = INTERVALS.drop_overlapping(conserved, cds)
                ces_dropped = len(conserved) - len(out_rows)
                log_stream.write(
                    f"CEs dropped for CDS overlap: {ces_dropped}; CNEEs remaining: {len(out_rows)}\n"
                )

                with open(output.cnees_bed, "w") as out:
                    for chrom, s, e in out_rows:
                        if e > s:
                            out.write(f"{chrom}\t{s}\t{e}\n")

                os.makedirs(os.path.dirname(output.filter_summary), exist_ok=True)
                with open(output.filter_summary, "w") as sf:
                    sf.write("metric\tvalue\n")
                    sf.write(f"ces_raw\t{len(conserved_raw)}\n")
                    sf.write(f"ces_merged\t{len(conserved)}\n")
                    sf.write(f"ces_dropped_cds_overlap\t{ces_dropped}\n")
                    sf.write(f"cnees_after_cds_drop\t{len(out_rows)}\n")
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule cnees_to_bed4_chr:
    input:
        cnees_bed = rules.cnees_from_conserved_chr.output.cnees_bed
    output:
        cnees_bed4 = os.path.join(CNEES_STAGE_DIR, "{source}", "bed", "{chromosome_group}", "{ref_chromosome}.cnees.bed4")
    log:
        job_log = os.path.join(LOG_DIR, "cnees_to_bed4_chr", "{source}", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "cnees_to_bed4_chr", "{source}", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("cnees_to_bed4_chr")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.cnees_bed4), exist_ok=True)
                rows = INTERVALS.parse_bed3(input.cnees_bed)
                bed4_rows = INTERVALS.filter_and_id_bed4(rows, CNEE_MIN_LEN_BP, maf_chrom(wildcards.ref_chromosome))
                n = len(bed4_rows)
                dropped = len(rows) - n
                with open(output.cnees_bed4, "w") as outf:
                    for chrom, s, e, cid in bed4_rows:
                        outf.write(f"{chrom}\t{s}\t{e}\t{cid}\n")
                log_stream.write(
                    f"Filtered CNEEs by length > {CNEE_MIN_LEN_BP} bp: kept={n}, dropped={dropped}\n"
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
        manifest = os.path.join(CNEES_STAGE_DIR, "{source}",
                                CNEE_OUTPUT_FORMAT if CNEE_OUTPUT_FORMAT == "fasta" else "maf",
                                "{chromosome_group}", "{ref_chromosome}", "manifest.txt")
    params:
        outdir = os.path.join(CNEES_STAGE_DIR, "{source}",
                              CNEE_OUTPUT_FORMAT if CNEE_OUTPUT_FORMAT == "fasta" else "maf",
                              "{chromosome_group}", "{ref_chromosome}"),
        rule_name = "cnee_alignments_chr"
    log:
        job_log = os.path.join(LOG_DIR, "cnee_alignments_chr", "{source}", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "cnee_alignments_chr", "{source}", "{chromosome_group}", "{ref_chromosome}.txt")
    resources:
        **getRuleResources("cnee_alignments_chr")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(params.outdir, exist_ok=True)
                mchrom = maf_chrom(wildcards.ref_chromosome)  # CNEE ids/files are MAF-named
                # Remove stale per-CNEE outputs so manifest reflects current filtering.
                stale = glob.glob(os.path.join(params.outdir, f"{mchrom}.cnee*.maf"))
                stale += glob.glob(os.path.join(params.outdir, f"{mchrom}.cnee*.fa"))
                for fp in stale:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass

                # Graceful short-circuit: no CNEEs to extract (e.g. an empty phyloP
                # source) -> write an empty manifest, don't invoke mafutils on empty input.
                n_cnees = sum(1 for line in open(input.cnees_bed4) if line.strip())
                if n_cnees == 0:
                    with open(output.manifest, "w") as out:
                        pass
                    log_stream.write("No CNEEs to extract; wrote empty manifest.\n")
                else:
                    # Fail loud on a chromosome-name mismatch before mafutils silently drops it.
                    INTERVALS.assert_bed_chroms_in_index(input.cnees_bed4, input.maf_index_block, "cnee_alignments_chr")
                    p = min(int(resources.cpus_per_task), 4)
                    cmd = [
                        "mafutils", "fetch",
                        input.maf,
                        input.cnees_bed4,
                        "-i", input.maf_index_block,
                        "-o", params.outdir,
                        "-p", str(p),
                        "-m", "block",
                        "-b", "id",
                    ]
                    if CNEE_OUTPUT_FORMAT == "fasta":
                        cmd += ["-f", "-fh", CNEE_FASTA_HEADER]
                        cmd += ["--fasta-dedupe", CNEE_FASTA_DEDUPE]
                        if CNEE_EXPECTED_SPECIES:
                            cmd += ["--expected-species", ",".join(CNEE_EXPECTED_SPECIES)]
                    COMMON.runCommand(
                        cmd, log_stream, log_stream, params.rule_name,
                        wc=f"{wildcards.source}.{wildcards.chromosome_group}.{wildcards.ref_chromosome}"
                    )

                    ext = "fa" if CNEE_OUTPUT_FORMAT == "fasta" else "maf"
                    outs = sorted(glob.glob(os.path.join(params.outdir, f"{mchrom}.cnee*.{ext}")))

                    # Duplicate species (paralogous / both-strand alignments) are collapsed
                    # upstream by `mafutils fetch --fasta-dedupe` (cnee_fasta_dedupe), so every
                    # fetched FASTA already carries at most one record per species - manifest
                    # them all as-is.
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
