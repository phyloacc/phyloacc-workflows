#############################################################################
# Clean single-copy ortholog filtering for per-CNEE alignment FASTAs.
#
# `cnee_alignments_chr` extracts one aligned FASTA per CNEE (one row per species)
# with `mafutils fetch`. Because MAF blocks are contiguous in the REFERENCE only,
# a non-reference species' stitched row is not necessarily a contiguous single-
# copy locus. mafutils classifies every (element, species) as one of:
#     single | contiguous | split | multi_scaffold | multi_strand | no_bases
# and reports per-element how concentrated any splits are (whether one shared
# reference boundary explains them). This module turns those two mafutils TSVs
# (maf_fetch_summary.tsv, maf_fetch_loci.tsv) into a per-element decision:
#
#   - keep `single`/`contiguous` rows (a genuine single locus in that species);
#   - keep short `split` rows (<= split_max_gap_bp, a benign indel);
#   - keep ALL split rows when the element is a shared-boundary / reference-
#     deletion event (never drop these - the excluded sequence is non-reference
#     and unrecoverable, and the rows are still clean orthologs of the ref halves)
#     - the element is only FLAGGED;
#   - mask (-> N missing data) `multi_strand`/`multi_scaffold`/`no_bases` rows and
#     long lineage-specific `split` rows;
#   - the reference species is never masked and must be present;
#   - retain the element if >= min_species real (non-N) species remain, else the
#     caller moves it out of the PhyloAcc-input folder to a rejects tree.
#
# The logic here is pure (stdlib only) so it is unit-tested without Snakemake or
# mafutils; workflow/cnees.smk calls it after the fetch.
#############################################################################

from __future__ import annotations

import glob
import os
from collections import Counter
from typing import NamedTuple, Optional

#############################################################################
# mafutils TSV column names (matched by header, not position).

# maf_fetch_summary.tsv: one row per element. The element id is the `basename`
# column (mafutils fetch was run with `-b id`, so basename == CNEE id).
SUMMARY_ID_COL = "basename"
SUMMARY_N_SPLIT = "n.split"
SUMMARY_BOUNDARY_N_SPECIES = "split.max.boundary.n.species"
SUMMARY_BOUNDARY_GAP_SPREAD = "split.max.boundary.gap.spread"

# maf_fetch_loci.tsv: one row per (element, species, contributing block).
LOCI_ID_COL = "region_id"
LOCI_SPECIES_COL = "species"
LOCI_CLASS_COL = "class"
LOCI_MAXGAP_COL = "max_gap"
LOCI_SRC_SCAFFOLD_COL = "src_scaffold"
LOCI_SRC_SIZE_COL = "src_size"

# Class values.
CLASS_SINGLE = "single"
CLASS_CONTIGUOUS = "contiguous"
CLASS_SPLIT = "split"
CLASS_MULTI_SCAFFOLD = "multi_scaffold"
CLASS_MULTI_STRAND = "multi_strand"
CLASS_NO_BASES = "no_bases"

# Mask reasons (used for the per-chromosome funnel stats).
REASON_SPLIT_LONG = "split_long"
REASON_MULTI_SCAFFOLD = "multi_scaffold"
REASON_MULTI_STRAND = "multi_strand"
REASON_NO_BASES = "no_bases"
REASON_UNKNOWN = "unknown_class"

# Absolute floor on how many species must split at one boundary before it can be
# called a shared (ancestral / reference-deletion) event. 1-2 species agreeing is
# not evidence of a shared event. Documented constant, not a config knob.
SHARED_BOUNDARY_MIN_SPECIES_ABS = 3


class Thresholds(NamedTuple):
    split_max_gap_bp: int
    min_species: int
    shared_boundary_min_frac: float
    shared_boundary_max_gap_spread_bp: int


class SummaryRow(NamedTuple):
    region_id: str
    n_split: int
    boundary_n_species: Optional[int]
    boundary_gap_spread: Optional[int]


class LociInfo(NamedTuple):
    species: str
    cls: str
    max_gap: Optional[int]
    src_scaffold: Optional[str]
    src_size: Optional[int]


class ElementDecision(NamedTuple):
    region_id: str
    masked_species: frozenset          # species whose row -> N
    mask_reasons: Counter              # reason -> count (masked rows only)
    n_aligned: int                     # species with real bases before masking
    n_kept_real: int                   # real (non-N) species after masking
    retained: bool
    shared_boundary: bool
    ref_present: bool


#############################################################################
# Parsing helpers


def _to_int(x) -> Optional[int]:
    """mafutils writes '.', '' or 'None' where a number is undefined."""
    if x is None:
        return None
    s = str(x).strip()
    if s in ("", ".", "None", "NA", "nan"):
        return None
    try:
        return int(s)
    except ValueError:
        try:
            return int(float(s))
        except ValueError:
            return None


def _header_index(header_line: str, required, path: str):
    cols = header_line.rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(cols)}
    missing = [c for c in required if c not in idx]
    if missing:
        raise ValueError(
            f"{os.path.basename(path)} is missing expected column(s) {missing}; "
            f"got columns {cols}. This needs a mafutils build with the extended "
            f"summary/loci-table (scaffold headers, no_bases class, block_index/gap_to_next)."
        )
    return idx


def load_summary(path: str) -> dict:
    """maf_fetch_summary.tsv -> {region_id: SummaryRow}."""
    out = {}
    with open(path) as fh:
        header = fh.readline()
        if not header:
            return out
        idx = _header_index(
            header,
            (SUMMARY_ID_COL, SUMMARY_N_SPLIT, SUMMARY_BOUNDARY_N_SPECIES, SUMMARY_BOUNDARY_GAP_SPREAD),
            path,
        )
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rid = f[idx[SUMMARY_ID_COL]]
            out[rid] = SummaryRow(
                region_id=rid,
                n_split=_to_int(f[idx[SUMMARY_N_SPLIT]]) or 0,
                boundary_n_species=_to_int(f[idx[SUMMARY_BOUNDARY_N_SPECIES]]),
                boundary_gap_spread=_to_int(f[idx[SUMMARY_BOUNDARY_GAP_SPREAD]]),
            )
    return out


def load_loci(path: str) -> dict:
    """maf_fetch_loci.tsv -> {region_id: {species: LociInfo}}.

    class/max_gap are per (element, species) - mafutils denormalizes them onto every
    block row for that species, so the first row is authoritative (this is all the
    ortholog filter needs). src_scaffold/src_size are taken from the first row and are
    NOT safe for BED coordinate construction: the loci table is one row per contributing
    block and can carry a zero-size chunk on a different scaffold, so a real per-species
    BED must merge (min start, max end) over chunk_size>0 rows and forward-convert via
    src_size (see the mafutils README recipe), not read a single row here.
    """
    out: dict = {}
    with open(path) as fh:
        header = fh.readline()
        if not header:
            return out
        idx = _header_index(
            header,
            (LOCI_ID_COL, LOCI_SPECIES_COL, LOCI_CLASS_COL, LOCI_MAXGAP_COL),
            path,
        )
        has_src = LOCI_SRC_SCAFFOLD_COL in idx and LOCI_SRC_SIZE_COL in idx
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rid = f[idx[LOCI_ID_COL]]
            sp = f[idx[LOCI_SPECIES_COL]]
            per_sp = out.setdefault(rid, {})
            if sp in per_sp:
                continue  # class/max_gap are per-species; first row is authoritative
            per_sp[sp] = LociInfo(
                species=sp,
                cls=f[idx[LOCI_CLASS_COL]].strip(),
                max_gap=_to_int(f[idx[LOCI_MAXGAP_COL]]),
                src_scaffold=(f[idx[LOCI_SRC_SCAFFOLD_COL]] if has_src else None),
                src_size=(_to_int(f[idx[LOCI_SRC_SIZE_COL]]) if has_src else None),
            )
    return out


#############################################################################
# Policy


def is_shared_boundary(row: Optional[SummaryRow], thr: Thresholds) -> bool:
    """True when one reference boundary explains the element's splits (an ancestral
    indel / reference-specific deletion): enough species split, most of them at the
    same boundary, with tightly-agreeing gap sizes."""
    if row is None or row.n_split < SHARED_BOUNDARY_MIN_SPECIES_ABS:
        return False
    b_n, spread = row.boundary_n_species, row.boundary_gap_spread
    if b_n is None or spread is None:
        return False
    if b_n < SHARED_BOUNDARY_MIN_SPECIES_ABS:
        return False
    if b_n < thr.shared_boundary_min_frac * row.n_split:
        return False
    return spread <= thr.shared_boundary_max_gap_spread_bp


def classify_species(info: LociInfo, thr: Thresholds, shared_boundary: bool):
    """Return ('keep', None) or ('mask', reason) for one (element, species)."""
    cls = info.cls
    if cls in (CLASS_SINGLE, CLASS_CONTIGUOUS):
        return "keep", None
    if cls == CLASS_NO_BASES:
        return "mask", REASON_NO_BASES
    if cls == CLASS_MULTI_SCAFFOLD:
        return "mask", REASON_MULTI_SCAFFOLD
    if cls == CLASS_MULTI_STRAND:
        return "mask", REASON_MULTI_STRAND
    if cls == CLASS_SPLIT:
        if shared_boundary:
            return "keep", None
        if info.max_gap is not None and abs(info.max_gap) <= thr.split_max_gap_bp:
            return "keep", None
        return "mask", REASON_SPLIT_LONG
    # Unknown class -> mask (fail toward a clean set); surfaced in stats.
    return "mask", REASON_UNKNOWN


def decide_element(region_id: str, summary: dict, loci: dict,
                   ref_species: str, thr: Thresholds) -> ElementDecision:
    per_sp = loci.get(region_id, {})
    shared = is_shared_boundary(summary.get(region_id), thr)
    masked = set()
    reasons: Counter = Counter()
    for sp, info in per_sp.items():
        if sp == ref_species:
            continue  # reference is the anchor; never masked
        decision, reason = classify_species(info, thr, shared)
        if decision == "mask":
            masked.add(sp)
            reasons[reason] += 1
    n_aligned = len(per_sp)
    n_kept_real = n_aligned - len(masked)
    ref_present = ref_species in per_sp
    retained = ref_present and n_kept_real >= thr.min_species
    return ElementDecision(
        region_id=region_id,
        masked_species=frozenset(masked),
        mask_reasons=reasons,
        n_aligned=n_aligned,
        n_kept_real=n_kept_real,
        retained=retained,
        shared_boundary=shared,
        ref_present=ref_present,
    )


#############################################################################
# FASTA rewriting


def species_from_header(header: str) -> str:
    """Species token from a mafutils FASTA header, e.g.
    '>galGal6.chr1:97715-97807(+) id:...' or '>HLcolLiv2' -> the part before the
    first '.' or ':' of the first whitespace token."""
    name = header[1:].strip().split()[0] if header.startswith(">") else header.strip().split()[0]
    for i, ch in enumerate(name):
        if ch in ".:":
            return name[:i]
    return name


def filter_cnee_directory(fasta_paths, rejects_dir: str, summary: dict, loci: dict,
                          ref_species: str, thr: Thresholds):
    """Apply the ortholog filter to a directory of per-CNEE FASTAs.

    For each path: decide the element, mask retained files in place (bad rows -> N),
    move rejected files into rejects_dir. Returns (retained_basenames,
    rejected_basenames, stats) where stats has the `cnee_*` metric keys the report
    consumes. Does not write manifests - the caller owns that IO.
    """
    os.makedirs(rejects_dir, exist_ok=True)
    retained, rejected = [], []
    reason_totals: Counter = Counter()
    n_all_clean = n_masked_elems = n_shared = 0
    n_rej_minsp = n_rej_refabsent = 0
    for path in fasta_paths:
        base = os.path.basename(path)
        region_id = os.path.splitext(base)[0]
        d = decide_element(region_id, summary, loci, ref_species, thr)
        if d.retained:
            # stats below describe the retained clean set only: a rejected file is
            # moved wholesale, so its would-be-masked rows are never actually blanked.
            reason_totals.update(d.mask_reasons)
            if d.shared_boundary:
                n_shared += 1
            if d.masked_species:
                mask_fasta_records(path, path, d.masked_species)
                n_masked_elems += 1
            elif not d.shared_boundary:
                n_all_clean += 1
            retained.append(base)
        else:
            os.replace(path, os.path.join(rejects_dir, base))
            rejected.append(base)
            if not d.ref_present:
                n_rej_refabsent += 1
            else:
                n_rej_minsp += 1
    stats = {
        "cnee_fetched": len(fasta_paths),
        "cnee_retained": len(retained),
        "cnee_rejected": len(rejected),
        "cnee_elements_all_clean": n_all_clean,
        "cnee_elements_masked": n_masked_elems,
        "cnee_elements_shared_boundary": n_shared,
        "cnee_rows_masked_split_long": reason_totals.get(REASON_SPLIT_LONG, 0),
        "cnee_rows_masked_multi_scaffold": reason_totals.get(REASON_MULTI_SCAFFOLD, 0),
        "cnee_rows_masked_multi_strand": reason_totals.get(REASON_MULTI_STRAND, 0),
        "cnee_rows_masked_no_bases": reason_totals.get(REASON_NO_BASES, 0),
        "cnee_rows_masked_unknown": reason_totals.get(REASON_UNKNOWN, 0),
        "cnee_rejected_min_species": n_rej_minsp,
        "cnee_rejected_ref_absent": n_rej_refabsent,
    }
    return retained, rejected, stats


# The metric keys written to {chrom}.cnee-ortho-filter.tsv, in report order.
STATS_KEYS = (
    "cnee_fetched", "cnee_retained", "cnee_rejected",
    "cnee_elements_all_clean", "cnee_elements_masked", "cnee_elements_shared_boundary",
    "cnee_rows_masked_split_long", "cnee_rows_masked_multi_scaffold",
    "cnee_rows_masked_multi_strand", "cnee_rows_masked_no_bases", "cnee_rows_masked_unknown",
    "cnee_rejected_min_species", "cnee_rejected_ref_absent",
)


def stow_sidecars(outdir: str, keep: bool, dest_dir: Optional[str] = None,
                  prefix: str = "", patterns=("*.tsv", "*.log")):
    """Handle mafutils' fetch sidecars (maf_fetch_summary.tsv, maf_fetch_loci.tsv,
    maf_fetch.log) sitting beside the per-CNEE FASTAs, so the FASTA dir stays a clean
    `*.fa` + manifest.txt set for PhyloAcc.

    keep=True  -> MOVE them into dest_dir (default <outdir>/run-info), renamed to
                  <prefix><basename> so per-chromosome sidecars don't collide when several
                  chromosomes share one dest_dir (e.g. the source summary folder).
    keep=False -> delete them.
    Non-recursive: matches only the top level of outdir. Returns the destination basenames.
    """
    files = [f for pat in patterns for f in glob.glob(os.path.join(outdir, pat))
             if os.path.isfile(f)]
    handled = []
    if keep:
        if files:
            dest = dest_dir or os.path.join(outdir, "run-info")
            os.makedirs(dest, exist_ok=True)
            for fp in files:
                name = prefix + os.path.basename(fp)
                os.replace(fp, os.path.join(dest, name))
                handled.append(name)
    else:
        for fp in files:
            try:
                os.remove(fp)
            except OSError:
                pass
            handled.append(os.path.basename(fp))
    return handled


def mask_fasta_records(in_path: str, out_path: str, masked_species) -> int:
    """Rewrite a per-CNEE FASTA, replacing each masked species' sequence with
    'N' x (alignment length) so the row becomes missing data while the alignment
    stays rectangular. Returns the number of records masked. Writes atomically via
    a temp file even when out_path == in_path.
    """
    masked = set(masked_species)
    if not masked:
        if out_path != in_path:
            with open(in_path) as fi, open(out_path, "w") as fo:
                fo.write(fi.read())
        return 0

    # Parse into (header, seq) records (seq may span multiple lines).
    records = []
    with open(in_path) as fi:
        header = None
        seq_parts: list = []
        for line in fi:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line.rstrip("\n")
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_parts)))

    n_masked = 0
    tmp = out_path + ".tmp"
    with open(tmp, "w") as fo:
        for header, seq in records:
            if species_from_header(header) in masked:
                seq = "N" * len(seq)
                n_masked += 1
            fo.write(f"{header}\n{seq}\n")
    os.replace(tmp, out_path)
    return n_masked
