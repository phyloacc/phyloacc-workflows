#############################################################################
# Pure interval-list algorithms and small positional-data-format readers,
# shared across workflow/phastcons_cnees.smk, workflow/phylop_regions.smk, and
# utils/summary_report.py. Extracted from those rules' inline `run:` blocks so
# the logic is importable and unit-testable outside of a live Snakemake job
# (see tests/test_intervals.py) - a pure refactor, not a behavior change.
#############################################################################

import os

#############################################################################
# Generic interval-list algorithms

def merge_intervals(intervals, max_gap_bp=0):
    # intervals: iterable of (chrom, start, end) tuples.
    # Merges overlapping intervals, and near-adjacent ones within max_gap_bp.
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[1], x[2]))
    merged = [list(intervals[0])]
    for c, s, e in intervals[1:]:
        if s <= merged[-1][2] + max_gap_bp:
            if e > merged[-1][2]:
                merged[-1][2] = e
        else:
            merged.append([c, s, e])
    return [(c, s, e) for c, s, e in merged]


def drop_overlapping(a, b):
    # a, b: sorted, non-overlapping (chrom, start, end) lists (e.g. already run through
    # merge_intervals). Returns the subset of `a` with zero overlap in `b` - any interval
    # in `a` that overlaps an interval in `b` at all (fully containing it, splitting it,
    # or clipping either end) is dropped entirely, no fragments kept.
    a = sorted(a, key=lambda x: (x[1], x[2]))
    b = sorted(b, key=lambda x: (x[1], x[2]))
    out = []
    j = 0
    for chrom, s, e in a:
        # Advance j past any b-intervals that end at or before this a-interval starts -
        # they can't overlap this one or any later one (both lists are sorted by start).
        while j < len(b) and b[j][2] <= s:
            j += 1
        if j < len(b) and b[j][1] < e:
            continue
        out.append((chrom, s, e))
    return out


def complement_gaps(blocks, max_num_seqs_for_gap, min_gap_bp, min_keep_region_len):
    # blocks: iterable of (ref_start, ref_len, num_seqs) tuples (from a mafutils
    # .block.idx file - see read_chrom_length below for the same file format).
    # Merges consecutive low-coverage blocks (num_seqs <= max_num_seqs_for_gap) into
    # "gap runs", drops gap runs shorter than min_gap_bp, takes the complement of the
    # surviving gaps as candidate chunks, then drops chunks shorter than
    # min_keep_region_len. Returns a list of (start, end) chunk tuples.
    rows = sorted(blocks, key=lambda r: r[0])
    if not rows:
        return []

    chrom_start = rows[0][0]
    chrom_end = rows[-1][0] + rows[-1][1]

    gaps = []
    gap_start = None
    gap_end = None
    for ref_start, ref_len, num_seqs in rows:
        if num_seqs <= max_num_seqs_for_gap:
            if gap_start is None:
                gap_start = ref_start
            gap_end = ref_start + ref_len
        else:
            if gap_start is not None:
                gaps.append((gap_start, gap_end))
                gap_start = None
    if gap_start is not None:
        gaps.append((gap_start, gap_end))

    gaps = [(gs, ge) for gs, ge in gaps if ge - gs >= min_gap_bp]

    chunks = []
    cur = chrom_start
    for gs, ge in gaps:
        if gs > cur:
            chunks.append((cur, gs))
        cur = max(cur, ge)
    if cur < chrom_end:
        chunks.append((cur, chrom_end))

    return [(s, e) for s, e in chunks if e - s >= min_keep_region_len]


def cluster_sites(sites, max_gap_bp, min_count, min_len_bp):
    # sites: iterable of (chrom, start, end) tuples, assumed sorted/consecutive.
    # Greedily merges consecutive sites sharing a chromosome into regions when the gap
    # to the next site is <= max_gap_bp, tracking a count of merged sites per region;
    # drops regions with count < min_count or length < min_len_bp. Returns
    # (kept_regions, total_regions) - kept_regions is a list of (chrom, start, end,
    # count, idx) for surviving regions, where idx is the region's position in the
    # full (pre-filter) sequence, matching the original numbering so dropped regions
    # leave gaps in the id sequence rather than renumbering survivors; total_regions is
    # the pre-filter region count, for computing how many were dropped.
    regions = []
    current = None
    for chrom, start, end in sites:
        if current is None:
            current = {"chrom": chrom, "start": start, "end": end, "count": 1}
            continue
        if chrom == current["chrom"] and start - current["end"] <= max_gap_bp:
            current["end"] = max(current["end"], end)
            current["count"] += 1
        else:
            regions.append(current)
            current = {"chrom": chrom, "start": start, "end": end, "count": 1}
    if current is not None:
        regions.append(current)

    out = []
    for idx, region in enumerate(regions, start=1):
        region_len = region["end"] - region["start"]
        if region["count"] < min_count or region_len < min_len_bp:
            continue
        out.append((region["chrom"], region["start"], region["end"], region["count"], idx))
    return out, len(regions)

#############################################################################
# Small positional-data-format readers/writers

def parse_bed3(path, normalize_to=None):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            chrom = p[0]
            if normalize_to is not None:
                # Per-chromosome file - normalize aliases like NC_085107 -> NC_085107.1
                base_norm = normalize_to.rsplit(".", 1)[0]
                base_chrom = chrom.rsplit(".", 1)[0]
                if chrom == normalize_to or base_chrom == base_norm:
                    chrom = normalize_to
                else:
                    # Skip any unexpected chromosome labels in this per-chromosome file.
                    continue
            try:
                s = int(p[1])
                e = int(p[2])
            except ValueError:
                continue
            if e > s:
                rows.append((chrom, s, e))
    return rows


def filter_and_id_bed4(rows, min_len_bp, id_prefix):
    # rows: iterable of (chrom, start, end) tuples (e.g. from parse_bed3). Drops rows
    # with length <= min_len_bp, assigns sequential IDs "{id_prefix}.cnee{n:07d}" to
    # survivors. Returns a list of (chrom, start, end, id) 4-tuples.
    out = []
    n = 0
    for chrom, s, e in rows:
        if (e - s) <= min_len_bp:
            continue
        n += 1
        out.append((chrom, s, e, f"{id_prefix}.cnee{n:07d}"))
    return out


def gff_to_cds_bed(gff_lines, chrom):
    # gff_lines: iterable of raw GFF text lines. Filters to CDS features (case-
    # insensitive) on the given chromosome, converts 1-based inclusive GFF coordinates
    # to 0-based half-open BED, and drops malformed/non-positive-length rows. Returns a
    # list of (chrom, start, end) tuples.
    out = []
    for line in gff_lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        line_chrom, feature, start_s, end_s = parts[0], parts[2], parts[3], parts[4]
        if line_chrom != chrom or feature.upper() != "CDS":
            continue
        try:
            start = int(start_s) - 1
            end = int(end_s)
        except ValueError:
            continue
        if end > start >= 0:
            out.append((chrom, start, end))
    return out


def picard_interval_list_to_bed(lines, chrom):
    # lines: iterable of raw Picard interval-list text lines (may include "@" header
    # lines). Filters to the given chromosome, returns a list of "chrom:start-end"
    # strings (this pipeline's own on-disk convention for this intermediate format).
    out = []
    for line in lines:
        if line.startswith("@"):
            continue
        parts = line.strip().split("\t")
        line_chrom, start, end = parts[0], parts[1], parts[2]
        if line_chrom == chrom:
            out.append(f"{line_chrom}:{start}-{end}")
    return out


def tile_fixed_windows(chrom_len, window_size_bp, window_step_bp):
    # Tiles [0, chrom_len) into fixed-size windows of window_size_bp, stepping by
    # window_step_bp, with the last window clipped to chrom_len. Returns a list of
    # (start, end) tuples.
    windows = []
    start = 0
    while start < chrom_len:
        end = min(start + window_size_bp, chrom_len)
        windows.append((start, end))
        if end >= chrom_len:
            break
        start += window_step_bp
    return windows


def read_chrom_length(path):
    # Parses a mafutils-index .maf.block.idx file (columns: ref_scaff, ref_start,
    # ref_len, seq_len, line_len, num_seqs, byte_start, byte_end - see
    # complement_gaps above for the same file format). Coverage is gapless from
    # position 0 to the true chromosome end, so max(ref_start + ref_len) is the
    # chromosome length - no ref_fasta/.fai needed (not always available).
    if not path or not os.path.isfile(path):
        return None
    max_end = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                end = int(fields[1]) + int(fields[2])
            except ValueError:
                continue
            max_end = max(max_end, end)
    return max_end if max_end > 0 else None


def read_block_index_chroms(path):
    # Set of distinct reference chromosome names (column 1) in a mafutils .block.idx.
    chroms = set()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            name = line.split("\t", 1)[0].strip()
            if name:
                chroms.add(name)
    return chroms


def assert_bed_chroms_in_index(bed_path, block_index_path, rule_name):
    # Raise if any chromosome named in bed_path is absent from the MAF block index.
    # Catches a chromosome-name (prefix) mismatch loudly, before mafutils silently
    # extracts nothing for a name it can't find in the alignment.
    idx = read_block_index_chroms(block_index_path)
    bed = set()
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            bed.add(line.split("\t", 1)[0])
    missing = bed - idx
    if missing:
        raise ValueError(
            f"{rule_name}: chromosome(s) {sorted(missing)} in {bed_path} are not in the MAF "
            f"block index {block_index_path} (index has: {sorted(idx)[:8]}). This is a "
            f"chromosome-name mismatch between the bed and the MAF - check maf_prefix / gff_prefix."
        )


def maf_symlink_target_for_group(scaffold_index_path, group_bed_path):
    # Decide whether maf_split_chr_by_group can just SYMLINK the input MAF for this group
    # instead of re-running `mafutils fetch` (which, for a single-scaffold input, merely
    # copies the whole file). Returns the output basename to link to (the group bed's col4,
    # e.g. "chr1" -> link at "chr1.maf") IFF the input MAF is exactly one UNCOMPRESSED
    # scaffold that equals the group's single requested scaffold; otherwise None (must fetch).
    #
    # scaffold_index_path: the input MAF's .scaffold.idx (mafutils index; header line carries
    # "compression=none|gzip", data rows are "<scaffold>\t<offset>\t<size>").
    # group_bed_path: make_group_beds output (col1 = MAF scaffold name, col4 = output basename).
    compression = "none"
    scaffolds = set()
    with open(scaffold_index_path) as fh:
        for line in fh:
            if line.startswith("#"):
                if "compression=" in line:
                    compression = line.split("compression=", 1)[1].split()[0].strip()
                continue
            if not line.strip():
                continue
            scaffolds.add(line.split("\t", 1)[0].strip())
    # Only safe for a single UNCOMPRESSED scaffold: a symlink named "<x>.maf" must be a
    # plain MAF, and one input file can only map to one output file.
    if compression != "none" or len(scaffolds) != 1:
        return None
    bed_rows = []
    with open(group_bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                bed_rows.append((fields[0].strip(), fields[3].strip()))  # (scaffold, out_base)
    if len(bed_rows) != 1:
        return None
    scaffold, out_base = bed_rows[0]
    if scaffolds == {scaffold}:
        return out_base
    return None
