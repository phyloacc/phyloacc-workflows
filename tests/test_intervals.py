#############################################################################
# Unit tests for lib/intervals.py - pure interval-list algorithms and small
# positional-data-format readers extracted from workflow/*.smk rules'
# inline `run:` blocks. Run with `pytest tests/` from the repo root.
#############################################################################

import pytest

import lib.intervals as intervals

#############################################################################
# merge_intervals

def test_merge_intervals_overlapping():
    assert intervals.merge_intervals([("c", 100, 200), ("c", 150, 250)]) == [("c", 100, 250)]


def test_merge_intervals_adjacent_within_gap():
    assert intervals.merge_intervals([("c", 100, 200), ("c", 203, 300)], max_gap_bp=5) == [("c", 100, 300)]


def test_merge_intervals_beyond_gap_stays_separate():
    result = intervals.merge_intervals([("c", 100, 200), ("c", 210, 300)], max_gap_bp=5)
    assert result == [("c", 100, 200), ("c", 210, 300)]


def test_merge_intervals_empty():
    assert intervals.merge_intervals([]) == []

#############################################################################
# drop_overlapping (CDS-overlap-drop sweep)

def test_drop_overlapping_cds_fully_inside_ce():
    # A CDS fully inside a CE - the whole CE is dropped, no flanking fragments kept.
    result = intervals.drop_overlapping([("c", 100, 200)], [("c", 140, 160)])
    assert result == []


def test_drop_overlapping_cds_clips_one_end():
    result = intervals.drop_overlapping([("c", 100, 200)], [("c", 180, 220)])
    assert result == []


def test_drop_overlapping_no_overlap():
    result = intervals.drop_overlapping([("c", 100, 200)], [("c", 300, 400)])
    assert result == [("c", 100, 200)]


def test_drop_overlapping_touching_boundary_is_not_overlap():
    # CDS starts exactly where the CE ends - touching, not overlapping.
    result = intervals.drop_overlapping([("c", 100, 200)], [("c", 200, 300)])
    assert result == [("c", 100, 200)]


def test_drop_overlapping_mixed_multiple_ces():
    conserved = [("c", 100, 200), ("c", 300, 400), ("c", 500, 600)]
    cds = [("c", 150, 160), ("c", 550, 560)]
    assert intervals.drop_overlapping(conserved, cds) == [("c", 300, 400)]

#############################################################################
# complement_gaps

def test_complement_gaps_basic_split():
    # blocks: (ref_start, ref_len, num_seqs). A low-coverage run in the middle splits
    # the chromosome into two chunks.
    blocks = [(0, 100, 10), (100, 100, 1), (200, 100, 10)]
    chunks = intervals.complement_gaps(blocks, max_num_seqs_for_gap=3, min_gap_bp=50, min_keep_region_len=1)
    assert chunks == [(0, 100), (200, 300)]


def test_complement_gaps_short_gap_below_min_gap_bp_does_not_split():
    blocks = [(0, 100, 10), (100, 10, 1), (110, 100, 10)]
    chunks = intervals.complement_gaps(blocks, max_num_seqs_for_gap=3, min_gap_bp=50, min_keep_region_len=1)
    assert chunks == [(0, 210)]


def test_complement_gaps_drops_short_resulting_chunk():
    blocks = [(0, 5, 10), (5, 100, 1), (105, 100, 10)]
    chunks = intervals.complement_gaps(blocks, max_num_seqs_for_gap=3, min_gap_bp=50, min_keep_region_len=10)
    assert chunks == [(105, 205)]


def test_complement_gaps_no_blocks():
    assert intervals.complement_gaps([], max_num_seqs_for_gap=3, min_gap_bp=50, min_keep_region_len=1) == []

#############################################################################
# cluster_sites

def test_cluster_sites_merges_within_gap():
    sites = [("c", 100, 110), ("c", 115, 120), ("c", 200, 210)]
    kept, total = intervals.cluster_sites(sites, max_gap_bp=10, min_count=1, min_len_bp=1)
    assert total == 2
    assert kept == [("c", 100, 120, 2, 1), ("c", 200, 210, 1, 2)]


def test_cluster_sites_chromosome_boundary_forces_split():
    sites = [("c1", 100, 110), ("c2", 111, 120)]
    kept, total = intervals.cluster_sites(sites, max_gap_bp=100, min_count=1, min_len_bp=1)
    assert total == 2
    assert [k[0] for k in kept] == ["c1", "c2"]


def test_cluster_sites_drops_below_min_count_and_min_len():
    sites = [("c", 100, 105), ("c", 200, 201)]
    kept, total = intervals.cluster_sites(sites, max_gap_bp=1, min_count=2, min_len_bp=1)
    assert total == 2
    assert kept == []  # neither region has count >= 2


def test_cluster_sites_dropped_region_leaves_gap_in_idx():
    # Middle region (idx=2, length 2) is below min_len_bp=5 and dropped; surviving
    # idx values (1, 3) should not be renumbered.
    sites = [("c", 0, 10), ("c", 100, 102), ("c", 200, 210)]
    kept, total = intervals.cluster_sites(sites, max_gap_bp=1, min_count=1, min_len_bp=5)
    assert total == 3
    idxs = [k[4] for k in kept]
    assert idxs == [1, 3]

#############################################################################
# parse_bed3

def test_parse_bed3_exact_and_version_suffix_match(tmp_path):
    bed = tmp_path / "in.bed"
    bed.write_text("CM000994.3\t10\t20\nCM000994\t30\t40\nchrX\t50\t60\n")
    rows = intervals.parse_bed3(str(bed), normalize_to="CM000994.3")
    assert rows == [("CM000994.3", 10, 20), ("CM000994.3", 30, 40)]


def test_parse_bed3_no_normalize_keeps_all_chroms(tmp_path):
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t10\t20\nchr2\t30\t40\n")
    rows = intervals.parse_bed3(str(bed))
    assert rows == [("chr1", 10, 20), ("chr2", 30, 40)]


def test_parse_bed3_drops_malformed_and_zero_length(tmp_path):
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t10\t20\nchr1\tnotanumber\t20\nchr1\t50\t50\n")
    rows = intervals.parse_bed3(str(bed))
    assert rows == [("chr1", 10, 20)]

#############################################################################
# filter_and_id_bed4

def test_filter_and_id_bed4_length_threshold_and_ids():
    rows = [("c", 0, 100), ("c", 200, 210), ("c", 300, 500)]
    result = intervals.filter_and_id_bed4(rows, min_len_bp=50, id_prefix="chr1")
    assert result == [("c", 0, 100, "chr1.cnee0000001"), ("c", 300, 500, "chr1.cnee0000002")]

#############################################################################
# gff_to_cds_bed

def test_gff_to_cds_bed_filters_feature_and_chrom():
    lines = [
        "chr1\tsrc\tCDS\t11\t20\t.\t+\t.\tid=1\n",
        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tid=2\n",
        "chr2\tsrc\tCDS\t11\t20\t.\t+\t.\tid=3\n",
        "chr1\tsrc\tcds\t31\t40\t.\t+\t.\tid=4\n",  # lowercase feature
    ]
    rows = intervals.gff_to_cds_bed(lines, "chr1")
    assert rows == [("chr1", 10, 20), ("chr1", 30, 40)]


def test_gff_to_cds_bed_drops_malformed_coordinates():
    lines = ["chr1\tsrc\tCDS\tnotanumber\t20\t.\t+\t.\tid=1\n"]
    assert intervals.gff_to_cds_bed(lines, "chr1") == []

#############################################################################
# picard_interval_list_to_bed

def test_picard_interval_list_to_bed_skips_header_and_other_chroms():
    lines = [
        "@HD\tVN:1.6\n",
        "chr1\t10\t20\t+\tinterval1\n",
        "chr2\t30\t40\t+\tinterval2\n",
    ]
    assert intervals.picard_interval_list_to_bed(lines, "chr1") == ["chr1:10-20"]

#############################################################################
# tile_fixed_windows

def test_tile_fixed_windows_clips_last_window():
    windows = intervals.tile_fixed_windows(250, window_size_bp=100, window_step_bp=100)
    assert windows == [(0, 100), (100, 200), (200, 250)]


def test_tile_fixed_windows_overlap_step():
    windows = intervals.tile_fixed_windows(150, window_size_bp=100, window_step_bp=50)
    assert windows == [(0, 100), (50, 150)]

#############################################################################
# read_chrom_length

def test_read_chrom_length_gapless_block_index(tmp_path):
    idx = tmp_path / "chr1.maf.block.idx"
    idx.write_text(
        "# mafutils-index format=2 maf=chr1.maf compression=none size=1 mtime=1.0 hash=md5:x\n"
        "\n"
        "chr1\t0\t100\t100\t101\t5\t0\t100\n"
        "chr1\t100\t50\t50\t51\t5\t100\t150\n"
    )
    assert intervals.read_chrom_length(str(idx)) == 150


def test_read_chrom_length_missing_file_returns_none(tmp_path):
    assert intervals.read_chrom_length(str(tmp_path / "does-not-exist.idx")) is None


#############################################################################
# block-index chromosome set + fail-loud name check

def _write_idx(path):
    path.write_text(
        "# mafutils-index format=2 maf=x.maf compression=none size=1 mtime=1.0 hash=md5:x\n"
        "chr1\t0\t100\t100\t101\t5\t0\t100\n"
        "chr1\t100\t50\t50\t51\t5\t100\t150\n"
        "chr2\t0\t80\t80\t81\t5\t150\t230\n"
    )


def test_read_block_index_chroms(tmp_path):
    idx = tmp_path / "x.maf.block.idx"
    _write_idx(idx)
    assert intervals.read_block_index_chroms(str(idx)) == {"chr1", "chr2"}


def test_assert_bed_chroms_in_index_passes_when_names_match(tmp_path):
    idx = tmp_path / "x.maf.block.idx"; _write_idx(idx)
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t0\t50\nchr2\t10\t20\n")  # MAF names -> present in index
    intervals.assert_bed_chroms_in_index(str(bed), str(idx), "test_rule")  # no raise


def test_assert_bed_chroms_in_index_raises_on_prefix_mismatch(tmp_path):
    idx = tmp_path / "x.maf.block.idx"; _write_idx(idx)
    bed = tmp_path / "regions.bed"
    bed.write_text("1\t0\t50\n")  # bare name vs "chr1" in the index -> mismatch
    with pytest.raises(ValueError) as exc:
        intervals.assert_bed_chroms_in_index(str(bed), str(idx), "test_rule")
    assert "chromosome-name mismatch" in str(exc.value)
    assert "maf_prefix" in str(exc.value)


#############################################################################
# maf_symlink_target_for_group

def _write_scaffold_idx(path, scaffolds, compression="none"):
    # scaffolds: list of (name, size). Mirrors the mafutils .scaffold.idx format.
    lines = [f"# mafutils-index format=2 maf=x.maf compression={compression} size=999\n"]
    off = 16
    for name, size in scaffolds:
        lines.append(f"{name}\t{off}\t{size}\n")
        off += size
    path.write_text("".join(lines))


def test_maf_symlink_single_uncompressed_scaffold_matches(tmp_path):
    idx = tmp_path / "chr1.maf.scaffold.idx"
    _write_scaffold_idx(idx, [("chr1", 5000)])
    bed = tmp_path / "autosomes.bed"
    bed.write_text("chr1\t0\t4999\tchr1\n")  # col4 = output basename
    assert intervals.maf_symlink_target_for_group(str(idx), str(bed)) == "chr1"


def test_maf_symlink_none_when_multi_scaffold_input(tmp_path):
    idx = tmp_path / "genome.maf.scaffold.idx"
    _write_scaffold_idx(idx, [("chr1", 5000), ("chr2", 4000)])
    bed = tmp_path / "autosomes.bed"
    bed.write_text("chr1\t0\t4999\tchr1\n")
    assert intervals.maf_symlink_target_for_group(str(idx), str(bed)) is None


def test_maf_symlink_none_when_compressed(tmp_path):
    idx = tmp_path / "chr1.maf.gz.scaffold.idx"
    _write_scaffold_idx(idx, [("chr1", 5000)], compression="gzip")
    bed = tmp_path / "autosomes.bed"
    bed.write_text("chr1\t0\t4999\tchr1\n")
    assert intervals.maf_symlink_target_for_group(str(idx), str(bed)) is None


def test_maf_symlink_none_when_scaffold_name_mismatch(tmp_path):
    idx = tmp_path / "chr1.maf.scaffold.idx"
    _write_scaffold_idx(idx, [("chr1", 5000)])
    bed = tmp_path / "autosomes.bed"
    bed.write_text("chr2\t0\t4999\tchr2\n")  # group wants chr2, input has chr1
    assert intervals.maf_symlink_target_for_group(str(idx), str(bed)) is None


def test_maf_symlink_none_when_group_has_multiple_chroms(tmp_path):
    idx = tmp_path / "chr1.maf.scaffold.idx"
    _write_scaffold_idx(idx, [("chr1", 5000)])
    bed = tmp_path / "autosomes.bed"
    bed.write_text("chr1\t0\t4999\tchr1\nchr2\t0\t100\tchr2\n")  # >1 requested
    assert intervals.maf_symlink_target_for_group(str(idx), str(bed)) is None
