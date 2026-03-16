import os
import subprocess
import pytest

TEST_DIR = os.path.dirname(__file__)
SCRIPT_FILE = os.path.join(TEST_DIR, "..", "maf_fetch.py")
BED_FILE = os.path.join(TEST_DIR, "example.bed")
#BED_FILE = os.path.join(TEST_DIR, "crossblocks-missing-species.bed")
MAF_FILE = os.path.join(TEST_DIR, "example.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")
EXPECTED_MAF_DIR = os.path.join(TEST_DIR, "expected-maf")
EXPECTED_FASTA_DIR = os.path.join(TEST_DIR, "expected-fasta")

def get_regions():
    """Yields (chrom, start, end, id) for each region in the BED file."""
    with open(BED_FILE) as bed:
        for line in bed:
            if line.strip() and not line.startswith("#"):
                fields = line.strip().split()
                yield fields[0], int(fields[1]), int(fields[2]), fields[3]

@pytest.mark.parametrize("chrom,start,end,region_id", list(get_regions()))
def test_maf_fetch_region(chrom, start, end, region_id, tmp_path):
    """Run maf_fetch for just one region and compare output to expected."""
    out_dir = str(tmp_path)
    # Write a one-line BED for this region
    per_test_bed = os.path.join(out_dir, f"{region_id}.bed")
    os.makedirs(out_dir, exist_ok=True)
    with open(per_test_bed, "w") as out:
        out.write(f"{chrom}\t{start}\t{end}\t{region_id}\n")

    # Run maf_fetch
    output_file = os.path.join(out_dir, f"{region_id}.maf")
    cmd = [
        "python", SCRIPT_FILE,
        MAF_FILE,
        INDEX_FILE,
        per_test_bed,
        "-b", "id",
        "-o", out_dir
    ]
    expected_file = os.path.join(EXPECTED_MAF_DIR, f"{region_id}.maf")
    result = subprocess.run(cmd, check=False)

    if os.path.getsize(expected_file) == 0:
        assert result.returncode != 0, f"Expected non-zero exit for region {region_id}"
        assert not os.path.exists(output_file), f"Did not expect output file for region {region_id}"
        return

    assert result.returncode == 0, f"Unexpected non-zero exit for region {region_id}: {result.returncode}"

    got = open(output_file).read()
    expected = open(expected_file).read()
    
    # Strip trailing blank lines from both outputs before comparison
    got = got.rstrip('\n')
    expected = expected.rstrip('\n')

    assert got == expected, f"Output for region {region_id} does not match expected"

@pytest.mark.parametrize("chrom,start,end,region_id", list(get_regions()))
def test_maf_fetch_fasta_region(chrom, start, end, region_id, tmp_path):
    """Run maf_fetch for just one region (FASTA output) and compare to expected."""
    out_dir = str(tmp_path)
    per_test_bed = os.path.join(out_dir, f"{region_id}.bed")
    os.makedirs(out_dir, exist_ok=True)
    with open(per_test_bed, "w") as out:
        out.write(f"{chrom}\t{start}\t{end}\t{region_id}\n")

    # Run maf_fetch with --fasta/-f
    output_file = os.path.join(out_dir, f"{region_id}.fa")
    cmd = [
        "python", SCRIPT_FILE,
        MAF_FILE,
        INDEX_FILE,
        per_test_bed,
        "-b", "id",
        "-o", out_dir,
        "-f",
        "-fh", "species-coords-id"
    ]
    expected_file = os.path.join(EXPECTED_FASTA_DIR, f"{region_id}.fa")
    result = subprocess.run(cmd, check=False)

    if os.path.getsize(expected_file) == 0:
        assert result.returncode != 0, f"Expected non-zero FASTA exit for region {region_id}"
        assert not os.path.exists(output_file), f"Did not expect FASTA output file for region {region_id}"
        return

    assert result.returncode == 0, f"Unexpected non-zero FASTA exit for region {region_id}: {result.returncode}"

    got = open(output_file).read()
    expected = open(expected_file).read()

    # Strip trailing blank lines from both outputs before comparison
    got = got.rstrip('\n')
    expected = expected.rstrip('\n')

    assert got == expected, f"FASTA output for region {region_id} does not match expected"
