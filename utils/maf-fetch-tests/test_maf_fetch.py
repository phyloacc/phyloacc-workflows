import os
import subprocess
import pytest

TEST_DIR = os.path.dirname(__file__)
BED_FILE = os.path.join(TEST_DIR, "example.bed")
#BED_FILE = os.path.join(TEST_DIR, "crossblocks-missing-species.bed")
MAF_FILE = os.path.join(TEST_DIR, "example.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")
OUT_DIR = os.path.join(TEST_DIR, "results")
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
def test_maf_fetch_region(chrom, start, end, region_id):
    """Run maf_fetch for just one region and compare output to expected."""
    # Write a one-line BED for this region
    per_test_bed = os.path.join(OUT_DIR, f"{region_id}.bed")
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(per_test_bed, "w") as out:
        out.write(f"{chrom}\t{start}\t{end}\t{region_id}\n")

    # Run maf_fetch
    output_file = os.path.join(OUT_DIR, f"{region_id}.maf")
    cmd = [
        "python", "../maf_fetch.py",
        MAF_FILE,
        INDEX_FILE,
        per_test_bed,
        "-b", "id",
        "-o", OUT_DIR
    ]
    subprocess.run(cmd, check=True)

    got = open(output_file).read()
    expected_file = os.path.join(EXPECTED_MAF_DIR, f"{region_id}.maf")
    expected = open(expected_file).read()
    
    # Strip trailing blank lines from both outputs before comparison
    got = got.rstrip('\n')
    expected = expected.rstrip('\n')

    assert got == expected, f"Output for region {region_id} does not match expected"

@pytest.mark.parametrize("chrom,start,end,region_id", list(get_regions()))
def test_maf_fetch_fasta_region(chrom, start, end, region_id):
    """Run maf_fetch for just one region (FASTA output) and compare to expected."""
    per_test_bed = os.path.join(OUT_DIR, f"{region_id}.bed")
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(per_test_bed, "w") as out:
        out.write(f"{chrom}\t{start}\t{end}\t{region_id}\n")

    # Run maf_fetch with --fasta/-f
    output_file = os.path.join(OUT_DIR, f"{region_id}.fa")
    cmd = [
        "python", "../maf_fetch.py",
        MAF_FILE,
        INDEX_FILE,
        per_test_bed,
        "-b", "id",
        "-o", OUT_DIR,
        "-f"
    ]
    subprocess.run(cmd, check=True)

    got = open(output_file).read()
    expected_file = os.path.join(EXPECTED_FASTA_DIR, f"{region_id}.fa")
    expected = open(expected_file).read()

    # Strip trailing blank lines from both outputs before comparison
    got = got.rstrip('\n')
    expected = expected.rstrip('\n')

    assert got == expected, f"FASTA output for region {region_id} does not match expected"
