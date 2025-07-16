#!/usr/bin/env python3
"""
make_group_beds.py

Generate a BED file that defines the coordinates (0-based) of each chromosome
in a specified group. Each group is defined by a list of chromosome names
provided as command-line arguments. The length of each chromosome is obtained
from the reference index file.

Usage:
    python make_group_beds.py ref_index.tsv output.bed chrom1 chrom2 ...

Arguments:
    ref_index.tsv : A tab-delimited file with at least two columns:
                    chromosome name and chromosome length.
    output.bed    : Path to the BED file to write.
    chrom1 ...    : List of chromosome names to include (passed from Snakemake config).
"""

import sys

# Parse command-line arguments
ref_index = sys.argv[1]             # Path to reference index file (e.g. ref_genome.fa.fai)
chrom_prefix = sys.argv[2]
output_bed = sys.argv[3]            # Output BED file to write
chrom_list = [str(c) for c in sys.argv[4:]]  # Remaining args are chromosome names (strings)

# Check that chromosome names were passed
if not chrom_list:
    sys.stderr.write("Error: No chromosomes provided.\n")
    sys.exit(1)

# Convert list to a set for fast lookup
chrom_set = set(chrom_list)

# Open output BED file for writing
with open(output_bed, "w") as out_f:
    # Open reference index and process each line
    with open(ref_index) as in_f:
        for line in in_f:
            # Reference index expected to have at least: chrom <tab> length
            chrom, length = line.strip().split("\t")[:2]

            # Only write chroms that belong to the specified group
            if chrom in chrom_set:
                out_chrom = chrom_prefix + chrom
                end = int(length) - 1  # BED format: 0-based, half-open
                out_f.write(f"{out_chrom}\t0\t{end}\t{out_chrom}\n")

# Optional: write summary message to stderr/log
sys.stderr.write(f"Wrote {output_bed}\n")# for chromosome group: {', '.join(sorted(chrom_set))}\n")