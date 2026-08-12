#!/usr/bin/env python3
"""
make_group_beds.py

Generate a BED file that defines the coordinates (0-based) of each chromosome
in a specified group. Each group is defined by a list of chromosome names
provided as command-line arguments. The length of each chromosome is derived
from the MAF's own block index, so no reference FASTA/.fai is needed.

Usage:
    python make_group_beds.py maf_block_index.tsv output.bed chrom1 chrom2 ...

Arguments:
    maf_block_index.tsv : The .block.idx file produced by `mafutils index`
                           (columns: ref_scaff, ref_start, ref_len, seq_len,
                           line_len, num_seqs, byte_start, byte_end).
    output.bed           : Path to the BED file to write.
    chrom1 ...           : List of chromosome names to include (passed from Snakemake config).
"""

import sys

# Parse command-line arguments
maf_block_index = sys.argv[1]       # Path to the MAF's .block.idx (from mafutils index)
chrom_prefix = sys.argv[2]
output_bed = sys.argv[3]            # Output BED file to write
chrom_list = [str(c) for c in sys.argv[4:]]  # Remaining args are chromosome names (strings)

# Check that chromosome names were passed
if not chrom_list:
    sys.stderr.write("Error: No chromosomes provided.\n")
    sys.exit(1)

# Convert list to a set for fast lookup
chrom_set = set(chrom_list)

# A chromosome spans many block-index rows, so accumulate the furthest extent
# seen per scaffold in one pass rather than assuming one row = one chromosome.
chrom_max_end = {}
with open(maf_block_index) as in_f:
    for line in in_f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        chrom = parts[0]
        # A FASTA .fai lists bare chromosome names ("1"), but a MAF block index lists
        # them as they appear in the MAF - carrying maf_chr_prefix ("chr1"). Strip the
        # prefix so both input types match the (bare) config chromosome set; the output
        # below re-adds the prefix.
        if chrom_prefix and chrom.startswith(chrom_prefix):
            chrom = chrom[len(chrom_prefix):]
        if chrom not in chrom_set:
            continue
        ref_start = int(parts[1])
        ref_len = int(parts[2])
        end = ref_start + ref_len
        if end > chrom_max_end.get(chrom, 0):
            chrom_max_end[chrom] = end

missing = chrom_set - chrom_max_end.keys()
if missing:
    sys.stderr.write(f"Error: no blocks found for chromosome(s): {', '.join(sorted(missing))}\n")
    sys.exit(1)

# Open output BED file for writing
with open(output_bed, "w") as out_f:
    for chrom in chrom_list:
        out_chrom = chrom_prefix + chrom
        end = chrom_max_end[chrom] - 1  # BED format: 0-based, half-open
        out_f.write(f"{out_chrom}\t0\t{end}\t{out_chrom}\n")

# Optional: write summary message to stderr/log
sys.stderr.write(f"Wrote {output_bed}\n")# for chromosome group: {', '.join(sorted(chrom_set))}\n")