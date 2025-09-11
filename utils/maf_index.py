#############################################################################
# Given a MAF file, this script will create an index telling the location of
# each alignment block in the file
#
# Gregg Thomas, December 2023
#
# Output format for the mdx file is:
# 1. Reference scaffold
# 2. Reference position (0-based)
# 3. Reference interval length
# 4. Alignment block sequence length
# 5. Alignment block line length
# 6. Number of sequences in the alignment block
# 7. Start byte position of the alignment block in the maf file
# 8. End byte position of the alignment block in the maf file
#############################################################################

import sys
import os
import gzip

import lib.common as COMMON

#############################################################################

def processMAFBlock(block):
    ref_seq = block[1].split()  # Reference line (second line in block)
    ref_scaff = ref_seq[1].split(".", 1)[1]
    ref_len = block[1].split()[3]
    line_len = str(len(block[1]))
    num_seqs = str(len(block) - 1)
    seq_len = str(len(ref_seq[6]))  # Aligned sequence

    return [ref_scaff, ref_seq[2], ref_len, seq_len, line_len, num_seqs]

#############################################################################

if len(sys.argv) != 4:
    print("Usage: python index_maf.py <input.maf> <output.block.idx> <output.scaffold.idx>")
    sys.exit(1)

maf_file, mdx_block_file, mdx_scaffold_file = sys.argv[1:]

maf_compression = COMMON.detectCompression(maf_file)

if maf_compression == "gz":
    maf_stream = gzip.open(maf_file, "rt")
elif maf_compression == "none":
    maf_stream = open(maf_file, "r")
else:
    raise ValueError(f"Unsupported compression type: {maf_compression}")

#maf_file = "/n/holyscratch01/informatics/gwct/test.maf";
#mdx_file = "test.mdx";
#maf_file = "/n/holyscratch01/informatics/gwct/241-mammalian-2020v2b-mafSplit/all-chromosomes/chr22.00.maf";
#mdx_file = "test.mdx";
# 533.92user 28.90system 9:23.90elapsed 99%CPU (0avgtext+0avgdata 16680maxresident)k
# 499720inputs+786136outputs (0major+2548minor)pagefaults 0swaps
# Test inputs

with open(mdx_block_file, "w") as block_stream, open(mdx_scaffold_file, "w") as scaffold_stream:
    line = "init"
    block = []
    block_num = 0

    # Track region info for rdx
    current_scaffold = None
    region_start_byte = None
    region_end_byte = None

    while line != "":
        line = maf_stream.readline()
        if line.startswith("#") or line.strip() == "":
            continue

        if line.startswith("a"):
            header_line_len = len(line)

            if block:
                block_num += 1
                block_end = maf_stream.tell() - header_line_len

                # Process block info
                block_info = processMAFBlock(block)
                ref_scaffold = block_info[0]

                # Write block index
                mdx_line = block_info + [str(block_start), str(block_end)]
                block_stream.write("\t".join(mdx_line) + "\n")

                # Handle scaffold grouping for region index
                if current_scaffold is None:
                    # First block
                    current_scaffold = ref_scaffold
                    region_start_byte = block_start
                    region_end_byte = block_end
                elif ref_scaffold != current_scaffold:
                    # Scaffold changed — flush the previous region
                    scaffold_stream.write(f"{current_scaffold}\t{region_start_byte}\t{region_end_byte}\n")
                    current_scaffold = ref_scaffold
                    region_start_byte = block_start
                    region_end_byte = block_end
                else:
                    # Same scaffold — extend
                    region_end_byte = block_end

            block_start = maf_stream.tell() - header_line_len
            block = [line.strip()]
        else:
            block.append(line.strip())

    # Handle the last block
    if block:
        block_end = maf_stream.tell()
        block_info = processMAFBlock(block)
        ref_scaffold = block_info[0]

        mdx_line = block_info + [str(block_start), str(block_end)]
        block_stream.write("\t".join(mdx_line) + "\n")

        if current_scaffold != ref_scaffold:
            # First flush previous scaffold
            if current_scaffold is not None:
                scaffold_stream.write(f"{current_scaffold}\t{region_start_byte}\t{region_end_byte}\n")
            # Start new region
            current_scaffold = ref_scaffold
            region_start_byte = block_start
            region_end_byte = block_end
        else:
            # Same scaffold, extend it
            region_end_byte = block_end

        # Final block — always write
        scaffold_stream.write(f"{current_scaffold}\t{region_start_byte}\t{region_end_byte}\n")

    maf_stream.close()
    # Close the input maf file
## Close the output file


