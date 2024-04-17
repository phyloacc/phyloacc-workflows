#############################################################################
# Given a MAF file, this script will create an index telling the location of
# each alignment block in the file
#
# Gregg Thomas, December 2023
#
# Output format for the mdx file is:
# 1. Reference scaffold
# 2. Reference position (0-based)
# 3. Alignment block sequence length
# 4. Alignment block line length
# 5. Number of sequences in the alignment block (though this is incorrect right now -- need to filter empty lines)
# 6. Start byte position of the alignment block in the maf file
# 7. End byte position of the alignment block in the maf file
#############################################################################

import sys
import os

#############################################################################

def processMAFBlock(block):
    ref_seq = block[1].split();
    # The reference sequence is always the second line in a block

    ref_scaff = ref_seq[1].split(".")[1];
    line_len = str(len(block[1]));
    num_seqs = str(len(block) - 1);
    seq_len = str(len(ref_seq[6]));

    return [ref_scaff, ref_seq[2], seq_len, line_len, num_seqs];

#############################################################################

#maf_file, mdx_file = sys.argv[1:];
# Get inputs from command line

#maf_file = "/n/holyscratch01/informatics/gwct/test.maf";
#mdx_file = "test.mdx";
maf_file = "/n/holyscratch01/informatics/gwct/241-mammalian-2020v2b-mafSplit/chr22.00.maf";
mdx_file = "chr22.00.mdx";
# 533.92user 28.90system 9:23.90elapsed 99%CPU (0avgtext+0avgdata 16680maxresident)k
# 499720inputs+786136outputs (0major+2548minor)pagefaults 0swaps
# Test inputs

with open(mdx_file, "w") as mdx_stream:
    maf_stream = open(maf_file, "r");
    # Open the input maf file

    line = "init";
    # Initialize the line variable so the loop starts

    #first_block = True;
    block = [];

    while line != "":
        line = maf_stream.readline();
        # Read a line from the maf file
        
        if line.startswith("#"):
            continue;
        # Skip the blocks containing comments

        if line.startswith("a"):
        # If the line is a header line, process the previous block and start the next one

            header_line_len = len(line);

            if block:
                block_end = maf_stream.tell() - header_line_len;
                # Since the current header line is the one for the next block, we need to subtract its length to get the correct end
                # position for the previous block

                outline = processMAFBlock(block) + [str(block_start), str(block_end)];
                # Process the info in the maf block and add the start and end positions                

                mdx_stream.write("\t".join(outline) + "\n");
                # Write the index line for the current block to the mdx file
            # For all cases but the first block, process the previous block

            block_start = maf_stream.tell() - header_line_len;
            block = [line.strip()];
            # Start the next block and get its location in the file

        else:
            block.append(line.strip());
        # For lines that are not headers, add them to the current block
    ## End of the loop over lines in the maf file

    block_end = maf_stream.tell();
    # Since the last block doesn't have a header line, we need to get the end position of the last line in the file
    
    outline = processMAFBlock(block) + [str(block_start), str(block_end)];
    mdx_stream.write("\t".join(outline) + "\n"); 
    # Process the last block  

    maf_stream.close();
    # Close the input maf file
## Close the output file


