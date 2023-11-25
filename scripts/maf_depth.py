#############################################################################
# Given a MAF file, this script will create a wig file of depth per site
#
# Gregg Thomas, October 2023
#############################################################################

import sys
import os
import re
from itertools import groupby

#############################################################################

maf_file, chrome, wig_file = sys.argv[1:];
# Get inputs from command line

# maf_file = "/n/holyscratch01/informatics/gwct/241-mammalian-2020v2b-mafSplit-1Mb/chr1.69.maf";
# chrome = "chr17";
# wig_file = "test.wig";
# Test inputs

with open(maf_file, "r") as maf_stream, open(wig_file, "w") as out_stream:

    maf_iter = (x[1] for x in groupby(maf_stream, lambda line: line.startswith("a")));
    # Open the input MAF file and create an iterator that will return blocks of lines
    # This reads one block in a MAF file at a time

    readstr = lambda s : s.strip();
    # A function to strip whitespace from strings

    # num_blocks = 0;
    # num_neg_blocks = 0;
    for header_obj in maf_iter:
        header = readstr(header_obj.__next__());
        # The header object is an iterator. This gets the string.

        if header.startswith("#"):
           continue;
        # Skip the blocks containing comments

        block = [readstr(s) for s in maf_iter.__next__()];
        #num_blocks += 1;
        # The current header should correspond to the current iterator in maf_iter. This gets all those
        # lines and combines them as a string.
        
        # print(header);
        # print("\n".join(block));

        cur_pos = int(block[0].split()[2]);
        # Get the current position in the alignment. MAF files are 0-based, so it matches the wig/bed format we're after already
        # Negative strand MAF starts may also need to be adjusted (https://genome.ucsc.edu/FAQ/FAQformat.html#format5), but I'm not
        # seeing any negative strand reference blocks in the MAF files I'm working with

        specs, seqs = [], [];
        for row in block:
            if row.startswith("s"):
                spec = row.split()[1].split(".")[0];
                # QUESTION: Is "." always the delimiter between species id and scaffold id?
                # Get the species id of the current sequence

                if spec not in specs:
                    seqs.append(row.split()[6]);
                # Add the sequence to the list if that species hasn't been seen yet in this block
                ## TODO: Better way to deal with paralogs?

                specs.append(spec);
                # Add the species to the list of species seen in this block
        # Get the sequences for each species in the block

        sites = [''.join(s) for s in zip(*seqs)];
        # Transpose the sequences to get the sites

        for site in sites:
            depth = len(re.sub(r'[-nN]', '', site));
            # Get the depth for the current site. This is the length of the site minus the number of gaps

            outline = [chrome, str(cur_pos), str(cur_pos+1), str(depth)];
            out_stream.write("\t".join(outline) + "\n");
            cur_pos += 1;
            # Write the depth for each site to the output file
        # Write the depth for each site to the output file

        
        # if num_blocks > 10000:
        #     break;
        # Some debug code
## Close the input and output files

# print(num_blocks);
# print(num_neg_blocks);