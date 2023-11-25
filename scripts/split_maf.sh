#!/bin/bash
#############################################################################
# Split a MAF file into multiple MAF files by chromosome
# Assumes the first sequence listed in each block is the reference and names
# output files by the reference name
#
# Usage: ./split_maf.sh <maf_file> <output_dir>
#
# Gregg Thomas, October 2023
#############################################################################

filename=$1
outdir=$2
# Inputs, input file should be a MAF file
# Note that for now, output files must be erased between runs if the
# run fails...

mkdir -p $outdir
# Make the output directory if it doesn't exist

first=true
# Boolean to indicate whether or not we are on the first block

lines=()
# Stores the lines of the current block

block_line="a"
# The line that indicates the start of a new block

#i=0
# A counter for testing

while read line; do
    
    if [[ "$line" == "$block_line" ]]; then
    # Check if the line is the start of a new block
        
        if [[ "$first" == true ]]; then
        # Check if this is the first block

            first=false
            # All subsequent blocks will not be the first

            first_block=true
            # Boolean to indicate we've reached the first line in a block and need to parse the output name

            lines+=("$line")
            # Add the current line ("a") to the lines array

        else
        # If this is not the first block, we need to output the lines from the current block to a file

            for outline in "${lines[@]}"; do
                echo "$outline" >> "$outdir/$output_name"
            done
            # Output the lines to a file based on the reference name in the block

            first_block=true
            # Boolean to indicate we've reached the first line in the next block and need to parse the output name

            lines=("$line")
            # Reset the lines array for the next block, adding the current line ("a")

            #i=$((i+1))
            # if [[ $i -eq 20 ]]; then
            #     break
            # fi
            # For testing, only run the first 20 blocks

        fi
        continue
        # Skip to the next line since we've already added the "a" line to the lines array
    fi

    if [[ "$first_block" == true ]]; then
        output_name=$(echo "$line" | cut -f2)
        output_name="$output_name.maf"
        first_block=false
    fi
    # If this is the first line in a block, parse the output file name from the reference name

    lines+=("$line")
    # Add the current line to the lines array

done < "$filename"
# Read the input file line by line

#############################################################################
