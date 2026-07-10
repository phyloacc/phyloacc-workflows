#############################################################################
# After running maf_fetch.py, this script will remove duplicate sequences from
# each fasta file in the output directory.
#
# Gregg Thomas, January 2024
#############################################################################

import sys
import os
from collections import defaultdict

#############################################################################

def fastaGetDict(i_name):
# fastaGetDict reads a FASTA file and returns a dictionary containing all sequences in the file with 
# the key:value format as title:sequence.

    seqdict = {};
    for line in open(i_name, "r"):
        if line == "\n":
            continue;
        line = line.replace("\n", "");
        if line[0] == '>':
            curkey = line;
            seqdict[curkey] = "";
        else:
            seqdict[curkey] = seqdict[curkey] + line;

    return seqdict;

#############################################################################

#indir, outdir = sys.argv[1:];

chrom = "chrX";

indir = f"/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/08-conserved-elements-chr-fasta/{chrom}/80-40/";
outdir = f"/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/08-conserved-elements-chr-fasta/{chrom}/80-40-dedup-241/";
outbed = f"/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/08-conserved-elements-chr-fasta/{chrom}/241-mammalian-2020v2b.phylop-conserved-windows.{chrom}.0.05.80.40.241.bed";
overwrite = True;
# Input options

##########

spec_delim = ".";
# The delimiter in the sequence name that separates the species name from the rest of the header

min_spec = 241;
# The minimum number of species that must be present in the alignment after removing duplicates in order for it to be written to the output directory
# If this is set to 0, all alignments will be written to the output directory.

only_spec_name = True;
# Whether or not to only use the species name in the output sequence header

stats = True;
# Whether or not to print the distribution of number of species per alignment after deduplication

## Processing options - these should be added as options to the script

##########

if not os.path.isdir(indir):
    print("ERROR: Input directory does not exist");
    sys.exit(1);
# Make sure the input directory exists

if os.path.isdir(outdir) and not overwrite:
    print("ERROR: Overwite not set");
    sys.exit(1);
# If the output directory exists and the overwrite flag isn't set, exit here

if not os.path.isdir(outdir):
    os.makedirs(outdir);
# If the output directory doesn't exist, create it

valid_extensions = [".fa", ".fasta", ".fna"]
fasta_files = [ f for f in os.listdir(indir) if any(f.endswith(ext) for ext in valid_extensions) ]
# Get a list of all fasta files in the input directory

if not fasta_files:
    print("ERROR: No fasta files found in input directory");
    sys.exit(1);
# If no fasta files were found in the input directory, exit here

##########

num_seqs = defaultdict(int);
# A dictionary to keep track of the number of sequences in each alignment

with open(outbed, "w") as outbed_stream:
    for fasta_file in fasta_files:
        aln_file = os.path.join(indir, fasta_file);
        # Get the full path to the alignment file

        base, ext = os.path.splitext(fasta_file);
        file_chrom, start, end = base.split("-")[1:4];
        end = end.replace(".fa", "");
        new_filename = f"{base}-dedup{ext}";
        output_file = os.path.join(outdir, new_filename);
        # Get the full path to the output file by adding the dedup suffix to the base filename

        aln = fastaGetDict(aln_file);
        # Read the alignment into a dictionary

        spec_seqs = defaultdict(dict);
        for seq in aln:
            spec = seq.split(spec_delim)[0];
            spec_seqs[spec][seq] = aln[seq];
        # Split the sequences into species-specific dictionaries

        final_aln = {};
        for spec in spec_seqs:
            if len(spec_seqs[spec]) == 1:
                final_aln.update(spec_seqs[spec]);
            # If there is only one sequence for this species, add it to the final alignment

            else:
                min_gaps_seq = None;
                min_gaps = 999999;
                # Initialize min variables

                for seq in spec_seqs[spec]:
                    gaps = spec_seqs[spec][seq].count("-");
                    if gaps < min_gaps:
                        min_gaps = gaps;
                        min_gaps_seq = seq;
                ## Loop through each sequence for this species and find the one with the fewest gaps

                final_aln[min_gaps_seq] = spec_seqs[spec][min_gaps_seq];
            # If there are multiple sequences for this species, add the one with the fewest gaps to the final alignment
        ## End species loop

        if min_spec and len(final_aln) != min_spec:
            continue;
        # If the minimum number of species is set and the final alignment doesn't have that many species, skip this alignment
        else:
            outbed_stream.write(f"{file_chrom}\t{start}\t{end}\n");
            with open(output_file, "w") as out:
                for seq in final_aln:
                    if only_spec_name:
                        out.write(f"{seq.split(spec_delim)[0]}\n{final_aln[seq]}\n");
                    else:
                        out.write(f"{seq}\n{final_aln[seq]}\n");
            # Write the final alignment to the output file
        ## End minimum species check

        num_seqs[len(final_aln)] += 1;
        # Add the number of sequences in this alignment to the dictionary

    ## End fasta file loop
## Close bed file

##########

if stats:
    for num in sorted(num_seqs):
        print(f"{num}\t{num_seqs[num]}");
# Print the number of alignments with each number of sequences
