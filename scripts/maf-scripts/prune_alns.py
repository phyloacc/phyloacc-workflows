#############################################################################
# A script to prune alignments based on a list of species to remove
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

indir = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/08-conserved-elements-chr-fasta/chr1/80-40-dedup-241/";
outdir = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/08-conserved-elements-chr-fasta/chr1/80-40-dedup-241-pruned-97/";
prune_spec_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/summary-data/prune-tips.txt";
overwrite = True;
# Input options

##########

species_to_prune = open(prune_spec_file, "r").read().splitlines();
# Read the list of species to prune from the file

species_to_prune = [ ">" + spec for spec in species_to_prune ];
# Add the ">" to the beginning of each species name

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

for fasta_file in fasta_files:
    aln_file = os.path.join(indir, fasta_file);
    # Get the full path to the alignment file

    base, ext = os.path.splitext(fasta_file);
    new_filename = f"{base}-pruned{ext}";
    output_file = os.path.join(outdir, new_filename);
    # Get the full path to the output file by adding the dedup suffix to the base filename

    aln = fastaGetDict(aln_file);
    # Read the alignment into a dictionary

    aln = {spec: aln[spec] for spec in aln.keys() - species_to_prune}
    # Remove the species to prune from the alignment

    with open(output_file, "w") as out:
        for seq in aln:
            out.write(f"{seq}\n{aln[seq]}\n");
    # Write the final alignment to the output file