#############################################################################
# Given a bed file, maf file, and maf index, this script will retrieve the
# alignments from the regions in the bed file from the maf file.
#
# Gregg Thomas, December 2023
#############################################################################

## DO MAF FILES EVER HAVE GAPS IN THE REFERENCE SEQ?
## DO MAF FILES CONTAIN EVERY POSITION OF THE REFERENCE SEQ? OTHER SEQS?

## MAF FILES ARE 0-BASED
## https://genome.ucsc.edu/FAQ/FAQformat.html#format5

import sys
import os
import re
#from itertools import groupby
from collections import defaultdict
import subprocess

#############################################################################

def seqTest(genome_fasta, maf_fetch_seq, chrom, start, end, debug=False):
    ref_chrom = chrom.replace("chr", "");
    command = f"samtools faidx {genome_fasta} {ref_chrom}:{start+1}-{end}"
    output = "init";

    try:
        # Run the command and capture the output
        output = subprocess.check_output(command, shell=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools faidx: {e}")

    output_list = output.split("\n")[1:];
    genome_seq = "".join(output_list);

    if genome_seq.upper() == maf_fetch_seq.upper():
        if debug:
            print("\n-----");
            print("SUCCESS: MAF_FETCH SEQ MATCHES REFERENCE GENOME SEQ");
            print(output.split("\n")[0]);
            print(genome_seq.upper());
            print(">maf_fetch_seq " + chrom + ":" + str(start) + "-" + str(end));
            print(maf_fetch_seq.upper());            
            print("-----");
        return True;
    else:
        print("\n-----");
        print("ERROR: MAF_FETCH SEQ DOES NOT MATCH REFERENCE GENOME SEQ");
        print(output.split("\n")[0]);
        print(genome_seq.upper());
        print(">maf_fetch_seq " + chrom + ":" + str(start) + "-" + str(end));
        print(maf_fetch_seq.upper());
        print("-----");
        return False;

#############################################################################

maf_file, mdx_file, ref_spec, bed_file, outdir = sys.argv[1:];
# Get inputs from command line

# ref_genome_fasta = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/00-human-ref-ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# maf_file = "/n/holyscratch01/informatics/gwct/241-mammalian-2020v2b-mafSplit/chr1.00.maf";
# mdx_file = "/n/holyscratch01/informatics/gwct/241-mammalian-2020v2b-mafSplit/chr1.00.mdx"
# bed_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/07-conserved-elements-chr/chr1/80/241-mammalian-2020v2b.phylop-conserved-windows.chr1.0.05.80.40.bed";
# outdir = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/08-conserved-elements-chr-fasta/chr1/80-40/"


#bed_file = "test-chr22.bed";
#mdx_file = "chr22.00.mdx";
#outdir = "test-out-100";
overwrite = True;

# ref_spec = "Homo_sapiens";
# ref_chr = "chr1";
#query_spec = "Castor_canadensis";
#query_scaff = "CasCan_scaffold_576";
# query_spec = "Solenodon_paradoxus";
# query_scaff = "SolPar_flattened_line_9386";

# bed_sorted = True;
# Assume the bed file is sorted... i think this should be default since this is easy to do with bedtools

debug_level = 0;
test = False;
# Test inputs

if os.path.isdir(outdir) and not overwrite:
    print("ERROR: Overwite not set");
    sys.exit(1);

if not os.path.isdir(outdir):
    os.makedirs(outdir);

maf_stream = open(maf_file, "r");
# Open the maf file

block_iter = iter(open(mdx_file));
current_block = next(block_iter);
# Open the maf index (mdx) file and read the first block
# Doing this here ensures we don't need a nested loop to read the mdx file
# from the beginning for each bed region

for line in open(bed_file):

    if debug_level >= 1 or debug_level == -1:
        print("\n--------------------");
        print("REGION: ", line.strip());
    # Debug stuff

    chrom, region_start, region_end = line.strip().split()[:3];
    region_start = int(region_start);
    region_end = int(region_end);
    region_len = region_end - region_start;
    # Read the info from the line in the bed file

    if debug_level >= 2:
        print("REGION START:", region_start);
        print("REGION END:", region_end);
        print("REGION LEN:", region_len);
    # Debug stuff

    maf_start = region_start;
    # This is the starting position relative to what we've read in the maf file
    # This will be updated to be the end position of the last block we read so 
    # that we can continue getting sequence from that position

    remaining_region_len = region_end - region_start;
    # This is the length of the region we still need to get sequence for

    seqs = defaultdict(dict);
    # This will be a dictionary of species and their sequences for this region

    first_block = True;
    # This is used to keep track of whether we're on the first block of the region

    outfilename = os.path.join(outdir, ref_spec + "-" + chrom + "-" + str(region_start) + "-" + str(region_end) + ".fa");

    ## Bed region stuff
    ####################

    while remaining_region_len > 0:
    # Continue getting sequence until we've gotten all the sequence for this region

        block_scaff, block_start_pos, block_aln_len, block_line_len, block_num_seqs, block_start_byte, block_end_byte = current_block.strip().split();
        # Read the info from the current block from the mdx file

        if chrom != block_scaff:
            continue;
        # If the block is not on the same chromosome as the region, skip it

        block_start_pos = int(block_start_pos);
        block_aln_len = int(block_aln_len);
        block_start_byte = int(block_start_byte);
        block_end_byte = int(block_end_byte);
        block_end_pos = block_start_pos + block_aln_len;
        # Convert the info from the mdx file to integers

        if block_start_pos <= maf_start and block_end_pos > maf_start:
        # If the block we've read from the mdx file overlaps with the region we're getting sequence for

            if debug_level >= 2:
                print("BLOCK: ", block_start_pos, block_end_pos, region_start, region_end);
                print("BLOCK START POS:", block_start_pos);
                print("BLOCK END POS:", block_end_pos);
                print("BLOCK ALN LEN:", block_aln_len);                    
                # print(block_start_byte, block_end_byte);
            # Debug stuff

            maf_stream.seek(block_start_byte);
            block = maf_stream.read(block_end_byte - block_start_byte).split("\n");
            block = list(filter(None, block));
            # Read the block from the maf file

            if debug_level == -1:
                print("\n".join(block));

            first_spec = True;
            block_specs = [];
            # Variables for readin through blocks:
            # A tracker for whether we're on the first species in the block (which should be the reference species)
            # A list of the species in the block

            for spec_entry in block[1:]:

                spec_entry = spec_entry.split();
                # Parse the maf entry for this species
                
                spec, scaff = spec_entry[1].split(".", 1);
                block_specs.append(spec_entry[1]);
                # Read the current species and add it to the list of species in the block

                block_seq = spec_entry[6];
                block_pos = int(spec_entry[2]);
                block_len = int(spec_entry[3]);
                block_strand = spec_entry[4];
                # Get the sequence and length of the block for this species

                # if block_strand == "-":
                #     block_pos = block_pos - block_len;

                last_pos = int(block_pos) + len(block_seq.replace("-", ""));
                # Get the last position of the sequence for this species, excluding alignment gaps

                if debug_level == 2 and (spec == ref_spec or spec == query_spec and scaff == query_scaff) or debug_level > 2:
                    print("ENTRY: ", spec_entry);
                    print("ENTRY SPEC:", spec);
                    print("ENTRY SCAFF:", scaff);
                    print("ENTRY POS:", block_pos);
                    print("ENTRY LEN:", block_len);
                    print("ENTRY LAST POS:", last_pos);
                    print("REMAINING REGION LEN:", remaining_region_len);
                    #print("BLOCK SEQ:", block_seq);
                # Debug stuff

                if first_spec:
                    if test:
                        assert block_seq.count("-") == 0, "FOUND GAPS IN REFERENCE MAF SEQUENCE: " + block[0];

                    seq_start_index = maf_start - block_start_pos;
                    # The starting index for the sequence to get

                    if block_len < remaining_region_len:
                        seq_end_index = block_len;
                        #remaining_region_len -= seq_end_index - seq_start_index;
                    # If the length of the current block is less than what we need to get                        

                    elif block_len == remaining_region_len:
                        seq_end_index = block_len;
                        #remaining_region_len -= block_len;
                    # If the length of the block exactly equals what we need to get

                    elif block_len > remaining_region_len:
                        seq_end_index = seq_start_index + remaining_region_len;
                        if seq_end_index > block_len:
                            seq_end_index = block_len;
                        #remaining_region_len -= seq_end_index - seq_start_index;
                    # If the length of the block is greater than what we need to get
                    
                    first_spec = False;
                    # We're no longer on the first species in the block
                # For the first species in the block, get the start and end index for the sequence
                # and calculate the remaining region length based on them
                ## End first species (reference) check block        

                seq_to_add = block_seq[seq_start_index:seq_end_index];
                # Get the sequence to add based on the start and end indexes calculated from the reference sequence

                if debug_level == 2 and (spec == query_spec and scaff == query_scaff) or debug_level > 2:
                    print("MAF START:", maf_start);
                    print("SEQ START INDEX:", seq_start_index);
                    print("SEQ END INDEX:", seq_end_index);
                    print("LEN SEQ TO ADD:", len(seq_to_add));
                # Debug stuff       

                if spec not in seqs:
                    # print("ADDING SPECIES");
                    seqs[spec] = defaultdict(list);
                # If the species isn't in our dict for this region, add it

                if scaff not in seqs[spec]:
                    if debug_level == 2 and (spec == query_spec and scaff == query_scaff) or debug_level > 2:
                        print("SCAFF NOT FOUND - ADDING SCAFFOLD AND SEQ");
                    seqs[spec][scaff].append({'scaff' : scaff, 'start' : int(spec_entry[2]) + seq_start_index, 'seq' : "", 'last-pos' : int(block_pos), 'found' : False});
                    seqs[spec][scaff][0]['seq'] = "-" * (region_len - remaining_region_len);
                # If the scaffold isn't in our dict for this species, add it and fill in the sequence up to this point with gaps
                # If this is the first block for the region, no gaps should be added (region_len - remaining_region_len = 0)

                seq_added = False;
                if debug_level == 2 and (spec == query_spec and scaff == query_scaff) or debug_level > 2:
                    print("SEQ ADDED?", seq_added);
                    print(seq_to_add);
                    print(seqs[spec][scaff]);
                    print("-----");
                # Debug stuff

                for i in range(len(seqs[spec][scaff])):

                    if abs(seqs[spec][scaff][i]['last-pos'] - block_pos) < 100 and not seqs[spec][scaff][i]['found']:
                        if debug_level == 2 and (spec == query_spec and scaff == query_scaff) or debug_level > 2:
                            print("ADDING TO PREVIOUS SEQ");
                            print("LENGTH OF SEQUENCE BEFORE ADDING:", len(seqs[spec][scaff][i]['seq']));
                            print("PREV LAST POS", seqs[spec][scaff][i]['last-pos']);
                            print(abs(seqs[spec][scaff][i]['last-pos'] - block_pos));
                        # Debug stuff

                        seqs[spec][scaff][i]['seq'] += seq_to_add;
                        seqs[spec][scaff][i]['last-pos'] = last_pos;
                        seqs[spec][scaff][i]['found'] = True;
                        seq_added = True;
                        break;
                    # If there exists a sequence on the same scaffold for this species that is within 100bp of the last added position,
                    # add the sequence to that sequence and update the last position of the sequence
                ## End sequence adding loop

                if debug_level == 2 and (spec == query_spec and scaff == query_scaff) or debug_level > 2:
                    print("SEQ ADDED?", seq_added);
                    print(seq_to_add);
                    print("-----");
                # Debug stuff

                if not seq_added:
                    if debug_level == 2 and (spec == query_spec and scaff == query_scaff) or debug_level > 2:
                        print("SEQ NOT FOUND - ADDING NEW SEQ");
                    seqs[spec][scaff].append({'scaff' : scaff, 'start' : int(spec_entry[2]) + seq_start_index, 'seq' : "", 'last-pos' : int(block_pos), 'found' : True});
                    seqs[spec][scaff][-1]['seq'] = "-" * (region_len - remaining_region_len);
                    seqs[spec][scaff][-1]['seq'] += seq_to_add;
                # If the sequence hasn't been added to any existing sequence for this species, add a new sequence for this species
                # Fill in with gaps up to this point, then add the sequence
            ## End block loop

            first_block = False;

            remaining_region_len -= seq_end_index - seq_start_index;

            if debug_level >= 2:
                print("REMAINING REGION LEN:", remaining_region_len);

            #if debug_level >= 2:
            if debug_level >= 2:
            #     print(seqs[spec][scaff]);
            #     print("SEQ:", seqs[spec][scaff][0]['seq']);
            # elif debug_level > 2:
                for spec in seqs:
                    for scaff in seqs[spec]:
                        if spec == query_spec and scaff == query_scaff:
                            print(seqs[spec][scaff]);
                            for i in range(len(seqs[spec][scaff])):
                                print("SPEC:", spec);
                                print("SCAFF:", scaff);
                                print("SEQ:", seqs[spec][scaff][i]['seq']);
                                print("=====");
                print("-----------------");
                #sys.exit();
                #print(remaining_region_len);
            # Debug stuff

            for spec in seqs:
                for scaff in seqs[spec]:
                    for i in range(len(seqs[spec][scaff])):
                        if not seqs[spec][scaff][i]['found']:
                            seqs[spec][scaff][i]['seq'] += "-" * (seq_end_index - seq_start_index);
                        seqs[spec][scaff][i]['found'] = False;
            # For any species in our dict that wasn't in this block, add gaps

            maf_start = block_end_pos;
            # Update the starting index of the maf to be the end position of this block

            if remaining_region_len == 0:
                break;
            # End the MAF index loop if we've gotten all the sequence for this region
        ## End region overlap block

        try:
            current_block = next(block_iter);
        except StopIteration:
            break;
        # Read the next block from the mdx file, or break if its this was the last block
        # TODO: probably some check to make sure the regions in the bed file are in the mdx file

    ## End MAF index loop

    if test and not seqTest(ref_genome_fasta, seqs["Homo_sapiens"][ref_chr][0]['seq'], chrom, region_start, region_end, debug=True):
        print(outfilename);
        sys.exit(1);
    # A test to make sure the reference sequence we got matches what samtools gets from the genome fasta file


    with open(outfilename, "w") as out_stream:
        length_error = False;
        for spec in seqs:
            for scaff in seqs[spec]:
                for i in range(len(seqs[spec][scaff])):
                    if test:
                        if len(seqs[spec][scaff][i]['seq']) != region_len:
                            length_error = True;

                    spec_end = str(int(seqs[spec][scaff][i]['start']) + len(seqs[spec][scaff][i]['seq'].replace("-", "")));
                    header = ">" + spec + "." + scaff + ":" + str(seqs[spec][scaff][i]['start']) + "-" + spec_end;
                    out_stream.write(header + "\n");
                    out_stream.write(seqs[spec][scaff][i]['seq'] + "\n");

    if length_error:
        print("\n-----");
        print("ERROR: SOME RETRIEVED SEQUENCE LENGTHS DO NOT MATCH REGION LENGTH");
        print(outfilename);
        print("-----");
        sys.exit(1);

    if test:
        ref_gap_count = seqs[ref_spec][ref_chr][0]['seq'].count("-");
        if ref_gap_count != 0:
            print("\n-----");
            print("CHECK: REFERENCE SEQUENCE CONTAINS GAPS");
            print(outfilename);
            print("-----");
            sys.exit(1);
    ## End FASTA output loop
## End bed region loop

maf_stream.close();

#############################################################################
## STASH

## This was code for an unsorted bed file. incomplete.
# else:
#     for line in open(bed_file):
#         if line.startswith("#"):
#             continue;

#         #if debug:
#         print("REGION: ", line);

#         chrom, region_start, region_end = line.strip().split()[:3];
#         region_start = int(region_start);
#         maf_start = region_start;
#         region_end = int(region_end);
#         remaining_region_len = region_end - region_start;
#         seq = "";
#         get_seq = False;

#         for block in open(mdx_file):
#             block_scaff, block_start_pos, block_aln_len, block_line_len, block_num_seqs, block_start_byte, block_end_byte = block.strip().split();
#             if chrom != block_scaff:
#                 continue;

#             block_start_pos = int(block_start_pos);
#             block_aln_len = int(block_aln_len);
#             block_start_byte = int(block_start_byte);
#             block_end_byte = int(block_end_byte);
#             block_end_pos = block_start_pos + block_aln_len;

#             if block_start_pos <= maf_start and block_end_pos > maf_start or get_seq == True:
#                 if debug:
#                     print("BLOCK: ", block_start_pos, block_end_pos, region_start, region_end);
#                 # print(block_start_byte, block_end_byte);

#                 maf_stream.seek(int(block_start_byte));
#                 block = maf_stream.read(int(block_end_byte) - int(block_start_byte)).split("\n");

#                 block = list(filter(None, block))

#                 cur_spec_entry = block[1].split();
                
#                 if debug:
#                     print("ENTRY: ", cur_spec_entry);

#                 block_seq = cur_spec_entry[6];
#                 block_len = int(cur_spec_entry[3]);

#                 if block_len < remaining_region_len:
#                     seq_start_index = maf_start - block_start_pos;
#                     #print(seq_start_index);
#                     seq_to_add = block_seq[seq_start_index:];
#                     remaining_region_len -= len(seq_to_add);
#                     seq += seq_to_add;

#                 elif block_len == remaining_region_len:
#                     seq += block_seq;
#                     remaining_region_len -= block_len;

#                 elif block_len > remaining_region_len:
#                     seq_start_index = maf_start - block_start_pos;
#                     seq_to_add = block_seq[seq_start_index:seq_start_index + remaining_region_len];
#                     remaining_region_len -= len(seq_to_add);
#                     seq += seq_to_add;

#                 if debug:
#                     print("SEQ:", seq);
#                     #print(remaining_region_len);
#                 maf_start = block_end_pos;
#                 get_seq = True;

#                 #print(line);

#                 if remaining_region_len == 0:
#                     break;

#         if test and not seqTest(ref_genome_fasta, seq, chrom, region_start, region_end):
#             sys.exit(1);


