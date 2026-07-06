#############################################################################
# An attempt at a simple HMM to predict conserved elements from phyloP 
# conserved sites
#
# Gregg Thomas, April 2024
#############################################################################

import sys
import os
import numpy as np
import logging
import warnings
import traceback

#############################################################################

class OnlineHMM:
# A simple HMM class for online prediction of state sequences

    def __init__(self, n_states, n_symbols, trans_prob, emit_prob, start_prob):
        self.n_states = n_states;
        self.n_symbols = n_symbols;
        self.trans_prob = trans_prob;
        self.emit_prob = emit_prob;
        self.start_prob = start_prob;
    # Initialize the HMM with the number of states, number of symbols, 
    # transition probabilities, emission probabilities, and start probabilities

    ##########

    def predictLog(self, symbol):
    # Predict the most likely state given the current symbol, using log probabilities
        #epsilon = 1e-10;
        # Set a small value to prevent log(0) errors

        # try:
        #     with warnings.catch_warnings():
        #         warnings.simplefilter("error", RuntimeWarning)  
        #         # Turn RuntimeWarning into an error

        if not hasattr(self, 'log_state_prob'):
            self.log_state_prob = np.log(self.start_prob) + np.log(self.emit_prob[:, symbol]);
        # Initialize the log state probabilities for the first symbol

        else:
            #self.log_state_prob = np.logaddexp(np.log(self.trans_prob.T) + self.log_state_prob, axis=1) + np.log(self.emit_prob[:, symbol]);
            # print(self.trans_prob.T, self.log_state_prob, self.emit_prob[:, symbol]);
            # print(1)
            # print(np.log(self.trans_prob.T));
            # print(2)
            # print(np.logaddexp.reduce(np.log(self.trans_prob.T) + self.log_state_prob, axis=1));
            # print(3)
            # print(np.log(self.emit_prob[:, symbol] + epsilon));                    
            self.log_state_prob = np.logaddexp.reduce(np.log(self.trans_prob.T) + self.log_state_prob, axis=1) + np.log(self.emit_prob[:, symbol]);
        # Update the log state probabilities for the next symbol
            
        self.log_state_prob -= np.logaddexp.reduce(self.log_state_prob);
        # Normalize the log state probabilities

        # except RuntimeWarning as e:
        #     print(f"\nCaught the following Warning in predictLog:");
        #     traceback.print_exc();
        #     sys.exit(1)
        # # Catch RuntimeWarnings and exit with an error message

        return self.log_state_prob.argmax();
        # Return the most likely state
    ## End predictLog

    ##########

    def predict(self, symbol):
    # Predict the most likely state given the current symbol
        if not hasattr(self, 'state_prob'):
            self.state_prob = self.start_prob * self.emit_prob[:, symbol];
        # Initialize the state probabilities for the first symbol

        else:
            self.state_prob = (self.trans_prob.T @ self.state_prob) * self.emit_prob[:, symbol];
        # Update the state probabilities for the next symbol          
        
        self.state_prob /= self.state_prob.sum();
        # Normalize the state probabilities
        
        return self.state_prob.argmax();
        # Return the most likely state
    ## End predict

    ##########

    # def update(self, prev_state, curr_state, symbol):
    # # Update the transition and emission probabilities given the previous state, current state, and symbol
        
    #     self.trans_prob[prev_state, curr_state] += 1;
    #     self.trans_prob[prev_state] /= self.trans_prob[prev_state].sum();
    #     # Update transition probabilities
        
    #     self.emit_prob[curr_state, symbol] += 1;
    #     self.emit_prob[curr_state] /= self.emit_prob[curr_state].sum();
    #     # Update emission probabilities
    # ## End update

############################################################################# 

t0_0, t1_1, e0_0, e1_1, s0 = map(float, sys.argv[1:])
params_str = '-'.join(map(str, [t0_0, t1_1, e0_0, e1_1, s0]));
# Get the transition and emission probabilities from the command line

min_prob = 1e-10;
max_prob = 0.9999999999;
# Set the minimum and maximum probabilities

clamp = lambda n, minn, maxn: max(min(maxn, n), minn);
# Define a lambda function to apply the min and max constraints

t0_0, t1_1, e0_0, e1_1, s0 = [ clamp(v, min_prob, max_prob) for v in [t0_0, t1_1, e0_0, e1_1, s0] ];
# Apply the lambda function to each variable using a list comprehension

# t0_0 = max(t0_0, 1e-10)
# t1_1 = max(t1_1, 1e-10)
# # Set a minimum value for the transition probabilities

# e0_0 = max(e0_0, 1e-10);
# e1_1 = max(e1_1, 1e-10);
# # Set a minimum value for the emission probabilities

# s0 = min(s0, 0.9999999999)
# # Set a maximum value for the start probability

####################

n_states = 2;  
# Inside a conserved element (1), outside a conserved element (0)

n_symbols = 2;  
# Conserved site (1), non-conserved site (0)

#t0_0 = 0.9;         # state 0 -> state 0
t0_1 = 1.0 - t0_0;    # state 0 -> state 1

#t1_1 = 0.8;         # state 1 -> state 1
t1_0 = 1.0 - t1_1;    # state 1 -> state 0

trans_prob = np.array([[t0_0, t0_1], 
                       [t1_0, t1_1]]);
# [[stay outside conserved element, transition into conserved element ],
# [transition out of conserved element, stay inside conserved element]]

#e0_0 = 0.8;         # state 0 -> emit 0
e0_1 = 1 - e0_0;    # state 0 -> emit 1

#e1_1 = 0.5;         # state 1 -> emit 1
e1_0 = 1 - e1_1;    # state 1 -> emit 0

emit_prob = np.array([[e0_0, e0_1], 
                      [e1_0, e1_1]]);
# [[observe non-conserved site outside conserved element, observe conserved site outside conserved element],
# [observe non-conserved site inside conserved element, observe conserved site inside conserved element]]

#s0 = 0.9;           # start outside conserved element
s1 = 1 - s0;        # start inside conserved element

start_prob = np.array([s0, s1]);
# [start outside conserved element, start inside conserved element]

hmm = OnlineHMM(n_states, n_symbols, trans_prob, emit_prob, start_prob);
# Initialize the HMM

####################

test = False;
if test:
    test_data = np.random.choice(n_symbols, size=100);
    # Generate some test data

    state_sequence = [hmm.predictLog(symbol) for symbol in test_data];
    # Use the HMM to predict the state sequence for the test data

    print("Test data: ", test_data)
    print("Predicted state sequence: ", state_sequence)
    print(test_data == state_sequence)
    sys.exit();

####################

logging_file = "conserved_elements_hmm.debug";
# with open(logging_file, 'w'):
#     pass;
# Clear the logging file

logging.basicConfig(filename=logging_file, level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s');
# Set the logging level

####################

ref_index_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/00-human-ref-ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai";
index_prefix = "chr";
cons_sites_bed = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/04-phylop-bedgraph/chr19-conserved-sites.bed";
# Input files
## The ref_index_file and the cons_sites_bed must be sorted in the same way w.r.t chromosomes/scaffolds
## The cons_sites_bed must also be sorted by coordinates

chr_list = ["chr19"];
min_elem_len = 20;
max_elem_len = 10000;
# Options

params_str = str(min_elem_len) + "-" + params_str;

summary_file = "conserved-elements-" + params_str + ".log";
cons_elem_bed = "conserved-elements-" + params_str + ".bed";
state_seq_file = "state-seqs-" + params_str + ".txt";
write_state_seqs = False;
# Output files

####################

with open(ref_index_file, 'r') as ref_index_stream:
    ref_index = { index_prefix + line.strip().split('\t')[0]: int(line.strip().split('\t')[1]) for line in ref_index_stream };
# Read the reference index

if chr_list == []:
    chr_list = list(ref_index.keys());
else:
    if not all(chrome in ref_index for chrome in chr_list):
        print(f'Error: Some provided chromosomes in chr_list are not in the reference index', file=sys.stderr);
        sys.exit(1);
# If a list of chromosomes is not provided, set to all chromosomes in the index

####################

if write_state_seqs:
    state_stream = open(state_seq_file, "w");

with open(cons_sites_bed) as cons_stream, open(cons_elem_bed, "w") as out_stream, open(summary_file, "w") as summary_stream:
    log_headers = ['chromosome', 'chromosome.length', 'conserved.sites', 'total.conserved.elems', 'elems.too.short', 'elems.too.long', 'final.cons.elems',
                   'avg.cons.elem.len', 'avg.inter.elem.dist', 'cons.sites.in.cons.elems', 'avg.perc.cons.sites.per.elem',
                   't0_0', 't1_1', 'e0_0', 'e1_1', 's0', 'min.elem.len.cutoff'];
    print("# Command: " + ' '.join(sys.argv), file=summary_stream);
    print("# " + '\t'.join(log_headers), file=summary_stream);

    line = cons_stream.readline();
    cons_site, end_site = map(int, line.strip().split('\t')[1:3]);
    # Read the first line of the conserved sites bed

    prev_state = None
    # Initialize the previous state

    prev_state_end = None;
    # Initialize the previous state end coordinate

    cons_elem_start, cons_elem_end, cons_elem_len, cons_elem_len_sum, inter_elem_dist_sum = None, None, None, 0, 0;
    cons_sites_in_cons_elems, perc_cons_sites_per_elem_sum = 0, 0;
    # Initialize the conserved element coordinates

    symbol, state = None, None;
    # Initialize the symbol and state

    for chrome in chr_list:
        if chrome != line.strip().split('\t')[0]:
            print(f'Error: Mismatch between ref_index and cons_sites_bed at {chrome}', file=sys.stderr);
            sys.exit(1);
        # At the beginning of each chrome, make sure the next conserved site is on the same chrome
        # These could mis-match if they aren't sorted correctly

        num_cons_elems = 0;
        cons_sites_in_cur_elem = 0;
        total_cons_sites_in_chr = 0;
        num_elems_too_short = 0;
        num_elems_too_long = 0;
        # Initialize counters

        for i in range(ref_index[chrome]):
        # Loop over every site in the chromosome, using 0-based coordinates (as in the input bed file)

            #logging.debug(f"a {i} {symbol} {state} {prev_state} {cons_site} {end_site}");
            # DEBUG

            cons_str = "OUT";
            trans_str = "";

            if i == cons_site:
            # If the current site is the next conserved site

                #print("SITE", line);

                symbol = 1;
                # Update the current symbol to 1, conserved site

                cons_site += 1;
                # Increment the conserved site position

                #num_cons_sites += 1;
                total_cons_sites_in_chr += 1;
                # Increment the counters

                if cons_site >= end_site:
                # If the next site is past the end coordinate of this chunk of conserved sites

                    line = cons_stream.readline();
                    if line:
                        cons_site, end_site = map(int, line.strip().split('\t')[1:3]);
                    else:
                        cons_site = None;
                        end_site = None;
                    # Read the next conserved site chunk, or set to None if there are no more
            else:
                symbol = 0;

            # If the current site is not a conserved site, set the symbol to 0

            state = hmm.predictLog(symbol)
            # Predict the current state (in conserved element (1) or not in conserved element (0)) 
            # given the current symbol (conserved site (1) or non-conserved site (0)

            if state == 1:
                cons_str = "IN";
                if symbol == 1:
                    cons_sites_in_cur_elem += 1;
            # Increment the conserved site counter if the current state is inside a conserved element
            
            if prev_state is not None:
            # If this isn't the first conserved element state

                if prev_state == 0 and state == 1:

                    #logging.debug("START ELEMENT");
                    # DEBUG

                    trans_str = "START";

                    cons_elem_start = i;
                    num_cons_elems += 1;
                # Check if the previous state was outside a conserved element and the current state is inside a conserved element
                # If so, set the start of the conserved element

                elif prev_state == 1 and state == 0:
                # Check if the previous state was inside a conserved element and the current state is outside a conserved element

                    #logging.debug("END ELEMENT");
                    # DEBUG
                    
                    trans_str = "END";
                    elem_str = "";
                    
                    cons_elem_end = i;
                    cons_elem_len = cons_elem_end - cons_elem_start;
                    # Set the end of the conserved element and calculate the length

                    if cons_elem_len <= min_elem_len:
                        elem_str = "S";
                        num_elems_too_short += 1;
                    elif cons_elem_len >= max_elem_len:
                        elem_str = "L";
                        num_elems_too_long += 1;
                    else:
                        cons_elem_len_sum += cons_elem_len;

                        if prev_state_end is not None:
                            inter_elem_dist_sum += cons_elem_start - prev_state_end;
                        # Calculate the distance between this conserved element and the previous one

                        prev_state_end = cons_elem_end;
                        # Update the previous state end coordinate

                        perc_cons_sites_in_elem = cons_sites_in_cur_elem / cons_elem_len;
                        perc_cons_sites_per_elem_sum += perc_cons_sites_in_elem;
                        # Calculate the percentage of conserved sites in the conserved element

                        cons_sites_in_cons_elems += cons_sites_in_cur_elem;
                        # Increment the conserved site counter for the conserved element

                        outline = f'{chrome}\t{cons_elem_start}\t{cons_elem_end}\t{cons_elem_len}\t{cons_sites_in_cur_elem}\t{perc_cons_sites_in_elem}';
                        print(outline, file=out_stream);
                        # Write the conserved element to the output file

                        #logging.debug(outline);
                        # DEBUG

                        if cons_sites_in_cur_elem > cons_elem_len:
                            print("ERROR: More conserved sites than conserved element length. See debug file.", file=sys.stderr);
                            logging.error('Command: ' + ' '.join(sys.argv));
                            logging.error(f"b {i} {symbol} {state} {prev_state} {cons_site} {end_site}")
                            sys.exit(1);
                    # If the conserved element is within the length limits, write it to the output file

                    cons_sites_in_cur_elem = 0;
                    # Reset the conserved site counter for the next element
            ## End state transition block

            prev_state = state;
            # Update the previous state for the next iteration

            if write_state_seqs:
                stateline = f'{chrome}\t{i}\t{symbol}\t{state}\t{cons_str}\t{trans_str}';
                if trans_str == "END":
                    stateline += f'\t{cons_elem_start}\t{cons_elem_end}\t{cons_elem_len}';
                    if elem_str not in ["S", "L"]:
                        stateline += f'\t{perc_cons_sites_in_elem}';

                print(stateline, file=state_stream);
            # Write the state sequence to a file if that option is set

            #logging.debug(f"b {i} {symbol} {state} {prev_state} {cons_site} {end_site}");
            # DEBUG
        ## End site loop

        final_cons_elems = num_cons_elems - num_elems_too_short - num_elems_too_long;
        avg_elem_len = cons_elem_len_sum / final_cons_elems if final_cons_elems > 0 else 0;
        avg_inter_elem_dist = inter_elem_dist_sum / final_cons_elems if final_cons_elems > 0 else 0;
        avg_perc_cons_sites_per_elem = perc_cons_sites_per_elem_sum / final_cons_elems if final_cons_elems > 0 else 0;
        # Calculate the summary statistics for the chromosome

        print(
            f'{chrome}\t{ref_index[chrome]}\t{total_cons_sites_in_chr}\t{num_cons_elems}\t{num_elems_too_short}\t{num_elems_too_long}\t{final_cons_elems}\t'
            f'{avg_elem_len}\t{avg_inter_elem_dist}\t{cons_sites_in_cons_elems}\t{avg_perc_cons_sites_per_elem}\t'
            f'{t0_0}\t{t1_1}\t{e0_0}\t{e1_1}\t{s0}\t{min_elem_len}', 
            file=summary_stream
        )        
        # Write the summary statistics for the chromosome to the log file
          
    ## End chrome loop 
## Close conserved files

if write_state_seqs:
    state_stream.close();
# Close the state sequence file if it was opened

####################
