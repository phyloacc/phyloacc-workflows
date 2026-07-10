#############################################################################
# Given a file with NCBI accessions, this script uses the NCBI datasets tool
# to look up the GC content of each assembly and add that column to the file.
#
# Gregg Thomas, October 2023
#############################################################################

import sys
import os
import datetime
import subprocess
import json

#############################################################################

now = datetime.datetime.now()
date_time_string = now.strftime("%Y-%m-%d %H:%M:%S")

# Get the input file names and optional accession header name from the command line
if len(sys.argv) < 4 or len(sys.argv) > 5:
    print("Usage: {} sample_file output_file avg_gc_file [accession_header]".format(sys.argv[0]))
    sys.exit(1)
sample_file, output_file, avg_gc_file = sys.argv[1:4]
accession_header = sys.argv[4] if len(sys.argv) == 5 else False
# Get the input file name from the command line

#file_extension = os.path.splitext(sample_file)[1];
#output_file = sample_file.replace(file_extension, "-gc" + file_extension);
# Add "-gc" to the file name before the extension

#avg_gc_file = sample_file.replace(file_extension, "-avg-gc" + file_extension);
# A file to write just the average GC percent to

with open(output_file, "w") as out_stream:
    first = True;
    num_retrievals_failed, num_asm_no_gc = 0, 0;
    num_gc_retrieved, gc_sum = 0.0, 0.0;
    for line in open(sample_file):
        if line.startswith("#"):
            out_stream.write(line);
            continue;
        # Skip comment lines in the input file

        line_list = line.strip().split(",");
        # Parse the line

        if first:
            out_stream.write("# Updated with get_gc_content.py, " + date_time_string + "\n");
            # Write a new comment to the file with some upated info

            headers = [ col.lower() for col in line_list ];
            # The first line without a comment is the header line

            if accession_header:
                if accession_header.lower() not in headers:
                    print("Error: Specified accession header \"" + accession_header + "\" not found in headers");
                    sys.exit(1);
                else:
                    accession_index = headers.index(accession_header.lower());

            else:
                accession_columns = [i for i, col in enumerate(headers) if "accession" in col.lower()]
                if len(accession_columns) == 0:
                    print("Error: No column with 'accession' found in headers")
                    sys.exit(1)
                elif len(accession_columns) > 1:
                    print("Error: Multiple columns with 'accession' found in headers")
                    sys.exit(1)
                else:
                    accession_index = accession_columns[0]

            # if "accession" not in headers:
            #     print("Error: \"accession\" not found in headers");
            #     sys.exit(1);
            # else:
            #     accession_index = headers.index("accession");
            # Look for the accession column in the headers and error out if it doesn't exist

            has_gc_col = False;
            gc_index = False;
            if "gc" in headers:
                print("# input file has gc column");
                gc_index = headers.index("gc");
                has_gc_col = True;
            else:
                print("# input file does not have gc column...adding it");
                gc_index = len(line_list);
                line_list.append("gc");
                has_gc_col = False;
            # Check to see if a gc header exists, and if not add one

            out_stream.write(",".join(line_list) + "\n");
            # Write the headers back to the output file

            first = False;
            continue;

        ##########

        accession = line_list[accession_index];
        print("accession " + accession);
        # Get the accession from the line

        get_gc = True;
        if has_gc_col:
            gc = line_list[gc_index];
            

            try:
                gc = float(gc);
                print("gc from input file: " + str(gc));
                get_gc = False;

                num_gc_retrieved += 1.0;
                gc_sum += gc;

            except ValueError:
                print("gc from input file is not a number...retrieving it (" + gc + ")");
                get_gc = True;
        else:
            line_list.append("NA");
        # Try and look up the gc percent from the line, and if it's not a number then run the datasets command below

        ##########
        
        if get_gc:
            cmd = ["datasets", "summary", "genome", "accession", accession]
            print("Running command:", " ".join(cmd))
            # Build the datasets command

            result = subprocess.run(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
            # Run the command and capture its output

            if result.returncode == 0:
                output = result.stdout.decode('utf-8');
                # Store the output in a string

                asm_data = json.loads(output);
                # Parse the JSON string

                try:
                    gc = asm_data['reports'][0]['assembly_stats']['gc_percent'];
                    # Get the GC percent from the JSON object

                    try:
                        gc = float(gc);
                        print("gc from datasets command: " + str(gc));
                        num_gc_retrieved += 1.0;
                        gc_sum += gc;
                    except ValueError:
                        print("gc from datasets command is not a number...setting to NA (" + gc + ")");
                        gc = "NA";
                        num_asm_no_gc += 1;
                    # Make sure the GC percent is a number

                except KeyError:
                    print("gc key not found...setting to NA")
                    gc = "NA";
                    num_asm_no_gc += 1;
                # Get the GC percent from the JSON object
                # TODO: check if there are multiple reports? print(len(asm_data['reports']));

            else:
                error_message = result.stderr.decode('utf-8');
                print("error during datasets command: " + error_message);
                gc = "NA";
                # Print the error message
                # TODO: Re-run if the error is "Too many requests" or "Connection reset by peer", or maybe pause between running

                num_retrievals_failed += 1;
            # Check if the command completed without errors

        ##########

        line_list[gc_index] = str(gc);
        # Put the gc percent back into the line

        out_stream.write(",".join(line_list) + "\n");
        # Write the line back to the output file

print("num accessions failed lookup:", num_retrievals_failed);
print("num accessions with no gc value:", num_asm_no_gc);
# TODO: Check if too many failed and exit with an error

avg_gc = (gc_sum / num_gc_retrieved) / 100;
print("average gc percent for all accessions with gc value: ", avg_gc);
with open(avg_gc_file, "w") as out_stream:
    out_stream.write(str(avg_gc));



#############################################################################