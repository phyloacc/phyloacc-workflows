#!/bin/bash

#SBATCH --job-name=hmm_run
#SBATCH --array=0-9999
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=hmm_run_%A_%a.out
#SBATCH --error=hmm_run_%A_%a.err

# Define the lists of values
values=($(seq 0 0.1 1.0))

# Calculate the indices for each value
index1=$((SLURM_ARRAY_TASK_ID / 1000))
index2=$((SLURM_ARRAY_TASK_ID % 1000 / 100))
index3=$((SLURM_ARRAY_TASK_ID % 100 / 10))
index4=$((SLURM_ARRAY_TASK_ID % 10))

# Get the values
val1=${values[$index1]}
val2=${values[$index2]}
val3=${values[$index3]}
val4=${values[$index4]}

# Run the Python script
python ../conserved_elements_hmm.py $val1 $val2 $val3 $val4 0.9

# # Define the lists of values
# t0_0=(0.1 0.25 0.5 0.75 0.9)
# t1_1=(0.1 0.25 0.5 0.75 0.9)
# e0_0=(0.1 0.25 0.5 0.75 0.9)
# e1_1=(0.1 0.25 0.5 0.75 0.9)
# s0=(0.1 0.25 0.5 0.75 0.9)

# # Nested loops to iterate over the lists
# for val1 in "${t0_0[@]}"; do
#     for val2 in "${t1_1[@]}"; do
#         for val3 in "${e0_0[@]}"; do
#             for val4 in "${e1_1[@]}"; do
#                 echo "Running with values: $val1 $val2 $val3 $val4"

#                 python ../conserved_elements_hmm.py $val1 $val2 $val3 $val4 0.9


#             done
#         done
#     done
# done

# Define the lists of values
# t0_0=(0.1 0.25 0.5 0.75 0.9)
# t1_1=(0.1 0.25 0.5 0.75 0.9)
# e0_0=(0.1 0.25 0.5 0.75 0.9)
# e1_1=(0.1 0.25 0.5 0.75 0.9)
# s0=(0.1 0.25 0.5 0.75 0.9)

# # Export the Python script path as a function for parallel to use
# do_work() {
#     val1=$1
#     val2=$2
#     val3=$3
#     val4=$4
#     echo "Running with values: $val1 $val2 $val3 $val4"
#     python ../conserved_elements_hmm.py $val1 $val2 $val3 $val4 0.9
# }
# export -f do_work

# # Use parallel to run the function
# parallel -j 10 do_work ::: "${t0_0[@]}" ::: "${t1_1[@]}" ::: "${e0_0[@]}" ::: "${e1_1[@]}"