#!/bin/bash
#SBATCH --job-name={rule}.{wildcards}
#SBATCH --output=logs/{rule}.{wildcards}.out
#SBATCH --error=logs/{rule}.{wildcards}.err
#SBATCH --mem={resources.mem}
#SBATCH --cpus-per-task={resources.cpus}
#SBATCH --time={resources.time}
#SBATCH --partition={resources.partition}

echo "SLURM job ID: $SLURM_JOB_ID"
echo "Job running on $(hostname)"
echo "Using $(nproc) CPU cores"
echo "Started at $(date)"
echo "Node: $(hostname)"
echo "Working in: $(pwd)"
{exec_job}
echo "Ended at $(date)"
