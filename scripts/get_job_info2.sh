#!/bin/bash

# Get the date 14 days ago
startdate=$(date -d "14 days ago" +%Y-%m-%d)

# Get a list of job IDs
jobids=$(sacct --noheader --format=JobID --starttime $startdate | grep -oP '^\d+' | sort | uniq)

# For each job ID
for jobid in $jobids; do
  # Get the job name
  jobname=$(sacct -j $jobid --format=JobName%30 --noheader | head -n 1)

  # Get the elapsed time and MaxRSS
  #elapsed_maxrss_kb=$(sacct -j $jobid --format=JobID%30,Elapsed,MaxRSS --noheader | grep 'batch' | awk '{print $2, $3}')

  # Convert MaxRSS to GB
  #elapsed_maxrss_gb=$(echo $elapsed_maxrss_kb | awk '{printf "%s %.2f\n", $1, $2 / 1048576}')

  job_info=$(sacct -j $jobid --format=JobID%30,Elapsed,MaxRSS,End,State --noheader | grep 'batch' | awk '{printf "%s %.2f %s %s\n", $2, $3 / 1048576, $4, $5}')

  # Print the job ID, job name, elapsed time, and MaxRSS
  echo "$jobid $jobname $job_info"
done