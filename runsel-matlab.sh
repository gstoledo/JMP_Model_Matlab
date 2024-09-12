#!/bin/bash

# Define an array of unique job names
job_names=("costs" "costsq06")

# Loop over each job name
for job_name in "${job_names[@]}"
do
    # Submit the job with the current job name
    sbatch --job-name=$job_name runjobs-matlab.sbatch
done
