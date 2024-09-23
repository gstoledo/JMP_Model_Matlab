#!/bin/bash

# Define an array of unique job names
job_names=("costsq05" "costsq04"  "costsq06" "costsq055" "alpha_A" "alpha_rho")

# Loop over each job name
for job_name in "${job_names[@]}"
do
    # Submit the job with the current job name
    sbatch --job-name=$job_name runjobs-matlab.sbatch
done
