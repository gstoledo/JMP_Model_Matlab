#!/bin/bash

#SBATCH --job-name=calibrationv0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60GB
#SBATCH --time=05:00:00
#SBATCH --output=Matlab-%j.out
#SBATCH --error=Matlab-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gst247@nyu.edu

module purge
module load matlab/2023b

# Create temporary directories in your home directory
mkdir -p /home/gst247/tmp
export MATLAB_PREFDIR=$(mktemp -d /home/gst247/tmp/matlab-XXXX)
export MATLAB_LOG_DIR=/home/gst247/tmp



cat<<EOF | matlab -nodisplay

% Add your working directory to MATLAB's path
addpath('/home/gst247/HPC_Model_Matlab/JMP_Model_Matlab');

c = parcluster;
c.AdditionalProperties.WallTime = '05:00:00';
pool = c.parpool(20);

% Run your MATLAB script
calibrationv0

exit    

EOF
