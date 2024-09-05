#!/bin/bash

#SBATCH --job-name=calibrationv0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB
#SBATCH --time=07:00:00
#SBATCH --output=Matlab-%j.out
#SBATCH --error=Matlab-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gst247@nyu.edu

module purge
module load matlab/2023b

# Create temporary directories in your home directory
export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)
export MATLAB_LOG_DIR=$SLURM_JOBTMP

cp -rp /share/apps/matlab-slurm/20201018/slurm/shared/parallel.mlsettings $MATLAB_PREFDIR

cat<<EOF | matlab -nodisplay

% Add your working directory to MATLAB's path
addpath('/home/gst247/HPC_Model_Matlab/JMP_Model_Matlab');



% Run your MATLAB script
calibrationv0

exit

EOF
