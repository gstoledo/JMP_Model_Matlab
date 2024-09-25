#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB
#SBATCH --time=20:00:00
#SBATCH --output=/home/gst247/HPC_Model_Matlab/Calibration_outcomes/%x/%x_%j.out   # Redirect standard output to Calibration_outcomes/filename/MATLAB_jobid.out
#SBATCH --error=/home/gst247/HPC_Model_Matlab/Calibration_outcomes/%x/%x_%j.err    # Redirect standard error to Calibration_outcomes/filename/MATLAB_jobid.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gst247@nyu.edu

module purge
module load matlab/2023b

# Create temporary directories in your home directory
export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)
export MATLAB_LOG_DIR=$SLURM_JOBTMP

cp -rp /share/apps/matlab-slurm/20201018/slurm/shared/parallel.mlsettings $MATLAB_PREFDIR

srun matlab -nodisplay -r "calibration_jobs('$SLURM_JOB_NAME')"

exit