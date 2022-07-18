#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0079
#SBATCH --job-name=Run_2d_unsteady_serial
# Outputs of the job
#SBATCH --output=out_medium.%j
#SBATCH --error=err_medium.%j
# Wall clock limit
#SBATCH --time=01:00:00

#SBATCH --exclusive
# run the process
./serial ./settings.medium.in
