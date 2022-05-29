#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0079
#SBATCH --job-name=Run_2d_unsteady
# Outputs of the job
#SBATCH --output=out_finest.%j
#SBATCH --error=err_finest.%j
# Wall clock limit
#SBATCH --time=01:00:00

# run the process
./2d_Unsteady settings.finest.in
