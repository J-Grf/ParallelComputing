#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0079
# Outputs of the job
#SBATCH --output=out.%j
#SBATCH --error=err.%j
# Wall clock limit
#SBATCH --time=0:01:00

# run the process
./2d_Unsteady settings.coarse.in
