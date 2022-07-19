#!/usr/local_rwth/bin/zsh
#SBATCH --account=lect0079
#SBATCH --exclusive
#SBATCH --ntasks=1
# name the job
#SBATCH --job-name=mpi1
# declare the merged STDOUT/STDERR file
#SBATCH --output=output_medium_1CPU.%J.txt
#SBATCH --time=1:00:00

$MPIEXEC $FLAGS_MPI_BATCH ./mpi ./settings.medium.in
