#!/usr/local_rwth/bin/zsh
#SBATCH --account=lect0079
#SBATCH --exclusive
#SBATCH --ntasks=12
# name the job
#SBATCH --job-name=mpi12
# declare the merged STDOUT/STDERR file
#SBATCH --output=output_fine_12CPU.%J.txt
#SBATCH --time=1:00:00

$MPIEXEC $FLAGS_MPI_BATCH ./mpi ./settings.fine.in
