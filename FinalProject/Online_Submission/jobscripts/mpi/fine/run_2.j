#!/usr/local_rwth/bin/zsh
#SBATCH --account=lect0079
#SBATCH --exclusive
#SBATCH --ntasks=2
# name the job
#SBATCH --job-name=mpi2
# declare the merged STDOUT/STDERR file
#SBATCH --output=output_fine_2CPU.%J.txt
#SBATCH --time=1:00:00

$MPIEXEC $FLAGS_MPI_BATCH ./mpi ./settings.fine.in
