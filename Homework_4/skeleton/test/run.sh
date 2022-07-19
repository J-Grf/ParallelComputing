#!/usr/local_rwth/bin/zsh
#SBATCH --account=lect0079
## #SBATCH --exclusive
#SBATCH --ntasks=<1,2,4,6,8,12,16>
# name the job
#SBATCH --job-name=homework4
# declare the merged STDOUT/STDERR file
#SBATCH --output=output.coarse-bis.cpu4.%J.txt
#SBATCH --time=1:00:00

$MPIEXEC $FLAGS_MPI_BATCH ../build/2d_Unsteady ./settings.<fine,finest,finest-bis>.in
