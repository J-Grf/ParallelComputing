#!/usr/local_rwth/bin/zsh

#SBATCH --nodes=3
#SBATCH --ntasks=144
#SBATCH --ntasks-per-node=48

#################
# ATTENTION !!! #
#################
# Divide the needed memory per task through the cpus-per-task, as slurm requests memory per cpu, not per task !
# Example:
# You need 24 GB memory per task, you have 24 cpus per task ordered
# order 24GB/24 -> 1G memory per cpu (i.e., per thread)
#SEBATCH --mem-per-cpu=1G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)
 
# name the job
#SBATCH --job-name=Run
 
# declare the merged STDOUT/STDERR file
#SBATCH --output=output.txt
#SEBATCH --error=err_init.txt

#wall Clock limit
#SBATCH --time=20:00:00 
### beginning of executable commands
# Note: the OMP_NUM_THREADS envvar is set automatically - do not overwrite!

echo "Unloading modules:"
module unload intel
module unload intelmpi
echo "Loading modules:"
module load gcc/9
module load openmpi


$MPIEXEC $FLAGS_MPI_BATCH ./maia properties_run.toml 
