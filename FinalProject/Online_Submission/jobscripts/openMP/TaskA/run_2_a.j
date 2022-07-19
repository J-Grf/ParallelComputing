#!/usr/bin/env bash
# Lecture queue
#SBATCH --account=lect0079

### Ask for exactly one node -> all allocated cpus must be on this one
#SBATCH --nodes=1

### Ask for <1 or 2 or 4 or 6 or 8 or 12> cpus.
#SBATCH --cpus-per-task=2

### Divide the needed memory per task through the cpus-per-task,
### as slurm requests memory per cpu, not per task!
### Example:
### You need 2 GB memory per task, you have 8 cpus per task ordered
### order 2048/ <1 or 2 or 4 or 6 or 8 or 12> ->
### <2048M or 1024M or 512M or 342M or 256M or 171M> memory per task.
### M is the default and can therefore be omitted,
### but could also be K(ilo)|G(iga)|T(era).
#SBATCH --mem-per-cpu=1024M

### Name the job.
#SBATCH --job-name=2d_Unsteady_OpenMP_A

# Outputs of the job.
#SBATCH --output=out_2CPU_A.%J.txt
#SBATCH --error=err_2CPU_A.%J.txt

# Wall clock limit.
#SBATCH --time=1:00:00

#SBATCH --exclusive

# run the process
# export OMP_SCHEDULE="static,128"
./openmpA ./settings.coarse.in
