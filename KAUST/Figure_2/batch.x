#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=1
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "Kaust Figure 2"   # job name
#SBATCH --qos=test


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

./plt_nfields.py 0 -n 1 -f
./plt_nvel.py 0 -n v -f
./plt_a.py 78 90 100 -n 1 -f
./plt_b.py 78 90 100 -n v -f
./Vel_mix.py
