#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --gpus=1
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "Kaust Figure 3"   # job name
#SBATCH --qos=test


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

cd 11 

./plt_fields.py 23 -n 1 -f
./plt_a.py 50 -n 1 -f
./plt_vel.py 23 -n v -f
./plt_b.py 50 -n v -f

cd ../12

./plt_fields.py 7 -n 1 -f
./plt_a.py 48 -n 1 -f
./plt_vel.py 7 -n v -f
./plt_b.py 48 -n v -f

cd ..

./Vel_mix.py
