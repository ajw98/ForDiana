#!/bin/bash
#SBATCH --job-name=plot_times
#SBATCH --output=slurm-%A.out
#SBATCH --time=1:00:00            # Time limit
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks=1                # Number of tasks (processes)
#SBATCH --mem=10G                 # Total memory
#SBATCH --gpus=1                   # Number of GPUs
#SBATCH --constraint=pascal        # GPU constraint

module load python

for t in $(seq 8 17); do
    echo "Running t = $t"
    ./plt_nfields.py $t
done
