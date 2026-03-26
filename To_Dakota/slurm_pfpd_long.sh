#!/bin/bash
#SBATCH --job-name=extended_time          # Optional: name of your job
#SBATCH --time=72:00:00            # Time limit
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks=1                # Number of tasks (processes)
#SBATCH --mem=42G                 # Total memory
#SBATCH --gpus=1                   # Number of GPUs
#SBATCH --constraint=pascal        # GPU constraint

# Run your executable
srun ./rungpu.x
