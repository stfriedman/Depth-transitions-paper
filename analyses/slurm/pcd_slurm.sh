#!/bin/bash 
#SBATCH --partition=general
#SBATCH --job-name=musscrat
#SBATCH --mem-per-cpu=2000
#SBATCH --time=30-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=sarahtfried@gmail.com
#SBATCH --output=output/%x_output_%a_%j.txt
#SBATCH --error=error/%x_err_%a_%j.txt

echo "ARRAY ID: " $SLURM_ARRAY_TASK_ID
rb musscrat_UCLN_state_dependent_mrm.Rev $SLURM_ARRAY_TASK_ID 