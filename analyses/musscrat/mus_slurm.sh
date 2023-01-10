#!/bin/bash 
#SBATCH --partition=general
#SBATCH --job-name=musscrat
#SBATCH --mem-per-cpu=2000
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=30-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=sarahtfried@gmail.com
#SBATCH --output=output/%x_output_%a_%j.txt
#SBATCH --error=error/%x_err_%a_%j.txt

module load revbayes/1.0.11-foss-2018b
rb musscrat_UCLN_state_dependent_mrm.Rev