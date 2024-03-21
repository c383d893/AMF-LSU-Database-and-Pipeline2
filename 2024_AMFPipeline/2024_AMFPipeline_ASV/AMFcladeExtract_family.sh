#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --time=00-06:00:00
#SBATCH --output=./slurmOutputs/AMFcladeExtract_family.out
#SBATCH --error=./slurmOutputs/AMFcladeExtract_family.out
#SBATCH --job-name=AMFcladeExtract_family.out
#SBATCH --partition=sixhour

### Activate conda
. ~/.bashrc
conda activate $C_ENV

Rscript AMFcladeExtract_family.R
