#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --time=00-06:00:00
#SBATCH --output=./slurmOutputs/AMFcladeExtract.out
#SBATCH --error=./slurmOutputs/AMFcladeExtract.out
#SBATCH --job-name=AMFcladesExtract.out

### Activate conda
. ~/.bashrc
conda activate $C_ENV

Rscript AMFcladesExtract.R
