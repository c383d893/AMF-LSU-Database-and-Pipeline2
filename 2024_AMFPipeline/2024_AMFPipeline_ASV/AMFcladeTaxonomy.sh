#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --time=00-06:00:00
#SBATCH --output=./slurmOutputs/AMFcladeTaxonomy.out
#SBATCH --error=./slurmOutputs/AMFcladeTaxonomy.out
#SBATCH --job-name=AMFcladeTaxonomy.out

### Activate conda
. ~/.bashrc
conda activate $C_ENV

Rscript AMFcladeTaxonomy.R
