#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --partition=sixhour
#SBATCH --output=./slurmOutputs/Alignseqs
#SBATCH --error=./slurmOutputs/Alignseqs
#SBATCH --job-name=Alignseqs

module load qiime2
module load R/3.6

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat CutLSUdb.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > CutLSUdbwithCutoffs.R

echo;echo "Beginning cutdb pipeline..."
Rscript CutLSUdbwithCutoffs.R
echo;echo “cutdb pipeline complete”

rm CutLSUdbwithCutoffs.R 

# join cut LSU DB and study seqs

cat V15_LSUDB_3.23.21_cut.fasta otu97repseqs_clean.fasta > otu97plusV15_3.23.21_cut.fasta

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat Alignseqs.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > AlignseqswithCutoffs.R

echo;echo "Beginning split R1-R2 pipeline..."
Rscript AlignseqswithCutoffs.R
echo;echo “split R1-R2 pipeline complete”

# Align R1 and R2 independently

# R1
module load qiime2
qiime tools import \
  --input-path R1.otu97plusV15_3.23.21_cut.fasta \
  --output-path R1.otu97plusV15_3.23.21_cut \
  --type 'FeatureData[Sequence]'
qiime alignment mafft \
  --i-sequences R1.otu97plusV15_3.23.21_cut.qza \
  --o-alignment aligned_R1.otu97plusV15_3.23.21_cut.qza
qiime tools export \
  --input-path aligned_R1.otu97plusV15_3.23.21_cut.qza  \
  --output-path R1.otu97plusV15_3.23.21_cut_out

qiime tools import \
  --input-path R2.otu97plusV15_3.23.21_cut.fasta \
  --output-path R2.otu97plusV15_3.23.21_cut \
  --type 'FeatureData[Sequence]'
qiime alignment mafft \
  --i-sequences R2.otu97plusV15_3.23.21_cut.qza \
  --o-alignment aligned_R2.otu97plusV15_3.23.21_cut.qza
qiime tools export \
  --input-path aligned_R2.otu97plusV15_3.23.21_cut.qza  \
  --output-path R2.otu97plusV15_3.23.21_cut_out

# Concatenate alignments (R1 and R2)
echo;echo "Beginning concatenate R1-R2 pipeline..."
Rscript Concatseqs.R
echo;echo “concatenate R1-R2 pipeline complete”


