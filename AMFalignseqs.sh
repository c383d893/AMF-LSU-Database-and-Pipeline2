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
module load R

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat AMFcutLSUdb.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > AMFcutLSUdbwithCutoffs.R

echo;echo "Beginning cutdb pipeline..."
Rscript AMFcutLSUdbwithCutoffs.R
echo;echo “cutdb pipeline complete”
rm AMFcutLSUdbwithCutoffs.R 

# join cut LSU DB and study seqs

cat V15_LSUDB_3.23.21_cut.fasta otu97repseqs_clean_BLAST.fasta > BLAST_otu97plusV15_3.23.21_cut.fasta

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat AMFalignseqs.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > AMFalignseqswithCutoffs.R

echo;echo "Beginning split R1-R2 pipeline..."
Rscript AMFalignseqswithCutoffs.R
echo;echo “split R1-R2 pipeline complete”
rm AMFalignseqswithCutoffs.R

# Align R1 and R2 independently: import .fasta, align via mafft, export back to .fasta

# R1
module load qiime2
qiime tools import --input-path R1.BLAST_otu97plusV15_3.23.21_cut.fasta --output-path R1.BLAST_otu97plusV15_3.23.21_cut --type 'FeatureData[Sequence]'
qiime alignment mafft --i-sequences R1.BLAST_otu97plusV15_3.23.21_cut.qza --o-alignment aligned_R1.BLAST_otu97plusV15_3.23.21_cut.qza
qiime tools export --input-path aligned_R1.BLAST_otu97plusV15_3.23.21_cut.qza --output-path R1.BLAST_otu97plusV15_3.23.21_cut_out

# R2
qiime tools import --input-path R2.BLAST_otu97plusV15_3.23.21_cut.fasta --output-path R2.BLAST_otu97plusV15_3.23.21_cut --type 'FeatureData[Sequence]'
qiime alignment mafft --i-sequences R2.BLAST_otu97plusV15_3.23.21_cut.qza --o-alignment aligned_R2.BLAST_otu97plusV15_3.23.21_cut.qza
qiime tools export --input-path aligned_R2.BLAST_otu97plusV15_3.23.21_cut.qza --output-path R2.BLAST_otu97plusV15_3.23.21_cut_out
echo;echo “Split R1-R2 alignment complete”

# Concatenate alignments (R1 and R2)
echo;echo "Beginning concatenate R1-R2 pipeline..."
Rscript AMFconcatseqs.R
echo;echo “Concatenate R1-R2 pipeline complete”


