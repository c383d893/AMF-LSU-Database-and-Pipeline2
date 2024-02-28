#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4gb
#SBATCH --output=./slurmOutputs/Alignseqs
#SBATCH --error=./slurmOutputs/Alignseqs
#SBATCH --job-name=Alignseqs

module load gcc/8.2.0 r/4.2.2

# Get the working directory 
SCRIPT_DIR=$1

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat AMFcutLSUdb.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > AMFcutLSUdbwithCutoffs.R

echo;echo "Beginning cutdb pipeline..."
Rscript AMFcutLSUdbwithCutoffs.R
echo;echo “cutdb pipeline complete”
rm AMFcutLSUdbwithCutoffs.R 

# join cut LSU DB and study seqs

cat v16_LSUDB_2024_cut.fasta otu97repseqs_clean_BLAST.fasta > BLAST_otu97plusv16_2024_cut.fasta

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat AMFalignseqs.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > AMFalignseqswithCutoffs.R

echo;echo "Beginning split R1-R2 pipeline..."
Rscript AMFalignseqswithCutoffs.R
echo;echo “split R1-R2 pipeline complete”
rm AMFalignseqswithCutoffs.R

# Align R1 and R2 independently: import .fasta, align via mafft, export back to .fasta

### Activate conda 
shift
. ~/.bashrc
conda activate /cluster/project/crowther/miniconda3/envs/microbio

# Define a temporary folder 
# mkdir $SCRIPT_DIR/tmp/
export TMPDIR=$SCRIPT_DIR/tmp/

# R1
qiime tools import --input-path R1.BLAST_otu97plusv16_2024_cut.fasta --output-path R1.BLAST_otu97plusv16_2024_cut --type 'FeatureData[Sequence]'
qiime alignment mafft --i-sequences R1.BLAST_otu97plusv16_2024_cut.qza --o-alignment aligned_R1.BLAST_otu97plusv16_2024_cut.qza
qiime tools export --input-path aligned_R1.BLAST_otu97plusv16_2024_cut.qza --output-path R1.BLAST_otu97plusv16_2024_cut_out

# R2
qiime tools import --input-path R2.BLAST_otu97plusv16_2024_cut.fasta --output-path R2.BLAST_otu97plusv16_2024_cut --type 'FeatureData[Sequence]'
qiime alignment mafft --i-sequences R2.BLAST_otu97plusv16_2024_cut.qza --o-alignment aligned_R2.BLAST_otu97plusv16_2024_cut.qza
qiime tools export --input-path aligned_R2.BLAST_otu97plusv16_2024_cut.qza --output-path R2.BLAST_otu97plusv16_2024_cut_out
echo;echo “Split R1-R2 alignment complete”

### Deactivate conda 
conda deactivate

# Load the modules again 
module load gcc/8.2.0 r/4.2.2

# Concatenate alignments (R1 and R2)
echo;echo "Beginning concatenate R1-R2 pipeline..."
Rscript AMFconcatseqs.R
echo;echo “Concatenate R1-R2 pipeline complete”

# remove once no longer needed:
rm R1.BLAST_otu97plusv16_2024_cut.fasta
rm R1.BLAST_otu97plusv16_2024_cut_out
rm R1.BLAST_otu97plusv16_2024_cut.qza
rm R2.BLAST_otu97plusv16_2024_cut.fasta
rm R2.BLAST_otu97plusv16_2024_cut_out
rm R2.BLAST_otu97plusv16_2024_cut.qza
rm aligned_R1.BLAST_otu97plusv16_2024_cut.qza
rm aligned_R2.BLAST_otu97plusv16_2024_cut.qza
rm R1.BLAST_otu97plusv16_2024_cut_out
rm R2.BLAST_otu97plusv16_2024_cut_out

rm -r tmp/
