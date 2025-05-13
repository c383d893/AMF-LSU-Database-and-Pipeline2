#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4gb
#SBATCH --output=./slurmOutputs/Alignseqs
#SBATCH --error=./slurmOutputs/Alignseqs
#SBATCH --job-name=Alignseqs

# Get the working directory 
SCRIPT_DIR=$1

### Activate conda 
. ~/.bashrc
conda activate $C_ENV

# Define a temporary folder
mkdir $SCRIPT_DIR/tmp/
export TMPDIR=$SCRIPT_DIR/tmp/

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat AMFcutLSUdb.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > AMFcutLSUdbwithCutoffs.R

echo;echo "Beginning cutdb pipeline..."
Rscript AMFcutLSUdbwithCutoffs.R
echo;echo “cutdb pipeline complete”
rm AMFcutLSUdbwithCutoffs.R 

# join cut LSU DB and study seqs

cat V18_LSUDB_052025_cut.fasta otu97repseqs_clean_BLAST.fasta > BLAST_otu97plusV18_052025_cut.fasta

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2
cat AMFalignseqs.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > AMFalignseqswithCutoffs.R

echo;echo "Beginning split R1-R2 pipeline..."
Rscript AMFalignseqswithCutoffs.R
echo;echo “split R1-R2 pipeline complete”
rm AMFalignseqswithCutoffs.R

# Align R1 and R2 independently: import .fasta, align via mafft, export back to .fasta

# R1
qiime tools import --input-path R1.BLAST_otu97plusV18_052025_cut.fasta --output-path R1.BLAST_otu97plusV18_052025_cut --type 'FeatureData[Sequence]'
qiime alignment mafft --i-sequences R1.BLAST_otu97plusV18_052025_cut.qza --o-alignment aligned_R1.BLAST_otu97plusV18_052025_cut.qza
qiime tools export --input-path aligned_R1.BLAST_otu97plusV18_052025_cut.qza --output-path R1.BLAST_otu97plusV18_052025_cut_out

# R2
qiime tools import --input-path R2.BLAST_otu97plusV18_052025_cut.fasta --output-path R2.BLAST_otu97plusV18_052025_cut --type 'FeatureData[Sequence]'
qiime alignment mafft --i-sequences R2.BLAST_otu97plusV18_052025_cut.qza --o-alignment aligned_R2.BLAST_otu97plusV18_052025_cut.qza
qiime tools export --input-path aligned_R2.BLAST_otu97plusV18_052025_cut.qza --output-path R2.BLAST_otu97plusV18_052025_cut_out
echo;echo “Split R1-R2 alignment complete”

# Concatenate alignments (R1 and R2)
echo;echo "Beginning concatenate R1-R2 pipeline..."
Rscript AMFconcatseqs.R
echo;echo “Concatenate R1-R2 pipeline complete”

# remove once no longer needed:
rm R1.BLAST_otu97plusV18_052025_cut.fasta
rm R1.BLAST_otu97plusV18_052025_cut.qza
rm -r R1.BLAST_otu97plusV18_052025_cut_out
rm R2.BLAST_otu97plusV18_052025_cut.fasta
rm R2.BLAST_otu97plusV18_052025_cut.qza
rm -r R2.BLAST_otu97plusV18_052025_cut_out
rm aligned_R1.BLAST_otu97plusV18_052025_cut.qza
rm aligned_R2.BLAST_otu97plusV18_052025_cut.qza

# Delete temporary folder
rm -r $SCRIPT_DIR/tmp/
