
### *** README *** ###

##IMPORTANT NOTES:
#1. AMF LSU Pipeline
#2. primers : LROR (ACCCGCTGAACTTAAGC); FLR2 (TCGTTTAAAGCCATTACGTC)
#3. assumes demultiplexed files with names *R1_001.fastq.gz and *R2_001.fastq.gz
#4. Adapted to run on the ETH Euler cluster (tested with 10 samples from the STRI project)

###################################
####### Required programs #########
###################################

# QIIME2/2019.10.0
# cutadapt/4.4 (python/3.11.2)
# fastqc/0.11.9
# R/4.2.2
# raxml/8.2.12
# blast-plus/2.12.0

###################################
####### Required R packages #######
###################################

# dada2
# ShortRead
# Biostrings
# stringr
# ape
# TreeTools

###################################
####### Directory structure #######
###################################

# Main working directory:
# could be named anything. All scripts use paths relative to this directory
mkdir myWD/ 
cd myWD/

# Subdirectory for raw sequence data:
mkdir ./raw/

# Subdirectory for primer-trimmed sequence data:
mkdir ./trimmed/

# Subdirectory for fastQC output files:
mkdir ./fastQCoutput

# Subdirectory for truncated + quality-filtered + concatenated sequence data:
mkdir ./filtered/

# Subdirectory for DADA2 output
mkdir ./dada2output

# Subdirectory for Qiime2 files
mkdir ./q2files

# Subdirectory for RAxML files
mkdir ./RAxMLfiles

# Subdirectory for all Output/Error files
mkdir ./slurmOutputs

###################################
#### Required files & locations ###
###################################

# *** All paths are relative to the working directory (.)

# Raw demultiplexed study sequences (paired-end):
./raw/*R1_001.fastq.gz
./raw/*R2_001.fastq.gz

# Reference database:
./v16_LSUDB_2024.fasta

# AMF only reference database:
./v16_LSUDB_2024_AMFONLY.fasta

# Bash scripts to run on SLURM
./AMFseqTrimParallel.sh   # this one will trim primers from the raw sequences
./AMFtrimmedToASVs.sh   # this one will quality-filter the primer-trimmed sequences, truncate reads to fixed lengths based on user input, denoise, concatenate FWD and REV reads into ASVs. Then, this separated FWD and REV to BLAST against AMF in the reference database and subsets the ASV seqs and table to include only BLAST positive ASVS
./AMFalignseqs.sh # align sequence to database sequences (FWD and REV separately), then concatenate
./AMFbuildTree.sh   # this one places environmental (OTU) sequences onto the backbone tree using RAxML
./AMFlaunchTrees.sh   # this one optimizes AMFbuildTree.sh for parallel processing for a particular dataset (just decides how many parallel jobs need to be run for a given set of ASVs)

# R scripts:
./AMFdada2.R    # this gets called by AMFtrimmedToASVs.sh
./AMFcladeExtract.R    # outputs ASV table and rep-seqs file subsetted to include AMF only
./AMFcladeExtract_family.R # outputs ASV table and rep-seqs file subsetted to include AMF for each family
./AMFconcatseqs.R # concatenates aligned seqs
./AMFalignseqs.R # aligns R1 and R2 separately 
./AMFcutLSUdb.R # Cuts database to match study seq cutoffs
./AMFextractBLAST.R # extract blast ASVs
./AMFsplitR1R2 # split R1 and R2 before blasting

###################################
####### fastQC on raw data ########
###################################

module load gcc/8.2.0 fastqc/0.11.9

zcat ./raw/*R1_001.fastq.gz | fastqc -o ./fastQCoutput stdin:R1_raw 
zcat ./raw/*R2_001.fastq.gz | fastqc -o ./fastQCoutput stdin:R2_raw

###################################
####### Initial processing: #######
####### Remove primers and ########
####### visualize quality #########
###################################

# IMPORTANT: update this to specify your number of samples: count: ls | wc -l
# i.e. --array=1-10%8 becomes --array=1-120%8 for 120 sampmles

sbatch --array=1-10%8 AMFseqTrimParallel.sh		# Samples will be trimmed in parallel (up to 8 simultaneously)

### After AMFseqTrimParallel.sh is done running: 

# 1. Fire up fastQC:
module load gcc/8.2.0 fastqc/0.11.9

# 2. Generate fastQC summaries for R1 and R2:

zcat ./trimmed/*R1_trimmed.fastq.gz | fastqc -o ./fastQCoutput stdin:R1 
zcat ./trimmed/*R2_trimmed.fastq.gz | fastqc -o ./fastQCoutput stdin:R2

# 3. Download fastqc outputs (./fastQCoutput/*.html)
# 4. View on web browser: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# 5. Choose truncation lengths for R1 and for R2 -- 
#    where does sequence quality start to drop off?
#    Specify your cutoffs in the next sbatch command (below)

###################################
######### Quality filter,##########
######### denoise, BLAST,##########
######### make ASV table ##########
###################################

### Specify truncation lengths for R1 and R2:
# These (110 and 105) are default values; you should change them
# Also change the working directory to current directory
# to fit your dataset (based on FastQC results)

sbatch --export=R1cutoff=150,R2cutoff=110 AMFtrimmedToASVs.sh $(cd -P -- "$(dirname -- "$0")" && pwd -P)

# ^ this bash script will do everything from dada2 pipeline (make ASV table) to 
# Determining ASVS, to
# filtering the ASV table, to
# subsetting only BLAST positive ASVS and filtering seqs and table in preparation for RAXML.

###################################
##### Align DB and study seqs #####
##### R1, R2, concatenate #########
###################################

# This script will align R1 and R2 reads with the reference database and then concatenate resulting alignments

sbatch --export=R1cutoff=150,R2cutoff=110 AMFalignseqs.sh $(cd -P -- "$(dirname -- "$0")" && pwd -P)

###################################
####### Place ASV sequences #######
####### on the backbone tree ######
###################################

# Place ASV representative sequences onto the backbone/reference tree using RAxML:
# This script will determine how many ASVs need to be placed, and launch a batch
# script with an appropriate number of tasks to be processed in parallel.

bash AMFlaunchTrees.sh

#Check that all trees finished

cd slurmOutputs/
find . -name "buildTree*.out" | xargs grep -E 'DUE TO TIME LIMIT'

###################################
#### Separate AMF from non-AMF ####
###################################

### *** NOTE: *** ###
## before proceeding, the following packages 
## must be installed in your R environment:
# ape
# TreeTools

# After all tasks are finished running for AMFbuildTree.sh,
# run this R script to determine which ASVs fall within the
# AMF clade, and make subsets of the rep-seqs and ASV table
# for AMF-ASVs and non-AMF-ASVs

module load gcc/8.2.0 r/4.2.2

Rscript AMFcladeExtract.R

##################################
#### Separate AMF by family ######
##################################
  
Rscript AMFcladeExtract_family.R
