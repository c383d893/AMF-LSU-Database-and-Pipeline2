
### *** README *** ###

##IMPORTANT NOTES:
#1. AMF LSU Pipeline
#2. primers : LROR (ACCCGCTGAACTTAAGC); FLR2 (TCGTTTAAAGCCATTACGTC)
#3. assumes demultiplexed files with names *R1_001.fastq.gz and *R2_001.fastq.gz
#4. Adapted to run on the ETH Euler cluster (tested with 10 samples from the STRI project)

###################################
####### Required programs #########
###################################

# conda/miniconda/miniforge
# QIIME2/2024.2


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

################################
#### Setup conda environment ###
################################

# 1. You must have conda installed on your cluster
#    To check if its installed correctly run the
#    following command.
conda env list

# If conda is installed correctly then you will
# see a list of environments and thier locations
# at a minimum you should have a "base" environment

# 2. Install the conda pipeline into a new
#    environment. We use a custom yml file
#    that sets up a qiime2 (2024.2) environment
#    with additional packages (r-treetools).
#    The name of the environment is by default
#    "amf-pipeline" but can be changes in the
#    command (-n "amf-pipeline").
conda env create -n "amf_pipeline" --file "amf_pipeline_requirements.yml"

#########################################
####### Activate Qiime2, this will ######
####### be our active environment #######
####### for most subsequent steps #######
#########################################

# 1. Set name/location of conda environment
CONDA_PIPELINE_ENV="amf_pipeline"

# 2. Activate the amf pipeline environment
conda activate ${CONDA_PIPELINE_ENV}

#########################################
####### View quality of raw data ########
#########################################

# 1.  Create qza object of raw reads
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path raw \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

# 2. Generate summaries for R1 and R2:
qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qzv

# 3. View demux-paired-end.qzv either on https://view.qiime2.org/ 
#    or using a local install of qiime2 and the command:
#    qiime tools view demux-paired-end.qzv

###################################
####### Initial processing: #######
####### Remove primers and ########
####### visualize quality #########
###################################

# IMPORTANT: update this to specify your number of samples: count: ls | wc -l
# i.e. --array=1-10%8 becomes --array=1-120%8 for 120 sampmles

sbatch --array=1-10%8 AMFseqTrimParallel.sh		# Samples will be trimmed in parallel (up to 8 simultaneously)

### After AMFseqTrimParallel.sh is done running:

# 1. Create qza object of trimmed reads
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path trimmed \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path trimmed-paired-end.qza

# 2. Generate summaries for R1 and R2:
qiime demux summarize \
  --i-data trimmed-paired-end.qza \
  --o-visualization trimmed-paired-end.qzv

# 3. Download  outputs (./trimmed-paired-end.qzv)

# 4. View trimmed-paired-end.qzv either on https://view.qiime2.org/ 
#    or using a local install of qiime2 and the command 
#    "qiime tools view trimmed-paired-end.qzv"

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
# to fit your dataset

sbatch --export=R1cutoff=150,R2cutoff=110,CONDA_PIPELINE_ENV AMFtrimmedToASVs.sh $(cd -P -- "$(dirname -- "$0")" && pwd -P)

# ^ this bash script will do everything from dada2 pipeline (make ASV table) to 
# Determining ASVS, to
# filtering the ASV table, to
# subsetting only BLAST positive ASVS and filtering seqs and table in preparation for RAXML.

###################################
##### Align DB and study seqs #####
##### R1, R2, concatenate #########
###################################

# This script will align R1 and R2 reads with the reference database and then concatenate resulting alignments

sbatch --export=R1cutoff=150,R2cutoff=110,CONDA_PIPELINE_ENV AMFalignseqs.sh $(cd -P -- "$(dirname -- "$0")" && pwd -P)

###################################
####### Place ASV sequences #######
####### on the backbone tree ######
###################################

# Place ASV representative sequences onto the backbone/reference tree using RAxML:
# This script will determine how many ASVs need to be placed, and launch a batch
# script with an appropriate number of tasks to be processed in parallel.

bash --export=CONDA_PIPELINE_ENV AMFlaunchTrees.sh

# Check that all trees finished

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

Rscript AMFcladeExtract.R

##################################
#### Separate AMF by family ######
##################################
  
Rscript AMFcladeExtract_family.R
