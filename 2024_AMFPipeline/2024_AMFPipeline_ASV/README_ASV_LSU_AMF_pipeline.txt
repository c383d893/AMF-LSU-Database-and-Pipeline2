######### *** README *** ##########

###################################
######## IMPORTANT NOTES ##########
###################################

##IMPORTANT NOTES:
# 1. AMF LSU Pipeline
# 2. Primers : LROR (ACCCGCTGAACTTAAGC); FLR2 (TCGTTTAAAGCCATTACGTC)
# 3. Assumes demultiplexed files with names *R1_001.fastq.gz and *R2_001.fastq.gz
# 4. All paths are relative to the working directory (.); follow directory structure

###################################
####### REQUIRED PROGRAMS #########
###################################

# 1. conda/miniconda/miniforge
# 2. AMF conda environment: contains all programs and packages required

##################################
####### INSTALLING CONDA #########
##################################

# These instructions will create a standard conda installation in your home folder 
# i.e. it is fully contained within your home directory and should not affect any 
# other users of your cluster. 

# 1. Download the miniforge installer

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

# 2. Run the installer: accept all defaults

bash Miniforge3-$(uname)-$(uname -m).sh

# IMPORTANT: Do you wish to update your shell profile to automatically initialize conda? respond 'no'

# 3. Register conda commands

~/miniforge3/condabin/conda init

# 4. Source your bash.rc file to load the conda settings

source ~/.bashrc

# 5. Turn of conda base environment: not doing this can cause software version conflicts

conda config --set auto_activate_base false

###########################################
#### SET UP AMF LSU CONDA ENVIRONMENT #####
###########################################

# 1. Check that conda is installed correctly: if installed correctly,
# you will see a list of environments and their locations

conda env list

# 2. Install the conda AMF LSU pipeline into a new environment: done only once.  
#    We use a custom .yml file that sets up a qiime2 (2024.2) environment with 
#    additional packages (r-treetools).

conda env create -n "amf_pipeline" --file "amf_pipeline_requirements.yml"

###########################################
########### DIRECTORY STRUCTURE ###########
###########################################

# Main working directory:
mkdir myWD/ 
cd myWD/

# Subdirectory for raw sequence data:
mkdir ./raw/

# Subdirectory for primer-trimmed sequence data:
mkdir ./trimmed/

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

###########################################
####### REQUIRED FILES & LOCATIONS ########
###########################################

# Raw demultiplexed study sequences (paired-end):
./raw/*R1_001.fastq.gz
./raw/*R2_001.fastq.gz

# Reference database:
./V16_LSUDB_2024.fasta

# AMF only reference database:
./V16_LSUDB_2024_AMFONLY.fasta

# Bash scripts to run on SLURM
./AMFvisualizeRaw.sh
./AMFvisualizeTrimmed.sh
./AMFseqTrimParallel.sh
./AMFtrimmedToASVs.sh
./AMFalignseqs.sh
./AMFbuildTree.sh
./AMFlaunchTrees.sh

# R scripts
./AMFdada2.R
./AMFcladeExtract.R
./AMFcladeExtract_family.R
./AMFconcatseqs.R
./AMFalignseqs.R
./AMFcutLSUdb.R
./AMFextractBLAST.R
./AMFsplitR1R2 

###########################################
######### SET CONDA ENVIRONMENT ###########
###########################################

# 1. Set name/location of conda environment
#    IMPORTANT: run each time you log in to the cluster

C_ENV="amf_pipeline"

#########################################
########### VIZUALIZE RAW DATA ##########
#########################################

# 1. Vizualize raw reads
    
sbatch --export=C_ENV=$C_ENV \
    AMFvisualizeRaw.sh  \
    $(readlink -f $(dirname AMFvisualizeRaw.sh))

# 2. Download outputs (./visualize_raw.qzv)

# 3. View visualize_raw.qzv on https://view.qiime2.org/ 

#########################################
########### INITIAL PROCESSING: #########
###### REMOVE PRIMERS & VISUALIZE #######
#########################################

# 1. Trim adapters and primers
#    IMPORTANT: update this to specify your number of samples: count: ls | wc -l
#    i.e. --array=1-10%8 becomes --array=1-120%8 for 120 samples

sbatch --array=1-10%8 --export=C_ENV=$C_ENV AMFseqTrimParallel.sh

# 2. Vizualize trimmed reads

sbatch --export=C_ENV=$C_ENV \
    AMFvisualizeTrimmed.sh \
    $(readlink -f $(dirname AMFvisualizeTrimmed.sh))

# 3. Download outputs (./visualize_trimmed.qzv)

# 4. View visualize_trimmed.qzv on https://view.qiime2.org/ 

# 5. Choose truncation lengths for R1 and for R2 -- 
#    where does sequence quality start to drop off?
#    Specify your cutoffs in the next sbatch command (below)

#########################################
########### QUALITY FILTER: #############
####### DENOISE, PRE-BLAST FILTER #######
######### CHIMERA REMOVAL ###############
#########################################

# 1. Specify truncation lengths for R1 and R2:
#    Determined in step 5 above from your dataset

sbatch --export=R1cutoff=170,R2cutoff=140,C_ENV=$C_ENV \
    AMFtrimmedToASVs.sh \
    $(readlink -f $(dirname AMFtrimmedToASVs.sh))

#########################################
########### PREPARE FOR RAXML: ##########
######## ALIGN DB & STUDY SEQS ##########
#########################################

# 1. Align R1 and R2 reads with the reference database and
#    then concatenate resulting alignments

sbatch --export=R1cutoff=170,R2cutoff=140,C_ENV=$C_ENV \
    AMFalignseqs.sh \
    $(readlink -f $(dirname AMFalignseqs.sh))

#########################################
###### PLACE SEQS ON BACKBONE TREE ######
#########################################

# 1. Place representative sequences onto the backbone/reference tree using
#    RAxML: This script will determine how many seqs need to be placed, and
#    launch a batch script with an appropriate number of tasks to be processed
#    in parallel.

bash AMFlaunchTrees.sh \
    $C_ENV \
    $(readlink -f $(dirname AMFlaunchTrees.sh))

# 2. Check that all trees finished

cd slurmOutputs/
find . -name "buildTree*.out" | xargs grep -E 'DUE TO TIME LIMIT'
cd ..

#########################################
###### EXTRACT SEQS IN AMF PHYLUM #######
#########################################

# 1. Determine which seqs fall within the AMF clade, 
#    and make subsets of the rep-seqs and seqs table
#    for AMF-seqs and non-AMF-seqs

sbatch --export=C_ENV=$C_ENV AMFcladeExtract.sh

#########################################
###### EXTRACT SEQS IN AMF FAMILIES #####
#########################################

# 1. Create subsets for each AMF clade:

sbatch --export=C_ENV=$C_ENV AMFcladeExtract_family.sh