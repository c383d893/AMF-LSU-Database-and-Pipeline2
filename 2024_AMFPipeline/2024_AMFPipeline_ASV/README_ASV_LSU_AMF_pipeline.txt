
### *** README *** ###

##IMPORTANT NOTES:
# 1. AMF LSU Pipeline
# 2. primers : LROR (ACCCGCTGAACTTAAGC); FLR2 (TCGTTTAAAGCCATTACGTC)
# 3. assumes demultiplexed files with names *R1_001.fastq.gz and 
#    *R2_001.fastq.gz
# 4. Adapted to run on the ETH Euler cluster (tested with 10 samples from the
#    STRI project)

###################################
####### Required programs #########
###################################

# 1. conda/miniconda/miniforge
# 2. QIIME2/2024.2
# 3. R package TreeTools

##################################
####### Installing Conda #########
##################################

# 1. Any version of the conda installer should work, but this pipeline
#    is tested with the miniforge installer.
#    Simply download the installer

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

# 2. Then run the installer. Generally running the installer as a normal user
#    and accepting all the defaults should work fine.

bash Miniforge3-$(uname)-$(uname -m).sh

# 3. You will likely want to run conda init so that conda commands will work in
#    your normal shell. You need to specify the full path to the conda
#    executable this one time. Afterwards the conda commands should be
#    registered.

~/miniforge3/condabin/conda init

# 4. You can then source your bash.rc file to load the conda settings. 
#    You can also just log out and log back in. Until you perform one of those
#    two steps, you will unable to run conda commands without specifiying thier
#    path

source ~/.bashrc

# 5. Finally, by default conda loads its "base" environment on login if you ran 
#    conda init. You may not want the base environemnt loaded by default,
#    at it could potentially cause software version conflicts. You can turn of
#    this feature by changing the conda config settings.

conda config --set auto_activate_base false

# 6. You can turn off the base environment in your current session simply by
#    typing.

conda deactivate

# Note, these instructions will create a standard conda installation in you 
# home folder. This installation is fully contained within your home directory 
# and should not affect any other users of your cluster. If you are using a 
# version of conda managed by your cluster administrator you may need to adjust 
# some of the subsequent commands. 
# i.e. conda activate ... >> module load conda && conda activate ...

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

# 1 You must have conda successfully installed on your cluster
#    To check if its installed correctly run the following command.

conda env list

# If conda is installed correctly then you will
# see a list of environments and thier locations
# at a minimum you should have a "base" environment

# 2. If you have not already, install the conda 
#    pipeline into a new environment. This
#    only need to be done once per user.
#    We use a custom yml file that sets up 
#    a qiime2 (2024.2) environment with 
#    additional packages (r-treetools).
#    The name of the environment is by default
#    "amf-pipeline" but can be changed in the
#    argument (-n "amf-pipeline").

conda env create -n "amf_pipeline" --file "amf_pipeline_requirements.yml"

###########################################
####### Save name of pipeline, this #######
####### be our active environment   #######
####### for most subsequent jobs    #######
###########################################

# 1. Set name/location of conda environment
C_ENV="amf_pipeline"

# 2. You may need to set a default partition if your cluster
#    does not specify one
# export SBATCH_PARTITION=sixhour

#########################################
####### View quality of raw data ########
#########################################

# 1. Vizialize raw reads
sbatch --export=C_ENV=$C_ENV AMFvisualizeRaw.sh

# 2. Download outputs (./visualize_raw.qzv)

# 3. View visualize_raw.qzv either on https://view.qiime2.org/ 
#    or using a local install of qiime2 and the command:
#    qiime tools view visualize_raw.qzv

###################################
####### Initial processing: #######
####### Remove primers and ########
####### visualize quality #########
###################################

# 1. Trim adapters and primers
#    IMPORTANT: update this to specify your number of samples: count: ls | wc -l
#    i.e. --array=1-10%8 becomes --array=1-120%8 for 120 sampmles

sbatch --array=1-10%8 --export=C_ENV=$C_ENV AMFseqTrimParallel.sh
# Samples will be trimmed in parallel (up to 8 simultaneously)

# 2. Visualization of trimmed sequence quality
#    After AMFseqTrimParallel.sh is done run:
sbatch --export=C_ENV=$C_ENV AMFvisualizeTrimmed.sh

# 3. Download outputs (./visualize_trimmed.qzv)

# 4. View visualize_trimmed.qzv either on https://view.qiime2.org/ 
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

# 1. Specify truncation lengths for R1 and R2:
#    These (110 and 105) are default values; you should change them
#    Also change the working directory to current directory
#    to fit your dataset

sbatch --export=R1cutoff=110,R2cutoff=105,C_ENV=$C_ENV \
    AMFtrimmedToASVs.sh \
    $(cd -P -- "$(dirname -- "$0")" && pwd -P)

# ^ this bash script will do everything from dada2 pipeline (make ASV table) to 
# Determining ASVS, to
# filtering the ASV table, to
# subsetting only BLAST positive ASVS and filtering seqs and table in 
# preparation for RAXML.

###################################
##### Align DB and study seqs #####
##### R1, R2, concatenate #########
###################################

# 1. This script will align R1 and R2 reads with the reference database and
#    then concatenate resulting alignments

sbatch --export=R1cutoff=110,R2cutoff=105,C_ENV=$C_ENV \
    AMFalignseqs.sh \
    $(readlink -f $(dirname AMFalignseqs.sh))

###################################
####### Place ASV sequences #######
####### on the backbone tree ######
###################################

# 1. Place ASV representative sequences onto the backbone/reference tree using
#    RAxML: This script will determine how many ASVs need to be placed, and
#    launch a batch script with an appropriate number of tasks to be processed
#    in parallel.

bash AMFlaunchTrees.sh \
    $C_ENV \
    $(readlink -f $(dirname AMFalignseqs.sh))

# 2. Check that all trees finished

cd slurmOutputs/
find . -name "buildTree*.out" | xargs grep -E 'DUE TO TIME LIMIT'
cd ..

###################################
#### Separate AMF from non-AMF ####
###################################

# 1. After all tasks are finished running for AMFbuildTree.sh,
#    run this R script to determine which ASVs fall within the
#    AMF clade, and make subsets of the rep-seqs and ASV table
#    for AMF-ASVs and non-AMF-ASVs

sbatch --export=C_ENV=$C_ENV AMFcladeExtract.sh

##################################
#### Separate AMF by family ######
##################################

# 1. Create subsets for each AMF clade:

sbatch --export=C_ENV=$C_ENV AMFcladeExtract_family.sh

