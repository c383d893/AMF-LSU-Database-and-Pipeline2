#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --time=00-06:00:00
#SBATCH --output=./slurmOutputs/AMFvisualizeRaw.out
#SBATCH --error=./slurmOutputs/AMFvisualizeRaw.out
#SBATCH --job-name=AMFvisualizeRaw

### Activate conda
. ~/.bashrc
conda activate $C_ENV

# 1.  Create qza object of raw reads
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./raw \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ./q2files/raw_files.qza

# 2. Generate summaries for R1 and R2:
qiime demux summarize \
  --i-data ./q2files/raw_files.qza \
  --o-visualization ./visualize_raw.qzv