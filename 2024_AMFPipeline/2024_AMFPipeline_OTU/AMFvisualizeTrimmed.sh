#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --time=00-06:00:00
#SBATCH --output=./slurmOutputs/AMFvisualizeTrimmed.out
#SBATCH --error=./slurmOutputs/AMFvisualizeTrimmed.out
#SBATCH --job-name=AMFvisualizeTrimmed

### Activate conda
. ~/.bashrc
conda activate $C_ENV

# 1.  Create qza object of cutadapt trimmed reads
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./trimmed/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ./q2files/trimmed-paired-end.qza

# 2. Generate summaries for R1 and R2:
qiime demux summarize \
  --i-data ./q2files/trimmed-paired-end.qza \
  --o-visualization ./visualize_trimmed.qzv
