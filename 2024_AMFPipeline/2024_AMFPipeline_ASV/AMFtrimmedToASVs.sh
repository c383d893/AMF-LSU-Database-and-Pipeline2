#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=00-24:00:00
#SBATCH --output=./slurmOutputs/AMFtrimmedToASVs.out
#SBATCH --error=./slurmOutputs/AMFtrimmedToASVs.out
#SBATCH --job-name=trm2asv

# Get the working directory
SCRIPT_DIR=$1

### Activate conda
. ~/.bashrc
conda activate $C_ENV

# Define a temporary folder
mkdir $SCRIPT_DIR/tmp/
export TMPDIR=$SCRIPT_DIR/tmp/

# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2: DADA2
cat AMFdada2.R | sed "s/R1trunclen.value/$R1cutoff/" | sed "s/R2trunclen.value/$R2cutoff/" > AMFdada2withCutoffs.R

# Run DADA2 pipeline
echo;echo "Beginning DADA2 pipeline..."
Rscript AMFdada2withCutoffs.R
echo;echo "DADA2 pipeline complete. Converting output files to Qiime format..."
rm AMFdada2withCutoffs.R # remove temporary script once it's done running

# Convert ASV table from .tsv format to .biom format
biom convert -i ./dada2output/ASVtable.tsv -o ./q2files/ASVtable.biom --to-json --table-type="OTU table"

# Convert ASV table from .biom to .qza, ASV sequences from .fasta to .qza
qiime tools import --input-path ./q2files/ASVtable.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path ./q2files/ASVtable.qza

qiime tools import --input-path ./dada2output/ASVs.fasta --type 'FeatureData[Sequence]' --output-path ./q2files/ASVseqs.qza

# Convert reference database from .fasta to .qza
qiime tools import --input-path ./V18_LSUDB_052025.fasta --type 'FeatureData[Sequence]' --output-path ./q2files/AMFreferenceSeqs.qza
echo; echo "ASV table, ASV sequences, and reference sequences have been converted to .qza"

# Additional chimera removal
echo;echo "Beginning Additional chimera removal step..."
qiime vsearch uchime-denovo --i-table ./q2files/ASVtable.qza --i-sequences ./q2files/ASVseqs.qza --output-dir uchime-dn-out
qiime feature-table filter-features --i-table ./q2files/ASVtable.qza --m-metadata-file uchime-dn-out/nonchimeras.qza --o-filtered-table ./q2files/ASVtable.chim.qza
echo; echo "Additional uchime chimera removal step complete"

# ASV table cleanup 
# Remove ASVs that occur fewer than 5 times
qiime feature-table filter-features --i-table ./q2files/ASVtable.chim.qza --p-min-frequency 5 --o-filtered-table ./q2files/ASVtable_minusRareOTUs_allSamples.qza --verbose
# Remove samples with fewer than 1000 reads
qiime feature-table filter-samples --i-table ./q2files/ASVtable_minusRareOTUs_allSamples.qza --p-min-frequency 1000 --o-filtered-table ./q2files/ASVtable_clean.qza --verbose
# Update ASV representative sequences to exclude OTUs that were removed:
qiime feature-table filter-seqs --i-data ./q2files/ASVseqs.qza --i-table ./q2files/ASVtable_clean.qza --o-filtered-data ./q2files/ASVrepseqs_clean.qza --verbose
echo; echo "Rare OTUs and low-coverage samples have been removed from the dataset."

# Export ASV representative sequences into .fasta format:
# Qiime2 forces output into a file called dna-sequences.fasta . We direct it into ./q2files/
qiime tools export --input-path ./q2files/ASVrepseqs_clean.qza --output-path ./q2files/
cat ./q2files/dna-sequences.fasta > ./ASVrepseqs_clean.fasta # rename file, relabel ASVs as OTUs, and any perfect matched to database append OTU and move to working directory.

# Export ASV table into .tsv format:
# Qiime2 forces output into a .biom file called feature-table.biom . We direct it into ./q2files/
qiime tools export --input-path ./q2files/ASVtable_clean.qza --output-path ./q2files/
# Convert from .biom to .tsv:
biom convert -i ./q2files/feature-table.biom -o ./q2files/feature-table.tsv --to-tsv
cat ./q2files/feature-table.tsv | sed 's/OTU/ASV/g' > ./ASVtable_clean.tsv # rename file, any perfect matched to database append OTU and move to working directory                                              
echo;echo "export of sequences and ASV table complete"


# BLAST selection for tree building
# Replace placeholder text in R script with user-provided truncation lengths for R1 and R2: SplitR1R2
cat AMFsplitR1R2.R | sed "s/ftrunc/$R1cutoff/" | sed "s/rtrunc/$R2cutoff/" > AMFsplitR1R2withCutoffs.R

echo;echo "Beginning split R1-R2 pipeline..."
Rscript AMFsplitR1R2withCutoffs.R
echo;echo “split R1-R2 pipeline complete”
rm AMFsplitR1R2withCutoffs.R # remove temporary script once it's done running

#Make BLAST reference database
makeblastdb -in V18_LSUDB_052025_AMFONLY.fasta -input_type fasta -dbtype nucl -out V18_LSUDB_052025_AMFONLY

#Blast R1, extract top hit based on bit score, extract col 1, replace header with feature-id
blastn -db V18_LSUDB_052025_AMFONLY -query R1.ASVrepseqs_clean.fasta -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge | cut -f1 > S.BLAST.R1.ASVrepseqs_clean.txt 

#Blast R2, extract top hit based on bit score, extract col 1, replace header with feature-id
blastn -db V18_LSUDB_052025_AMFONLY -query R2.ASVrepseqs_clean.fasta -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge | cut -f1 > S.BLAST.R2.ASVrepseqs_clean.txt

# Join files, get unique OTU and add header
cat S.BLAST.R1.ASVrepseqs_clean.txt S.BLAST.R2.ASVrepseqs_clean.txt |sort| uniq | sed -e '1i\asvs' > BLAST.R1_R2.ASVrepseqs_clean_cut.tsv
echo; echo "BLAST has been completed on forward and reverse reads; each OTU with a hit saved"

# Subset .fasta and .txt based on this name list:
Rscript AMFextractBLAST.R

# remove once no longer needed:
rm S.BLAST.R1.ASVrepseqs_clean.txt 
rm S.BLAST.R2.ASVrepseqs_clean.txt
rm BLAST.R1_R2.ASVrepseqs_clean_cut.tsv 
rm R1.ASVrepseqs_clean.fasta 
rm R2.ASVrepseqs_clean.fasta 
rm -r $SCRIPT_DIR/tmp/

echo;echo "Sequences and ASV table subset to BLAST positive ASVs"

echo; echo "AMFtrimmedToASV.sh has completed running. Be sure to check the output file for error messages."
