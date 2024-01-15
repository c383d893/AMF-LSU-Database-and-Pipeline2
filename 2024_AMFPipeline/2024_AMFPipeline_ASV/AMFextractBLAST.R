####### Load required packages #######

library(Biostrings)
library(stringr)
library(ape)

# Load ASV table and prune:
print('Loading qiime2 output to be subset - ASV table')
ASVs <- read.delim('./ASVtable_clean.tsv',sep='\t',skip=1,header=TRUE)
BLAST.ASVs <- read.delim('./BLAST.R1_R2.ASVrepseqs_clean_cut.tsv',sep='\t',header=TRUE)
BLASTASVs <- subset(ASVs,X.ASV.ID %in% BLAST.ASVs$asvs) # filter out non-AMF ASVs using the list we generated above
otherASVs <- subset(ASVs,!(X.ASV.ID %in% BLAST.ASVs$asvs)) # filter out non-AMF ASVs using the list we generated above

# Write to file:
write.table(BLASTASVs,file='./ASVtable_clean_BLASTonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
write.table(otherASVs,file='./ASVtable_clean_nonBLAST.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

# Load ASV rep seqs and prune:
print('Loading qiime2 output to be pruned - ASV rep. seqs')
ASVs <- readDNAStringSet('./ASVrepseqs_clean.fasta')
BLAST.ASVs <- read.delim('./BLAST.R1_R2.ASVrepseqs_clean_cut.tsv',sep='\t',header=TRUE)
BLASTseqs <- ASVs[BLAST.ASVs$asvs]
otherseqs <- ASVs[BLAST.ASVs$asvs]

# Write to file:
writeXStringSet(BLASTseqs, "./ASVrepseqs_clean_BLAST.fasta",width=10000)

print('BLAST - LSU pipeline complete.')
print('BLAST - only ASV table and sequences have been written to file.')
