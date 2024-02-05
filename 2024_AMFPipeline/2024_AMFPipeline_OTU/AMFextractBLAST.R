####### Load required packages #######

library(Biostrings)
library(stringr)
library(ape)

# Load OTU table and prune:
print('Loading qiime2 output to be subset - OTU table')
otus <- read.delim('./otu97table_clean.tsv',sep='\t',skip=1,header=TRUE)
BLAST.OTUs <- read.delim('./BLAST.R1_R2.otu97repseqs_clean_cut.tsv',sep='\t',header=TRUE)
BLASTotus <- subset(otus,X.OTU.ID %in% BLAST.OTUs$otus) # filter out non-AMF OTUs using the list we generated above
otherotus <- subset(otus,!(X.OTU.ID %in% BLAST.OTUs$otus)) # filter out non-AMF OTUs using the list we generated above

# Write to file:
write.table(BLASTotus,file='./otu97table_clean_BLASTonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
write.table(otherotus,file='./otu97table_clean_nonBLAST.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

# Load OTU rep seqs and prune:
print('Loading qiime2 output to be pruned - OTU rep. seqs')
otus <- readDNAStringSet('./otu97repseqs_clean.fasta')
BLAST.OTUs <- read.delim('./BLAST.R1_R2.otu97repseqs_clean_cut.tsv',sep='\t',header=TRUE)
BLASTseqs <- otus[BLAST.OTUs$otus]
otherseqs <- otus[BLAST.OTUs$otus]

# Write to file:
writeXStringSet(BLASTseqs, "./otu97repseqs_clean_BLAST.fasta",width=10000)

print('BLAST - LSU pipeline complete.')
print('BLAST - only OTU table and sequences have been written to file.')
