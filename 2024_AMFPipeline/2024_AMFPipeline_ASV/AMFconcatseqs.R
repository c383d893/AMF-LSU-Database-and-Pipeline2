####### Load required packages #######

library(Biostrings)
library(stringr)

#read in fasta to get names
SeqstoAlign <- readDNAStringSet('BLAST_ASVplusV15_3.23.21_cut.fasta')
SeqstoAlign <- as.vector(SeqstoAlign)

#extract forward and reverse reads, then join
names<-names(SeqstoAlign)

# bring alignments here and concatenate:

#read in database and format seqs based on study cutoffs
forward.aligned <- readDNAStringSet('R1.BLAST_ASVplusV15_3.23.21_cut_out/aligned-dna-sequences.fasta')
reverse.aligned <- readDNAStringSet('R2.BLAST_ASVplusV15_3.23.21_cut_out/aligned-dna-sequences.fasta')

#extract forward and reverse reads, then join
joined.aligned<-paste(forward.aligned,reverse.aligned, sep="NNNNNNNNNN") # don't need Ns anymore, but doesn't hurt

joined.aligned<-DNAStringSet(joined.aligned) #turn into DNAString
names(joined.aligned)<-names #add names back in
str(joined.aligned) #check looks ok

writeXStringSet(joined.aligned, "aligned_BLAST_ASVplusV15_3.23.21_cut.fasta",width=10000) #write out.
