####### Load required packages #######

library(Biostrings)
library(stringr)

#read in fasta and format seqs based on study cutoffs
SeqstoSplit <- readDNAStringSet('ASVrepseqs_clean.fasta')
SeqstoSplit <- as.vector(SeqstoSplit)

#extract forward and reverse reads, then join
names<-names(SeqstoSplit)
forwardread<-str_sub(SeqstoSplit, 1, ftrunc) 
reverseread<-str_sub(SeqstoSplit,-rtrunc)

# write out R1
forwardread<-DNAStringSet(forwardread) #turn into DNAString
names(forwardread)<-names #add names back in
str(forwardread) #check looks ok
writeXStringSet(forwardread, "R1.ASVrepseqs_clean.fasta",width=10000) #write out.

# write out R2
reverseread<-DNAStringSet(reverseread) #turn into DNAString
names(reverseread)<-names #add names back in
str(reverseread) #check looks ok
writeXStringSet(reverseread, "R2.ASVrepseqs_clean.fasta",width=10000) #write out.
