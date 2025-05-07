####### Load required packages #######

library(Biostrings)
library(stringr)

#read in fasta and format seqs based on study cutoffs
SeqstoAlign <- readDNAStringSet('BLAST_otu97plusV18_052025_cut.fasta')
SeqstoAlign <- as.vector(SeqstoAlign)

#extract forward and reverse reads, then join
names<-names(SeqstoAlign)
forwardread<-str_sub(SeqstoAlign, 1, ftrunc) 
reverseread<-str_sub(SeqstoAlign,-rtrunc)

# write out R1
forwardread<-DNAStringSet(forwardread) #turn into DNAString
names(forwardread)<-names #add names back in
str(forwardread) #check looks ok
writeXStringSet(forwardread, "R1.BLAST_otu97plusV18_052025_cut.fasta",width=10000) #write out.

# write out R2
reverseread<-DNAStringSet(reverseread) #turn into DNAString
names(reverseread)<-names #add names back in
str(reverseread) #check looks ok
writeXStringSet(reverseread, "R2.BLAST_otu97plusV18_052025_cut.fasta",width=10000) #write out.
