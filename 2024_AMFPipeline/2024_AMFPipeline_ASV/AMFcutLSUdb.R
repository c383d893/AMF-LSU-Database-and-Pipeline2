####### Load required packages #######

library(Biostrings)
library(stringr)

#read in database and format seqs based on study cutoffs
LSUdb <- readDNAStringSet('V18_LSUDB_052025.fasta')
LSUdb <- as.vector(LSUdb)

#extract forward and reverse reads, then join
names<-names(LSUdb)
forwardread<-str_sub(LSUdb, 1, ftrunc) 
reverseread<-str_sub(LSUdb,-rtrunc)
new.LSUdb<-paste(forwardread,reverseread, sep="NNNNNNNNNN")

new.LSUdb<-DNAStringSet(new.LSUdb) #turn into DNAString
names(new.LSUdb)<-names #add names back in
str(new.LSUdb) #check looks ok

writeXStringSet(new.LSUdb, "V18_LSUDB_052025_cut.fasta",width=10000) #write out.
