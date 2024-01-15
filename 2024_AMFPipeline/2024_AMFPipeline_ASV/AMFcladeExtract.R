####### Load required packages #######

library(ape)
library(TreeTools)

####### Determine which ASVs fall within the AMF clade #######
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from ASV subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define AMF clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='AM183920_Geosiphon_pyriformis')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832207_Acaulospora_tuberculata')$TipNumber)
  # Determine the node for the AMF clade:
  AMF.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF clade and identifying AMF ASVs:')
  # Extract the AMF clade:
  AMF.clade <- extract.clade(mytree.rooted.preorder,AMF.node)
  # plotTree(AMF.clade)
  # str(AMF.clade) # check how many tips
  ### Extract names of ASVs from the AMF clade:
  AMFnames<-c(AMF.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFnames[grep("ASV",AMFnames)]
  if (tree==treeFiles[1]) { AMF.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMF.ASVs <- append(AMF.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMF.node, AMF.clade, AMFnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-AMF ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune:
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFASVs <- subset(ASVs,X.ASV.ID %in% AMF.ASVs) # filter out non-AMF ASVs using the list we generated above
otherASVs <- subset(ASVs,!(X.ASV.ID %in% AMF.ASVs)) # filter out non-AMF ASVs using the list we generated above
# Write to file:
write.table(AMFASVs,file='./ASVtable_clean_AMFonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
write.table(otherASVs,file='./ASVtable_clean_nonAMF.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune:
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFseqs <- seqs[AMF.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMF.ASVs)]
# Write to file:
write.FASTA(AMFseqs,'./ASVrepseqs_clean_AMFonly.fasta')
write.FASTA(otherseqs,'./ASVrepseqs_clean_nonAMF.fasta')

print('AMF - LSU pipeline complete.')
print('AMF-only ASV table and sequences have been written to file.')
print('non-AMF ASV table and sequences have also been written to file.')