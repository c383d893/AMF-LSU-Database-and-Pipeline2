####### Load required packages #######

library(ape)
library(TreeTools)

#################################################################################
####### Determine which OTUs fall within the Acaulosporaceae family clade #######
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define ACA clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='FN547517_Acaulospora_laevis')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832207_Acaulospora_tuberculata')$TipNumber)
  # Determine the node for the ACA clade
  AMFACA.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF ACA clade and identifying OTUs:')
  # Extract the ACA clade
  AMFACA.clade <- extract.clade(mytree.rooted.preorder,AMFACA.node)
  ### Extract names of OTUs from the ACA clade
  AMFACAnames<-c(AMFACA.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFACAnames[grep("OTU",AMFACAnames)]
  if (tree==treeFiles[1]) { AMFACA.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFACA.OTUs <- append(AMFACA.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF ACA OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFACA.node, AMFACA.clade, AMFACAnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-ACA OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)

AMFACAotus <- subset(otus,X.OTU.ID %in% AMFACA.OTUs) # filter out non-ACA OTUs using the list we generated above
# Write to file
write.table(AMFACAotus,file='./otu97table_clean_AMFACAonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFACAseqs <- seqs[AMFACA.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFACA.OTUs)]
# Write to file
write.FASTA(AMFACAseqs,'./otu97repseqs_clean_AMFACAonly.fasta')

print('AMF-ACA LSU pipeline complete.')
print('AMF-ACA OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Ambisporaceae family clade ##########
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define AMB clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MT832166_Ambispora_leptoticha')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='JF439210_Ambispora_gerdemannii')$TipNumber)
  # Determine the node for the AMB clade
  AMFAMB.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF AMB clade and identifying OTUs:')
  # Extract the AMB clade
  AMFAMB.clade <- extract.clade(mytree.rooted.preorder,AMFAMB.node)
  ### Extract names of OTUs from the AMB clade
  AMFAMBnames<-c(AMFAMB.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFAMBnames[grep("OTU",AMFAMBnames)]
  if (tree==treeFiles[1]) { AMFAMB.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFAMB.OTUs <- append(AMFAMB.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF AMB OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFAMB.node, AMFAMB.clade, AMFAMBnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-AMB OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFAMBotus <- subset(otus,X.OTU.ID %in% AMFAMB.OTUs) # filter out non-AMB OTUs using the list we generated above
# Write to file
write.table(AMFAMBotus,file='./otu97table_clean_AMFAMBonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFAMBseqs <- seqs[AMFAMB.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFAMB.OTUs)]
# Write to file
write.FASTA(AMFAMBseqs,'./otu97repseqs_clean_AMFAMBonly.fasta')

print('AMF-AMB LSU pipeline complete.')
print('AMF-AMB OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Archaeosporaceae family clade #########
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define ARC clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MT832160_Archaeospora_trappei')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832165_Archaeospora_schencki')$TipNumber)
  # Determine the node for the ARC clade
  AMFARC.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF ARC clade and identifying OTUs:')
  # Extract the ARC clade
  AMFARC.clade <- extract.clade(mytree.rooted.preorder,AMFARC.node)
  ### Extract names of OTUs from the ARC clade
  AMFARCnames<-c(AMFARC.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFARCnames[grep("OTU",AMFARCnames)]
  if (tree==treeFiles[1]) { AMFARC.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFARC.OTUs <- append(AMFARC.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF ARC OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFARC.node, AMFARC.clade, AMFARCnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-ARC OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFARCotus <- subset(otus,X.OTU.ID %in% AMFARC.OTUs) # filter out non-ARC OTUs using the list we generated above
# Write to file
write.table(AMFARCotus,file='./otu97table_clean_AMFARConly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFARCseqs <- seqs[AMFARC.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFARC.OTUs)]
# Write to file
write.FASTA(AMFARCseqs,'./otu97repseqs_clean_AMFARConly.fasta')

print('AMF-ARC LSU pipeline complete.')
print('AMF-ARC OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Claroideoglomeraceae family clade #####
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define CLA clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='AJ972464.1_Albahypha_drummondii')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832172_Claroideoglomus_etunicatum')$TipNumber)
  # Determine the node for the CLA clade
  AMFCLA.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF CLA clade and identifying OTUs:')
  # Extract the CLA clade
  AMFCLA.clade <- extract.clade(mytree.rooted.preorder,AMFCLA.node)
  ### Extract names of OTUs from the CLA clade
  AMFCLAnames<-c(AMFCLA.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFCLAnames[grep("OTU",AMFCLAnames)]
  if (tree==treeFiles[1]) { AMFCLA.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFCLA.OTUs <- append(AMFCLA.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF CLA OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFCLA.node, AMFCLA.clade, AMFCLAnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-CLA OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFCLAotus <- subset(otus,X.OTU.ID %in% AMFCLA.OTUs) # filter out non-CLA OTUs using the list we generated above
# Write to file
write.table(AMFCLAotus,file='./otu97table_clean_AMFCLAonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFCLAseqs <- seqs[AMFCLA.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFCLA.OTUs)]
# Write to file
write.FASTA(AMFCLAseqs,'./otu97repseqs_clean_AMFCLAonly.fasta')

print('AMF-CLA LSU pipeline complete.')
print('AMF-CLA OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Diversisporaceae family clade #########
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define DIV clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MG459208.1_Desertispora_omaniana')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832219_Diversispora_trimurales')$TipNumber)
  # Determine the node for the DIV clade
  AMFDIV.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF DIV clade and identifying OTUs:')
  # Extract the DIV clade
  AMFDIV.clade <- extract.clade(mytree.rooted.preorder,AMFDIV.node)
  ### Extract names of OTUs from the DIV clade
  AMFDIVnames<-c(AMFDIV.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFDIVnames[grep("OTU",AMFDIVnames)]
  if (tree==treeFiles[1]) { AMFDIV.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFDIV.OTUs <- append(AMFDIV.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF DIV OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFDIV.node, AMFDIV.clade, AMFDIVnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-DIV OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFDIVotus <- subset(otus,X.OTU.ID %in% AMFDIV.OTUs) # filter out non-DIV OTUs using the list we generated above
# Write to file
write.table(AMFDIVotus,file='./otu97table_clean_AMFDIVonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFDIVseqs <- seqs[AMFDIV.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFDIV.OTUs)]
# Write to file
write.FASTA(AMFDIVseqs,'./otu97repseqs_clean_AMFDIVonly.fasta')

print('AMF-DIV LSU pipeline complete.')
print('AMF-DIV OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Gigasporaceae family clade ############
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define GIG clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MT832220_Scutellospora_calospora')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832238_Gigaspora_gigantea')$TipNumber)
  # Determine the node for the GIG clade
  AMFGIG.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF GIG clade and identifying OTUs:')
  # Extract the GIG clade
  AMFGIG.clade <- extract.clade(mytree.rooted.preorder,AMFGIG.node)
  ### Extract names of OTUs from the GIG clade
  AMFGIGnames<-c(AMFGIG.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFGIGnames[grep("OTU",AMFGIGnames)]
  if (tree==treeFiles[1]) { AMFGIG.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFGIG.OTUs <- append(AMFGIG.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF GIG OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFGIG.node, AMFGIG.clade, AMFGIGnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-GIG OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFGIGotus <- subset(otus,X.OTU.ID %in% AMFGIG.OTUs) # filter out non-GIG OTUs using the list we generated above
# Write to file
write.table(AMFGIGotus,file='./otu97table_clean_AMFGIGonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFGIGseqs <- seqs[AMFGIG.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFGIG.OTUs)]
# Write to file
write.FASTA(AMFGIGseqs,'./otu97repseqs_clean_AMFGIGonly.fasta')

print('AMF-GIG LSU pipeline complete.')
print('AMF-GIG OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Glomeraceae family clade ##############
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define GLO clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MK875635.1_Nanoglomus_plukenetia')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832193_Rhizophagus_dimorphicus')$TipNumber)
  # Determine the node for the GLO clade
  AMFGLO.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF GLO clade and identifying OTUs:')
  # Extract the GLO clade
  AMFGLO.clade <- extract.clade(mytree.rooted.preorder,AMFGLO.node)
  ### Extract names of OTUs from the GLO clade
  AMFGLOnames<-c(AMFGLO.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFGLOnames[grep("OTU",AMFGLOnames)]
  if (tree==treeFiles[1]) { AMFGLO.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFGLO.OTUs <- append(AMFGLO.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF GLO OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFGLO.node, AMFGLO.clade, AMFGLOnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-GLO OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFGLOotus <- subset(otus,X.OTU.ID %in% AMFGLO.OTUs) # filter out non-GLO OTUs using the list we generated above
# Write to file
write.table(AMFGLOotus,file='./otu97table_clean_AMFGLOonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFGLOseqs <- seqs[AMFGLO.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFGLO.OTUs)]
# Write to file
write.FASTA(AMFGLOseqs,'./otu97repseqs_clean_AMFGLOonly.fasta')

print('AMF-GLO LSU pipeline complete.')
print('AMF-GLO OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Pacisporaceae family clade ############
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define PAC clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='FM876831.1_Pacispora_scintillans')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='FM876832.1_Pacispora_scintillans')$TipNumber)
  # Determine the node for the PAC clade
  AMFPAC.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF PAC clade and identifying OTUs:')
  # Extract the PAC clade
  AMFPAC.clade <- extract.clade(mytree.rooted.preorder,AMFPAC.node)
  ### Extract names of OTUs from the PAC clade
  AMFPACnames<-c(AMFPAC.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFPACnames[grep("OTU",AMFPACnames)]
  if (tree==treeFiles[1]) { AMFPAC.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFPAC.OTUs <- append(AMFPAC.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF PAC OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFPAC.node, AMFPAC.clade, AMFPACnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-PAC OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFPACotus <- subset(otus,X.OTU.ID %in% AMFPAC.OTUs) # filter out non-PAC OTUs using the list we generated above
# Write to file
write.table(AMFPACotus,file='./otu97table_clean_AMFPAConly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFPACseqs <- seqs[AMFPAC.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFPAC.OTUs)]
# Write to file
write.FASTA(AMFPACseqs,'./otu97repseqs_clean_AMFPAConly.fasta')

print('AMF-PAC LSU pipeline complete.')
print('AMF-PAC OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Paraglomeraceae family clade ##########
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define PAR clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='KY630236.1_Pervetustus_simplex')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832155_Paraglomus_occultum')$TipNumber)
  # Determine the node for the PAR clade
  AMFPAR.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF PAR clade and identifying OTUs:')
  # Extract the PAR clade
  AMFPAR.clade <- extract.clade(mytree.rooted.preorder,AMFPAR.node)
  ### Extract names of OTUs from the PAR clade
  AMFPARnames<-c(AMFPAR.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFPARnames[grep("OTU",AMFPARnames)]
  if (tree==treeFiles[1]) { AMFPAR.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFPAR.OTUs <- append(AMFPAR.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF PAR OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFPAR.node, AMFPAR.clade, AMFPARnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-PAR OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFPARotus <- subset(otus,X.OTU.ID %in% AMFPAR.OTUs) # filter out non-PAR OTUs using the list we generated above
# Write to file
write.table(AMFPARotus,file='./otu97table_clean_AMFPARonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFPARseqs <- seqs[AMFPAR.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFPAR.OTUs)]
# Write to file
write.FASTA(AMFPARseqs,'./otu97repseqs_clean_AMFPARonly.fasta')

print('AMF-PAR LSU pipeline complete.')
print('AMF-PAR OTU table and sequences have been written to file.')

#################################################################################
#### Determine which OTUs fall within the Sacculosporaceae family clade #########
#################################################################################
print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*') # make list of the tree files output from RAxML
for (tree in treeFiles) {
  print(paste0('Processing tree from OTU subset: ',tree))
  mytree <- read.tree(paste0('RAxMLfiles/',tree)) # load one tree
  # reformat tree to Preorder
  mytree.rooted <- RootTree(mytree,"M11585.1_Oryza_sativa")
  mytree.rooted.preorder <- Preorder(mytree.rooted)
  # Make data frame with tip names and numbers
  tipDF <- data.frame('TipLabel'=TipLabels(mytree.rooted.preorder))
  tipDF$TipNumber <- rownames(tipDF) # copy row names (integer IDs) to column
  # get integer values of tips that define SAC clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='KX355820.1_Sacculospora_baltica')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='KX345939.1_Sacculospora_felinovii')$TipNumber)
  # Determine the node for the SAC clade
  AMFSAC.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF SAC clade and identifying OTUs:')
  # Extract the SAC clade
  AMFSAC.clade <- extract.clade(mytree.rooted.preorder,AMFSAC.node)
  ### Extract names of OTUs from the SAC clade
  AMFSACnames<-c(AMFSAC.clade$tip.label)
  ### Subset only those that start with 'OTU'
  hits<-AMFSACnames[grep("OTU",AMFSACnames)]
  if (tree==treeFiles[1]) { AMFSAC.OTUs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFSAC.OTUs <- append(AMFSAC.OTUs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF SAC OTU names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFSAC.node, AMFSAC.clade, AMFSACnames, hits) # clean up except for the AMF.OTUs vector
}

####### Load qiime2 output and prune out non-SAC OTUs #######
print('Loading qiime2 output to be pruned - OTU table')
# Load OTU table and prune
otus <- read.delim('./otu97table_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFSACotus <- subset(otus,X.OTU.ID %in% AMFSAC.OTUs) # filter out non-SAC OTUs using the list we generated above
# Write to file
write.table(AMFSACotus,file='./otu97table_clean_AMFSAConly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - OTU rep. seqs')
# Load OTU rep seqs and prune
seqs <- read.FASTA('./otu97repseqs_clean_BLAST.fasta', type = 'DNA')
AMFSACseqs <- seqs[AMFSAC.OTUs]
otherseqs <- seqs[setdiff(names(seqs),AMFSAC.OTUs)]
# Write to file
write.FASTA(AMFSACseqs,'./otu97repseqs_clean_AMFSAConly.fasta')

print('AMF-SAC LSU pipeline complete.')
print('AMF-SAC OTU table and sequences have been written to file.')