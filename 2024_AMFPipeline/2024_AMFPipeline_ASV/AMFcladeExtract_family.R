####### Load required packages #######

library(ape)
library(TreeTools)

#################################################################################
####### Determine which ASVs fall within the Acaulosporaceae family clade #######
#################################################################################
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
  # get integer values of tips that define ACA clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='FN547517_Acaulospora_laevis')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832207_Acaulospora_tuberculata')$TipNumber)
  # Determine the node for the ACA clade
  AMFACA.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF ACA clade and identifying ASVs:')
  # Extract the ACA clade
  AMFACA.clade <- extract.clade(mytree.rooted.preorder,AMFACA.node)
  ### Extract names of ASVs from the ACA clade
  AMFACAnames<-c(AMFACA.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFACAnames[grep("ASV",AMFACAnames)]
  if (tree==treeFiles[1]) { AMFACA.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFACA.ASVs <- append(AMFACA.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF ACA ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFACA.node, AMFACA.clade, AMFACAnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-ACA ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFACAASVs <- subset(ASVs,X.ASV.ID %in% AMFACA.ASVs) # filter out non-ACA ASVs using the list we generated above
# Write to file
write.table(AMFACAASVs,file='./ASVtable_clean_AMFACAonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFACAseqs <- seqs[AMFACA.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFACA.ASVs)]
# Write to file
write.FASTA(AMFACAseqs,'./ASVrepseqs_clean_AMFACAonly.fasta')

print('AMF-ACA LSU pipeline complete.')
print('AMF-ACA ASV table and sequences have been written to file.')


#################################################################################
#### Determine which ASVs fall within the Ambisporaceae family clade ##########
#################################################################################
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
  # get integer values of tips that define AMB clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MT832166_Ambispora_leptoticha')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='JF439210_Ambispora_gerdemannii')$TipNumber)
  # Determine the node for the AMB clade
  AMFAMB.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF AMB clade and identifying ASVs:')
  # Extract the AMB clade
  AMFAMB.clade <- extract.clade(mytree.rooted.preorder,AMFAMB.node)
  ### Extract names of ASVs from the AMB clade
  AMFAMBnames<-c(AMFAMB.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFAMBnames[grep("ASV",AMFAMBnames)]
  if (tree==treeFiles[1]) { AMFAMB.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFAMB.ASVs <- append(AMFAMB.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF AMB ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFAMB.node, AMFAMB.clade, AMFAMBnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-AMB ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFAMBASVs <- subset(ASVs,X.ASV.ID %in% AMFAMB.ASVs) # filter out non-AMB ASVs using the list we generated above
# Write to file
write.table(AMFAMBASVs,file='./ASVtable_clean_AMFAMBonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFAMBseqs <- seqs[AMFAMB.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFAMB.ASVs)]
# Write to file
write.FASTA(AMFAMBseqs,'./ASVrepseqs_clean_AMFAMBonly.fasta')

print('AMF-AMB LSU pipeline complete.')
print('AMF-AMB ASV table and sequences have been written to file.')

#################################################################################
#### Determine which ASVs fall within the Archaeosporaceae family clade #########
#################################################################################
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
  # get integer values of tips that define ARC clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MT832160_Archaeospora_trappei')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832165_Archaeospora_schencki')$TipNumber)
  # Determine the node for the ARC clade
  AMFARC.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF ARC clade and identifying ASVs:')
  # Extract the ARC clade
  AMFARC.clade <- extract.clade(mytree.rooted.preorder,AMFARC.node)
  ### Extract names of ASVs from the ARC clade
  AMFARCnames<-c(AMFARC.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFARCnames[grep("ASV",AMFARCnames)]
  if (tree==treeFiles[1]) { AMFARC.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFARC.ASVs <- append(AMFARC.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF ARC ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFARC.node, AMFARC.clade, AMFARCnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-ARC ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFARCASVs <- subset(ASVs,X.ASV.ID %in% AMFARC.ASVs) # filter out non-ARC ASVs using the list we generated above
# Write to file
write.table(AMFARCASVs,file='./ASVtable_clean_AMFARConly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFARCseqs <- seqs[AMFARC.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFARC.ASVs)]
# Write to file
write.FASTA(AMFARCseqs,'./ASVrepseqs_clean_AMFARConly.fasta')

print('AMF-ARC LSU pipeline complete.')
print('AMF-ARC ASV table and sequences have been written to file.')

#################################################################################
#### Determine which ASVs fall within the Claroideoglomeraceae family clade #####
#################################################################################
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
  # get integer values of tips that define CLA clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='AJ972464.1_Albahypha_drummondii')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832172_Claroideoglomus_etunicatum')$TipNumber)
  # Determine the node for the CLA clade
  AMFCLA.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF CLA clade and identifying ASVs:')
  # Extract the CLA clade
  AMFCLA.clade <- extract.clade(mytree.rooted.preorder,AMFCLA.node)
  ### Extract names of ASVs from the CLA clade
  AMFCLAnames<-c(AMFCLA.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFCLAnames[grep("ASV",AMFCLAnames)]
  if (tree==treeFiles[1]) { AMFCLA.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFCLA.ASVs <- append(AMFCLA.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF CLA ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFCLA.node, AMFCLA.clade, AMFCLAnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-CLA ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFCLAASVs <- subset(ASVs,X.ASV.ID %in% AMFCLA.ASVs) # filter out non-CLA ASVs using the list we generated above
# Write to file
write.table(AMFCLAASVs,file='./ASVtable_clean_AMFCLAonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFCLAseqs <- seqs[AMFCLA.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFCLA.ASVs)]
# Write to file
write.FASTA(AMFCLAseqs,'./ASVrepseqs_clean_AMFCLAonly.fasta')

print('AMF-CLA LSU pipeline complete.')
print('AMF-CLA ASV table and sequences have been written to file.')

#################################################################################
#### Determine which ASVs fall within the Diversisporaceae family clade #########
#################################################################################
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
  # get integer values of tips that define DIV clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MG459208.1_Desertispora_omaniana')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832219_Diversispora_trimurales')$TipNumber)
  # Determine the node for the DIV clade
  AMFDIV.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF DIV clade and identifying ASVs:')
  # Extract the DIV clade
  AMFDIV.clade <- extract.clade(mytree.rooted.preorder,AMFDIV.node)
  ### Extract names of ASVs from the DIV clade
  AMFDIVnames<-c(AMFDIV.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFDIVnames[grep("ASV",AMFDIVnames)]
  if (tree==treeFiles[1]) { AMFDIV.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFDIV.ASVs <- append(AMFDIV.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF DIV ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFDIV.node, AMFDIV.clade, AMFDIVnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-DIV ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFDIVASVs <- subset(ASVs,X.ASV.ID %in% AMFDIV.ASVs) # filter out non-DIV ASVs using the list we generated above
# Write to file
write.table(AMFDIVASVs,file='./ASVtable_clean_AMFDIVonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFDIVseqs <- seqs[AMFDIV.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFDIV.ASVs)]
# Write to file
write.FASTA(AMFDIVseqs,'./ASVrepseqs_clean_AMFDIVonly.fasta')

print('AMF-DIV LSU pipeline complete.')
print('AMF-DIV ASV table and sequences have been written to file.')

#################################################################################
#### Determine which ASVs fall within the Gigasporaceae family clade ############
#################################################################################
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
  # get integer values of tips that define GIG clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MT832220_Scutellospora_calospora')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832238_Gigaspora_gigantea')$TipNumber)
  # Determine the node for the GIG clade
  AMFGIG.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF GIG clade and identifying ASVs:')
  # Extract the GIG clade
  AMFGIG.clade <- extract.clade(mytree.rooted.preorder,AMFGIG.node)
  ### Extract names of ASVs from the GIG clade
  AMFGIGnames<-c(AMFGIG.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFGIGnames[grep("ASV",AMFGIGnames)]
  if (tree==treeFiles[1]) { AMFGIG.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFGIG.ASVs <- append(AMFGIG.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF GIG ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFGIG.node, AMFGIG.clade, AMFGIGnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-GIG ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_AMFonly.tsv',sep='\t',header=TRUE)
AMFGIGASVs <- subset(ASVs,X.ASV.ID %in% AMFGIG.ASVs) # filter out non-GIG ASVs using the list we generated above
# Write to file
write.table(AMFGIGASVs,file='./ASVtable_clean_AMFGIGonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFGIGseqs <- seqs[AMFGIG.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFGIG.ASVs)]
# Write to file
write.FASTA(AMFGIGseqs,'./ASVrepseqs_clean_AMFGIGonly.fasta')

print('AMF-GIG LSU pipeline complete.')
print('AMF-GIG ASV table and sequences have been written to file.')

#################################################################################
#### Determine which ASVs fall within the Glomeraceae family clade ##############
#################################################################################
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
  # get integer values of tips that define GLO clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='MK875635.1_Nanoglomus_plukenetia')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832193_Rhizophagus_dimorphicus')$TipNumber)
  # Determine the node for the GLO clade
  AMFGLO.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF GLO clade and identifying ASVs:')
  # Extract the GLO clade
  AMFGLO.clade <- extract.clade(mytree.rooted.preorder,AMFGLO.node)
  ### Extract names of ASVs from the GLO clade
  AMFGLOnames<-c(AMFGLO.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFGLOnames[grep("ASV",AMFGLOnames)]
  if (tree==treeFiles[1]) { AMFGLO.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFGLO.ASVs <- append(AMFGLO.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF GLO ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFGLO.node, AMFGLO.clade, AMFGLOnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-GLO ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFGLOASVs <- subset(ASVs,X.ASV.ID %in% AMFGLO.ASVs) # filter out non-GLO ASVs using the list we generated above
# Write to file
write.table(AMFGLOASVs,file='./ASVtable_clean_AMFGLOonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFGLOseqs <- seqs[AMFGLO.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFGLO.ASVs)]
# Write to file
write.FASTA(AMFGLOseqs,'./ASVrepseqs_clean_AMFGLOonly.fasta')

print('AMF-GLO LSU pipeline complete.')
print('AMF-GLO ASV table and sequences have been written to file.')

#################################################################################
#### Determine which ASVs fall within the Pacisporaceae family clade ############
#################################################################################
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
  # get integer values of tips that define PAC clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='FM876831.1_Pacispora_scintillans')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='FM876832.1_Pacispora_scintillans')$TipNumber)
  # Determine the node for the PAC clade
  AMFPAC.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF PAC clade and identifying ASVs:')
  # Extract the PAC clade
  AMFPAC.clade <- extract.clade(mytree.rooted.preorder,AMFPAC.node)
  ### Extract names of ASVs from the PAC clade
  AMFPACnames<-c(AMFPAC.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFPACnames[grep("ASV",AMFPACnames)]
  if (tree==treeFiles[1]) { AMFPAC.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFPAC.ASVs <- append(AMFPAC.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF PAC ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFPAC.node, AMFPAC.clade, AMFPACnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-PAC ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_AMFonly.tsv',sep='\t',header=TRUE)
AMFPACASVs <- subset(ASVs,X.ASV.ID %in% AMFPAC.ASVs) # filter out non-PAC ASVs using the list we generated above
# Write to file
write.table(AMFPACASVs,file='./ASVtable_clean_AMFPAConly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFPACseqs <- seqs[AMFPAC.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFPAC.ASVs)]
# Write to file
write.FASTA(AMFPACseqs,'./ASVrepseqs_clean_AMFPAConly.fasta')

print('AMF-PAC LSU pipeline complete.')
print('AMF-PAC ASV table and sequences have been written to file.')

#################################################################################
#### Determine which ASVs fall within the Paraglomeraceae family clade ##########
#################################################################################
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
  # get integer values of tips that define PAR clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='KY630232.1_Innospora_majewskii')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='MT832155_Paraglomus_occultum')$TipNumber)
  # Determine the node for the PAR clade
  AMFPAR.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF PAR clade and identifying ASVs:')
  # Extract the PAR clade
  AMFPAR.clade <- extract.clade(mytree.rooted.preorder,AMFPAR.node)
  ### Extract names of ASVs from the PAR clade
  AMFPARnames<-c(AMFPAR.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFPARnames[grep("ASV",AMFPARnames)]
  if (tree==treeFiles[1]) { AMFPAR.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFPAR.ASVs <- append(AMFPAR.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF PAR ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFPAR.node, AMFPAR.clade, AMFPARnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-PAR ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFPARASVs <- subset(ASVs,X.ASV.ID %in% AMFPAR.ASVs) # filter out non-PAR ASVs using the list we generated above
# Write to file
write.table(AMFPARASVs,file='./ASVtable_clean_AMFPARonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFPARseqs <- seqs[AMFPAR.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFPAR.ASVs)]
# Write to file
write.FASTA(AMFPARseqs,'./ASVrepseqs_clean_AMFPARonly.fasta')

print('AMF-PAR LSU pipeline complete.')
print('AMF-PAR ASV table and sequences have been written to file.')

#################################################################################
#### Determine which ASVs fall within the Pervestustaceae family clade ##########
#################################################################################
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
  # get integer values of tips that define PER clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='KY630236.1_Pervetustus_simplex')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='KY630244.1_Pervetustus_simplex')$TipNumber)
  # Determine the node for the PER clade
  AMFPER.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF PER clade and identifying ASVs:')
  # Extract the PER clade
  AMFPER.clade <- extract.clade(mytree.rooted.preorder,AMFPER.node)
  ### Extract names of ASVs from the PER clade
  AMFPERnames<-c(AMFPER.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFPERnames[grep("ASV",AMFPERnames)]
  if (tree==treeFiles[1]) { AMFPER.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFPER.ASVs <- append(AMFPER.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF PER ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFPER.node, AMFPER.clade, AMFPERnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-PER ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_AMFonly.tsv',sep='\t',header=TRUE)
AMFPERASVs <- subset(ASVs,X.ASV.ID %in% AMFPER.ASVs) # filter out non-PER ASVs using the list we generated above
# Write to file
write.table(AMFPERASVs,file='./ASVtable_clean_AMFPERonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFPERseqs <- seqs[AMFPER.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFPER.ASVs)]
# Write to file
write.FASTA(AMFPERseqs,'./ASVrepseqs_clean_AMFPERonly.fasta')

print('AMF-PER LSU pipeline complete.')
print('AMF-PER ASV table and sequences have been written to file.')


#################################################################################
#### Determine which ASVs fall within the Sacculosporaceae family clade #########
#################################################################################
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
  # get integer values of tips that define SAC clade edges
  cladeStart <- as.integer(subset(tipDF,TipLabel=='KX355820.1_Sacculospora_baltica')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='KX345939.1_Sacculospora_felinovii')$TipNumber)
  # Determine the node for the SAC clade
  AMFSAC.node <- MRCA(cladeStart,cladeStop,AllAncestors(parent=mytree.rooted.preorder$edge[,1],child=mytree.rooted.preorder$edge[,2]))
  print('---Extracting the AMF SAC clade and identifying ASVs:')
  # Extract the SAC clade
  AMFSAC.clade <- extract.clade(mytree.rooted.preorder,AMFSAC.node)
  ### Extract names of ASVs from the SAC clade
  AMFSACnames<-c(AMFSAC.clade$tip.label)
  ### Subset only those that start with 'ASV'
  hits<-AMFSACnames[grep("ASV",AMFSACnames)]
  if (tree==treeFiles[1]) { AMFSAC.ASVs <- hits } # if this is the first tree in the list, make a new vector to store the hits
  else {AMFSAC.ASVs <- append(AMFSAC.ASVs, hits)} # otherwise, add the hits from this tree to the growing vector of AMF SAC ASV names
  rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMFSAC.node, AMFSAC.clade, AMFSACnames, hits) # clean up except for the AMF.ASVs vector
}

####### Load qiime2 output and prune out non-SAC ASVs #######
print('Loading qiime2 output to be pruned - ASV table')
# Load ASV table and prune
ASVs <- read.delim('./ASVtable_clean_BLASTonly.tsv',sep='\t',header=TRUE)
AMFSACASVs <- subset(ASVs,X.ASV.ID %in% AMFSAC.ASVs) # filter out non-SAC ASVs using the list we generated above
# Write to file
write.table(AMFSACASVs,file='./ASVtable_clean_AMFSAConly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFSACseqs <- seqs[AMFSAC.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMFSAC.ASVs)]
# Write to file
write.FASTA(AMFSACseqs,'./ASVrepseqs_clean_AMFSAConly.fasta')

print('AMF-SAC LSU pipeline complete.')
print('AMF-SAC ASV table and sequences have been written to file.')