####### Load required packages #######

library(ape)
library(TreeTools)

##############################################################
####### Determine which ASVs fall within the AMF clade #######
##############################################################

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
  cladeStart <- as.integer(subset(tipDF,TipLabel=='FM876840_Geosiphon_pyriformis')$TipNumber)
  cladeStop <- as.integer(subset(tipDF,TipLabel=='FR750083.1_Paraglomeraceae_Paraglomus_laccatum')$TipNumber)
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
write.table(AMFASVs,file='./cladeOutputs/ASVtable_clean_AMFonly.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
write.table(otherASVs,file='./cladeOutputs/ASVtable_clean_nonAMF.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)

print('Loading qiime2 output to be pruned - ASV rep. seqs')
# Load ASV rep seqs and prune:
seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
AMFseqs <- seqs[AMF.ASVs]
otherseqs <- seqs[setdiff(names(seqs),AMF.ASVs)]
# Write to file:
write.FASTA(AMFseqs,'./cladeOutputs/ASVrepseqs_clean_AMFonly.fasta')
write.FASTA(otherseqs,'./cladeOutputs/ASVrepseqs_clean_nonAMF.fasta')

print('AMF - LSU pipeline complete.')
print('AMF-only ASV table and sequences have been written to file.')
print('non-AMF ASV table and sequences have also been written to file.')

######################################################################
####### Determine which ASVs fall within the AMF family clades #######
######################################################################

print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*')

cladeDefs <- list(
  POL = c("MZ359659.1_Polonosporaceae_Polonospora_polonica", "MZ359658.1_Polonosporaceae_Polonospora_polonica"),
  AMB = c("MT832166.1_Ambisporaceae_Ambispora_leptoticha", "FN547539.1_Ambisporaceae_Ambispora_fennica"),
  ARC = c("MT832161.1_Archaeosporaceae_Archaeospora_trappei", "HG977202.1_Archaeosporaceae_Archaeospora_spainiae"),
  ENT = c("AJ972466.1_Entrophosporaceae_Entrophospora_Albahypha_drummondi", "MT722022.1_Entrophosporaceae_Entrophospora_argentinensis"),
  GLO = c("MK036783.1_Glomeraceae_Sclerocarpum_amazonicum", "PP753776.1_Glomeraceae_Glomus_mongioiense"),
  SEP = c("FN547496.1_Septoglomeraceae_Funneliformis_caledonium", "HF674439.1_Septoglomeraceae_Septoglomus_altomontanum"),
  DOM = c("MK875637.1_Dominikiaceae_Nanoglomus_plukenetiae", "KJ564165.2_Dominikiaceae_Dominikia_minuta"),
  SCL_KAM = c("MT832183.1_Sclerocystaceae_Oehlia_diaphana", "MN130956.1_Sclerocystaceae_Rhizophagus_maiae"),
  GIG = c("HQ871519.1_Gigasporaceae_Scutellospora_pernambucana", "OP205496.1_Gigasporaceae_Gigaspora_polymorphira"),
  PAC = c("FM876832.1_Pacisporaceae_Pacispora_scintillans", "FM876831.1_Pacisporaceae_Pacispora_scintillans"),
  SAC = c("KX345943.1_Sclerocystaceae_Sacculospora_felinovii", "KX355820.1_Sclerocystaceae_Sacculospora_baltica"),
  DIV = c("KF154769.1_Diversisporaceae_Desertispora_omaniana", "AY639236.1_Diversisporaceae_Diversispora_celata"),
  ACA = c("MT832197.1_Acaulosporaceae_Acaulospora_colombiana", "OK356207.1_Acaulosporaceae_Acaulospora_excavata"),
  PER = c("KY630236.1_Pervetustaceae_Pervetustus_simplex", "MT832159.1_Pervetustaceae_Pervetustus_simplex"),
  PAR = c("KY630232.1_Paraglomeraceae_Innospora_majewskii", "FR750083.1_Paraglomeraceae_Paraglomus_laccatum"),
)

for (cladeName in names(cladeDefs)) {
  print(paste('Processing clade:', cladeName))
  cladeTips <- cladeDefs[[cladeName]]
  allHits <- c()
  
  for (tree in treeFiles) {
    print(paste0('Processing tree from ASV subset: ', tree))
    mytree <- read.tree(paste0('RAxMLfiles/', tree))
    mytree.rooted <- RootTree(mytree, "M11585.1_Oryza_sativa")
    mytree.rooted.preorder <- Preorder(mytree.rooted)
    
    tipDF <- data.frame('TipLabel' = TipLabels(mytree.rooted.preorder))
    tipDF$TipNumber <- rownames(tipDF)
    
    cladeStart <- as.integer(subset(tipDF, TipLabel == cladeTips[1])$TipNumber)
    cladeStop <- as.integer(subset(tipDF, TipLabel == cladeTips[2])$TipNumber)
    
    AMF.node <- MRCA(cladeStart, cladeStop, AllAncestors(parent=mytree.rooted.preorder$edge[,1], child=mytree.rooted.preorder$edge[,2]))
    AMF.clade <- extract.clade(mytree.rooted.preorder, AMF.node)
    
    tipNames <- AMF.clade$tip.label
    hits <- tipNames[grep("ASV", tipNames)]
    allHits <- append(allHits, hits)
    
    rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMF.node, AMF.clade, tipNames, hits)
  }
  
  # Prune ASV table
  print(paste('Pruning ASV table for', cladeName))
  ASVs <- read.delim('./cladeOutputs/ASVtable_clean_AMFonly.tsv', sep='\t', header=TRUE)
  cladeASVs <- subset(ASVs, X.ASV.ID %in% allHits)
  write.table(cladeASVs, file=paste0('./cladeOutputs/ASVtable_clean_AMF', cladeName, 'only.tsv'), sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  # Prune FASTA
  print(paste('Pruning FASTA for', cladeName))
  seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
  cladeSeqs <- seqs[allHits]
  write.FASTA(cladeSeqs, paste0('./cladeOutputs/ASVrepseqs_clean_AMF', cladeName, 'only.fasta'))
}

print('Multi-clade AMF family pipeline complete.')

#####################################################################
####### Determine which ASVs fall within the AMF genus clades #######
#####################################################################

print('Loading RAxML best trees...')
treeFiles <- dir('./RAxMLfiles/', pattern ='RAxML_bestTree*')

cladeDefs <- list(
  POL_POL = c("MZ359659.1_Polonosporaceae_Polonospora_polonica", "MZ359658.1_Polonosporaceae_Polonospora_polonica"),
  AMB_AMB = c("MT832166.1_Ambisporaceae_Ambispora_leptoticha", "FN547539.1_Ambisporaceae_Ambispora_fennica"),
  ARC_ARC = c("MT832161.1_Archaeosporaceae_Archaeospora_trappei", "HG977202.1_Archaeosporaceae_Archaeospora_spainiae"),
  ENT_ENT = c("AJ972466.1_Entrophosporaceae_Entrophospora_Albahypha_drummondi", "MT722022.1_Entrophosporaceae_Entrophospora_argentinensis"),
  GLO_SCL = c("MK036783.1_Glomeraceae_Sclerocarpum_amazonicum", "MK036781.1_Glomeraceae_Sclerocarpum_amazonicum"),
  GLO_COM = c("OQ437305.1_Glomeraceae_Complexispora_mediterranea", "OQ437299.1_Glomeraceae_Complexispora_multistratosa"),
  GLO_GLO = c("MH560608.1_Glomeraceae_Glomus_bareae", "PP753776.1_Glomeraceae_Glomus_mongioiense"),
  SEP_SEP = c("JQ312667.1_Septoglomeraceae_Septoglomus_titan", "HF674439.1_Septoglomeraceae_Septoglomus_altomontanum"),
  SEP_FUNG = c("MK348930.1_Septoglomeraceae_Funneliglomus_sanmartinense", "MK348935.1_Septoglomeraceae_Funneliglomus_sanmartinense"),
  SEP_FUNF = c("FN547496.1_Septoglomeraceae_Funneliformis_caledonium", "FJ461828.1_Septoglomeraceae_Funneliformis_coronatus"),
  DOM_NAN = c("MK875637.1_Dominikiaceae_Nanoglomus_plukenetiae", "MK875635.1_Dominikiaceae_Nanoglomus_plukenetiae"),
  DOM_ORI = c("KY555051.1_Dominikiaceae_Orientoglomus_emiratium", "KY555052.1_Dominikiaceae_Orientoglomus_emiratium"),
  DOM_MIC = c("MG710519.1_Dominikiaceae_Microdominikia_litorea", "MG710518.1_Dominikiaceae_Microdominikia_litorea"),
  DOM_MAC = c("HG798897.1_Dominikiaceae_Macrodominikia_compressa", "HG798898.1_Dominikiaceae_Macrodominikia_compressa"),
  DOM_DOM = c("KM056664.1_Dominikiaceae_Dominikia_aurea", "KJ564165.2_Dominikiaceae_Dominikia_minuta"),
  SCL_KAM_OEH = c("MT832183.1_Sclerocystaceae_Oehlia_diaphana", "MT832184.1_Sclerocystaceae_Oehlia_diaphana"),
  SCL_KAM_KAM = c("KJ564133.1_Kamienskiaceae_Kamienskia_bistrata", "KJ564137.1_Kamienskiaceae_Kamienskia_bistrata"),
  SCL_KAM_EPI = c("OQ437315.1_Kamienskiaceae_Epigeocarpum_japonicum", "MW507156.1_Kamienskiaceae_Epigeocarpum_crypticum"),
  SCL_KAM_HAL = c("MH560602.1_Sclerocystaceae_Halonatospora_pansihalos", "MH560604.1_Sclerocystaceae_Halonatospora_pansihalos"),
  SCL_KAM_MIC = c("KX758123.1_Kamienskiaceae_Microkamienskia_divaricata", "KJ564144.1_Kamienskiaceae_Microkamienskia_perpusilla"),
  SCL_KAM_SCL = c("MT832185.1_Sclerocystaceae_Sclerocystis_sinuosa", "MT832186.1_Sclerocystaceae_Sclerocystis_sinuosa"),
  SCL_KAM_SIL = c("KY362438.1_Sclerocystaceae_Silvaspora_neocaledonica", "KY362436.1_Sclerocystaceae_Silvaspora_neocaledonica"),
  SCL_KAM_RHI = c("HG964396.1_Sclerocystaceae_Rhizophagus_melanus", "MN130956.1_Sclerocystaceae_Rhizophagus_maiae"),
  GIG_SCU = c("HQ871519.1_Gigasporaceae_Scutellospora_pernambucana", "OR669037.1_Gigasporaceae_Scutellospora_intraundulata"),
  GIG_BUL = c("KJ944325.1_Gigasporaceae_Bulbospora_minima", "KJ944323.1_Gigasporaceae_Bulbospora_minima"),
  GIG_CET = c("MT832222.1_Gigasporaceae_Cetraspora_pellucida", "OP004059.1_Gigasporaceae_Cetraspora_huaxica"),
  GIG_RAC = c("FR750135.1_Gigasporaceae_Racocetra_weresubiae", "MT967288.1_Gigasporaceae_Racocetra_gregaria"),
  GIG_DEN = c("MT832225.1_Gigasporaceae_Dentiscutata_nigra", "JN971066.1_Gigasporaceae_Dentiscutata_Fuscutata_aurea"),
  GIG_INT = c("JN971072.1_Gigasporaceae_Intraornatospora_intraornata", "JN971073.1_Gigasporaceae_Intraornatospora_intraornata"),
  GIG_PAR = c("JN971081.1_Gigasporaceae_Paradentiscutata_maritima", "JN971070.1_Gigasporaceae_Paradentiscutata_baiana"),
  GIG_GIG = c("OQ680684.1_Gigasporaceae_Gigaspora_siqueirae", "OP205496.1_Gigasporaceae_Gigaspora_polymorphira"),
  PAC_PAC = c("FM876832.1_Pacisporaceae_Pacispora_scintillans", "FM876831.1_Pacisporaceae_Pacispora_scintillans"),
  SAC_SAC = c("KX345943.1_Sclerocystaceae_Sacculospora_felinovii", "KX355820.1_Sclerocystaceae_Sacculospora_baltica"),
  DIV_DES = c("KF154769.1_Diversisporaceae_Desertispora_omaniana", "MG459208.1_Diversisporaceae_Desertispora_omaniana"),
  DIV_COR = c("FJ461836.1_Diversisporaceae_Corymbiglomus_globiferum", "KF060295.1_Diversisporaceae_Corymbiglomus_corymbiforme"),
  DIV_RED = c("HG518628.1_Diversisporaceae_Redeckera_megalocarpum", "MT832215.1_Diversisporaceae_Redeckera_megalocarpum"),
  DIV_DIV = c("KJ850201.1_Diversisporaceae_Diversispora_gibbosa", "AY639236.1_Diversisporaceae_Diversispora_celata"),
  ACA_ACA = c("MT832197.1_Acaulosporaceae_Acaulospora_colombiana", "OK356207.1_Acaulosporaceae_Acaulospora_excavata"),
  PER_PER = c("KY630236.1_Pervetustaceae_Pervetustus_simplex", "MT832159.1_Pervetustaceae_Pervetustus_simplex"),
  PAR_INN = c("KY630232.1_Paraglomeraceae_Innospora_majewskii", "JN131593.1_Paraglomeraceae_Innospora_majewskii"),
  PAR_PAR = c("FR750049.1_Paraglomeraceae_Paraglomus_brasilianum", "FR750083.1_Paraglomeraceae_Paraglomus_laccatum")
)

for (cladeName in names(cladeDefs)) {
  print(paste('Processing clade:', cladeName))
  cladeTips <- cladeDefs[[cladeName]]
  allHits <- c()
  
  for (tree in treeFiles) {
    print(paste0('Processing tree from ASV subset: ', tree))
    mytree <- read.tree(paste0('RAxMLfiles/', tree))
    mytree.rooted <- RootTree(mytree, "M11585.1_Oryza_sativa")
    mytree.rooted.preorder <- Preorder(mytree.rooted)
    
    tipDF <- data.frame('TipLabel' = TipLabels(mytree.rooted.preorder))
    tipDF$TipNumber <- rownames(tipDF)
    
    cladeStart <- as.integer(subset(tipDF, TipLabel == cladeTips[1])$TipNumber)
    cladeStop <- as.integer(subset(tipDF, TipLabel == cladeTips[2])$TipNumber)
    
    AMF.node <- MRCA(cladeStart, cladeStop, AllAncestors(parent=mytree.rooted.preorder$edge[,1], child=mytree.rooted.preorder$edge[,2]))
    AMF.clade <- extract.clade(mytree.rooted.preorder, AMF.node)
    
    tipNames <- AMF.clade$tip.label
    hits <- tipNames[grep("ASV", tipNames)]
    allHits <- append(allHits, hits)
    
    rm(mytree, mytree.rooted, mytree.rooted.preorder, tipDF, cladeStart, cladeStop, AMF.node, AMF.clade, tipNames, hits)
  }
  
  # Prune ASV table
  print(paste('Pruning ASV table for', cladeName))
  ASVs <- read.delim('./cladeOutputs/ASVtable_clean_AMFonly.tsv', sep='\t', header=TRUE)
  cladeASVs <- subset(ASVs, X.ASV.ID %in% allHits)
  write.table(cladeASVs, file=paste0('./cladeOutputs/ASVtable_clean_AMF', cladeName, 'only.tsv'), sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  # Prune FASTA
  print(paste('Pruning FASTA for', cladeName))
  seqs <- read.FASTA('./ASVrepseqs_clean_BLAST.fasta', type = 'DNA')
  cladeSeqs <- seqs[allHits]
  write.FASTA(cladeSeqs, paste0('./cladeOutputs/ASVrepseqs_clean_AMF', cladeName, 'only.fasta'))
}

print('Multi-clade AMF genus pipeline complete.')
