####### Load required packages #######

library(tidyverse)

###### Bring in Glomeromycota (AMF) table ####

AMF <- read.table(file='./ASVtable_clean_AMFonly.tsv', header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(phylum = "Glomeromycota")

###### Bring in Family tables #####

POL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPOLonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Polonosporaceae")
AMB <- read.table(file='./cladeOutputs/ASVtable_clean_AMFAMBonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Ambisporaceae")
ARC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFARConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Archaeosporaceae")
ENT <- read.table(file='./cladeOutputs/ASVtable_clean_AMFENTonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Entrophosporaceae")
GLO <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGLOonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Glomeraceae")
SEP <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSEPonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Septoglomeraceae")
DOM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Dominikiaceae")
SCL_KAM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Sclerocystaceae_Kamienskiaceae")
GIG <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIGonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Gigasporaceae")
PAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Pacisporaceae")
SAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Sacculosporaceae")
DIV <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIVonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Diversisporaceae")
ACA <- read.table(file='./cladeOutputs/ASVtable_clean_AMFACAonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Acaulosporaceae")
PER <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPERonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Pervetustaceae")
PAR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPARonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(family = "Paraglomeraceae")

###### Bring in Genera tables #####

POL_POL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPOL_POLonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Polonospora")
AMB_AMB <- read.table(file='./cladeOutputs/ASVtable_clean_AMFAMB_AMBonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Ambispora")
ARC_ARC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFARC_ARConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Archaeospora")
ENT_ENT <- read.table(file='./cladeOutputs/ASVtable_clean_AMFENT_ENTonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Entrophospora")
GLO_SCL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGLO_SCLonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Sclerocarpum")
GLO_COM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGLO_COMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Complexispora")
GLO_GLO <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGLO_GLOonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Glomus")
SEP_SEP <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSEP_SEPonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Septoglomus")
SEP_FUNG <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSEP_FUNGonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Funneliformis")
SEP_FUNF <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSEP_FUNFonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Funneliformis")
DOM_NAN <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_NANonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Nanoglomus")
DOM_ORI <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_ORIonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Orientoglomus")
DOM_MIC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_MIConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Microdominikia")
DOM_MAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_MAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Macrodominikia")
DOM_DOM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_DOMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Dominikia")
SCL_KAM_OEH <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_OEHonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Kamienskia")
SCL_KAM_KAM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_KAMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Kamienskia")
SCL_KAM_EPI <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_EPIonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Kamienskia")
SCL_KAM_HAL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_HALonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Kamienskia")
SCL_KAM_MIC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_SCLonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Kamienskia")
SCL_KAM_SCL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_KAMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Kamienskia")
SCL_KAM_SIL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_SILonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Kamienskia")
SCL_KAM_RHI <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_RHIonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Kamienskia")
GIG_SCU <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_SCUonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Scutellospora")
GIG_BUL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_BULonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Bulbospora")
GIG_CET <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_CETonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Cetraspora") 
GIG_RAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_RAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Racocetra") 
GIG_DEN <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_DENonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Dentiscutata")
GIG_INT <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_INTonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Intraornatospora")
GIG_PAR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_PARonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Paradentiscutata")
GIG_GIG <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_GIGonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Gigaspora")
PAC_PAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPAC_PAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Pacispora") 
SAC_SAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSAC_SAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Sacculospora")  
DIV_DES <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_DESonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Desertispora") 
DIV_COR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_CORonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Redeckera") 
DIV_RED <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_REDonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Desertispora") 
DIV_DIV <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_DIVonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Redeckera") 
ACA_ACA <- read.table(file='./cladeOutputs/ASVtable_clean_AMFACA_ACAonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Acaulospora") 
PER_PER <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPER_PERonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Pervetustus") 
PAR_INN <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPAR_INNonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Innospora") 
PAR_PAR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPAR_PARonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(genus = "Paraglomus") 

###### Assign taxonomy to AMF #####

AMF.tax <- AMF %>%
  #join with family level taxonomy
  left_join(POL, by = "X.ASV.ID") %>%
  left_join(AMB, by = "X.ASV.ID") %>%
  left_join(ARC, by = "X.ASV.ID") %>%
  left_join(ENT, by = "X.ASV.ID") %>%
  left_join(GLO, by = "X.ASV.ID") %>%
  left_join(SEP, by = "X.ASV.ID") %>%
  left_join(DOM, by = "X.ASV.ID") %>%
  left_join(SCL_KAM, by = "X.ASV.ID") %>%
  left_join(GIG, by = "X.ASV.ID") %>%
  left_join(PAC, by = "X.ASV.ID") %>%
  left_join(SAC, by = "X.ASV.ID") %>%
  left_join(DIV, by = "X.ASV.ID") %>%
  left_join(ACA, by = "X.ASV.ID") %>%
  left_join(PER, by = "X.ASV.ID") %>%
  left_join(PAR, by = "X.ASV.ID") %>%
  #join with genus level taxonomy 
  left_join(POL_POL, by = "X.ASV.ID") %>%
  left_join(AMB_AMB, by = "X.ASV.ID") %>%
  left_join(ARC_ARC, by = "X.ASV.ID") %>%
  left_join(ENT_ENT, by = "X.ASV.ID") %>%
  left_join(GLO_SCL, by = "X.ASV.ID") %>%
  left_join(GLO_COM, by = "X.ASV.ID") %>%
  left_join(GLO_GLO, by = "X.ASV.ID") %>%
  left_join(SEP_SEP, by = "X.ASV.ID") %>%
  left_join(SEP_FUNG, by = "X.ASV.ID") %>%
  left_join(SEP_FUNF, by = "X.ASV.ID") %>%
  left_join(DOM_NAN, by = "X.ASV.ID") %>%
  left_join(DOM_ORI, by = "X.ASV.ID") %>%
  left_join(DOM_MIC, by = "X.ASV.ID") %>%
  left_join(DOM_MAC, by = "X.ASV.ID") %>%
  left_join(DOM_DOM, by = "X.ASV.ID") %>%
  left_join(SCL_KAM_OEH, by = "X.ASV.ID") %>%
  left_join(SCL_KAM_KAM, by = "X.ASV.ID") %>%
  left_join(SCL_KAM_EPI, by = "X.ASV.ID") %>%
  left_join(SCL_KAM_HAL, by = "X.ASV.ID") %>%
  left_join(SCL_KAM_MIC, by = "X.ASV.ID") %>%
  left_join(SCL_KAM_SCL, by = "X.ASV.ID") %>%
  left_join(SCL_KAM_SIL, by = "X.ASV.ID") %>%
  left_join(SCL_KAM_RHI, by = "X.ASV.ID") %>%
  left_join(GIG_SCU, by = "X.ASV.ID") %>%
  left_join(GIG_BUL, by = "X.ASV.ID") %>%
  left_join(GIG_CET, by = "X.ASV.ID") %>%
  left_join(GIG_RAC, by = "X.ASV.ID") %>%
  left_join(GIG_DEN, by = "X.ASV.ID") %>%
  left_join(GIG_INT, by = "X.ASV.ID") %>%
  left_join(GIG_PAR, by = "X.ASV.ID") %>%
  left_join(GIG_GIG, by = "X.ASV.ID") %>%
  left_join(PAC_PAC, by = "X.ASV.ID") %>%
  left_join(SAC_SAC, by = "X.ASV.ID") %>%
  left_join(DIV_DES, by = "X.ASV.ID") %>%
  left_join(DIV_COR, by = "X.ASV.ID") %>%
  left_join(DIV_RED, by = "X.ASV.ID") %>%
  left_join(DIV_DIV, by = "X.ASV.ID") %>%
  left_join(ACA_ACA, by = "X.ASV.ID") %>%
  left_join(PER_PER, by = "X.ASV.ID") %>%
  left_join(PAR_INN, by = "X.ASV.ID") %>%
  left_join(PAR_PAR, by = "X.ASV.ID") 

# Write to file:
write.table(AMF.tax,file='./ASVtable_clean_AMFonly_taxonomy.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
  