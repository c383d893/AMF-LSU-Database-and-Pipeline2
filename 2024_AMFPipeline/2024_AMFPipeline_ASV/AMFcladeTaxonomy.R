####### Load required packages #######

library(tidyverse)

###### Bring in Glomeromycota (AMF) table ####

AMF <- read.table(file='./cladeOutputs/ASVtable_clean_AMFonly.tsv', header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),phylum = "Glomeromycota")

###### Bring in Family tables #####

POL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPOLonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Polonosporaceae")
AMB <- read.table(file='./cladeOutputs/ASVtable_clean_AMFAMBonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Ambisporaceae")
ARC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFARConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Archaeosporaceae")
ENT <- read.table(file='./cladeOutputs/ASVtable_clean_AMFENTonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Entrophosporaceae")
GLO <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGLOonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Glomeraceae")
SEP <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSEPonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Septoglomeraceae")
DOM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Dominikiaceae")
SCL_KAM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Sclerocystaceae_Kamienskiaceae")
GIG <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIGonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Gigasporaceae")
PAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Pacisporaceae")
SAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Sacculosporaceae") 
DIV <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIVonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Diversisporaceae")
ACA <- read.table(file='./cladeOutputs/ASVtable_clean_AMFACAonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Acaulosporaceae")
PER <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPERonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Pervetustaceae")
PAR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPARonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),family = "Paraglomeraceae")

###### Bring in Genera tables #####

POL_POL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPOL_POLonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Polonospora")
AMB_AMB <- read.table(file='./cladeOutputs/ASVtable_clean_AMFAMB_AMBonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Ambispora")
ARC_ARC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFARC_ARConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Archaeospora")
ENT_ENT <- read.table(file='./cladeOutputs/ASVtable_clean_AMFENT_ENTonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Entrophospora")
GLO_SCL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGLO_SCLonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Sclerocarpum")
GLO_COM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGLO_COMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Complexispora")
GLO_GLO <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGLO_GLOonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Glomus")
SEP_SEP <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSEP_SEPonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Septoglomus")
SEP_FUNG <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSEP_FUNGonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Funneliglomus")
SEP_FUNF <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSEP_FUNFonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Funneliformis")
DOM_NAN <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_NANonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Nanoglomus")
DOM_ORI <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_ORIonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Orientoglomus")
DOM_MIC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_MIConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Microdominikia")
DOM_MAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_MAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Macrodominikia")
DOM_DOM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDOM_DOMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Dominikia")
SCL_KAM_OEH <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_OEHonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Oehlia")
SCL_KAM_KAM <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_KAMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Kamienskia")
SCL_KAM_EPI <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_EPIonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Epigeocarpum")
SCL_KAM_HAL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_HALonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Halonatospora")
SCL_KAM_MIC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_SCLonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Microkamienskia")
SCL_KAM_SCL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_KAMonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Sclerocystis")
SCL_KAM_SIL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_SILonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Silvaspora")
SCL_KAM_RHI <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSCL_KAM_RHIonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Rhizophagus")
GIG_SCU <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_SCUonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Scutellospora")
GIG_BUL <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_BULonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Bulbospora")
GIG_CET <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_CETonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Cetraspora") 
GIG_RAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_RAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Racocetra") 
GIG_DEN <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_DENonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Dentiscutata")
GIG_INT <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_INTonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Intraornatospora")
GIG_PAR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_PARonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Paradentiscutata")
GIG_GIG <- read.table(file='./cladeOutputs/ASVtable_clean_AMFGIG_GIGonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Gigaspora")
PAC_PAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPAC_PAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Pacispora") 
SAC_SAC <- read.table(file='./cladeOutputs/ASVtable_clean_AMFSAC_SAConly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Sacculospora")  
DIV_DES <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_DESonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Desertispora") 
DIV_PAR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_PARonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Paracorymbiglomus") 
DIV_COR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_CORonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Corymbiglomus") 
DIV_RED <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_REDonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Redeckera") 
DIV_DIV <- read.table(file='./cladeOutputs/ASVtable_clean_AMFDIV_DIVonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Diversispora") 
ACA_ACA <- read.table(file='./cladeOutputs/ASVtable_clean_AMFACA_ACAonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Acaulospora") 
PER_PER <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPER_PERonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Pervetustus") 
PAR_INN <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPAR_INNonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Innospora") 
PAR_PAR <- read.table(file='./cladeOutputs/ASVtable_clean_AMFPAR_PARonly.tsv',header = TRUE, sep = '\t') %>% select(X.ASV.ID) %>% mutate(X.ASV.ID = as.character(X.ASV.ID),genus = "Paraglomus") 

###### Assign taxonomy to AMF #####

AMF.tax <- merge(
  AMF,
  bind_rows(POL, AMB, ARC, ENT, GLO, SEP, DOM, SCL_KAM, GIG, PAC,SAC, DIV, ACA, PER, PAR),
  by = "X.ASV.ID", all.x = TRUE
)

AMF.tax <- merge(
  AMF.tax,
  bind_rows(
    POL_POL, AMB_AMB, ARC_ARC, ENT_ENT, GLO_SCL, GLO_COM, GLO_GLO, SEP_SEP, SEP_FUNG, SEP_FUNF, DOM_NAN, 
    DOM_ORI, DOM_MIC, DOM_MAC, DOM_DOM, SCL_KAM_OEH, SCL_KAM_KAM, SCL_KAM_EPI, SCL_KAM_HAL, SCL_KAM_MIC, 
    SCL_KAM_SCL, SCL_KAM_SIL, SCL_KAM_RHI, GIG_SCU, GIG_BUL, GIG_CET, GIG_RAC, GIG_DEN, GIG_INT, GIG_PAR, 
    GIG_GIG, PAC_PAC, SAC_SAC, DIV_DES, DIV_PAR, DIV_COR, DIV_RED, DIV_DIV, ACA_ACA, PER_PER, PAR_INN, PAR_PAR
  ),
  by = "X.ASV.ID", all.x = TRUE
)

# Write to file:
write.table(AMF.tax,file='./ASVtable_clean_AMFonly_taxonomy.tsv',sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
  